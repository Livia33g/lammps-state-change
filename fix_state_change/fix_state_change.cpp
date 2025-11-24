/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_state_change.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "random_park.h"
#include "update.h"
#include "utils.h"

#include <mpi.h>

#include <cmath>
#include <cstring>
#include <map>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixStateChange::FixStateChange(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), last_change(nullptr), effective_type(nullptr),
    prev_coord_A(nullptr), prev_coord_B(nullptr),
    c_coord_A(nullptr), c_coord_B(nullptr)
{
  if (narg < 9) error->all(FLERR, "Illegal fix state/change command");

  // Parse arguments:
  // fix ID group state/change check_every cooldown_steps probability cutoff group_A group_B
  
  check_every = utils::inumeric(FLERR, arg[3], false, lmp);
  cooldown_steps = utils::inumeric(FLERR, arg[4], false, lmp);
  probability = utils::numeric(FLERR, arg[5], false, lmp);
  cutoff = utils::numeric(FLERR, arg[6], false, lmp);
  
  // Get group IDs
  group_A = group->find(arg[7]);
  group_B = group->find(arg[8]);
  
  if (group_A < 0) error->all(FLERR, "Fix state/change group_A not found");
  if (group_B < 0) error->all(FLERR, "Fix state/change group_B not found");

  if (check_every <= 0) error->all(FLERR, "Illegal fix state/change check_every value");
  if (cooldown_steps < 0) error->all(FLERR, "Illegal fix state/change cooldown_steps value");
  if (probability < 0.0 || probability > 1.0)
    error->all(FLERR, "Illegal fix state/change probability value");
  if (cutoff <= 0.0) error->all(FLERR, "Illegal fix state/change cutoff value");

  // Per-atom arrays
  peratom_flag = 1;
  size_peratom_cols = 0;
  peratom_freq = 1;
  
  // Restart support
  restart_global = 1;
  restart_peratom = 1;
  
  // Initialize
  nlevels_respa = 0;
  next_check = update->ntimestep + check_every;
  nchanges = 0;
  nattempts = 0;
  seed = 12345;  // Default seed, can be made configurable

  // Provide scalar output for thermo (number of successful state changes)
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;
  
  // Allocate per-atom arrays
  grow_arrays(atom->nmax);
  
  // Initialize effective_type to match current atom types
  // Only patches (type 2 or 3) will have effective_type set
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (type[i] == 2) {
        effective_type[i] = 2;
        last_change[i] = -1;  // -1 means ready for a shot (not attached or detached)
        prev_coord_A[i] = 0.0;
        prev_coord_B[i] = 0.0;
      } else if (type[i] == 3) {
        effective_type[i] = 3;
        last_change[i] = -1;  // -1 means ready for a shot (not attached or detached)
        prev_coord_A[i] = 0.0;
        prev_coord_B[i] = 0.0;
      } else {
        effective_type[i] = type[i];  // Body atoms keep their type
        last_change[i] = -1;
        prev_coord_A[i] = 0.0;
        prev_coord_B[i] = 0.0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

FixStateChange::~FixStateChange()
{
  memory->destroy(last_change);
  memory->destroy(effective_type);
  memory->destroy(prev_coord_A);
  memory->destroy(prev_coord_B);
}

/* ---------------------------------------------------------------------- */

int FixStateChange::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStateChange::init()
{
  // Find computes for coordination numbers
  // These should be created in the input script as:
  // compute cA all coord/atom cutoff <cutoff> group patches_A
  // compute cB all coord/atom cutoff <cutoff> group patches_B
  
  int icompute = modify->find_compute("cA");
  if (icompute < 0) error->all(FLERR, "Fix state/change requires compute cA");
  c_coord_A = modify->compute[icompute];
  
  icompute = modify->find_compute("cB");
  if (icompute < 0) error->all(FLERR, "Fix state/change requires compute cB");
  c_coord_B = modify->compute[icompute];
  
  if (c_coord_A->peratom_flag == 0)
    error->all(FLERR, "Fix state/change compute cA does not calculate per-atom values");
  if (c_coord_B->peratom_flag == 0)
    error->all(FLERR, "Fix state/change compute cB does not calculate per-atom values");
  
  if (strcmp(c_coord_A->style, "coord/atom") != 0)
    error->all(FLERR, "Fix state/change compute cA is not coord/atom style");
  if (strcmp(c_coord_B->style, "coord/atom") != 0)
    error->all(FLERR, "Fix state/change compute cB is not coord/atom style");
  
  // Initialize next_check
  if (next_check < update->ntimestep) next_check = update->ntimestep + check_every;
}

/* ---------------------------------------------------------------------- */

void FixStateChange::post_force(int /*vflag*/)
{
  // Check if it's time to evaluate state changes
  if (update->ntimestep < next_check) return;
  
  next_check = update->ntimestep + check_every;
  nchanges = 0;
  nattempts = 0;
  
  // Invoke computes to get coordination numbers
  c_coord_A->compute_peratom();
  c_coord_B->compute_peratom();
  
  // Check and perform state changes (but don't update types yet)
  // We'll update types at end_of_step to avoid breaking rigid bodies
  check_and_change();
  
  // Don't update types here - do it at end_of_step instead
}

/* ---------------------------------------------------------------------- */

void FixStateChange::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa - 1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixStateChange::end_of_step()
{
  // Update atom types at end of step (after forces are computed)
  // This is safer for rigid bodies - avoids breaking structure during force calculation
  if (nchanges > 0) {
    update_atom_types();
  }
}

/* ---------------------------------------------------------------------- */

void FixStateChange::check_and_change()
{
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  bigint timestep = update->ntimestep;
  
  double *coord_A = c_coord_A->vector_atom;
  double *coord_B = c_coord_B->vector_atom;
  
  // CRITICAL: Synchronize effective_type BEFORE processing
  // This ensures we start with consistent state across all processors
  comm->forward_comm(this);
  
  // Random number generator - use timestep and MPI rank for unique seeds
  int my_seed = seed + static_cast<int>(timestep) + comm->me;
  RanPark *random = new RanPark(lmp, my_seed);
  
  // First pass: Group patches by molecule to ensure unified behavior
  // This ensures all patches in a molecule are processed together as one unit
  std::map<int, std::vector<int> > mol_patches;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if (type[i] != 2 && type[i] != 3) continue;
    int mol_id = molecule[i];
    mol_patches[mol_id].push_back(i);
  }
  
  // Process each molecule's patches together as a single unit
  for (auto &mol_pair : mol_patches) {
    int mol_id = mol_pair.first;
    std::vector<int> &patch_indices = mol_pair.second;
    
    if (patch_indices.empty()) continue;
    
    // Representative patch index (use first patch for status checks)
    int rep_idx = patch_indices[0];
    
    // CRITICAL: First, ensure all patches in molecule have same effective_type
    // This prevents drift from previous steps
    int count_type2 = 0;
    int count_type3 = 0;
    for (int idx : patch_indices) {
      if (effective_type[idx] == 2) count_type2++;
      else if (effective_type[idx] == 3) count_type3++;
    }
    int mol_effective_type = (count_type3 > count_type2) ? 3 : 2;
    
    // Apply unified type to ALL patches in molecule FIRST
    for (int idx : patch_indices) {
      effective_type[idx] = mol_effective_type;
    }
    
    // Find MAXIMUM coordination across all patches in molecule
    // If ANY patch is attached, the whole molecule is considered attached
    double max_coord_A = 0.0;
    double max_coord_B = 0.0;
    double max_prev_coord_A = 0.0;
    double max_prev_coord_B = 0.0;
    
    for (int idx : patch_indices) {
      if (coord_A[idx] > max_coord_A) max_coord_A = coord_A[idx];
      if (coord_B[idx] > max_coord_B) max_coord_B = coord_B[idx];
      if (prev_coord_A[idx] > max_prev_coord_A) max_prev_coord_A = prev_coord_A[idx];
      if (prev_coord_B[idx] > max_prev_coord_B) max_prev_coord_B = prev_coord_B[idx];
    }
    
    // Determine if molecule is currently attached (use maximum coordination)
    // HIGHER threshold (0.8) ensures only FULL dimer attachment triggers state change
    // This requires patches to be very close/overlapping, not just nearby
    bool currently_attached = false;
    if (mol_effective_type == 2 && max_coord_A >= 0.8) {
      currently_attached = true;
    } else if (mol_effective_type == 3 && max_coord_B >= 0.8) {
      currently_attached = true;
    }
    
    // ADDITIONAL CHECK: Only allow state change if molecules are in dimer (COM distance < 3.0)
    // This prevents state changes when patches are close but molecules are far apart
    const double DIMER_COM_CUTOFF = 3.0;  // COM distance threshold for dimers
    bool in_dimer = false;
    
    if (currently_attached) {
      // Find the patch with maximum coordination (the one that triggered attachment)
      int max_coord_patch_idx = patch_indices[0];
      double max_coord_val = (max_coord_A > max_coord_B) ? max_coord_A : max_coord_B;
      for (int idx : patch_indices) {
        double coord_val = (coord_A[idx] > coord_B[idx]) ? coord_A[idx] : coord_B[idx];
        if (coord_val > max_coord_val) {
          max_coord_val = coord_val;
          max_coord_patch_idx = idx;
        }
      }
      
      // Calculate COM for this molecule (using all atoms in molecule)
      // CRITICAL: For rigid bodies, all atoms in a molecule are on the same processor
      // So using nlocal is safe - we don't need ghost atoms
      double com_x = 0.0, com_y = 0.0, com_z = 0.0;
      int com_count = 0;
      for (int i = 0; i < nlocal; i++) {
        if (molecule[i] == mol_id) {
          // Safety check: ensure atom index is valid
          if (i >= atom->nmax) {
            error->one(FLERR, "Fix state/change: Invalid atom index in COM calculation");
          }
          com_x += atom->x[i][0];
          com_y += atom->x[i][1];
          com_z += atom->x[i][2];
          com_count++;
        }
      }
      if (com_count == 0) {
        // No atoms found for this molecule - skip COM check
        in_dimer = false;
        currently_attached = false;
        // Skip rest of processing for this molecule by breaking out of the if block
        // The rest of the logic will handle currently_attached = false correctly
      } else {
        com_x /= com_count;
        com_y /= com_count;
        com_z /= com_count;
        
        // Find coordinated patches from other molecules and check COM distance
        double *prd = domain->prd;
        double min_com_distance = 1e10;
        
        // Check all patches to find ones coordinated with this molecule's patches
        for (int j = 0; j < nlocal; j++) {
        if (molecule[j] == mol_id) continue;  // Skip same molecule
        if (!(mask[j] & groupbit)) continue;
        if (type[j] != 2 && type[j] != 3) continue;
        
        // Check if this patch is coordinated with any of our patches
        bool is_coordinated = false;
        if (mol_effective_type == 2 && type[j] == 2 && coord_A[j] >= 0.8) {
          is_coordinated = true;
        } else if (mol_effective_type == 3 && type[j] == 3 && coord_B[j] >= 0.8) {
          is_coordinated = true;
        }
        
        if (is_coordinated) {
          int other_mol_id = molecule[j];
          
          // Calculate COM for other molecule
          // CRITICAL: For rigid bodies, all atoms in a molecule are on the same processor
          double other_com_x = 0.0, other_com_y = 0.0, other_com_z = 0.0;
          int other_com_count = 0;
          for (int k = 0; k < nlocal; k++) {
            if (molecule[k] == other_mol_id) {
              // Safety check: ensure atom index is valid
              if (k >= atom->nmax) {
                error->one(FLERR, "Fix state/change: Invalid atom index in other COM calculation");
              }
              other_com_x += atom->x[k][0];
              other_com_y += atom->x[k][1];
              other_com_z += atom->x[k][2];
              other_com_count++;
            }
          }
          if (other_com_count == 0) {
            // Other molecule not found on this processor - skip this coordination check
            continue;
          }
          other_com_x /= other_com_count;
          other_com_y /= other_com_count;
          other_com_z /= other_com_count;
          
          // Calculate minimum image distance (PBC)
          double dx = other_com_x - com_x;
          double dy = other_com_y - com_y;
          double dz = other_com_z - com_z;
          
          dx -= prd[0] * std::round(dx / prd[0]);
          dy -= prd[1] * std::round(dy / prd[1]);
          dz -= prd[2] * std::round(dz / prd[2]);
          
          double dist = sqrt(dx*dx + dy*dy + dz*dz);
          if (dist < min_com_distance) {
            min_com_distance = dist;
          }
        }  // End of if (is_coordinated)
      }  // End of for (int j = 0; j < nlocal; j++)
        
        // Only consider attached if COM distance < 3.0 (dimer threshold)
        in_dimer = (min_com_distance < DIMER_COM_CUTOFF);
        currently_attached = currently_attached && in_dimer;
      }  // End of else block (com_count > 0)
    }  // End of if (currently_attached)
    
    // Check previous attachment state (use maximum previous coordination)
    bool was_attached = (max_prev_coord_A >= 0.8) || (max_prev_coord_B >= 0.8);
    bool new_attachment = !was_attached && currently_attached;
    bool detached = was_attached && !currently_attached;
    
    double current_coord = (max_coord_A > max_coord_B) ? max_coord_A : max_coord_B;
    
    // Reset logic for all patches in molecule
    if (detached) {
      for (int idx : patch_indices) {
        last_change[idx] = -1;
      }
    }
    
    if (last_change[rep_idx] != -1 && current_coord < 0.35) {
      for (int idx : patch_indices) {
        last_change[idx] = -1;
      }
    }
    
    // Check if molecule should try to switch (use representative patch's status)
    bool ready_for_shot = (last_change[rep_idx] == -1);
    bool still_together_after_change = (last_change[rep_idx] != -1) && (current_coord >= 0.35);
    bool should_try_switch = new_attachment && ready_for_shot && !still_together_after_change;
    
    bool should_change = false;
    int new_effective_type = mol_effective_type;
    
    if (should_try_switch) {
      // Give molecule ONE chance to switch (all patches switch together)
      if (mol_effective_type == 2) {
        double rand_val = random->uniform();
        if (rand_val < probability) {
          should_change = true;
          new_effective_type = 3;
        }
      } else if (mol_effective_type == 3) {
        double rand_val = random->uniform();
        if (rand_val < probability) {
          should_change = true;
          new_effective_type = 2;
        }
      }
      
      // Mark ALL patches in molecule as having had their shot
      for (int idx : patch_indices) {
        last_change[idx] = timestep;
      }
      nattempts++;
    }
    
    // Apply change to ALL patches in molecule simultaneously
    if (should_change) {
      for (int idx : patch_indices) {
        effective_type[idx] = new_effective_type;
        nchanges++;
      }
    }
    
    // Update previous coordination for all patches in molecule (use maximum)
    for (int idx : patch_indices) {
      prev_coord_A[idx] = max_coord_A;
      prev_coord_B[idx] = max_coord_B;
    }
  }
  
  delete random;
  
  // CRITICAL: Synchronize effective_type across all processors BEFORE summing
  // This ensures patches on different processors have consistent types
  comm->forward_comm(this);
  
  // Sum changes across all MPI ranks
  MPI_Allreduce(MPI_IN_PLACE, &nchanges, 1, MPI_INT, MPI_SUM, world);
}

/* ---------------------------------------------------------------------- */

void FixStateChange::update_atom_types()
{
  // Update actual atom types based on effective_type
  // This is done directly without unfixing rigid bodies
  // All patches in the same molecule change together (handled in check_and_change)
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  
  // CRITICAL: Store original atom count before any changes
  bigint natoms_before = atom->natoms;
  
  // CRITICAL: First, synchronize effective_type across all processors
  // This ensures patches on different processors get the same effective_type
  comm->forward_comm(this);
  
  // First pass: ensure all patches in same molecule have same effective_type
  // Group by molecule and use the most common type
  std::map<int, std::vector<int> > mol_patches;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if (type[i] != 2 && type[i] != 3) continue;
    int mol_id = molecule[i];
    mol_patches[mol_id].push_back(i);
  }
  
  // For each molecule, find the most common effective_type and apply to all patches
  for (auto &mol_pair : mol_patches) {
    int mol_id = mol_pair.first;
    std::vector<int> &patch_indices = mol_pair.second;
    
    if (patch_indices.empty()) continue;
    
    // Count types in this molecule (local atoms only)
    int count_type2 = 0;
    int count_type3 = 0;
    for (int idx : patch_indices) {
      if (effective_type[idx] == 2) count_type2++;
      else if (effective_type[idx] == 3) count_type3++;
    }
    
    // Determine unified type (use majority, or type 2 if tie)
    // For rigid bodies, all patches should be on same processor, so local count should work
    int unified_type = (count_type3 > count_type2) ? 3 : 2;
    
    // Apply unified type to ALL patches in molecule (local only)
    for (int idx : patch_indices) {
      effective_type[idx] = unified_type;
    }
  }
  
  // Synchronize again after unification to ensure consistency across processors
  comm->forward_comm(this);
  
  // Second pass: update actual types
  // CRITICAL: Group by molecule and update all patches together
  std::map<int, std::vector<int> > mol_patches_update;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if (type[i] != 2 && type[i] != 3) continue;
    int mol_id = molecule[i];
    mol_patches_update[mol_id].push_back(i);
  }
  
  // For each molecule, ensure all patches have same effective_type, then update atom types
  for (auto &mol_pair : mol_patches_update) {
    int mol_id = mol_pair.first;
    std::vector<int> &patch_indices = mol_pair.second;
    
    if (patch_indices.empty()) continue;
    
    // Find most common effective_type in this molecule
    int count_type2 = 0;
    int count_type3 = 0;
    for (int idx : patch_indices) {
      if (effective_type[idx] == 2) count_type2++;
      else if (effective_type[idx] == 3) count_type3++;
    }
    int unified_type = (count_type3 > count_type2) ? 3 : 2;
    
    // Apply unified type to ALL patches in molecule
    for (int idx : patch_indices) {
      effective_type[idx] = unified_type;
      // Update atom type to match effective_type
      // CRITICAL: Only change if type is actually 2 or 3 (safety check)
      if (type[idx] == 2 || type[idx] == 3) {
        if (type[idx] != unified_type) {
          type[idx] = unified_type;
        }
      }
    }
  }
  
  // CRITICAL: Verify atom count before communication
  if (atom->natoms != natoms_before) {
    error->all(FLERR, "Fix state/change: Atom count changed during type update - aborting to prevent corruption");
  }
  
  // Don't call domain->pbc() or domain->reset_box() - these can break rigid bodies!
  // The rigid fix handles PBC automatically
  
  // CRITICAL: Update ghost atoms (for MPI) - this is critical for parallel simulations
  // Use forward_comm() to update all atom properties including types
  comm->forward_comm();
  
  // CRITICAL: Verify atom count after communication
  if (atom->natoms != natoms_before) {
    error->all(FLERR, "Fix state/change: Atom count changed after communication - possible atom loss");
  }
  
  // Force neighbor list rebuild on next step
  // This ensures pair interactions are recalculated with new types
  neighbor->decide();
  
  // Final safety check: verify atoms are still present
  bigint natoms_after = atom->natoms;
  int natoms_local = atom->nlocal;
  int natoms_ghost = atom->nghost;
  
  if (natoms_after != natoms_before) {
    char msg[256];
    sprintf(msg, "Fix state/change: Atom loss detected! Before: %ld, After: %ld, Local: %d, Ghost: %d",
            natoms_before, natoms_after, natoms_local, natoms_ghost);
    error->all(FLERR, msg);
  }
}

/* ---------------------------------------------------------------------- */

double FixStateChange::compute_scalar()
{
  return static_cast<double>(nchanges);
}

/* ---------------------------------------------------------------------- */

void FixStateChange::grow_arrays(int nmax)
{
  memory->grow(last_change, nmax, "fix_state_change:last_change");
  memory->grow(effective_type, nmax, "fix_state_change:effective_type");
  memory->grow(prev_coord_A, nmax, "fix_state_change:prev_coord_A");
  memory->grow(prev_coord_B, nmax, "fix_state_change:prev_coord_B");
}

/* ---------------------------------------------------------------------- */

void FixStateChange::copy_arrays(int i, int j, int /*delflag*/)
{
  last_change[j] = last_change[i];
  effective_type[j] = effective_type[i];
  prev_coord_A[j] = prev_coord_A[i];
  prev_coord_B[j] = prev_coord_B[i];
}

/* ---------------------------------------------------------------------- */

int FixStateChange::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i, j, m = 0;
  
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = static_cast<double>(last_change[j]);
    buf[m++] = static_cast<double>(effective_type[j]);
    buf[m++] = prev_coord_A[j];
    buf[m++] = prev_coord_B[j];
  }
  
  return m;
}

/* ---------------------------------------------------------------------- */

void FixStateChange::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m = 0, last = first + n;
  
  for (i = first; i < last; i++) {
    last_change[i] = static_cast<int>(buf[m++]);
    effective_type[i] = static_cast<int>(buf[m++]);
    prev_coord_A[i] = buf[m++];
    prev_coord_B[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int FixStateChange::pack_exchange(int i, double *buf)
{
  int m = 0;
  buf[m++] = static_cast<double>(last_change[i]);
  buf[m++] = static_cast<double>(effective_type[i]);
  buf[m++] = prev_coord_A[i];
  buf[m++] = prev_coord_B[i];
  return m;
}

/* ---------------------------------------------------------------------- */

int FixStateChange::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  last_change[nlocal] = static_cast<int>(buf[m++]);
  effective_type[nlocal] = static_cast<int>(buf[m++]);
  prev_coord_A[nlocal] = buf[m++];
  prev_coord_B[nlocal] = buf[m++];
  return m;
}

/* ---------------------------------------------------------------------- */

int FixStateChange::pack_restart(int i, double *buf)
{
  int m = 0;
  buf[m++] = 4;  // Number of values (last_change, effective_type, prev_coord_A, prev_coord_B)
  buf[m++] = static_cast<double>(last_change[i]);
  buf[m++] = static_cast<double>(effective_type[i]);
  buf[m++] = prev_coord_A[i];
  buf[m++] = prev_coord_B[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void FixStateChange::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;
  
  int m = 0;
  int n = static_cast<int>(extra[nlocal][nth]);
  last_change[nlocal] = static_cast<int>(extra[nlocal][m++]);
  effective_type[nlocal] = static_cast<int>(extra[nlocal][m++]);
  prev_coord_A[nlocal] = extra[nlocal][m++];
  prev_coord_B[nlocal] = extra[nlocal][m++];
}

/* ---------------------------------------------------------------------- */

int FixStateChange::size_restart(int /*nlocal*/)
{
  return 5;  // count + 4 values (last_change, effective_type, prev_coord_A, prev_coord_B)
}

/* ---------------------------------------------------------------------- */

int FixStateChange::maxsize_restart()
{
  return 5;
}

