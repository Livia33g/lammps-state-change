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

#include "fix_state_change_ksat.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "random_park.h"
#include "update.h"

#include <mpi.h>
#include <cmath>
#include <cstring>
#include <map>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixStateChangeKsat::FixStateChangeKsat(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), last_change(nullptr), effective_type(nullptr),
    prev_coord_type2(nullptr), patch_mask(0)
{
  if (narg < 7) error->all(FLERR, "Illegal fix state/change/ksat command");

  // Parse arguments:
  // fix ID group state/change/ksat check_every cooldown_steps probability cutoff group_patches
  
  check_every = utils::inumeric(FLERR, arg[3], false, lmp);
  cooldown_steps = utils::inumeric(FLERR, arg[4], false, lmp);
  probability = utils::numeric(FLERR, arg[5], false, lmp);
  cutoff = utils::numeric(FLERR, arg[6], false, lmp);
  
  // Get group ID for patches
  group_patches = group->find(arg[7]);
  if (group_patches < 0) error->all(FLERR, "Fix state/change/ksat group_patches not found");

  if (check_every <= 0) error->all(FLERR, "Illegal fix state/change/ksat check_every value");
  if (cooldown_steps < 0) error->all(FLERR, "Illegal fix state/change/ksat cooldown_steps value");
  if (probability < 0.0 || probability > 1.0)
    error->all(FLERR, "Illegal fix state/change/ksat probability value");
  if (cutoff <= 0.0) error->all(FLERR, "Illegal fix state/change/ksat cutoff value");

  // Per-atom arrays
  peratom_flag = 1;
  size_peratom_cols = 0;
  peratom_freq = 1;
  
  // Restart support
  restart_global = 1;
  restart_peratom = 1;
  
  // Initialize
  next_check = update->ntimestep + check_every;
  nchanges = 0;
  nattempts = 0;
  seed = 12345;

  // Provide scalar output for thermo (number of successful state changes)
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;
  
  // Allocate per-atom arrays (match octahedron pattern exactly)
  grow_arrays(atom->nmax);
  
  // Initialize arrays to safe defaults - actual initialization happens in init()
  // This avoids accessing group->bitmask in constructor when groups might not be ready
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    effective_type[i] = atom->type[i];  // Default to current type
    last_change[i] = -1;
    prev_coord_type2[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

FixStateChangeKsat::~FixStateChangeKsat()
{
  memory->destroy(last_change);
  memory->destroy(effective_type);
  memory->destroy(prev_coord_type2);
}

/* ---------------------------------------------------------------------- */

int FixStateChangeKsat::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStateChangeKsat::init()
{
  // No compute dependency - we use direct distance calculation
  // Initialize next_check (like octahedron)
  if (next_check < update->ntimestep) next_check = update->ntimestep + check_every;
  
  // CRITICAL: Ensure arrays are properly sized
  if (atom->nmax > 0) {
    grow_arrays(atom->nmax);
  }
  
  // CRITICAL: Now initialize effective_type properly - groups are ready in init()
  // Type 2: never changes (monomer A)
  // Type 3: initial state (monomer B) - can change to 4
  // Type 4: final state (monomer B) - never changes
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  // Safety check: ensure group_patches is valid
  if (group_patches < 0 || group_patches >= group->ngroup) {
    error->all(FLERR, "Fix state/change/ksat: Invalid group_patches in init()");
  }
  
  // Cache the patch bitmask for safe access later
  patch_mask = group->bitmask[group_patches];
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (mask[i] & patch_mask) {
        // This is a patch
        if (type[i] == 2 || type[i] == 3 || type[i] == 4) {
          effective_type[i] = type[i];  // Direct mapping for types 2, 3, 4
        } else {
          // Unknown patch type - default to 3 (can change to 4)
          effective_type[i] = 3;
        }
        last_change[i] = -1;  // -1 means ready for change (not in cooldown)
        prev_coord_type2[i] = 0.0;
      } else {
        // Body atom - keep its type
        effective_type[i] = type[i];
        last_change[i] = -1;
        prev_coord_type2[i] = 0.0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixStateChangeKsat::post_force(int /*vflag*/)
{
  // Check if it's time to evaluate state changes
  if (update->ntimestep < next_check) return;
  
  next_check = update->ntimestep + check_every;
  nchanges = 0;
  nattempts = 0;
  
  // Check and perform state changes (but don't update types yet)
  check_and_change();
}

/* ---------------------------------------------------------------------- */

void FixStateChangeKsat::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  // Not implemented for rRESPA
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixStateChangeKsat::end_of_step()
{
  // Update atom types at end of step (after forces are computed)
  if (nchanges > 0) {
    update_atom_types();
  }
}

/* ---------------------------------------------------------------------- */

double FixStateChangeKsat::get_coordination(int i, int target_type)
{
  // Count neighbors of target_type within cutoff
  // Safety checks: ensure atom index is valid and domain is initialized
  if (i < 0 || i >= atom->nmax) return 0.0;
  if (!domain || !atom->x || !atom->x[i]) return 0.0;
  
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  if (!type || !mask || !molecule) return 0.0;
  
  int my_mol = molecule[i];
  int nlocal = atom->nlocal;
  double *prd = domain->prd;
  
  if (!prd) return 0.0;  // Domain not initialized
  
  double coord = 0.0;
  double cutoffsq = cutoff * cutoff;
  
  // Check local atoms for coordination
  for (int j = 0; j < nlocal; j++) {
    if (i == j) continue;
    if (j >= atom->nmax || !atom->x[j]) continue;  // Safety check
    if (molecule[j] == my_mol) continue;  // Skip same molecule
    if (!(mask[j] & patch_mask)) continue;  // Must be a patch
    if (type[j] != target_type) continue;  // Must be target type
    
    // Calculate distance with PBC
    double dx = atom->x[j][0] - atom->x[i][0];
    double dy = atom->x[j][1] - atom->x[i][1];
    double dz = atom->x[j][2] - atom->x[i][2];
    
    // PBC with safety check for zero periodicity
    if (prd[0] > 0.0) dx -= prd[0] * std::round(dx / prd[0]);
    if (prd[1] > 0.0) dy -= prd[1] * std::round(dy / prd[1]);
    if (prd[2] > 0.0) dz -= prd[2] * std::round(dz / prd[2]);
    
    double rsq = dx*dx + dy*dy + dz*dz;
    if (rsq < cutoffsq) {
      coord += 1.0;
    }
  }
  
  return coord;
}

/* ---------------------------------------------------------------------- */

void FixStateChangeKsat::check_and_change()
{
  // Safety check: ensure arrays are allocated
  if (!last_change || !effective_type || !prev_coord_type2) return;
  if (!atom || !atom->type || !atom->mask || !atom->molecule) return;
  if (!group || !group->bitmask) return;
  
  // Safety check: ensure group_patches is valid
  if (group_patches < 0 || group_patches >= group->ngroup) return;
  
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  
  // Safety check: ensure nlocal is valid
  if (nlocal <= 0 || nlocal > atom->nmax) return;
  
  bigint timestep = update->ntimestep;
  
  // Random number generator - use timestep and MPI rank for unique seeds
  int my_seed = seed + static_cast<int>(timestep) + comm->me;
  RanPark *random = new RanPark(lmp, my_seed);
  
  // Group patches by molecule to ensure all patches in molecule change together
  std::map<int, std::vector<int> > mol_patches;
  
  // Use cached patch_mask (set in init())
  // Safety check: ensure patch_mask is valid
  if (patch_mask == 0) {
    // If mask is 0, something went wrong - skip processing
    return;
  }
  
  for (int i = 0; i < nlocal; i++) {
    // Safety check: ensure index is valid
    if (i >= atom->nmax || i < 0) break;
    if (!(mask[i] & groupbit)) continue;
    if (!(mask[i] & patch_mask)) continue;  // Must be a patch
    // Safety check: ensure arrays are valid
    if (!effective_type) continue;
    // Only process type 3 patches (can change to 4)
    // Type 2 and 4 never change
    if (type[i] == 3 || effective_type[i] == 3) {
      int mol_id = molecule[i];
      mol_patches[mol_id].push_back(i);
    }
  }
  
  // Process each molecule's patches together as a single unit
  for (auto &mol_pair : mol_patches) {
    int mol_id = mol_pair.first;
    std::vector<int> &patch_indices = mol_pair.second;
    
    if (patch_indices.empty()) continue;
    
    // Representative patch index (use first patch for status checks)
    int rep_idx = patch_indices[0];
    
    // Safety check: ensure rep_idx is valid
    if (rep_idx < 0 || rep_idx >= atom->nmax) continue;
    
    // Ensure all patches in molecule have same effective_type
    int mol_effective_type = effective_type[rep_idx];
    for (int idx : patch_indices) {
      effective_type[idx] = mol_effective_type;
    }
    
    // Only process type 3 patches (type 2 and 4 never change)
    if (mol_effective_type != 3) continue;
    
    // Check cooldown period
    bool in_cooldown = false;
    if (last_change[rep_idx] >= 0) {
      bigint steps_since_change = timestep - last_change[rep_idx];
      if (steps_since_change < cooldown_steps) {
        in_cooldown = true;
      } else {
        // Cooldown expired
        last_change[rep_idx] = -1;
        for (int idx : patch_indices) {
          last_change[idx] = -1;
        }
      }
    }
    
    if (in_cooldown) continue;
    
    // Calculate coordination with type 2 patches (use maximum across all patches)
    double max_coord_type2 = 0.0;
    for (int idx : patch_indices) {
      double coord_val = get_coordination(idx, 2);  // Get coordination with type 2
      if (coord_val > max_coord_type2) max_coord_type2 = coord_val;
    }
    
    // Get previous coordination
    double max_prev_coord_type2 = prev_coord_type2[rep_idx];
    
    // Check if currently attached to type 2 (threshold >= 0.8 means strong attachment)
    bool currently_attached = (max_coord_type2 >= 0.8);
    
    // Check if was previously attached
    bool was_attached = (max_prev_coord_type2 >= 0.8);
    
    // Detect new attachment (transition from not attached to attached)
    bool new_attachment = !was_attached && currently_attached;
    
    // Update previous coordination
    for (int idx : patch_indices) {
      prev_coord_type2[idx] = max_coord_type2;
    }
    
    // State change logic: Type 3 → 4 when attached to type 2
    bool should_change = false;
    
    if (new_attachment && mol_effective_type == 3) {
      // Type 3 patches attach to type 2 → change to type 4
      double rand_val = random->uniform();
      if (rand_val < probability) {
        should_change = true;
        nattempts++;
      }
    }
    
    // Apply change to ALL patches in molecule simultaneously
    if (should_change) {
      for (int idx : patch_indices) {
        effective_type[idx] = 4;  // Change to type 4
        last_change[idx] = timestep;  // Mark as changed
        nchanges++;
      }
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

void FixStateChangeKsat::update_atom_types()
{
  // Safety checks
  if (!atom || !atom->type || !atom->mask || !atom->molecule) return;
  if (!group || !group->bitmask) return;
  if (group_patches < 0 || group_patches >= group->ngroup) return;
  
  // CRITICAL: Store original atom count before any changes
  bigint natoms_before = atom->natoms;
  
  // CRITICAL: First, synchronize effective_type across all processors
  // This ensures patches on different processors get the same effective_type
  comm->forward_comm(this);
  
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  
  // Safety check: ensure nlocal is valid
  if (nlocal <= 0 || nlocal > atom->nmax) return;
  
  // First pass: ensure all patches in same molecule have same effective_type
  // Group by molecule and use the most common type
  std::map<int, std::vector<int> > mol_patches;
  // Use cached patch_mask (set in init())
  if (patch_mask == 0) return;  // Safety check
  for (int i = 0; i < nlocal; i++) {
    if (i >= atom->nmax || i < 0) break;
    if (!(mask[i] & groupbit)) continue;
    if (!(mask[i] & patch_mask)) continue;
    if (type[i] == 2 || type[i] == 3 || type[i] == 4) {
      int mol_id = molecule[i];
      mol_patches[mol_id].push_back(i);
    }
  }
  
  // For each molecule, find the most common effective_type and apply to all patches
  for (auto &mol_pair : mol_patches) {
    int mol_id = mol_pair.first;
    std::vector<int> &patch_indices = mol_pair.second;
    
    if (patch_indices.empty()) continue;
    
    // Count types in this molecule (local atoms only)
    int count_type2 = 0;
    int count_type3 = 0;
    int count_type4 = 0;
    for (int idx : patch_indices) {
      if (effective_type[idx] == 2) count_type2++;
      else if (effective_type[idx] == 3) count_type3++;
      else if (effective_type[idx] == 4) count_type4++;
    }
    
    // Determine unified type (use majority)
    int unified_type = 2;  // Default
    if (count_type4 > count_type3 && count_type4 > count_type2) unified_type = 4;
    else if (count_type3 > count_type2) unified_type = 3;
    
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
    if (i >= atom->nmax || i < 0) break;
    if (!(mask[i] & groupbit)) continue;
    if (!(mask[i] & patch_mask)) continue;
    if (type[i] == 2 || type[i] == 3 || type[i] == 4) {
      int mol_id = molecule[i];
      mol_patches_update[mol_id].push_back(i);
    }
  }
  
  // For each molecule, ensure all patches have same effective_type, then update atom types
  for (auto &mol_pair : mol_patches_update) {
    int mol_id = mol_pair.first;
    std::vector<int> &patch_indices = mol_pair.second;
    
    if (patch_indices.empty()) continue;
    
    // Find most common effective_type in this molecule
    int count_type2 = 0;
    int count_type3 = 0;
    int count_type4 = 0;
    for (int idx : patch_indices) {
      if (effective_type[idx] == 2) count_type2++;
      else if (effective_type[idx] == 3) count_type3++;
      else if (effective_type[idx] == 4) count_type4++;
    }
    int unified_type = 2;  // Default
    if (count_type4 > count_type3 && count_type4 > count_type2) unified_type = 4;
    else if (count_type3 > count_type2) unified_type = 3;
    
    // Apply unified type to ALL patches in molecule
    for (int idx : patch_indices) {
      effective_type[idx] = unified_type;
      // Update atom type to match effective_type
      // CRITICAL: Only change if type is actually 2, 3, or 4 (safety check)
      if (type[idx] == 2 || type[idx] == 3 || type[idx] == 4) {
        if (type[idx] != unified_type) {
          type[idx] = unified_type;
        }
      }
    }
  }
  
  // CRITICAL: Verify atom count before communication
  if (atom->natoms != natoms_before) {
    error->all(FLERR, "Fix state/change/ksat: Atom count changed during type update - aborting to prevent corruption");
  }
  
  // CRITICAL: Update ghost atoms (for MPI) - this is critical for parallel simulations
  // Use forward_comm() to update all atom properties including types
  comm->forward_comm();
  
  // CRITICAL: Verify atom count after communication
  if (atom->natoms != natoms_before) {
    error->all(FLERR, "Fix state/change/ksat: Atom count changed after communication - possible atom loss");
  }
  
  // Force neighbor list rebuild on next step
  // This ensures pair interactions are recalculated with new types
  neighbor->decide();
}

/* ---------------------------------------------------------------------- */

double FixStateChangeKsat::compute_scalar()
{
  return static_cast<double>(nchanges);
}

/* ---------------------------------------------------------------------- */

void FixStateChangeKsat::grow_arrays(int nmax)
{
  memory->grow(last_change, nmax, "fix_state_change_ksat:last_change");
  memory->grow(effective_type, nmax, "fix_state_change_ksat:effective_type");
  memory->grow(prev_coord_type2, nmax, "fix_state_change_ksat:prev_coord_type2");
}

/* ---------------------------------------------------------------------- */

void FixStateChangeKsat::copy_arrays(int i, int j, int /*delflag*/)
{
  last_change[j] = last_change[i];
  effective_type[j] = effective_type[i];
  prev_coord_type2[j] = prev_coord_type2[i];
}

/* ---------------------------------------------------------------------- */

int FixStateChangeKsat::pack_exchange(int i, double *buf)
{
  int n = 0;
  buf[n++] = static_cast<double>(last_change[i]);
  buf[n++] = static_cast<double>(effective_type[i]);
  buf[n++] = prev_coord_type2[i];
  return n;
}

/* ---------------------------------------------------------------------- */

int FixStateChangeKsat::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  last_change[nlocal] = static_cast<int>(buf[n++]);
  effective_type[nlocal] = static_cast<int>(buf[n++]);
  prev_coord_type2[nlocal] = buf[n++];
  return n;
}

/* ---------------------------------------------------------------------- */

int FixStateChangeKsat::pack_restart(int i, double *buf)
{
  int n = 0;
  buf[n++] = static_cast<double>(last_change[i]);
  buf[n++] = static_cast<double>(effective_type[i]);
  buf[n++] = prev_coord_type2[i];
  return n;
}

/* ---------------------------------------------------------------------- */

void FixStateChangeKsat::unpack_restart(int nlocal, int nth)
{
  // Restart data will be unpacked by restart system
}

/* ---------------------------------------------------------------------- */

int FixStateChangeKsat::size_restart(int /*nlocal*/)
{
  return 3;
}

/* ---------------------------------------------------------------------- */

int FixStateChangeKsat::maxsize_restart()
{
  return 3;
}

/* ---------------------------------------------------------------------- */

int FixStateChangeKsat::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i, j, m = 0;
  
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = static_cast<double>(last_change[j]);
    buf[m++] = static_cast<double>(effective_type[j]);
    buf[m++] = prev_coord_type2[j];
  }
  
  return m;
}

/* ---------------------------------------------------------------------- */

void FixStateChangeKsat::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m = 0, last = first + n;
  
  for (i = first; i < last; i++) {
    last_change[i] = static_cast<int>(buf[m++]);
    effective_type[i] = static_cast<int>(buf[m++]);
    prev_coord_type2[i] = buf[m++];
  }
}

