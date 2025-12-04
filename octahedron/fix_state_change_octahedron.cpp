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

#include "fix_state_change_octahedron.h"

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

FixStateChangeOctahedron::FixStateChangeOctahedron(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), last_change(nullptr), effective_type(nullptr),
    prev_coord(nullptr)
{
  if (narg < 7) error->all(FLERR, "Illegal fix state/change/octahedron command");

  // Parse arguments:
  // fix ID group state/change/octahedron check_every cooldown_steps probability cutoff group_patches
  
  check_every = utils::inumeric(FLERR, arg[3], false, lmp);
  cooldown_steps = utils::inumeric(FLERR, arg[4], false, lmp);
  probability = utils::numeric(FLERR, arg[5], false, lmp);
  cutoff = utils::numeric(FLERR, arg[6], false, lmp);
  
  // Get group ID for patches
  group_patches = group->find(arg[7]);
  if (group_patches < 0) error->all(FLERR, "Fix state/change/octahedron group_patches not found");

  if (check_every <= 0) error->all(FLERR, "Illegal fix state/change/octahedron check_every value");
  if (cooldown_steps < 0) error->all(FLERR, "Illegal fix state/change/octahedron cooldown_steps value");
  if (probability < 0.0 || probability > 1.0)
    error->all(FLERR, "Illegal fix state/change/octahedron probability value");
  if (cutoff <= 0.0) error->all(FLERR, "Illegal fix state/change/octahedron cutoff value");

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
  
  // Allocate per-atom arrays
  grow_arrays(atom->nmax);
  
  // Initialize effective_type to match current atom types
  // Patches should start as type 1
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (mask[i] & group->bitmask[group_patches]) {
        // This is a patch - initialize based on current type
        // Type 2 = initial "patch type 1", Types 3, 4, 5 are evolved states
        if (type[i] == 2 || type[i] == 3 || type[i] == 4 || type[i] == 5) {
          // Map type 2 to effective_type 1 (for state change logic)
          // Types 3, 4, 5 stay as-is
          if (type[i] == 2) {
            effective_type[i] = 1;  // Map LAMMPS type 2 to effective "patch type 1"
          } else {
            effective_type[i] = type[i];
          }
        } else {
          // If not a recognized patch type, default to 1 (initial state)
          effective_type[i] = 1;
        }
        last_change[i] = -1;  // -1 means ready for change (not in cooldown)
        prev_coord[i] = 0.0;
      } else {
        // Body atom - keep its type
        effective_type[i] = type[i];
        last_change[i] = -1;
        prev_coord[i] = 0.0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

FixStateChangeOctahedron::~FixStateChangeOctahedron()
{
  memory->destroy(last_change);
  memory->destroy(effective_type);
  memory->destroy(prev_coord);
}

/* ---------------------------------------------------------------------- */

int FixStateChangeOctahedron::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStateChangeOctahedron::init()
{
  // No compute dependency - we'll use neighbor lists directly
  // Initialize next_check
  if (next_check < update->ntimestep) next_check = update->ntimestep + check_every;
}

/* ---------------------------------------------------------------------- */

void FixStateChangeOctahedron::post_force(int /*vflag*/)
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

void FixStateChangeOctahedron::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  // Not implemented for rRESPA - would need nlevels_respa member
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixStateChangeOctahedron::end_of_step()
{
  // Update atom types at end of step (after forces are computed)
  if (nchanges > 0) {
    update_atom_types();
  }
}

/* ---------------------------------------------------------------------- */

double FixStateChangeOctahedron::get_coordination(int i, int effective_patch_type)
{
  // Count neighbors of effective type within cutoff
  // effective_patch_type: 1 = initial "patch type 1" (LAMMPS type 2), 3/4/5 = evolved (same types)
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int my_mol = molecule[i];
  int nlocal = atom->nlocal;
  double *prd = domain->prd;
  
  double coord = 0.0;
  double cutoffsq = cutoff * cutoff;
  
  // Map effective type to LAMMPS type
  int lammps_target_type;
  if (effective_patch_type == 1) {
    lammps_target_type = 2;  // Map effective type 1 to LAMMPS type 2
  } else {
    lammps_target_type = effective_patch_type;  // Types 3, 4, 5 map directly
  }
  
  // Check local atoms for coordination (use neighbor list or check all local atoms)
  // For rigid bodies with MPI, all atoms of a molecule should be on same processor
  for (int j = 0; j < nlocal; j++) {
    if (i == j) continue;
    if (molecule[j] == my_mol) continue;  // Skip same molecule
    if (!(mask[j] & group->bitmask[group_patches])) continue;  // Must be a patch
    if (type[j] != lammps_target_type) continue;  // Must be target type
    
    // Calculate distance with PBC
    double dx = atom->x[j][0] - atom->x[i][0];
    double dy = atom->x[j][1] - atom->x[i][1];
    double dz = atom->x[j][2] - atom->x[i][2];
    
    dx -= prd[0] * std::round(dx / prd[0]);
    dy -= prd[1] * std::round(dy / prd[1]);
    dz -= prd[2] * std::round(dz / prd[2]);
    
    double rsq = dx*dx + dy*dy + dz*dz;
    if (rsq < cutoffsq) {
      coord += 1.0;
    }
  }
  
  return coord;
}

/* ---------------------------------------------------------------------- */

bool FixStateChangeOctahedron::check_same_type_coordination(int i, int effective_patch_type)
{
  // Check if patch i (of effective type) is coordinated with another patch of same type
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int my_mol = molecule[i];
  int nlocal = atom->nlocal;
  double *prd = domain->prd;
  
  double cutoffsq = cutoff * cutoff;
  
  // Map effective type to LAMMPS type
  int lammps_target_type;
  if (effective_patch_type == 1) {
    lammps_target_type = 2;  // Map effective type 1 to LAMMPS type 2
  } else {
    lammps_target_type = effective_patch_type;  // Types 3, 4, 5 map directly
  }
  
  // Check local atoms for same-type coordination
  for (int j = 0; j < nlocal; j++) {
    if (i == j) continue;
    if (molecule[j] == my_mol) continue;  // Skip same molecule
    if (!(mask[j] & group->bitmask[group_patches])) continue;  // Must be a patch
    if (type[j] != lammps_target_type) continue;  // Must be same type
    
    // Calculate distance with PBC
    double dx = atom->x[j][0] - atom->x[i][0];
    double dy = atom->x[j][1] - atom->x[i][1];
    double dz = atom->x[j][2] - atom->x[i][2];
    
    dx -= prd[0] * std::round(dx / prd[0]);
    dy -= prd[1] * std::round(dy / prd[1]);
    dz -= prd[2] * std::round(dz / prd[2]);
    
    double rsq = dx*dx + dy*dy + dz*dz;
    if (rsq < cutoffsq) {
      return true;  // Found same-type neighbor within cutoff
    }
  }
  
  return false;
}

/* ---------------------------------------------------------------------- */

void FixStateChangeOctahedron::check_and_change()
{
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  bigint timestep = update->ntimestep;
  
  // Random number generator - use timestep and MPI rank for unique seeds
  int my_seed = seed + static_cast<int>(timestep) + comm->me;
  RanPark random(lmp, my_seed);
  
  // Group patches by molecule to ensure all patches in molecule change together
  std::map<int, std::vector<int> > mol_patches;
    for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if (!(mask[i] & group->bitmask[group_patches])) continue;  // Must be a patch
    // Check if it's a patch type (2 = initial "patch type 1", 3, 4, 5 = evolved)
    if (type[i] != 2 && type[i] != 3 && type[i] != 4 && type[i] != 5) continue;
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
    
    // Ensure all patches in molecule have same effective_type
    // Use majority vote or first patch's type
    int mol_effective_type = effective_type[rep_idx];
    for (int idx : patch_indices) {
      effective_type[idx] = mol_effective_type;
    }
    
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
    
    // Calculate coordination for this molecule (check all patches, use maximum)
    double max_coord = 0.0;
    bool has_same_type_coord = false;
    
    for (int idx : patch_indices) {
      // Get coordination with all patch types
      // effective type 1 = initial "patch type 1" (LAMMPS type 2)
      double coord_1 = get_coordination(idx, 1);  // Maps to LAMMPS type 2
      double coord_3 = get_coordination(idx, 3);
      double coord_4 = get_coordination(idx, 4);
      double coord_5 = get_coordination(idx, 5);
      double total_coord = coord_1 + coord_3 + coord_4 + coord_5;
      
      if (total_coord > max_coord) {
        max_coord = total_coord;
      }
      
      // Check for same-type coordination if patch is type 3, 4, or 5
      if (mol_effective_type == 3 || mol_effective_type == 4 || mol_effective_type == 5) {
        if (check_same_type_coordination(idx, mol_effective_type)) {
          has_same_type_coord = true;
        }
      }
    }
    
    // Determine if attached (coordination > 0.5 threshold)
    bool is_attached = (max_coord >= 0.5);
    bool was_attached = (prev_coord[rep_idx] >= 0.5);
    bool new_attachment = !was_attached && is_attached;
    
    bool should_change = false;
    int new_type = mol_effective_type;
    
    // Logic for state changes:
    // 1. Type 1 -> when attaches to any patch, randomly change to 3, 4, or 5
    // 2. Type 3/4/5 -> when attaches to same type, randomly change to 3, 4, or 5
    
    if (mol_effective_type == 1 && new_attachment) {
      // Type 1 attaching - randomly choose 3, 4, or 5
      double rand_val = random.uniform();
      if (rand_val < probability) {
        should_change = true;
        // Randomly choose from {3, 4, 5}
        double type_rand = random.uniform();
        if (type_rand < 0.333333) {
          new_type = 3;
        } else if (type_rand < 0.666666) {
          new_type = 4;
        } else {
          new_type = 5;
        }
        nattempts++;
      }
    } else if ((mol_effective_type == 3 || mol_effective_type == 4 || mol_effective_type == 5) 
               && has_same_type_coord && is_attached) {
      // Type 3/4/5 attaching to same type - randomly change to 3, 4, or 5
      double rand_val = random.uniform();
      if (rand_val < probability) {
        should_change = true;
        // Randomly choose from {3, 4, 5}
        double type_rand = random.uniform();
        if (type_rand < 0.333333) {
          new_type = 3;
        } else if (type_rand < 0.666666) {
          new_type = 4;
        } else {
          new_type = 5;
        }
        nattempts++;
      }
    }
    
    // Apply change to ALL patches in molecule simultaneously
    if (should_change) {
      for (int idx : patch_indices) {
        effective_type[idx] = new_type;
        last_change[idx] = timestep;
        nchanges++;
      }
    }
    
    // Update previous coordination for all patches
    for (int idx : patch_indices) {
      prev_coord[idx] = max_coord;
    }
  }
  
  // Sum changes across all MPI ranks
  MPI_Allreduce(MPI_IN_PLACE, &nchanges, 1, MPI_INT, MPI_SUM, world);
}

/* ---------------------------------------------------------------------- */

void FixStateChangeOctahedron::update_atom_types()
{
  // Update actual atom types based on effective_type
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  
  // Group by molecule and update all patches together
  std::map<int, std::vector<int> > mol_patches;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if (!(mask[i] & group->bitmask[group_patches])) continue;
    if (type[i] != 2 && type[i] != 3 && type[i] != 4 && type[i] != 5) continue;
    int mol_id = molecule[i];
    mol_patches[mol_id].push_back(i);
  }
  
  // For each molecule, ensure all patches have same effective_type, then update
  for (auto &mol_pair : mol_patches) {
    int mol_id = mol_pair.first;
    std::vector<int> &patch_indices = mol_pair.second;
    
    if (patch_indices.empty()) continue;
    
    // Use first patch's effective_type for all patches in molecule
    int unified_effective_type = effective_type[patch_indices[0]];
    for (int idx : patch_indices) {
      effective_type[idx] = unified_effective_type;
      // Map effective type to LAMMPS type
      if (unified_effective_type == 1) {
        type[idx] = 2;  // Map effective type 1 to LAMMPS type 2
      } else {
        type[idx] = unified_effective_type;  // Types 3, 4, 5 map directly
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixStateChangeOctahedron::compute_scalar()
{
  return static_cast<double>(nchanges);
}

/* ---------------------------------------------------------------------- */

void FixStateChangeOctahedron::grow_arrays(int nmax)
{
  memory->grow(last_change, nmax, "fix_state_change_octahedron:last_change");
  memory->grow(effective_type, nmax, "fix_state_change_octahedron:effective_type");
  memory->grow(prev_coord, nmax, "fix_state_change_octahedron:prev_coord");
}

/* ---------------------------------------------------------------------- */

void FixStateChangeOctahedron::copy_arrays(int i, int j, int /*delflag*/)
{
  last_change[j] = last_change[i];
  effective_type[j] = effective_type[i];
  prev_coord[j] = prev_coord[i];
}

/* ---------------------------------------------------------------------- */

int FixStateChangeOctahedron::pack_exchange(int i, double *buf)
{
  int n = 0;
  buf[n++] = static_cast<double>(last_change[i]);
  buf[n++] = static_cast<double>(effective_type[i]);
  buf[n++] = prev_coord[i];
  return n;
}

/* ---------------------------------------------------------------------- */

int FixStateChangeOctahedron::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  last_change[nlocal] = static_cast<int>(buf[n++]);
  effective_type[nlocal] = static_cast<int>(buf[n++]);
  prev_coord[nlocal] = buf[n++];
  return n;
}

/* ---------------------------------------------------------------------- */

int FixStateChangeOctahedron::pack_restart(int i, double *buf)
{
  int n = 0;
  buf[n++] = static_cast<double>(last_change[i]);
  buf[n++] = static_cast<double>(effective_type[i]);
  buf[n++] = prev_coord[i];
  return n;
}

/* ---------------------------------------------------------------------- */

void FixStateChangeOctahedron::unpack_restart(int nlocal, int nth)
{
  // Restart data will be unpacked by restart system
}

/* ---------------------------------------------------------------------- */

int FixStateChangeOctahedron::size_restart(int /*nlocal*/)
{
  return 3;
}

/* ---------------------------------------------------------------------- */

int FixStateChangeOctahedron::maxsize_restart()
{
  return 3;
}

