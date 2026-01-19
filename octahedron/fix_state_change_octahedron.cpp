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
#include <set>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixStateChangeOctahedron::FixStateChangeOctahedron(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), last_change(nullptr), cooldown_duration(nullptr),
    effective_type(nullptr), prev_coord(nullptr), contact_timer(nullptr)
{
  if (narg < 7) error->all(FLERR, "Illegal fix state/change/octahedron command");

  // Parse arguments:
  // fix ID group state/change/octahedron check_every cooldown_steps probability cutoff group_patches [hysteresis_threshold]
  
  check_every = utils::inumeric(FLERR, arg[3], false, lmp);
  cooldown_steps = utils::inumeric(FLERR, arg[4], false, lmp);
  probability = utils::numeric(FLERR, arg[5], false, lmp);
  cutoff = utils::numeric(FLERR, arg[6], false, lmp);
  
  // Get group ID for patches
  group_patches = group->find(arg[7]);
  if (group_patches < 0) error->all(FLERR, "Fix state/change/octahedron group_patches not found");
  
  // Parse hysteresis threshold (optional, default to 1000 if not provided)
  if (narg > 8) {
    hysteresis_threshold = utils::inumeric(FLERR, arg[8], false, lmp);
  } else {
    hysteresis_threshold = 1000;  // Default wait time for hysteresis
  }
  if (hysteresis_threshold < 0) error->all(FLERR, "Illegal fix state/change/octahedron hysteresis_threshold value");

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
  
  // Communication: We send 4 values per atom (last_change, cooldown_duration, effective_type, prev_coord)
  comm_forward = 4;
  
  // Initialize
  next_check = update->ntimestep + check_every;
  nchanges = 0;
  nattempts = 0;
  seed = 12345;
  max_changes_per_step = 1;  // CRITICAL: Global rate limiter - REDUCED to max 1 molecule per timestep (prevents mass behavior more aggressively)
  
  // Debug counters
  nconfident_contacts = 0;
  ntrigger_attempts = 0;
  ncooldown_blocked = 0;
  
  // Consistency sweep diagnostic counters
  nsweep_attempts = 0;
  nsweep_blocked_cooldown = 0;
  nsweep_blocked_timestamp = 0;
  nsweep_applied = 0;
  
  // Provide scalar output for thermo (number of successful state changes)
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;
  
  // Allocate per-atom arrays
  grow_arrays(atom->nmax);
  
  // Initialize arrays to safe defaults - actual initialization happens in init()
  // This avoids accessing group->bitmask in constructor when groups might not be ready
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    effective_type[i] = atom->type[i];  // Default to current type
    last_change[i] = -1;
    prev_coord[i] = 0.0;
    contact_timer[i] = 0;  // Initialize contact timer to 0
  }
}

/* ---------------------------------------------------------------------- */

FixStateChangeOctahedron::~FixStateChangeOctahedron()
{
  memory->destroy(last_change);
  memory->destroy(cooldown_duration);
  memory->destroy(effective_type);
  memory->destroy(prev_coord);
  memory->destroy(contact_timer);
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
  
  // CRITICAL: Ensure arrays are properly sized
  if (atom->nmax > 0) {
    grow_arrays(atom->nmax);
  }
  
  // CRITICAL: Now initialize effective_type properly - groups are ready in init()
  // Patches should start as type 1
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  // Safety check: ensure group_patches is valid
  if (group_patches < 0 || group_patches >= group->ngroup) {
    error->all(FLERR, "Fix state/change/octahedron: Invalid group_patches in init()");
  }
  
  int patch_mask = group->bitmask[group_patches];
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (mask[i] & patch_mask) {
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
        cooldown_duration[i] = cooldown_steps;  // Default, will be randomized on first change
        prev_coord[i] = 0.0;
        contact_timer[i] = 0;  // Initialize contact timer to 0
      } else {
        // Body atom - keep its type
        effective_type[i] = type[i];
        last_change[i] = -1;
        cooldown_duration[i] = cooldown_steps;  // Default, will be randomized on first change
        prev_coord[i] = 0.0;
        contact_timer[i] = 0;
      }
    }
  }
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
    
    // Safety check: ensure atom index is valid
    if (j >= atom->nmax || !atom->x[j]) continue;
    
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

bool FixStateChangeOctahedron::check_same_type_coordination(int i, int effective_patch_type)
{
  // Check if patch i (of effective type) is coordinated with another patch of same type
  // Safety checks
  if (i < 0 || i >= atom->nmax) return false;
  if (!domain || !atom->x || !atom->x[i]) return false;
  
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  if (!type || !mask || !molecule) return false;
  
  int my_mol = molecule[i];
  int nlocal = atom->nlocal;
  double *prd = domain->prd;
  
  if (!prd) return false;  // Domain not initialized
  
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
    
    // Safety check: ensure atom index is valid
    if (j >= atom->nmax || !atom->x[j]) continue;
    
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
      return true;  // Found same-type neighbor within cutoff
    }
  }
  
  return false;
}

/* ---------------------------------------------------------------------- */

void FixStateChangeOctahedron::check_and_change()
{
  // CRITICAL FIX: Update ghost atoms immediately before conflict checks
  // This ensures all processors have the latest effective_type and last_change values
  // from the previous step, preventing "phantom" conflicts from stale ghost data
  comm->forward_comm(this);
  
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule; // Need this to find "brothers" (all patches in same monomer)
  int nlocal = atom->nlocal;
  bigint timestep = update->ntimestep;

  // Map to store molecule changes to avoid double-processing
  // key = mol_id, value = new_type
  std::map<int, int> mol_changes;

  // 1. SCAN LOOP: Look for triggers (independent per-patch detection)
  for (int i = 0; i < nlocal; i++) {
    // CRITICAL FIX: Create random generator per atom with molecule ID in seed
    // This ensures molecules with different IDs get different random sequences,
    // breaking correlation even if processed on the same timestep
    int my_seed = seed + static_cast<int>(timestep) + comm->me + molecule[i];
    RanPark random(lmp, my_seed);
    // Basic filters: Must be in group, must be a patch type
    if (!(mask[i] & groupbit)) continue;
    if (!(mask[i] & group->bitmask[group_patches])) continue;
    if (type[i] != 2 && type[i] != 3 && type[i] != 4 && type[i] != 5) continue;
    
    // Skip if this molecule is already scheduled to change in this step
    if (mol_changes.find(molecule[i]) != mol_changes.end()) continue;

    // Cooldown check using PER-ATOM randomized cooldown duration (jitter breaks synchronization)
    if (last_change[i] >= 0) {
        bigint steps_since_change = timestep - last_change[i];
        // Use per-atom randomized cooldown instead of fixed cooldown_steps
        if (steps_since_change < cooldown_duration[i]) continue;
    }

    int my_eff_type = effective_type[i];
    
    // Compute coordination for THIS patch only
    double coord_1 = get_coordination(i, 1);
    double coord_3 = get_coordination(i, 3);
    double coord_4 = get_coordination(i, 4);
    double coord_5 = get_coordination(i, 5);
    double total_coord = coord_1 + coord_3 + coord_4 + coord_5;

    bool is_attached = (total_coord >= 1.0);
    
    // --- HYSTERESIS LOGIC: Update contact timer ---
    // Increment timer if attached, reset if detached
    if (is_attached) {
        contact_timer[i]++;  // Bond is holding - count consecutive steps
    } else {
        contact_timer[i] = 0;  // Bond broke - reset timer
    }

    // Only proceed if we have "confidence" (Timer > Threshold)
    // ADJUSTED: Trigger when timer reaches threshold+1, but also allow triggering every 50 steps after that
    // This makes it less restrictive - molecules won't miss the exact trigger moment
    bool confident = (contact_timer[i] > hysteresis_threshold);
    // RANDOMIZED TRIGGER: Instead of deterministic intervals (every 50 steps), use probabilistic triggering
    // This desynchronizes molecules - they trigger independently based on probability, not fixed timesteps
    bool trigger_moment = (contact_timer[i] == (hysteresis_threshold + 1)) || 
                          (contact_timer[i] > hysteresis_threshold && random.uniform() < (1.0 / 50.0));  // 2% chance per check after threshold
    
    // Debug: track confident contacts
    if (confident && is_attached) {
        nconfident_contacts++;
    }

    bool should_change = false;
    int new_type = my_eff_type;

    // LOGIC A: Unreacted Monomer (Type 1) -> Changes when touching SAME-TYPE (Type 1) patches
    // Same rule as reacted monomers: HIGHER ID monomer must change
    if (my_eff_type == 1 && confident && trigger_moment) {
        bool touching_same_type_1 = (coord_1 >= 1.0);  // Type 1 touching other Type 1 patches
        
        if (touching_same_type_1) {
            // Check if I have higher ID than all conflicted neighbors
            int my_mol = molecule[i];
            bool all_neighbors_lower_id = true;
            double cutoffsq = cutoff * cutoff;
            
            // Find all type 1 neighbors I'm touching
            int nall = atom->nlocal + atom->nghost;
            for (int j = 0; j < nall; j++) {
                if (j == i) continue;
                if (molecule[j] == my_mol) continue;
                if (!(mask[j] & group->bitmask[group_patches])) continue;
                if (effective_type[j] != 1) continue;  // Must be type 1
                
                // Check distance
                double dx = atom->x[j][0] - atom->x[i][0];
                double dy = atom->x[j][1] - atom->x[i][1];
                double dz = atom->x[j][2] - atom->x[i][2];
                double *prd = domain->prd;
                if (prd[0] > 0.0) dx -= prd[0] * std::round(dx / prd[0]);
                if (prd[1] > 0.0) dy -= prd[1] * std::round(dy / prd[1]);
                if (prd[2] > 0.0) dz -= prd[2] * std::round(dz / prd[2]);
                double rsq = dx*dx + dy*dy + dz*dz;
                
                if (rsq < cutoffsq) {
                    int neighbor_mol = molecule[j];
                    if (neighbor_mol > my_mol) {
                        all_neighbors_lower_id = false;  // Neighbor has higher ID
                    }
                }
            }
            
            // If I have higher ID than all neighbors, I must change
            if (all_neighbors_lower_id) {
             should_change = true;
                // Pick new type randomly from 3, 4, 5 (33% each)
                // Note: Type 1 can evolve to any of {3, 4, 5} with equal probability
             double r = random.uniform();
             if (r < 0.333333) new_type = 3;
             else if (r < 0.666666) new_type = 4;
             else new_type = 5;
             nattempts++;
                // Reset timer with jitter
             int jitter_magnitude = 10;
             int random_jitter = static_cast<int>(random.uniform() * jitter_magnitude);
                contact_timer[i] = -random_jitter;
            }
        }
    }
    // LOGIC B: Reacted Monomer touches SAME type -> SYMMETRY-BREAKING LOGIC
    // This implements the "Seniority + Cooldown" mechanism to prevent synchronized state changes
    else if ((my_eff_type == 3 || my_eff_type == 4 || my_eff_type == 5) && confident && trigger_moment) {
        bool touching_same = false;
        if (my_eff_type == 3 && coord_3 >= 1.0) touching_same = true;
        if (my_eff_type == 4 && coord_4 >= 1.0) touching_same = true;
        if (my_eff_type == 5 && coord_5 >= 1.0) touching_same = true;

        if (touching_same) {
            // Check if self is in cooldown (using per-atom randomized duration)
            bool self_in_cooldown = false;
            if (last_change[i] >= 0) {
                bigint steps_since_change = timestep - last_change[i];
                if (steps_since_change < cooldown_duration[i]) {
                    self_in_cooldown = true;
                }
            }
            
            // If self is in cooldown, do nothing (I just changed, so I'm stable)
            if (self_in_cooldown) {
                // Do nothing - stay locked
                ncooldown_blocked++;
            } else {
                // Self is NOT in cooldown. Check neighbors to see if I should change.
                // CORRECTED LOGIC: When same-type patches touch, the monomer with HIGHER ID must change
                // The change is deterministic (if condition met, must attempt), but type selection is random
                int my_mol = molecule[i];
                bool should_i_change = false;
                bool found_conflict = false;
                bool all_neighbors_lower_id = true;  // True if my ID > all neighbor IDs (I'm higher)
                double cutoffsq = cutoff * cutoff;
                
                // PASS 1: Check all neighbors to find same-type conflicts
                int nall = atom->nlocal + atom->nghost;  // Check ghosts too for MPI
                for (int j = 0; j < nall; j++) {
                    if (j == i) continue;
                    if (molecule[j] == my_mol) continue;  // Skip same molecule
                    if (!(mask[j] & group->bitmask[group_patches])) continue;  // Must be a patch
                    if (type[j] < 2 || type[j] > 5) continue;
                    
                    // Check if j is of the same effective type (SAME TYPE CONFLICT)
                    int j_eff_type = effective_type[j];
                    if (j_eff_type != my_eff_type) continue;
                    
                    // Check distance to see if we're attached (touching/interacting)
                    double dx = atom->x[j][0] - atom->x[i][0];
                    double dy = atom->x[j][1] - atom->x[i][1];
                    double dz = atom->x[j][2] - atom->x[i][2];
                    double *prd = domain->prd;
                    if (prd[0] > 0.0) dx -= prd[0] * std::round(dx / prd[0]);
                    if (prd[1] > 0.0) dy -= prd[1] * std::round(dy / prd[1]);
                    if (prd[2] > 0.0) dz -= prd[2] * std::round(dz / prd[2]);
                    double rsq = dx*dx + dy*dy + dz*dz;
                    
                    if (rsq < cutoffsq) {
                        found_conflict = true;
                        int neighbor_mol = molecule[j];
                        
                        // Check if neighbor has HIGHER ID than me
                        if (neighbor_mol > my_mol) {
                            all_neighbors_lower_id = false;  // At least one neighbor has higher ID
                        }
                    }
                }
                
                // PASS 2: Apply rule - HIGHER ID monomer must change
                // If I have higher ID than ALL conflicted neighbors, I must change
                if (found_conflict && all_neighbors_lower_id) {
                    should_i_change = true;
                }
                
                // PASS 3: If condition is met, ALWAYS attempt change (no probability threshold)
                // The randomness comes from type selection (50% chance each for the OTHER two types)
                if (should_i_change && found_conflict) {
                    ntrigger_attempts++;
                        should_change = true;

                    // CRITICAL FIX: TYPE SELECTION must EXCLUDE current type to resolve conflict
                    // If we're type 3, pick from {4, 5} with 50% each
                    // If we're type 4, pick from {3, 5} with 50% each
                    // If we're type 5, pick from {3, 4} with 50% each
                            double r = random.uniform();
                    if (my_eff_type == 3) {
                        // Currently type 3, touching type 3 -> change to 4 or 5
                        new_type = (r < 0.5) ? 4 : 5;
                    } else if (my_eff_type == 4) {
                        // Currently type 4, touching type 4 -> change to 3 or 5
                        new_type = (r < 0.5) ? 3 : 5;
                    } else if (my_eff_type == 5) {
                        // Currently type 5, touching type 5 -> change to 3 or 4
                        new_type = (r < 0.5) ? 3 : 4;
                    } else {
                        // Fallback (should not happen) - pick from all 3
                        if (r < 0.333333) new_type = 3;
                        else if (r < 0.666666) new_type = 4;
                        else new_type = 5;
                        }

                        nattempts++;
                    // Reset timer with jitter to desynchronize future triggers
                        int jitter_magnitude = 10;
                        int random_jitter = static_cast<int>(random.uniform() * jitter_magnitude);
                        contact_timer[i] = -random_jitter;  // Negative jitter delays next trigger
                }
            }
        }
    }

    // Schedule the change for the WHOLE molecule (monochromatic rule)
    if (should_change) {
        mol_changes[molecule[i]] = new_type;
    }
    
    // Store history
    prev_coord[i] = total_coord;
  }

  // 2. GLOBAL RATE LIMITER: Prevent mass synchronized changes
  // CRITICAL FIX: Limit must be enforced GLOBALLY across all processors
  // Problem: Each processor was limiting to 2 changes locally, so with N processors = 2N total changes!
  // Solution: Check global total first, then proportionally reduce per processor
  
  int local_change_count = static_cast<int>(mol_changes.size());
  int total_changes_global = 0;
  MPI_Allreduce(&local_change_count, &total_changes_global, 1, MPI_INT, MPI_SUM, world);
  
  if (total_changes_global > max_changes_per_step) {
      // Create random generator for rate limiter (separate from loop to avoid scope issues)
      int rate_limiter_seed = seed + static_cast<int>(timestep) + comm->me + 7777;  // Different offset
      RanPark rate_limiter_random(lmp, rate_limiter_seed);
      
      // Too many want to change globally - proportionally reduce each processor's quota
      // This ensures total global changes <= max_changes_per_step
      double reduction_factor = static_cast<double>(max_changes_per_step) / static_cast<double>(total_changes_global);
      int allowed_local = static_cast<int>(local_change_count * reduction_factor + 0.5);
      if (allowed_local > local_change_count) allowed_local = local_change_count;
      if (allowed_local < 0) allowed_local = 0;
      
      if (allowed_local < local_change_count) {
          // Randomly select which molecules to keep
          std::vector<int> mol_ids;
          for (const auto& pair : mol_changes) {
              mol_ids.push_back(pair.first);
          }
          
          // Shuffle
          for (size_t i = 0; i < mol_ids.size(); i++) {
              size_t j = static_cast<size_t>(rate_limiter_random.uniform() * (mol_ids.size() - i)) + i;
              std::swap(mol_ids[i], mol_ids[j]);
          }
          
          // Keep only allowed_local
          std::map<int, int> limited_mol_changes;
          for (int i = 0; i < allowed_local && i < static_cast<int>(mol_ids.size()); i++) {
              limited_mol_changes[mol_ids[i]] = mol_changes[mol_ids[i]];
          }
          mol_changes = limited_mol_changes;
      }
  }

  // 3. UPDATE LOOP: Apply changes to ALL atoms in the affected molecules
  // This enforces "monochromatic" rule: all patches in a monomer get the same type
  // CRITICAL: Randomize cooldown PER MOLECULE (not per atom) to ensure all patches in molecule wake up together
  std::map<int, int> mol_cooldown_durations;  // Store randomized cooldown per molecule
  if (!mol_changes.empty()) {
      // Create random generator for cooldown assignment (seed includes timestep and processor)
      int update_seed = seed + static_cast<int>(timestep) + comm->me + 9999;  // Offset to avoid correlation
      RanPark update_random(lmp, update_seed);
      
      // First pass: Assign randomized cooldown to each molecule
      for (const auto& pair : mol_changes) {
          int mol_id = pair.first;
          // RANDOMIZED COOLDOWN: INCREASED VARIANCE - 0.3x to 2.0x of base cooldown_steps (wider spread breaks synchronization better)
          double cooldown_factor = 0.3 + update_random.uniform() * 1.7;  // Range: 0.3 to 2.0 (wider spread prevents mass behavior)
          mol_cooldown_durations[mol_id] = static_cast<int>(cooldown_steps * cooldown_factor);
      }
      
      // Second pass: Apply changes and cooldowns to all atoms in affected molecules
      for (int i = 0; i < nlocal; i++) {
          auto it = mol_changes.find(molecule[i]);
          if (it != mol_changes.end()) {
              effective_type[i] = it->second; // UNIFY - all patches in molecule get same type
              last_change[i] = timestep;
              // Apply same randomized cooldown to all atoms in this molecule
              cooldown_duration[i] = mol_cooldown_durations[molecule[i]];
              
              // CRITICAL FIX: Reset timer with jitter to desynchronize future triggers
              int jitter_magnitude = 10;
              int mol_id = molecule[i];
              int timer_seed = seed + static_cast<int>(timestep) + comm->me + mol_id + 8888;  // Different offset for jitter
              RanPark timer_random(lmp, timer_seed);
              int random_jitter = static_cast<int>(timer_random.uniform() * jitter_magnitude);
              contact_timer[i] = -random_jitter;  // Negative jitter delays next trigger
              nchanges++;
          }
      }
  }

  // Sync results across processors
  comm->forward_comm(this);
  
  // Debug output (every 10000 steps to avoid spam)
  if (timestep % 10000 == 0 && comm->me == 0) {
    int total_confident = 0, total_attempts = 0, total_blocked = 0, total_changes = 0;
    MPI_Allreduce(&nconfident_contacts, &total_confident, 1, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(&ntrigger_attempts, &total_attempts, 1, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(&ncooldown_blocked, &total_blocked, 1, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(&nchanges, &total_changes, 1, MPI_INT, MPI_SUM, world);
    if (update->ntimestep > 30000) {  // Only print after equilibration
      fprintf(stderr, "DEBUG[step %ld]: Confident=%d, TriggerAttempts=%d, CooldownBlocked=%d, Changes=%d\n",
              timestep, total_confident, total_attempts, total_blocked, total_changes);
    }
  }
  // Reset debug counters for next check
  nconfident_contacts = 0;
  ntrigger_attempts = 0;
  ncooldown_blocked = 0;

  // --- SYNC FIX: Consistency Sweep ---
  // This fixes the "Fruit Salad" bug where split monomers have different colors across MPI boundaries.
  // Logic: "If my brother (ghost or local) has a newer timestamp, copy him."
  int nall = atom->nlocal + atom->nghost; // Look at ghosts too
  std::map<int, int> best_type_per_mol;
  std::map<int, bigint> newest_time_per_mol;

  // Pass 1: Find the most recent change for each molecule (checking all atoms including ghosts)
  for (int i = 0; i < nall; i++) {
      if (molecule[i] <= 0) continue;
      
      // Only consider patches (types 2-5)
      if (type[i] < 2 || type[i] > 5) continue;
      
      int mol_id = molecule[i];
      
      // If this atom has no history, skip
      if (last_change[i] < 0) continue;

      // --- SANITY FILTER: Ignore "Ghost from the Future" ---
      // If a ghost atom claims to be from the future, it's garbage memory - ignore it.
      // This prevents impossible timestamps (e.g., 9,999,999,999) from overwriting real state changes.
      if (last_change[i] > timestep) continue;

      // Check if this atom has a newer change than what we've seen for this mol
      if (newest_time_per_mol.find(mol_id) == newest_time_per_mol.end() ||
          last_change[i] > newest_time_per_mol[mol_id]) {
          newest_time_per_mol[mol_id] = last_change[i];
          best_type_per_mol[mol_id] = effective_type[i];
      }
  }

  // Pass 2: Force local atoms to match the newest type found
  // CRITICAL FIX (Issues 1 & 8): Prevent consistency sweep from thrashing molecules in cooldown
  // This prevents the sweep from forcing state changes that bypass cooldown periods,
  // which was causing cascading conflicts and mass synchronized changes
  int consistency_fixes = 0;
  // Initialize sweep diagnostic counters
  nsweep_attempts = 0;
  nsweep_blocked_cooldown = 0;
  nsweep_blocked_timestamp = 0;
  nsweep_applied = 0;
  for (int i = 0; i < nlocal; i++) {
      if (molecule[i] <= 0) continue;
      if (type[i] < 2 || type[i] > 5) continue;
      
          int mol_id = molecule[i];
      if (best_type_per_mol.find(mol_id) != best_type_per_mol.end()) {
          nsweep_attempts++;  // Track every molecule we consider for consistency enforcement
          int target_type = best_type_per_mol[mol_id];
          
          // CRITICAL FIX: Multiple gates to prevent consistency sweep from causing cascades
          
          // GATE 1: Cooldown check - skip if molecule is currently in cooldown
          bool in_cooldown = false;
          if (last_change[i] >= 0) {
              bigint steps_since_change = timestep - last_change[i];
              // Ensure we don't overwrite a local change that is valid and cooling down
              // with an "older" or "conflicting" update from a confused ghost
              if (steps_since_change < cooldown_duration[i]) {
                  in_cooldown = true;
                  nsweep_blocked_cooldown++;  // Track cooldown blocks
              }
          }
          
          // GATE 2: Only update if target timestamp is SIGNIFICANTLY newer (not just slightly newer)
          // This prevents the sweep from forcing updates based on tiny timestamp differences
          // that might occur due to MPI communication timing or race conditions
          // Require at least 100 steps difference to avoid forcing updates from recent changes
          bool target_is_significantly_newer = false;
          bigint target_timestamp = newest_time_per_mol[mol_id];
          if (last_change[i] < 0) {
              // No local history - always accept (molecule never changed before)
              target_is_significantly_newer = true;
          } else {
              // Additional safety: ignore if target timestamp is in the future (shouldn't happen, but double-check)
              if (target_timestamp > timestep) {
                  target_is_significantly_newer = false;
                  nsweep_blocked_timestamp++;  // DIAGNOSTIC: Track future timestamp blocks
              } else {
                  // Require target to be at least 100 steps newer to avoid race conditions
                  bigint timestamp_diff = target_timestamp - last_change[i];
                  if (timestamp_diff > 100) {
                      target_is_significantly_newer = true;
                  } else {
                      nsweep_blocked_timestamp++;  // DIAGNOSTIC: Track timestamp blocks (difference too small)
                  }
              }
          }
          
          // GATE 3: Only apply consistency enforcement if all gates pass and types differ
          // CRITICAL: Set last_change to CURRENT timestep, not target_timestamp
          // This prevents the sweep from "backdating" changes which could reset cooldowns incorrectly
          if (!in_cooldown && target_is_significantly_newer && effective_type[i] != target_type) {
               effective_type[i] = target_type;
               last_change[i] = timestep;  // Use current timestep, not target_timestamp
               consistency_fixes++;
               nsweep_applied++;  // Track successful sweep updates
          }
      }
  }
  
  // Final Sync to ensure ghosts are updated with the corrections
  comm->forward_comm(this);
  MPI_Allreduce(MPI_IN_PLACE, &nchanges, 1, MPI_INT, MPI_SUM, world);
  
  // Debug output (every 10000 steps to avoid spam)
  if (timestep % 10000 == 0 && comm->me == 0) {
    int total_confident = 0, total_attempts = 0, total_blocked = 0, total_changes_glob = 0;
    int total_sweep_attempts = 0, total_sweep_blocked_cooldown = 0, total_sweep_blocked_timestamp = 0, total_sweep_applied = 0;
    MPI_Allreduce(&nconfident_contacts, &total_confident, 1, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(&ntrigger_attempts, &total_attempts, 1, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(&ncooldown_blocked, &total_blocked, 1, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(&nsweep_attempts, &total_sweep_attempts, 1, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(&nsweep_blocked_cooldown, &total_sweep_blocked_cooldown, 1, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(&nsweep_blocked_timestamp, &total_sweep_blocked_timestamp, 1, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(&nsweep_applied, &total_sweep_applied, 1, MPI_INT, MPI_SUM, world);
    total_changes_glob = nchanges;  // Already reduced above
    if (timestep > 30000) {  // Only print after equilibration
      fprintf(stderr, "DEBUG[step %ld]: Confident=%d, TriggerAttempts=%d, CooldownBlocked=%d, Changes=%d | SWEEP: Attempts=%d, Blocked(Cooldown=%d,Timestamp=%d), Applied=%d\n",
              timestep, total_confident, total_attempts, total_blocked, total_changes_glob,
              total_sweep_attempts, total_sweep_blocked_cooldown, total_sweep_blocked_timestamp, total_sweep_applied);
    }
  }
  // Reset debug counters for next check
  nconfident_contacts = 0;
  ntrigger_attempts = 0;
  ncooldown_blocked = 0;
  nsweep_attempts = 0;
  nsweep_blocked_cooldown = 0;
  nsweep_blocked_timestamp = 0;
  nsweep_applied = 0;
}

/* ---------------------------------------------------------------------- */

void FixStateChangeOctahedron::update_atom_types()
{
  // CRITICAL: Store original atom count before any changes
  bigint natoms_before = atom->natoms;
  
  // CRITICAL: First, synchronize effective_type across all processors
  // This ensures patches on different processors get the same effective_type
  comm->forward_comm(this);
  
  // Update actual atom types based on effective_type
  // IMPORTANT CHANGE:
  // - We no longer perform any "majority vote" or unification across a molecule.
  // - Each patch keeps the effective_type it acquired from local contacts
  //   in check_and_change(), and we simply map that to the LAMMPS atom type.
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  
  // Map per-patch effective_type to actual LAMMPS atom type, patch by patch
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if (!(mask[i] & group->bitmask[group_patches])) continue;
    if (type[i] != 2 && type[i] != 3 && type[i] != 4 && type[i] != 5) continue;

    int eff = effective_type[i];
    int new_lammps_type;

    // Effective type 1 corresponds to LAMMPS type 2 (initial patch state)
    if (eff == 1) {
      new_lammps_type = 2;
    // Effective types 3,4,5 map directly to LAMMPS types 3,4,5
    } else if (eff == 3 || eff == 4 || eff == 5) {
      new_lammps_type = eff;
    } else {
      // Unknown effective type: leave atom type unchanged
      continue;
    }

    if (type[i] != new_lammps_type) {
      type[i] = new_lammps_type;
    }
  }
  
  // CRITICAL: Verify atom count before communication
  if (atom->natoms != natoms_before) {
    error->all(FLERR, "Fix state/change/octahedron: Atom count changed during type update - aborting to prevent corruption");
  }
  
  // CRITICAL: Update ghost atoms (for MPI) - this is critical for parallel simulations
  // Use forward_comm() to update all atom properties including types
  comm->forward_comm();
  
  // CRITICAL: Verify atom count after communication
  if (atom->natoms != natoms_before) {
    error->all(FLERR, "Fix state/change/octahedron: Atom count changed after communication - possible atom loss");
  }
  
  // Force neighbor list rebuild on next step
  // This ensures pair interactions are recalculated with new types
  neighbor->decide();
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
  memory->grow(cooldown_duration, nmax, "fix_state_change_octahedron:cooldown_duration");
  memory->grow(effective_type, nmax, "fix_state_change_octahedron:effective_type");
  memory->grow(prev_coord, nmax, "fix_state_change_octahedron:prev_coord");
  memory->grow(contact_timer, nmax, "fix_state_change_octahedron:contact_timer");
  
  // Initialize new memory to safe defaults
  for (int i = atom->nmax; i < nmax; i++) {
    last_change[i] = -1;
    cooldown_duration[i] = cooldown_steps;  // Default to base cooldown, will be randomized on first change
    effective_type[i] = 1;  // Default initial patch type
    prev_coord[i] = 0.0;
    contact_timer[i] = 0;
  }
}

/* ---------------------------------------------------------------------- */

void FixStateChangeOctahedron::copy_arrays(int i, int j, int /*delflag*/)
{
  last_change[j] = last_change[i];
  cooldown_duration[j] = cooldown_duration[i];
  effective_type[j] = effective_type[i];
  prev_coord[j] = prev_coord[i];
  contact_timer[j] = contact_timer[i];
}

/* ---------------------------------------------------------------------- */

int FixStateChangeOctahedron::pack_exchange(int i, double *buf)
{
  int n = 0;
  buf[n++] = static_cast<double>(last_change[i]);
  buf[n++] = static_cast<double>(cooldown_duration[i]);
  buf[n++] = static_cast<double>(effective_type[i]);
  buf[n++] = prev_coord[i];
  return n;
}

/* ---------------------------------------------------------------------- */

int FixStateChangeOctahedron::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  last_change[nlocal] = static_cast<int>(buf[n++]);
  cooldown_duration[nlocal] = static_cast<int>(buf[n++]);
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

/* ---------------------------------------------------------------------- */

int FixStateChangeOctahedron::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i, j, m = 0;
  
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = static_cast<double>(last_change[j]);
    buf[m++] = static_cast<double>(cooldown_duration[j]);
    buf[m++] = static_cast<double>(effective_type[j]);
    buf[m++] = prev_coord[j];
  }
  
  return m;
}

/* ---------------------------------------------------------------------- */

void FixStateChangeOctahedron::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m = 0, last = first + n;
  
  for (i = first; i < last; i++) {
    last_change[i] = static_cast<int>(buf[m++]);
    cooldown_duration[i] = static_cast<int>(buf[m++]);
    effective_type[i] = static_cast<int>(buf[m++]);
    prev_coord[i] = buf[m++];
  }
}

