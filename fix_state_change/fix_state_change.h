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

#ifdef FIX_CLASS
// clang-format off
FixStyle(state/change,FixStateChange);
// clang-format on
#else

#ifndef LMP_FIX_STATE_CHANGE_H
#define LMP_FIX_STATE_CHANGE_H

#include "fix.h"

namespace LAMMPS_NS {

class RanPark;

class FixStateChange : public Fix {
 public:
  FixStateChange(class LAMMPS *, int, char **);
  ~FixStateChange() override;
  int setmask() override;
  void init() override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void end_of_step() override;
  double compute_scalar() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;

 private:
  int nlevels_respa;
  int check_every;              // Check state changes every N steps
  int cooldown_steps;           // Cooldown period after state change
  double probability;           // Probability of state change when conditions met
  double cutoff;                // Cutoff distance for coordination check
  int group_A;                  // Group ID for type 2 patches
  int group_B;                  // Group ID for type 3 patches
  
  int *last_change;             // Per-atom: -1 if ready for shot, timestep if had shot while attached
  int *effective_type;          // Per-atom: effective type (2 or 3) for patches
  double *prev_coord_A;         // Per-atom: previous coordination with type 2 (to detect attachment)
  double *prev_coord_B;         // Per-atom: previous coordination with type 3 (to detect attachment)
  
  class Compute *c_coord_A;      // Compute for coordination with type 2
  class Compute *c_coord_B;      // Compute for coordination with type 3
  
  bigint next_check;            // Next timestep to check for state changes
  int nchanges;                 // Number of state changes this check (local)
  int nattempts;                // Number of change attempts this check (local)
  int seed;                     // Random seed
  
  void check_and_change();
  void update_atom_types();
  int get_coordination(int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif

