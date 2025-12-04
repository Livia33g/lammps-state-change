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
FixStyle(state/change/octahedron,FixStateChangeOctahedron);
// clang-format on
#else

#ifndef LMP_FIX_STATE_CHANGE_OCTAHEDRON_H
#define LMP_FIX_STATE_CHANGE_OCTAHEDRON_H

#include "fix.h"

namespace LAMMPS_NS {

class RanPark;

class FixStateChangeOctahedron : public Fix {
 public:
  FixStateChangeOctahedron(class LAMMPS *, int, char **);
  ~FixStateChangeOctahedron() override;
  int setmask() override;
  void init() override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void end_of_step() override;
  double compute_scalar() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;

 private:
  int check_every;              // Check state changes every N steps
  int cooldown_steps;           // Cooldown period after state change
  double probability;           // Probability of state change when conditions met
  double cutoff;                // Cutoff distance for coordination check
  int group_patches;            // Group ID for all patches (types 1, 3, 4, 5)
  
  int *last_change;             // Per-atom: timestep of last change (for cooldown)
  int *effective_type;          // Per-atom: effective type (1, 3, 4, or 5 for patches)
  double *prev_coord;           // Per-atom: previous coordination number
  
  bigint next_check;            // Next timestep to check for state changes
  int nchanges;                 // Number of state changes this check (local)
  int nattempts;                // Number of change attempts this check (local)
  int seed;                     // Random seed
  
  void check_and_change();
  void update_atom_types();
  double get_coordination(int, int);
  bool check_same_type_coordination(int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif

