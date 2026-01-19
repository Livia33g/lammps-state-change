/* -*- c++ -*- ----------------------------------------------------------
   Auto-generated LAMMPS fix header from policy YAML
   Policy: avoid_same_type_octahedron
   Problem: octahedron
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
 FixStyle(state/change/octahedron/avoid_same_type_octahedron,FixStateChangeOctahedronAvoidSameTypeOctahedron);
// clang-format on
#else

#ifndef LMP_FIX_STATE_CHANGE_OCTAHEDRON_AVOID_SAME_TYPE_OCTAHEDRON_H
#define LMP_FIX_STATE_CHANGE_OCTAHEDRON_AVOID_SAME_TYPE_OCTAHEDRON_H

#include "fix.h"

namespace LAMMPS_NS {

class RanPark;

class FixStateChangeOctahedronAvoidSameTypeOctahedron : public Fix {
 public:
  FixStateChangeOctahedronAvoidSameTypeOctahedron(class LAMMPS *, int, char **);
  ~FixStateChangeOctahedronAvoidSameTypeOctahedron() override;
  int setmask() override;
  void init() override;
  void post_force(int) override;
  void end_of_step() override;
  double compute_scalar() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;

 private:
  int check_every;
  int cooldown_steps;
  double probability;
  double cutoff;
  int group_patches;
  int hysteresis_threshold;
  int max_changes_per_step;

  int *last_change;
  int *cooldown_duration;
  int *effective_type;
  double *prev_coord;
  int *contact_timer;

  bigint next_check;
  int nchanges;
  int nattempts;
  int seed;

  int nconfident_contacts;
  int ntrigger_attempts;
  int ncooldown_blocked;

  void check_and_change();
  void update_atom_types();
  double get_coordination(int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
