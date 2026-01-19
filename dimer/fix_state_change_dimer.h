/* -*- c++ -*- ----------------------------------------------------------*/
/* A minimal state-change fix for dimer system:                                  */
/* - Detects any contact between Red patches (type 2) and Blue patches (type 3)  */
/*   within a cutoff every `nevery` steps.                                       */
/* - With probability `pflip`, flips the entire Red molecule's patches 2 -> 3.   */
/* - Core atoms (type 1) remain unchanged.                                       */
/*                                                                              */
/* This is intentionally simpler than the octahedron fix:                       */
/* no cooldowns, no consistency sweep—just a deterministic molecule flip.       */
/* ------------------------------------------------------------------------- */
#ifdef FIX_CLASS
// clang-format off
FixStyle(state/change/dimer,FixStateChangeDimer);
// clang-format on
#else

#ifndef LMP_FIX_STATE_CHANGE_DIMER_H
#define LMP_FIX_STATE_CHANGE_DIMER_H

#include "fix.h"
#include <vector>

namespace LAMMPS_NS {

class FixStateChangeDimer : public Fix {
 public:
  FixStateChangeDimer(class LAMMPS *, int, char **);
  int setmask() override;
  void init() override;
  void post_integrate() override;

 private:
  int nevery;
  double cutoff, cutoffsq;
  double pflip;
  int seed;
  int group_patches;  // group id for patches
  int hysteresis_checks;  // require contact for N consecutive checks before flip

  int *contact_counter;   // per-atom contact counter (for patch atoms)

  class RanPark *random;
  class Pair *pair;  // pointer to pair potential for energy calculations

  void detect_and_schedule(int, std::vector<int> &);
  double calculate_flip_energy(int mol_id, int old_type, int new_type);
};

}  // namespace LAMMPS_NS

#endif

#endif

