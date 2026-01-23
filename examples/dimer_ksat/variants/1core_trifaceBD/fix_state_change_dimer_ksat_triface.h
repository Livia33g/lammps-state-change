/* -*- c++ -*- ----------------------------------------------------------*/
/* State-change fix for "dimer_ksat_triface":
 *
 * System conventions (by atom type):
 * - A patches: type 2
 * - C patches: type 4
 * - B/D patches live on 3 faces, encoded as 3 patch types:
 *     face +x : type 3
 *     face mid: type 6   (between +x and -x face)
 *     face -x : type 5
 * - B vs D is distinguished by core type on the molecule:
 *     B core: type 1
 *     D core: type 7
 *
 * Rules:
 * - For B molecules (core type 1):
 *     If B has >=2 DISTINCT A molecules attached across ANY faces (3/6/5),
 *     then after hysteresis, flip the lowest-ID A among the attached set: 2->4.
 * - For D molecules (core type 7):
 *     If D has ALL 3 faces occupied by A simultaneously (faces 3,6,5 each have an A),
 *     and those correspond to 3 DISTINCT A molecules, then after hysteresis flip
 *     the lowest-ID A among the three: 2->4.
 *
 * Arguments:
 *   fix ID group-ID state/change/dimer_ksat_triface nevery cutoff pflip group_patches [hysteresis_checks]
 * ------------------------------------------------------------------------- */
#ifdef FIX_CLASS
// clang-format off
FixStyle(state/change/dimer_ksat_triface,FixStateChangeDimerKsatTriFace);
// clang-format on
#else

#ifndef LMP_FIX_STATE_CHANGE_DIMER_KSAT_TRIFACE_H
#define LMP_FIX_STATE_CHANGE_DIMER_KSAT_TRIFACE_H

#include "fix.h"

#include <unordered_map>

namespace LAMMPS_NS {

class FixStateChangeDimerKsatTriFace : public Fix {
 public:
  FixStateChangeDimerKsatTriFace(class LAMMPS *, int, char **);
  int setmask() override;
  void post_integrate() override;

 private:
  int nevery;
  double cutoff, cutoffsq;
  double pflip;
  int seed;
  int group_patches;
  int hysteresis_checks;

  class RanPark *random;

  // per-(B or D) molecule consecutive-check counter for trigger condition
  std::unordered_map<int, int> mol_counter;
};

}  // namespace LAMMPS_NS

#endif

#endif



