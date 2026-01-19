/* -*- c++ -*- ----------------------------------------------------------*/
/* State-change fix for "dimer_ksat_twoside_twins".
 *
 * Two independent "twin" channels:
 *
 * Channel 1 (A/B/C):
 * - A patches: type 2
 * - C patches: type 4
 * - B has two patch faces:
 *     face 1 patches: type 3
 *     face 2 patches: type 5
 * - Rule: if a B molecule has A bound on BOTH faces simultaneously,
 *         flip the lower molecule-ID A (among those two) from 2->4.
 *
 * Channel 2 (E/D/F):
 * - E patches: type 8
 * - F patches: type 10
 * - D has two patch faces:
 *     face 1 patches: type 9
 *     face 2 patches: type 11
 * - Rule: if a D molecule has E bound on BOTH faces simultaneously,
 *         flip the lower molecule-ID E (among those two) from 8->10.
 *
 * Both rules use the same nevery/cutoff/pflip/hysteresis parameters.
 *
 * Arguments:
 *   fix ID group-ID state/change/dimer_ksat_twoside_twins nevery cutoff pflip group_patches [hysteresis_checks]
 * ------------------------------------------------------------------------- */
#ifdef FIX_CLASS
// clang-format off
FixStyle(state/change/dimer_ksat_twoside_twins,FixStateChangeDimerKsatTwoSideTwins);
// clang-format on
#else

#ifndef LMP_FIX_STATE_CHANGE_DIMER_KSAT_TWOSIDE_TWINS_H
#define LMP_FIX_STATE_CHANGE_DIMER_KSAT_TWOSIDE_TWINS_H

#include "fix.h"

#include <unordered_map>

namespace LAMMPS_NS {

class FixStateChangeDimerKsatTwoSideTwins : public Fix {
 public:
  FixStateChangeDimerKsatTwoSideTwins(class LAMMPS *, int, char **);
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

  // per-molecule consecutive-check counter for B and D trigger conditions
  std::unordered_map<int, int> mol_counter;
};

}  // namespace LAMMPS_NS

#endif

#endif



