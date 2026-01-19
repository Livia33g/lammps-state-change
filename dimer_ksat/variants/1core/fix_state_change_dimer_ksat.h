/* -*- c++ -*- ----------------------------------------------------------*/
/* State-change fix for "dimer_ksat" A->C catalyzed by B contact.
 *
 * Goal behavior:
 * - Only A can switch to C.
 * - Switching is triggered by A-patch (type 2) being within `cutoff` of any
 *   B-patch (type 3) for `hysteresis_checks` consecutive checks.
 * - When a molecule switches, all its patch atoms change type 2 -> 4.
 * - C-patches (type 4) never switch.
 *
 * Arguments:
 *   fix ID group-ID state/change/dimer_ksat nevery cutoff pflip group_patches [hysteresis_checks]
 *
 * Notes:
 * - This is intentionally minimal (like state/change/dimer) and is meant to be
 *   paired with a LAMMPS input where pair_coeff(2,3) == pair_coeff(3,4).
 * ------------------------------------------------------------------------- */
#ifdef FIX_CLASS
// clang-format off
FixStyle(state/change/dimer_ksat,FixStateChangeDimerKsat);
// clang-format on
#else

#ifndef LMP_FIX_STATE_CHANGE_DIMER_KSAT_H
#define LMP_FIX_STATE_CHANGE_DIMER_KSAT_H

#include "fix.h"
#include <vector>

namespace LAMMPS_NS {

class FixStateChangeDimerKsat : public Fix {
 public:
  FixStateChangeDimerKsat(class LAMMPS *, int, char **);
  int setmask() override;
  void post_integrate() override;

 private:
  int nevery;
  double cutoff, cutoffsq;
  double pflip;
  int seed;
  int group_patches;
  int hysteresis_checks;

  int *contact_counter;

  class RanPark *random;

  void detect_and_schedule(int i, std::vector<int> &mol_list);
};

}  // namespace LAMMPS_NS

#endif

#endif





