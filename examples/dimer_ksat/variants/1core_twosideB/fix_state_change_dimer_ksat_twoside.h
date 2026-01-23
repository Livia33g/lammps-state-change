/* -*- c++ -*- ----------------------------------------------------------*/
/* State-change fix for "dimer_ksat_twoside":
 *
 * Goal behavior:
 * - A patches (type 2) can flip to C patches (type 4).
 * - B has TWO patch faces, encoded as two patch types on the same molecule:
 *     - B face 1 patches: type 3
 *     - B face 2 patches: type 5
 * - A flips ONLY if a B molecule simultaneously has:
 *     - at least one A patch (type 2) within cutoff of any type-3 patch AND
 *     - at least one (other) A patch (type 2) within cutoff of any type-5 patch
 * - When triggered, flip the A molecule with the LOWEST molecule-ID among the two
 *   As attached on the two faces.
 * - C patches (type 4) never switch back.
 *
 * Arguments:
 *   fix ID group-ID state/change/dimer_ksat_twoside nevery cutoff pflip group_patches [hysteresis_checks]
 *
 * Notes:
 * - Intended to be paired with identical interactions between A and both B faces:
 *     pair_coeff(2,3) == pair_coeff(2,5) and pair_coeff(3,4) == pair_coeff(4,5)
 * ------------------------------------------------------------------------- */
#ifdef FIX_CLASS
// clang-format off
FixStyle(state/change/dimer_ksat_twoside,FixStateChangeDimerKsatTwoSide);
// clang-format on
#else

#ifndef LMP_FIX_STATE_CHANGE_DIMER_KSAT_TWOSIDE_H
#define LMP_FIX_STATE_CHANGE_DIMER_KSAT_TWOSIDE_H

#include "fix.h"

#include <unordered_map>
#include <vector>

namespace LAMMPS_NS {

class FixStateChangeDimerKsatTwoSide : public Fix {
 public:
  FixStateChangeDimerKsatTwoSide(class LAMMPS *, int, char **);
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

  // per-B-molecule consecutive-check counter for "A on both sides" condition
  std::unordered_map<int, int> b_counter;
};

}  // namespace LAMMPS_NS

#endif

#endif

/* -*- c++ -*- ----------------------------------------------------------*/
/* State-change fix for a "two-sided B" variant of dimer_ksat.
 *
 * Behavior:
 * - A patches are type 2 (switchable), B patches are type 3, C patches are type 4 (non-switchable).
 * - A can switch to C (2 -> 4) ONLY if a B molecule has at least one A bound on EACH side
 *   (right-side patch set and left-side patch set) at the same check time.
 * - When triggered, flip the A molecule with the LOWER molecule-ID among the two side-attached A's.
 * - C never switches back.
 *
 * Side detection for B:
 * - For a B patch atom, compute dx = x_patch - x_core for the same B molecule.
 * - dx >= 0 => "right" side, dx < 0 => "left" side.
 *
 * Arguments:
 *   fix ID group-ID state/change/dimer_ksat_twoside nevery cutoff pflip group_patches [hysteresis_checks]
 *
 * Notes:
 * - This is intentionally minimal (like state/change/dimer_ksat).
 * - Hysteresis is applied per B molecule (via its core atom index).
 * ------------------------------------------------------------------------- */
#ifdef FIX_CLASS
// clang-format off
FixStyle(state/change/dimer_ksat_twoside,FixStateChangeDimerKsatTwoSide);
// clang-format on
#else

#ifndef LMP_FIX_STATE_CHANGE_DIMER_KSAT_TWOSIDE_H
#define LMP_FIX_STATE_CHANGE_DIMER_KSAT_TWOSIDE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixStateChangeDimerKsatTwoSide : public Fix {
 public:
  FixStateChangeDimerKsatTwoSide(class LAMMPS *, int, char **);
  int setmask() override;
  void post_integrate() override;

 private:
  int nevery;
  double cutoff, cutoffsq;
  double pflip;
  int seed;
  int group_patches;
  int hysteresis_checks;

  int *contact_counter;   // per-atom counter used on B core atom indices

  class RanPark *random;
};

}  // namespace LAMMPS_NS

#endif

#endif


