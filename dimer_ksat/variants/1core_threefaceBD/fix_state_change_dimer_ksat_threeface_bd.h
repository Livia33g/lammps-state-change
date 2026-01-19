/* -*- c++ -*- ----------------------------------------------------------*/
/* State-change fix for "dimer_ksat_threeface_bd".
 *
 * Geometry / types (expected by the generator):
 * - A patches: type 2 (switchable)
 * - C patches: type 4 (sink, non-switchable)
 *
 * - B monomer: core type CORE_B (default 7), patches on 3 faces:
 *     face0: type 3
 *     face1: type 5
 *     face2: type 6
 * - D monomer: core type CORE_D (default 8), patches on the same 3 faces (types 3/5/6)
 *
 * Rule:
 * - For B molecules: if >= 2 distinct A molecules are attached anywhere (any face),
 *   then after hysteresis, flip the lowest-ID attached A molecule 2->4.
 * - For D molecules: if ALL 3 faces have at least one A attached (face0/1/2),
 *   then after hysteresis, flip the lowest-ID attached A molecule 2->4.
 *
 * Arguments:
 *   fix ID group-ID state/change/dimer_ksat_threeface_bd nevery cutoff pflip group_patches [hysteresis_checks]
 *
 * Notes:
 * - Patches of types 3/5/6 should have identical energetics to A/C (pair_coeff equality).
 * ------------------------------------------------------------------------- */
#ifdef FIX_CLASS
// clang-format off
FixStyle(state/change/dimer_ksat_threeface_bd,FixStateChangeDimerKsatThreeFaceBD);
// clang-format on
#else

#ifndef LMP_FIX_STATE_CHANGE_DIMER_KSAT_THREEFACE_BD_H
#define LMP_FIX_STATE_CHANGE_DIMER_KSAT_THREEFACE_BD_H

#include "fix.h"

#include <unordered_map>

namespace LAMMPS_NS {

class FixStateChangeDimerKsatThreeFaceBD : public Fix {
 public:
  FixStateChangeDimerKsatThreeFaceBD(class LAMMPS *, int, char **);
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

  // per-(B or D)-molecule consecutive-check counter for its rule condition
  std::unordered_map<int, int> mol_counter;
};

}  // namespace LAMMPS_NS

#endif

#endif


