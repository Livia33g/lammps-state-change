/* -*- c++ -*- ----------------------------------------------------------*/
/* Minimal state-change fix for "dimer_ksat_twoside" system:
 *
 * - Every `nevery` steps:
 *   For each B molecule, detect if it has an A bound on face-1 (B patches type 3)
 *   AND an A bound on face-2 (B patches type 5) simultaneously (within cutoff).
 *   If so for `hysteresis_checks` consecutive checks, flip the A molecule with the
 *   LOWEST molecule-ID among the two attached As: type-2 patches -> type-4 patches.
 *
 * Arguments:
 *   fix ID group-ID state/change/dimer_ksat_twoside nevery cutoff pflip group_patches [hysteresis_checks]
 *
 * Example:
 *   fix sc patches state/change/dimer_ksat_twoside 100 0.7 1.0 patches 5
 * ---------------------------------------------------------------------- */

#include "fix_state_change_dimer_ksat_twoside.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "random_park.h"
#include "update.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>
#include <unordered_map>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;

namespace {
constexpr int TYPE_A_PATCH = 2;
constexpr int TYPE_B_FACE1 = 3;
constexpr int TYPE_C_PATCH = 4;
constexpr int TYPE_B_FACE2 = 5;

struct MinPair {
  int minA_face1 = 0;
  int minA_face2 = 0;
};
}  // namespace

FixStateChangeDimerKsatTwoSide::FixStateChangeDimerKsatTwoSide(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp, narg, arg), random(nullptr) {
  if (narg < 7)
    error->all(FLERR,
               "Illegal fix state/change/dimer_ksat_twoside command "
               "(expected: ID group-ID state/change/dimer_ksat_twoside nevery cutoff pflip "
               "group_patches)");

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  cutoff = utils::numeric(FLERR, arg[4], false, lmp);
  pflip = utils::numeric(FLERR, arg[5], false, lmp);

  group_patches = group->find(arg[6]);
  if (group_patches < 0)
    error->all(FLERR, "Fix state/change/dimer_ksat_twoside group_patches not found");

  hysteresis_checks = 1;
  if (narg >= 8) {
    hysteresis_checks = utils::inumeric(FLERR, arg[7], false, lmp);
    if (hysteresis_checks < 1)
      error->all(FLERR, "Illegal hysteresis_checks in state/change/dimer_ksat_twoside");
  }

  if (nevery <= 0) error->all(FLERR, "Illegal nevery in state/change/dimer_ksat_twoside");
  if (cutoff <= 0.0) error->all(FLERR, "Illegal cutoff in state/change/dimer_ksat_twoside");
  if (pflip < 0.0 || pflip > 1.0)
    error->all(FLERR, "Illegal pflip in state/change/dimer_ksat_twoside");

  cutoffsq = cutoff * cutoff;
  seed = 12345;
  peratom_flag = 0;
  restart_peratom = 0;
  random = new RanPark(lmp, seed);
}

int FixStateChangeDimerKsatTwoSide::setmask() {
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

void FixStateChangeDimerKsatTwoSide::post_integrate() {
  if (update->ntimestep % nevery) return;

  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  double **x = atom->x;
  double *prd = domain->prd;

  // Map B-molecule -> (min A on face1, min A on face2)
  std::unordered_map<int, MinPair> mins;
  mins.reserve(64);

  auto in_patch_group = [&](int i) -> bool {
    return (mask[i] & group->bitmask[group_patches]);
  };

  auto min_image_rsq = [&](int i, int j) -> double {
    double dx = x[j][0] - x[i][0];
    double dy = x[j][1] - x[i][1];
    double dz = x[j][2] - x[i][2];
    if (prd) {
      if (prd[0] > 0.0) dx -= prd[0] * std::round(dx / prd[0]);
      if (prd[1] > 0.0) dy -= prd[1] * std::round(dy / prd[1]);
      if (prd[2] > 0.0) dz -= prd[2] * std::round(dz / prd[2]);
    }
    return dx * dx + dy * dy + dz * dz;
  };

  // 1) Collect lowest-ID A molecule on each B face, per B molecule.
  for (int i = 0; i < nlocal; i++) {
    if (!in_patch_group(i)) continue;
    if (type[i] != TYPE_B_FACE1 && type[i] != TYPE_B_FACE2) continue;

    const int molB = molecule[i];
    if (molB <= 0) continue;

    MinPair &mp = mins[molB];

    for (int j = 0; j < nall; j++) {
      if (j == i) continue;
      if (!in_patch_group(j)) continue;
      if (type[j] != TYPE_A_PATCH) continue;

      if (min_image_rsq(i, j) < cutoffsq) {
        const int molA = molecule[j];
        if (molA <= 0) continue;

        if (type[i] == TYPE_B_FACE1) {
          mp.minA_face1 = (mp.minA_face1 == 0) ? molA : std::min(mp.minA_face1, molA);
        } else {
          mp.minA_face2 = (mp.minA_face2 == 0) ? molA : std::min(mp.minA_face2, molA);
        }
      }
    }
  }

  if (mins.empty()) return;

  // 2) Apply hysteresis per B molecule and schedule flips of A molecules.
  std::set<int> to_flip_A;
  for (const auto &kv : mins) {
    const int molB = kv.first;
    const MinPair &mp = kv.second;

    const bool has_both = (mp.minA_face1 > 0 && mp.minA_face2 > 0 && mp.minA_face1 != mp.minA_face2);
    int &count = b_counter[molB];
    if (has_both) count += 1;
    else count = 0;

    if (has_both && count >= hysteresis_checks) {
      const int molA_target = std::min(mp.minA_face1, mp.minA_face2);
      if (molA_target > 0) to_flip_A.insert(molA_target);
      count = 0;
    }
  }

  if (to_flip_A.empty()) return;

  // 3) Apply flips (2->4) for scheduled A molecules (local atoms only).
  for (int molA : to_flip_A) {
    double r = random->uniform();
    if (r > pflip) continue;

    for (int i = 0; i < nlocal; i++) {
      if (molecule[i] != molA) continue;
      if (!in_patch_group(i)) continue;
      if (type[i] == TYPE_A_PATCH) type[i] = TYPE_C_PATCH;
    }

    if (comm->me == 0) {
      fprintf(stderr,
              "STATECHANGE dimer_ksat_twoside: step %ld molA %d flipped 2->4\n",
              update->ntimestep, molA);
    }
  }
}

/* ---------------------------------------------------------------------- */



