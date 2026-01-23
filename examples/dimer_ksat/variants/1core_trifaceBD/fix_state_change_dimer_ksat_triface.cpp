/* -*- c++ -*- ----------------------------------------------------------*/
/* Minimal state-change fix for "dimer_ksat_triface".
 * See header for full description.
 * ---------------------------------------------------------------------- */

#include "fix_state_change_dimer_ksat_triface.h"

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
constexpr int TYPE_FACE_POS = 3;
constexpr int TYPE_C_PATCH = 4;
constexpr int TYPE_FACE_NEG = 5;
constexpr int TYPE_FACE_MID = 6;
constexpr int TYPE_CORE_B = 1;
constexpr int TYPE_CORE_D = 7;

struct FaceMins {
  int minA_pos = 0;
  int minA_mid = 0;
  int minA_neg = 0;
  int core_type = 0;  // 1 for B, 7 for D
};
}  // namespace

FixStateChangeDimerKsatTriFace::FixStateChangeDimerKsatTriFace(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp, narg, arg), random(nullptr) {
  if (narg < 7)
    error->all(FLERR,
               "Illegal fix state/change/dimer_ksat_triface command "
               "(expected: ID group-ID state/change/dimer_ksat_triface nevery cutoff pflip "
               "group_patches)");

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  cutoff = utils::numeric(FLERR, arg[4], false, lmp);
  pflip = utils::numeric(FLERR, arg[5], false, lmp);

  group_patches = group->find(arg[6]);
  if (group_patches < 0)
    error->all(FLERR, "Fix state/change/dimer_ksat_triface group_patches not found");

  hysteresis_checks = 1;
  if (narg >= 8) {
    hysteresis_checks = utils::inumeric(FLERR, arg[7], false, lmp);
    if (hysteresis_checks < 1)
      error->all(FLERR, "Illegal hysteresis_checks in state/change/dimer_ksat_triface");
  }

  if (nevery <= 0) error->all(FLERR, "Illegal nevery in state/change/dimer_ksat_triface");
  if (cutoff <= 0.0) error->all(FLERR, "Illegal cutoff in state/change/dimer_ksat_triface");
  if (pflip < 0.0 || pflip > 1.0)
    error->all(FLERR, "Illegal pflip in state/change/dimer_ksat_triface");

  cutoffsq = cutoff * cutoff;
  seed = 12345;
  peratom_flag = 0;
  restart_peratom = 0;
  random = new RanPark(lmp, seed);
}

int FixStateChangeDimerKsatTriFace::setmask() {
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

void FixStateChangeDimerKsatTriFace::post_integrate() {
  if (update->ntimestep % nevery) return;

  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  double **x = atom->x;
  double *prd = domain->prd;

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

  // Collect per-molecule face minima for all molecules that have face patches.
  std::unordered_map<int, FaceMins> mins;
  mins.reserve(128);

  // 1) Determine core type per molecule (local only)
  for (int i = 0; i < nlocal; i++) {
    int mol = molecule[i];
    if (mol <= 0) continue;
    if (type[i] == TYPE_CORE_B || type[i] == TYPE_CORE_D) {
      mins[mol].core_type = type[i];
    }
  }

  // 2) For each face patch (types 3,6,5), find lowest-ID A molecule within cutoff.
  for (int i = 0; i < nlocal; i++) {
    if (!in_patch_group(i)) continue;

    const int t = type[i];
    if (t != TYPE_FACE_POS && t != TYPE_FACE_MID && t != TYPE_FACE_NEG) continue;

    const int molX = molecule[i];  // molecule owning this face patch (B or D)
    if (molX <= 0) continue;

    FaceMins &fm = mins[molX];
    // If core type unknown (not found locally), skip; avoids misclassifying A/C molecules.
    if (fm.core_type != TYPE_CORE_B && fm.core_type != TYPE_CORE_D) continue;

    for (int j = 0; j < nall; j++) {
      if (j == i) continue;
      if (!in_patch_group(j)) continue;
      if (type[j] != TYPE_A_PATCH) continue;

      if (min_image_rsq(i, j) < cutoffsq) {
        const int molA = molecule[j];
        if (molA <= 0) continue;
        if (t == TYPE_FACE_POS) {
          fm.minA_pos = (fm.minA_pos == 0) ? molA : std::min(fm.minA_pos, molA);
        } else if (t == TYPE_FACE_MID) {
          fm.minA_mid = (fm.minA_mid == 0) ? molA : std::min(fm.minA_mid, molA);
        } else {
          fm.minA_neg = (fm.minA_neg == 0) ? molA : std::min(fm.minA_neg, molA);
        }
      }
    }
  }

  if (mins.empty()) return;

  // 3) Evaluate triggers (with hysteresis) and schedule flips.
  std::set<int> to_flip_A;
  std::set<int> triggered_by_B;
  std::set<int> triggered_by_D;

  for (const auto &kv : mins) {
    const int molX = kv.first;
    const FaceMins &fm = kv.second;

    if (fm.core_type != TYPE_CORE_B && fm.core_type != TYPE_CORE_D) continue;

    // Unique A molecule IDs attached to this molecule (across faces)
    std::vector<int> attached;
    attached.reserve(3);
    if (fm.minA_pos > 0) attached.push_back(fm.minA_pos);
    if (fm.minA_mid > 0) attached.push_back(fm.minA_mid);
    if (fm.minA_neg > 0) attached.push_back(fm.minA_neg);
    std::sort(attached.begin(), attached.end());
    attached.erase(std::unique(attached.begin(), attached.end()), attached.end());

    bool trigger = false;
    if (fm.core_type == TYPE_CORE_B) {
      // B: trigger if >=2 distinct A attached
      trigger = (static_cast<int>(attached.size()) >= 2);
    } else {
      // D: trigger only if all 3 faces have A and they are 3 distinct molecules
      trigger = (fm.minA_pos > 0 && fm.minA_mid > 0 && fm.minA_neg > 0 &&
                 static_cast<int>(attached.size()) >= 3);
    }

    int &count = mol_counter[molX];
    if (trigger) count += 1;
    else count = 0;

    if (trigger && count >= hysteresis_checks) {
      const int molA_target = attached.empty() ? 0 : attached.front();
      if (molA_target > 0) {
        to_flip_A.insert(molA_target);
        if (fm.core_type == TYPE_CORE_B) triggered_by_B.insert(molX);
        else triggered_by_D.insert(molX);
      }
      count = 0;
    }
  }

  if (to_flip_A.empty()) return;

  // 4) Apply flips (2->4) for scheduled A molecules (local atoms only).
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
              "STATECHANGE dimer_ksat_triface: step %ld molA %d flipped 2->4\n",
              update->ntimestep, molA);
    }
  }
}

/* ---------------------------------------------------------------------- */



