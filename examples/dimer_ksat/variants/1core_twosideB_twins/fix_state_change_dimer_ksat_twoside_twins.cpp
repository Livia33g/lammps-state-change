/* -*- c++ -*- ----------------------------------------------------------*/
/* Minimal state-change fix for "dimer_ksat_twoside_twins". */

#include "fix_state_change_dimer_ksat_twoside_twins.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "pair.h"
#include "random_park.h"
#include "update.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;

namespace {
// Channel 1
constexpr int TYPE_A = 2;
constexpr int TYPE_B_F1 = 3;
constexpr int TYPE_C = 4;
constexpr int TYPE_B_F2 = 5;

// Channel 2
constexpr int TYPE_E = 8;
constexpr int TYPE_D_F1 = 9;
constexpr int TYPE_F = 10;
constexpr int TYPE_D_F2 = 11;

struct MinPair {
  int m1 = 0;
  int m2 = 0;
};

inline void min_set(int &cur, int val) {
  if (val <= 0) return;
  cur = (cur == 0) ? val : std::min(cur, val);
}

inline double min_image_rsq(double **x, int i, int j, double *prd) {
  double dx = x[j][0] - x[i][0];
  double dy = x[j][1] - x[i][1];
  double dz = x[j][2] - x[i][2];
  if (prd) {
    if (prd[0] > 0.0) dx -= prd[0] * std::round(dx / prd[0]);
    if (prd[1] > 0.0) dy -= prd[1] * std::round(dy / prd[1]);
    if (prd[2] > 0.0) dz -= prd[2] * std::round(dz / prd[2]);
  }
  return dx * dx + dy * dy + dz * dz;
}

double delta_energy_for_molecule_typeflip(
    Pair *pair, double **x, double *prd, int nall,
    int *type, int *molecule, int *mask, int bit_patches,
    int mol_id, int old_patch_type, int new_patch_type) {
  // Collect all patch atoms (local+ghost) in this molecule that have old_patch_type.
  std::vector<int> changed;
  changed.reserve(8);
  for (int i = 0; i < nall; i++) {
    if (molecule[i] != mol_id) continue;
    if (!(mask[i] & bit_patches)) continue;
    if (type[i] != old_patch_type) continue;
    changed.push_back(i);
  }
  if (changed.empty()) return 0.0;

  std::unordered_set<int> changed_set;
  changed_set.reserve(changed.size() * 2);
  for (int i : changed) changed_set.insert(i);

  double dE = 0.0;
  double fforce = 0.0;

  for (int i : changed) {
    for (int j = 0; j < nall; j++) {
      if (j == i) continue;

      const bool j_changed = (changed_set.find(j) != changed_set.end());
      if (j_changed && j < i) continue;  // count within-changed pairs once

      // Determine old/new types for both endpoints
      const int it_old = old_patch_type;
      const int it_new = new_patch_type;
      int jt_old = type[j];
      int jt_new = type[j];
      if (j_changed) {
        jt_old = old_patch_type;
        jt_new = new_patch_type;
      }

      const double rsq = min_image_rsq(x, i, j, prd);
      const double cut_old = pair->cutsq[it_old][jt_old];
      const double cut_new = pair->cutsq[it_new][jt_new];
      if (rsq >= cut_old && rsq >= cut_new) continue;

      const double e_old =
          (rsq < cut_old) ? pair->single(i, j, it_old, jt_old, rsq, 1.0, 1.0, fforce) : 0.0;
      const double e_new =
          (rsq < cut_new) ? pair->single(i, j, it_new, jt_new, rsq, 1.0, 1.0, fforce) : 0.0;

      dE += (e_new - e_old);
    }
  }

  return dE;
}
}  // namespace

FixStateChangeDimerKsatTwoSideTwins::FixStateChangeDimerKsatTwoSideTwins(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp, narg, arg), random(nullptr) {
  if (narg < 7)
    error->all(FLERR,
               "Illegal fix state/change/dimer_ksat_twoside_twins command "
               "(expected: ID group-ID state/change/dimer_ksat_twoside_twins nevery cutoff pflip "
               "group_patches)");

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  cutoff = utils::numeric(FLERR, arg[4], false, lmp);
  pflip = utils::numeric(FLERR, arg[5], false, lmp);

  group_patches = group->find(arg[6]);
  if (group_patches < 0)
    error->all(FLERR, "Fix state/change/dimer_ksat_twoside_twins group_patches not found");

  hysteresis_checks = 1;
  if (narg >= 8) {
    hysteresis_checks = utils::inumeric(FLERR, arg[7], false, lmp);
    if (hysteresis_checks < 1)
      error->all(FLERR, "Illegal hysteresis_checks in state/change/dimer_ksat_twoside_twins");
  }

  if (nevery <= 0) error->all(FLERR, "Illegal nevery in state/change/dimer_ksat_twoside_twins");
  if (cutoff <= 0.0) error->all(FLERR, "Illegal cutoff in state/change/dimer_ksat_twoside_twins");
  if (pflip < 0.0 || pflip > 1.0)
    error->all(FLERR, "Illegal pflip in state/change/dimer_ksat_twoside_twins");

  cutoffsq = cutoff * cutoff;
  seed = 12345;
  peratom_flag = 0;
  restart_peratom = 0;
  random = new RanPark(lmp, seed);
}

int FixStateChangeDimerKsatTwoSideTwins::setmask() {
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

void FixStateChangeDimerKsatTwoSideTwins::post_integrate() {
  if (update->ntimestep % nevery) return;

  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  double **x = atom->x;
  double *prd = domain->prd;
  Pair *pair = force->pair;
  if (!pair) error->all(FLERR, "Fix state/change/dimer_ksat_twoside_twins requires a pair style");

  const int bit_patches = group->bitmask[group_patches];

  auto in_patch_group = [&](int i) -> bool {
    return (mask[i] & bit_patches);
  };

  auto min_image_rsq_local = [&](int i, int j) -> double { return min_image_rsq(x, i, j, prd); };

  // Collect per-B and per-D minima
  std::unordered_map<int, MinPair> bmins;
  std::unordered_map<int, MinPair> dmins;
  bmins.reserve(64);
  dmins.reserve(64);

  // 1) scan face patches, record min attached molecule id for each face
  for (int i = 0; i < nlocal; i++) {
    if (!in_patch_group(i)) continue;
    const int t = type[i];
    const int molX = molecule[i];
    if (molX <= 0) continue;

    // B faces: look for A patches nearby
    if (t == TYPE_B_F1 || t == TYPE_B_F2) {
      MinPair &mp = bmins[molX];
      for (int j = 0; j < nall; j++) {
        if (j == i) continue;
        if (!in_patch_group(j)) continue;
        if (type[j] != TYPE_A) continue;
        if (min_image_rsq_local(i, j) < cutoffsq) {
          int molA = molecule[j];
          if (molA <= 0) continue;
          if (t == TYPE_B_F1) min_set(mp.m1, molA);
          else min_set(mp.m2, molA);
        }
      }
      continue;
    }

    // D faces: look for E patches nearby
    if (t == TYPE_D_F1 || t == TYPE_D_F2) {
      MinPair &mp = dmins[molX];
      for (int j = 0; j < nall; j++) {
        if (j == i) continue;
        if (!in_patch_group(j)) continue;
        if (type[j] != TYPE_E) continue;
        if (min_image_rsq_local(i, j) < cutoffsq) {
          int molE = molecule[j];
          if (molE <= 0) continue;
          if (t == TYPE_D_F1) min_set(mp.m1, molE);
          else min_set(mp.m2, molE);
        }
      }
      continue;
    }
  }

  // 2) hysteresis and scheduling
  std::set<int> flip_A_to_C;
  std::set<int> flip_E_to_F;

  for (const auto &kv : bmins) {
    const int molB = kv.first;
    const MinPair &mp = kv.second;
    const bool has_both = (mp.m1 > 0 && mp.m2 > 0 && mp.m1 != mp.m2);
    int &count = mol_counter[molB];
    if (has_both) count += 1;
    else count = 0;
    if (has_both && count >= hysteresis_checks) {
      flip_A_to_C.insert(std::min(mp.m1, mp.m2));
      count = 0;
    }
  }

  for (const auto &kv : dmins) {
    const int molD = kv.first;
    const MinPair &mp = kv.second;
    const bool has_both = (mp.m1 > 0 && mp.m2 > 0 && mp.m1 != mp.m2);
    int &count = mol_counter[molD];
    if (has_both) count += 1;
    else count = 0;
    if (has_both && count >= hysteresis_checks) {
      flip_E_to_F.insert(std::min(mp.m1, mp.m2));
      count = 0;
    }
  }

  if (flip_A_to_C.empty() && flip_E_to_F.empty()) return;

  // 3) apply flips (local atoms only)
  for (int molA : flip_A_to_C) {
    double r = random->uniform();
    if (r > pflip) continue;

    const double dE = delta_energy_for_molecule_typeflip(
        pair, x, prd, nall, type, molecule, mask, bit_patches, molA, TYPE_A, TYPE_C);
    for (int i = 0; i < nlocal; i++) {
      if (molecule[i] != molA) continue;
      if (!in_patch_group(i)) continue;
      if (type[i] == TYPE_A) type[i] = TYPE_C;
    }
    if (comm->me == 0) {
      fprintf(stderr,
              "STATECHANGE dimer_ksat_twoside_twins: step %ld molA %d flipped 2->4 dE %.15g\n",
              update->ntimestep, molA, dE);
    }
  }

  for (int molE : flip_E_to_F) {
    double r = random->uniform();
    if (r > pflip) continue;

    const double dE = delta_energy_for_molecule_typeflip(
        pair, x, prd, nall, type, molecule, mask, bit_patches, molE, TYPE_E, TYPE_F);
    for (int i = 0; i < nlocal; i++) {
      if (molecule[i] != molE) continue;
      if (!in_patch_group(i)) continue;
      if (type[i] == TYPE_E) type[i] = TYPE_F;
    }
    if (comm->me == 0) {
      fprintf(stderr,
              "STATECHANGE dimer_ksat_twoside_twins: step %ld molE %d flipped 8->10 dE %.15g\n",
              update->ntimestep, molE, dE);
    }
  }
}

/* ---------------------------------------------------------------------- */



