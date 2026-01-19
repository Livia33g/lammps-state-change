/* -*- c++ -*- ----------------------------------------------------------*/
/* Minimal state-change fix for "dimer_ksat_threeface_bd".
 *
 * See header for the rule details.
 * ---------------------------------------------------------------------- */

#include "fix_state_change_dimer_ksat_threeface_bd.h"

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
#include <unordered_map>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;

namespace {
constexpr int TYPE_A_PATCH = 2;
constexpr int TYPE_C_PATCH = 4;

constexpr int TYPE_FACE0 = 3;
constexpr int TYPE_FACE1 = 5;
constexpr int TYPE_FACE2 = 6;

constexpr int TYPE_CORE_B = 7;
constexpr int TYPE_CORE_D = 8;

struct MolInfo {
  // two smallest distinct A molecule IDs attached anywhere
  int a_min1 = 0;
  int a_min2 = 0;
  // face occupancy (A present)
  bool face_has[3] = {false, false, false};
  // what kind of monomer this is
  bool is_B = false;
  bool is_D = false;
};

inline void update_two_mins(int molA, int &min1, int &min2) {
  if (molA <= 0) return;
  if (min1 == 0 || molA < min1) {
    if (molA != min1) min2 = min1;
    min1 = molA;
    return;
  }
  if (molA == min1) return;
  if (min2 == 0 || molA < min2) min2 = molA;
}
}  // namespace

FixStateChangeDimerKsatThreeFaceBD::FixStateChangeDimerKsatThreeFaceBD(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp, narg, arg), random(nullptr) {
  if (narg < 7)
    error->all(FLERR,
               "Illegal fix state/change/dimer_ksat_threeface_bd command "
               "(expected: ID group-ID state/change/dimer_ksat_threeface_bd nevery cutoff pflip "
               "group_patches)");

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  cutoff = utils::numeric(FLERR, arg[4], false, lmp);
  pflip = utils::numeric(FLERR, arg[5], false, lmp);

  group_patches = group->find(arg[6]);
  if (group_patches < 0)
    error->all(FLERR, "Fix state/change/dimer_ksat_threeface_bd group_patches not found");

  hysteresis_checks = 1;
  if (narg >= 8) {
    hysteresis_checks = utils::inumeric(FLERR, arg[7], false, lmp);
    if (hysteresis_checks < 1)
      error->all(FLERR, "Illegal hysteresis_checks in state/change/dimer_ksat_threeface_bd");
  }

  if (nevery <= 0) error->all(FLERR, "Illegal nevery in state/change/dimer_ksat_threeface_bd");
  if (cutoff <= 0.0) error->all(FLERR, "Illegal cutoff in state/change/dimer_ksat_threeface_bd");
  if (pflip < 0.0 || pflip > 1.0)
    error->all(FLERR, "Illegal pflip in state/change/dimer_ksat_threeface_bd");

  cutoffsq = cutoff * cutoff;
  seed = 12345;
  peratom_flag = 0;
  restart_peratom = 0;
  random = new RanPark(lmp, seed);
}

int FixStateChangeDimerKsatThreeFaceBD::setmask() {
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

void FixStateChangeDimerKsatThreeFaceBD::post_integrate() {
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

  // Map (B/D molecule id) -> collected info
  std::unordered_map<int, MolInfo> info;
  info.reserve(128);

  // 0) Identify which molecules are B vs D by core type.
  // Scan local+ghost so we still recognize B/D if the core is a ghost on this rank.
  for (int i = 0; i < nall; i++) {
    int mol = molecule[i];
    if (mol <= 0) continue;
    if (type[i] == TYPE_CORE_B) {
      info[mol].is_B = true;
    } else if (type[i] == TYPE_CORE_D) {
      info[mol].is_D = true;
    }
  }

  // 1) For each B/D face patch (types 3/5/6), search A patches within cutoff.
  for (int i = 0; i < nlocal; i++) {
    if (!in_patch_group(i)) continue;
    int t = type[i];
    int face = -1;
    if (t == TYPE_FACE0) face = 0;
    else if (t == TYPE_FACE1) face = 1;
    else if (t == TYPE_FACE2) face = 2;
    else continue;

    int molBD = molecule[i];
    if (molBD <= 0) continue;
    auto it = info.find(molBD);
    if (it == info.end()) continue;  // not a B/D molecule (or core not local)
    MolInfo &mi = it->second;
    if (!(mi.is_B || mi.is_D)) continue;

    for (int j = 0; j < nall; j++) {
      if (j == i) continue;
      if (!in_patch_group(j)) continue;
      if (type[j] != TYPE_A_PATCH) continue;
      if (min_image_rsq(i, j) >= cutoffsq) continue;

      int molA = molecule[j];
      if (molA <= 0) continue;

      mi.face_has[face] = true;
      update_two_mins(molA, mi.a_min1, mi.a_min2);
    }
  }

  // 2) Apply hysteresis and schedule flips
  std::vector<int> to_flip;
  to_flip.reserve(64);

  for (auto &kv : info) {
    int molBD = kv.first;
    MolInfo &mi = kv.second;
    if (!(mi.is_B || mi.is_D)) continue;

    bool cond = false;
    if (mi.is_B) {
      // B flips if >=2 A molecules attached anywhere
      cond = (mi.a_min2 > 0);
    } else if (mi.is_D) {
      // D flips if all 3 faces have at least one A attached
      cond = (mi.face_has[0] && mi.face_has[1] && mi.face_has[2]);
    }

    int &count = mol_counter[molBD];
    if (cond) count += 1;
    else count = 0;

    if (cond && count >= hysteresis_checks) {
      if (mi.a_min1 > 0) to_flip.push_back(mi.a_min1);
      count = 0;
    }
  }

  if (to_flip.empty()) return;

  // 3) Apply flips (2->4) for scheduled A molecules (local atoms only)
  for (int molA : to_flip) {
    double r = random->uniform();
    if (r > pflip) continue;

    for (int i = 0; i < nlocal; i++) {
      if (molecule[i] != molA) continue;
      if (!in_patch_group(i)) continue;
      if (type[i] == TYPE_A_PATCH) type[i] = TYPE_C_PATCH;
    }

    if (comm->me == 0) {
      fprintf(stderr,
              "STATECHANGE dimer_ksat_threeface_bd: step %ld molA %d flipped 2->4\n",
              update->ntimestep, molA);
    }
  }
}

/* ---------------------------------------------------------------------- */


