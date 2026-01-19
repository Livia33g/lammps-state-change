/* -*- c++ -*- ----------------------------------------------------------*/
/* Minimal state-change fix for "dimer_ksat" system:
 * - Every `nevery` steps, if any A patch (type 2) is within `cutoff`
 *   of any B patch (type 3), flip that A molecule's patches 2 -> 4 (C).
 * - Flip probability `pflip` (0-1) applied per molecule.
 * - Core atoms (type 1) unchanged.
 *
 * Arguments:
 *   fix ID group-ID state/change/dimer_ksat nevery cutoff pflip group_patches [hysteresis_checks]
 *
 * Example:
 *   fix sc patches state/change/dimer_ksat 100 1.6 1.0 patches 5
 * ---------------------------------------------------------------------- */

#include "fix_state_change_dimer_ksat.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "random_park.h"
#include "update.h"

#include <cmath>
#include <cstdio>
#include <set>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;

FixStateChangeDimerKsat::FixStateChangeDimerKsat(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp, narg, arg), contact_counter(nullptr), random(nullptr) {
  // Args:
  //   fix ID group-ID state/change/dimer_ksat nevery cutoff pflip group_patches [hysteresis_checks]
  if (narg < 7)
    error->all(FLERR,
               "Illegal fix state/change/dimer_ksat command "
               "(expected: ID group-ID state/change/dimer_ksat nevery cutoff pflip "
               "group_patches)");

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  cutoff = utils::numeric(FLERR, arg[4], false, lmp);
  pflip = utils::numeric(FLERR, arg[5], false, lmp);

  group_patches = group->find(arg[6]);
  if (group_patches < 0)
    error->all(FLERR, "Fix state/change/dimer_ksat group_patches not found");

  hysteresis_checks = 1;
  if (narg >= 8) {
    hysteresis_checks = utils::inumeric(FLERR, arg[7], false, lmp);
    if (hysteresis_checks < 1)
      error->all(FLERR, "Illegal hysteresis_checks in state/change/dimer_ksat");
  }

  if (nevery <= 0) error->all(FLERR, "Illegal nevery in state/change/dimer_ksat");
  if (cutoff <= 0.0) error->all(FLERR, "Illegal cutoff in state/change/dimer_ksat");
  if (pflip < 0.0 || pflip > 1.0)
    error->all(FLERR, "Illegal pflip in state/change/dimer_ksat");

  cutoffsq = cutoff * cutoff;
  seed = 12345;

  peratom_flag = 0;
  restart_peratom = 0;
  random = new RanPark(lmp, seed);

  // allocate contact counter for all atoms (we'll use it for patch atoms only)
  memory->create(contact_counter, atom->nmax, "state/change/dimer_ksat:contact_counter");
  for (int i = 0; i < atom->nmax; i++) contact_counter[i] = 0;
}

int FixStateChangeDimerKsat::setmask() {
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

void FixStateChangeDimerKsat::detect_and_schedule(int i, std::vector<int> &mol_list) {
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nall = atom->nlocal + atom->nghost;

  // Only A patches (type 2) in patch group can switch
  if (!(mask[i] & group->bitmask[group_patches])) return;
  if (type[i] != 2) return;

  // check neighbors for B patch (type 3)
  double **x = atom->x;
  double *prd = domain->prd;
  bool found = false;
  for (int j = 0; j < nall; j++) {
    if (j == i) continue;
    if (!(mask[j] & group->bitmask[group_patches])) continue;
    if (type[j] != 3) continue;

    // distance (minimum image)
    double dx = x[j][0] - x[i][0];
    double dy = x[j][1] - x[i][1];
    double dz = x[j][2] - x[i][2];
    if (prd) {
      if (prd[0] > 0.0) dx -= prd[0] * std::round(dx / prd[0]);
      if (prd[1] > 0.0) dy -= prd[1] * std::round(dy / prd[1]);
      if (prd[2] > 0.0) dz -= prd[2] * std::round(dz / prd[2]);
    }
    double rsq = dx * dx + dy * dy + dz * dz;
    if (rsq < cutoffsq) {
      found = true;
      break;
    }
  }

  // Update hysteresis counter (for this patch atom)
  if (found) {
    if (contact_counter) contact_counter[i] += 1;
  } else {
    if (contact_counter) contact_counter[i] = 0;
  }

  // Only schedule if we've held contact for long enough
  if (found && contact_counter && contact_counter[i] >= hysteresis_checks) {
    int mol = molecule[i];
    if (mol > 0) mol_list.push_back(mol);
    // Reset so we don't immediately reschedule on next check
    contact_counter[i] = 0;
  }
}

void FixStateChangeDimerKsat::post_integrate() {
  // Only act every nevery steps
  if (update->ntimestep % nevery) return;

  int nlocal = atom->nlocal;
  int *molecule = atom->molecule;
  int *type = atom->type;
  int *mask = atom->mask;

  // 1) detect mols to flip
  std::vector<int> to_flip;
  to_flip.reserve(64);

  for (int i = 0; i < nlocal; i++) {
    detect_and_schedule(i, to_flip);
  }

  if (to_flip.empty()) return;

  // Unique mol list
  std::set<int> unique_mols(to_flip.begin(), to_flip.end());

  // 2) apply flips (only to local atoms of those molecules)
  for (int mol : unique_mols) {
    // probability check per molecule
    double r = random->uniform();
    if (r > pflip) continue;

    for (int i = 0; i < nlocal; i++) {
      if (molecule[i] != mol) continue;
      if (!(mask[i] & group->bitmask[group_patches])) continue;
      if (type[i] == 2) type[i] = 4;  // A -> C
    }
    if (comm->me == 0) {
      fprintf(stderr, "STATECHANGE dimer_ksat: step %ld mol %d flipped 2->4\n",
              update->ntimestep, mol);
    }
  }
}

/* ---------------------------------------------------------------------- */





