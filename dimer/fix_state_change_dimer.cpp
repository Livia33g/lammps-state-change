/* -*- c++ -*- ----------------------------------------------------------*/
/* Minimal state/change fix for dimer system                              */
/* - Every `nevery` steps, if any Red patch (type 2) is within `cutoff`   */
/*   of any Blue patch (type 3), flip that Red molecule's patches 2 -> 3. */
/* - Flip probability `pflip` (0-1).                                      */
/* - Core atoms (type 1) unchanged.                                       */
/*                                                                        */
/* Arguments:                                                             */
/*   fix ID group-ID state/change/dimer nevery cutoff pflip group_patches */
/* Example:                                                               */
/*   fix sc all state/change/dimer 100 1.6 1.0 patches                    */
/*                                                                        */
/* ---------------------------------------------------------------------- */

#include "fix_state_change_dimer.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "math_extra.h"
#include "memory.h"
#include "pair.h"
#include "random_park.h"
#include "update.h"
#include "memory.h"

#include <set>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <cstdio>

using namespace LAMMPS_NS;
using namespace FixConst;

FixStateChangeDimer::FixStateChangeDimer(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp, narg, arg), contact_counter(nullptr), random(nullptr), pair(nullptr) {
  // Args:
  //   fix ID group-ID state/change/dimer nevery cutoff pflip group_patches [hysteresis_checks]
  // where hysteresis_checks requires sustained contact across N consecutive checks.
  if (narg < 7)
    error->all(FLERR,
               "Illegal fix state/change/dimer command "
               "(expected: ID group-ID state/change/dimer nevery cutoff pflip "
               "group_patches)");

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  cutoff = utils::numeric(FLERR, arg[4], false, lmp);
  pflip = utils::numeric(FLERR, arg[5], false, lmp);

  group_patches = group->find(arg[6]);
  if (group_patches < 0)
    error->all(FLERR, "Fix state/change/dimer group_patches not found");

  hysteresis_checks = 1;
  if (narg >= 8) {
    hysteresis_checks = utils::inumeric(FLERR, arg[7], false, lmp);
    if (hysteresis_checks < 1)
      error->all(FLERR, "Illegal hysteresis_checks in state/change/dimer");
  }

  if (nevery <= 0) error->all(FLERR, "Illegal nevery in state/change/dimer");
  if (cutoff <= 0.0) error->all(FLERR, "Illegal cutoff in state/change/dimer");
  if (pflip < 0.0 || pflip > 1.0)
    error->all(FLERR, "Illegal pflip in state/change/dimer");

  cutoffsq = cutoff * cutoff;
  seed = 12345;

  peratom_flag = 0;
  restart_peratom = 0;
  random = new RanPark(lmp, seed);

  // allocate contact counter for all atoms (we'll use it for patch atoms only)
  memory->create(contact_counter, atom->nmax, "state/change/dimer:contact_counter");
  for (int i = 0; i < atom->nmax; i++) contact_counter[i] = 0;
}

int FixStateChangeDimer::setmask() {
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

void FixStateChangeDimer::init() {
  // Get pointer to pair potential for energy calculations
  if (!force->pair)
    error->all(FLERR, "Fix state/change/dimer requires a pair style");
  pair = force->pair;
}

double FixStateChangeDimer::calculate_flip_energy(int mol_id, int old_type, int new_type) {
  // Calculate the potential energy change when flipping a molecule's patches
  // from old_type to new_type.
  //
  // Strategy:
  // 1. Find all patch atoms in this molecule with old_type
  // 2. For each such atom i, iterate over all other atoms j
  // 3. Calculate pair energy before flip: E_old(i,j) with types (old_type, type[j])
  // 4. Calculate pair energy after flip: E_new(i,j) with types (new_type, type[j])
  // 5. Accumulate dE = sum(E_new - E_old)

  int *type = atom->type;
  int *molecule = atom->molecule;
  int *mask = atom->mask;
  double **x = atom->x;
  double *prd = domain->prd;
  int nall = atom->nlocal + atom->nghost;

  // Collect all patch atoms in this molecule that will be changed
  std::vector<int> changed_atoms;
  changed_atoms.reserve(4);

  for (int i = 0; i < nall; i++) {
    if (molecule[i] != mol_id) continue;
    if (!(mask[i] & group->bitmask[group_patches])) continue;
    if (type[i] != old_type) continue;
    changed_atoms.push_back(i);
  }

  if (changed_atoms.empty()) return 0.0;

  // Build a set for fast lookup
  std::unordered_set<int> changed_set(changed_atoms.begin(), changed_atoms.end());

  double dE = 0.0;
  double fforce = 0.0;  // dummy variable for pair->single()

  for (int i : changed_atoms) {
    for (int j = 0; j < nall; j++) {
      if (j == i) continue;

      // If both i and j are in changed_atoms, only count once (when j > i)
      bool j_changed = (changed_set.find(j) != changed_set.end());
      if (j_changed && j < i) continue;

      // Determine old and new types for both atoms
      int it_old = old_type;
      int it_new = new_type;
      int jt_old = type[j];
      int jt_new = type[j];

      if (j_changed) {
        jt_old = old_type;
        jt_new = new_type;
      }

      // Calculate distance with minimum image convention
      double dx = x[j][0] - x[i][0];
      double dy = x[j][1] - x[i][1];
      double dz = x[j][2] - x[i][2];

      if (prd) {
        if (prd[0] > 0.0) dx -= prd[0] * std::round(dx / prd[0]);
        if (prd[1] > 0.0) dy -= prd[1] * std::round(dy / prd[1]);
        if (prd[2] > 0.0) dz -= prd[2] * std::round(dz / prd[2]);
      }

      double rsq = dx * dx + dy * dy + dz * dz;

      // Check cutoffs for both old and new type pairs
      double cut_old = pair->cutsq[it_old][jt_old];
      double cut_new = pair->cutsq[it_new][jt_new];

      if (rsq >= cut_old && rsq >= cut_new) continue;

      // Calculate pair energies
      double e_old = (rsq < cut_old)
          ? pair->single(i, j, it_old, jt_old, rsq, 1.0, 1.0, fforce)
          : 0.0;

      double e_new = (rsq < cut_new)
          ? pair->single(i, j, it_new, jt_new, rsq, 1.0, 1.0, fforce)
          : 0.0;

      dE += (e_new - e_old);
    }
  }

  return dE;
}

void FixStateChangeDimer::detect_and_schedule(int i, std::vector<int> &mol_list) {
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nall = atom->nlocal + atom->nghost;

  // Only Red patches (type 2) in patch group
  if (!(mask[i] & group->bitmask[group_patches])) return;
  if (type[i] != 2) return;

  // check neighbors for Blue patch (type 3)
  double **x = atom->x;
  double *prd = domain->prd;
  bool found = false;
  for (int j = 0; j < nall; j++) {
    if (j == i) continue;
    if (!(mask[j] & group->bitmask[group_patches])) continue;
    if (type[j] != 3) continue;

    // distance
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

void FixStateChangeDimer::post_integrate() {
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

  double **x = atom->x;  // (unused, but keep for consistency)
  int nall = atom->nlocal + atom->nghost;

  // 2) apply flips (only to local atoms of those molecules)
  for (int mol : unique_mols) {
    // probability check per molecule
    double r = random->uniform();
    if (r > pflip) continue;

    // Calculate energy change BEFORE applying the flip
    const double dE = calculate_flip_energy(mol, 2, 3);

    // Apply the flip
    for (int i = 0; i < nlocal; i++) {
      if (molecule[i] != mol) continue;
      if (!(mask[i] & group->bitmask[group_patches])) continue;
      if (type[i] == 2) type[i] = 3;  // flip Red patch to Blue
    }

    // Print state change with instantaneous energy change
    if (comm->me == 0) {
      fprintf(stderr, "STATECHANGE dimer: step %ld mol %d flipped 2->3 dE %.15g\n",
              update->ntimestep, mol, dE);
    }
  }
}

/* ---------------------------------------------------------------------- */