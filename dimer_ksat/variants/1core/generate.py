#!/usr/bin/env python3
"""
Generate data/in files for the rigid "dimer_ksat" state-change simulation,
but with a *single* repulsive core site per monomer (instead of 3 cores).

Geometry is taken from ksat's `part1_3patch.mol` (core + 3 patches):
  core:   (0.0, 0.0, 0.0)
  patch1: (0.5, 0.0, 0.1)
  patch2: (0.5, 0.0866025404, -0.05)
  patch3: (0.5, -0.0866025404, -0.05)

Species / patch types:
- A monomers: patch type 2 (switchable)
- B monomers: patch type 3 (catalyst)
- C monomers: patch type 4 (non-switchable; only appears via flips)

Desired physics:
- ONLY attractions are between B-A (2-3) and B-C (3-4)
- B-A and B-C attractions are IDENTICAL (same Morse parameters)
- Patches are allowed to "overlap": no patch repulsion; Morse r0=0
- A->C flip triggered by A-B patch contact (custom fix state/change/dimer_ksat)

Key tuning knobs:
- CORE_DIAMETER controls excluded volume of the monomer body (pair_style soft cutoff)
- MORSE_RCUT controls attraction range (pair_style morse cutoff)
- CONTACT_CUTOFF should match MORSE_RCUT for state-change detection
"""

import math
import random
from pathlib import Path


def make_positions(n_mol, box_l, seed=1234):
    """
    Place molecules on a cubic lattice spanning the whole box.

    Important: choose lattice sites spread out across the box (not the first N sites),
    otherwise you visually get "everything in one corner".

    We preferentially pick a checkerboard sublattice (parity constraint) to increase the
    initial nearest-neighbor spacing.
    """
    rng = random.Random(seed)
    per_side = math.ceil(n_mol ** (1 / 3))
    spacing = box_l / per_side
    start = -box_l / 2 + 0.5 * spacing

    all_sites = [(ix, iy, iz) for ix in range(per_side) for iy in range(per_side) for iz in range(per_side)]
    checker_sites = [s for s in all_sites if ((s[0] + s[1] + s[2]) % 2 == 0)]

    if len(checker_sites) >= n_mol:
        sites = checker_sites
    else:
        sites = all_sites

    rng.shuffle(sites)
    sites = sites[:n_mol]

    return [(start + spacing * ix, start + spacing * iy, start + spacing * iz) for (ix, iy, iz) in sites]


def main():
    out_dir = Path("simulation_cpp")
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- composition (requested) ---
    n_A = 20
    n_B = 10
    n_C = 0  # appears only via flips
    n_mol = n_A + n_B + n_C

    # --- masses ---
    # Old 3-core version had total mass 0.9/monomer (3*0.2 + 3*0.1).
    # Keep the same total mass while using 1 core + 3 patches.
    mass_patch = 0.1
    mass_core = 0.9 - 3 * mass_patch  # = 0.6

    # --- geometry (relative to COM) ---
    core = [(0.0, 0.0, 0.0)]
    patches = [
        (0.5, 0.0, 0.1),
        (0.5, 0.0866025404, -0.05),
        (0.5, -0.0866025404, -0.05),
    ]

    atoms_per_mol = 4
    natoms = n_mol * atoms_per_mol

    # --- interaction tuning ---
    # Smaller monomer body: patch anchor is at x=0.5, so choose core radius ~0.5 -> diameter 1.0
    CORE_DIAMETER = 1.0
    SOFT_A = 500.0

    # Patch radius (user-chosen length scale). Patches are allowed to overlap (r0=0),
    # but the interaction range should scale with this radius.
    PATCH_RADIUS = 0.1
    # Attraction/contact range: keep short-range, but long enough for encounters at this density.
    MORSE_RCUT = 7.0 * PATCH_RADIUS  # 0.7
    CONTACT_CUTOFF = MORSE_RCUT

    # Binding depth: keep equal for A–B and C–B (you've been tuning this)
    D0_BIND = 1.0
    MORSE_ALPHA = 5.0
    MORSE_R0 = 0.0

    # concentration -> box size (same scheme as before; monomer concentration)
    # NOTE: "concentration" here is a number density in LJ units.
    # With the smaller 1-core monomer, the same number density can feel too dilute.
    # We therefore allow shrinking the box by a scale factor while keeping the nominal c_total.
    c_total = 0.001
    volume = n_mol / c_total
    box_l_nominal = math.pow(volume, 1 / 3)
    BOX_SHRINK = 0.5  # <1.0 increases effective density; tune as needed
    box_l = box_l_nominal * BOX_SHRINK

    positions = make_positions(n_mol, box_l, seed=1234)
    assert len(positions) == n_mol

    # --- write data file ---
    data_path = out_dir / "data.dimer_ksat_1core"
    with data_path.open("w") as f:
        f.write("LAMMPS data file for rigid dimer_ksat (1 core) state-change\n\n")
        f.write(f"{natoms} atoms\n")
        f.write("0 bonds\n0 angles\n0 dihedrals\n0 impropers\n\n")
        f.write("4 atom types\n\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} xlo xhi\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} ylo yhi\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} zlo zhi\n\n")
        f.write("Masses\n\n")
        f.write(f"1 {mass_core:.6f}\n")
        f.write(f"2 {mass_patch:.6f}\n")
        f.write(f"3 {mass_patch:.6f}\n")
        f.write(f"4 {mass_patch:.6f}\n\n")
        f.write("Atoms # full\n\n")

        atom_id = 1
        mol_id = 1

        # ordering: all A, then all B, then all C
        for idx, (cx, cy, cz) in enumerate(positions):
            if idx < n_A:
                patch_type = 2  # A
            elif idx < n_A + n_B:
                patch_type = 3  # B
            else:
                patch_type = 4  # C

            # core (type 1)
            for (dx, dy, dz) in core:
                f.write(
                    f"{atom_id} {mol_id} 1 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n"
                )
                atom_id += 1

            # patches (type 2/3/4)
            for (dx, dy, dz) in patches:
                f.write(
                    f"{atom_id} {mol_id} {patch_type} 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n"
                )
                atom_id += 1

            mol_id += 1

    # --- write input script ---
    in_path = out_dir / "in.dimer_ksat_1core"
    with in_path.open("w") as f:
        f.write("# dimer_ksat (1 core) state-change simulation (A->C catalyzed by B)\n")
        f.write("units           lj\n")
        f.write("atom_style      full\n")
        f.write("boundary        p p p\n\n")
        f.write("read_data       data.dimer_ksat_1core\n\n")

        f.write("# Groups\n")
        f.write("group cores type 1\n")
        f.write("group patches type 2 3 4\n")
        f.write("group patches_A type 2\n")
        f.write("group patches_B type 3\n")
        f.write("group patches_C type 4\n\n")

        f.write("# Pair potentials\n")
        f.write(f"pair_style hybrid morse {MORSE_RCUT:.6f} soft {CORE_DIAMETER:.6f}\n")
        f.write("# 1-1 soft repulsion (excluded volume)\n")
        f.write(f"pair_coeff 1 1 soft {SOFT_A:.6f} {CORE_DIAMETER:.6f}\n")

        # Neutral interactions (morse D0=0)
        f.write(f"pair_coeff 1 2 morse 0.0 {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
        f.write(f"pair_coeff 1 3 morse 0.0 {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
        f.write(f"pair_coeff 1 4 morse 0.0 {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
        f.write(f"pair_coeff 2 2 morse 0.0 {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
        f.write(f"pair_coeff 2 4 morse 0.0 {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
        f.write(f"pair_coeff 3 3 morse 0.0 {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
        f.write(f"pair_coeff 4 4 morse 0.0 {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")

        # Only attractions: B-A and B-C, identical strength
        f.write(f"pair_coeff 2 3 morse {D0_BIND:.6f} {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
        f.write(f"pair_coeff 3 4 morse {D0_BIND:.6f} {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n\n")

        f.write("# Neighbor settings\n")
        # Use a reasonable skin for the larger Morse cutoff
        f.write("neighbor        1.0 bin\n")
        f.write("neigh_modify    every 1 delay 0 check yes\n\n")

        f.write("# Rigid body integration + Langevin thermostat\n")
        f.write("group all_monomers molecule > 0\n")
        f.write("fix fx_nve all_monomers rigid/nve molecule\n")
        f.write("fix fx_langevin all_monomers langevin 0.5 0.5 0.5 12345\n\n")

        f.write("# State change: A (2) -> C (4) upon A-B contact\n")
        f.write("# Hysteresis: require 5 consecutive checks in contact (5 * 100 steps)\n")
        f.write(f"fix sc patches state/change/dimer_ksat 100 {CONTACT_CUTOFF:.6f} 1.0 patches 5\n\n")

        f.write("# Output\n")
        f.write("thermo_style    custom step temp pe ke etotal\n")
        f.write("thermo          1000\n")
        f.write("dump            d1 all custom 100 dump.dimer_ksat_1core.lammpstrj id mol type x y z\n")
        f.write("dump_modify     d1 sort id\n\n")

        f.write("# Run\n")
        f.write("timestep        0.005\n")
        f.write("run             2000000\n")

    print(f"Wrote {data_path} and {in_path}")


if __name__ == "__main__":
    main()


