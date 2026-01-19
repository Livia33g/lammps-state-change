#!/usr/bin/env python3
"""
Generate data/in files for a "two-side-B" variant of the 1-core dimer_ksat simulation.

Monomer geometries:
- A and C: identical to the working 1-core monomer
    core (type 1): (0,0,0)
    3 patches on +x side:
      (0.5, 0.0,           0.1)
      (0.5, 0.0866025404, -0.05)
      (0.5,-0.0866025404, -0.05)

- B: same core (type 1) but TWO identical patch faces:
    face +x: 3 patches of type 3 at the same coordinates as above
    face -x: 3 patches of type 5 at the coordinates rotated 180° about y:
      (-0.5, 0.0,          -0.1)
      (-0.5, 0.0866025404,  0.05)
      (-0.5,-0.0866025404,  0.05)

State-change rule (implemented by C++ fix):
- A (type-2 patches) flips to C (type-4 patches) iff a B molecule has A bound
  on BOTH faces simultaneously (type-3 and type-5 contacts present).
- Flip the lowest-ID A among the two As (one on each face).

Potentials / parameters:
- Keep identical to generate_dimer_ksat_1core_cpp.py, except we add type 5.
- Ensure interactions are identical for both B faces:
    (2,3) == (2,5) and (4,3) == (4,5)
"""

import math
import random
from pathlib import Path


def make_positions(n_mol, box_l, seed=1234):
    rng = random.Random(seed)
    per_side = math.ceil(n_mol ** (1 / 3))
    spacing = box_l / per_side
    start = -box_l / 2 + 0.5 * spacing

    all_sites = [(ix, iy, iz) for ix in range(per_side) for iy in range(per_side) for iz in range(per_side)]
    checker_sites = [s for s in all_sites if ((s[0] + s[1] + s[2]) % 2 == 0)]
    sites = checker_sites if len(checker_sites) >= n_mol else all_sites
    rng.shuffle(sites)
    sites = sites[:n_mol]
    return [(start + spacing * ix, start + spacing * iy, start + spacing * iz) for (ix, iy, iz) in sites]


def main():
    out_dir = Path("simulation_cpp")
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- composition (same as your standard 1-core runs) ---
    n_A = 20
    n_B = 10
    n_C = 0
    n_mol = n_A + n_B + n_C

    # --- masses (same as 1-core) ---
    mass_patch = 0.1
    mass_core = 0.9 - 3 * mass_patch  # 0.6

    # --- interaction parameters (same as 1-core) ---
    CORE_DIAMETER = 1.0
    SOFT_A = 500.0

    PATCH_RADIUS = 0.1
    MORSE_RCUT = 7.0 * PATCH_RADIUS  # 0.7
    CONTACT_CUTOFF = MORSE_RCUT

    D0_BIND = 1.0
    MORSE_ALPHA = 5.0
    MORSE_R0 = 0.0

    # --- concentration / box size (same as 1-core) ---
    c_total = 0.001
    volume = n_mol / c_total
    box_l_nominal = math.pow(volume, 1 / 3)
    BOX_SHRINK = 0.5
    box_l = box_l_nominal * BOX_SHRINK

    # --- geometry ---
    core = [(0.0, 0.0, 0.0)]
    patches_pos = [
        (0.5, 0.0, 0.1),
        (0.5, 0.0866025404, -0.05),
        (0.5, -0.0866025404, -0.05),
    ]
    patches_neg = [(-x, y, -z) for (x, y, z) in patches_pos]  # rotate 180° about y

    # atoms per molecule differ by species:
    atoms_A = 1 + 3
    atoms_C = 1 + 3
    atoms_B = 1 + 6
    natoms = n_A * atoms_A + n_B * atoms_B + n_C * atoms_C

    positions = make_positions(n_mol, box_l, seed=1234)
    assert len(positions) == n_mol
    random.Random(1234).shuffle(positions)

    # --- write data file ---
    data_path = out_dir / "data.dimer_ksat_1core_twosideB"
    with data_path.open("w") as f:
        f.write("LAMMPS data file for rigid dimer_ksat (1 core, two-side B) state-change\n\n")
        f.write(f"{natoms} atoms\n")
        f.write("0 bonds\n0 angles\n0 dihedrals\n0 impropers\n\n")
        f.write("5 atom types\n\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} xlo xhi\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} ylo yhi\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} zlo zhi\n\n")
        f.write("Masses\n\n")
        f.write(f"1 {mass_core:.6f}\n")
        f.write(f"2 {mass_patch:.6f}\n")
        f.write(f"3 {mass_patch:.6f}\n")
        f.write(f"4 {mass_patch:.6f}\n")
        f.write(f"5 {mass_patch:.6f}\n\n")
        f.write("Atoms # full\n\n")

        atom_id = 1
        mol_id = 1

        # assign first n_A molecules as A, next n_B as B, last n_C as C (positions are shuffled)
        for idx, (cx, cy, cz) in enumerate(positions):
            if idx < n_A:
                kind = "A"
            elif idx < n_A + n_B:
                kind = "B"
            else:
                kind = "C"

            # core
            for (dx, dy, dz) in core:
                f.write(f"{atom_id} {mol_id} 1 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                atom_id += 1

            if kind in ("A", "C"):
                patch_type = 2 if kind == "A" else 4
                for (dx, dy, dz) in patches_pos:
                    f.write(f"{atom_id} {mol_id} {patch_type} 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
            else:
                # B: 3 patches on +x face are type 3, 3 patches on -x face are type 5
                for (dx, dy, dz) in patches_pos:
                    f.write(f"{atom_id} {mol_id} 3 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
                for (dx, dy, dz) in patches_neg:
                    f.write(f"{atom_id} {mol_id} 5 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1

            mol_id += 1

    # --- write input script ---
    in_path = out_dir / "in.dimer_ksat_1core_twosideB"
    with in_path.open("w") as f:
        f.write("# dimer_ksat (1 core, two-side B) state-change simulation\n")
        f.write("units           lj\n")
        f.write("atom_style      full\n")
        f.write("boundary        p p p\n\n")
        f.write("read_data       data.dimer_ksat_1core_twosideB\n\n")

        f.write("# Groups\n")
        f.write("group cores type 1\n")
        f.write("group patches type 2 3 4 5\n")
        f.write("group patches_A type 2\n")
        f.write("group patches_B type 3 5\n")
        f.write("group patches_C type 4\n\n")

        f.write("# Pair potentials (same parameters as 1-core)\n")
        f.write(f"pair_style hybrid morse {MORSE_RCUT:.6f} soft {CORE_DIAMETER:.6f}\n")
        f.write(f"pair_coeff 1 1 soft {SOFT_A:.6f} {CORE_DIAMETER:.6f}\n")

        # Neutral interactions (morse D0=0)
        for t in (2, 3, 4, 5):
            f.write(f"pair_coeff 1 {t} morse 0.0 {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
        for (t1, t2) in ((2, 2), (2, 4), (2, 5), (3, 3), (3, 4), (3, 5), (4, 4), (4, 5), (5, 5)):
            f.write(f"pair_coeff {t1} {t2} morse 0.0 {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")

        # Only attractions: A-B and C-B, identical for both B faces
        f.write(f"pair_coeff 2 3 morse {D0_BIND:.6f} {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
        f.write(f"pair_coeff 2 5 morse {D0_BIND:.6f} {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
        f.write(f"pair_coeff 3 4 morse {D0_BIND:.6f} {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
        f.write(f"pair_coeff 4 5 morse {D0_BIND:.6f} {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n\n")

        f.write("# Neighbor settings\n")
        f.write("neighbor        1.0 bin\n")
        f.write("neigh_modify    every 1 delay 0 check yes\n\n")

        f.write("# Rigid body integration + Langevin thermostat\n")
        f.write("group all_monomers molecule > 0\n")
        f.write("fix fx_nve all_monomers rigid/nve molecule\n")
        f.write("fix fx_langevin all_monomers langevin 0.5 0.5 0.5 12345\n\n")

        f.write("# State change: requires A on both B faces, flips lowest-ID A\n")
        f.write(f"fix sc patches state/change/dimer_ksat_twoside 100 {CONTACT_CUTOFF:.6f} 1.0 patches 5\n\n")

        f.write("# Output\n")
        f.write("thermo_style    custom step temp pe ke etotal\n")
        f.write("thermo          1000\n")
        f.write("dump            d1 all custom 100 dump.dimer_ksat_1core_twosideB.lammpstrj id mol type x y z\n")
        f.write("dump_modify     d1 sort id\n\n")

        f.write("# Run\n")
        f.write("timestep        0.005\n")
        f.write("run             2000000\n")

    print(f"Wrote {data_path} and {in_path}")


if __name__ == "__main__":
    main()
