#!/usr/bin/env python3
"""
Generate data/in files for a 1-core "tri-face" B/D simulation.

Monomer types:
- A: 1 core (type 1) + 3 patches on +x (type 2)
- C: same geometry as A, patches type 4
- B: 1 core (type 1) + 9 patches (3 faces × 3 patches), faces are:
    +x face patches: type 3
    -x face patches: type 5
    "mid" face patches: type 6 (placed between +x and -x faces)
- D: identical to B in patch layout/types, but core is type 7 (so the fix can apply a different rule)

State-change logic:
Use `fix state/change/dimer_ksat_triface`:
- B triggers when >=2 distinct A are attached (any faces).
- D triggers only when all 3 faces are occupied by A (and 3 distinct A).
- Flip lowest-ID A to C (2->4).

Potentials: keep identical to the existing 1-core parameterization.
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


def rot_y_minus_90(v):
    # Rotate around y by -90 degrees: (x,y,z) -> (-z, y, x)
    x, y, z = v
    return (-z, y, x)


def main():
    out_dir = Path("simulation_cpp")
    out_dir.mkdir(parents=True, exist_ok=True)

    # counts: equal A=C and B=D
    n_A = 10
    n_C = 10
    n_B = 5
    n_D = 5
    n_mol = n_A + n_B + n_C + n_D

    # masses (same as 1-core baseline)
    mass_patch = 0.1
    mass_core = 0.9 - 3 * mass_patch  # 0.6

    # interaction parameters (same as 1-core baseline)
    CORE_DIAMETER = 1.0
    SOFT_A = 500.0

    PATCH_RADIUS = 0.1
    MORSE_RCUT = 7.0 * PATCH_RADIUS  # 0.7
    CONTACT_CUTOFF = MORSE_RCUT

    # Slightly stronger attraction (keep identical across faces and for A/C)
    D0_BIND = 1.2
    MORSE_ALPHA = 5.0
    MORSE_R0 = 0.0

    # concentration / box size (same as 1-core baseline)
    c_total = 0.001
    volume = n_mol / c_total
    box_l_nominal = math.pow(volume, 1 / 3)
    BOX_SHRINK = 0.5
    box_l = box_l_nominal * BOX_SHRINK

    # geometry
    core = [(0.0, 0.0, 0.0)]
    patches_pos = [
        (0.5, 0.0, 0.1),
        (0.5, 0.0866025404, -0.05),
        (0.5, -0.0866025404, -0.05),
    ]
    patches_neg = [(-x, y, -z) for (x, y, z) in patches_pos]  # 180 about y
    # "mid" face: rotate +x face by -90deg about y -> points along +z
    patches_mid = [rot_y_minus_90(p) for p in patches_pos]

    atoms_A = 1 + 3
    atoms_C = 1 + 3
    atoms_B = 1 + 9
    atoms_D = 1 + 9
    natoms = n_A * atoms_A + n_C * atoms_C + n_B * atoms_B + n_D * atoms_D

    positions = make_positions(n_mol, box_l, seed=1234)
    random.Random(1234).shuffle(positions)

    # assign molecule kinds across shuffled positions
    kinds = (["A"] * n_A) + (["B"] * n_B) + (["C"] * n_C) + (["D"] * n_D)
    random.Random(4321).shuffle(kinds)
    assert len(kinds) == n_mol

    data_path = out_dir / "data.dimer_ksat_1core_trifaceBD"
    with data_path.open("w") as f:
        f.write("LAMMPS data file for dimer_ksat 1-core with tri-face B/D\n\n")
        f.write(f"{natoms} atoms\n")
        f.write("0 bonds\n0 angles\n0 dihedrals\n0 impropers\n\n")
        f.write("7 atom types\n\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} xlo xhi\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} ylo yhi\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} zlo zhi\n\n")
        f.write("Masses\n\n")
        f.write(f"1 {mass_core:.6f}\n")  # core (A/B/C)
        f.write(f"2 {mass_patch:.6f}\n")  # A patches
        f.write(f"3 {mass_patch:.6f}\n")  # B/D face +x
        f.write(f"4 {mass_patch:.6f}\n")  # C patches
        f.write(f"5 {mass_patch:.6f}\n")  # B/D face -x
        f.write(f"6 {mass_patch:.6f}\n")  # B/D face mid
        f.write(f"7 {mass_core:.6f}\n\n")  # D core
        f.write("Atoms # full\n\n")

        atom_id = 1
        mol_id = 1
        for (cx, cy, cz), kind in zip(positions, kinds):
            # core type
            core_type = 7 if kind == "D" else 1
            f.write(f"{atom_id} {mol_id} {core_type} 0.0 {cx:.6f} {cy:.6f} {cz:.6f}\n")
            atom_id += 1

            if kind == "A":
                for (dx, dy, dz) in patches_pos:
                    f.write(f"{atom_id} {mol_id} 2 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
            elif kind == "C":
                for (dx, dy, dz) in patches_pos:
                    f.write(f"{atom_id} {mol_id} 4 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
            else:
                # B or D: 3 faces (types 3,6,5)
                for (dx, dy, dz) in patches_pos:
                    f.write(f"{atom_id} {mol_id} 3 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
                for (dx, dy, dz) in patches_mid:
                    f.write(f"{atom_id} {mol_id} 6 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
                for (dx, dy, dz) in patches_neg:
                    f.write(f"{atom_id} {mol_id} 5 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1

            mol_id += 1

    in_path = out_dir / "in.dimer_ksat_1core_trifaceBD"
    with in_path.open("w") as f:
        f.write("# dimer_ksat 1-core tri-face B/D simulation\n")
        f.write("units           lj\n")
        f.write("atom_style      full\n")
        f.write("boundary        p p p\n\n")
        f.write("read_data       data.dimer_ksat_1core_trifaceBD\n\n")

        f.write("# Groups\n")
        f.write("group cores type 1 7\n")
        f.write("group patches type 2 3 4 5 6\n\n")

        f.write("# Pair potentials\n")
        f.write(f"pair_style hybrid morse {MORSE_RCUT:.6f} soft {CORE_DIAMETER:.6f}\n")
        # core-core repulsion for core types 1 and 7
        f.write(f"pair_coeff 1 1 soft {SOFT_A:.6f} {CORE_DIAMETER:.6f}\n")
        f.write(f"pair_coeff 1 7 soft {SOFT_A:.6f} {CORE_DIAMETER:.6f}\n")
        f.write(f"pair_coeff 7 7 soft {SOFT_A:.6f} {CORE_DIAMETER:.6f}\n")

        # Neutral morse (D0=0) for everything else by default
        for c in (1, 7):
            for p in (2, 3, 4, 5, 6):
                f.write(f"pair_coeff {c} {p} morse 0.0 {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
        for (t1, t2) in (
            (2, 2), (2, 4), (2, 3), (2, 5), (2, 6),
            (3, 3), (3, 4), (3, 5), (3, 6),
            (4, 4), (4, 5), (4, 6),
            (5, 5), (5, 6),
            (6, 6),
        ):
            f.write(f"pair_coeff {t1} {t2} morse 0.0 {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")

        # Attractions: A/C to any face patch type (3,5,6) are identical
        for face in (3, 5, 6):
            f.write(f"pair_coeff 2 {face} morse {D0_BIND:.6f} {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
            f.write(f"pair_coeff 4 {face} morse {D0_BIND:.6f} {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
        f.write("\n")

        f.write("# Neighbor settings\n")
        f.write("neighbor        1.0 bin\n")
        f.write("neigh_modify    every 1 delay 0 check yes\n\n")

        f.write("# Rigid body integration + Langevin thermostat\n")
        f.write("group all_monomers molecule > 0\n")
        f.write("fix fx_nve all_monomers rigid/nve molecule\n")
        f.write("fix fx_langevin all_monomers langevin 0.5 0.5 0.5 12345\n\n")

        f.write("# State change (B vs D logic in fix)\n")
        f.write(f"fix sc patches state/change/dimer_ksat_triface 100 {CONTACT_CUTOFF:.6f} 1.0 patches 5\n\n")

        f.write("# Output\n")
        f.write("thermo_style    custom step temp pe ke etotal\n")
        f.write("thermo          1000\n")
        f.write("dump            d1 all custom 100 dump.dimer_ksat_1core_trifaceBD.lammpstrj id mol type x y z\n")
        f.write("dump_modify     d1 sort id\n\n")

        f.write("timestep        0.005\n")
        f.write("run             2000000\n")

    print(f"Wrote {data_path} and {in_path}")


if __name__ == "__main__":
    main()


