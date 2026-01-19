#!/usr/bin/env python3
"""
Generate data/in files for a 1-core "three-face" B/D simulation.

Monomer geometries:
- A and C: same as the working 1-core monomer
    core type 1 at (0,0,0)
    3 patches on +x face (type 2 for A, type 4 for C)

- B and D: same core size, but THREE patch faces (3 patches each = 9 patches total)
    core type 7 for B, core type 8 for D  (cores are energetically identical)
    face types (energetically identical):
      +x face: type 3
      -x face: type 5
      mid face (~45deg): type 6

Rules (implemented in C++ fix state/change/dimer_ksat_threeface_bd):
- B: if >=2 distinct A molecules attached anywhere, flip the lowest-ID A to C.
- D: if all 3 faces have at least one A attached, flip the lowest-ID A to C.

Potentials / parameters: match generate_dimer_ksat_1core_cpp.py
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


def rotate_y(point, theta_rad):
    x, y, z = point
    c = math.cos(theta_rad)
    s = math.sin(theta_rad)
    return (x * c + z * s, y, -x * s + z * c)


def main():
    out_dir = Path("simulation_cpp")
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- composition (requested: equal B and D; equal A and C initially) ---
    # Default keeps total monomers = 30 like prior runs.
    n_A = 10
    n_C = 10
    n_B = 5
    n_D = 5
    n_mol = n_A + n_B + n_C + n_D

    # --- masses (same totals as 1-core baseline) ---
    mass_patch = 0.1
    mass_core_AC = 0.6  # A/C core (type 1)
    mass_core_BD = 0.6  # B core (type 7) and D core (type 8)

    # --- interaction parameters (match 1-core baseline) ---
    CORE_DIAMETER = 1.0
    SOFT_A = 500.0

    PATCH_RADIUS = 0.1
    MORSE_RCUT = 7.0 * PATCH_RADIUS  # 0.7
    CONTACT_CUTOFF = MORSE_RCUT

    D0_BIND = 1.0
    MORSE_ALPHA = 5.0
    MORSE_R0 = 0.0

    # --- concentration / box size (same as 1-core baseline) ---
    c_total = 0.001
    volume = n_mol / c_total
    box_l_nominal = math.pow(volume, 1 / 3)
    BOX_SHRINK = 0.5
    box_l = box_l_nominal * BOX_SHRINK

    # --- geometry ---
    core_AC = [(0.0, 0.0, 0.0)]
    core_BD = [(0.0, 0.0, 0.0)]

    patches_pos = [
        (0.5, 0.0, 0.1),
        (0.5, 0.0866025404, -0.05),
        (0.5, -0.0866025404, -0.05),
    ]
    patches_neg = [(-x, y, -z) for (x, y, z) in patches_pos]  # 180° about y
    patches_mid = [rotate_y(p, math.pi / 4.0) for p in patches_pos]  # +45° about y

    atoms_A = 1 + 3
    atoms_C = 1 + 3
    atoms_B = 1 + 9
    atoms_D = 1 + 9
    natoms = n_A * atoms_A + n_B * atoms_B + n_C * atoms_C + n_D * atoms_D

    positions = make_positions(n_mol, box_l, seed=1234)
    rng = random.Random(1234)
    rng.shuffle(positions)

    # --- write data file ---
    data_path = out_dir / "data.dimer_ksat_1core_threefaceBD"
    with data_path.open("w") as f:
        f.write("LAMMPS data file for rigid dimer_ksat 1-core three-face B/D\n\n")
        f.write(f"{natoms} atoms\n")
        f.write("0 bonds\n0 angles\n0 dihedrals\n0 impropers\n\n")
        f.write("8 atom types\n\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} xlo xhi\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} ylo yhi\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} zlo zhi\n\n")
        f.write("Masses\n\n")
        f.write(f"1 {mass_core_AC:.6f}\n")
        f.write(f"2 {mass_patch:.6f}\n")
        f.write(f"3 {mass_patch:.6f}\n")
        f.write(f"4 {mass_patch:.6f}\n")
        f.write(f"5 {mass_patch:.6f}\n")
        f.write(f"6 {mass_patch:.6f}\n")
        f.write(f"7 {mass_core_BD:.6f}\n")
        f.write(f"8 {mass_core_BD:.6f}\n\n")
        f.write("Atoms # full\n\n")

        atom_id = 1
        mol_id = 1

        # Assign shuffled positions: A, B, C, D blocks (counts fixed)
        # (We don't care about spatial correlation because positions were shuffled.)
        kinds = (["A"] * n_A) + (["B"] * n_B) + (["C"] * n_C) + (["D"] * n_D)
        assert len(kinds) == n_mol
        rng.shuffle(kinds)

        for (kind, (cx, cy, cz)) in zip(kinds, positions):
            if kind == "A":
                core_type = 1
                patch_type = 2
                # core
                for (dx, dy, dz) in core_AC:
                    f.write(f"{atom_id} {mol_id} {core_type} 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
                # patches (+x only)
                for (dx, dy, dz) in patches_pos:
                    f.write(f"{atom_id} {mol_id} {patch_type} 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1

            elif kind == "C":
                core_type = 1
                patch_type = 4
                for (dx, dy, dz) in core_AC:
                    f.write(f"{atom_id} {mol_id} {core_type} 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
                for (dx, dy, dz) in patches_pos:
                    f.write(f"{atom_id} {mol_id} {patch_type} 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1

            elif kind == "B":
                core_type = 7
                for (dx, dy, dz) in core_BD:
                    f.write(f"{atom_id} {mol_id} {core_type} 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
                # face +x: type 3
                for (dx, dy, dz) in patches_pos:
                    f.write(f"{atom_id} {mol_id} 3 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
                # face -x: type 5
                for (dx, dy, dz) in patches_neg:
                    f.write(f"{atom_id} {mol_id} 5 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
                # mid face: type 6
                for (dx, dy, dz) in patches_mid:
                    f.write(f"{atom_id} {mol_id} 6 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1

            else:  # "D"
                core_type = 8
                for (dx, dy, dz) in core_BD:
                    f.write(f"{atom_id} {mol_id} {core_type} 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
                for (dx, dy, dz) in patches_pos:
                    f.write(f"{atom_id} {mol_id} 3 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
                for (dx, dy, dz) in patches_neg:
                    f.write(f"{atom_id} {mol_id} 5 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
                for (dx, dy, dz) in patches_mid:
                    f.write(f"{atom_id} {mol_id} 6 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1

            mol_id += 1

    # --- write input script ---
    in_path = out_dir / "in.dimer_ksat_1core_threefaceBD"
    with in_path.open("w") as f:
        f.write("# dimer_ksat 1-core three-face B/D simulation\n")
        f.write("units           lj\n")
        f.write("atom_style      full\n")
        f.write("boundary        p p p\n\n")
        f.write("read_data       data.dimer_ksat_1core_threefaceBD\n\n")

        f.write("# Groups\n")
        f.write("group cores_AC type 1\n")
        f.write("group cores_BD type 7 8\n")
        f.write("group cores union cores_AC cores_BD\n")
        f.write("group patches type 2 3 4 5 6\n\n")

        f.write("# Pair potentials (match 1-core baseline)\n")
        f.write(f"pair_style hybrid morse {MORSE_RCUT:.6f} soft {CORE_DIAMETER:.6f}\n")

        # Core-core repulsion for all core type combinations
        core_types = [1, 7, 8]
        for i, t1 in enumerate(core_types):
            for t2 in core_types[i:]:
                f.write(f"pair_coeff {t1} {t2} soft {SOFT_A:.6f} {CORE_DIAMETER:.6f}\n")

        # Neutral morse for all non-core-core pairs (D0=0)
        # Core-patch (1/7/8 with 2/3/4/5/6)
        for tcore in core_types:
            for tpatch in [2, 3, 4, 5, 6]:
                f.write(f"pair_coeff {tcore} {tpatch} morse 0.0 {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")

        # Patch-patch neutral by default
        patch_types = [2, 3, 4, 5, 6]
        for i, t1 in enumerate(patch_types):
            for t2 in patch_types[i:]:
                f.write(f"pair_coeff {t1} {t2} morse 0.0 {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")

        # Attractions: A/C to B/D faces (all faces identical)
        for tf in [3, 5, 6]:
            f.write(f"pair_coeff 2 {tf} morse {D0_BIND:.6f} {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
            f.write(f"pair_coeff 4 {tf} morse {D0_BIND:.6f} {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
        f.write("\n")

        f.write("# Neighbor settings\n")
        f.write("neighbor        1.0 bin\n")
        f.write("neigh_modify    every 1 delay 0 check yes\n\n")

        f.write("# Rigid body integration + Langevin thermostat\n")
        f.write("group all_monomers molecule > 0\n")
        f.write("fix fx_nve all_monomers rigid/nve molecule\n")
        f.write("fix fx_langevin all_monomers langevin 0.5 0.5 0.5 12345\n\n")

        f.write("# State change\n")
        f.write(f"fix sc patches state/change/dimer_ksat_threeface_bd 100 {CONTACT_CUTOFF:.6f} 1.0 patches 5\n\n")

        f.write("# Output\n")
        f.write("thermo_style    custom step temp pe ke etotal\n")
        f.write("thermo          1000\n")
        f.write("dump            d1 all custom 100 dump.dimer_ksat_1core_threefaceBD.lammpstrj id mol type x y z\n")
        f.write("dump_modify     d1 sort id\n\n")

        f.write("timestep        0.005\n")
        f.write("run             2000000\n")

    print(f"Wrote {data_path} and {in_path}")


if __name__ == "__main__":
    main()


