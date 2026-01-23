#!/usr/bin/env python3
"""
Generate data/in files for a "twins" version of the 1-core two-side-B simulation.

We create two independent channels that do NOT cross-interact:

Channel ABC:
- A patches: type 2
- C patches: type 4
- B two faces: type 3 (+x face) and type 5 (-x face)
- Rule (in fix): if B has A on both faces => flip lower-ID A: 2->4

Channel EFD:
- E patches: type 8
- F patches: type 10
- D two faces: type 9 (+x face) and type 11 (-x face)
- Rule (in fix): if D has E on both faces => flip lower-ID E: 8->10

Important: All attractions are only within-channel, with identical parameters to the baseline.
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

    # --- composition (requested structure: equal pairs) ---
    n_A = 10
    n_C = 10
    n_E = 10
    n_F = 10
    n_B = 5
    n_D = 5
    n_mol = n_A + n_B + n_C + n_D + n_E + n_F

    # --- masses (same as 1-core baseline) ---
    mass_patch = 0.1
    mass_core = 0.9 - 3 * mass_patch  # 0.6

    # --- interaction parameters (same as 1-core baseline) ---
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
    core = [(0.0, 0.0, 0.0)]
    patches_pos = [
        (0.5, 0.0, 0.1),
        (0.5, 0.0866025404, -0.05),
        (0.5, -0.0866025404, -0.05),
    ]
    patches_neg = [(-x, y, -z) for (x, y, z) in patches_pos]  # 180° about y

    atoms_simple = 1 + 3          # A,C,E,F
    atoms_twoface = 1 + 6         # B,D
    natoms = (n_A + n_C + n_E + n_F) * atoms_simple + (n_B + n_D) * atoms_twoface

    positions = make_positions(n_mol, box_l, seed=1234)
    random.Random(1234).shuffle(positions)

    kinds = (["A"] * n_A) + (["C"] * n_C) + (["E"] * n_E) + (["F"] * n_F) + (["B"] * n_B) + (["D"] * n_D)
    random.Random(4321).shuffle(kinds)
    assert len(kinds) == n_mol

    data_path = out_dir / "data.dimer_ksat_1core_twosideB_twins"
    with data_path.open("w") as f:
        f.write("LAMMPS data file for dimer_ksat 1-core two-side B twins (ABC + EFD)\n\n")
        f.write(f"{natoms} atoms\n")
        f.write("0 bonds\n0 angles\n0 dihedrals\n0 impropers\n\n")
        # Types:
        # - Channel ABC cores: type 1
        # - Channel EFD cores: type 12
        # - Patches use types {2,3,4,5,8,9,10,11}
        f.write("12 atom types\n\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} xlo xhi\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} ylo yhi\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} zlo zhi\n\n")

        f.write("Masses\n\n")
        # core types
        f.write(f"1 {mass_core:.6f}\n")
        f.write(f"12 {mass_core:.6f}\n")
        # patch-like types (set same mass)
        for t in range(2, 12):  # 2..11
            f.write(f"{t} {mass_patch:.6f}\n")
        f.write("\n")

        f.write("Atoms # full\n\n")

        atom_id = 1
        mol_id = 1
        for (cx, cy, cz), kind in zip(positions, kinds):
            # Core type: ABC channel uses 1, EFD channel uses 12
            core_type = 12 if kind in ("E", "F", "D") else 1
            f.write(f"{atom_id} {mol_id} {core_type} 0.0 {cx:.6f} {cy:.6f} {cz:.6f}\n")
            atom_id += 1

            if kind == "A":
                patch_type = 2
                for (dx, dy, dz) in patches_pos:
                    f.write(f"{atom_id} {mol_id} {patch_type} 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
            elif kind == "C":
                patch_type = 4
                for (dx, dy, dz) in patches_pos:
                    f.write(f"{atom_id} {mol_id} {patch_type} 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
            elif kind == "E":
                patch_type = 8
                for (dx, dy, dz) in patches_pos:
                    f.write(f"{atom_id} {mol_id} {patch_type} 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
            elif kind == "F":
                patch_type = 10
                for (dx, dy, dz) in patches_pos:
                    f.write(f"{atom_id} {mol_id} {patch_type} 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
            elif kind == "B":
                for (dx, dy, dz) in patches_pos:
                    f.write(f"{atom_id} {mol_id} 3 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
                for (dx, dy, dz) in patches_neg:
                    f.write(f"{atom_id} {mol_id} 5 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
            else:  # D
                for (dx, dy, dz) in patches_pos:
                    f.write(f"{atom_id} {mol_id} 9 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1
                for (dx, dy, dz) in patches_neg:
                    f.write(f"{atom_id} {mol_id} 11 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1

            mol_id += 1

    in_path = out_dir / "in.dimer_ksat_1core_twosideB_twins"
    with in_path.open("w") as f:
        f.write("# dimer_ksat 1-core two-side B twins (ABC + EFD)\n")
        f.write("units           lj\n")
        f.write("atom_style      full\n")
        f.write("boundary        p p p\n\n")
        f.write("read_data       data.dimer_ksat_1core_twosideB_twins\n\n")

        f.write("# Groups\n")
        f.write("group cores type 1 12\n")
        f.write("group patches type 2 3 4 5 8 9 10 11\n\n")

        f.write("# Pair potentials\n")
        f.write(f"pair_style hybrid morse {MORSE_RCUT:.6f} soft {CORE_DIAMETER:.6f}\n")

        # IMPORTANT: In LAMMPS, pair coefficients must be defined for ALL type pairs 1..Ntypes,
        # even if some types have zero atoms. We'll set a neutral Morse for everything, and then
        # override the specific attractive pairs. Core-core (1/12) uses soft repulsion.
        f.write(f"pair_coeff 1 1 soft {SOFT_A:.6f} {CORE_DIAMETER:.6f}\n")
        f.write(f"pair_coeff 1 12 soft {SOFT_A:.6f} {CORE_DIAMETER:.6f}\n")
        f.write(f"pair_coeff 12 12 soft {SOFT_A:.6f} {CORE_DIAMETER:.6f}\n")

        all_types = list(range(1, 13))  # 1..12
        for i, t1 in enumerate(all_types):
            for t2 in all_types[i:]:
                if (t1, t2) in ((1, 1), (1, 12), (12, 12)):
                    continue  # already set to soft
                f.write(f"pair_coeff {t1} {t2} morse 0.0 {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")

        # Channel ABC attractions (A/C to B faces)
        for bf in (3, 5):
            f.write(f"pair_coeff 2 {bf} morse {D0_BIND:.6f} {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
            f.write(f"pair_coeff 4 {bf} morse {D0_BIND:.6f} {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")

        # Channel EFD attractions (E/F to D faces)
        for df in (9, 11):
            f.write(f"pair_coeff 8 {df} morse {D0_BIND:.6f} {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
            f.write(f"pair_coeff 10 {df} morse {D0_BIND:.6f} {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
        f.write("\n")

        f.write("# Neighbor settings\n")
        f.write("neighbor        1.0 bin\n")
        f.write("neigh_modify    every 1 delay 0 check yes\n\n")

        f.write("# Rigid body integration + Langevin thermostat\n")
        f.write("group all_monomers molecule > 0\n")
        f.write("fix fx_nve all_monomers rigid/nve molecule\n")
        f.write("fix fx_langevin all_monomers langevin 0.5 0.5 0.5 12345\n\n")

        f.write("# State change (both channels)\n")
        f.write(f"fix sc patches state/change/dimer_ksat_twoside_twins 100 {CONTACT_CUTOFF:.6f} 1.0 patches 5\n\n")

        f.write("# Output\n")
        f.write("thermo_style    custom step temp pe ke etotal\n")
        f.write("thermo          1000\n")
        f.write("dump            d1 all custom 100 dump.dimer_ksat_1core_twosideB_twins.lammpstrj id mol type x y z\n")
        f.write("dump_modify     d1 sort id\n\n")

        f.write("timestep        0.005\n")
        f.write("run             2000000\n")

    print(f"Wrote {data_path} and {in_path}")


if __name__ == "__main__":
    main()


