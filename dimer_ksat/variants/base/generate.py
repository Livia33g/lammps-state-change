#!/usr/bin/env python3
"""
Generate data/in files for the rigid "dimer_ksat" state-change simulation.

System:
- A monomers: patch type 2 (switchable)
- B monomers: patch type 3 (catalyst)
- C monomers: patch type 4 (non-switchable)

Desired physics:
- ONLY attractions are between B-A (2-3) and B-C (3-4)
- B-A and B-C attractions are IDENTICAL (same Morse parameters)
- A can switch to C only when A is in contact with B (handled by custom fix)
- C never switches

Rigid body:
- 6 atoms per monomer: 3 cores (type 1) + 3 patches (type 2/3/4)
- Same geometry for all monomers; only patch type differs

Potentials (matched to the working dimer_cpp setup):
- 1-1: soft repulsion A=500, rmax=2.0 (pair_style soft)
- 2-3 and 3-4: Morse attraction (tunable) D0=D0_BIND, alpha=5.0, r0=0.0, rcut=1.6
- all other morse pairs: D0=0

State change:
- `fix state/change/dimer_ksat` flips patch types 2->4 for an entire molecule
  when any type-2 patch is within cutoff of any type-3 patch.
"""

import math
import random
from pathlib import Path


def make_positions(n_mol, box_l, spacing=5.0, seed=1234):
    random.seed(seed)
    coords = []
    per_side = math.ceil(n_mol ** (1 / 3))
    idx = 0
    spacing = min(spacing, box_l / per_side * 0.9)
    start = -box_l / 2 + 0.5 * spacing

    for ix in range(per_side):
        for iy in range(per_side):
            for iz in range(per_side):
                if idx >= n_mol:
                    return coords
                coords.append(
                    (
                        start + spacing * ix,
                        start + spacing * iy,
                        start + spacing * iz,
                    )
                )
                idx += 1
    return coords


def main():
    out_dir = Path("simulation_cpp")
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- composition (edit freely) ---
    # Requested default: 30 monomers total, 20 A and 10 B, and no C initially.
    # C (type 4) appears only via A->C flips.
    n_A = 20
    n_B = 10
    n_C = 0
    n_mol = n_A + n_B + n_C

    # --- interaction strength tuning ---
    # B–A and B–C must match exactly; reduce D0 if dimers are too stable.
    D0_BIND = 1.0
    MORSE_ALPHA = 5.0
    MORSE_R0 = 0.0

    atoms_per_mol = 6
    natoms = n_mol * atoms_per_mol

    # geometry (relative to COM) -- same as working dimer
    cores = [
        (0.0, 0.0, 1.0),
        (0.0, 0.866, -0.5),
        (0.0, -0.866, -0.5),
    ]
    patches = [
        (1.0, 0.0, 0.3),
        (1.0, 0.260, -0.150),
        (1.0, -0.260, -0.150),
    ]

    # masses
    mass_core = 0.2
    mass_patch = 0.1

    # concentration -> box size (same as working dimer)
    c_total = 0.001
    volume = n_mol / c_total
    box_l = math.pow(volume, 1 / 3)

    positions = make_positions(n_mol, box_l, spacing=3.0)
    assert len(positions) == n_mol

    # --- write data file ---
    data_path = out_dir / "data.dimer_ksat"
    with data_path.open("w") as f:
        f.write("LAMMPS data file for rigid dimer_ksat state-change\n\n")
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

            # cores
            for (dx, dy, dz) in cores:
                f.write(
                    f"{atom_id} {mol_id} 1 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n"
                )
                atom_id += 1
            # patches
            for (dx, dy, dz) in patches:
                f.write(
                    f"{atom_id} {mol_id} {patch_type} 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n"
                )
                atom_id += 1

            mol_id += 1

    # --- write input script ---
    in_path = out_dir / "in.dimer_ksat"
    with in_path.open("w") as f:
        f.write("# dimer_ksat state-change simulation (A->C catalyzed by B)\n")
        f.write("units           lj\n")
        f.write("atom_style      full\n")
        f.write("boundary        p p p\n\n")
        f.write("read_data       data.dimer_ksat\n\n")

        f.write("# Groups\n")
        f.write("group cores type 1\n")
        f.write("group patches type 2 3 4\n")
        f.write("group patches_A type 2\n")
        f.write("group patches_B type 3\n")
        f.write("group patches_C type 4\n\n")

        f.write("# Pair potentials\n")
        f.write("pair_style hybrid morse 1.6 soft 2.0\n")
        f.write("# 1-1 soft repulsion (excluded volume)\n")
        f.write("pair_coeff 1 1 soft 500.0 2.0\n")

        # Neutral interactions (morse D0=0)
        f.write("pair_coeff 1 2 morse 0.0 5.0 0.0\n")
        f.write("pair_coeff 1 3 morse 0.0 5.0 0.0\n")
        f.write("pair_coeff 1 4 morse 0.0 5.0 0.0\n")
        f.write("pair_coeff 2 2 morse 0.0 5.0 0.0\n")
        f.write("pair_coeff 2 4 morse 0.0 5.0 0.0\n")
        f.write("pair_coeff 3 3 morse 0.0 5.0 0.0\n")
        f.write("pair_coeff 4 4 morse 0.0 5.0 0.0\n")

        # Only attractions: B-A and B-C, identical strength
        f.write(f"pair_coeff 2 3 morse {D0_BIND:.6f} {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n")
        f.write(f"pair_coeff 3 4 morse {D0_BIND:.6f} {MORSE_ALPHA:.6f} {MORSE_R0:.6f}\n\n")

        f.write("# Neighbor settings\n")
        f.write("neighbor        1.0 bin\n")
        f.write("neigh_modify    every 1 delay 0 check yes\n\n")

        f.write("# Rigid body integration + Langevin thermostat\n")
        f.write("group all_monomers molecule > 0\n")
        f.write("fix fx_nve all_monomers rigid/nve molecule\n")
        f.write("fix fx_langevin all_monomers langevin 0.5 0.5 0.5 12345\n\n")

        f.write("# State change: A (2) -> C (4) upon A-B contact\n")
        f.write("# Hysteresis: require 5 consecutive checks in contact (5 * 100 steps)\n")
        f.write("fix sc patches state/change/dimer_ksat 100 1.6 1.0 patches 5\n\n")

        f.write("# Output\n")
        f.write("thermo_style    custom step temp pe ke etotal\n")
        f.write("thermo          1000\n")
        f.write("dump            d1 all custom 100 dump.dimer_ksat.lammpstrj id mol type x y z\n")
        f.write("dump_modify     d1 sort id\n\n")

        f.write("# Run\n")
        f.write("timestep        0.005\n")
        f.write("run             2000000\n")

    print(f"Wrote {data_path} and {in_path}")


if __name__ == "__main__":
    main()


