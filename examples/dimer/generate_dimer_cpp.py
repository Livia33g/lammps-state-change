#!/usr/bin/env python3
"""
Generate `data.dimer` + `in.dimer` for the rigid dimer state-change simulation.

This simulation has two monomer "colors":
- Red monomers: patch atom type 2
- Blue monomers: patch atom type 3

Both colors have identical geometry; only the patch atom type differs.

Common usage in this repo:
- kT = 0.5 (Langevin thermostat target)
- total monomer concentration c_total = 0.001 (monomers / sigma^3)
- 50/50 mixture Red/Blue (requires an even number of monomers)

Potentials used in the generated input:
- core-core (1-1): soft repulsion
- patch attractions: Morse for (2-3) and (3-3)
- everything else: Morse with D0=0 (neutral)
"""

import os
import math
import random
import argparse
from pathlib import Path


def make_positions(n_mol, box_l, spacing=5.0, seed=1234):
    random.seed(seed)
    coords = []
    per_side = math.ceil(n_mol ** (1 / 3))
    idx = 0
    # Place on a simple cubic lattice inside the box.
    # Choose spacing from box length so we don't place outside when box is small.
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
    parser = argparse.ArgumentParser(
        description="Generate LAMMPS data/input for dimer state-change simulation"
    )
    parser.add_argument(
        "--out_dir",
        type=Path,
        default=Path("dimer_simulation_cpp"),
        help="Output directory (default: dimer_simulation_cpp)",
    )
    parser.add_argument(
        "--c_total",
        type=float,
        default=0.001,
        help="Total monomer concentration in 1/sigma^3 (default: 0.001)",
    )
    parser.add_argument(
        "--n_total",
        type=int,
        default=16,
        help="Total number of monomers (default: 16; pick even for 50/50)",
    )
    parser.add_argument(
        "--red_fraction",
        type=float,
        default=0.5,
        help="Fraction of Red monomers (default: 0.5)",
    )
    parser.add_argument(
        "--kt",
        type=float,
        default=0.5,
        help="Thermostat target temperature kT (default: 0.5)",
    )
    parser.add_argument(
        "--morse_rcut",
        type=float,
        default=1.6,
        help="Morse cutoff for patch interactions (default: 1.6)",
    )
    parser.add_argument(
        "--morse_alpha",
        type=float,
        default=5.0,
        help="Morse alpha for patch interactions (default: 5.0)",
    )
    parser.add_argument(
        "--morse_r0",
        type=float,
        default=0.0,
        help="Morse r0 (location of minimum) for patch interactions (default: 0.0)",
    )
    parser.add_argument(
        "--D23",
        type=float,
        default=6.0,
        help="Morse well depth D0 for 2-3 (Red-Blue) attraction (default: 6.0)",
    )
    parser.add_argument(
        "--D33",
        type=float,
        default=12.0,
        help="Morse well depth D0 for 3-3 (Blue-Blue) attraction (default: 12.0)",
    )
    parser.add_argument(
        "--dump_every",
        type=int,
        default=1000,
        help="Dump cadence in timesteps (default: 1000)",
    )
    parser.add_argument(
        "--thermo_every",
        type=int,
        default=1000,
        help="Thermo print cadence in timesteps (default: 1000)",
    )
    parser.add_argument(
        "--run_steps",
        type=int,
        default=20000000,
        help="Number of MD timesteps to run (default: 20000000)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1234,
        help="Seed for initial placement (default: 1234)",
    )
    args = parser.parse_args()

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.n_total <= 0:
        raise SystemExit("--n_total must be > 0")
    if not (0.0 <= args.red_fraction <= 1.0):
        raise SystemExit("--red_fraction must be in [0, 1]")
    n_red = int(round(args.n_total * args.red_fraction))
    n_blue = args.n_total - n_red
    if args.red_fraction == 0.5 and (args.n_total % 2 != 0):
        raise SystemExit(
            "50/50 requested (red_fraction=0.5) but n_total is odd; pick an even --n_total."
        )
    n_mol = n_red + n_blue
    atoms_per_mol = 6
    natoms = n_mol * atoms_per_mol

    # geometry (relative to COM)
    cores = [
        (0.0, 0.0, 1.0),
        (0.0, 0.866, -0.5),
        (0.0, -0.866, -0.5),
    ]
    patches_red = [
        (1.0, 0.0, 0.3),
        (1.0, 0.260, -0.150),
        (1.0, -0.260, -0.150),
    ]
    # Same geometry for Blue as Red (only type differs)
    patches_blue = list(patches_red)

    # masses
    mass_core = 0.2
    mass_patch = 0.1

    # Box size from total monomer concentration (monomers / sigma^3):
    #   V = N_monomers / c_total
    volume = n_mol / args.c_total  # sigma^3
    box_l = math.pow(volume, 1 / 3)

    # Smaller initial spacing helps increase encounter rate without changing concentration
    positions = make_positions(n_mol, box_l, spacing=3.0, seed=args.seed)
    assert len(positions) == n_mol

    # Write data file
    data_path = out_dir / "data.dimer"
    with data_path.open("w") as f:
        f.write("LAMMPS data file for rigid dimer state-change\n\n")
        f.write(f"{natoms} atoms\n")
        f.write("0 bonds\n0 angles\n0 dihedrals\n0 impropers\n\n")
        f.write("3 atom types\n\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} xlo xhi\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} ylo yhi\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} zlo zhi\n\n")
        f.write("Masses\n\n")
        f.write(f"1 {mass_core:.6f}\n")
        f.write(f"2 {mass_patch:.6f}\n")
        f.write(f"3 {mass_patch:.6f}\n\n")
        f.write("Atoms # full\n\n")

        atom_id = 1
        mol_id = 1
        for idx, (cx, cy, cz) in enumerate(positions):
            is_red = idx < n_red
            patch_coords = patches_red if is_red else patches_blue
            patch_type = 2 if is_red else 3
            # cores
            for (dx, dy, dz) in cores:
                f.write(f"{atom_id} {mol_id} 1 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                atom_id += 1
            # patches
            for (dx, dy, dz) in patch_coords:
                f.write(f"{atom_id} {mol_id} {patch_type} 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                atom_id += 1
            mol_id += 1

    # Write input script
    in_path = out_dir / "in.dimer"
    with in_path.open("w") as f:
        f.write("# Dimer state-change NESS simulation\n")
        f.write("units           lj\n")
        f.write("atom_style      full\n")
        f.write("boundary        p p p\n\n")
        f.write("read_data       data.dimer\n\n")
        f.write("# Groups\n")
        f.write("group cores type 1\n")
        f.write("group patches type 2 3\n\n")
        f.write("# Pair potentials\n")
        # NOTE: for LAMMPS morse, r0 sets the location of the energy minimum.
        # We want patch "overlap"/contact to be favorable, so use r0 = 0.0.
        f.write(f"pair_style hybrid morse {args.morse_rcut:.6f} soft 2.0\n")
        f.write("# 1-1 soft repulsion\n")
        f.write("pair_coeff 1 1 soft 500.0 2.0\n")
        f.write("# Neutral interactions\n")
        f.write(f"pair_coeff 1 2 morse 0.0 {args.morse_alpha:.6f} {args.morse_r0:.6f}\n")
        f.write(f"pair_coeff 1 3 morse 0.0 {args.morse_alpha:.6f} {args.morse_r0:.6f}\n")
        f.write(f"pair_coeff 2 2 morse 0.0 {args.morse_alpha:.6f} {args.morse_r0:.6f}\n")
        f.write("# Transition and sink\n")
        # Requirement: 3-3 is 2x stronger than 2-3
        f.write(
            f"pair_coeff 2 3 morse {args.D23:.6f} {args.morse_alpha:.6f} {args.morse_r0:.6f}\n"
        )
        f.write(
            f"pair_coeff 3 3 morse {args.D33:.6f} {args.morse_alpha:.6f} {args.morse_r0:.6f}\n\n"
        )
        f.write("# Neighbor settings\n")
        f.write("neighbor        1.0 bin\n")
        f.write("neigh_modify    every 1 delay 0 check yes\n")
        # Rigid bodies should not have intra-molecule pair forces; they can create large
        # internal forces/torques and can be confusing in visualization.
        f.write("neigh_modify    exclude molecule/intra all\n\n")
        f.write("# Rigid body integration + Langevin thermostat (more stable rotations)\n")
        f.write("group all_monomers molecule > 0\n")
        f.write("fix fx_nve all_monomers rigid/nve molecule\n")
        # Stronger damping reduces excess spinning (smaller damp => stronger friction)
        f.write(
            f"fix fx_langevin all_monomers langevin {args.kt:.6f} {args.kt:.6f} 0.5 12345\n\n"
        )
        f.write("# State change: check every 100 steps, pflip=1.0\n")
        # Hysteresis: require 5 consecutive checks in contact (500 steps) to ensure a real RB dimer first
        f.write("fix sc patches state/change/dimer 100 1.6 1.0 patches 5\n\n")
        f.write("# Output\n")
        f.write("thermo_style    custom step temp pe ke etotal\n")
        f.write(f"thermo          {args.thermo_every}\n")
        # Dump every 100 steps (same as nevery) so you can see RB contact before flip
        f.write(
            f"dump            d1 all custom {args.dump_every} dump.dimer.lammpstrj id mol type x y z\n"
        )
        f.write("dump_modify     d1 sort id\n\n")
        f.write("# Run\n")
        f.write("timestep        0.005\n")
        # Longer run to allow encounters at low concentration
        f.write(f"run             {args.run_steps}\n")

    print(
        f"Wrote {data_path} and {in_path} | n_total={n_mol} (n_red={n_red}, n_blue={n_blue}) "
        f"| c_total={args.c_total} | box_l={box_l:.6f} | kT={args.kt}"
    )


if __name__ == "__main__":
    main()




