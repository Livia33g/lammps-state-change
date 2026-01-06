#!/usr/bin/env python3
"""
Generate data/in files for the rigid dimer state-change simulation.

Specs:
- 10 Red monomers (patch type 2), 5 Blue monomers (patch type 3)
- 6 atoms per monomer (3 cores type 1, 3 patches)
- IMPORTANT: Red and Blue have the SAME geometry; only patch type differs.
- Masses: cores 0.2 each (total 0.6), patches 0.1 each
- Density ~0.001 -> box ~44.8 sigma side for 90 atoms
- Potentials:
    1-1 soft repulsion A=500, rmax=2.0, alpha=2.5 (pair_style soft)
    2-3 Morse D0=4.0, alpha=5.0, r0=0.0, rcut=1.6
    3-3 Morse D0=10.0, alpha=5.0, r0=0.0, rcut=1.6
    2-2 and all others: D0=0
- State change handled by fix_state_change_dimer (2->3 on 2-3 contact)
"""

import os
import math
import random
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
    out_dir = Path("dimer_simulation_cpp")
    out_dir.mkdir(parents=True, exist_ok=True)

    n_red = 10
    n_blue = 5
    n_mol = n_red + n_blue
    atoms_per_mol = 6
    natoms = n_mol * atoms_per_mol  # 90

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

    # Box size from total monomer concentration c_total ~ 0.001 (monomers / sigma^3)
    # i.e. V = N_monomers / c_total
    c_total = 0.001
    volume = n_mol / c_total  # sigma^3
    box_l = math.pow(volume, 1 / 3)  # ~44.8

    # Smaller initial spacing helps increase encounter rate without changing concentration
    positions = make_positions(n_mol, box_l, spacing=3.0)
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
        f.write("pair_style hybrid morse 1.6 soft 2.0\n")
        f.write("# 1-1 soft repulsion\n")
        f.write("pair_coeff 1 1 soft 500.0 2.0\n")
        f.write("# Neutral interactions\n")
        f.write("pair_coeff 1 2 morse 0.0 5.0 0.0\n")
        f.write("pair_coeff 1 3 morse 0.0 5.0 0.0\n")
        f.write("pair_coeff 2 2 morse 0.0 5.0 0.0\n")
        f.write("# Transition and sink\n")
        # Requirement: 3-3 is 2x stronger than 2-3
        f.write("pair_coeff 2 3 morse 6.0 5.0 0.0\n")
        f.write("pair_coeff 3 3 morse 12.0 5.0 0.0\n\n")
        f.write("# Neighbor settings\n")
        f.write("neighbor        1.0 bin\n")
        f.write("neigh_modify    every 1 delay 0 check yes\n\n")
        f.write("# Rigid body integration + Langevin thermostat (more stable rotations)\n")
        f.write("group all_monomers molecule > 0\n")
        f.write("fix fx_nve all_monomers rigid/nve molecule\n")
        # Stronger damping reduces excess spinning (smaller damp => stronger friction)
        f.write("fix fx_langevin all_monomers langevin 0.5 0.5 0.5 12345\n\n")
        f.write("# State change: check every 100 steps, pflip=1.0\n")
        # Hysteresis: require 5 consecutive checks in contact (500 steps) to ensure a real RB dimer first
        f.write("fix sc patches state/change/dimer 100 1.6 1.0 patches 5\n\n")
        f.write("# Output\n")
        f.write("thermo_style    custom step temp pe ke etotal\n")
        f.write("thermo          1000\n")
        # Dump every 100 steps (same as nevery) so you can see RB contact before flip
        f.write("dump            d1 all custom 100 dump.dimer.lammpstrj id mol type x y z\n")
        f.write("dump_modify     d1 sort id\n\n")
        f.write("# Run\n")
        f.write("timestep        0.005\n")
        # Longer run to allow encounters at low concentration
        f.write("run             2000000\n")

    print(f"Wrote {data_path} and {in_path}")


if __name__ == "__main__":
    main()




