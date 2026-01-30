#!/usr/bin/env python3
"""
Generate LAMMPS simulation files from problem.json + params.json.

This is the "system generator" that creates the actual LAMMPS data and input
files needed to run a simulation for a given problem with user-specified parameters.

Usage:
    python generate_system_from_problem.py \
        --problem problems/problem-001-dimer-ksat/problem.json \
        --params submissions/problem-001/user/params.json \
        --policy submissions/problem-001/user/policy.json \
        --output sim_output/

Generates:
    sim_output/data.{system_name}    # LAMMPS data file
    sim_output/in.{system_name}      # LAMMPS input script
"""

import json
import sys
import math
import random
from pathlib import Path
from typing import Dict, Any, List, Tuple
import argparse


def load_json(path: str) -> Dict[str, Any]:
    """Load JSON file."""
    with open(path) as f:
        return json.load(f)


def merge_params(problem: Dict[str, Any], params: Dict[str, Any] = None) -> Dict[str, Any]:
    """Merge user params with defaults from problem definition."""
    # Extract default values from problem.json
    defaults = {}
    for param in problem['constraints'].get('tunable_parameters', []):
        defaults[param['name']] = param['default']

    # Override with user-provided params
    if params:
        for key, value in params.items():
            if key in defaults:
                defaults[key] = value
            elif key.startswith('morse_') or key.startswith('contact_'):
                # Allow params not explicitly in tunable list (for flexibility)
                defaults[key] = value

    return defaults


def make_positions(n_mol: int, box_l: float, seed: int = 1234) -> List[Tuple[float, float, float]]:
    """
    Place molecules on a cubic lattice spanning the whole box.
    Uses checkerboard sublattice to increase initial spacing.

    Based on examples/dimer_ksat/variants/1core/generate.py
    """
    rng = random.Random(seed)
    per_side = math.ceil(n_mol ** (1 / 3))
    spacing = box_l / per_side
    start = -box_l / 2 + 0.5 * spacing

    all_sites = [(ix, iy, iz)
                 for ix in range(per_side)
                 for iy in range(per_side)
                 for iz in range(per_side)]

    checker_sites = [s for s in all_sites if ((s[0] + s[1] + s[2]) % 2 == 0)]

    sites = checker_sites if len(checker_sites) >= n_mol else all_sites
    rng.shuffle(sites)
    sites = sites[:n_mol]

    return [(start + spacing * ix, start + spacing * iy, start + spacing * iz)
            for (ix, iy, iz) in sites]


def generate_data_file(problem: Dict[str, Any], params: Dict[str, Any],
                       output_path: Path) -> None:
    """Generate LAMMPS data file."""

    encoding = problem['encoding']
    system_params = encoding['system_parameters']

    # Get species info
    species = {s['label']: s for s in encoding['species']}
    composition = encoding['initial_composition']

    # Count molecules and atoms
    n_mol = sum(composition.values())

    # Determine atoms per molecule from geometry
    if encoding['geometry'] == '1core_3patch':
        atoms_per_mol = 4  # 1 core + 3 patches
        core_pos = [(0.0, 0.0, 0.0)]
        patch_pos = [
            (0.5, 0.0, 0.1),
            (0.5, 0.0866025404, -0.05),
            (0.5, -0.0866025404, -0.05),
        ]
        patch_neg = None
        is_twoside = False
    elif encoding['geometry'] == '1core_twosideB_twins':
        # Two independent channels: ABC (core type 1) and EFD (core type 12)
        # Simple molecules (A,C,E,F): 1 core + 3 patches on +x face
        # Two-face molecules (B,D): 1 core + 6 patches (3 on +x, 3 on -x)
        core_pos = [(0.0, 0.0, 0.0)]
        patch_pos = [
            (0.5, 0.0, 0.1),
            (0.5, 0.0866025404, -0.05),
            (0.5, -0.0866025404, -0.05),
        ]
        patch_neg = [(-x, y, -z) for (x, y, z) in patch_pos]  # 180° rotation about y
        is_twoside = True
        # Calculate total atoms (varies by species)
        atoms_simple = 1 + 3  # A,C,E,F
        atoms_twoface = 1 + 6  # B,D
        n_simple = sum(composition.get(label, 0) for label in ['A', 'C', 'E', 'F'])
        n_twoface = sum(composition.get(label, 0) for label in ['B', 'D'])
        natoms = n_simple * atoms_simple + n_twoface * atoms_twoface
    else:
        raise NotImplementedError(f"Geometry {encoding['geometry']} not yet supported")

    if not is_twoside:
        natoms = n_mol * atoms_per_mol

    # Box size
    c_total = system_params['density']
    volume = n_mol / c_total
    box_l_nominal = math.pow(volume, 1 / 3)
    box_shrink = system_params.get('box_shrink_factor', 1.0)
    box_l = box_l_nominal * box_shrink

    # Generate positions
    positions = make_positions(n_mol, box_l, seed=1234)

    # Determine species order (A, B, C, ...)
    species_order = list(composition.keys())

    # Write data file
    with output_path.open('w') as f:
        f.write(f"LAMMPS data file for {problem['problem_id']}\n\n")
        f.write(f"{natoms} atoms\n")
        f.write("0 bonds\n0 angles\n0 dihedrals\n0 impropers\n\n")

        # Count atom types (for twins: 12 types total: 1,2,3,4,5,8,9,10,11,12)
        if encoding['geometry'] == '1core_twosideB_twins':
            ntypes = 12
        else:
            ntypes = len(species)
        f.write(f"{ntypes} atom types\n\n")

        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} xlo xhi\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} ylo yhi\n")
        f.write(f"{-box_l/2:.4f} {box_l/2:.4f} zlo zhi\n\n")

        f.write("Masses\n\n")
        if encoding['geometry'] == '1core_twosideB_twins':
            # Write masses for all types 1-12
            mass_core = species.get('core_ABC', {}).get('mass', 0.6)
            mass_patch = species.get('A', {}).get('mass', 0.1)
            f.write(f"1 {mass_core:.6f}\n")  # ABC core
            for t in range(2, 6):  # Types 2,3,4,5
                f.write(f"{t} {mass_patch:.6f}\n")
            for t in range(8, 12):  # Types 8,9,10,11
                f.write(f"{t} {mass_patch:.6f}\n")
            f.write(f"12 {mass_core:.6f}\n")  # EFD core
        else:
            for label, spec in species.items():
                f.write(f"{spec['lammps_type']} {spec['mass']:.6f}\n")
        f.write("\n")

        f.write("Atoms # full\n\n")

        atom_id = 1
        mol_id = 1

        if encoding['geometry'] == '1core_twosideB_twins':
            # Create list of molecule kinds (A, B, C, D, E, F)
            kinds = []
            for label in ['A', 'B', 'C', 'D', 'E', 'F']:
                kinds.extend([label] * composition.get(label, 0))
            random.Random(4321).shuffle(kinds)
            
            # Place molecules
            for (cx, cy, cz), kind in zip(positions, kinds):
                # Determine core type: ABC channel uses 1, EFD channel uses 12
                core_type = 12 if kind in ("E", "F", "D") else 1
                f.write(f"{atom_id} {mol_id} {core_type} 0.0 {cx:.6f} {cy:.6f} {cz:.6f}\n")
                atom_id += 1
                
                if kind == "A":
                    patch_type = 2
                    for (dx, dy, dz) in patch_pos:
                        f.write(f"{atom_id} {mol_id} {patch_type} 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                        atom_id += 1
                elif kind == "C":
                    patch_type = 4
                    for (dx, dy, dz) in patch_pos:
                        f.write(f"{atom_id} {mol_id} {patch_type} 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                        atom_id += 1
                elif kind == "E":
                    patch_type = 8
                    for (dx, dy, dz) in patch_pos:
                        f.write(f"{atom_id} {mol_id} {patch_type} 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                        atom_id += 1
                elif kind == "F":
                    patch_type = 10
                    for (dx, dy, dz) in patch_pos:
                        f.write(f"{atom_id} {mol_id} {patch_type} 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                        atom_id += 1
                elif kind == "B":
                    # Two faces: +x face (type 3) and -x face (type 5)
                    for (dx, dy, dz) in patch_pos:
                        f.write(f"{atom_id} {mol_id} 3 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                        atom_id += 1
                    for (dx, dy, dz) in patch_neg:
                        f.write(f"{atom_id} {mol_id} 5 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                        atom_id += 1
                else:  # D
                    # Two faces: +x face (type 9) and -x face (type 11)
                    for (dx, dy, dz) in patch_pos:
                        f.write(f"{atom_id} {mol_id} 9 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                        atom_id += 1
                    for (dx, dy, dz) in patch_neg:
                        f.write(f"{atom_id} {mol_id} 11 0.0 {cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                        atom_id += 1
                
                mol_id += 1
        else:
            # Place molecules according to composition (original logic)
            for idx, (cx, cy, cz) in enumerate(positions):
                # Determine which species this molecule is
                species_label = None
                cumulative = 0
                for label in species_order:
                    cumulative += composition[label]
                    if idx < cumulative:
                        species_label = label
                        break

                if species_label is None:
                    species_label = species_order[-1]

                spec = species[species_label]
                patch_type = spec['lammps_type']

                # Write core atom (type 1)
                for (dx, dy, dz) in core_pos:
                    core_type = species['core']['lammps_type']
                    f.write(f"{atom_id} {mol_id} {core_type} 0.0 "
                           f"{cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1

                # Write patch atoms (type depends on species)
                for (dx, dy, dz) in patch_pos:
                    f.write(f"{atom_id} {mol_id} {patch_type} 0.0 "
                           f"{cx+dx:.6f} {cy+dy:.6f} {cz+dz:.6f}\n")
                    atom_id += 1

                mol_id += 1

    print(f"✓ Generated {output_path}")


def generate_input_script(problem: Dict[str, Any], params: Dict[str, Any],
                          policy: Dict[str, Any], output_path: Path,
                          data_file: str) -> None:
    """Generate LAMMPS input script."""

    encoding = problem['encoding']
    system_params = encoding['system_parameters']
    interactions = encoding['interactions']
    constraints = problem['constraints']

    # Get fix name from policy
    fix_name = policy['policy_name'].lower().replace(' ', '_').replace('-', '_')
    check_freq = policy['check_frequency']

    # Get contact cutoff (from policy or params)
    contact_cutoff = params.get('contact_cutoff', 0.7)

    # Get Morse parameters
    morse_depth_AB = params.get('morse_depth_AB', 1.0)
    morse_depth_BC = params.get('morse_depth_BC', 1.0)
    morse_alpha = params.get('morse_alpha', 5.0)
    morse_r0 = params.get('morse_r0', 0.0)

    # Calculate Morse cutoff
    patch_radius = system_params['patch_radius']
    morse_rcut = 7.0 * patch_radius

    # Core parameters
    core_diameter = system_params['core_diameter']
    soft_a = 500.0
    temperature = system_params['temperature']
    timestep = system_params['timestep']
    max_steps = constraints['max_timesteps']

    with output_path.open('w') as f:
        f.write(f"# {problem['title']}\n")
        f.write(f"# Generated from {problem['problem_id']}\n")
        f.write(f"# Policy: {policy['policy_name']}\n\n")

        f.write("units           lj\n")
        f.write("atom_style      full\n")
        f.write("boundary        p p p\n\n")

        f.write(f"read_data       {data_file}\n\n")

        f.write("# Groups\n")
        if encoding['geometry'] == '1core_twosideB_twins':
            f.write("group cores type 1 12\n")
            f.write("group patches type 2 3 4 5 8 9 10 11\n")
        else:
            f.write("group cores type 1\n")
            # Create groups for each species
            species = {s['label']: s for s in encoding['species']}
            patch_types = [s['lammps_type'] for s in species.values() if s['role'] == 'patch']
            f.write(f"group patches type {' '.join(map(str, patch_types))}\n")

            for label, spec in species.items():
                if spec['role'] == 'patch':
                    f.write(f"group patches_{label} type {spec['lammps_type']}\n")
        f.write("\n")

        f.write("# Pair potentials\n")
        f.write(f"pair_style hybrid morse {morse_rcut:.6f} soft {core_diameter:.6f}\n")

        if encoding['geometry'] == '1core_twosideB_twins':
            # Core-core soft repulsion (1-1, 1-12, 12-12)
            f.write(f"pair_coeff 1 1 soft {soft_a:.6f} {core_diameter:.6f}\n")
            f.write(f"pair_coeff 1 12 soft {soft_a:.6f} {core_diameter:.6f}\n")
            f.write(f"pair_coeff 12 12 soft {soft_a:.6f} {core_diameter:.6f}\n")
            
            # Set neutral Morse for all other pairs first
            all_types = list(range(1, 13))  # 1..12
            for i, t1 in enumerate(all_types):
                for t2 in all_types[i:]:
                    if (t1, t2) in ((1, 1), (1, 12), (12, 12)):
                        continue  # Already set to soft
                    f.write(f"pair_coeff {t1} {t2} morse 0.0 {morse_alpha:.6f} {morse_r0:.6f}\n")
            
            # Channel ABC attractions (A/C to B faces: types 2,4 with 3,5)
            morse_depth = params.get('morse_depth', 1.0)
            for bf in (3, 5):  # B_face1 and B_face2
                f.write(f"pair_coeff 2 {bf} morse {morse_depth:.6f} {morse_alpha:.6f} {morse_r0:.6f}\n")
                f.write(f"pair_coeff 4 {bf} morse {morse_depth:.6f} {morse_alpha:.6f} {morse_r0:.6f}\n")
            
            # Channel EFD attractions (E/F to D faces: types 8,10 with 9,11)
            for df in (9, 11):  # D_face1 and D_face2
                f.write(f"pair_coeff 8 {df} morse {morse_depth:.6f} {morse_alpha:.6f} {morse_r0:.6f}\n")
                f.write(f"pair_coeff 10 {df} morse {morse_depth:.6f} {morse_alpha:.6f} {morse_r0:.6f}\n")
        else:
            # Core-core soft repulsion
            f.write(f"pair_coeff 1 1 soft {soft_a:.6f} {core_diameter:.6f}\n")

            # Neutral interactions (morse D0=0)
            ntypes = len(species)
            for i in range(1, ntypes + 1):
                for j in range(i, ntypes + 1):
                    if i == 1 and j == 1:
                        continue  # Already defined core-core
                    f.write(f"pair_coeff {i} {j} morse 0.0 {morse_alpha:.6f} {morse_r0:.6f}\n")

            # Attractive interactions (from problem definition)
            for pair_info in interactions.get('attractive_pairs', []):
                types = pair_info['types'].split('-')
                type1, type2 = int(types[0]), int(types[1])

                # Determine which depth to use
                if '2' in types and '3' in types:  # A-B
                    depth = morse_depth_AB
                elif '3' in types and '4' in types:  # B-C
                    depth = morse_depth_BC
                else:
                    depth = 1.0  # default

                f.write(f"pair_coeff {type1} {type2} morse {depth:.6f} {morse_alpha:.6f} {morse_r0:.6f}\n")
        f.write("\n")

        f.write("# Neighbor settings\n")
        f.write("neighbor        1.0 bin\n")
        f.write("neigh_modify    every 1 delay 0 check yes\n\n")

        f.write("# Rigid body integration + Langevin thermostat\n")
        f.write("group all_monomers molecule > 0\n")
        f.write("fix fx_nve all_monomers rigid/nve molecule\n")
        f.write(f"fix fx_langevin all_monomers langevin {temperature:.3f} {temperature:.3f} 0.5 12345\n\n")

        f.write("# State change fix (auto-generated from policy)\n")
        f.write(f"# Policy: {policy['policy_name']}\n")
        f.write(f"# Check frequency: every {check_freq} steps\n")
        f.write(f"fix sc patches state/change/{fix_name} {check_freq} {contact_cutoff:.6f} 1.0 patches\n\n")

        f.write("# Output\n")
        f.write("thermo_style    custom step temp pe ke etotal\n")
        f.write("thermo          1000\n")

        # Dump file name from problem
        dump_file = f"dump.{problem['problem_id'].replace('problem-', '').replace('-', '_')}.lammpstrj"
        f.write(f"dump            d1 all custom 100 {dump_file} id mol type x y z\n")
        f.write("dump_modify     d1 sort id\n\n")

        f.write("# Run\n")
        f.write(f"timestep        {timestep:.6f}\n")
        f.write(f"run             {max_steps}\n")

    print(f"✓ Generated {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate LAMMPS simulation files from problem + policy + params'
    )
    parser.add_argument('--problem', required=True, help='Path to problem.json')
    parser.add_argument('--policy', required=True, help='Path to policy.json')
    parser.add_argument('--params', help='Path to params.json (optional)')
    parser.add_argument('--output', required=True, help='Output directory')

    args = parser.parse_args()

    # Load files
    problem = load_json(args.problem)
    policy = load_json(args.policy)
    params = load_json(args.params) if args.params else {}

    # Merge params with defaults
    merged_params = merge_params(problem, params)

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate system name
    system_name = problem['problem_id'].replace('problem-', '').replace('-', '_')

    # Generate data file
    data_file = f"data.{system_name}"
    data_path = output_dir / data_file
    generate_data_file(problem, merged_params, data_path)

    # Generate input script
    input_file = f"in.{system_name}"
    input_path = output_dir / input_file
    generate_input_script(problem, merged_params, policy, input_path, data_file)

    print(f"\n✓ System generation complete!")
    print(f"  Data file: {data_path}")
    print(f"  Input script: {input_path}")
    print(f"\nNext steps:")
    print(f"  1. Generate C++ fix from policy:")
    print(f"     python core/generators/generate_fix_from_policy.py {args.policy} {output_dir}/")
    print(f"  2. Compile LAMMPS with the fix")
    print(f"  3. Run: lmp_mpi -in {input_file}")


if __name__ == "__main__":
    main()
