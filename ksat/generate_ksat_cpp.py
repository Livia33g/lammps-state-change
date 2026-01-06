"""
Generate LAMMPS input script for ksat patchy particles with state changes
using the custom C++ fix_state_change_ksat.

Structure:
- Each monomer has 1 core (type 1) + 3 patches = 4 atoms total
- Patches on ONE side only (+x side, from part1tripatch_mid.mol)
- Monomer A: ALL patches are type 2 (never change)
- Monomer B: ALL patches start as type 3, can change to type 4 when attached to type 2
- Only attractive interactions: 2-3 and 2-4 (Morse)
"""

import numpy as np
import os
from datetime import datetime

def create_lammps_ksat_script_cpp(
    num_monomers_A=200,  # Number of type A monomers (patches type 2, never change)
    num_monomers_B=400,  # Number of type B monomers (patches type 3, can change to 4)
    box_size=None,
    patch_coordination_cutoff=0.28,  # Coordination cutoff for state change detection
    state_change_probability=0.7,
    timesteps=2700000,  # Match original: 900k at T=1.0 + 1800k at T=0.4
    thermo_freq=20000,
    dump_freq=10000,
    output_dir="ksat_simulation_cpp",
    seed=12345,
    state_change_freq=100,
    cooldown_steps=1000,
    timestep=0.0002,
    temperature_initial=1.0,
    temperature_final=0.4,
    steps_T1=900000,
    steps_T2=1800000,
):
    """
    Generates LAMMPS input script for ksat monomers using the custom C++ fix_state_change_ksat.
    """
    
    num_monomers = num_monomers_A + num_monomers_B
    
    # Calculate box size if not provided
    if box_size is None:
        # Use similar density to original (600 monomers in 20x20x20 box)
        volume = num_monomers * (20.0**3) / 600.0
        box_size = volume ** (1.0/3.0)
        # Add padding
        box_size *= 1.1
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    data_filename = os.path.join(output_dir, "data.ksat_monomers")
    lammps_input_filename = os.path.join(output_dir, "in.ksat_monomers")
    
    # Read molecule template (7-atom structure: 1 core + 6 patches)
    mol_template_path = os.path.join(os.path.dirname(__file__), "part1tripatch_mid.mol")
    if not os.path.exists(mol_template_path):
        raise FileNotFoundError(f"Molecule template not found: {mol_template_path}")
    
    # Parse molecule template to get relative coordinates
    with open(mol_template_path, 'r') as f:
        lines = f.readlines()
    
    # Extract coordinates from template
    coords_start = False
    types_start = False
    coords = []
    types = []
    for line in lines:
        if line.strip().startswith("Coords"):
            coords_start = True
            continue
        if line.strip().startswith("Types"):
            coords_start = False
            types_start = True
            continue
        if line.strip().startswith("Masses") or line.strip().startswith("#"):
            types_start = False
            continue
        if coords_start and line.strip() and not line.strip().startswith("#"):
            parts = line.split()
            if len(parts) >= 4:
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        if types_start and line.strip() and not line.strip().startswith("#"):
            parts = line.split()
            if len(parts) >= 2:
                types.append(int(parts[1]))
    
    # Separate body (type 1) and patches - use only ONE side (3 patches instead of 6)
    # Template has: type 1 = core, type 2 = patches on +x side, type 3 = patches on -x side
    # We'll use only the +x side patches (type 2 in template)
    body_coords = []
    patch_coords = []
    for i, (coord, atom_type) in enumerate(zip(coords, types)):
        if atom_type == 1:
            body_coords.append(coord)
        elif atom_type == 2:  # Only use +x side patches (type 2 in template)
            patch_coords.append(coord)
        # Skip type 3 patches (-x side) - we only want patches on one side
    
    if len(body_coords) != 1:
        raise ValueError(f"Expected 1 body particle, got {len(body_coords)}")
    if len(patch_coords) != 3:
        raise ValueError(f"Expected 3 patches (one side only), got {len(patch_coords)}")
    
    # Place monomers in box
    atom_data_lines = []
    current_atom_id = 1
    existing_centers = []
    min_center_distance = 2.0  # Minimum distance between monomer centers
    
    np.random.seed(seed)
    max_attempts = 10000
    
    # Place monomers with grid-based + random fallback
    grid_size = int(np.ceil(num_monomers**(1.0/3.0)))
    grid_positions = []
    for i in range(grid_size):
        for j in range(grid_size):
            for k in range(grid_size):
                grid_positions.append((i, j, k))
    
    np.random.shuffle(grid_positions)
    
    placed_count = 0
    spacing = box_size / (grid_size + 1)
    
    # Place first batch in grid
    for idx, (gi, gj, gk) in enumerate(grid_positions):
        if placed_count >= num_monomers:
            break
        com_x = spacing * (gi + 1) + np.random.uniform(-spacing*0.3, spacing*0.3)
        com_y = spacing * (gj + 1) + np.random.uniform(-spacing*0.3, spacing*0.3)
        com_z = spacing * (gk + 1) + np.random.uniform(-spacing*0.3, spacing*0.3)
        
        # Apply PBC
        com_x = com_x % box_size
        com_y = com_y % box_size
        com_z = com_z % box_size
        
        existing_centers.append((com_x, com_y, com_z))
        placed_count += 1
    
    # Random placement for remaining
    for mol_id in range(placed_count + 1, num_monomers + 1):
        placed = False
        for attempt in range(max_attempts):
            com_x = np.random.uniform(0, box_size)
            com_y = np.random.uniform(0, box_size)
            com_z = np.random.uniform(0, box_size)
            
            too_close = False
            for prev_center in existing_centers:
                dx = com_x - prev_center[0]
                dy = com_y - prev_center[1]
                dz = com_z - prev_center[2]
                dx -= box_size * np.round(dx / box_size)
                dy -= box_size * np.round(dy / box_size)
                dz -= box_size * np.round(dz / box_size)
                dist = np.sqrt(dx * dx + dy * dy + dz * dz)
                if dist < min_center_distance:
                    too_close = True
                    break
            
            if not too_close:
                existing_centers.append((com_x, com_y, com_z))
                placed = True
                break
        
        if not placed:
            raise RuntimeError(f"Failed to place monomer {mol_id} after {max_attempts} attempts")
    
    # Generate atoms
    for mol_id, (com_x, com_y, com_z) in enumerate(existing_centers, start=1):
        # Determine monomer type (A or B)
        if mol_id <= num_monomers_A:
            patch_type = 2  # Monomer A: patches are type 2 (never change)
        else:
            patch_type = 3  # Monomer B: patches start as type 3 (can change to 4)
        
        # Add body particle (type 1)
        body_coord = body_coords[0]
        x = com_x + body_coord[0]
        y = com_y + body_coord[1]
        z = com_z + body_coord[2]
        atom_data_lines.append(f"{current_atom_id} {mol_id} 1 {x} {y} {z}")
        current_atom_id += 1
        
        # Add ALL patches - all 3 patches (one side) have the same type on each monomer
        for patch_coord in patch_coords:
            x = com_x + patch_coord[0]
            y = com_y + patch_coord[1]
            z = com_z + patch_coord[2]
            atom_data_lines.append(f"{current_atom_id} {mol_id} {patch_type} {x} {y} {z}")
            current_atom_id += 1
    
    num_atoms = len(atom_data_lines)
    num_atom_types = 4  # Types: 1=body, 2=patches (monomer A), 3=patches (monomer B initial), 4=patches (monomer B final)
    
    # Write data file
    with open(data_filename, 'w') as f:
        f.write(f"LAMMPS data file for ksat patchy monomers\n\n")
        f.write(f"{num_atoms} atoms\n")
        f.write(f"{num_atom_types} atom types\n\n")
        f.write(f"0.0 {box_size} xlo xhi\n")
        f.write(f"0.0 {box_size} ylo yhi\n")
        f.write(f"0.0 {box_size} zlo zhi\n\n")
        
        f.write("Masses\n\n")
        f.write(f"1 1.0\n")      # Body (size=1.0, mass = size^3)
        f.write(f"2 1e-08\n")    # Patches - essentially massless (matches octahedron)
        f.write(f"3 1e-08\n")
        f.write(f"4 1e-08\n")
        f.write("\n")
        
        f.write("Atoms # molecular\n\n")
        f.write("\n".join(atom_data_lines))
        f.write("\n")
    
    print(f"Generated {data_filename}")
    print(f"Total atoms: {num_atoms}, Atoms per monomer: {num_atoms // num_monomers} (1 core + 3 patches on one side)")
    print(f"Monomer A: {num_monomers_A} (all 3 patches type 2), Monomer B: {num_monomers_B} (all 3 patches type 3 → 4)")
    
    # Write LAMMPS input script
    with open(lammps_input_filename, 'w') as f:
        f.write("# LAMMPS input script for ksat patchy particles with state change\n")
        f.write(f"# Generated: {datetime.now()}\n")
        f.write(f"# Monomers: {num_monomers_A} type A + {num_monomers_B} type B = {num_monomers} total\n")
        f.write(f"# State change: Type 3 → 4 when attached to type 2\n\n")
        
        # Setup
        f.write("units lj\n")
        f.write("dimension 3\n")
        f.write("atom_style molecular\n")
        f.write("boundary p p p\n\n")
        
        f.write("# Read data file\n")
        f.write(f"read_data {os.path.basename(data_filename)}\n\n")
        
        # Groups
        f.write("# Define groups\n")
        f.write("group body type 1\n")
        f.write("group patches type 2 3 4\n")
        f.write("group type2 type 2\n")
        f.write("group all type 1 2 3 4\n\n")
        
        # Soft relaxation
        f.write("# ===================== Soft relaxation =====================\n")
        f.write("pair_style soft 2.0\n")
        f.write("pair_coeff * * 300.0\n\n")
        
        f.write("timestep 0.0001\n")
        f.write("neighbor 1.0 bin\n")
        f.write("neigh_modify delay 0 every 1 check yes\n")
        f.write("comm_modify cutoff 4.5\n\n")
        
        f.write("fix R all rigid/nve/small molecule langevin 0.5 0.5 0.1 98765\n")
        f.write("run 4000\n")
        f.write("unfix R\n\n")
        
        # Physical interactions
        f.write("# ===================== Physical interactions: WCA + Morse =====================\n")
        f.write("pair_style hybrid/overlay lj/cut 2.5 morse 0.6\n")
        f.write("pair_modify shift yes\n\n")
        
        # WCA core-core
        f.write("# WCA core-core (sigma=1.0, rcut=1.1224620483)\n")
        f.write("pair_coeff 1 1 lj/cut 1.0 1.0 1.1224620483\n\n")
        
        # WCA core-patch
        f.write("# WCA core-patch (sigma=0.6, rcut=0.67347722899)\n")
        f.write("pair_coeff 1 2 lj/cut 1.0 0.6 0.67347722899\n")
        f.write("pair_coeff 1 3 lj/cut 1.0 0.6 0.67347722899\n")
        f.write("pair_coeff 1 4 lj/cut 1.0 0.6 0.67347722899\n\n")
        
        # WCA patch-patch
        f.write("# WCA patch-patch (sigma=0.20, rcut=0.22449240966)\n")
        f.write("variable rcPP equal 0.22449240966\n")
        for t1 in [2, 3, 4]:
            for t2 in [2, 3, 4]:
                if t2 >= t1:  # Avoid duplicates
                    f.write(f"pair_coeff {t1} {t2} lj/cut 1.0 0.20 ${{rcPP}}\n")
        f.write("\n")
        
        # Stronger steric clash for same-type patches (3-3 and 4-4)
        f.write("# Stronger steric clash for same-type patches\n")
        f.write("variable fac equal 1.122462048309373\n")
        f.write("variable sig33 equal 0.30\n")
        f.write("variable rc33 equal ${fac}*${sig33}\n")
        f.write("pair_coeff 3 3 lj/cut 1.0 ${sig33} ${rc33}\n")
        f.write("pair_coeff 4 4 lj/cut 1.0 ${sig33} ${rc33}\n\n")
        
        # Morse attractions (2-3 and 2-4 only)
        f.write("# Morse attractions: 2-3 and 2-4 only\n")
        f.write("pair_coeff 2 3 morse 8.0 20.0 0.225\n")
        f.write("pair_coeff 2 4 morse 8.0 20.0 0.225\n\n")
        
        # State change setup
        f.write("# ===================== State change setup =====================\n")
        f.write(f"fix STATE all state/change/ksat {state_change_freq} {cooldown_steps} {state_change_probability} {patch_coordination_cutoff} patches\n\n")
        
        # MD run
        f.write("# ===================== Rigid-body MD =====================\n")
        f.write(f"timestep {timestep}\n")
        f.write("neighbor 1.0 bin\n")
        f.write("neigh_modify delay 0 every 1 check yes\n")
        f.write("neigh_modify exclude molecule/intra all\n\n")
        
        f.write(f"thermo {thermo_freq}\n")
        f.write("thermo_style custom step temp press pe ke f_STATE\n")
        f.write("thermo_modify lost warn\n")
        f.write(f"dump mydmp all custom {dump_freq} dump.ksat_monomers.lammpstrj id mol type x y z\n\n")
        
        f.write("velocity all create 1.0 987654 dist gaussian rot yes\n\n")
        
        f.write(f"# Run at T={temperature_initial}\n")
        f.write(f"fix R all rigid/nve/small molecule langevin {temperature_initial} {temperature_initial} 0.1 1234\n")
        f.write(f"run {steps_T1}\n")
        f.write("unfix R\n\n")
        
        f.write(f"# Run at T={temperature_final}\n")
        f.write(f"fix R all rigid/nve/small molecule langevin {temperature_initial} {temperature_final} 0.1 4321\n")
        f.write(f"run {steps_T2}\n")
        f.write("unfix R\n\n")
        
        f.write("print 'Simulation completed successfully'\n")
    
    print(f"Generated {lammps_input_filename}")
    print("✅ Script generated using C++ fix - NO unfix/refix cycles needed!")
    print("   State changes happen automatically during the run.")

if __name__ == "__main__":
    create_lammps_ksat_script_cpp(
        num_monomers_A=200,
        num_monomers_B=400,
        box_size=None,
        patch_coordination_cutoff=0.28,
        state_change_probability=0.7,
    )

