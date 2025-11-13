import numpy as np
import os
from datetime import datetime # Corrected: added 'import' keyword

def create_lammps_rigid_patchy_monomers_script(
    num_monomers,
    box_size=None,  # If None, calculated from volume = N/0.001
    rigid_sphere_radius=1.0,
    patch_radius=0.33,  # Changed from 0.4
    patch_offset_dist=1.0,
    patch_interaction_cutoff_morse=2.5,  # Increased from 1.5 to 2.5 for longer-range attraction
    state_change_probability=0.5,
    timesteps=500000,
    thermo_freq=100,
    dump_freq=1000,
    output_dir="rigid_patchy_monomer_simulation",
    seed=12345
):
    """
    Generates LAMMPS data and input scripts for rigid patchy monomers with reversible state change.

    Args:
        num_monomers (int): Number of rigid monomer-patch pairs.
        box_size (float, optional): Size of the simulation box. 
            If None, calculated from volume formula: volume = N/0.001, box_size = volume^(1/3)
        rigid_sphere_radius (float): Radius of the 3 spheres forming the rigid body (Type 1).
        patch_radius (float): Radius of the patch (Type 2/3). Default 0.33.
        patch_offset_dist (float): Distance of the patch from the rigid body's center of mass.
        patch_interaction_cutoff_morse (float): Cutoff for Morse potential and state change detection.
        state_change_probability (float): Probability (0-1) for a state change.
        timesteps (int): Total simulation timesteps.
        thermo_freq (int): Frequency for thermodynamic output.
        dump_freq (int): Frequency for trajectory output.
        output_dir (str): Directory to save output files.
        seed (int): Random seed for reproducibility.
    """
    # Calculate box_size from volume formula: volume = N/0.001, box_size = volume^(1/3)
    if box_size is None:
        volume = num_monomers / 0.001
        box_size = volume ** (1.0/3.0)
        print(f"Calculated box_size from volume formula (N/0.001 = {volume}): {box_size:.4f}")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    data_filename = os.path.join(output_dir, "data.rigid_patchy_monomers")
    lammps_input_filename = os.path.join(output_dir, "in.rigid_patchy_monomers")

    # --- Atom Type Definitions ---
    # Type 1: Main Body Sphere (soft sphere potential)
    # Type 2: Patch A (Morse attractive to other Type 2, changes to Type 3)
    # Type 3: Patch B (Morse attractive to other Type 3, changes to Type 2)

    main_body_mass = 1.0  # Each main sphere has mass 1.0
    patch_mass = 1e-13  # Patches are essentially massless

    # --- Rigid Body Relative Positions ---
    # Centered at (0,0,0) for definition, will be translated later
    # The three spheres of radius 'rigid_sphere_radius'
    a = rigid_sphere_radius
    c6 = np.cos(np.pi/6) # cos(30 deg)
    s6 = np.sin(np.pi/6) # sin(30 deg)

    rigid_body_coords_relative = [
        [0.0, 0.0, a],
        [0.0, a * c6, -a * s6],
        [0.0, -a * c6, -a * s6]
    ]
    
    # Patch at the center of the three body spheres (which is at the COM)
    # The patch is part of the rigid body.
    patch_coord_relative = [0.75, 0.0, 0.0]

    # --- Generate Initial Atom Data ---
    atom_data_lines = []
    current_atom_id = 1
    
    np.random.seed(seed) # Set seed for reproducible initial placement

    # Calculate maximum extent of rigid body to ensure atoms stay in box
    max_x_extent = max([abs(c[0]) for c in rigid_body_coords_relative] + [abs(patch_coord_relative[0])]) + rigid_sphere_radius
    max_y_extent = max([abs(c[1]) for c in rigid_body_coords_relative] + [abs(patch_coord_relative[1])]) + rigid_sphere_radius
    max_z_extent = max([abs(c[2]) for c in rigid_body_coords_relative] + [abs(patch_coord_relative[2])]) + rigid_sphere_radius
    
    # Ensure center of mass placement leaves room for rigid body extent
    padding_x = max_x_extent + 0.1
    padding_y = max_y_extent + 0.1
    padding_z = max_z_extent + 0.1
    
    for mol_id in range(1, num_monomers + 1):
        # Random initial position for the center of mass, ensuring all atoms stay in box
        # Use stricter bounds to ensure all atoms are well inside the box
        com_x = np.random.uniform(padding_x, box_size - padding_x)
        com_y = np.random.uniform(padding_y, box_size - padding_y)
        com_z = np.random.uniform(padding_z, box_size - padding_z)

        # Calculate all atom positions first
        atom_positions = []
        for rel_coord in rigid_body_coords_relative:
            x = com_x + rel_coord[0]
            y = com_y + rel_coord[1]
            z = com_z + rel_coord[2]
            atom_positions.append((x, y, z))
        
        # Add patch position
        x_p = com_x + patch_coord_relative[0]
        y_p = com_y + patch_coord_relative[1]
        z_p = com_z + patch_coord_relative[2]
        atom_positions.append((x_p, y_p, z_p))
        
        # Check if any atom is outside the box and adjust COM if needed
        min_x = min([pos[0] for pos in atom_positions])
        max_x = max([pos[0] for pos in atom_positions])
        min_y = min([pos[1] for pos in atom_positions])
        max_y = max([pos[1] for pos in atom_positions])
        min_z = min([pos[2] for pos in atom_positions])
        max_z = max([pos[2] for pos in atom_positions])
        
        # If any atom is outside, shift COM to ensure all atoms are inside
        # Add small margin (0.01) to keep atoms well inside the box
        margin = 0.01
        if min_x < margin:
            com_x += (margin - min_x)
        if max_x > box_size - margin:
            com_x -= (max_x - (box_size - margin))
        if min_y < margin:
            com_y += (margin - min_y)
        if max_y > box_size - margin:
            com_y -= (max_y - (box_size - margin))
        if min_z < margin:
            com_z += (margin - min_z)
        if max_z > box_size - margin:
            com_z -= (max_z - (box_size - margin))
        
        # Recalculate all positions with adjusted COM (preserves rigid body structure)
        atom_positions = []
        for rel_coord in rigid_body_coords_relative:
            x = com_x + rel_coord[0]
            y = com_y + rel_coord[1]
            z = com_z + rel_coord[2]
            atom_positions.append((x, y, z))
        
        x_p = com_x + patch_coord_relative[0]
        y_p = com_y + patch_coord_relative[1]
        z_p = com_z + patch_coord_relative[2]
        atom_positions.append((x_p, y_p, z_p))
        
        # Write atoms (rigid body structure is preserved)
        for i, (x, y, z) in enumerate(atom_positions[:-1]):  # All but last (patch)
            atom_data_lines.append(f"{current_atom_id} {mol_id} 1 {x} {y} {z}")
            current_atom_id += 1
        
        # Write patch (last position) - assign types: 9 monomers with type 2, 1 with type 3
        x_p, y_p, z_p = atom_positions[-1]
        # First monomer gets type 3, all others get type 2
        if mol_id == 1:
            patch_type = 3  # First monomer has type 3 (patch B)
        else:
            patch_type = 2  # All other monomers have type 2 (patch A)
        atom_data_lines.append(f"{current_atom_id} {mol_id} {patch_type} {x_p} {y_p} {z_p}")
        current_atom_id += 1
        
    num_atoms = len(atom_data_lines)
    num_atom_types = 3 # Type 1 (body), Type 2 (patch A), Type 3 (patch B)
    
    # --- Write LAMMPS Data File ---
    with open(data_filename, 'w') as f:
        f.write(f"LAMMPS data file for rigid patchy monomers\n\n")
        f.write(f"{num_atoms} atoms\n")
        f.write(f"{num_atom_types} atom types\n\n")
        f.write(f"0.0 {box_size} xlo xhi\n")
        f.write(f"0.0 {box_size} ylo yhi\n")
        f.write(f"0.0 {box_size} zlo zhi\n\n")

        f.write("Masses\n\n")
        f.write(f"1 {main_body_mass}\n")
        f.write(f"2 {patch_mass}\n")
        f.write(f"3 {patch_mass}\n\n")

        f.write("Atoms # molecular\n\n")
        f.write("\n".join(atom_data_lines))
        f.write("\n")

    print(f"Generated {data_filename}")

    # Calculate half steps for splitting runs in the loop
    steps_half = thermo_freq // 2

    # --- Write LAMMPS Input Script (.in file) ---
    # Start of the LAMMPS input script string
    lammps_script_content = f"""# LAMMPS Input Script for Rigid Patchy Monomers with Reversible State Change
# Generated by Python script on {datetime.now().isoformat(timespec='seconds')}

# 1. Initialization
units           lj
atom_style      molecular
boundary        p p p
newton          on

# 2. Atom Definition
read_data       data.rigid_patchy_monomers

# 2.5. Expand box by 5% and remap atoms to ensure all atoms are inside
# Use remap BEFORE setting up rigid bodies to move atoms into box
# This is safe because rigid bodies haven't been created yet
change_box      all x scale 1.05 remap
change_box      all y scale 1.05 remap
change_box      all z scale 1.05 remap

# 3. Force Field Parameters
pair_style      hybrid/overlay soft 1.0 morse {patch_interaction_cutoff_morse}
pair_coeff      * * soft 1.0
pair_coeff      1 1 soft 0.3  # Reduced body-body repulsion for more attraction
pair_coeff      2 2 morse 10.0 2.0 0.66  # Strong attraction: D0=10.0 (well depth), r0=0.66 (equilibrium distance)
pair_coeff      3 3 morse 10.0 2.0 0.66  # Strong attraction: D0=10.0 (well depth), r0=0.66 (equilibrium distance)
pair_coeff      2 3 soft 1.0 # Soft repulsion between different patch types
pair_coeff      1 2 soft 1.0 # Soft repulsion between body and patch
pair_coeff      1 3 soft 1.0 # Soft repulsion between body and patch
# Set all morse pair coefficients (required for hybrid/overlay)
pair_coeff      1 1 morse 0.0 1.0 1.0 # No morse interaction for body-body
pair_coeff      1 2 morse 0.0 1.0 1.0 # No morse interaction for body-patch_A
pair_coeff      1 3 morse 0.0 1.0 1.0 # No morse interaction for body-patch_B
pair_coeff      2 3 morse 0.0 1.0 1.0 # No morse interaction for patch_A-patch_B
"""

    # Add rigid body definition - Use SINGLE fix for ALL molecules
    lammps_script_content += f"""
# 4. Rigid Body Definition (temporary, for minimization only)
# Create a temporary rigid fix for minimization (without thermostat)
fix             rigid_temp all rigid molecule

# 5. Simulation Settings
# Communication cutoff for rigid bodies (should be large enough for rigid body extent)
neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes page 200000 one 20000
comm_style      brick
comm_modify     vel yes cutoff 5.0

# 5.5. Minimize initial configuration to avoid overlaps
minimize        1.0e-4 1.0e-6 100 1000

# 5.6. Reset timestep after minimization
reset_timestep  0

# 5.7. Remove temporary rigid fix and create the NVT version
unfix           rigid_temp

# 5.8. Initialize velocities BEFORE creating the NVT fix
# This ensures the system has kinetic energy to start dynamics
# Use rot no to minimize rotational momentum (only translational)
velocity        all create 1.0 87287 mom yes rot no

# 6. Dynamics - use very small time step for stability with rigid bodies
timestep        0.001  # Increased timestep for faster dynamics

# 6.5. Create NVT fix AFTER velocity initialization
# Use fix rigid/nvt which combines rigid body constraints with NVT thermostat
fix             rigid_nvt all rigid/nvt molecule temp 1.0 1.0 0.01

# 6.6. Verify all atoms are in box before dynamics
variable        natoms_before equal count(all)
print          "DEBUG: Atoms before dynamics: ${{natoms_before}}"

# 7. State Change Logic - Define groups
group           patches_A type 2
group           patches_B type 3

# 8. Output - MUST be defined before loop
thermo          100
thermo_style    custom step temp pe ke etotal press
dump            1 all custom 1000 dump.rigid_patchy_monomers.lammpstrj id mol type x y z
dump_modify     1 sort id

# 9. Run single steps first to check stability
run             1
variable        natoms_after1 equal count(all)
print          "DEBUG: Atoms after step 1: ${{natoms_after1}}"

run             1
variable        natoms_after2 equal count(all)
print          "DEBUG: Atoms after step 2: ${{natoms_after2}}"

# If we get here, continue with more steps
run             8

# 10. NOW define computes for coordination numbers (after system is stable)
# Use 'all' so computes work for all atoms (returns 0.0 for atoms not in the group)
compute         cA all coord/atom cutoff {patch_interaction_cutoff_morse} group patches_A
compute         cB all coord/atom cutoff {patch_interaction_cutoff_morse} group patches_B

# Variable for random seed (will be incremented in loop)
variable        rand_seed equal {seed}
variable        rand_val atom random(0,1,${{rand_seed}})

# Define variables for type changes (will be redefined in loop)
# Strategy: Use a single variable that works for all atoms
# For patches_A (type 2): if coordinated and random<prob, change to 3, else stay 2
# For patches_B (type 3): if coordinated and random<prob, change to 2, else stay 3
# For body atoms (type 1): always stay 1
# Use the 'type' keyword to get current type, then modify based on conditions
# Type 1: type = 1, so 1 + 0 - 0 = 1 (no change)
# Type 2: type = 2, so 2 + (condition) - 0 = 2 or 3
# Type 3: type = 3, so 3 + 0 - (condition) = 2 or 3
variable        new_type atom "type + (type==2) * (c_cA>0.1) * (v_rand_val<{state_change_probability}) - (type==3) * (c_cB>0.1) * (v_rand_val<{state_change_probability})"

# 11. Main simulation loop with state changes
variable        total_steps equal {timesteps}
variable        steps_per_iteration equal {thermo_freq}
variable        iteration equal 0

label           loop

  # Check if we've reached total steps
  variable      current_step equal ${{iteration}}*${{steps_per_iteration}}
  if "${{current_step}} >= ${{total_steps}}" then "jump SELF end_loop"
  
  # Increment iteration counter
  variable      iteration equal ${{iteration}}+1
  
  # Delete old computes first (must be done before redefining groups they reference)
  uncompute       cA
  uncompute       cB
  
  # Update groups based on current types (redefining automatically replaces old definition)
  group           patches_A type 2
  group           patches_B type 3
  
  # Redefine computes for updated groups (use 'all' so computes work for all atoms)
  compute         cA all coord/atom cutoff {patch_interaction_cutoff_morse} group patches_A
  compute         cB all coord/atom cutoff {patch_interaction_cutoff_morse} group patches_B
  
  # Regenerate random values with new seed
  variable        rand_seed equal ${{rand_seed}}+1000
  variable        rand_val atom random(0,1,${{rand_seed}})
  
  # Redefine variables for type changes after computes are recreated
  # Use same formula as initial definition to ensure consistency
  variable        new_type atom "type + (type==2) * (c_cA>0.1) * (v_rand_val<{state_change_probability}) - (type==3) * (c_cB>0.1) * (v_rand_val<{state_change_probability})"
  
  # Run simulation steps - split into smaller chunks for stability
  run             {steps_half}
  run             {steps_half}
  
  # CRITICAL: Delete rigid fix before changing atom types
  unfix           rigid_nvt
  
  # Remap atoms into box before type changes (safety measure)
  change_box      all remap
  
  # Apply type changes using single variable for all atoms
  # The variable formula ensures:
  # - Type 1 atoms stay 1
  # - Type 2 atoms become 2 or 3 based on coordination and random
  # - Type 3 atoms become 2 or 3 based on coordination and random
  set             atom * type v_new_type
  
  # Remap atoms again after type changes
  change_box      all remap
  
  # CRITICAL: Recreate rigid fix after type change
  fix             rigid_nvt all rigid/nvt molecule temp 1.0 1.0 0.01
  
  # Reinitialize velocities after recreating rigid bodies (to maintain temperature)
  velocity        all scale 1.0
  
  # Reinitialize system after recreating rigid bodies
  run             0
  
jump            SELF loop

label           end_loop

print "Simulation finished. Total steps: ${{current_step}}"
"""

    # Write the complete script to the file in one go
    with open(lammps_input_filename, 'w') as f:
        f.write(lammps_script_content)

    print(f"Generated {lammps_input_filename}")

# Example usage (parameters can be adjusted as needed):
create_lammps_rigid_patchy_monomers_script(
    num_monomers=10,
    box_size=None,  # Auto-calculated: (10/0.001)^(1/3) ≈ 21.54
    rigid_sphere_radius=1.0,
    patch_radius=0.33,
    patch_offset_dist=1.0,
    patch_interaction_cutoff_morse=2.5,  # Increased from 1.5 to 2.5 for longer-range attraction
    state_change_probability=0.7,  # 70% probability of state change when coordinated
    timesteps=500000000,  # 500 million timesteps (1000x original)
    thermo_freq=5000,     # Thermo output every 5000 steps (fewer outputs for long sim)
    dump_freq=50000,      # Trajectory every 50000 steps (10000 frames for 500M steps)
    output_dir="rigid_patchy_simulation",
    seed=12345
)
