"""
Generate LAMMPS input script for rigid octahedron monomers with state changes
using the custom C++ fix_state_change.

Octahedron structure:
- 6 vertices (body particles, type 1)
- Each vertex has 4 patches (alternating type 2 and 3)
- Total: 30 particles per monomer (6 vertices × 5 particles each)
"""

import numpy as np
import os
from datetime import datetime

def create_lammps_octahedron_script_cpp(
    num_monomers,
    box_size=None,
    vertex_radius=2.0,
    patch_radius=0.5,
    patch_coordination_cutoff=0.34,
    state_change_probability=0.7,
    timesteps=500000000,
    thermo_freq=5000,
    dump_freq=50000,
    output_dir="octahedron_simulation_cpp",
    seed=12345,
    state_change_freq=100,
    cooldown_steps=1000,
    timestep=0.001,
    morse_D0_22=None,  # Morse D0 for type 2-2 interactions
    morse_D0_33=None,  # Morse D0 for type 3-3 interactions
    morse_D0_23=None,  # Morse D0 for type 2-3 interactions (cross-interaction)
    rep_epsilon=None,  # Repulsion epsilon for body-body
    temperature=None,
    morse_cutoff=None  # Morse potential cutoff
):
    """
    Generates LAMMPS input script for octahedron monomers using the custom C++ fix_state_change.
    """
    
    # Calculate max extent first (needed for box size calculation)
    # We'll calculate this after building geometry, but set a reasonable default
    if box_size is None:
        # Use density formula similar to JAX code: box = (300 / 0.001) ** (1 / 3)
        density = 0.001
        volume = num_monomers / density
        box_size = volume ** (1.0/3.0)
    
    # Store original box_size for later validation
    original_box_size = box_size
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    data_filename = os.path.join(output_dir, "data.octahedron_monomers")
    lammps_input_filename = os.path.join(output_dir, "in.octahedron_monomers")

    # --- Atom Type Definitions ---
    # Type 1: Vertex centers (body particles)
    # Type 2: Patch type A (can change to type 3)
    # Type 3: Patch type B (can change to type 2)
    vertex_mass = 0.5  # From JAX code: [0.5, 1e-8, 1e-8, 1e-8, 1e-8]
    patch_mass = 1.0e-8  # Small but non-zero

    # --- Monomer Geometry ---
    # Each monomer = 1 vertex (body particle) + 4 patches = 5 particles total
    # The full octahedron is formed by assembling 6 such monomers
    # 
    # From JAX code structure:
    # - Each vertex has 1 body particle at the vertex center
    # - 4 patches positioned around the vertex, pointing toward nearest neighbors
    # - The patches are positioned at vertex_radius (2.0) inward from the vertex
    
    # For a single monomer (vertex), the body particle is at the origin (COM)
    # Patches are positioned relative to the body/vertex
    
    # Get patch positions for a vertex monomer
    # Patches point toward the 4 nearest neighbors (other octahedron vertices)
    # For an octahedron, each vertex has 4 neighbors
    
    # Octahedron vertex positions (for reference - to calculate patch directions)
    scale = (2.0 / np.sqrt(2.0)) * vertex_radius
    octahedron_vertex_vectors = np.array([
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ])
    octahedron_vertex_positions = scale * octahedron_vertex_vectors
    
    # For a generic monomer, we'll position patches to point toward 4 neighbor directions
    # Use the first vertex as reference: its neighbors are vertices 2, 3, 4, 5 (indices 1, 2, 3, 4)
    reference_vertex_pos = octahedron_vertex_positions[0]  # [1, 0, 0] scaled
    neighbor_positions = octahedron_vertex_positions[1:5]  # 4 nearest neighbors
    
    # Build rigid body structure for one monomer
    # Body particle is at origin (COM)
    rigid_body_coords_relative = [[0.0, 0.0, 0.0]]  # Single body particle at COM
    
    # 4 patches positioned around the body, pointing toward neighbors
    patch_coords_relative = []
    for neighbor_pos in neighbor_positions:
        # Vector from body to neighbor
        vec_to_neighbor = neighbor_pos - reference_vertex_pos
        vec_norm = np.linalg.norm(vec_to_neighbor)
        if vec_norm > 1e-6:
            vec_to_neighbor_unit = vec_to_neighbor / vec_norm
            # Patch position is vertex_radius inward from body (toward neighbors)
            # Patches should point outward from the body, but be positioned slightly inward
            patch_pos = -vertex_radius * vec_to_neighbor_unit  # Negative because patches point toward center
            patch_coords_relative.append(patch_pos.tolist())
    
    # --- Generate Initial Atom Data ---
    atom_data_lines = []
    current_atom_id = 1
    
    np.random.seed(seed)

    # Calculate max extent for placement
    all_coords = rigid_body_coords_relative + patch_coords_relative
    max_extent = max([np.linalg.norm(c) for c in all_coords]) + vertex_radius
    padding = max_extent + 0.5
    
    # Minimum distance between monomer centers
    # Use a smaller multiplier - 1.5x instead of 2.0x for better packing
    min_center_distance = 1.5 * max_extent + 0.5
    
    # Calculate minimum box size needed for num_monomers
    # For N monomers, we need space for N monomers with min_center_distance spacing
    # Rough estimate: box should be large enough for a 3D grid of monomers
    # Using a simple scaling: each monomer needs ~min_center_distance^3 of volume
    min_volume_per_monomer = min_center_distance ** 3
    min_total_volume = num_monomers * min_volume_per_monomer * 1.5  # 1.5x for safety
    min_box_size_from_volume = (min_total_volume) ** (1.0/3.0)
    
    # Also ensure minimum for at least 2 monomers at edges
    min_box_size_for_placement = 2 * padding + min_center_distance
    
    # Use the larger of the two minimums
    min_box_size = max(min_box_size_from_volume, min_box_size_for_placement)
    
    if box_size < min_box_size:
        print(f"Warning: Calculated box_size {box_size:.4f} too small for {num_monomers} monomers.")
        print(f"  Minimum required: {min_box_size:.4f}")
        print(f"  Adjusting box_size to minimum required size.")
        box_size = min_box_size
    
    existing_centers = []
    max_attempts = 10000  # Increased attempts
    
    # Try a more systematic approach: use a grid-based initial placement, then random attempts
    # Calculate grid dimensions
    grid_size = int(np.ceil(num_monomers ** (1.0/3.0))) + 1
    if grid_size < 2:
        grid_size = 2
    grid_spacing = (box_size - 2 * padding) / max(1, grid_size - 1)
    
    # First try grid-based placement
    grid_positions = []
    for i in range(grid_size):
        for j in range(grid_size):
            for k in range(grid_size):
                com_x = padding + i * grid_spacing
                com_y = padding + j * grid_spacing
                com_z = padding + k * grid_spacing
                grid_positions.append((com_x, com_y, com_z))
    
    # Shuffle grid positions for randomness
    np.random.shuffle(grid_positions)
    
    placed_count = 0
    
    # Try placing monomers on grid points first
    for grid_pos in grid_positions:
        if placed_count >= num_monomers:
            break
            
        com_x, com_y, com_z = grid_pos
        
        # Check distance to previously placed monomers (with PBC)
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
            placed_count += 1
    
    # Fill remaining monomers with random placement
    for mol_id in range(placed_count + 1, num_monomers + 1):
        placed = False
        for attempt in range(max_attempts):
            com_x = np.random.uniform(padding, box_size - padding)
            com_y = np.random.uniform(padding, box_size - padding)
            com_z = np.random.uniform(padding, box_size - padding)

            # Check distance to previously placed monomers (with PBC)
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

            if too_close:
                continue

            existing_centers.append((com_x, com_y, com_z))
            placed = True
            break

        if not placed:
            raise RuntimeError(
                f"Failed to place monomer {mol_id} without overlap after {max_attempts} attempts. "
                f"Placed {len(existing_centers)}/{num_monomers} monomers. "
                f"Box size: {box_size:.4f}, min_distance: {min_center_distance:.4f}. "
                f"Increase box size or decrease num_monomers."
            )

    # Now generate atoms for all placed monomers
    for mol_id, (com_x, com_y, com_z) in enumerate(existing_centers, start=1):
        # Add body particle (single vertex, type 1)
        # Should only be one body particle per monomer
        if len(rigid_body_coords_relative) != 1:
            raise RuntimeError(f"Expected 1 body particle per monomer, got {len(rigid_body_coords_relative)}")
        
        rel_coord = rigid_body_coords_relative[0]
        x = com_x + rel_coord[0]
        y = com_y + rel_coord[1]
        z = com_z + rel_coord[2]
        atom_data_lines.append(f"{current_atom_id} {mol_id} 1 {x} {y} {z}")
        current_atom_id += 1

        # Add patches - ALL start as type 2 (user's "patch type 1", will change to 3/4/5 during simulation)
        # Note: We use type 2 in LAMMPS because type 1 is used for vertices (body particles)
        # The fix will treat type 2 as the initial patch type ("patch type 1")
        for i in range(len(patch_coords_relative)):
            patch_coord = patch_coords_relative[i]
            x = com_x + patch_coord[0]
            y = com_y + patch_coord[1]
            z = com_z + patch_coord[2]
            
            # All patches start as type 2 (initial "patch type 1", state changes will change to 3, 4, or 5)
            patch_type = 2
            
            atom_data_lines.append(f"{current_atom_id} {mol_id} {patch_type} {x} {y} {z}")
            current_atom_id += 1

    num_atoms = len(atom_data_lines)
    num_atom_types = 5  # Types: 1=vertices (body), 2=patches (initial "patch type 1"), 3,4,5=patches (evolved)
    
    # --- Write LAMMPS Data File ---
    with open(data_filename, 'w') as f:
        f.write(f"LAMMPS data file for rigid octahedron monomers\n\n")
        f.write(f"{num_atoms} atoms\n")
        f.write(f"{num_atom_types} atom types\n\n")
        f.write(f"0.0 {box_size} xlo xhi\n")
        f.write(f"0.0 {box_size} ylo yhi\n")
        f.write(f"0.0 {box_size} zlo zhi\n\n")

        f.write("Masses\n\n")
        f.write(f"1 {vertex_mass}\n")  # Type 1: vertices (body particles)
        f.write(f"2 {patch_mass}\n")   # Type 2: patches (initial "patch type 1")
        f.write(f"3 {patch_mass}\n")   # Type 3: patches (evolved state)
        f.write(f"4 {patch_mass}\n")   # Type 4: patches (evolved state)
        f.write(f"5 {patch_mass}\n")   # Type 5: patches (evolved state)
        f.write("\n")

        f.write("Atoms # molecular\n\n")
        f.write("\n".join(atom_data_lines))
        f.write("\n")

    print(f"Generated {data_filename}")
    print(f"Total atoms: {num_atoms}, Atoms per monomer: {num_atoms // num_monomers}")

    # --- Parameters for potentials ---
    # All patch-patch interactions are attractive with similar epsilon values
    morse_rcut = morse_cutoff if morse_cutoff is not None else 4.0
    morse_alpha = 2.0  # From JAX code
    morse_r0 = 0.0
    
    # Base epsilon for patch-patch attractions (all similar, varying by ~0.5-1.0)
    base_epsilon = 10.0
    morse_D0_11 = base_epsilon  # Type 1-1
    morse_D0_13 = base_epsilon + 0.5  # Type 1-3
    morse_D0_14 = base_epsilon + 0.3  # Type 1-4
    morse_D0_15 = base_epsilon + 0.7  # Type 1-5
    morse_D0_33 = base_epsilon + 1.0  # Type 3-3
    morse_D0_34 = base_epsilon + 0.6  # Type 3-4
    morse_D0_35 = base_epsilon + 0.2  # Type 3-5
    morse_D0_44 = base_epsilon + 0.8  # Type 4-4
    morse_D0_45 = base_epsilon + 0.4  # Type 4-5
    morse_D0_55 = base_epsilon + 0.9  # Type 5-5
    
    # Body-body repulsion
    rep_epsilon = rep_epsilon if rep_epsilon is not None else 10000.0
    rep_sigma = 2.0 * vertex_radius  # Contact distance
    rep_rmax_body = rep_sigma + 1.0
    
    # Temperature
    temp_value = temperature if temperature is not None else 1.0
    init_temp_value = 0.3 * temp_value

    # --- Write LAMMPS Input Script ---
    lammps_script_content = f"""# LAMMPS Input Script for Rigid Octahedron Monomers with State Changes
# Generated by Python script on {datetime.now().isoformat(timespec='seconds')}
# Uses custom C++ fix_state_change (NO unfix/refix required!)
# Octahedron monomer: 1 body (vertex) + 4 patches = 5 particles per monomer

# 1. Initialization
units           lj
atom_style      molecular
boundary        p p p
newton          on

# 2. Atom Definition
read_data       data.octahedron_monomers

# 2.5. Expand box by 5% and remap atoms
change_box      all x scale 1.05 remap
change_box      all y scale 1.05 remap
change_box      all z scale 1.05 remap

# 3. Force Field Parameters
# Using hybrid/overlay: morse for patch-patch, lj/cut for body-body repulsion
pair_style      hybrid/overlay morse {morse_rcut} lj/cut {rep_rmax_body}

# Initialize all interactions to zero
pair_coeff      * * morse 0.0 1.0 1.0 1.0
pair_coeff      * * lj/cut 0.0 1.0 1.0

# Patch-patch Morse attractions (all pairs can attract)
# Type 2 (initial patches, user's "patch type 1") interactions
pair_coeff      2 2 morse {morse_D0_11} {morse_alpha} {morse_r0} {morse_rcut}
pair_coeff      2 3 morse {morse_D0_13} {morse_alpha} {morse_r0} {morse_rcut}
pair_coeff      2 4 morse {morse_D0_14} {morse_alpha} {morse_r0} {morse_rcut}
pair_coeff      2 5 morse {morse_D0_15} {morse_alpha} {morse_r0} {morse_rcut}
# Type 3 interactions
pair_coeff      3 3 morse {morse_D0_33} {morse_alpha} {morse_r0} {morse_rcut}
pair_coeff      3 4 morse {morse_D0_34} {morse_alpha} {morse_r0} {morse_rcut}
pair_coeff      3 5 morse {morse_D0_35} {morse_alpha} {morse_r0} {morse_rcut}
# Type 4 interactions
pair_coeff      4 4 morse {morse_D0_44} {morse_alpha} {morse_r0} {morse_rcut}
pair_coeff      4 5 morse {morse_D0_45} {morse_alpha} {morse_r0} {morse_rcut}
# Type 5 interactions
pair_coeff      5 5 morse {morse_D0_55} {morse_alpha} {morse_r0} {morse_rcut}

# Body-body LJ repulsion (only vertex-vertex, type 1)
pair_coeff      1 1 lj/cut {rep_epsilon} {rep_sigma:.4f} {rep_rmax_body}

# Body-patch repulsion (weaker, between vertices and patches)
pair_coeff      1 2 lj/cut {rep_epsilon * 0.1} {rep_sigma * 0.8:.4f} {rep_rmax_body * 0.8}
pair_coeff      1 3 lj/cut {rep_epsilon * 0.1} {rep_sigma * 0.8:.4f} {rep_rmax_body * 0.8}
pair_coeff      1 4 lj/cut {rep_epsilon * 0.1} {rep_sigma * 0.8:.4f} {rep_rmax_body * 0.8}
pair_coeff      1 5 lj/cut {rep_epsilon * 0.1} {rep_sigma * 0.8:.4f} {rep_rmax_body * 0.8}

# 4. Rigid Body Definition (temporary, for minimization)
fix             rigid_temp all rigid molecule

# 5. Simulation Settings
neighbor        {max(morse_rcut, rep_rmax_body) + 0.5:.2f} bin
neigh_modify    delay 0 every 1 check yes page 200000 one 20000
comm_style      brick
comm_modify     vel yes cutoff {max(morse_rcut, rep_rmax_body) + 0.5:.2f}

# 5.5. Minimize
minimize        1.0e-4 1.0e-6 100 1000

# 5.6. Reset timestep
reset_timestep  0

# 5.7. Remove temporary rigid fix
unfix           rigid_temp

# 5.8. Initialize velocities
velocity        all create {init_temp_value} {seed} mom yes rot no

# 6. Dynamics
timestep        {timestep:.6f}

# 6.5. Rigid-body integrator and thermostat
fix             rigid_nvt all rigid/nvt molecule temp {temp_value} {temp_value} 200.0

# 7. Define group for all patches (types 2, 3, 4, 5)
group           patches type 2 3 4 5

# 8. Add state change fix (C++ fix - NO unfix/refix needed!)
# Syntax: fix ID group state/change/octahedron check_every cooldown_steps probability cutoff group_patches
# Note: Type 2 = initial "patch type 1", Types 3,4,5 = evolved states
fix             state_change all state/change/octahedron {state_change_freq} {cooldown_steps} {state_change_probability} {patch_coordination_cutoff} patches

# 10. Output
thermo          {thermo_freq}
thermo_style    custom step temp pe ke etotal press f_state_change
dump            1 all custom {dump_freq} dump.octahedron_monomers.lammpstrj id mol type x y z
dump_modify     1 sort id

# 11. Run simulation
# State changes happen automatically during the run - no loop needed!
variable        final_step equal step
run             {timesteps}
variable        final_step equal step

print "Simulation finished. Total steps: ${{final_step}}"
"""

    with open(lammps_input_filename, 'w') as f:
        f.write(lammps_script_content)

    print(f"Generated {lammps_input_filename}")
    print("\n✅ Script generated using C++ fix - NO unfix/refix cycles needed!")
    print("   State changes happen automatically during the run.")
    print(f"   Octahedron monomers: {num_monomers} monomers × 5 particles = {num_atoms} total atoms")

if __name__ == "__main__":
    create_lammps_octahedron_script_cpp(
        num_monomers=50,
        box_size=None,
        vertex_radius=2.0,
        patch_radius=0.5,
        patch_coordination_cutoff=0.34,
        state_change_probability=0.7,
        timesteps=500000000,
        thermo_freq=5000,
        dump_freq=10000,
        output_dir="octahedron_simulation_cpp",
        seed=12345,
        state_change_freq=100,
        cooldown_steps=1000,
        timestep=0.001,
        morse_D0_22=10.0,
        morse_D0_33=10.0,
        morse_D0_23=0.0,
        rep_epsilon=10000.0,
        temperature=1.0,
        morse_cutoff=4.0
    )

