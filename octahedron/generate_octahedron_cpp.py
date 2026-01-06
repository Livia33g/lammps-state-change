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
    vertex_radius=1.0,  # Body particle radius (matches reference geometry, though now hardcoded to 1.0)
    patch_radius=0.5,  # Not used - patches have fixed positions from reference
    patch_coordination_cutoff=0.34,
    state_change_probability=0.7,
    timesteps=500000000,
    thermo_freq=5000,
    dump_freq=50000,
    output_dir="octahedron_simulation_cpp",
    seed=12345,
    state_change_freq=100,
    cooldown_steps=1000,
    hysteresis_threshold=5000,  # New: Number of consecutive steps of contact required before state change
    timestep=0.0001,  # Further reduced from 0.0002 for stability (prevents NaNs from large forces)
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
        # INCREASED DENSITY for better cluster formation (collision frequency ∝ [concentration]²)
        # Previous: 0.0005 (too sparse, particles spend too much time in void)
        # New: 0.005 (10x increase for better collision rate, still safe for stability)
        density = 0.005  # Increased from 0.0005 for better cluster formation
        volume = num_monomers / density
        box_size = volume ** (1.0/3.0)
    
    # REDUCED box multiplier for higher effective concentration
    # Previous: 3.5x (too large, particles too spread out)
    # New: 2.0x (higher collision frequency for cluster formation)
    box_size = box_size * 2.0
    
    # Store original box_size for later validation
    original_box_size = box_size
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    data_filename = os.path.join(output_dir, "data.octahedron_monomers")
    lammps_input_filename = os.path.join(output_dir, "in.octahedron_monomers")

    # --- Atom Type Definitions ---
    # Type 1: Vertex centers (body particles)
    # Types 2–5: Patch types (can change state via custom fix)
    #
    # CRITICAL PHYSICS: Mass ≠ Size in LAMMPS
    # - Mass (m): Determines inertia and rotational stability (F=ma, I=Σm_i*r_i^2)
    # - Diameter/Potential (σ): Determines collision/overlap behavior
    # - Patches interact ONLY via Morse potential (soft, allows overlap regardless of mass)
    # - There is NO patch-patch LJ repulsion, so patches can overlap freely
    #
    # Mass Configuration for Rotational Stability:
    # - Body mass: 0.6 (provides base rotational inertia)
    # - Patch mass: 0.1 each (4 patches × 0.1 = 0.4 total)
    # - Total mass per monomer: 0.6 + 0.4 = 1.0
    # - This creates high moment of inertia (I = Σ m_i * r_i^2) for stable rotation
    # - Patches act as "flywheel" weights at distance ~2.0 from COM
    #
    # Overlap Behavior:
    # - Patches are point particles with NO finite-size repulsion
    # - Only Morse attraction between different patch types
    # - Heavy patches can still overlap because mass doesn't prevent it - only repulsion does
    # - Attracting patches will "sink into" each other creating deep energy wells
    vertex_mass = 0.6        # Body particle mass
    patch_mass = 0.1         # Patch mass (4 patches × 0.1 = 0.4 total, for rotational stability)

    # --- Monomer Geometry ---
    # Each monomer = 1 body particle + 4 patches = 5 particles total
    # Geometry matches reference files in octahedron_shape/:
    # - Body particle: at origin, mass = 1.0
    # - 4 patches: at distance 2.0 from center, mass = 1e-08
    # Interaction radii: body = 1.5, patch = 0.33 (for repulsion interactions)
    
    # Build rigid body structure for one monomer
    # Body particle is at origin (COM)
    rigid_body_coords_relative = [[0.0, 0.0, 0.0]]  # Single body particle at COM
    
    # 4 patches positioned at distance 2.0 from center
    # Patch positions from reference geometry (vertex_shape_points.npy):
    patch_coords_relative = [
        [-1.41421356,  1.39285561,  0.24485356],  # Patch 1: distance = 2.0
        [-1.41421356, -1.39285560, -0.24485356],  # Patch 2: distance = 2.0
        [-1.41421356,  0.24485356, -1.39285561],  # Patch 3: distance = 2.0
        [-1.41421356, -0.24485356,  1.39285560],  # Patch 4: distance = 2.0
    ]
    
    # --- Generate Initial Atom Data ---
    atom_data_lines = []
    current_atom_id = 1
    
    np.random.seed(seed)

    # Calculate max extent for placement
    # Patches are at distance 2.0 from center
    # Use body interaction radius = 1.5 for placement calculations
    all_coords = rigid_body_coords_relative + patch_coords_relative
    max_extent = max([np.linalg.norm(c) for c in all_coords]) + 1.5  # Use body radius for extent
    padding = max_extent + 0.5
    
    # Minimum distance between monomer centers
    # LARGER SPACING: Make monomers more dispersed to reduce bond breaking
    min_center_distance = 3.0 * max_extent + 1.0  # Increased from 1.5x to 3.0x for much more dispersion
    
    # Calculate minimum box size needed for num_monomers
    # For N monomers, we need space for N monomers with min_center_distance spacing
    # Rough estimate: box should be large enough for a 3D grid of monomers
    # Using a simple scaling: each monomer needs ~min_center_distance^3 of volume
    min_volume_per_monomer = min_center_distance ** 3
    min_total_volume = num_monomers * min_volume_per_monomer * 3.0  # 3.0x for much larger box (more dispersed)
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
        # For atom_style molecular: format is id mol type x y z
        # Diameter will be set via 'set' command in input script
        atom_data_lines.append(f"{current_atom_id} {mol_id} 1 {x} {y} {z}")
        current_atom_id += 1

        # Add patches - ALL start as type 2 (all patches on monomer are same type)
        # Note: We use type 2 in LAMMPS because type 1 is used for vertices (body particles)
        # The fix will treat type 2 as the initial patch type ("patch type 1")
        # All 4 patches on a monomer start as type 2 and change together as a unit
        for i in range(len(patch_coords_relative)):
            patch_coord = patch_coords_relative[i]
            x = com_x + patch_coord[0]
            y = com_y + patch_coord[1]
            z = com_z + patch_coord[2]
            
            # All patches start as type 2 (all patches on monomer start same, state changes change to 3, 4, or 5 together)
            patch_type = 2
            # For atom_style molecular: format is id mol type x y z
            # Diameter will be set via 'set' command in input script
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
        f.write(f"1 {vertex_mass}\n")  # Type 1: body (central particle)
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
    # STABILITY FIXES: Increased cutoff for better capture, decreased alpha for wider well
    # morse_rcut: Increased from 1.2 to 2.5 - acts like a funnel, guiding patches in from further away
    # morse_alpha: Decreased from 2.0 to 1.2 - widens the attractive basin, making bonds more forgiving
    morse_rcut = morse_cutoff if morse_cutoff is not None else 2.5  # Increased from 1.2 for better capture radius
    morse_alpha = 1.2  # Decreased from 2.0 to widen the well (was 2.0 from JAX code)
    morse_r0 = 0.0  # Allows deep overlap for stable bonding
    
    # Base epsilon for patch-patch attractions (ONLY for different types)
    # Same-type patches (3-3, 4-4, 5-5) have NO interaction (no attraction, no repulsion)
    # BOOSTED from 3.0 to 5.0 to increase bond strength and sticking probability
    base_epsilon = 5.0
    morse_D0_11 = base_epsilon  # Type 1-1 (initial state) - can attract to itself
    morse_D0_13 = base_epsilon + 0.5  # Type 1-3 (different types - attraction)
    morse_D0_14 = base_epsilon + 0.3  # Type 1-4 (different types - attraction)
    morse_D0_15 = base_epsilon + 0.7  # Type 1-5 (different types - attraction)
    # Same-type patches have NO interaction
    morse_D0_33 = 0.0  # Type 3-3 - NO interaction (no attraction, no repulsion)
    morse_D0_34 = base_epsilon + 0.6  # Type 3-4 (different types - attraction)
    morse_D0_35 = base_epsilon + 0.2  # Type 3-5 (different types - attraction)
    morse_D0_44 = 0.0  # Type 4-4 - NO interaction (no attraction, no repulsion)
    morse_D0_45 = base_epsilon + 0.4  # Type 4-5 (different types - attraction)
    morse_D0_55 = 0.0  # Type 5-5 - NO interaction (no attraction, no repulsion)
    
    # Body-body repulsion - use softer parameters to match JAX soft_sphere behavior
    # Body interaction radius = 1.5 (for repulsion sigma), patch radius = 0.5 (user: 0.33 or 0.5 both work)
    body_radius = 1.5  # Interaction radius for body particles
    patch_radius = 0.5  # Interaction radius for patch particles (0.33 or 0.5 both work)
    # Body-body repulsion strength - MATCH JAX MD soft_eps=10000.0 style
    # JAX uses soft_eps=10000.0 to prevent unwanted large clusterization
    # We use 5000.0 as a balance (stronger than 800.0 but not as extreme as 10000)
    # Strong repulsion ensures bodies bounce apart, preventing large unwanted clusters
    rep_epsilon = rep_epsilon if rep_epsilon is not None else 5000.0  # Increased from 800.0 to match JAX-style strong repulsion
    rep_sigma = 2.0 * body_radius  # Contact distance = 3.0 (2 * 1.5)
    # VERY SHORT RANGE: Almost only when touching
    rep_rmax_body = 1.1 * rep_sigma  # Very short range (3.3 instead of 4.5) - almost only when touching
    
    # Body-patch repulsion parameters (much weaker, short range)
    # Body radius = 1.5, patch radius = 0.33
    body_patch_sigma = body_radius + patch_radius  # Contact distance ~1.83 (body radius + patch radius)
    body_patch_rmax = 1.1 * body_patch_sigma  # Very short range (~2.0 instead of ~2.75) - almost only when touching
    body_patch_epsilon = rep_epsilon * 0.05  # Reduced from 0.1 to 0.05
    
    # Temperature (REDUCED for stability - working sims use 0.3)
    temp_value = temperature if temperature is not None else 0.3  # Reduced from 1.0 to 0.3
    init_temp_value = temp_value  # Use same temp as initial (not 0.3*temp)

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

# 2.1. Moment of inertia for rigid bodies is calculated automatically by LAMMPS
# from atom positions and masses relative to center of mass: I = Σ m_i * r_i^2
# No need to set diameters - atom_style molecular doesn't support it anyway
# Heavy rotational damping (angmom 100.0) and zero initial rotation (rot no) 
# will minimize rotation regardless of moment of inertia

# 2.5. Expand box by 20% and remap atoms (larger expansion for stability, prevents NaNs)
change_box      all x scale 1.20 remap
change_box      all y scale 1.20 remap
change_box      all z scale 1.20 remap

# 3. Force Field Parameters
# Using hybrid/overlay: morse for patch-patch, lj/cut (WCA-style) for body-body repulsion
pair_style      hybrid/overlay morse {morse_rcut} lj/cut {rep_rmax_body}

# Initialize all interactions to zero
pair_coeff      * * morse 0.0 1.0 1.0 1.0
pair_coeff      * * lj/cut 0.0 1.0 1.0

# Patch-patch interactions: ONLY Morse potential (NO LJ repulsion)
# CRITICAL: Patches have NO finite-size repulsion, so they can overlap freely
# Mass doesn't prevent overlap - only repulsion does. Heavy patches (mass=0.1) can still overlap.
# Morse potential is soft and allows deep overlap for stable bonding.
# Type 2 (initial patches, user's "patch type 1") interactions
pair_coeff      2 2 morse {morse_D0_11} {morse_alpha} {morse_r0} {morse_rcut}
pair_coeff      2 3 morse {morse_D0_13} {morse_alpha} {morse_r0} {morse_rcut}  # Different types - attraction
pair_coeff      2 4 morse {morse_D0_14} {morse_alpha} {morse_r0} {morse_rcut}  # Different types - attraction
pair_coeff      2 5 morse {morse_D0_15} {morse_alpha} {morse_r0} {morse_rcut}  # Different types - attraction
# Type 3 interactions
# Same-type (3-3) has NO interaction (no attraction, no repulsion)
pair_coeff      3 3 morse 0.0 {morse_alpha} {morse_r0} {morse_rcut}  # No interaction
pair_coeff      3 4 morse {morse_D0_34} {morse_alpha} {morse_r0} {morse_rcut}  # Different types - attraction
pair_coeff      3 5 morse {morse_D0_35} {morse_alpha} {morse_r0} {morse_rcut}  # Different types - attraction
# Type 4 interactions
# Same-type (4-4) has NO interaction (no attraction, no repulsion)
pair_coeff      4 4 morse 0.0 {morse_alpha} {morse_r0} {morse_rcut}  # No interaction
pair_coeff      4 5 morse {morse_D0_45} {morse_alpha} {morse_r0} {morse_rcut}  # Different types - attraction
# Type 5 interactions
# Same-type (5-5) has NO interaction (no attraction, no repulsion)
pair_coeff      5 5 morse 0.0 {morse_alpha} {morse_r0} {morse_rcut}  # No interaction
# NOTE: There is NO patch-patch LJ/cut interaction, so patches have zero repulsion and can overlap

# Body-body WCA repulsion (LJ cut at minimum, prevents NaNs)
# WCA: repulsive only, cut at r = 2^(1/6) * sigma
pair_coeff      1 1 lj/cut {rep_epsilon} {rep_sigma:.4f} {rep_rmax_body:.4f}

# Body-patch WCA repulsion (weaker, smaller sigma)
pair_coeff      1 2 lj/cut {body_patch_epsilon} {body_patch_sigma:.4f} {body_patch_rmax:.4f}
pair_coeff      1 3 lj/cut {body_patch_epsilon} {body_patch_sigma:.4f} {body_patch_rmax:.4f}
pair_coeff      1 4 lj/cut {body_patch_epsilon} {body_patch_sigma:.4f} {body_patch_rmax:.4f}
pair_coeff      1 5 lj/cut {body_patch_epsilon} {body_patch_sigma:.4f} {body_patch_rmax:.4f}

# 4. Simulation Settings (NO MINIMIZATION - LAMMPS doesn't support minimization with rigid bodies)
# Skipping minimization prevents patches from collapsing during energy minimization
# The initial configuration from data file is already correct
neighbor        {max(morse_rcut, rep_rmax_body) + 0.5:.2f} bin
neigh_modify    delay 0 every 1 check yes page 200000 one 20000
comm_style      brick
comm_modify     vel yes cutoff {max(morse_rcut, rep_rmax_body) + 0.5:.2f}

# 5. Initialize velocities at VERY LOW temperature for gentle start (prevents NaNs from overlaps)
# Start at 0.01 temperature to gently relax any overlaps before ramping up
# rot no = NO rotational velocity - monomers move translationally only (like point particles)
velocity        all create 0.01 {seed} mom yes rot no

# 5.5. Use smaller timestep initially for stability (prevents NaNs from large forces)
timestep        0.000050

# 6. Rigid-body integrator and thermostat - START AT VERY LOW TEMPERATURE
# Use rigid/nvt/small so both translation and rotation are thermostatted.
# Mass configuration: Body (0.6) + Patches (0.1 each) = 1.0 total per monomer
# This provides rotational stability through high moment of inertia (I = Σ m_i * r_i^2)
# Patches can still overlap because they interact only via Morse (soft, no repulsion)
# Tdamp=1.0 provides strong "shock absorber" to instantly absorb energy spikes from state changes
# This prevents spinning when strong Morse attractions (D0≈3.5) create sudden energy drops
# Start at 0.01 temperature to gently relax overlaps, then ramp up.
fix             rigid_nvt_equil all rigid/nvt/small molecule temp 0.01 0.01 1.0

# 6.5. Gentle equilibration run (10000 steps at low T to relax overlaps)
# This prevents NaNs from initial overlaps without using minimization
run             10000

# 6.6. Ramp up temperature gradually (prevents sudden force spikes)
# MUST unfix the equilibration fix before adding the ramp fix
# Tdamp=1.0 for strong energy absorption during state changes
unfix           rigid_nvt_equil
fix             rigid_nvt_ramp all rigid/nvt/small molecule temp 0.01 {temp_value} 1.0
run             20000

# 6.7. Now run at target temperature
# Tdamp=1.0 provides strong "shock absorber" to instantly absorb energy spikes from state changes
# When patches change state (e.g., Type 2→3), strong Morse attraction (D0≈3.5) creates huge energy drop
# This converts to kinetic energy → angular velocity (spin). Strong thermostat (Tdamp=1.0) immediately
# applies drag to remove excess KE, preventing runaway rotation. Patches bond smoothly without spinning.
unfix           rigid_nvt_ramp
fix             rigid_nvt all rigid/nvt/small molecule temp {temp_value} {temp_value} 1.0

# 6.8. Increase timestep to normal value after equilibration
timestep        {timestep:.6f}

# 7. Define group for all patches (types 2, 3, 4, 5)
group           patches type 2 3 4 5

# 8. Add state change fix (C++ fix - NO unfix/refix needed!)
# Syntax: fix ID group state/change/octahedron check_every cooldown_steps probability cutoff group_patches [hysteresis_threshold]
# Note: Type 2 = initial "patch type 1", Types 3,4,5 = evolved states
# Hysteresis: Requires consecutive contact for hysteresis_threshold steps before triggering state change (filters noise)
fix             state_change all state/change/octahedron {state_change_freq} {cooldown_steps} {state_change_probability} {patch_coordination_cutoff} patches {hysteresis_threshold}

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
        timestep=0.0002,  # Further reduced for stability
        morse_D0_22=10.0,
        morse_D0_33=10.0,
        morse_D0_23=0.0,
        rep_epsilon=10000.0,
        temperature=1.0,
        morse_cutoff=4.0
    )

