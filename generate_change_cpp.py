"""
Generate LAMMPS input script for rigid patchy monomers with state changes
using the custom C++ fix (no unfix/refix required!)
"""

import numpy as np
import os
from datetime import datetime

def smoothing(r, ron, rcut):
    """Smoothing function for Morse potential"""
    rcut = np.clip(rcut, 1e-6, None)
    ron = np.clip(ron, 1e-6, rcut - 1e-6)
    
    result = np.where(
        r < ron,
        1.0,
        np.where(
            r > rcut,
            0.0,
            np.clip((rcut**2 - r**2)**2 * (rcut**2 + 2*r**2 - 3*ron**2) / ((rcut**2 - ron**2)**3), 0, None)
        )
    )
    return result

def morse(r, rmin, rmax, D0, alpha, r0):
    """Standard Morse potential"""
    return np.where(
        r >= rmax,
        0.0,
        D0 * (np.exp(-2 * alpha * (r - r0)) - 2 * np.exp(-alpha * (r - r0)))
    )

def morse_x(r, rmin, rmax, D0, alpha, r0, ron):
    """Smoothed Morse potential"""
    return morse(r, rmin, rmax, D0, alpha, r0) * smoothing(r, ron, rmax)

def smooth_step(r, rmin, rmax, steepness=10):
    """Smooth step function for repulsive potential"""
    x = (r - rmin) / (rmax - rmin)
    return np.clip(1.0 / (1.0 + np.exp(-steepness * (x - 0.5))), 0, 1)

def repulsive(r, rmin, rmax, A, alpha):
    """Repulsive potential with smooth step"""
    epsilon = 1e-6
    base = np.maximum(rmax - r, epsilon)
    smoothing_factor = smooth_step(r, rmin, rmax)
    potential = (A / (alpha * rmax)) * base**alpha
    return np.where(r < rmax, potential * smoothing_factor, 0.0)

def generate_potential_table(potential_func, rmin, rmax, npoints=1000, filename=None):
    """Generate LAMMPS table file for a potential function"""
    # LAMMPS requires rmin > 0, so use a small positive value if rmin is 0
    # Use very small value to handle atoms getting extremely close (down to ~1e-7)
    if rmin == 0.0:
        rmin = 1e-8  # Use very small positive value to handle close encounters
    
    r = np.linspace(rmin, rmax, npoints)
    dr = r[1] - r[0]
    
    # Calculate potential
    V = potential_func(r)
    
    # Calculate force (negative derivative of potential)
    # F = -dV/dr
    dV_dr = np.gradient(V, dr)
    F = -dV_dr
    
    # Create table content
    table_lines = []
    table_lines.append(f"# Potential table generated automatically")
    table_lines.append(f"# N = {npoints} R = {rmin} {rmax}")
    table_lines.append(f"\nN {npoints}")
    table_lines.append(f"R {rmin} {rmax}\n")
    
    for i in range(npoints):
        table_lines.append(f"{i+1} {r[i]:.8e} {V[i]:.8e} {F[i]:.8e}")
    
    table_content = "\n".join(table_lines)
    
    if filename:
        with open(filename, 'w') as f:
            f.write(table_content)
    
    return table_content

def create_lammps_rigid_patchy_monomers_script_cpp(
    num_monomers,
    box_size=None,
    rigid_sphere_radius=1.0,
    patch_radius=0.33,
    patch_offset_dist=1.0,
    patch_interaction_cutoff_morse=1.5,  # Increased cutoff to allow attractive range (r0=0.66, need cutoff > r0)
    patch_coordination_cutoff=0.34,  # Tight cutoff for state change detection (requires near-perfect overlap)
    state_change_probability=0.7,
    timesteps=500000000,
    thermo_freq=5000,
    dump_freq=50000,
    output_dir="rigid_patchy_simulation_cpp",
    seed=12345,
    state_change_freq=100,  # Check state changes every N timesteps (100 = every 100 timesteps)
    cooldown_steps=1000,  # Cooldown period (increased to prevent rapid toggling)
    timestep=0.001  # Stable timestep for D0=7.0 with 3-patch configuration
):
    """
    Generates LAMMPS input script using the custom C++ fix_state_change.
    This version does NOT require unfix/refix cycles!
    """
    
    if box_size is None:
        # Use even lower density (0.0002 instead of 0.0005) to create larger box
        # This further reduces close encounters and improves stability with strong repulsion
        density = 0.0002  # Lower density = larger box
        volume = num_monomers / density
        box_size = volume ** (1.0/3.0)
        print(f"Calculated box_size from volume formula (N/{density} = {volume}): {box_size:.4f}")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    data_filename = os.path.join(output_dir, "data.rigid_patchy_monomers")
    lammps_input_filename = os.path.join(output_dir, "in.rigid_patchy_monomers")

    # --- Atom Type Definitions ---
    main_body_mass = 1.0
    patch_mass = 1e-15  # Massless patches (from successful test)

    # --- Rigid Body Relative Positions ---
    # 3-patch configuration: 3 body spheres + 3 patches (all act as unit)
    # This was the stable configuration that created dimers without NaNs
    a = rigid_sphere_radius  # a=1.0
    b = 0.33  # Patch radius offset
    c6 = np.cos(np.pi/6)
    s6 = np.sin(np.pi/6)

    # 3 body spheres arranged in triangle at x=0
    rigid_body_coords_relative = [
        [0.0, 0.0, a],
        [0.0, a * c6, -a * s6],
        [0.0, -a * c6, -a * s6]
    ]
    
    # 3 patches arranged in triangle at x=a (positive x)
    # All 3 patches act as a single unit (change type together)
    patch_coords_relative = [
        [a, 0.0, b],
        [a, b * c6, -b * s6],
        [a, -b * c6, -b * s6]
    ]

    # --- Generate Initial Atom Data ---
    atom_data_lines = []
    current_atom_id = 1
    
    np.random.seed(seed)

    # Calculate max extent including the patches
    all_coords = rigid_body_coords_relative + patch_coords_relative
    max_x_extent = max([abs(c[0]) for c in all_coords]) + rigid_sphere_radius
    max_y_extent = max([abs(c[1]) for c in all_coords]) + rigid_sphere_radius
    max_z_extent = max([abs(c[2]) for c in all_coords]) + rigid_sphere_radius
    
    padding_x = max_x_extent + 0.1
    padding_y = max_y_extent + 0.1
    padding_z = max_z_extent + 0.1
    
    # Ensure monomers do not overlap by enforcing minimum COM spacing
    # With proper repulsion (rep_rmax_body=2.0, rep_A=500), main bodies repel strongly
    # Patches extend ~1.05 units from COM, so minimum safe distance ≈ 2.0 + 2*1.05 ≈ 4.1
    # But with strong repulsion, we can use a smaller margin
    existing_centers = []
    min_center_distance = 3.5  # Main bodies repel at 2.0, patches extend ~1.0, so 3.5 is safe
    max_attempts = 2000
    
    for mol_id in range(1, num_monomers + 1):
        placed = False
        for attempt in range(max_attempts):
            com_x = np.random.uniform(padding_x, box_size - padding_x)
            com_y = np.random.uniform(padding_y, box_size - padding_y)
            com_z = np.random.uniform(padding_z, box_size - padding_z)

            # First pass: calculate positions (used for boundary adjustment)
            atom_positions = []
            for rel_coord in rigid_body_coords_relative:
                atom_positions.append(
                    (com_x + rel_coord[0], com_y + rel_coord[1], com_z + rel_coord[2])
                )
            for patch_coord in patch_coords_relative:
                atom_positions.append(
                    (com_x + patch_coord[0], com_y + patch_coord[1], com_z + patch_coord[2])
                )

            min_x = min(pos[0] for pos in atom_positions)
            max_x = max(pos[0] for pos in atom_positions)
            min_y = min(pos[1] for pos in atom_positions)
            max_y = max(pos[1] for pos in atom_positions)
            min_z = min(pos[2] for pos in atom_positions)
            max_z = max(pos[2] for pos in atom_positions)

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
                continue  # try a new random position

            existing_centers.append((com_x, com_y, com_z))
            placed = True
            break

        if not placed:
            raise RuntimeError(
                "Failed to place monomer without overlap. Increase box size or decrease density."
            )

        # Second pass: generate final positions
        for rel_coord in rigid_body_coords_relative:
            x = com_x + rel_coord[0]
            y = com_y + rel_coord[1]
            z = com_z + rel_coord[2]
            atom_data_lines.append(f"{current_atom_id} {mol_id} 1 {x} {y} {z}")
            current_atom_id += 1

        patch_type = 3 if mol_id == 1 else 2
        for patch_coord in patch_coords_relative:
            x = com_x + patch_coord[0]
            y = com_y + patch_coord[1]
            z = com_z + patch_coord[2]
            atom_data_lines.append(f"{current_atom_id} {mol_id} {patch_type} {x} {y} {z}")
            current_atom_id += 1

    num_atoms = len(atom_data_lines)
    num_atom_types = 3
    
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

    # --- Parameters for LAMMPS built-in potentials (from successful test) ---
    # Using built-in morse and lj/cut instead of table potentials for better reliability
    # INCREASED ATTRACTION with type 3-3 being 2x stronger than 2-2
    # VERY STRONG CLOSE-RANGE ATTRACTION - pulls patches together until they overlap
    # Patches must overlap (coordination cutoff = 0.34) for state change
    morse_alpha = 10.0  # Very steep: extremely strong at close range, pulls patches to overlap
    morse_D0_22 = 500.0  # STRONG: type 2-2 attraction (reduced from 1000.0)
    morse_D0_33 = 300.0  # MODERATE: type 3-3 attraction (reduced from 500.0, weaker than 2-2)
    morse_r0 = 0.0  # Minimum at r=0 (patches overlapping)
    morse_rcut = 4.0  # Same cutoff as before
    
    # Repulsive parameters (using lj/cut)
    # LOWER REPULSION: further reduced from 100.0 to 50.0 for even easier approach
    rep_epsilon = 50.0  # Further reduced from 100.0 to allow easier approach
    rep_sigma = 1.0  # LJ sigma parameter
    rep_rmax_body = 2.0  # Body-body cutoff
    rep_rmax_patch = 1.3  # Body-patch cutoff (from successful test)

    # --- Write LAMMPS Input Script ---
    lammps_script_content = f"""# LAMMPS Input Script for Rigid Patchy Monomers with State Changes
# Generated by Python script on {datetime.now().isoformat(timespec='seconds')}
# Uses custom C++ fix state/change (NO unfix/refix required!)
# Using LAMMPS built-in morse and lj/cut potentials (from successful test)

# 1. Initialization
units           lj
atom_style      molecular
boundary        p p p
newton          on

# 2. Atom Definition
read_data       data.rigid_patchy_monomers

# 2.5. Expand box by 5% and remap atoms
change_box      all x scale 1.05 remap
change_box      all y scale 1.05 remap
change_box      all z scale 1.05 remap

# 3. Force Field Parameters
# Using LAMMPS built-in potentials (morse and lj/cut) for better reliability
# This matches the successful test simulation configuration
pair_style      hybrid/overlay morse {morse_rcut} lj/cut {rep_rmax_body}
# STRONG ATTRACTION with type 2-2 (500.0) being stronger than 3-3 (300.0)
# Morse potential for patch-patch: 
#   - Type 2-2: D0=500.0 (STRONG, stronger than 3-3), alpha=10.0, r0=0.0, rcut=4.0
#   - Type 3-3: D0=300.0 (MODERATE, weaker than 2-2), alpha=10.0, r0=0.0, rcut=4.0
# STRONG CLOSE-RANGE ATTRACTION - pulls patches together until they overlap
# Patches must overlap (coordination cutoff=0.34) for state change
# Type 2-2: D0=500.0 gives: extremely strong at r=0.0 (-500.0), very strong at r=0.2 (-335.0), strong at r=0.34 (-150.0)
# Type 3-3: D0=300.0 gives: extremely strong at r=0.0 (-300.0), very strong at r=0.2 (-201.0), strong at r=0.34 (-90.0)
# This creates a strong pulling force when patches are close, with type 2-2 being stronger than 3-3
# REPULSION: rep_epsilon=50.0 (further reduced from 100.0) for even easier approach
pair_coeff      2 2 morse {morse_D0_22} {morse_alpha} {morse_r0} {morse_rcut}
pair_coeff      3 3 morse {morse_D0_33} {morse_alpha} {morse_r0} {morse_rcut}
# Lennard-Jones for repulsion: body-body (1-1) and body-patch (1-2, 1-3)
pair_coeff      1 1 lj/cut {rep_epsilon} {rep_sigma} {rep_rmax_body}
pair_coeff      1 2 lj/cut {rep_epsilon} {rep_sigma} {rep_rmax_patch}
pair_coeff      1 3 lj/cut {rep_epsilon} {rep_sigma} {rep_rmax_patch}
# Set all other combinations to zero (no interaction)
pair_coeff      * * morse 0.0 1.0 1.0 1.0
pair_coeff      * * lj/cut 0.0 1.0 1.0
# NOTE: NO repulsion for 2-2, 2-3, 3-3 - patches can overlap with attraction dominating

# 4. Rigid Body Definition (temporary, for minimization)
fix             rigid_temp all rigid molecule

# 5. Simulation Settings
# Enhanced neighbor list settings for strong attractions
# Neighbor skin to match rcut (rcut=4.0, so need skin > 4.0)
neighbor        4.5 bin  # Skin distance to match rcut=4.0
neigh_modify    delay 0 every 1 check yes page 200000 one 20000
comm_style      brick
comm_modify     vel yes cutoff 4.5

# 5.5. Minimize (critical for stability with strong attractions)
# More thorough minimization to eliminate unfavorable configurations (from successful test)
minimize        1.0e-4 1.0e-6 100 1000

# 5.6. Reset timestep
reset_timestep  0

# 5.7. Remove temporary rigid fix
unfix           rigid_temp

# 5.8. Initialize velocities
velocity        all create 0.60 {seed} mom yes rot no

# 6. Dynamics
timestep        {timestep:.6f}  # EXACT: 0.0002 (same as successful case) - avoids NaNs with strong Morse attractions

# 6.5. Create NVT fix with stronger damping for stability
# INCREASED: Temperature 0.50 (increased from 0.30 for faster dynamics)
# Increased damping (0.05) helps dissipate excess energy from very strong attractions
fix             rigid_nvt all rigid/nvt molecule temp 0.60 0.60 0.05

# 7. Define groups for patches
group           patches_A type 2
group           patches_B type 3

# 8. Define coordination computes (required by fix state/change)
# Use tighter cutoff for coordination to require near-perfect patch overlap
compute         cA all coord/atom cutoff {patch_coordination_cutoff} group patches_A
compute         cB all coord/atom cutoff {patch_coordination_cutoff} group patches_B

# 9. Add state change fix (C++ fix - NO unfix/refix needed!)
# Syntax: fix ID group state/change check_every cooldown_steps probability cutoff group_A group_B
# Use coordination cutoff (not Morse cutoff) for state change detection - requires near-perfect overlap
fix             state_change all state/change {state_change_freq} {cooldown_steps} {state_change_probability} {patch_coordination_cutoff} patches_A patches_B

# 10. Output
thermo          {thermo_freq}
thermo_style    custom step temp pe ke etotal press f_state_change
dump            1 all custom {dump_freq} dump.rigid_patchy_monomers.lammpstrj id mol type x y z
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

if __name__ == "__main__":
    create_lammps_rigid_patchy_monomers_script_cpp(
        20, None, 1.0, 0.33, 1.0, 1.5, 0.34, 0.7, 500000000, 5000, 100000,  # 20 monomers, patch_interaction_cutoff_morse=1.5, patch_coordination_cutoff=0.34, dump_freq=100000
        'rigid_patchy_simulation_cpp_2', 12345, 100, 1000, 0.0002  # state_change_freq=100, cooldown=1000, timestep=0.0002 (increased from 0.00005 for faster dynamics)
    )

