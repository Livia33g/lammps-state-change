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
    # LAMMPS pair_style table expects:
    #   <section_name>
    #   N <npoints> R <rmin> <rmax>
    #   i  r  V  F
    table_lines = []
    table_lines.append(f"# Potential table generated automatically")
    table_lines.append(f"# N = {npoints} R = {rmin} {rmax}")
    # Section label required by LAMMPS table pair style
    table_lines.append("REPULSIVE_11")
    table_lines.append(f"N {npoints} R {rmin} {rmax}")
    
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
    patch_offset_dist=None,
    patch_interaction_cutoff_morse=None,
    patch_coordination_cutoff=0.34,  # Tight cutoff for state change detection (requires near-perfect overlap)
    state_change_probability=0.7,
    timesteps=500000000,
    thermo_freq=5000,
    dump_freq=50000,
    output_dir="rigid_patchy_simulation_cpp",
    seed=12345,
    state_change_freq=100,  # Check state changes every N timesteps (100 = every 100 timesteps)
    cooldown_steps=1000,  # Cooldown period (increased to prevent rapid toggling)
    timestep=0.001,  # Stable timestep for D0=7.0 with 3-body configuration (can be overridden)
    morse_D0_22=None,  # Optional: override morse D0 for type 2-2 (default: 13.0)
    morse_D0_33=None,  # Optional: override morse D0 for type 3-3 (default: 13.0)
    rep_epsilon=None,  # Optional: override repulsion epsilon (default: 500.0)
    temperature=None  # Optional: override temperature (default: 1.0)
):
    """
    Generates LAMMPS input script using the custom C++ fix_state_change.
    This version does NOT require unfix/refix cycles!
    """
    
    if box_size is None:
        # Use higher density (0.001) for higher concentration
        # Higher concentration = smaller box = more frequent collisions
        density = 0.001  # Higher density = smaller box = higher concentration
        volume = num_monomers / density
        box_size = volume ** (1.0/3.0)
        print(f"Calculated box_size from volume formula (N/{density} = {volume}): {box_size:.4f}")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    data_filename = os.path.join(output_dir, "data.rigid_patchy_monomers")
    lammps_input_filename = os.path.join(output_dir, "in.rigid_patchy_monomers")

    # --- Atom Type Definitions ---
    # Back to original 3-body geometry:
    #   - Three main spheres of radius 1.0 (type 1), each with mass 1.0
    #   - One patch per monomer (type 2 or 3), essentially massless
    main_body_mass = 1.0
    patch_mass = 1.0e-3  # Small but non-zero to avoid singular inertia

    # --- Rigid Body Relative Positions ---
    # Geometry MATCHED to the JAX optimize trajectory (mon_shape1 in optimize.py):
    #   - Three main spheres (body atoms) at:
    #       [0, 0, a]
    #       [0,  a*cos(pi/6), -a*sin(pi/6)]
    #       [0, -a*cos(pi/6), -a*sin(pi/6)]
    #   - Three patches at:
    #       [a, 0, b]
    #       [a,  b*cos(pi/6), -b*sin(pi/6)]
    #       [a, -b*cos(pi/6), -b*sin(pi/6)]
    # where a = 1.0 (body radius) and b = patch_radius = 0.33.
    a = rigid_sphere_radius
    b = patch_radius
    c6 = np.cos(np.pi / 6.0)  # cos(30 deg)
    s6 = np.sin(np.pi / 6.0)  # sin(30 deg)

    # Three body spheres
    rigid_body_coords_relative = [
        [0.0, 0.0, a],
        [0.0,  a * c6, -a * s6],
        [0.0, -a * c6, -a * s6],
    ]

    # Three patches on the +x side, identical to optimize.py
    # (we ignore patch_offset_dist here to keep 1–1 matching with the JAX geometry)
    patch_coords_relative = [
        [a, 0.0, b],
        [a,  b * c6, -b * s6],
        [a, -b * c6, -b * s6],
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
    # Approximate collision distance between bodies: 2 * sphere radius.
    # (3-sphere geometry is anisotropic, but this gives a reasonable isotropic scale.)
    body_body_collision_distance = 2.0 * rigid_sphere_radius
    # Total radial extent including patches along the most extended axis.
    # Use max_x_extent (which already includes +rigid_sphere_radius) as a conservative
    # estimate for the "radius" of each rigid body, then space COMs by ~2x that.
    patch_total_extent = max_x_extent
    # Require COM spacing that keeps both body spheres and patches from overlapping initially
    existing_centers = []
    min_center_distance = 2.0 * patch_total_extent + 0.5
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

        # Initialize ALL patches as type 2; type changes to 3 will be driven
        # solely by the state_change fix during the dynamics.
        patch_type = 2
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

    # --- Parameters for potentials ---
    # Goal: very short-range attraction that only overcomes repulsion when
    # patches are essentially overlapping, and even then only slightly.

    # Morse patch-patch:
    #   - center at r0 = 0 (maximum at full overlap)
    #   - very short cutoff: rcut = 0.33
    #   - steep alpha so the well is narrow
    morse_rcut_default = 0.33
    morse_rcut = patch_interaction_cutoff_morse if patch_interaction_cutoff_morse is not None else morse_rcut_default
    morse_alpha = 10.0
    morse_r0 = 0.0
    # Moderate well depth so that 3 bonds modestly beat repulsion at full overlap
    morse_D0_22 = morse_D0_22 if morse_D0_22 is not None else 3.0
    morse_D0_33 = morse_D0_33 if morse_D0_33 is not None else 3.0

    # Body-body repulsion: LJ tuned so minimum is ~at contact (2R),
    # with epsilon small compared to 3*D0 so that at full patch overlap
    # the net attraction is dominated by patches, but still providing
    # a short-range core when patches are not engaged.
    #   epsilon ≈ 2
    #   sigma chosen so r_min = 2^(1/6) * sigma ≈ 2R
    #   cutoff slightly beyond contact for short-range repulsion.
    rep_epsilon = rep_epsilon if rep_epsilon is not None else 2.0
    rep_sigma = body_body_collision_distance / (2.0 ** (1.0 / 6.0))
    rep_rmax_body = body_body_collision_distance + 0.5
    
    # Temperature parameter
    temp_value = temperature if temperature is not None else 1.0  # Target kT ~ 1
    init_temp_value = 0.3 * temp_value  # Cold start to minimize initial KE

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
# Using hybrid/overlay:
#   - morse for patch-patch attraction (types 2 and 3)
#   - lj/cut for body-body repulsion (steep, short-ranged)
pair_style      hybrid/overlay morse {morse_rcut} lj/cut {rep_rmax_body}
# First zero all interactions, then override the pairs we actually use.
pair_coeff      * * morse 0.0 1.0 1.0 1.0
pair_coeff      * * lj/cut 0.0 1.0 1.0
pair_coeff      2 2 morse {morse_D0_22} {morse_alpha} {morse_r0} {morse_rcut}
pair_coeff      3 3 morse {morse_D0_33} {morse_alpha} {morse_r0} {morse_rcut}
# Body-body LJ repulsion only (no body-patch interactions)
pair_coeff      1 1 lj/cut {rep_epsilon} {rep_sigma:.4f} {rep_rmax_body}
# NOTE: No body-patch interactions (1-2, 1-3); patches only attract via short-range Morse.

# 4. Rigid Body Definition (temporary, for minimization)
fix             rigid_temp all rigid molecule

# 5. Simulation Settings
# Neighbor list tuned to the short-ranged interactions:
#   - Morse cutoff ~{morse_rcut}
#   - LJ cutoff ~{rep_rmax_body}
neighbor        {max(morse_rcut, rep_rmax_body) + 0.5:.2f} bin
neigh_modify    delay 0 every 1 check yes page 200000 one 20000
comm_style      brick
comm_modify     vel yes cutoff {max(morse_rcut, rep_rmax_body) + 0.5:.2f}

# 5.5. Minimize (critical for stability with strong attractions)
# More thorough minimization to eliminate unfavorable configurations (from successful test)
minimize        1.0e-4 1.0e-6 100 1000

# 5.6. Reset timestep
reset_timestep  0

# 5.7. Remove temporary rigid fix
unfix           rigid_temp

# 5.8. Initialize velocities
# Initialize at a lower temperature (init_temp_value) to avoid a hot start;
# thermostat will then heat/cool towards target temp_value.
velocity        all create {init_temp_value} {seed} mom yes rot no

# 6. Dynamics
timestep        {timestep:.6f}

# 6.5. Rigid-body integrator and thermostat
# Use fix rigid/nvt with strong (but not extreme) damping to keep temperature
# near target while allowing some relaxation of torques.
fix             rigid_nvt all rigid/nvt molecule temp {temp_value} {temp_value} 200.0

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
        20, None, 2.0, 0.33, 2.0, 0.66, 0.34, 0.7, 500000000, 5000, 100000,  # Example standalone run with updated geometry
        'rigid_patchy_simulation_cpp_2', 12345, 100, 1000, 0.0002  # state_change_freq=100, cooldown=1000, timestep=0.0002
    )

