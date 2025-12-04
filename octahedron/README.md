# Octahedron Monomers with State Changes

This directory contains LAMMPS simulation scripts for rigid octahedron monomers with dynamic state changes using the custom C++ `fix_state_change`.

## Overview

The octahedron monomers have:
- **1 body (vertex)**: Single body particle (type 1) - this monomer represents one vertex of the octahedron
- **4 patches**: Patches positioned around the body (types 2 initially, can evolve to 3, 4, or 5)
- **Total**: 5 particles per monomer (1 body + 4 patches)

The full octahedron structure is formed by assembling 6 such monomers:
- Each monomer = 1 vertex of the octahedron
- 6 monomers × 5 particles = 30 total particles in a complete octahedron

The structure matches the JAX MD implementation:
- Each monomer has 1 body particle at its center (COM)
- 4 patches are positioned around the body, pointing toward nearest neighbors
- All 4 patches in a monomer are always the same type (change together as a unit)

## State Change Mechanism

Uses the same C++ `fix_state_change` as the rigid patchy simulation:
- Patches of type 2 can change to type 3 when coordinated (within cutoff)
- Patches of type 3 can change to type 2 when coordinated (within cutoff)
- Only occurs when patches are in coordination (cutoff = 0.34)
- Probabilistic with cooldown period to prevent rapid toggling

## Files

- `generate_octahedron_cpp.py` - Python script to generate LAMMPS input files
- `submit_octahedron.slurm` - SLURM submission script
- `octahedron_simulation_cpp/` - Output directory (created on first run)

## Quick Start

### 1. Ensure LAMMPS is Built with fix_state_change

The C++ fix must be installed in your custom LAMMPS build. See `../fix_state_change/INSTALL_STEP_BY_STEP.md` for details.

Verify it's available:
```bash
/work/nvme/bewl/lguttieres/lammps_build/lammps/src/lmp_mpi -help | grep state/change
```

### 2. Generate and Run Simulation

From the `octahedron/` directory:

```bash
# Submit job (will generate input files if needed)
sbatch submit_octahedron.slurm

# Or generate manually first
python3 generate_octahedron_cpp.py
```

### 3. Monitor Progress

```bash
# Check job status
squeue -u $USER

# Watch output
tail -f slurm_octahedron-*.out

# Check simulation output
tail -f octahedron_simulation_cpp/lammps_stdout.log
```

## Key Parameters

### Geometry Parameters
- `num_monomers`: Number of octahedron monomers (default: 50)
- `vertex_radius`: Radius of vertex spheres (default: 2.0)
- `patch_radius`: Radius of patch spheres (default: 0.5, for reference)
- `box_size`: Simulation box size (auto-calculated from density if None)

### State Change Parameters
- `state_change_probability`: Probability when conditions met (default: 0.7)
- `patch_coordination_cutoff`: Distance cutoff for coordination detection (default: 0.34)
- `state_change_freq`: Check for state changes every N steps (default: 100)
- `cooldown_steps`: Minimum steps between state changes for same patch (default: 1000)

### Interaction Parameters
- `morse_D0_22`: Morse well depth for type 2-2 interactions (default: 10.0)
- `morse_D0_33`: Morse well depth for type 3-3 interactions (default: 10.0)
- `morse_D0_23`: Cross-interaction (default: 0.0, no interaction)
- `rep_epsilon`: Repulsion strength for body-body interactions (default: 10000.0)
- `temperature`: Simulation temperature (default: 1.0)
- `morse_cutoff`: Morse potential cutoff distance (default: 4.0)

### Simulation Parameters
- `timesteps`: Total simulation steps (default: 500000000)
- `timestep`: Integration timestep (default: 0.001)
- `dump_freq`: Trajectory output frequency (default: 10000)
- `thermo_freq`: Thermodynamic output frequency (default: 5000)

## Differences from Rigid Patchy Simulation

1. **Geometry**: Octahedron monomer (1 body, 4 patches) vs 3-body spheres (3 bodies, 3 patches)
2. **Particles per monomer**: 5 vs 6
3. **Patch arrangement**: 4 patches (all change together) vs 3 patches (can change individually)
4. **Assembly**: 6 monomers form full octahedron structure
4. **Masses**: Vertex mass = 0.5 (from JAX code) vs body mass = 1.0

## Output Files

After running, `octahedron_simulation_cpp/` contains:
- `data.octahedron_monomers` - Initial configuration
- `in.octahedron_monomers` - LAMMPS input script
- `dump.octahedron_monomers.lammpstrj` - Trajectory file
- `log.lammps` - LAMMPS log file
- `lammps_stdout.log` - Standard output

## Customization

To modify parameters, edit the function call in `submit_octahedron.slurm` or run `generate_octahedron_cpp.py` directly with custom arguments:

```python
from generate_octahedron_cpp import create_lammps_octahedron_script_cpp

create_lammps_octahedron_script_cpp(
    num_monomers=100,
    temperature=0.8,
    morse_D0_22=15.0,
    morse_D0_33=15.0,
    # ... other parameters
)
```

## Troubleshooting

### "Unknown fix style: state/change"
- Rebuild LAMMPS with the fix installed (see installation docs)

### State changes not occurring
- Check that patches are coordinating (cutoff may be too tight)
- Verify `state_change_probability` is not too low
- Check `cooldown_steps` is not too long

### Simulation instability
- Try reducing `timestep`
- Increase `rep_epsilon` for stronger repulsion
- Reduce Morse well depths (`morse_D0_22`, `morse_D0_33`)

## Notes

- The octahedron geometry matches the JAX MD implementation structure
- State changes preserve the octahedron rigid body structure (no unfix/refix needed)
- The C++ fix handles all state changes automatically during dynamics

