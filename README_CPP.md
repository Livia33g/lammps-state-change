# C++ State Change Implementation for LAMMPS

This directory contains the C++ implementation of the state change fix for rigid patchy monomers in LAMMPS.

## Overview

The C++ implementation (`fix_state_change`) provides a more robust and efficient way to handle state changes compared to the Python-based approach. Key advantages:

- **No unfix/refix cycles**: Changes atom types directly without breaking rigid body constraints
- **Better performance**: Compiled C++ code is much faster than Python loops
- **More stable**: Proper MPI communication and atom tracking prevents "lost atoms" errors
- **Dimer-aware**: Only allows state changes when molecules are in dimer configuration (COM distance < 3.0)

## Files

### Core C++ Files
- `fix_state_change/fix_state_change.cpp` - Main implementation
- `fix_state_change/fix_state_change.h` - Header file

### Generation Scripts
- `generate_change_cpp.py` - Python script to generate LAMMPS input files for C++ fix
- `submit_cpp_fix.slurm` - SLURM submission script

### Installation Documentation
- `fix_state_change/INSTALL_STEP_BY_STEP.md` - Detailed installation guide
- `fix_state_change/INSTALL.md` - Quick installation reference
- `fix_state_change/README.md` - General documentation

## Quick Start

### 1. Install the C++ Fix in LAMMPS

See `fix_state_change/INSTALL_STEP_BY_STEP.md` for detailed instructions. Summary:

1. Copy `fix_state_change.cpp` and `fix_state_change.h` to your LAMMPS `src/` directory
2. Add `fix_state_change.cpp` to the Makefile
3. Rebuild LAMMPS: `make mpi`

### 2. Generate Input Files

```bash
python generate_change_cpp.py
```

Or use the SLURM script:
```bash
sbatch submit_cpp_fix.slurm
```

### 3. Run Simulation

The SLURM script will automatically:
- Generate input files if needed
- Run LAMMPS with the custom fix
- Save trajectory and log files

## Key Parameters

The `generate_change_cpp.py` script accepts various parameters:

- `num_monomers`: Number of rigid patchy monomers (default: 40)
- `morse_D0_22`: Morse well depth for type 2-2 interactions (default: 500.0)
- `morse_D0_33`: Morse well depth for type 3-3 interactions (default: 300.0)
- `temperature`: Simulation temperature (default: 0.60)
- `state_change_probability`: Probability of state change when conditions met (default: 0.7)
- `dump_freq`: Trajectory output frequency (default: 1000 steps)

## Fix Syntax

```
fix ID group-ID state/change check_every cooldown_steps probability cutoff group_A group_B
```

Example:
```
compute cA all coord/atom cutoff 0.34 group patches_A
compute cB all coord/atom cutoff 0.34 group patches_B
fix state_change all state/change 100 1000 0.7 0.34 patches_A patches_B
```

## Differences from Python Implementation

1. **No unfix/refix**: The C++ fix changes atom types directly, so rigid bodies remain intact
2. **Dimer check**: Only allows state changes when molecules are in dimer (COM distance < 3.0)
3. **Better coordination**: Uses `coord/atom` compute for more accurate coordination detection
4. **MPI-safe**: Proper communication ensures no "lost atoms" errors

## Troubleshooting

### "Lost atoms" error
- This should NOT happen with the C++ fix
- If it does, check that atom types are correctly initialized (2 or 3 for patches)
- Verify groups are correctly defined

### State changes not occurring
- Check that patches are overlapping (coordination cutoff = 0.34)
- Verify molecules are in dimer (COM distance < 3.0)
- Check that `state_change_probability` is not too low
- Ensure `cooldown_steps` is not too long

### Compilation errors
- Make sure you're using LAMMPS 2023-08-17 or later
- Check that all required packages are enabled (`make yes-rigid yes-molecule`)
- Verify MPI is properly configured

## Output Files

- `dump.rigid_patchy_monomers.lammpstrj` - Trajectory file
- `log.lammps` - LAMMPS log file
- `lammps_stdout.log` - Standard output

## Analysis

Use `analyze_state_changes_with_dimers.py` to analyze state changes and dimer formation from the trajectory.

