# LAMMPS Fix: State Change for Rigid Patchy Monomers

This custom LAMMPS fix enables seamless state changes for patch atoms in rigid body simulations without requiring unfix/refix cycles.

## Features

- **No unfix/refix required**: Changes atom types directly without breaking rigid body constraints
- **Coordination-based**: State changes occur when patches are coordinated with same-type patches
- **Cooldown mechanism**: Prevents rapid toggling with configurable cooldown period
- **Probabilistic**: State changes occur with a configurable probability
- **MPI-safe**: Properly handles parallel simulations

## Installation

### Option 1: Build as LAMMPS package (Recommended)

1. Copy the fix files to your LAMMPS source:
   ```bash
   cp fix_state_change.h /path/to/lammps/src/
   cp fix_state_change.cpp /path/to/lammps/src/
   ```

2. Edit `/path/to/lammps/src/Makefile` and add `fix_state_change` to the list of source files.

3. Rebuild LAMMPS:
   ```bash
   cd /path/to/lammps/src
   make mpi
   ```

### Option 2: Build as shared library

1. Set `LAMMPS_DIR` in the Makefile to point to your LAMMPS installation
2. Run `make` to build `libfix_state_change.a`
3. Link this library when building LAMMPS

## Usage

### Syntax

```
fix ID group-ID state/change check_every cooldown_steps probability cutoff group_A group_B
```

- `ID`: Fix identifier
- `group-ID`: Group of atoms to apply fix to (typically `all`)
- `check_every`: Check for state changes every N timesteps (1 = every timestep)
- `cooldown_steps`: Cooldown period - patches can't change again for this many steps after changing
- `probability`: Probability (0-1) of state change when conditions are met
- `cutoff`: Cutoff distance for coordination check (should match Morse potential cutoff)
- `group_A`: Group name for type 2 patches (e.g., `patches_A`)
- `group_B`: Group name for type 3 patches (e.g., `patches_B`)

### Required Computes

The fix requires two coordination computes to be defined:

```
compute cA all coord/atom cutoff <cutoff> group patches_A
compute cB all coord/atom cutoff <cutoff> group patches_B
```

### Example

```
# Define groups
group patches_A type 2
group patches_B type 3

# Define coordination computes
compute cA all coord/atom cutoff 2.5 group patches_A
compute cB all coord/atom cutoff 2.5 group patches_B

# Create rigid fix (no need to unfix/refix!)
fix rigid_nvt all rigid/nvt molecule temp 1.0 1.0 0.01

# Add state change fix
fix state_change all state/change 1 100 0.7 2.5 patches_A patches_B

# Run simulation - state changes happen automatically!
run 1000000
```

## How It Works

1. The fix monitors coordination numbers for each patch atom
2. Every `check_every` steps, it evaluates state change conditions:
   - Type 2 patches: Change to type 3 if coordinated with type 2 AND random < probability AND cooldown passed
   - Type 3 patches: Change to type 2 if coordinated with type 3 AND random < probability AND cooldown passed
3. When conditions are met, atom types are changed directly (no unfix/refix needed)
4. The fix maintains per-atom properties for cooldown tracking

## Advantages Over Script-Based Approach

- **No atom loss**: Avoids the "Lost atoms" error by not unfixing rigid bodies
- **Efficient**: No expensive unfix/refix cycles
- **Seamless**: State changes happen automatically during dynamics
- **MPI-safe**: Properly handles parallel communication

## Notes

- The fix changes atom types directly, which updates pair interactions automatically
- The neighbor list is flagged for rebuild after type changes
- Works seamlessly with `fix rigid/nvt` and other rigid body fixes
- Per-atom data is preserved across restarts

