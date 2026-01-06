# Octahedron Monomers with State Changes - Complete Guide

This directory contains LAMMPS simulation scripts for rigid octahedron monomers with dynamic state changes using the custom C++ `fix_state_change/octahedron`.

## Overview

The octahedron monomers have:
- **1 body particle** (type 1): Central repulsive core
- **4 patches** (types 2-5): Attractive patches positioned around the body
- **Total**: 5 particles per monomer (1 body + 4 patches)

**Key Physics Principle**: In LAMMPS, **Mass ≠ Size**. Patches can have mass for rotational stability while still allowing overlap (no repulsion between patches).

## The "Holy Trinity" of Stable State Change Simulations

All state change simulations require three critical components:

### 1. Memory Safety (MPI Communication)
**Problem**: Ghost atoms on different processors carry garbage memory, causing phantom state changes.

**Solution**: Set `comm_forward = 3` in the fix constructor to properly transmit:
- `last_change` (timestep of last state change)
- `effective_type` (current effective type)
- `prev_coord` (previous coordination number)

**Code Location**: `fix_state_change_octahedron.cpp` constructor
```cpp
comm_forward = 3;  // Must match number of values in pack_forward_comm
```

### 2. Physical Stability (Rotation & Energy)
**Problem**: Light patches cause infinite spin (moment of inertia → 0).

**Solution**: 
- **Mass Distribution**: Body = 0.6, Patches = 0.1 each (total = 1.0)
- **Strong Damping**: `Tdamp = 1.0` (smaller = stronger damping)
- **Wide Wells**: `morse_alpha = 1.2` (wider, more forgiving bonds)

**Why It Works**: Heavy patches create high moment of inertia (I = Σ m_i × r_i²), preventing wild rotation. Mass doesn't prevent overlap - only repulsion does.

### 3. Kinetic Control (Reaction Rate)
**Problem**: High probability + frequent checks = explosive cascades ("flash mob" effect).

**Solution**:
- **Low Probability**: `probability = 0.05` (5% chance, not 70%)
- **Slow Checking**: `check_every = 5000` steps (not 100)
- **Long Cooldown**: `cooldown_steps = 10000` (not 1000)

**Why It Works**: System has time to relax between state changes. Each contact is a "suggestion" that's ignored 95% of the time.

## State Change Logic

### Monochromatic Rule
**All patches in a monomer must have the same type** (enforced by consistency sweep).

**Implementation**:
1. **Trigger Phase**: Each patch is checked independently for contacts
2. **Update Phase**: If one patch triggers, ALL patches in that molecule change together
3. **Consistency Sweep**: After state changes, check all molecules (including ghosts) and enforce monochromatic rule

### State Change Rules
- **Type 1 (initial) → Types 3/4/5**: When patch touches ANY other patch
- **Types 3/4/5 → Types 3/4/5**: When patch touches SAME type (symmetry breaking via probability)

## Critical Fixes Implemented

### Fix 1: MPI Consistency Sweep
**Location**: `check_and_change()` function, end of function

Prevents "fruit salad" monomers (mixed colors) caused by MPI split-brain:
- Checks all atoms including ghosts
- Finds most recent state change per molecule
- Forces all patches to match newest type

### Fix 2: Atom-Based Triggering
**Location**: `check_and_change()` function, scan loop

Each patch is evaluated independently for contacts, but changes are applied to the whole molecule.

### Fix 3: Communication Setup
**Location**: Constructor and `pack_forward_comm`/`unpack_forward_comm`

Ensures ghost atoms receive correct memory:
```cpp
// Constructor
comm_forward = 3;

// pack_forward_comm - must send in same order as unpack
buf[m++] = static_cast<double>(last_change[j]);
buf[m++] = static_cast<double>(effective_type[j]);
buf[m++] = prev_coord[j];
```

## Common Issues and Solutions

### Issue 1: NaNs (Not a Number)
**Symptoms**: Simulation crashes with NaN values in energy/forces

**Causes & Solutions**:
1. **Timestep too large**: Reduce `timestep` from 0.001 to 0.0001 or 0.0002
2. **Morse attraction too strong**: Reduce `morse_D0` from 10.0 to 5.0
3. **Repulsion too strong**: Reduce `rep_epsilon` from 800.0 to 400.0
4. **Density too high**: Reduce density or increase box size

**Current Stable Values**:
- `timestep = 0.0001`
- `morse_D0 = 5.0` (base_epsilon)
- `rep_epsilon = 800.0`
- `density = 0.005`

### Issue 2: Segmentation Fault
**Symptoms**: Simulation crashes immediately with segfault

**Causes & Solutions**:
1. **Uninitialized arrays**: Ensure all arrays are allocated in constructor, not `init()`
2. **Missing comm_forward**: Must set `comm_forward = 3` in constructor
3. **Array access out of bounds**: Check `i < atom->nmax` before accessing arrays

**Fix**: Move all array allocation to constructor, initialize to safe defaults (-1 for last_change, 1 for effective_type).

### Issue 3: "Lost Atoms" Error
**Symptoms**: `Lost atoms: original 100 current 90`

**Causes & Solutions**:
1. **Temperature too high**: Reduce from 1.0 to 0.3
2. **Initial velocities too large**: Scale initial velocities by 0.5
3. **Thermostat too weak**: Reduce `Tdamp` from 100.0 to 1.0 (stronger damping)
4. **Minimization incompatible**: Remove minimization steps (LAMMPS doesn't support with rigid bodies)

**Current Stable Values**:
- `temperature = 0.3`
- `Tdamp = 1.0`
- No minimization

### Issue 4: Wild Spinning ("Soccer Ball" Effect)
**Symptoms**: Monomers spin wildly during state changes

**Causes & Solutions**:
1. **Patches too light**: Increase patch mass from 0.000001 to 0.1
2. **Body too light**: Increase body mass to 0.6 (total = 1.0)
3. **Thermostat too slow**: Reduce `Tdamp` to 1.0
4. **Energy spike too large**: Reduce Morse attractions or increase cooldown

**Physics**: Moment of Inertia I = Σ m_i × r_i². If patches are massless, I → 0, causing infinite rotation speed.

**Current Stable Values**:
- `body_mass = 0.6`
- `patch_mass = 0.1` each
- `Tdamp = 1.0`

### Issue 5: "Fruit Salad" Monomers (Mixed Colors)
**Symptoms**: Single monomer has patches of different colors

**Causes & Solutions**:
1. **MPI split-brain**: Add consistency sweep at end of `check_and_change()`
2. **Missing comm_forward**: Set `comm_forward = 3` in constructor
3. **Incorrect pack/unpack order**: Ensure pack and unpack use same order

**Fix**: Consistency sweep checks all atoms (including ghosts) and forces all patches in a molecule to match the newest type.

### Issue 6: Phantom State Changes (No Contact)
**Symptoms**: Monomers change type without touching anything

**Causes & Solutions**:
1. **Ghost atom garbage memory**: Set `comm_forward = 3` and ensure pack/unpack transmit `last_change`
2. **Cutoff too large**: Reduce `morse_rcut` from 2.5 to 1.5 if needed
3. **Uninitialized last_change**: Initialize to -1 in constructor

**Fix**: Proper MPI communication ensures ghost atoms have correct `last_change` values.

### Issue 7: "Flash Mob" Cascades
**Symptoms**: All monomers change type simultaneously in a cascade

**Causes & Solutions**:
1. **Probability too high**: Reduce from 0.7 to 0.05
2. **Check interval too short**: Increase from 100 to 5000 steps
3. **Cooldown too short**: Increase from 1000 to 10000 steps

**Current Stable Values**:
- `probability = 0.05`
- `check_every = 5000`
- `cooldown_steps = 10000`

### Issue 8: No Cluster Formation
**Symptoms**: Monomers float around but never form clusters

**Causes & Solutions**:
1. **Density too low**: Increase from 0.0005 to 0.005 (10x)
2. **Box too large**: Reduce box multiplier from 3.5 to 2.0
3. **Cutoff too small**: Increase `morse_rcut` to 2.5
4. **Not enough particles**: Increase `num_monomers` from 20 to 50
5. **Run too short**: Self-assembly takes 10^7-10^8 steps

**Current Stable Values**:
- `density = 0.005`
- `box_multiplier = 2.0`
- `morse_rcut = 2.5`
- `num_monomers = 50`

### Issue 9: All Monomers Become Same Type
**Symptoms**: System converges to all type 3 (or 4, or 5)

**Causes & Solutions**:
1. **Unification logic bug**: Ensure each patch changes independently (atom-based triggering)
2. **Probability = 1.0**: Should be < 1.0 for symmetry breaking
3. **No cooldown**: Cooldown prevents rapid cascades

**Fix**: Remove any "majority vote" or unification logic that forces all patches to match.

## Production-Quality Parameters

### Current Stable Configuration
```python
# Physical Parameters
body_mass = 0.6
patch_mass = 0.1  # each (4 patches × 0.1 = 0.4 total)
temperature = 0.3
timestep = 0.0001
Tdamp = 1.0  # Strong damping

# Interaction Parameters
morse_alpha = 1.2  # Wide well
morse_rcut = 2.5   # Large capture radius
base_epsilon = 5.0  # Morse well depth
rep_epsilon = 800.0  # Body repulsion

# State Change Parameters
probability = 0.05      # 5% chance
check_every = 5000     # Check every 5000 steps
cooldown_steps = 10000 # 10k step cooldown

# System Parameters
density = 0.005
box_multiplier = 2.0
num_monomers = 50
```

## Quick Start

### 1. Build LAMMPS with Fix
```bash
cd /work/nvme/bewl/lguttieres/lammps_build/lammps/src
cp /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/octahedron/fix_state_change_octahedron.* .
make mpi -j 4
```

Verify installation:
```bash
./lmp_mpi -help | grep "state/change/octahedron"
```

### 2. Run Simulation
```bash
cd /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/octahedron
FORCE_REGENERATE=1 sbatch submit_octahedron.slurm
```

### 3. Monitor Progress
```bash
# Check job status
squeue -u $USER

# Watch output
tail -f octahedron_simulation_cpp/lammps_stdout.log | grep -E "Step|Temp|f_state"

# Check for errors
tail -f slurm_octahedron-*.err
```

## Troubleshooting Checklist

When starting a new state change simulation, verify:

- [ ] `comm_forward` is set correctly (matches number of values in pack_forward_comm)
- [ ] All arrays allocated in constructor, not `init()`
- [ ] Arrays initialized to safe defaults (-1 for last_change)
- [ ] Consistency sweep implemented for monochromatic rule
- [ ] Probability is low (0.05-0.1, not 0.7)
- [ ] Check interval is long (5000+, not 100)
- [ ] Cooldown is long (10000+, not 1000)
- [ ] Mass distribution provides rotational stability
- [ ] Tdamp is strong (1.0, not 100.0)
- [ ] Morse cutoff is reasonable (1.5-2.5)
- [ ] Density allows collisions but not instant gelation

## Expected Behavior

### Initial Phase (0-100k steps)
- Slow, gradual evolution
- Rare state changes (5% probability)
- Monomers jostle and explore

### Assembly Phase (100k-10M steps)
- Controlled dimerization
- Gradual cluster growth
- No instant cascades

### Visual Check (OVITO/VMD)
- ✅ Monochromatic monomers (all patches same color)
- ✅ No phantom changes
- ✅ Smooth, controlled evolution
- ❌ Mixed colors = MPI bug
- ❌ Instant cascades = probability too high

## Files

- `generate_octahedron_cpp.py` - Generates LAMMPS input files
- `submit_octahedron.slurm` - SLURM submission script
- `fix_state_change_octahedron.cpp` - Custom C++ fix (with all critical fixes)
- `fix_state_change_octahedron.h` - Header file
- `octahedron_simulation_cpp/` - Output directory

## Notes

- **Mass ≠ Size**: Patches can be heavy (for rotation) but still overlap (no repulsion)
- **Slow is Good**: Low probability + long intervals = stable, controlled reactions
- **MPI is Critical**: Always set `comm_forward` and implement consistency sweep
- **Patience Required**: Self-assembly takes 10^7-10^8 steps with slow reaction rates

## References

This implementation is based on extensive debugging and optimization. Key lessons:
1. Memory safety (MPI communication) is non-negotiable
2. Physical stability (mass + damping) prevents numerical disasters
3. Kinetic control (slow reactions) prevents logical cascades

For questions or issues, refer to this README's troubleshooting section first.
