# NaN Diagnosis for Octahedron Simulation

## Root Causes Identified

### 1. **CRITICAL: Patch Mass Too Small** ⚠️
- **Data file has**: `1e-06` (0.000001) 
- **Code expects**: `0.001`
- **Problem**: With mass = 1e-06, any force F causes acceleration a = F/m = F/1e-06 = 1,000,000×F
- **Result**: Huge accelerations → velocities → positions → NaNs
- **Fix**: Regenerate data file with correct mass (0.001) OR manually edit data file

### 2. **Repulsion Epsilon Too High**
- **Current**: `rep_epsilon = 400.0` for body-body
- **Problem**: If particles overlap even slightly, repulsion force = 400.0 / r^12 is HUGE
- **Fix**: Reduce to 100.0-200.0, or ensure better initial placement

### 3. **Box Size May Be Too Small**
- **Current**: Density = 0.0005 (reduced from 0.001)
- **Problem**: If monomers are too close, they overlap → huge repulsion forces
- **Fix**: Further reduce density to 0.0003 or increase box size manually

### 4. **Rigid Body Construction**
- **Check**: Ensure all atoms in a molecule are properly assigned
- **Problem**: If rigid body is malformed, forces can be wrong
- **Fix**: Verify molecule IDs in data file

## Recommended Fixes (in order of priority)

### Priority 1: Fix Patch Mass
```bash
# Regenerate data file with FORCE_REGENERATE=1
cd /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/octahedron
FORCE_REGENERATE=1 sbatch submit_octahedron.slurm
```

OR manually edit the data file:
```bash
cd octahedron_simulation_cpp
sed -i 's/2 1e-06/2 0.001/g' data.octahedron_monomers
sed -i 's/3 1e-06/3 0.001/g' data.octahedron_monomers
sed -i 's/4 1e-06/4 0.001/g' data.octahedron_monomers
sed -i 's/5 1e-06/5 0.001/g' data.octahedron_monomers
```

### Priority 2: Reduce Repulsion
In `submit_octahedron.slurm`, change:
```python
rep_epsilon=200.0,  # Reduced from 400.0
```

### Priority 3: Increase Box Size
In `submit_octahedron.slurm`, change:
```python
# Lower density = larger box
density = 0.0003  # In generate_octahedron_cpp.py
```

### Priority 4: Verify Rigid Bodies
Check that each molecule has exactly 5 atoms (1 body + 4 patches):
```bash
grep "^[0-9]* [0-9]* 1" data.octahedron_monomers | wc -l  # Should equal num_monomers
```

## Why This Happens

1. **Mass too small**: F = ma → a = F/m. If m is tiny, a is huge → v explodes → NaNs
2. **Repulsion too strong**: LJ repulsion ~1/r^12. If r is small, force is enormous
3. **Box too small**: Particles overlap → small r → huge forces
4. **Rigid body issues**: If not properly constrained, patches can "fly off" → huge forces

## Quick Test

After fixing patch mass, run a short test:
```bash
# Edit in.octahedron_monomers to run only 1000 steps
# Check if NaNs appear in first few steps
```

