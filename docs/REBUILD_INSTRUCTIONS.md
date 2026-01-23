# Rebuilding LAMMPS with All State Change Fixes

This guide explains how to rebuild LAMMPS to include all custom state change fixes:
- `state/change/dimer` (dimer flip fix; source: `dimer/fix_state_change_dimer.*`)
- `state/change/dimer_ksat` (A->C catalyzed by B contact; source: `dimer_ksat/variants/base/fix_state_change_dimer_ksat.*`)
- `state/change/dimer_ksat_twoside` (variant; source: `dimer_ksat/variants/1core_twosideB/fix_state_change_dimer_ksat_twoside.*`)
- `state/change/octahedron` (octahedron fix; source: `octahedron/fix_state_change_octahedron.*`)
- `state/change/ksat` (ksat fix; source: `ksat/fix_state_change_ksat.*`)

## Prerequisites

- Access to an interactive compute node
- LAMMPS source code at `/work/nvme/bewl/lguttieres/lammps_build/lammps/src`
- All fix source files in their respective directories

## Method: Interactive Compute Node

**The recommended approach is to use an interactive compute node** because:
- MPI compilers are available
- Environment modules work correctly
- Build completes without module loading issues

### Step 1: Get an Interactive Node

You have two options:

**Option A: SSH to gpuA node directly**
```bash
ssh gpuA
```

**Option B: Request an interactive SLURM session**
```bash
salloc --partition=gpuA100x4 --account=bewl-delta-gpu --nodes=1 --ntasks-per-node=1 --gres=gpu:1 --time=01:00:00
```

### Step 2: Set Up Environment

On the interactive node, run:

```bash
# Deactivate conda if active (can cause compiler conflicts)
conda deactivate 2>/dev/null || true

# Add Cray MPI compiler to PATH
export PATH=/opt/cray/pe/mpich/8.1.32/ofi/gnu/11.2/bin:$PATH

# Verify mpicxx is available
which mpicxx
mpicxx --version | head -1
```

### Step 3: Navigate to LAMMPS Source

```bash
cd /work/nvme/bewl/lguttieres/lammps_build/lammps/src
```

### Step 4: Copy All Fix Files

Copy the fix files you want compiled into LAMMPS:

```bash
# Dimer fix (current)
cp /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/dimer/fix_state_change_dimer.{cpp,h} .

# Dimer_ksat fix (A->C catalyzed by B contact)
cp /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/dimer_ksat/variants/base/fix_state_change_dimer_ksat.{cpp,h} .

# Dimer_ksat_twoside fix (A->C only if B has A on both faces)
cp /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/dimer_ksat/variants/1core_twosideB/fix_state_change_dimer_ksat_twoside.{cpp,h} .

# Octahedron fix
cp /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/octahedron/fix_state_change_octahedron.{cpp,h} .

# Ksat fix
cp /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/ksat/fix_state_change_ksat.{cpp,h} .

# Verify all files copied
ls -lh fix_state_change*.{cpp,h}
```

You should see 6 files:
- `fix_state_change_dimer.cpp` and `fix_state_change_dimer.h`
- `fix_state_change_dimer_ksat.cpp` and `fix_state_change_dimer_ksat.h`
- `fix_state_change_octahedron.cpp` and `fix_state_change_octahedron.h`
- `fix_state_change_ksat.cpp` and `fix_state_change_ksat.h`

### Step 5: Clean Old Build

Remove previous build artifacts:

```bash
rm -rf Obj_mpi
rm -f liblammps_mpi.a lmp_mpi
```

### Step 6: Build LAMMPS

Build with MPI support (takes 10-15 minutes):

```bash
make mpi -j 4
```

The `-j 4` flag uses 4 parallel jobs for faster compilation. Adjust based on available cores.

### Step 7: Verify Build

Check that the binary was created and all fixes are registered:

```bash
# Check binary exists
ls -lh lmp_mpi

# Verify all state change fixes are available
./lmp_mpi -help 2>&1 | grep "state/change"
```

You should see output like:
```
state/change
state/change/ksat
state/change/octahedron
```

## Alternative: Using the Script

A helper script is available at:
```
/work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/rebuild_manual.sh
```

**Important:** Still run it on an interactive node:

```bash
# Get interactive node
ssh gpuA

# Run the script
bash /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/rebuild_manual.sh
```

## Troubleshooting

### Error: mpicxx not found

If `mpicxx` is not found after setting PATH:
1. Verify you're on a compute node (not login node)
2. Try loading modules: `module load gcc-native/13.2 cray-mpich/8.1.32`
3. Check if mpicxx exists: `ls -la /opt/cray/pe/mpich/8.1.32/ofi/gnu/11.2/bin/mpicxx`

### Error: Conda compiler conflicts

If you get PIE object errors or linker issues:
1. Make sure conda is deactivated: `conda deactivate`
2. Check environment: `echo $CONDA_DEFAULT_ENV` (should be empty)
3. Use system compiler by ensuring conda bin is not in PATH

### Build takes too long

The build typically takes 10-15 minutes. If it's much longer:
- Check system load: `top` or `htop`
- Reduce parallel jobs: `make mpi -j 2` instead of `-j 4`

## Adding a New Fix

When adding a new state change fix:

1. Copy the new fix files to LAMMPS src directory
2. The Makefile should auto-detect `.cpp` files (if using wildcard)
3. If not auto-detected, add to Makefile manually
4. Rebuild following steps 5-7 above

## Verification Checklist

After rebuilding, verify:
- [ ] `lmp_mpi` binary exists and is executable
- [ ] All three fixes appear in `./lmp_mpi -help` output
- [ ] Binary size is reasonable (~100-150MB)
- [ ] Simulations using the fixes run without errors

## Notes

- The LAMMPS binary is at: `/work/nvme/bewl/lguttieres/lammps_build/lammps/src/lmp_mpi`
- All simulation scripts should use this custom binary (not the module version)
- Rebuild is only needed when adding new fixes or updating existing fix code
- The build includes all three fixes - you don't need separate builds for each

