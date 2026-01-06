# Step-by-Step Installation Guide for fix_state_change

## Your Situation
You're using a pre-built LAMMPS module (`lammps/2023.08.cuda.s11`). To add a custom fix, you need to build LAMMPS from source.

## Option 1: Build LAMMPS in Your User Directory (Recommended)

### Step 1: Download LAMMPS Source

```bash
cd ~
mkdir -p lammps_build
cd lammps_build

# Download LAMMPS (same version as your module: 2023-08-17)
wget https://github.com/lammps/lammps/archive/refs/tags/stable_23Aug2023.tar.gz
tar -xzf stable_23Aug2023.tar.gz
cd lammps-stable_23Aug2023
```

Or use git:
```bash
git clone -b stable_23Aug2023 https://github.com/lammps/lammps.git
cd lammps
```

### Step 2: Copy the Custom Fix

```bash
# Copy the fix files to LAMMPS src directory
cp /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/fix_state_change/fix_state_change.h src/
cp /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/fix_state_change/fix_state_change.cpp src/
```

### Step 3: Add Fix to Makefile

```bash
cd src

# Edit Makefile.mpi (or Makefile.serial if you want serial version)
# Find the section with fix_*.cpp files and add:
#   fix_state_change.cpp \

# Quick way to add it:
sed -i '/^fix_rigid_nvt\.cpp \\/a\        fix_state_change.cpp \\' Makefile.mpi
```

Or manually edit `src/Makefile.mpi` and add `fix_state_change.cpp` to the list of source files.

### Step 4: Build LAMMPS

```bash
cd ~/lammps_build/lammps-stable_23Aug2023/src

# Load required modules (same as your SLURM script)
module purge
module load gcc/11.4.0
module load openmpi/4.1.6
module load cuda/12.2.1  # If you want GPU support

# Build LAMMPS
make yes-rigid
make yes-molecule
make mpi -j 4

# This will create lmp_mpi in the src directory
```

### Step 5: Test Installation

```bash
cd ~/lammps_build/lammps-stable_23Aug2023/src
./lmp_mpi -help | grep state/change
```

You should see:
```
fix state/change
```

### Step 6: Update Your SLURM Script

Edit your SLURM script to use your custom LAMMPS:

```bash
# Instead of:
# module load lammps/2023.08.cuda.s11
# srun lmp -in in.rigid_patchy_monomers

# Use:
# srun ~/lammps_build/lammps-stable_23Aug2023/src/lmp_mpi -in in.rigid_patchy_monomers
```

Or add to PATH:
```bash
export PATH=~/lammps_build/lammps-stable_23Aug2023/src:$PATH
srun lmp_mpi -in in.rigid_patchy_monomers
```

## Option 2: Quick Test Build (Minimal)

If you just want to test quickly:

```bash
cd ~
mkdir lammps_test
cd lammps_test
git clone -b stable_23Aug2023 https://github.com/lammps/lammps.git
cd lammps/src

# Copy fix
cp /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/fix_state_change/fix_state_change.* .

# Add to Makefile
echo "fix_state_change.cpp \\" >> Makefile.mpi.tmp
# Then manually merge into Makefile.mpi

# Build
module load gcc/11.4.0 openmpi/4.1.6
make yes-rigid yes-molecule
make mpi -j 4
```

## Option 3: Contact System Administrator

If you don't have build tools or want the fix available system-wide:

1. Contact your HPC administrator
2. Provide them with:
   - The fix files (`fix_state_change.h` and `fix_state_change.cpp`)
   - Request to add it to the LAMMPS build
   - They can rebuild the module with your fix included

## Troubleshooting

### Compilation Errors

**Error: "utils.h: No such file"**
- Make sure you're using LAMMPS 2023-08-17 or later
- Check that you're in the `src/` directory when building

**Error: "MPI functions not found"**
- Make sure `module load openmpi/4.1.6` is loaded
- Check that MPI is in your PATH: `which mpicc`

**Error: "fix_state_change.cpp: undefined reference"**
- Make sure `fix_state_change.cpp` is listed in `Makefile.mpi`
- Try `make clean` then `make mpi` again

### Runtime Errors

**"Unknown fix style: state/change"**
- The fix wasn't compiled into LAMMPS
- Rebuild LAMMPS and verify with `lmp_mpi -help | grep state/change`

**"Fix state/change compute cA not found"**
- Make sure you define the computes before the fix in your input script:
  ```
  compute cA all coord/atom cutoff 2.5 group patches_A
  compute cB all coord/atom cutoff 2.5 group patches_B
  ```

## Quick Verification Script

Create `test_fix.sh`:

```bash
#!/bin/bash
module load gcc/11.4.0 openmpi/4.1.6

cd ~/lammps_build/lammps-stable_23Aug2023/src

if ./lmp_mpi -help | grep -q "state/change"; then
    echo "✅ Fix installed successfully!"
else
    echo "❌ Fix not found. Rebuild LAMMPS."
fi
```

Run: `bash test_fix.sh`

