# How to Run the Installation

## Quick Command

From the `fix_state_change` directory, run:

```bash
cd /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/fix_state_change
bash install_fix.sh
```

## What Happens

The script will:
1. ✅ Create `~/lammps_build` directory
2. ✅ Download LAMMPS source code (if not already present)
3. ✅ Copy `fix_state_change.h` and `fix_state_change.cpp` to LAMMPS src/
4. ✅ Add the fix to `Makefile.mpi`
5. ✅ Enable required packages (rigid, molecule)
6. ✅ Build LAMMPS (this takes ~10-15 minutes)
7. ✅ Verify the fix is installed

## Expected Output

You should see messages like:
```
=== Installing fix_state_change in LAMMPS ===
Step 1: Creating build directory...
Step 2: Downloading LAMMPS source...
Step 3: Copying fix files...
✅ Copied fix files to src/
Step 4: Adding fix to Makefile...
✅ Added fix_state_change.cpp to Makefile.mpi
Step 5: Enabling required LAMMPS packages...
Step 6: Building LAMMPS...
   (This may take several minutes...)
[compilation output...]
✅✅✅ LAMMPS built successfully! ✅✅✅
```

## After Installation

Your LAMMPS will be at:
```
~/lammps_build/lammps/src/lmp_mpi
```

Verify with:
```bash
~/lammps_build/lammps/src/lmp_mpi -help | grep state/change
```

Should output: `fix state/change`

## If Something Goes Wrong

1. **Build fails**: Check error messages. Common issues:
   - Missing modules: `module load gcc/11.4.0 openmpi/4.1.6`
   - Makefile issue: Manually check that `fix_state_change.cpp` is in `Makefile.mpi`

2. **Can't find git**: The script needs git to download LAMMPS. Install it or download manually.

3. **Permission errors**: Make sure you have write access to `~/lammps_build`

## Alternative: Manual Installation

If the script doesn't work, follow `INSTALL_STEP_BY_STEP.md` for manual instructions.

