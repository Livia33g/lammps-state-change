# Managing Multiple State Change Fixes in LAMMPS

## Current Architecture

**You have ONE LAMMPS build with MULTIPLE fixes compiled into it.**

### How It Works

1. **All fixes coexist**: Each fix has a unique class name and registers with LAMMPS using a unique identifier.
2. **No conflicts**: Different fixes are accessed by different names in LAMMPS input scripts.
3. **Single build**: All fixes are compiled into the same binary at `/work/nvme/bewl/lguttieres/lammps_build/lammps/src/lmp_mpi`

### Currently Available Fixes

Your LAMMPS build currently has:
- `fix state/change` → For dimer/rigid patchy monomer simulations
- `fix state/change/octahedron` → For octahedron monomer simulations

Verify with:
```bash
/work/nvme/bewl/lguttieres/lammps_build/lammps/src/lmp_mpi -help | grep "state/change"
```

### Why This Approach?

✅ **Advantages:**
- One build to maintain instead of many
- Easy to add new fixes (just add files and rebuild once)
- All fixes available in one binary
- No need to switch between different LAMMPS installations

❌ **Alternative (NOT recommended):**
- Separate LAMMPS build for each problem type
- Would require maintaining multiple installations
- More complex SLURM scripts (would need to select correct build)
- More disk space

## Adding a New State Change Fix

When you create a new simulation with different state change rules:

### Step 1: Create Your Fix Files

Create two files with a unique name:
- `fix_state_change_newproblem.h`
- `fix_state_change_newproblem.cpp`

**Important**: The class name and registration must be unique. In your `.cpp` file:
```cpp
FixStateChangeNewProblem::FixStateChangeNewProblem(...)
{
    ...
}

// Register with LAMMPS
FixStateChangeNewProblem::FixStateChangeNewProblem(LAMMPS *lmp, int narg, char **arg) 
    : Fix(lmp, narg, arg)
{
    ...
}
```

And in the constructor, register it with a unique name:
```cpp
// In the constructor or init()
// This registers as "state/change/newproblem"
```

### Step 2: Add to LAMMPS Build

**Option A: Use the helper script (recommended)**
```bash
cd /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change
./add_new_fix.sh newproblem \
    fix_state_change_newproblem.h \
    fix_state_change_newproblem.cpp
```

**Option B: Manual steps**
```bash
# 1. Copy files to LAMMPS src
cd /work/nvme/bewl/lguttieres/lammps_build/lammps/src
cp /path/to/fix_state_change_newproblem.h .
cp /path/to/fix_state_change_newproblem.cpp .

# 2. Add to Makefile
# Edit Makefile and add: fix_state_change_newproblem.cpp \
# (after other fix_state_change lines)

# 3. Rebuild
make mpi -j 4

# 4. Verify
./lmp_mpi -help | grep "state/change"
```

### Step 3: Use in Your Simulation

In your LAMMPS input script, use the new fix:
```lammps
fix 1 all state/change/newproblem <arguments>
```

## Workflow for New Problems

1. **Create fix files** in your problem directory (e.g., `newproblem/fix_state_change_newproblem.*`)
2. **Add to LAMMPS** using the helper script
3. **Create generation script** (e.g., `generate_newproblem_cpp.py`)
4. **Create SLURM script** that uses `fix state/change/newproblem`

## Directory Structure

```
state_change/
├── fix_state_change/              # Original dimer fix
│   ├── fix_state_change.h
│   └── fix_state_change.cpp
├── octahedron/                    # Octahedron fix
│   ├── fix_state_change_octahedron.h
│   └── fix_state_change_octahedron.cpp
├── newproblem/                    # Your next fix
│   ├── fix_state_change_newproblem.h
│   ├── fix_state_change_newproblem.cpp
│   ├── generate_newproblem_cpp.py
│   └── submit_newproblem.slurm
└── add_new_fix.sh                 # Helper script
```

## FAQ

**Q: Do I need to rebuild LAMMPS every time I modify a fix?**  
A: Yes, if you change the C++ code of a fix, you need to rebuild LAMMPS. But you can do incremental builds: `make mpi -j 4` will only recompile changed files.

**Q: Can I have fixes with the same logic but different parameters?**  
A: Usually no - each fix should have unique logic. If you just need different parameters, pass them as arguments to the same fix.

**Q: What if I want to test multiple versions of the same fix?**  
A: Create separate fixes with different names (e.g., `state/change/octahedron_v1`, `state/change/octahedron_v2`). This allows you to compare behavior.

**Q: How do I know which fixes are available?**  
A: Run: `/work/nvme/bewl/lguttieres/lammps_build/lammps/src/lmp_mpi -help | grep "state/change"`

**Q: What if I accidentally break a fix during development?**  
A: LAMMPS will fail to compile or the fix will be unavailable. Revert your changes and rebuild. The other fixes will still work once compilation succeeds.

## Current Build Location

- **Binary**: `/work/nvme/bewl/lguttieres/lammps_build/lammps/src/lmp_mpi`
- **Source**: `/work/nvme/bewl/lguttieres/lammps_build/lammps/src/`
- **Fix files**: Individual fix files are in each problem directory, with copies in LAMMPS src/

## Summary

✅ **You do NOT need separate LAMMPS builds for each problem**  
✅ **All fixes live in ONE build and are accessed by name**  
✅ **Adding a new fix is as simple as: copy files → add to Makefile → rebuild**  
✅ **Use the helper script to automate this process**

This is the standard LAMMPS way - one build, many fixes!

