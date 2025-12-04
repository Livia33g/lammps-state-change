# Installing fix_state_change/octahedron in LAMMPS

The new `fix_state_change/octahedron` needs to be compiled into your LAMMPS build before use.

## Quick Installation

If you already have a custom LAMMPS build at `/work/nvme/bewl/lguttieres/lammps_build/lammps/src`:

```bash
# 1. Copy the fix files to LAMMPS src directory
cd /work/nvme/bewl/lguttieres/lammps_build/lammps/src
cp /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/octahedron/fix_state_change_octahedron.h .
cp /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/octahedron/fix_state_change_octahedron.cpp .

# 2. Add to Makefile.mpi
# Find the section with fix_*.cpp files and add:
#   fix_state_change_octahedron.cpp \
# 
# Quick way (adds after fix_state_change.cpp if it exists, otherwise after fix_rigid_nvt.cpp):
if grep -q "fix_state_change.cpp" Makefile.mpi; then
  sed -i '/^fix_state_change\.cpp \\/a\        fix_state_change_octahedron.cpp \\' Makefile.mpi
else
  sed -i '/^fix_rigid_nvt\.cpp \\/a\        fix_state_change_octahedron.cpp \\' Makefile.mpi
fi

# 3. Rebuild LAMMPS
make mpi -j 4

# 4. Verify installation
./lmp_mpi -help | grep "state/change"
# Should show both "state/change" and "state/change/octahedron"
```

## Verification

After installation, verify the fix is available:

```bash
/work/nvme/bewl/lguttieres/lammps_build/lammps/src/lmp_mpi -help | grep state/change
```

You should see:
```
fix state/change
fix state/change/octahedron
```

## Troubleshooting

### "Unknown fix style: state/change/octahedron"
- The fix files weren't copied correctly
- Makefile wasn't updated
- LAMMPS wasn't rebuilt after adding the fix

### Compilation errors
- Make sure you have the required LAMMPS packages enabled
- Check that the header file is included correctly

