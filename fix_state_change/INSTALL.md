# Installation Guide for fix_state_change

## Quick Start

1. **Locate your LAMMPS source directory**
   ```bash
   # Typically something like:
   # /path/to/lammps/src
   # or
   # $HOME/lammps/src
   ```

2. **Copy the fix files**
   ```bash
   cd /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/fix_state_change
   cp fix_state_change.h /path/to/lammps/src/
   cp fix_state_change.cpp /path/to/lammps/src/
   ```

3. **Edit LAMMPS Makefile**
   
   Edit `/path/to/lammps/src/Makefile` and find the section that lists source files. Add:
   ```
   fix_state_change.cpp
   ```
   
   For example, if you see lines like:
   ```
   OBJ = \
   ...
   fix_rigid_nvt.o \
   ...
   ```
   
   Add:
   ```
   fix_state_change.o \
   ```

4. **Rebuild LAMMPS**
   ```bash
   cd /path/to/lammps/src
   make clean
   make mpi
   ```

5. **Verify installation**
   ```bash
   lmp_mpi -help | grep state/change
   ```
   
   You should see:
   ```
   fix state/change
   ```

## Alternative: Build as Package

If you prefer to keep the fix separate:

1. Create a package directory:
   ```bash
   mkdir -p /path/to/lammps/src/USER-STATE_CHANGE
   cp fix_state_change.* /path/to/lammps/src/USER-STATE_CHANGE/
   ```

2. Edit `/path/to/lammps/src/Makefile` to include the package:
   ```makefile
   # Add to the list of packages
   PACKAGES = ... USER-STATE_CHANGE ...
   ```

3. Rebuild LAMMPS

## Troubleshooting

### Compilation Errors

- **Missing includes**: Make sure you're using a recent version of LAMMPS (2020 or later)
- **MPI errors**: Ensure MPI is properly configured in your LAMMPS build
- **Utils errors**: The `utils.h` include requires LAMMPS 2018 or later

### Runtime Errors

- **"Fix state/change compute cA not found"**: Make sure you define the computes before the fix:
  ```
  compute cA all coord/atom cutoff 2.5 group patches_A
  compute cB all coord/atom cutoff 2.5 group patches_B
  ```

- **"Lost atoms"**: This should NOT happen with the C++ fix. If it does, check:
  - Atom types are correctly initialized (2 or 3 for patches)
  - Groups are correctly defined
  - Cutoff matches your Morse potential cutoff

## Testing

Create a simple test script `test_state_change.in`:

```
units lj
atom_style molecular
read_data test.data

group patches_A type 2
group patches_B type 3

compute cA all coord/atom cutoff 2.5 group patches_A
compute cB all coord/atom cutoff 2.5 group patches_B

fix rigid_nvt all rigid/nvt molecule temp 1.0 1.0 0.01
fix state_change all state/change 1 100 0.7 2.5 patches_A patches_B

run 1000
```

If this runs without errors, the fix is installed correctly!

