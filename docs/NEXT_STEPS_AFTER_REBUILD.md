# Next Steps After Rebuilding LAMMPS

## Step 1: Verify the Build

```bash
# Check that LAMMPS binary exists
ls -lh /work/nvme/bewl/lguttieres/lammps_build/lammps/src/lmp_mpi

# Verify all fixes are available
/work/nvme/bewl/lguttieres/lammps_build/lammps/src/lmp_mpi -help 2>&1 | grep "state/change"
```

You should see:
- `fix state/change` (dimer - original)
- `fix state/change/ksat` (ksat fix)
- `fix state/change/octahedron` (octahedron fix)

## Step 2: Test with a Short Simulation

### For Ksat:

```bash
cd /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/ksat

# Run a quick test (modify submit script for shorter run first, or run locally)
# Quick local test:
cd ksat_simulation_cpp
/work/nvme/bewl/lguttieres/lammps_build/lammps/src/lmp_mpi -in in.ksat_monomers -log test.log
```

### For Octahedron:

```bash
cd /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/octahedron

# Quick local test:
cd octahedron_simulation_cpp
/work/nvme/bewl/lguttieres/lammps_build/lammps/src/lmp_mpi -in in.octahedron_monomers -log test.log
```

## Step 3: Submit Full Simulations

### For Ksat:

```bash
cd /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/ksat
sbatch submit_ksat.slurm
```

### For Octahedron:

```bash
cd /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/octahedron
sbatch submit_octahedron.slurm
```

## Step 4: Monitor the Simulations

```bash
# Check job status
squeue -u $USER

# Watch ksat output
tail -f /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/ksat/slurm_ksat-*.out

# Watch octahedron output
tail -f /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/octahedron/slurm_octahedron-*.out

# Check for errors
grep -i "error\|segmentation\|lost atoms" slurm_ksat-*.err
grep -i "error\|segmentation\|lost atoms" slurm_octahedron-*.err
```

## Step 5: What to Look For

### ✅ Success Indicators:
- No "Segmentation fault" errors
- No "Lost atoms" errors
- Simulation runs without crashing
- State changes are occurring (check thermo output or trajectory)
- Atom count remains constant

### ❌ Failure Indicators:
- Segmentation faults (should be fixed now!)
- "Lost atoms" errors
- Simulation crashes immediately
- Atom count changes

## Step 6: If There Are Still Issues

1. **Check the error logs** in `.err` files
2. **Verify the fixes are in the binary**:
   ```bash
   /work/nvme/bewl/lguttieres/lammps_build/lammps/src/lmp_mpi -help 2>&1 | grep state/change
   ```
3. **Check compilation warnings** - look at the build output for any warnings
4. **Compare with working dimer case** - the fixes should now work the same way

## Quick Verification Commands

```bash
# All-in-one verification
cd /work/nvme/bewl/lguttieres/lammps_build/lammps/src
echo "=== Checking LAMMPS binary ==="
[ -f lmp_mpi ] && echo "✅ Binary exists" || echo "❌ Binary missing"

echo ""
echo "=== Checking available fixes ==="
./lmp_mpi -help 2>&1 | grep "state/change" || echo "❌ No state/change fixes found"

echo ""
echo "=== Binary size ==="
ls -lh lmp_mpi
```

