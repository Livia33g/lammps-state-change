# Quick Start: Install fix_state_change

## The Fastest Way

Just run the installation script:

```bash
cd /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/fix_state_change
bash install_fix.sh
```

This will:
1. Download LAMMPS source (if needed)
2. Copy the fix files
3. Add fix to Makefile
4. Build LAMMPS with the fix
5. Verify installation

**Time required:** ~10-15 minutes (mostly compilation)

## After Installation

Your custom LAMMPS will be at:
```
~/lammps_build/lammps/src/lmp_mpi
```

### Update Your SLURM Script

Edit `submit_continuous.slurm` (or create a new one):

```bash
# Comment out or remove:
# module load lammps/2023.08.cuda.s11

# Add before running LAMMPS:
export PATH=$HOME/lammps_build/lammps/src:$PATH

# Change the run command:
srun lmp_mpi -in in.rigid_patchy_monomers | tee lammps_stdout.log
```

### Or Use Full Path

```bash
srun ~/lammps_build/lammps/src/lmp_mpi -in in.rigid_patchy_monomers | tee lammps_stdout.log
```

## Verify It Works

```bash
~/lammps_build/lammps/src/lmp_mpi -help | grep state/change
```

Should output:
```
fix state/change
```

## Test with Your Simulation

```bash
cd /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change
python3 generate_change_cpp.py
cd rigid_patchy_simulation_cpp
~/lammps_build/lammps/src/lmp_mpi -in in.rigid_patchy_monomers
```

## Troubleshooting

**Script fails at "Downloading LAMMPS"**
- Check internet connection
- Or manually download and extract LAMMPS to `~/lammps_build/lammps`

**Build fails with "module not found"**
- Try: `module avail` to see available modules
- Or build without modules (uses system compiler)

**"fix_state_change.cpp not found" in Makefile**
- Manually edit `~/lammps_build/lammps/src/Makefile.mpi`
- Find the list of `.cpp` files and add `fix_state_change.cpp \`

For more details, see `INSTALL_STEP_BY_STEP.md`

