#!/bin/bash
#SBATCH --job-name=lammps_continuous_state
#SBATCH --account=bewl-delta-gpu
#SBATCH --partition=gpuA100x4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=slurm_continuous-%j.out
#SBATCH --error=slurm_continuous-%j.err

# Fail fast on any error or undefined variable
set -euo pipefail

echo "=== Starting LAMMPS continuous state change job (with cooldown) ==="
# The job starts in the directory where you run 'sbatch'
cd "$SLURM_SUBMIT_DIR"
echo "Working directory: $(pwd)"

# ------------------------------------------------------------
# 1. Load Environment
# ------------------------------------------------------------
echo "--- Loading environment ---"
module purge
module load lammps/2023.08.cuda.s11
echo "Loaded LAMMPS module"

set +u
source /u/lguttieres/miniforge3/bin/activate
conda activate py-env
set -u
echo "Activated conda environment"

# ------------------------------------------------------------
# 2. Generate Inputs if Missing
# ------------------------------------------------------------
SIM_DIR="rigid_patchy_simulation_continuous"
FORCE_REGENERATE="${FORCE_REGENERATE:-0}"  # Set to 1 to force regeneration

if [[ ! -d "$SIM_DIR" ]] || [[ "$FORCE_REGENERATE" == "1" ]]; then
    if [[ -d "$SIM_DIR" ]]; then
        echo "--- Force regeneration requested. Removing existing '$SIM_DIR' directory... ---"
        rm -rf "$SIM_DIR"
    fi
    echo "--- Generating input files with continuous state changes and cooldown... ---"
    # This script is in the same directory as generate_change.py
    # Generate directly to the new directory name
    python -c "
import sys
sys.path.insert(0, '.')
from generate_change import create_lammps_rigid_patchy_monomers_script

create_lammps_rigid_patchy_monomers_script(
    num_monomers=10,
    box_size=None,
    rigid_sphere_radius=1.0,
    patch_radius=0.33,
    patch_offset_dist=1.0,
    patch_interaction_cutoff_morse=2.5,
    state_change_probability=0.7,
    timesteps=500000000,
    thermo_freq=5000,
    dump_freq=50000,
    output_dir='$SIM_DIR',
    seed=12345,
    state_change_freq=1,
    cooldown_steps=100
)
"
    echo "✅ Input generation complete."
else
    echo "Found existing simulation directory '$SIM_DIR' — skipping generation."
    echo "  (Set FORCE_REGENERATE=1 to force regeneration)"
fi

# ------------------------------------------------------------
# 3. Run LAMMPS
# ------------------------------------------------------------
echo "--- Entering simulation directory and running LAMMPS ---"
# This is the crucial step: change into the subdirectory
cd "$SIM_DIR"

# Safety check for input files
if [[ ! -f "in.rigid_patchy_monomers" || ! -f "data.rigid_patchy_monomers" ]]; then
    echo "❌ ERROR: Missing LAMMPS input files in $(pwd). Exiting."
    exit 1
fi

# Check if hybrid/overlay pair style is used (not compatible with GPU)
if grep -q "pair_style.*hybrid/overlay" in.rigid_patchy_monomers; then
  echo "Detected hybrid/overlay pair style - GPU acceleration not supported, removing GPU commands if present..."
  # Remove GPU package commands if they exist
  sed -i '/^package gpu/d' in.rigid_patchy_monomers
  sed -i '/^suffix gpu/d' in.rigid_patchy_monomers
  # Remove empty lines at the start
  sed -i '/^$/N;/^\n$/d' in.rigid_patchy_monomers
elif ! grep -q "package gpu" in.rigid_patchy_monomers; then
  echo "Adding GPU package commands to LAMMPS input file..."
  # Use '1i' to insert at the first line
  sed -i '1i\
package gpu 1\
suffix gpu\
' in.rigid_patchy_monomers
fi

echo "Running LAMMPS in $(pwd)"
echo "This version uses continuous state change checks (every timestep) with cooldown mechanism"
# Run LAMMPS. The input script will find 'data.rigid_patchy_monomers' in the current directory.
srun lmp -in in.rigid_patchy_monomers | tee lammps_stdout.log

echo "✅ LAMMPS job completed successfully."

# ------------------------------------------------------------
# 4. Deactivate Environment
# ------------------------------------------------------------
conda deactivate
echo "Job finished."

