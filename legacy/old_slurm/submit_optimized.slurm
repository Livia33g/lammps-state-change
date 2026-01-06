#!/bin/bash
#SBATCH --job-name=lammps_optimized
#SBATCH --account=bewl-delta-gpu
#SBATCH --partition=gpuA100x4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --output=optimized_slurm-%j.out
#SBATCH --error=optimized_slurm-%j.err

# Fail fast on any error or undefined variable
set -euo pipefail

echo "=== Starting LAMMPS OPTIMIZED run with parameters from optimization ==="
echo "=== D0=5.137, epsilon_rep=500.0, kT=0.460, conc=0.001 ==="
cd "$SLURM_SUBMIT_DIR"
echo "Working directory: $(pwd)"

# ------------------------------------------------------------
# 1. Load Environment (NO LAMMPS MODULE - using custom build)
# ------------------------------------------------------------
echo "--- Loading environment ---"
module purge

# Load other required modules
module load gcc/11.4.0 2>/dev/null || echo "⚠️  gcc module not available"
module load openmpi/4.1.6 2>/dev/null || echo "⚠️  openmpi module not available"

set +u
source /u/lguttieres/miniforge3/bin/activate
conda activate py-env
set -u
echo "Activated conda environment"

# Set path to custom LAMMPS (in work directory for better I/O performance)
export PATH=/work/nvme/bewl/lguttieres/lammps_build/lammps/src:$PATH
echo "Using custom LAMMPS at: /work/nvme/bewl/lguttieres/lammps_build/lammps/src/lmp_mpi"

# ------------------------------------------------------------
# 2. Create optimized directory and generate inputs
# ------------------------------------------------------------
OPT_DIR="optimized"
mkdir -p "$OPT_DIR"
cd "$OPT_DIR"
echo "--- Working in optimized directory: $(pwd) ---"

SIM_DIR="rigid_patchy_simulation_cpp"
FORCE_REGENERATE="${FORCE_REGENERATE:-1}"  # Set to 1 to force regeneration with new parameters

# Optimized parameters from optimization results
OPT_D0=5.137350042125283
OPT_EPSILON_REP=500.0
OPT_KT=0.46043863918287997
OPT_CONC=0.001

# Calculate box size from concentration
# V = N / conc, box_size = V^(1/3)
NUM_MONOMERS=40
VOLUME=$(python3 -c "print($NUM_MONOMERS / $OPT_CONC)")
BOX_SIZE=$(python3 -c "import math; print(($NUM_MONOMERS / $OPT_CONC) ** (1.0/3.0))")
echo "--- Optimized parameters ---"
echo "   D0 (morse_D0_22 and morse_D0_33): $OPT_D0"
echo "   epsilon_rep: $OPT_EPSILON_REP"
echo "   kT: $OPT_KT"
echo "   concentration: $OPT_CONC"
echo "   num_monomers: $NUM_MONOMERS"
echo "   box_size: $BOX_SIZE"

if [[ ! -d "$SIM_DIR" ]] || [[ "$FORCE_REGENERATE" == "1" ]]; then
    if [[ -d "$SIM_DIR" ]]; then
        echo "--- Force regeneration requested. Removing existing '$SIM_DIR' directory... ---"
        rm -rf "$SIM_DIR"
    fi
    echo "--- Generating input files with OPTIMIZED parameters... ---"
    python3 -c "
import sys
import os
# Add parent directory to path to import generate_change_cpp
sys.path.insert(0, '..')
from generate_change_cpp import create_lammps_rigid_patchy_monomers_script_cpp

# Use optimized parameters directly in function call
create_lammps_rigid_patchy_monomers_script_cpp(
    num_monomers=$NUM_MONOMERS,
    box_size=$BOX_SIZE,  # Use calculated box size from concentration
    rigid_sphere_radius=2.0,
    patch_radius=0.33,
    patch_offset_dist=2.0,
    patch_interaction_cutoff_morse=0.66,
    patch_coordination_cutoff=0.34,
    state_change_probability=0.7,
    timesteps=500000000,
    thermo_freq=5000,
    dump_freq=1000,  # Save every 1k steps to capture state changes
    output_dir='$SIM_DIR',
    seed=12345,
    state_change_freq=100,  # Check every 100 steps
    cooldown_steps=1000,  # Cooldown to prevent rapid toggling
    timestep=0.0002,  # Timestep for stability
    morse_D0_22=$OPT_D0,  # Optimized D0 for type 2-2
    morse_D0_33=$OPT_D0,  # Optimized D0 for type 3-3
    rep_epsilon=$OPT_EPSILON_REP,  # Optimized repulsion epsilon
    temperature=$OPT_KT  # Optimized temperature
)
print('✅ Input generation complete in optimized/$SIM_DIR/')
"
    echo "✅ Input generation and parameter update complete."
else
    echo "Found existing simulation directory '$SIM_DIR' — skipping generation."
    echo "  (Set FORCE_REGENERATE=1 to force regeneration)"
fi

# ------------------------------------------------------------
# 3. Run LAMMPS
# ------------------------------------------------------------
echo "--- Entering simulation directory and running LAMMPS ---"
cd "$SIM_DIR"

# Safety check for input files
if [[ ! -f "in.rigid_patchy_monomers" || ! -f "data.rigid_patchy_monomers" ]]; then
    echo "❌ ERROR: Missing LAMMPS input files in $(pwd). Exiting."
    exit 1
fi

# Verify custom LAMMPS is available (in work directory)
LAMMPS_EXE="/work/nvme/bewl/lguttieres/lammps_build/lammps/src/lmp_mpi"
if [[ ! -f "$LAMMPS_EXE" ]]; then
    echo "❌ ERROR: Custom LAMMPS not found at $LAMMPS_EXE"
    echo "   Please rebuild LAMMPS first:"
    echo "   cd /work/nvme/bewl/lguttieres/lammps_build/lammps/src"
    echo "   make mpi -j 4"
    exit 1
fi

# Check if hybrid/overlay pair style is used (not compatible with GPU)
if grep -q "pair_style.*hybrid/overlay" in.rigid_patchy_monomers; then
  echo "Detected hybrid/overlay pair style - GPU acceleration not supported, removing GPU commands if present..."
  sed -i '/^package gpu/d' in.rigid_patchy_monomers
  sed -i '/^suffix gpu/d' in.rigid_patchy_monomers
  sed -i '/^$/N;/^\n$/d' in.rigid_patchy_monomers
fi

echo "Running LAMMPS with OPTIMIZED parameters:"
echo "  - D0 = $OPT_D0 (morse_D0_22 and morse_D0_33)"
echo "  - epsilon_rep = $OPT_EPSILON_REP"
echo "  - kT = $OPT_KT"
echo "  - concentration = $OPT_CONC"
echo "This version uses the custom fix - NO unfix/refix cycles needed!"

# Run custom LAMMPS (using work directory path)
# Unset conflicting SLURM memory variables if they exist
unset SLURM_MEM_PER_CPU SLURM_MEM_PER_GPU SLURM_MEM_PER_NODE 2>/dev/null || true
srun "$LAMMPS_EXE" -in in.rigid_patchy_monomers | tee lammps_stdout.log

echo "✅ LAMMPS job completed successfully."
echo "All output files are in: optimized/$SIM_DIR/"

# ------------------------------------------------------------
# 4. Deactivate Environment
# ------------------------------------------------------------
conda deactivate
echo "Job finished. Check optimized/rigid_patchy_simulation_cpp/ for results."

