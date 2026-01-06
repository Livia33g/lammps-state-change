#!/bin/bash
#SBATCH --job-name=lammps_debugging
#SBATCH --account=bewl-delta-gpu
#SBATCH --partition=gpuA100x4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --output=debugging_slurm-%j.out
#SBATCH --error=debugging_slurm-%j.err

# Fail fast on any error or undefined variable
set -euo pipefail

echo "=== Starting LAMMPS DEBUGGING run with updated potentials ==="
echo "=== D0=3.0, kT=1.0, rep_epsilon=2.0, mass=1.0, strong damping (200.0) ==="
cd "$SLURM_SUBMIT_DIR"
echo "Working directory: $(pwd)"

# ------------------------------------------------------------
# 1. Load Environment (NO LAMMPS MODULE - using custom build)
# ------------------------------------------------------------
echo "--- Loading environment ---"
module purge
# Don't load lammps module - we're using custom build

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
# 2. Create debugging directory and generate inputs
# ------------------------------------------------------------
DEBUG_DIR="debugging"
mkdir -p "$DEBUG_DIR"
cd "$DEBUG_DIR"
echo "--- Working in debugging directory: $(pwd) ---"

SIM_DIR="rigid_patchy_simulation_cpp"
FORCE_REGENERATE="${FORCE_REGENERATE:-1}"  # Set to 1 to force regeneration with new parameters

if [[ ! -d "$SIM_DIR" ]] || [[ "$FORCE_REGENERATE" == "1" ]]; then
    if [[ -d "$SIM_DIR" ]]; then
        echo "--- Force regeneration requested. Removing existing '$SIM_DIR' directory... ---"
        rm -rf "$SIM_DIR"
    fi
    echo "--- Generating input files with DEBUGGING parameters... ---"
    echo "   - D0 = 3.0 (Morse epsilon)"
    echo "   - kT = 1.0 (temperature)"
    echo "   - 3 main body spheres, R = 1.0, one patch per monomer at +x face (r = 0.33)"
    echo "   - rep_epsilon = 2.0 (short-range body-body repulsion)"
    echo "   - mass = 1.0 main body / ~0 mass patches"
    echo "   - damping = 200.0 (strong damping to reduce rotation)"
    echo "   - No body-patch interactions (only same-species interactions)"
    python3 -c "
import sys
import os
# Add parent directory to path to import generate_change_cpp
sys.path.insert(0, '..')
from generate_change_cpp import create_lammps_rigid_patchy_monomers_script_cpp

create_lammps_rigid_patchy_monomers_script_cpp(
    num_monomers=40,
    box_size=None,
    rigid_sphere_radius=1.0,
    patch_radius=0.33,
    patch_offset_dist=0.75,
    patch_interaction_cutoff_morse=0.33,
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
    morse_D0_22=3.0,
    morse_D0_33=3.0,
    rep_epsilon=2.0,
    temperature=1.0
)
print('✅ Input generation complete in debugging/$SIM_DIR/')
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

echo "Running LAMMPS with DEBUGGING parameters:"
echo "  - D0 = 3.0 (morse_D0_22 and morse_D0_33)"
echo "  - kT = 1.0 (temperature)"
echo "  - rep_epsilon = 20.0 (body-body repulsion)"
echo "  - mass = 1.0 (main body mass)"
echo "  - damping = 200.0 (strong damping to reduce rotation)"
echo "  - No body-patch interactions"
echo "  - Only same-type patch interactions (2-2 and 3-3)"
echo "This version uses the custom fix - NO unfix/refix cycles needed!"

# Run custom LAMMPS (using work directory path)
# Unset conflicting SLURM memory variables if they exist
unset SLURM_MEM_PER_CPU SLURM_MEM_PER_GPU SLURM_MEM_PER_NODE 2>/dev/null || true
srun "$LAMMPS_EXE" -in in.rigid_patchy_monomers | tee lammps_stdout.log

echo "✅ LAMMPS job completed successfully."
echo "All output files are in: debugging/$SIM_DIR/"

# ------------------------------------------------------------
# 4. Deactivate Environment
# ------------------------------------------------------------
conda deactivate
echo "Job finished. Check debugging/rigid_patchy_simulation_cpp/ for results."
