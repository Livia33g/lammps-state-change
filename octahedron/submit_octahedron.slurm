#!/bin/bash
#SBATCH --job-name=lammps_octahedron
#SBATCH --account=bewl-delta-gpu
#SBATCH --partition=gpuA100x4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --output=slurm_octahedron-%j.out
#SBATCH --error=slurm_octahedron-%j.err

# Fail fast on any error or undefined variable
set -euo pipefail

echo "=== Starting LAMMPS with C++ fix_state_change for Octahedron Monomers ==="
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
# Load libfabric (required by LAMMPS binary that was built with MPI)
# CRITICAL: libfabric must load successfully or LAMMPS will fail with "libfabric.so.1: cannot open shared object file"
module load libfabric/1.22.0 || { echo "❌ ERROR: Failed to load libfabric/1.22.0"; exit 1; }

set +u
source /u/lguttieres/miniforge3/bin/activate
conda activate py-env
set -u
echo "Activated conda environment"

# Set path to custom LAMMPS (in work directory for better I/O performance)
export PATH=/work/nvme/bewl/lguttieres/lammps_build/lammps/src:$PATH
echo "Using custom LAMMPS at: /work/nvme/bewl/lguttieres/lammps_build/lammps/src/lmp_mpi"

# ------------------------------------------------------------
# 2. Generate Inputs if Missing
# ------------------------------------------------------------
SIM_DIR="octahedron_simulation_cpp"
FORCE_REGENERATE="${FORCE_REGENERATE:-0}"  # Set to 1 to force regeneration with new parameters

if [[ ! -d "$SIM_DIR" ]] || [[ "$FORCE_REGENERATE" == "1" ]]; then
    if [[ -d "$SIM_DIR" ]]; then
        echo "--- Force regeneration requested. Removing existing '$SIM_DIR' directory... ---"
        rm -rf "$SIM_DIR"
    fi
    echo "--- Generating input files for octahedron monomers with C++ fix_state_change... ---"
    python -c "
import sys
import os
# Add parent directory to path to access the script
sys.path.insert(0, '$SLURM_SUBMIT_DIR/octahedron')
from generate_octahedron_cpp import create_lammps_octahedron_script_cpp

create_lammps_octahedron_script_cpp(
    num_monomers=50,  # Number of octahedron monomers (increased for better cluster formation)
    box_size=None,  # Auto-calculate from density
    vertex_radius=1.0,  # Body particle radius (matches reference geometry)
    patch_radius=0.5,  # Patch sphere radius (for visualization/planning)
    patch_coordination_cutoff=0.34,  # Tight cutoff for state change detection
    state_change_probability=0.25,  # INCREASED: 25% probability when conditions met (was 0.1, increased to allow more state changes while keeping anti-sync protection)
    timesteps=500000000,  # Total simulation steps
    thermo_freq=5000,  # Thermodynamic output frequency
    dump_freq=1000,  # Trajectory dump frequency (reduced from 10000 for smoother video - 10x more frames)
    output_dir='$SIM_DIR',
    seed=12345,
    state_change_freq=10,  # REDUCED: Check state changes every 10 steps for almost dynamic/asynchronous behavior (was 100, reduced for more responsive state changes)
    cooldown_steps=10000,  # Cooldown period between state changes (prevents rapid cascades, randomized 0.5x-1.5x in code)
    hysteresis_threshold=75,  # REDUCED: Require 75 consecutive steps of contact before state change (was 200, reduced to allow more state changes while still filtering noise)
    timestep=0.0001,  # Simulation timestep (further reduced for stability, prevents NaNs)
    morse_D0_22=10.0,  # Morse well depth for type 2-2 interactions
    morse_D0_33=10.0,  # Morse well depth for type 3-3 interactions
    morse_D0_23=0.0,  # No cross-interaction between 2-3
    rep_epsilon=5000.0,  # Repulsion strength for body-body (INCREASED to match JAX soft_eps=10000 style, prevents unwanted large clusters)
    temperature=0.3,  # Simulation temperature (reduced for stability)
    morse_cutoff=2.5  # Morse potential cutoff (INCREASED for better capture radius, stability fix)
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
cd "$SLURM_SUBMIT_DIR"
if [ ! -d "$SIM_DIR" ]; then
    echo "❌ ERROR: Simulation directory '$SIM_DIR' not found in $(pwd)"
    exit 1
fi
cd "$SIM_DIR"

# Safety check for input files
if [[ ! -f "in.octahedron_monomers" || ! -f "data.octahedron_monomers" ]]; then
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
if grep -q "pair_style.*hybrid/overlay" in.octahedron_monomers; then
  echo "Detected hybrid/overlay pair style - GPU acceleration not supported, removing GPU commands if present..."
  sed -i '/^package gpu/d' in.octahedron_monomers
  sed -i '/^suffix gpu/d' in.octahedron_monomers
  sed -i '/^$/N;/^\n$/d' in.octahedron_monomers
fi

echo "Running LAMMPS with C++ fix_state_change/octahedron for Octahedron Monomers"
echo "This version uses the custom fix - NO unfix/refix cycles needed!"
echo "Octahedron structure: 1 body + 4 patches = 5 particles per monomer"

# Run custom LAMMPS (using work directory path)
# Unset conflicting SLURM memory variables if they exist
unset SLURM_MEM_PER_CPU SLURM_MEM_PER_GPU SLURM_MEM_PER_NODE 2>/dev/null || true
srun "$LAMMPS_EXE" -in in.octahedron_monomers | tee lammps_stdout.log

echo "✅ LAMMPS job completed successfully."

# ------------------------------------------------------------
# 4. Deactivate Environment
# ------------------------------------------------------------
conda deactivate
echo "Job finished."

