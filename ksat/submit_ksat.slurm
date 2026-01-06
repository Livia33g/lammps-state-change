#!/bin/bash
#SBATCH --job-name=lammps_ksat
#SBATCH --account=bewl-delta-gpu
#SBATCH --partition=gpuA100x4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --output=slurm_ksat-%j.out
#SBATCH --error=slurm_ksat-%j.err

# Fail fast on any error or undefined variable
set -euo pipefail

echo "=== Starting LAMMPS with C++ fix_state_change/ksat ==="
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
module load libfabric/1.22.0 2>/dev/null || echo "⚠️  libfabric module not available"

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
SIM_DIR="ksat_simulation_cpp"
FORCE_REGENERATE="${FORCE_REGENERATE:-0}"  # Set to 1 to force regeneration with new parameters

if [[ ! -d "$SIM_DIR" ]] || [[ "$FORCE_REGENERATE" == "1" ]]; then
    if [[ -d "$SIM_DIR" ]]; then
        echo "--- Force regeneration requested. Removing existing '$SIM_DIR' directory... ---"
        rm -rf "$SIM_DIR"
    fi
    echo "--- Generating input files for ksat monomers with C++ fix_state_change... ---"
    python -c "
import sys
import os
# Add parent directory to path to access the script
sys.path.insert(0, '$SLURM_SUBMIT_DIR')
from generate_ksat_cpp import create_lammps_ksat_script_cpp

create_lammps_ksat_script_cpp(
    num_monomers_A=200,  # Type A: patches type 2 (never change)
    num_monomers_B=400,  # Type B: patches type 3 (can change to 4)
    box_size=None,  # Auto-calculate from density
    patch_coordination_cutoff=0.28,  # Coordination cutoff for state change detection
    state_change_probability=0.7,  # 70% probability when conditions met
    timesteps=2700000,  # Total simulation steps
    thermo_freq=20000,  # Thermodynamic output frequency
    dump_freq=10000,  # Trajectory dump frequency
    output_dir='$SIM_DIR',
    seed=12345,
    state_change_freq=100,  # Check state changes every 100 steps
    cooldown_steps=1000,  # Cooldown period between state changes
    timestep=0.0002,  # Simulation timestep
    temperature_initial=1.0,  # Initial temperature
    temperature_final=0.4,  # Final temperature
    steps_T1=900000,  # Steps at T=1.0
    steps_T2=1800000  # Steps at T=0.4
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
if [[ ! -f "in.ksat_monomers" || ! -f "data.ksat_monomers" ]]; then
    echo "❌ ERROR: Missing LAMMPS input files in $(pwd). Exiting."
    exit 1
fi

# Verify custom LAMMPS is available (in work directory)
LAMMPS_EXE="/work/nvme/bewl/lguttieres/lammps_build/lammps/src/lmp_mpi"
if [[ ! -f "$LAMMPS_EXE" ]]; then
    echo "❌ ERROR: Custom LAMMPS not found at $LAMMPS_EXE"
    echo "   Please rebuild LAMMPS first with fix_state_change_ksat:"
    echo "   cd /work/nvme/bewl/lguttieres/lammps_build/lammps/src"
    echo "   make mpi -j 4"
    exit 1
fi

# Check if hybrid/overlay pair style is used (not compatible with GPU)
if grep -q "pair_style.*hybrid/overlay" in.ksat_monomers; then
  echo "Detected hybrid/overlay pair style - GPU acceleration not supported, removing GPU commands if present..."
  sed -i '/^package gpu/d' in.ksat_monomers
  sed -i '/^suffix gpu/d' in.ksat_monomers
  sed -i '/^$/N;/^\n$/d' in.ksat_monomers
fi

echo "Running LAMMPS with C++ fix_state_change/ksat"
echo "This version uses the custom fix - NO unfix/refix cycles needed!"
echo "Ksat structure: 1 body + 3 patches = 4 particles per monomer"
echo "Monomer A: patches type 2 (never change)"
echo "Monomer B: patches type 3 → 4 when attached to type 2"

# Run custom LAMMPS (using work directory path)
# Unset conflicting SLURM memory variables if they exist
unset SLURM_MEM_PER_CPU SLURM_MEM_PER_GPU SLURM_MEM_PER_NODE 2>/dev/null || true
srun "$LAMMPS_EXE" -in in.ksat_monomers | tee lammps_stdout.log

echo "✅ LAMMPS job completed successfully."

# ------------------------------------------------------------
# 4. Deactivate Environment
# ------------------------------------------------------------
conda deactivate
echo "Job finished."

