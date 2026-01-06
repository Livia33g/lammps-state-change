#!/bin/bash
#SBATCH --job-name=rebuild_lammps_ksat
#SBATCH --account=bewl-delta-gpu
#SBATCH --partition=gpuA100x4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --output=rebuild_lammps_ksat-%j.out
#SBATCH --error=rebuild_lammps_ksat-%j.err
# Note: Requesting GPU to use GPU partition, but only using CPU for compilation

set -euo pipefail

echo "=== Rebuilding LAMMPS with fix_state_change_ksat ==="
echo "Job started at: $(date)"

LAMMPS_SRC="/work/nvme/bewl/lguttieres/lammps_build/lammps/src"
cd "$LAMMPS_SRC"

# Load modules
module purge
module load gcc/11.4.0 2>/dev/null || echo "⚠️  Using system gcc"
module load openmpi/4.1.6 2>/dev/null || echo "⚠️  Using system MPI"

# Verify fix files are present
if [[ ! -f "fix_state_change_ksat.cpp" ]] || [[ ! -f "fix_state_change_ksat.h" ]]; then
    echo "❌ ERROR: fix_state_change_ksat files not found!"
    exit 1
fi

echo "✅ Fix files found:"
ls -lh fix_state_change_ksat.*

# Clean old build (optional - remove if you want incremental build)
echo ""
echo "Cleaning old build artifacts..."
rm -rf Obj_mpi
rm -f liblammps_mpi.a lmp_mpi 2>/dev/null || true

# Build LAMMPS
echo ""
echo "Building LAMMPS (this may take 10-15 minutes)..."
make mpi -j 4

# Verify build
echo ""
echo "=== Verifying installation ==="
if [[ -f "lmp_mpi" ]]; then
    echo "✅ LAMMPS binary created successfully!"
    echo ""
    echo "Available state/change fixes:"
    ./lmp_mpi -help 2>&1 | grep "state/change" || echo "⚠️  No state/change fixes found in help"
    echo ""
    echo "Build completed at: $(date)"
else
    echo "❌ ERROR: lmp_mpi not found after build"
    exit 1
fi

