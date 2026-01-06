#!/bin/bash
#SBATCH --job-name=rebuild_lammps_all
#SBATCH --account=bewl-delta-gpu
#SBATCH --partition=gpuA100x4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --output=rebuild_all_fixes-%j.out
#SBATCH --error=rebuild_all_fixes-%j.err
# Note: Requesting GPU to use GPU partition, but only using CPU for compilation

set -euo pipefail

echo "=========================================="
echo "=== FULL REBUILD LAMMPS WITH ALL FIXES ==="
echo "=========================================="
echo "Job started at: $(date)"
echo ""

LAMMPS_SRC="/work/nvme/bewl/lguttieres/lammps_build/lammps/src"
STATE_CHANGE_DIR="/work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change"

cd "$LAMMPS_SRC"
echo "Working directory: $(pwd)"
echo ""

# ------------------------------------------------------------
# 1. Environment Setup
# ------------------------------------------------------------
echo "--- Environment Setup ---"

# Deactivate conda if active (conda compilers cause linker issues)
if [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then
    echo "⚠️  Conda environment '$CONDA_DEFAULT_ENV' is active - deactivating"
    conda deactivate
fi

# Load modules
module purge

# Try loading modules (similar to other working scripts)
module load gcc/11.4.0 2>/dev/null || echo "⚠️  gcc/11.4.0 module not available, trying system gcc"
module load openmpi/4.1.6 2>/dev/null || echo "⚠️  openmpi/4.1.6 module not available, trying system MPI"

# Try alternative Cray modules if above didn't work
if ! command -v mpicxx &>/dev/null; then
    echo "⚠️  mpicxx not found, trying Cray modules..."
    module load gcc-native/13.2 2>/dev/null || true
    module load cray-mpich/8.1.32 2>/dev/null || true
    module load craype-network-ofi 2>/dev/null || true
fi

# Verify mpicxx is available
if command -v mpicxx &>/dev/null; then
    echo "✅ Found mpicxx at: $(which mpicxx)"
    mpicxx --version 2>&1 | head -1 || true
else
    echo "❌ ERROR: mpicxx not found after module loading"
    echo "   Tried: gcc/11.4.0, openmpi/4.1.6, gcc-native/13.2, cray-mpich/8.1.32"
    echo "   PATH: $PATH"
    echo "   Cannot proceed with build"
    exit 1
fi

echo "✅ Environment ready"
echo ""

# ------------------------------------------------------------
# 2. Copy All Fix Files
# ------------------------------------------------------------
echo "--- Copying Fix Files ---"

# Original dimer fix
if [[ -f "$STATE_CHANGE_DIR/fix_state_change/fix_state_change.cpp" ]]; then
    cp "$STATE_CHANGE_DIR/fix_state_change/fix_state_change.cpp" .
    cp "$STATE_CHANGE_DIR/fix_state_change/fix_state_change.h" .
    echo "✅ Copied fix_state_change (original dimer)"
else
    echo "⚠️  WARNING: Original fix_state_change not found at $STATE_CHANGE_DIR/fix_state_change/"
fi

# Octahedron fix
if [[ -f "$STATE_CHANGE_DIR/octahedron/fix_state_change_octahedron.cpp" ]]; then
    cp "$STATE_CHANGE_DIR/octahedron/fix_state_change_octahedron.cpp" .
    cp "$STATE_CHANGE_DIR/octahedron/fix_state_change_octahedron.h" .
    echo "✅ Copied fix_state_change_octahedron"
else
    echo "⚠️  WARNING: Octahedron fix not found at $STATE_CHANGE_DIR/octahedron/"
fi

# Ksat fix
if [[ -f "$STATE_CHANGE_DIR/ksat/fix_state_change_ksat.cpp" ]]; then
    cp "$STATE_CHANGE_DIR/ksat/fix_state_change_ksat.cpp" .
    cp "$STATE_CHANGE_DIR/ksat/fix_state_change_ksat.h" .
    echo "✅ Copied fix_state_change_ksat"
else
    echo "⚠️  WARNING: Ksat fix not found at $STATE_CHANGE_DIR/ksat/"
fi

echo ""

# Verify all fix files are present
echo "--- Verifying Fix Files ---"
ls -lh fix_state_change*.{cpp,h} 2>/dev/null | grep -v "^total" || echo "⚠️  No fix files found!"
echo ""

# ------------------------------------------------------------
# 3. Update Makefile
# ------------------------------------------------------------
echo "--- Updating Makefile ---"

# Check if Makefile uses wildcard (auto-detects .cpp files)
if grep -q "wildcard.*fix.*cpp" Makefile 2>/dev/null || grep -q '\$.*wildcard.*\.cpp' Makefile 2>/dev/null; then
    echo "✅ Makefile uses wildcard - will auto-detect fix files"
else
    echo "⚠️  Makefile doesn't use wildcard - checking if fixes are listed..."
    
    # Add fixes to Makefile if not present
    for fix in fix_state_change fix_state_change_octahedron fix_state_change_ksat; do
        if [[ -f "${fix}.cpp" ]] && ! grep -q "^${fix}\.cpp" Makefile 2>/dev/null; then
            echo "Adding ${fix}.cpp to Makefile..."
            # Try to add after fix_state_change.cpp or similar fix
            if grep -q "^fix_state_change\.cpp" Makefile 2>/dev/null; then
                sed -i "/^fix_state_change\.cpp/a\\${fix}.cpp \\" Makefile
                echo "✅ Added ${fix}.cpp after fix_state_change.cpp"
            else
                echo "⚠️  Could not auto-add ${fix}.cpp - may need manual Makefile edit"
            fi
        fi
    done
fi
echo ""

# ------------------------------------------------------------
# 4. Full Clean Rebuild
# ------------------------------------------------------------
echo "--- Full Clean Rebuild ---"
echo "⚠️  This will delete all compiled objects and rebuild from scratch!"
echo ""

# Remove old build artifacts
echo "Cleaning old build artifacts..."
rm -rf Obj_mpi
rm -f liblammps_mpi.a lmp_mpi
echo "✅ Cleaned old build files"
echo ""

# Build LAMMPS
echo "Building LAMMPS (this will take 10-15 minutes)..."
echo "Using: make mpi -j 4"
make mpi -j 4

echo ""

# ------------------------------------------------------------
# 5. Verify Build
# ------------------------------------------------------------
echo "=== Verifying Installation ==="

if [[ -f "lmp_mpi" ]]; then
    echo "✅ LAMMPS binary created successfully!"
    echo ""
    echo "File size:"
    ls -lh lmp_mpi
    echo ""
    echo "Available state/change fixes:"
    ./lmp_mpi -help 2>&1 | grep "state/change" || echo "⚠️  No state/change fixes found in help"
    echo ""
    echo "Build completed at: $(date)"
    echo ""
    echo "✅ SUCCESS! All fixes should now be available."
else
    echo "❌ ERROR: lmp_mpi not found after build"
    echo "Build failed - check error messages above"
    exit 1
fi

