#!/bin/bash
# Manual rebuild commands for LAMMPS with all state change fixes
# 
# IMPORTANT: This script should be run on an interactive compute node
# To get an interactive node:
#   1. ssh to gpuA (or use: salloc --partition=gpuA100x4 --account=bewl-delta-gpu --nodes=1 --ntasks-per-node=1 --gres=gpu:1 --time=01:00:00)
#   2. Then run this script or the commands manually
#
# Run this script or copy commands one by one

set -e  # Exit on error

echo "=========================================="
echo "=== MANUAL REBUILD LAMMPS WITH ALL FIXES ==="
echo "=========================================="
echo "NOTE: This should be run on an interactive compute node (ssh to gpuA)"
echo ""

# ------------------------------------------------------------
# 1. Environment Setup
# ------------------------------------------------------------
echo "--- Environment Setup ---"

# Deactivate conda if active
if [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then
    echo "⚠️  Conda environment '$CONDA_DEFAULT_ENV' is active - deactivating"
    conda deactivate
fi

# Add Cray MPI to PATH (mpicxx is at this location)
export PATH=/opt/cray/pe/mpich/8.1.32/ofi/gnu/11.2/bin:$PATH

# Verify mpicxx is available
if command -v mpicxx &>/dev/null; then
    echo "✅ Found mpicxx at: $(which mpicxx)"
    mpicxx --version 2>&1 | head -1 || true
else
    echo "❌ ERROR: mpicxx not found"
    echo "   Check if PATH is correct or load MPI modules"
    exit 1
fi

echo ""

# ------------------------------------------------------------
# 2. Navigate to LAMMPS source
# ------------------------------------------------------------
LAMMPS_SRC="/work/nvme/bewl/lguttieres/lammps_build/lammps/src"
STATE_CHANGE_DIR="/work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change"

cd "$LAMMPS_SRC"
echo "Working directory: $(pwd)"
echo ""

# ------------------------------------------------------------
# 3. Copy All Fix Files
# ------------------------------------------------------------
echo "--- Copying Fix Files ---"

# Original dimer fix
if [[ -f "$STATE_CHANGE_DIR/fix_state_change/fix_state_change.cpp" ]]; then
    cp "$STATE_CHANGE_DIR/fix_state_change/fix_state_change.cpp" .
    cp "$STATE_CHANGE_DIR/fix_state_change/fix_state_change.h" .
    echo "✅ Copied fix_state_change (original dimer)"
else
    echo "⚠️  WARNING: Original fix_state_change not found"
fi

# Octahedron fix
if [[ -f "$STATE_CHANGE_DIR/octahedron/fix_state_change_octahedron.cpp" ]]; then
    cp "$STATE_CHANGE_DIR/octahedron/fix_state_change_octahedron.cpp" .
    cp "$STATE_CHANGE_DIR/octahedron/fix_state_change_octahedron.h" .
    echo "✅ Copied fix_state_change_octahedron"
else
    echo "⚠️  WARNING: Octahedron fix not found"
fi

# Ksat fix
if [[ -f "$STATE_CHANGE_DIR/ksat/fix_state_change_ksat.cpp" ]]; then
    cp "$STATE_CHANGE_DIR/ksat/fix_state_change_ksat.cpp" .
    cp "$STATE_CHANGE_DIR/ksat/fix_state_change_ksat.h" .
    echo "✅ Copied fix_state_change_ksat"
else
    echo "⚠️  WARNING: Ksat fix not found"
fi

echo ""

# Verify all fix files are present
echo "--- Verifying Fix Files ---"
ls -lh fix_state_change*.{cpp,h} 2>/dev/null | grep -v "^total" || echo "⚠️  No fix files found!"
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

