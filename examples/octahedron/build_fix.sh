#!/bin/bash
# Build LAMMPS with fix_state_change_octahedron
# IMPORTANT: Run this WITHOUT conda environment active!

set -euo pipefail

LAMMPS_SRC="/work/nvme/bewl/lguttieres/lammps_build/lammps/src"
cd "$LAMMPS_SRC"

echo "=== Building LAMMPS with fix_state_change_octahedron ==="
echo "⚠️  Make sure conda environment is DEACTIVATED!"
echo ""

# Check if we're in a conda environment
if [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then
    echo "❌ ERROR: Conda environment '$CONDA_DEFAULT_ENV' is active!"
    echo "   Please run: conda deactivate"
    echo "   Then run this script again"
    exit 1
fi

# Try to load modules, but continue if they don't exist
module purge 2>/dev/null || true
module load gcc/11.4.0 2>/dev/null || module load gcc-native/13.2 2>/dev/null || echo "⚠️  Using system gcc"
module load openmpi/4.1.6 2>/dev/null || echo "⚠️  Using system mpicxx"

echo "Using compiler: $(which mpicxx 2>/dev/null || which g++ || echo 'system compiler')"
if command -v mpicxx &>/dev/null; then
    echo "Compiler version: $(mpicxx --version 2>&1 | head -1)"
elif command -v g++ &>/dev/null; then
    echo "Compiler version: $(g++ --version 2>&1 | head -1)"
fi
echo ""

# Check if fix file exists
if [[ ! -f "fix_state_change_octahedron.cpp" ]]; then
    echo "❌ ERROR: fix_state_change_octahedron.cpp not found"
    echo "   Copying from project directory..."
    cp /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/octahedron/fix_state_change_octahedron.* .
fi

# Check if fix is in Makefile
if ! grep -q "fix_state_change_octahedron.cpp" Makefile 2>/dev/null; then
    echo "Adding fix_state_change_octahedron.cpp to Makefile..."
    # Try to add after fix_state_change.cpp
    if grep -q "fix_state_change.cpp" Makefile; then
        sed -i '/^fix_state_change\.cpp/a\fix_state_change_octahedron.cpp' Makefile
        echo "✅ Added to Makefile"
    else
        echo "⚠️  Could not find fix_state_change.cpp in Makefile"
        echo "   Please manually add 'fix_state_change_octahedron.cpp' to Makefile"
        exit 1
    fi
fi

echo ""
echo "=== Building LAMMPS (this may take several minutes) ==="
make clean || true
make mpi -j 4

echo ""
echo "=== Verifying installation ==="
if [[ -f "lmp_mpi" ]]; then
    if ./lmp_mpi -help 2>&1 | grep -q "state/change/octahedron"; then
        echo "✅ SUCCESS! fix_state_change/octahedron is installed!"
        echo ""
        ./lmp_mpi -help 2>&1 | grep "state/change"
    else
        echo "⚠️  WARNING: LAMMPS built but fix not found in help"
        echo "   The fix may still work, but verification failed"
    fi
else
    echo "❌ ERROR: lmp_mpi not found after build"
    exit 1
fi
