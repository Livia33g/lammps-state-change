#!/bin/bash
# Full clean rebuild of LAMMPS with fix_state_change_octahedron

set -euo pipefail

LAMMPS_SRC="/work/nvme/bewl/lguttieres/lammps_build/lammps/src"
cd "$LAMMPS_SRC"

echo "=== FULL CLEAN REBUILD of LAMMPS ==="
echo "⚠️  This will delete all compiled objects and rebuild from scratch!"
echo ""

# Check if we're in a conda environment
if [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then
    echo "❌ ERROR: Conda environment '$CONDA_DEFAULT_ENV' is active!"
    echo "   Please run: conda deactivate"
    exit 1
fi

# Full clean - remove object directories and libraries
echo "Removing old build artifacts..."
rm -rf Obj_mpi
rm -f liblammps_mpi.a lmp_mpi
echo "✅ Cleaned old build files"

# Make sure fix file exists
if [[ ! -f "fix_state_change_octahedron.cpp" ]]; then
    echo "Copying fix files..."
    cp /work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/octahedron/fix_state_change_octahedron.* .
fi

# Add fix to Makefile if needed
if ! grep -q "fix_state_change_octahedron.cpp" Makefile 2>/dev/null; then
    echo "Adding fix_state_change_octahedron.cpp to Makefile..."
    if grep -q "fix_state_change.cpp" Makefile; then
        sed -i '/^fix_state_change\.cpp/a\fix_state_change_octahedron.cpp' Makefile
        echo "✅ Added to Makefile"
    else
        echo "❌ Could not find fix_state_change.cpp in Makefile"
        exit 1
    fi
fi

echo ""
echo "=== Building LAMMPS from scratch (this will take 10-15 minutes) ==="
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
    fi
else
    echo "❌ ERROR: lmp_mpi not found after build"
    exit 1
fi

