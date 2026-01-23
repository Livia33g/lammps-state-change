#!/bin/bash
# Quick script to rebuild LAMMPS with the new octahedron fix

set -euo pipefail

LAMMPS_SRC="/work/nvme/bewl/lguttieres/lammps_build/lammps/src"
cd "$LAMMPS_SRC"

echo "=== Rebuilding LAMMPS with fix_state_change_octahedron ==="

# Check if fix file exists
if [[ ! -f "fix_state_change_octahedron.cpp" ]]; then
    echo "❌ ERROR: fix_state_change_octahedron.cpp not found in $LAMMPS_SRC"
    exit 1
fi

# Load modules
module purge
module load gcc/11.4.0 2>/dev/null || echo "⚠️  gcc module not available"
module load openmpi/4.1.6 2>/dev/null || echo "⚠️  openmpi module not available"

# Check build system
if [[ -f "Makefile" ]]; then
    echo "Using Makefile build system..."
    # Check if fix_state_change_octahedron.cpp is in Makefile
    if ! grep -q "fix_state_change_octahedron.cpp" Makefile; then
        echo "Adding fix_state_change_octahedron.cpp to Makefile..."
        # Try to add after fix_state_change.cpp
        if grep -q "fix_state_change.cpp" Makefile; then
            sed -i '/^fix_state_change\.cpp/a\fix_state_change_octahedron.cpp' Makefile
        else
            echo "⚠️  Could not find fix_state_change.cpp in Makefile. Please add manually:"
            echo "   Add 'fix_state_change_octahedron.cpp' to the list of fix source files"
            exit 1
        fi
    fi
    make clean || true
    make mpi -j 4
elif [[ -f "CMakeLists.txt" ]]; then
    echo "Using CMake build system..."
    echo "Note: CMake should auto-detect fix files in src/"
    echo "Rebuilding with CMake..."
    # You may need to adjust this based on your CMake setup
    cmake --build . -j 4 || {
        echo "❌ CMake build failed. You may need to reconfigure CMake."
        exit 1
    }
else
    echo "❌ Could not determine build system. Please rebuild manually."
    exit 1
fi

echo ""
echo "=== Verifying installation ==="
if ./lmp_mpi -help 2>&1 | grep -q "state/change/octahedron"; then
    echo "✅ SUCCESS! fix_state_change/octahedron is installed!"
    echo ""
    ./lmp_mpi -help 2>&1 | grep "state/change"
else
    echo "❌ WARNING: fix_state_change/octahedron not found after rebuild"
    echo "   It may need to be added manually to the build system"
    exit 1
fi

