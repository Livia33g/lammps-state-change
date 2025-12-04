#!/bin/bash
# Helper script to add a new state change fix to the existing LAMMPS build
# Usage: ./add_new_fix.sh <fix_name> <path_to_fix.h> <path_to_fix.cpp>
#
# Example:
#   ./add_new_fix.sh newproblem fix_state_change_newproblem.h fix_state_change_newproblem.cpp

set -euo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 <fix_name> <path_to_fix.h> <path_to_fix.cpp>"
    echo ""
    echo "Example:"
    echo "  $0 newproblem ./fix_state_change_newproblem.h ./fix_state_change_newproblem.cpp"
    exit 1
fi

FIX_NAME="$1"
FIX_H="$2"
FIX_CPP="$3"

# Validate inputs
if [ ! -f "$FIX_H" ]; then
    echo "❌ Error: Header file not found: $FIX_H"
    exit 1
fi

if [ ! -f "$FIX_CPP" ]; then
    echo "❌ Error: Source file not found: $FIX_CPP"
    exit 1
fi

# LAMMPS source directory
LAMMPS_SRC="/work/nvme/bewl/lguttieres/lammps_build/lammps/src"

if [ ! -d "$LAMMPS_SRC" ]; then
    echo "❌ Error: LAMMPS source directory not found: $LAMMPS_SRC"
    exit 1
fi

echo "=== Adding fix_state_change_${FIX_NAME} to LAMMPS build ==="
echo ""

# Step 1: Copy files
echo "Step 1: Copying fix files..."
cp "$FIX_H" "${LAMMPS_SRC}/fix_state_change_${FIX_NAME}.h"
cp "$FIX_CPP" "${LAMMPS_SRC}/fix_state_change_${FIX_NAME}.cpp"
echo "✅ Copied files to ${LAMMPS_SRC}/"
echo ""

# Step 2: Add to Makefile
echo "Step 2: Adding to Makefile..."
cd "$LAMMPS_SRC"

# Check if already added
if grep -q "fix_state_change_${FIX_NAME}.cpp" Makefile; then
    echo "⚠️  fix_state_change_${FIX_NAME}.cpp already in Makefile, skipping..."
else
    # Find insertion point (after other fix_state_change fixes, or after fix_rigid_nvt)
    if grep -q "fix_state_change_octahedron.cpp" Makefile; then
        sed -i "/fix_state_change_octahedron\.cpp \\/a\        fix_state_change_${FIX_NAME}.cpp \\" Makefile
        echo "✅ Added fix_state_change_${FIX_NAME}.cpp to Makefile (after fix_state_change_octahedron)"
    elif grep -q "fix_state_change.cpp" Makefile; then
        sed -i "/fix_state_change\.cpp \\/a\        fix_state_change_${FIX_NAME}.cpp \\" Makefile
        echo "✅ Added fix_state_change_${FIX_NAME}.cpp to Makefile (after fix_state_change)"
    elif grep -q "fix_rigid_nvt.cpp" Makefile; then
        sed -i "/fix_rigid_nvt\.cpp \\/a\        fix_state_change_${FIX_NAME}.cpp \\" Makefile
        echo "✅ Added fix_state_change_${FIX_NAME}.cpp to Makefile (after fix_rigid_nvt)"
    else
        echo "⚠️  Could not find insertion point in Makefile"
        echo "   Please manually add 'fix_state_change_${FIX_NAME}.cpp \\' to the source file list"
        exit 1
    fi
fi
echo ""

# Step 3: Rebuild
echo "Step 3: Rebuilding LAMMPS..."
echo "   (This may take a few minutes...)"
echo ""

# Check if we need to load modules
if ! command -v mpicxx &> /dev/null; then
    echo "⚠️  mpicxx not found in PATH. You may need to load modules:"
    echo "   module load gcc/11.4.0 openmpi/4.1.6"
    echo ""
fi

# Rebuild
if make mpi -j 4; then
    echo ""
    echo "✅ LAMMPS rebuilt successfully!"
else
    echo ""
    echo "❌ Build failed. Please check the error messages above."
    exit 1
fi
echo ""

# Step 4: Verify
echo "Step 4: Verifying installation..."
if ./lmp_mpi -help 2>&1 | grep -q "state/change"; then
    echo "✅ Fix registered successfully!"
    echo ""
    echo "Available state/change fixes:"
    ./lmp_mpi -help 2>&1 | grep "state/change"
    echo ""
    echo "To use this fix in your LAMMPS input script:"
    echo "  fix 1 all state/change/${FIX_NAME} <arguments>"
else
    echo "⚠️  Warning: Could not verify fix registration. Check manually with:"
    echo "   ./lmp_mpi -help | grep state/change"
fi

echo ""
echo "=== Done! ==="

