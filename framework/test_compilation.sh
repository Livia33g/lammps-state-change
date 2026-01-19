#!/bin/bash
# Script to test compilation of generated fix
# Usage: ./test_compilation.sh [path_to_lammps_src]
#
# Example:
#   ./test_compilation.sh /work/nvme/bewl/lguttieres/lammps_build/lammps/src

set -euo pipefail

# Configuration
GENERATED_CPP="framework/generated_test/fix_state_change_octahedron_avoid_same_type_octahedron.cpp"
GENERATED_H="framework/generated_test/fix_state_change_octahedron_avoid_same_type_octahedron.h"
FIX_NAME="state_change_octahedron_avoid_same_type_octahedron"

# Determine LAMMPS source directory
if [ $# -eq 1 ]; then
    LAMMPS_SRC="$1"
else
    # Default location
    LAMMPS_SRC="/work/nvme/bewl/lguttieres/lammps_build/lammps/src"
fi

if [ ! -d "$LAMMPS_SRC" ]; then
    echo "❌ Error: LAMMPS source directory not found: $LAMMPS_SRC"
    echo ""
    echo "Usage: $0 [path_to_lammps_src]"
    exit 1
fi

echo "==================================================="
echo "COMPILATION TEST: Generated Fix"
echo "==================================================="
echo ""
echo "Generated fix: fix_${FIX_NAME}"
echo "LAMMPS src:    $LAMMPS_SRC"
echo ""

# Check that generated files exist
if [ ! -f "$GENERATED_CPP" ]; then
    echo "❌ Error: Generated .cpp file not found: $GENERATED_CPP"
    exit 1
fi

if [ ! -f "$GENERATED_H" ]; then
    echo "❌ Error: Generated .h file not found: $GENERATED_H"
    exit 1
fi

echo "Step 1: Copying generated files to LAMMPS src/"
echo "---------------------------------------------------"
cp -v "$GENERATED_CPP" "${LAMMPS_SRC}/fix_${FIX_NAME}.cpp"
cp -v "$GENERATED_H" "${LAMMPS_SRC}/fix_${FIX_NAME}.h"
echo "✅ Files copied successfully"
echo ""

echo "Step 2: Adding to Makefile (if not already present)"
echo "---------------------------------------------------"
cd "$LAMMPS_SRC"

if grep -q "fix_${FIX_NAME}.cpp" Makefile; then
    echo "⚠️  fix_${FIX_NAME}.cpp already in Makefile"
else
    # Try to add after fix_state_change_octahedron.cpp
    if grep -q "fix_state_change_octahedron.cpp" Makefile; then
        sed -i "/fix_state_change_octahedron\.cpp \\/a\        fix_${FIX_NAME}.cpp \\" Makefile
        echo "✅ Added fix_${FIX_NAME}.cpp to Makefile (after fix_state_change_octahedron)"
    elif grep -q "fix_state_change.cpp" Makefile; then
        sed -i "/fix_state_change\.cpp \\/a\        fix_${FIX_NAME}.cpp \\" Makefile
        echo "✅ Added fix_${FIX_NAME}.cpp to Makefile (after fix_state_change)"
    else
        echo "⚠️  Could not find insertion point in Makefile"
        echo "   Please manually add 'fix_${FIX_NAME}.cpp \\' to the source file list"
        exit 1
    fi
fi
echo ""

echo "Step 3: Loading required modules"
echo "---------------------------------------------------"
if command -v module &> /dev/null; then
    echo "Loading gcc and openmpi modules..."
    module purge 2>/dev/null || true
    module load gcc/11.4.0 2>/dev/null || echo "⚠️  gcc module not available"
    module load openmpi/4.1.6 2>/dev/null || echo "⚠️  openmpi module not available"
else
    echo "⚠️  module command not found, assuming compilers are in PATH"
fi

# Verify compiler is available
if ! command -v mpicxx &> /dev/null; then
    echo "❌ Error: mpicxx not found in PATH"
    echo "   You may need to load modules or set up your environment"
    exit 1
fi
echo "✅ Compiler available: $(mpicxx --version | head -1)"
echo ""

echo "Step 4: Compiling LAMMPS"
echo "---------------------------------------------------"
echo "This may take several minutes..."
echo ""

# Clean and rebuild
make clean || true

if make mpi -j 4; then
    echo ""
    echo "✅ COMPILATION SUCCESSFUL!"
    echo ""
else
    echo ""
    echo "❌ COMPILATION FAILED!"
    echo ""
    echo "Common issues:"
    echo "  1. Missing dependencies (check error messages above)"
    echo "  2. Compiler warnings treated as errors (check CCFLAGS in Makefile)"
    echo "  3. Syntax errors (unlikely - code passed static analysis)"
    echo ""
    echo "Next steps:"
    echo "  - Review error messages above"
    echo "  - Check compiler version compatibility"
    echo "  - Try compiling with warnings disabled: CCFLAGS=\"-w\""
    exit 1
fi

echo "Step 5: Verifying fix registration"
echo "---------------------------------------------------"
if ./lmp_mpi -help 2>&1 | grep -q "state/change/octahedron/avoid_same_type_octahedron"; then
    echo "✅ Fix registered successfully!"
    echo ""
    echo "Available state/change fixes:"
    ./lmp_mpi -help 2>&1 | grep "state/change"
    echo ""
else
    echo "⚠️  Warning: Could not find fix in help output"
    echo "   Try manually: ./lmp_mpi -help | grep state/change"
fi
echo ""

echo "==================================================="
echo "COMPILATION TEST: SUCCESS! ✅"
echo "==================================================="
echo ""
echo "The generated fix has been successfully compiled!"
echo ""
echo "To use this fix in a LAMMPS input script:"
echo "  fix ID group state/change/octahedron/avoid_same_type_octahedron \\"
echo "      check_every cooldown probability cutoff group_patches hysteresis"
echo ""
echo "Example:"
echo "  fix 1 patches state/change/octahedron/avoid_same_type_octahedron \\"
echo "      100 10000 1.0 2.5 patches 1000"
echo ""
echo "Next step: Run a test simulation to verify behavior"
echo "See: COMPILATION_READINESS_REPORT.md (Phase 2: Runtime Test)"
echo ""
