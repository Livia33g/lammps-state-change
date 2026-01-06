#!/bin/bash
# Quick installation script for fix_state_change
# This builds LAMMPS with the custom fix in your home directory

set -e  # Exit on error

echo "=== Installing fix_state_change in LAMMPS ==="

# Configuration
LAMMPS_VERSION="stable_23Aug2023"
BUILD_DIR="$HOME/lammps_build"
FIX_DIR="/work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/fix_state_change"

# Step 1: Create build directory
echo "Step 1: Creating build directory..."
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Step 2: Download LAMMPS if not already present
if [ ! -d "lammps-$LAMMPS_VERSION" ] && [ ! -d "lammps" ]; then
    echo "Step 2: Downloading LAMMPS source..."
    if command -v git &> /dev/null; then
        # Try cloning with the tag (not branch)
        if git clone --branch $LAMMPS_VERSION --single-branch --depth 1 https://github.com/lammps/lammps.git 2>/dev/null; then
            cd lammps
        elif git clone https://github.com/lammps/lammps.git; then
            cd lammps
            git checkout $LAMMPS_VERSION 2>/dev/null || git checkout stable_23Aug2023_update1 2>/dev/null || echo "⚠️  Could not checkout specific version, using latest"
        else
            echo "Git clone failed. Trying wget/curl download..."
            if command -v wget &> /dev/null; then
                wget https://github.com/lammps/lammps/archive/refs/tags/${LAMMPS_VERSION}.tar.gz -O lammps.tar.gz
                tar -xzf lammps.tar.gz
                cd lammps-${LAMMPS_VERSION}
            elif command -v curl &> /dev/null; then
                curl -L https://github.com/lammps/lammps/archive/refs/tags/${LAMMPS_VERSION}.tar.gz -o lammps.tar.gz
                tar -xzf lammps.tar.gz
                cd lammps-${LAMMPS_VERSION}
            else
                echo "❌ Error: Need git, wget, or curl to download LAMMPS"
                exit 1
            fi
        fi
    else
        echo "Git not found. Trying wget/curl..."
        if command -v wget &> /dev/null; then
            wget https://github.com/lammps/lammps/archive/refs/tags/${LAMMPS_VERSION}.tar.gz -O lammps.tar.gz
            tar -xzf lammps.tar.gz
            cd lammps-${LAMMPS_VERSION}
        elif command -v curl &> /dev/null; then
            curl -L https://github.com/lammps/lammps/archive/refs/tags/${LAMMPS_VERSION}.tar.gz -o lammps.tar.gz
            tar -xzf lammps.tar.gz
            cd lammps-${LAMMPS_VERSION}
        else
            echo "❌ Error: Need git, wget, or curl to download LAMMPS"
            exit 1
        fi
    fi
else
    echo "Step 2: LAMMPS source already exists, skipping download..."
    if [ -d "lammps-$LAMMPS_VERSION" ]; then
        cd "lammps-$LAMMPS_VERSION"
    else
        cd lammps
    fi
fi

# Step 3: Copy fix files
echo "Step 3: Copying fix files..."
cp "$FIX_DIR/fix_state_change.h" src/
cp "$FIX_DIR/fix_state_change.cpp" src/
echo "✅ Copied fix files to src/"

# Step 4: Add fix to Makefile
echo "Step 4: Adding fix to Makefile..."
cd src

# Check if already added
if grep -q "fix_state_change.cpp" Makefile.mpi; then
    echo "⚠️  fix_state_change.cpp already in Makefile.mpi, skipping..."
else
    # Find a good insertion point (after fix_rigid_nvt.cpp)
    if grep -q "fix_rigid_nvt.cpp" Makefile.mpi; then
        sed -i '/fix_rigid_nvt\.cpp \\/a\        fix_state_change.cpp \\' Makefile.mpi
        echo "✅ Added fix_state_change.cpp to Makefile.mpi"
    else
        echo "⚠️  Could not find insertion point in Makefile.mpi"
        echo "   Please manually add 'fix_state_change.cpp \\' to the source file list"
        echo "   Look for lines like 'fix_rigid_nvt.o \\' and add nearby"
    fi
fi

# Step 5: Enable required packages
echo "Step 5: Enabling required LAMMPS packages..."
make yes-rigid
make yes-molecule
echo "✅ Enabled required packages"

# Step 6: Load modules and build
echo "Step 6: Building LAMMPS..."
echo "   (This may take several minutes...)"

# Load required modules
module purge 2>/dev/null || true
module load gcc/11.4.0 2>/dev/null || echo "⚠️  Could not load gcc module, using system gcc"
module load openmpi/4.1.6 2>/dev/null || echo "⚠️  Could not load openmpi module"

# Build
if make mpi -j 4; then
    echo ""
    echo "✅✅✅ LAMMPS built successfully! ✅✅✅"
    echo ""
    echo "Your custom LAMMPS is at: $(pwd)/lmp_mpi"
    echo ""
    echo "To use it, either:"
    echo "  1. Use full path: $(pwd)/lmp_mpi -in your_input.in"
    echo "  2. Add to PATH: export PATH=$(pwd):\$PATH"
    echo ""
    echo "To verify the fix is installed:"
    echo "  $(pwd)/lmp_mpi -help | grep state/change"
    echo ""
else
    echo ""
    echo "❌ Build failed. Check the error messages above."
    echo ""
    echo "Common issues:"
    echo "  - Missing modules: module load gcc/11.4.0 openmpi/4.1.6"
    echo "  - Makefile issue: Check that fix_state_change.cpp is in Makefile.mpi"
    echo "  - Compilation errors: Check fix_state_change.cpp for syntax errors"
    echo ""
    exit 1
fi

# Step 7: Verify installation
echo "Step 7: Verifying installation..."
if ./lmp_mpi -help 2>&1 | grep -q "state/change"; then
    echo "✅ Fix verified! 'state/change' is available."
else
    echo "⚠️  Warning: Fix not found in help output, but build completed."
    echo "   Try running: ./lmp_mpi -help | grep state"
fi

echo ""
echo "=== Installation complete! ==="

