#!/bin/bash
# Simplified installation script - uses wget/curl to download tarball directly

set -e

echo "=== Installing fix_state_change in LAMMPS ==="

BUILD_DIR="$HOME/lammps_build"
FIX_DIR="/work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/fix_state_change"
LAMMPS_VERSION="stable_23Aug2023"
LAMMPS_URL="https://github.com/lammps/lammps/archive/refs/tags/${LAMMPS_VERSION}.tar.gz"

# Step 1: Create build directory
echo "Step 1: Creating build directory..."
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Step 2: Download LAMMPS
if [ ! -d "lammps-${LAMMPS_VERSION}" ] && [ ! -d "lammps" ]; then
    echo "Step 2: Downloading LAMMPS source..."
    if [ -f "${LAMMPS_VERSION}.tar.gz" ]; then
        echo "   Found existing tarball, extracting..."
        tar -xzf "${LAMMPS_VERSION}.tar.gz"
    else
        if command -v wget &> /dev/null; then
            echo "   Downloading with wget..."
            wget "$LAMMPS_URL" -O "${LAMMPS_VERSION}.tar.gz"
        elif command -v curl &> /dev/null; then
            echo "   Downloading with curl..."
            curl -L "$LAMMPS_URL" -o "${LAMMPS_VERSION}.tar.gz"
        else
            echo "❌ Error: Need wget or curl to download LAMMPS"
            exit 1
        fi
        tar -xzf "${LAMMPS_VERSION}.tar.gz"
    fi
    cd "lammps-${LAMMPS_VERSION}"
else
    echo "Step 2: LAMMPS source already exists, skipping download..."
    if [ -d "lammps-${LAMMPS_VERSION}" ]; then
        cd "lammps-${LAMMPS_VERSION}"
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

if grep -q "fix_state_change.cpp" Makefile.mpi; then
    echo "⚠️  fix_state_change.cpp already in Makefile.mpi, skipping..."
else
    # Find insertion point after fix_rigid_nvt.cpp
    if grep -q "fix_rigid_nvt.cpp" Makefile.mpi; then
        sed -i '/fix_rigid_nvt\.cpp \\/a\        fix_state_change.cpp \\' Makefile.mpi
        echo "✅ Added fix_state_change.cpp to Makefile.mpi"
    else
        echo "⚠️  Could not find insertion point. Please manually add 'fix_state_change.cpp \\' to Makefile.mpi"
        echo "   Look for lines with 'fix_*.cpp' and add it there"
    fi
fi

# Step 5: Enable packages
echo "Step 5: Enabling required packages..."
make yes-rigid
make yes-molecule
echo "✅ Enabled packages"

# Step 6: Build
echo "Step 6: Building LAMMPS (this will take ~10-15 minutes)..."
module purge 2>/dev/null || true
module load gcc/11.4.0 2>/dev/null || echo "⚠️  Using system gcc"
module load openmpi/4.1.6 2>/dev/null || echo "⚠️  Using system MPI"

if make mpi -j 4; then
    echo ""
    echo "✅✅✅ LAMMPS built successfully! ✅✅✅"
    echo ""
    echo "Your LAMMPS is at: $(pwd)/lmp_mpi"
    echo ""
    if ./lmp_mpi -help 2>&1 | grep -q "state/change"; then
        echo "✅ Fix verified! 'state/change' is available."
    else
        echo "⚠️  Fix not found in help, but build completed."
    fi
else
    echo "❌ Build failed. Check errors above."
    exit 1
fi

