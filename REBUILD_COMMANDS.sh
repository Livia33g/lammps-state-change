#!/bin/bash
# Commands to rebuild LAMMPS with fixed ksat and octahedron fixes
# Copy and paste these commands into your terminal

# Set paths
LAMMPS_SRC="/work/nvme/bewl/lguttieres/lammps_build/lammps/src"
KSAT_DIR="/work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/ksat"
OCTAHEDRON_DIR="/work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/octahedron"

# Step 1: Navigate to LAMMPS source directory
cd "$LAMMPS_SRC"
echo "✅ Changed to LAMMPS source directory: $(pwd)"

# Step 2: Copy fixed ksat files
echo ""
echo "--- Copying ksat fix files ---"
cp "$KSAT_DIR/fix_state_change_ksat.cpp" .
cp "$KSAT_DIR/fix_state_change_ksat.h" .
echo "✅ Copied ksat fix files"

# Step 3: Copy fixed octahedron files
echo ""
echo "--- Copying octahedron fix files ---"
cp "$OCTAHEDRON_DIR/fix_state_change_octahedron.cpp" .
cp "$OCTAHEDRON_DIR/fix_state_change_octahedron.h" .
echo "✅ Copied octahedron fix files"

# Step 4: Verify files are present
echo ""
echo "--- Verifying fix files ---"
ls -lh fix_state_change_ksat.* fix_state_change_octahedron.* 2>/dev/null || echo "⚠️  Some files missing!"

# Step 5: Check if Makefile needs updating
echo ""
echo "--- Checking Makefile ---"
if grep -q "fix_state_change_ksat.cpp" Makefile 2>/dev/null; then
    echo "✅ fix_state_change_ksat.cpp already in Makefile"
else
    echo "⚠️  fix_state_change_ksat.cpp NOT in Makefile - adding it..."
    if grep -q "fix_state_change.cpp" Makefile; then
        sed -i '/^fix_state_change\.cpp/a\fix_state_change_ksat.cpp' Makefile
        echo "✅ Added fix_state_change_ksat.cpp to Makefile"
    else
        echo "❌ Could not find insertion point in Makefile"
        echo "   Please manually add 'fix_state_change_ksat.cpp' to Makefile"
    fi
fi

if grep -q "fix_state_change_octahedron.cpp" Makefile 2>/dev/null; then
    echo "✅ fix_state_change_octahedron.cpp already in Makefile"
else
    echo "⚠️  fix_state_change_octahedron.cpp NOT in Makefile - adding it..."
    if grep -q "fix_state_change_ksat.cpp" Makefile; then
        sed -i '/^fix_state_change_ksat\.cpp/a\fix_state_change_octahedron.cpp' Makefile
        echo "✅ Added fix_state_change_octahedron.cpp to Makefile"
    elif grep -q "fix_state_change.cpp" Makefile; then
        sed -i '/^fix_state_change\.cpp/a\fix_state_change_octahedron.cpp' Makefile
        echo "✅ Added fix_state_change_octahedron.cpp to Makefile"
    else
        echo "❌ Could not find insertion point in Makefile"
        echo "   Please manually add 'fix_state_change_octahedron.cpp' to Makefile"
    fi
fi

# Step 6: Load modules (if available)
echo ""
echo "--- Loading modules ---"
module purge 2>/dev/null || true
module load gcc/11.4.0 2>/dev/null || echo "⚠️  gcc/11.4.0 not available, using system gcc"
module load openmpi/4.1.6 2>/dev/null || echo "⚠️  openmpi/4.1.6 not available, using system MPI"

# Step 7: Clean old build (optional - comment out if you want incremental build)
echo ""
echo "--- Cleaning old build (optional) ---"
read -p "Do you want to clean old build? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    rm -rf Obj_mpi
    rm -f liblammps_mpi.a lmp_mpi 2>/dev/null || true
    echo "✅ Cleaned old build files"
else
    echo "⏭️  Skipping clean (incremental build)"
fi

# Step 8: Build LAMMPS
echo ""
echo "--- Building LAMMPS (this will take 10-15 minutes) ---"
echo "   Using: $(which mpicxx 2>/dev/null || which g++ || echo 'system compiler')"
make mpi -j 4

# Step 9: Verify build
echo ""
echo "--- Verifying installation ---"
if [[ -f "lmp_mpi" ]]; then
    echo "✅ LAMMPS binary created successfully!"
    echo ""
    echo "Available state/change fixes:"
    ./lmp_mpi -help 2>&1 | grep "state/change" || echo "⚠️  No state/change fixes found in help"
    echo ""
    echo "✅✅✅ BUILD COMPLETE! ✅✅✅"
    echo "LAMMPS binary location: $LAMMPS_SRC/lmp_mpi"
else
    echo "❌ ERROR: lmp_mpi not found after build"
    exit 1
fi

