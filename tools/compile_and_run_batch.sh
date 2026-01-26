#!/bin/bash
# Batch process submissions: generate, compile, and optionally run simulations
#
# This script efficiently processes multiple submissions with:
# - Isolated build directories (cleaned up after)
# - Shared LAMMPS base (if same policy, reuse compilation)
# - Batch processing options
# - Automatic cleanup
#
# Usage:
#   tools/compile_and_run_batch.sh problem-001-dimer-ksat [--compile-only] [--run] [--cleanup]
#   tools/compile_and_run_batch.sh problem-001-dimer-ksat --username alice --username bob

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$SCRIPT_DIR"

PROBLEM_ID="${1}"
shift || true

if [ -z "$PROBLEM_ID" ]; then
    echo "Usage: tools/compile_and_run_batch.sh problem-001-dimer-ksat [options]"
    echo ""
    echo "Options:"
    echo "  --compile-only    Only compile, don't run simulations"
    echo "  --run             Run simulations after compilation"
    echo "  --cleanup         Clean up build directories after (default: yes)"
    echo "  --username NAME   Process only specific username(s)"
    echo "  --keep-builds     Keep build directories (for debugging)"
    echo ""
    exit 1
fi

# Parse options
COMPILE_ONLY=false
RUN_SIMULATIONS=false
CLEANUP=true
KEEP_BUILDS=false
USERNAMES=()

while [[ $# -gt 0 ]]; do
    case $1 in
        --compile-only)
            COMPILE_ONLY=true
            shift
            ;;
        --run)
            RUN_SIMULATIONS=true
            shift
            ;;
        --cleanup)
            CLEANUP=true
            shift
            ;;
        --keep-builds)
            KEEP_BUILDS=true
            CLEANUP=false
            shift
            ;;
        --username)
            USERNAMES+=("$2")
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

SUBMISSIONS_DIR="submissions/$PROBLEM_ID"
BUILD_BASE_DIR=".builds/$PROBLEM_ID"  # Temporary build directory
mkdir -p "$BUILD_BASE_DIR"

if [ ! -d "$SUBMISSIONS_DIR" ]; then
    echo "Error: Submissions directory not found: $SUBMISSIONS_DIR"
    exit 1
fi

# Find submissions to process
if [ ${#USERNAMES[@]} -eq 0 ]; then
    # Process all submissions
    SUBMISSION_DIRS=$(find "$SUBMISSIONS_DIR" -mindepth 1 -maxdepth 1 -type d ! -name ".*" | sort)
else
    # Process only specified usernames
    SUBMISSION_DIRS=""
    for USERNAME in "${USERNAMES[@]}"; do
        DIR="$SUBMISSIONS_DIR/$USERNAME"
        if [ -d "$DIR" ]; then
            SUBMISSION_DIRS="$SUBMISSION_DIRS $DIR"
        else
            echo "Warning: Submission not found: $DIR"
        fi
    done
fi

if [ -z "$SUBMISSION_DIRS" ]; then
    echo "No submissions found to process"
    exit 0
fi

COUNT=$(echo "$SUBMISSION_DIRS" | wc -w)
echo "=========================================="
echo "Batch Processing Submissions"
echo "=========================================="
echo "Problem: $PROBLEM_ID"
echo "Submissions: $COUNT"
echo "Build directory: $BUILD_BASE_DIR"
echo ""

# Track policy hashes to reuse compilations
declare -A POLICY_HASHES
declare -A BUILD_DIRS

# Step 1: Process all submissions (generate C++ and LAMMPS files)
echo "Step 1: Processing submissions (generating files)..."
PROCESSED=0
FAILED=0
SKIPPED=0

for SUB_DIR in $SUBMISSION_DIRS; do
    USERNAME=$(basename "$SUB_DIR")
    echo ""
    echo "----------------------------------------"
    echo "Processing: $USERNAME"
    echo "----------------------------------------"
    
    # Check if already processed
    if [ -d "$SUB_DIR/generated" ] && [ -n "$(find "$SUB_DIR/generated" -name "*.h" -o -name "*.cpp" 2>/dev/null)" ] && \
       [ -d "$SUB_DIR/simulation" ] && [ -n "$(find "$SUB_DIR/simulation" -name "data.*" -o -name "in.*" 2>/dev/null)" ]; then
        echo "✓ Already processed - skipping"
        ((SKIPPED++))
        
        # Still calculate hash for compilation check
        if [ -f "$SUB_DIR/policy.json" ]; then
            POLICY_HASH=$(md5sum "$SUB_DIR/policy.json" | cut -d' ' -f1)
            POLICY_HASHES["$USERNAME"]="$POLICY_HASH"
        fi
        continue
    fi
    
    if tools/process_submission.sh "$SUB_DIR" 2>&1 | tee "$BUILD_BASE_DIR/${USERNAME}_process.log"; then
        ((PROCESSED++))
        
        # Calculate policy hash for reuse
        if [ -f "$SUB_DIR/policy.json" ]; then
            POLICY_HASH=$(md5sum "$SUB_DIR/policy.json" | cut -d' ' -f1)
            POLICY_HASHES["$USERNAME"]="$POLICY_HASH"
        fi
    else
        ((FAILED++))
        echo "✗ Failed to process $USERNAME"
    fi
done

echo ""
echo "Processing complete: $PROCESSED new, $SKIPPED already processed, $FAILED failed"
echo ""

# Step 2: Compile LAMMPS with fixes
if [ "$COMPILE_ONLY" = true ] || [ "$RUN_SIMULATIONS" = true ]; then
    echo "=========================================="
    echo "Step 2: Compiling LAMMPS with fixes"
    echo "=========================================="
    echo ""
    
    # Check if LAMMPS path is set
    if [ -z "$LAMMPS_DIR" ]; then
        echo "Error: LAMMPS_DIR environment variable not set"
        echo "Set it to your LAMMPS source directory, e.g.:"
        echo "  export LAMMPS_DIR=/path/to/lammps"
        exit 1
    fi
    
    # Group submissions by policy hash (reuse builds)
    declare -A HASH_TO_BUILD
    
    for SUB_DIR in $SUBMISSION_DIRS; do
        USERNAME=$(basename "$SUB_DIR")
        
        if [ ! -d "$SUB_DIR/generated" ]; then
            echo "Skipping $USERNAME (no generated files)"
            continue
        fi
        
        POLICY_HASH="${POLICY_HASHES[$USERNAME]}"
        
        if [ -z "$POLICY_HASH" ]; then
            # Unique build per submission
            BUILD_DIR="$BUILD_BASE_DIR/${USERNAME}_build"
        else
            # Reuse build if same policy
            if [ -z "${HASH_TO_BUILD[$POLICY_HASH]}" ]; then
                BUILD_DIR="$BUILD_BASE_DIR/policy_${POLICY_HASH:0:8}_build"
                HASH_TO_BUILD["$POLICY_HASH"]="$BUILD_DIR"
                BUILD_DIRS["$USERNAME"]="$BUILD_DIR"
            else
                BUILD_DIR="${HASH_TO_BUILD[$POLICY_HASH]}"
                BUILD_DIRS["$USERNAME"]="$BUILD_DIR"
                echo "Reusing build for $USERNAME (same policy as previous submission)"
                
                # Check if binary already exists (already compiled)
                if [ -f "$BUILD_DIR/lmp_mpi" ]; then
                    echo "  ✓ Already compiled - skipping"
                    continue
                fi
            fi
        fi
        
        BUILD_DIRS["$USERNAME"]="$BUILD_DIR"
        mkdir -p "$BUILD_DIR"
        
        # Check if already compiled
        if [ -f "$BUILD_DIR/lmp_mpi" ]; then
            echo ""
            echo "✓ $USERNAME: Already compiled (binary exists)"
            echo "  Build dir: $BUILD_DIR"
            continue
        fi
        
        echo ""
        echo "Compiling for $USERNAME..."
        echo "  Build dir: $BUILD_DIR"
        
        # Copy generated fix files to build directory
        cp "$SUB_DIR/generated"/*.{h,cpp} "$BUILD_DIR/" 2>/dev/null || {
            echo "Error: No generated files found for $USERNAME"
            continue
        }
        
        # Create isolated LAMMPS build
        LAMMPS_BUILD="$BUILD_DIR/lammps_build"
        if [ ! -d "$LAMMPS_BUILD" ]; then
            echo "  Creating isolated LAMMPS build..."
            cp -r "$LAMMPS_DIR" "$LAMMPS_BUILD"
        fi
        
        # Copy fix files to LAMMPS src
        FIX_NAME=$(basename "$SUB_DIR/generated"/*.h .h)
        cp "$BUILD_DIR"/*.{h,cpp} "$LAMMPS_BUILD/src/"
        
        # Compile
        echo "  Compiling LAMMPS..."
        cd "$LAMMPS_BUILD/src"
        make yes-RIGID yes-MOLECULE 2>&1 | tee "$BUILD_DIR/compile.log"
        make mpi -j$(nproc) 2>&1 | tee -a "$BUILD_DIR/compile.log" || {
            echo "✗ Compilation failed for $USERNAME"
            cd "$SCRIPT_DIR"
            continue
        }
        cd "$SCRIPT_DIR"
        
        # Copy compiled binary to build directory
        cp "$LAMMPS_BUILD/src/lmp_mpi" "$BUILD_DIR/"
        
        echo "✓ Compilation complete for $USERNAME"
    done
fi

# Step 3: Run simulations (if requested)
if [ "$RUN_SIMULATIONS" = true ]; then
    echo ""
    echo "=========================================="
    echo "Step 3: Running Simulations"
    echo "=========================================="
    echo ""
    echo "Note: This step requires manual execution or SLURM setup"
    echo "For each submission, run:"
    echo ""
    
    for SUB_DIR in $SUBMISSION_DIRS; do
        USERNAME=$(basename "$SUB_DIR")
        BUILD_DIR="${BUILD_DIRS[$USERNAME]}"
        
        if [ -z "$BUILD_DIR" ] || [ ! -f "$BUILD_DIR/lmp_mpi" ]; then
            continue
        fi
        
        echo "  $USERNAME:"
        echo "    cd $SUB_DIR/simulation/"
        echo "    $BUILD_DIR/lmp_mpi -in in.*"
        echo ""
    done
fi

# Step 4: Cleanup (if requested)
if [ "$CLEANUP" = true ] && [ "$KEEP_BUILDS" = false ]; then
    echo ""
    echo "=========================================="
    echo "Step 4: Cleanup"
    echo "=========================================="
    echo ""
    echo "Removing build directories..."
    rm -rf "$BUILD_BASE_DIR"
    echo "✓ Cleanup complete"
fi

echo ""
echo "=========================================="
echo "Batch Processing Complete!"
echo "=========================================="
echo "Processed: $PROCESSED"
echo "Failed: $FAILED"
if [ "$KEEP_BUILDS" = true ]; then
    echo "Builds kept in: $BUILD_BASE_DIR"
fi
echo ""

