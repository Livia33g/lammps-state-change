#!/bin/bash
# Efficient submission processing with smart resource management
#
# This is a wrapper that provides different processing strategies:
# 1. One-by-one (minimal memory, sequential)
# 2. Batch with shared builds (efficient, reuses compilations)
# 3. Isolated builds (safe, no conflicts)
#
# Usage:
#   tools/process_efficiently.sh problem-001-dimer-ksat [strategy]
#
# Strategies:
#   one-by-one    - Process sequentially, clean up after each
#   batch-shared  - Process all, reuse builds for same policies (default)
#   batch-isolated - Process all, separate build per submission

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$SCRIPT_DIR"

PROBLEM_ID="${1}"
STRATEGY="${2:-batch-shared}"

if [ -z "$PROBLEM_ID" ]; then
    echo "Usage: tools/process_efficiently.sh problem-001-dimer-ksat [strategy]"
    echo ""
    echo "Strategies:"
    echo "  one-by-one      - Process one at a time, clean up immediately"
    echo "  batch-shared    - Process all, reuse builds for identical policies"
    echo "  batch-isolated  - Process all, separate build per submission"
    echo ""
    exit 1
fi

SUBMISSIONS_DIR="submissions/$PROBLEM_ID"

if [ ! -d "$SUBMISSIONS_DIR" ]; then
    echo "Error: Submissions directory not found: $SUBMISSIONS_DIR"
    exit 1
fi

# Find all submissions
SUBMISSION_DIRS=$(find "$SUBMISSIONS_DIR" -mindepth 1 -maxdepth 1 -type d ! -name ".*" | sort)
COUNT=$(echo "$SUBMISSION_DIRS" | wc -w)

echo "=========================================="
echo "Efficient Processing"
echo "=========================================="
echo "Problem: $PROBLEM_ID"
echo "Strategy: $STRATEGY"
echo "Submissions: $COUNT"
echo ""

case "$STRATEGY" in
    one-by-one)
        echo "Strategy: One-by-one (minimal resource usage)"
        echo ""
        
        for SUB_DIR in $SUBMISSION_DIRS; do
            USERNAME=$(basename "$SUB_DIR")
            TEMP_BUILD=".temp_build_$$"
            
            echo "----------------------------------------"
            echo "Processing: $USERNAME"
            echo "----------------------------------------"
            
            # Process submission
            tools/process_submission.sh "$SUB_DIR"
            
            # Compile in temporary directory
            if [ -n "$LAMMPS_DIR" ] && [ -d "$SUB_DIR/generated" ]; then
                mkdir -p "$TEMP_BUILD"
                cp "$SUB_DIR/generated"/*.{h,cpp} "$TEMP_BUILD/"
                
                # Quick compile (simplified - you'd integrate with your LAMMPS build)
                echo "  Compiling (temporary build)..."
                # ... compilation logic ...
                
                # Clean up immediately
                rm -rf "$TEMP_BUILD"
            fi
            
            echo "✓ Complete: $USERNAME"
            echo ""
        done
        ;;
        
    batch-shared)
        echo "Strategy: Batch with shared builds (efficient)"
        echo "Using compile_and_run_batch.sh with policy hash reuse"
        echo ""
        
        tools/compile_and_run_batch.sh "$PROBLEM_ID" --compile-only --cleanup
        ;;
        
    batch-isolated)
        echo "Strategy: Batch with isolated builds (safe)"
        echo "Each submission gets its own build directory"
        echo ""
        
        tools/compile_and_run_batch.sh "$PROBLEM_ID" --compile-only --keep-builds
        ;;
        
    *)
        echo "Error: Unknown strategy: $STRATEGY"
        exit 1
        ;;
esac

echo ""
echo "=========================================="
echo "Processing Complete!"
echo "=========================================="
echo ""

