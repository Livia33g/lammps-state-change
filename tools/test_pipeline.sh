#!/bin/bash
# Test script for the dimer case set problem pipeline
# This tests: policy.json -> C++ fix generation + problem.json -> LAMMPS system generation
#
# Usage:
#   tools/test_pipeline.sh [problem_id] [test_name]
#   tools/test_pipeline.sh problem-001-dimer-ksat baseline_test
#
# Outputs go to: tests/{problem_id}/{test_name}/

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$SCRIPT_DIR"

# Parse arguments
PROBLEM_ID="${1:-problem-001-dimer-ksat}"
TEST_NAME="${2:-baseline_test}"

# Set up paths
PROBLEM_DIR="problems/$PROBLEM_ID"
POLICY_FILE="$PROBLEM_DIR/baseline_policy.json"
PARAMS_FILE="$PROBLEM_DIR/baseline_params.json"
PROBLEM_FILE="$PROBLEM_DIR/problem.json"
TEST_OUTPUT_DIR="tests/$PROBLEM_ID/$TEST_NAME"

# Validate inputs exist
if [ ! -f "$PROBLEM_FILE" ]; then
    echo "Error: Problem file not found: $PROBLEM_FILE"
    exit 1
fi

if [ ! -f "$POLICY_FILE" ]; then
    echo "Error: Policy file not found: $POLICY_FILE"
    exit 1
fi

if [ ! -f "$PARAMS_FILE" ]; then
    echo "Error: Params file not found: $PARAMS_FILE"
    exit 1
fi

echo "=========================================="
echo "Testing Pipeline for $PROBLEM_ID"
echo "Test Name: $TEST_NAME"
echo "=========================================="
echo ""

echo "Step 1: Creating organized test output directory..."
mkdir -p "$TEST_OUTPUT_DIR/generated"
mkdir -p "$TEST_OUTPUT_DIR/simulation"

echo ""
echo "Step 2: Generating C++ fix from policy..."
echo "  Policy: $POLICY_FILE"
echo "  Output: $TEST_OUTPUT_DIR/generated/"
python3 core/generators/generate_fix_from_policy.py \
    "$POLICY_FILE" \
    "$TEST_OUTPUT_DIR/generated/"

echo ""
echo "Step 3: Generating LAMMPS system files..."
echo "  Problem: $PROBLEM_FILE"
echo "  Policy: $POLICY_FILE"
echo "  Params: $PARAMS_FILE"
echo "  Output: $TEST_OUTPUT_DIR/simulation/"

# Extract params from baseline_params.json (it has nested structure)
# Create a simple params file with just the values
TEMP_PARAMS="$TEST_OUTPUT_DIR/temp_params.json"
python3 << EOF
import json
with open("$PARAMS_FILE") as f:
    baseline = json.load(f)
    
# Extract interaction parameters
params = {}
if "interaction_parameters" in baseline:
    for key, val in baseline["interaction_parameters"].items():
        if isinstance(val, dict) and "value" in val:
            params[key] = val["value"]
        else:
            params[key] = val

if "state_change_parameters" in baseline:
    for key, val in baseline["state_change_parameters"].items():
        if isinstance(val, dict) and "value" in val:
            params[key] = val["value"]
        else:
            params[key] = val

with open("$TEMP_PARAMS", "w") as f:
    json.dump(params, f, indent=2)
EOF

python3 core/generators/generate_system_from_problem.py \
    --problem "$PROBLEM_FILE" \
    --policy "$POLICY_FILE" \
    --params "$TEMP_PARAMS" \
    --output "$TEST_OUTPUT_DIR/simulation/"

# Clean up temp file
rm -f "$TEMP_PARAMS"

echo ""
echo "=========================================="
echo "Pipeline Test Complete!"
echo "=========================================="
echo ""
echo "Generated files organized in:"
echo "  $TEST_OUTPUT_DIR/"
echo "    ├── generated/          # C++ fix files"
echo "    │   ├── fix_state_change_*.h"
echo "    │   └── fix_state_change_*.cpp"
echo "    └── simulation/         # LAMMPS input files"
echo "        ├── data.*"
echo "        └── in.*"
echo ""
echo "Next steps:"
echo "  1. Review generated C++ files: $TEST_OUTPUT_DIR/generated/"
echo "  2. Review LAMMPS input files: $TEST_OUTPUT_DIR/simulation/"
echo "  3. If you have LAMMPS compiled with the fix:"
echo "     cd $TEST_OUTPUT_DIR/simulation/"
echo "     lmp_mpi -in in.*"
echo ""

