#!/bin/bash
# Process a competitor submission: validate, generate C++ fix, generate simulation files
#
# Usage:
#   tools/process_submission.sh submissions/problem-001-dimer-ksat/competitor_username/
#
# This script:
#   1. Validates the submission (JSON format, schema compliance)
#   2. Generates C++ fix from policy.json
#   3. Generates LAMMPS system files from problem.json + policy.json + params.json
#   4. Creates organized directory structure

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$SCRIPT_DIR"

# Parse arguments
SUBMISSION_DIR="${1}"

if [ -z "$SUBMISSION_DIR" ]; then
    echo "Usage: tools/process_submission.sh submissions/problem-XXX/competitor_username/"
    exit 1
fi

# Normalize path
SUBMISSION_DIR="$(cd "$(dirname "$SUBMISSION_DIR")/$(basename "$SUBMISSION_DIR")" && pwd)"
SUBMISSION_DIR="${SUBMISSION_DIR#$SCRIPT_DIR/}"

# Extract problem ID from path
PROBLEM_ID=$(echo "$SUBMISSION_DIR" | cut -d'/' -f2)
COMPETITOR_NAME=$(basename "$SUBMISSION_DIR")

# Set up paths
POLICY_FILE="$SUBMISSION_DIR/policy.json"
PARAMS_FILE="$SUBMISSION_DIR/params.json"
SUBMISSION_META="$SUBMISSION_DIR/submission.json"
PROBLEM_DIR="problems/$PROBLEM_ID"
PROBLEM_FILE="$PROBLEM_DIR/problem.json"

# Validate inputs exist
if [ ! -f "$POLICY_FILE" ]; then
    echo "Error: Policy file not found: $POLICY_FILE"
    exit 1
fi

if [ ! -f "$PROBLEM_FILE" ]; then
    echo "Error: Problem file not found: $PROBLEM_FILE"
    echo "Expected: $PROBLEM_FILE"
    exit 1
fi

# Check for params file (optional)
if [ ! -f "$PARAMS_FILE" ]; then
    echo "Warning: No params.json found, using defaults from problem.json"
    PARAMS_FILE="$PROBLEM_DIR/baseline_params.json"
fi

echo "=========================================="
echo "Processing Submission"
echo "=========================================="
echo "Problem:    $PROBLEM_ID"
echo "Competitor: $COMPETITOR_NAME"
echo "Directory:  $SUBMISSION_DIR"
echo ""

# Create output directories
mkdir -p "$SUBMISSION_DIR/generated"
mkdir -p "$SUBMISSION_DIR/simulation"
mkdir -p "$SUBMISSION_DIR/results"

echo "Step 1: Validating submission..."
if [ -f "tools/validate_submission.py" ]; then
    python3 tools/validate_submission.py "$SUBMISSION_DIR" || {
        echo "Warning: Validation failed or validator not available"
    }
else
    echo "  (Skipping - validator not found)"
fi

echo ""
echo "Step 2: Generating C++ fix from policy..."
echo "  Policy: $POLICY_FILE"
echo "  Output: $SUBMISSION_DIR/generated/"
python3 core/generators/generate_fix_from_policy.py \
    "$POLICY_FILE" \
    "$SUBMISSION_DIR/generated/"

echo ""
echo "Step 3: Generating LAMMPS system files..."
echo "  Problem: $PROBLEM_FILE"
echo "  Policy:  $POLICY_FILE"
echo "  Params:  $PARAMS_FILE"
echo "  Output:  $SUBMISSION_DIR/simulation/"

# Extract params if needed (handle nested structure)
TEMP_PARAMS="$SUBMISSION_DIR/temp_params.json"
python3 << EOF
import json
import sys

try:
    with open("$PARAMS_FILE") as f:
        baseline = json.load(f)
    
    # Extract interaction parameters (handle nested structure)
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
    
    # If it's already a flat params file, use as-is
    if not params and isinstance(baseline, dict):
        # Check if it looks like a flat params file
        if any(k.startswith("morse_") or k.startswith("contact_") for k in baseline.keys()):
            params = baseline
    
    with open("$TEMP_PARAMS", "w") as f:
        json.dump(params, f, indent=2)
except Exception as e:
    print(f"Error processing params: {e}", file=sys.stderr)
    # Create empty params file
    with open("$TEMP_PARAMS", "w") as f:
        json.dump({}, f)
EOF

python3 core/generators/generate_system_from_problem.py \
    --problem "$PROBLEM_FILE" \
    --policy "$POLICY_FILE" \
    --params "$TEMP_PARAMS" \
    --output "$SUBMISSION_DIR/simulation/"

# Clean up temp file
rm -f "$TEMP_PARAMS"

echo ""
echo "=========================================="
echo "Submission Processing Complete!"
echo "=========================================="
echo ""
echo "Generated files in: $SUBMISSION_DIR/"
echo "  ├── generated/          # C++ fix files"
echo "  ├── simulation/         # LAMMPS input files"
echo "  └── results/            # (empty, for simulation outputs)"
echo ""
echo "Next steps:"
echo "  1. Review generated files"
echo "  2. Compile LAMMPS with fix from: $SUBMISSION_DIR/generated/"
echo "  3. Run simulation from: $SUBMISSION_DIR/simulation/"
echo "  4. Analyze results and update leaderboard"
echo ""

