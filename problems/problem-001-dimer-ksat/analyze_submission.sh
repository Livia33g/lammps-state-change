#!/bin/bash
# Problem-specific analysis script for problem-001-dimer-ksat
#
# This script analyzes simulation results for the dimer k-SAT problem.
# It generates timeseries data and computes metrics specific to this problem.
#
# Usage: Called by tools/analyze_submission_results.sh
# Arguments are passed via environment or command line

set -e

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --submission-dir)
            SUBMISSION_DIR="$2"
            shift 2
            ;;
        --results-dir)
            RESULTS_DIR="$2"
            shift 2
            ;;
        --analysis-dir)
            ANALYSIS_DIR="$2"
            shift 2
            ;;
        --dump)
            DUMP_FILE="$2"
            shift 2
            ;;
        --thermo)
            THERMO_FILE="$2"
            shift 2
            ;;
        --events)
            EVENTS_FILE="$2"
            shift 2
            ;;
        --username)
            USERNAME="$2"
            shift 2
            ;;
        --problem-id)
            PROBLEM_ID="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# Get script directory to find problem.json
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROBLEM_FILE="$SCRIPT_DIR/problem.json"

if [ ! -f "$PROBLEM_FILE" ]; then
    echo "Error: problem.json not found in $SCRIPT_DIR"
    exit 1
fi

# Load analysis parameters from problem.json
ANALYSIS_PARAMS=$(python3 << EOF
import json
with open("$PROBLEM_FILE") as f:
    problem = json.load(f)

decoding = problem.get("decoding", {})
analysis_args = decoding.get("analysis_args", {})

# Extract parameters
site_types = analysis_args.get("site_types", [2, 3, 4])
bond_cutoff = analysis_args.get("bond_cutoff", 0.7)
yield_mode = analysis_args.get("yield_mode", "species_fraction")
species_label = analysis_args.get("species_label", "4")
yield_threshold = decoding.get("threshold", 0.6)

print(f"SITE_TYPES={' '.join(map(str, site_types))}")
print(f"BOND_CUTOFF={bond_cutoff}")
print(f"YIELD_MODE={yield_mode}")
print(f"SPECIES_LABEL={species_label}")
print(f"YIELD_THRESHOLD={yield_threshold}")
EOF
)

# Export as environment variables
eval "$ANALYSIS_PARAMS"

echo "Analysis parameters:"
echo "  Site types: $SITE_TYPES"
echo "  Bond cutoff: $BOND_CUTOFF"
echo "  Yield mode: $YIELD_MODE"
echo "  Species label: $SPECIES_LABEL"
echo "  Yield threshold: $YIELD_THRESHOLD"
echo ""

# Step 1: Generate timeseries CSV
echo "Step 1: Generating timeseries CSV from trajectory..."
TIMESERIES_CSV="$ANALYSIS_DIR/timeseries.csv"

python3 "$SCRIPT_DIR/../../core/analysis/analyze_trajectory_target_yield_and_work.py" \
    --dump "$DUMP_FILE" \
    --thermo "${THERMO_FILE:-$RESULTS_DIR/lammps_stdout.log}" \
    --events "${EVENTS_FILE:-$RESULTS_DIR/stderr.err}" \
    --site-types $SITE_TYPES \
    --bond-cutoff "$BOND_CUTOFF" \
    --yield-mode "$YIELD_MODE" \
    --species-label "$SPECIES_LABEL" \
    --out "$ANALYSIS_DIR/timeseries" \
    > "$ANALYSIS_DIR/analysis.log" 2>&1 || {
    echo "Warning: Check analysis.log for details"
    if [ ! -f "$ANALYSIS_DIR/timeseries.csv" ]; then
        echo "Error: Failed to generate timeseries CSV"
        exit 1
    fi
}

# Use the generated CSV
if [ -f "$ANALYSIS_DIR/timeseries.csv" ]; then
    TIMESERIES_CSV="$ANALYSIS_DIR/timeseries.csv"
fi

# Step 2: Compute metrics and generate leaderboard row
echo ""
echo "Step 2: Computing metrics and generating leaderboard row..."
LEADERBOARD_CSV="$ANALYSIS_DIR/leaderboard_row.csv"

python3 "$SCRIPT_DIR/../../core/benchmark/score_policy_from_timeseries.py" \
    --csv "$TIMESERIES_CSV" \
    --yield-threshold "$YIELD_THRESHOLD" \
    --out "$LEADERBOARD_CSV"

# Step 3: Add metadata to leaderboard row
echo ""
echo "Step 3: Adding metadata to leaderboard row..."
python3 << EOF
import csv
import json
from pathlib import Path

# Read existing leaderboard row
leaderboard_path = Path("$LEADERBOARD_CSV")
with open(leaderboard_path) as f:
    reader = csv.DictReader(f)
    row = next(reader)

# Add metadata
row["problem_id"] = "$PROBLEM_ID"
row["username"] = "$USERNAME"

# Try to get info from submission.json
submission_json = Path("$SUBMISSION_DIR/submission.json")
if submission_json.exists():
    with open(submission_json) as f:
        sub_data = json.load(f)
        row["team_name"] = sub_data.get("team_name", "$USERNAME")
        row["submission_date"] = sub_data.get("submission_date", "")
        row["policy_version"] = sub_data.get("policy_version", "")

# Try to get policy name from policy.json
policy_json = Path("$SUBMISSION_DIR/policy.json")
if policy_json.exists():
    with open(policy_json) as f:
        policy_data = json.load(f)
        row["policy_name"] = policy_data.get("policy_name", "")

# Write updated row with metadata first
fieldnames = ["problem_id", "username", "team_name", "policy_name", "submission_date", "policy_version"] + \
             [k for k in row.keys() if k not in ["problem_id", "username", "team_name", "policy_name", "submission_date", "policy_version"]]

with open(leaderboard_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerow(row)
EOF

echo "✓ Analysis complete!"
echo "  Timeseries: $TIMESERIES_CSV"
echo "  Leaderboard row: $LEADERBOARD_CSV"

