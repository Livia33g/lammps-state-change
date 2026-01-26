#!/bin/bash
# Analyze simulation results for a submission and generate metrics
#
# Usage:
#   tools/analyze_submission_results.sh submissions/problem-001-dimer-ksat/username/ [--force]
#
# This analyzes the trajectory and generates:
#   - Timeseries CSV
#   - Leaderboard row CSV
#   - Results summary
#
# Options:
#   --force    Force re-analysis even if results already exist

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$SCRIPT_DIR"

SUBMISSION_DIR="${1}"
FORCE=false

# Parse options
shift || true
while [[ $# -gt 0 ]]; do
    case $1 in
        --force)
            FORCE=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

if [ -z "$SUBMISSION_DIR" ]; then
    echo "Usage: tools/analyze_submission_results.sh submissions/problem-XXX/username/ [--force]"
    exit 1
fi

# Normalize path
SUBMISSION_DIR="$(cd "$(dirname "$SUBMISSION_DIR")/$(basename "$SUBMISSION_DIR")" && pwd)"
SUBMISSION_DIR="${SUBMISSION_DIR#$SCRIPT_DIR/}"

# Extract problem ID and username
PROBLEM_ID=$(echo "$SUBMISSION_DIR" | cut -d'/' -f2)
USERNAME=$(basename "$SUBMISSION_DIR")

# Set up paths
RESULTS_DIR="$SUBMISSION_DIR/results"
ANALYSIS_DIR="$RESULTS_DIR/analysis"
DUMP_FILE=$(find "$RESULTS_DIR" -name "*.lammpstrj" | head -1)
THERMO_FILE=$(find "$RESULTS_DIR" -name "lammps_stdout.log" -o -name "*.log" | head -1)
EVENTS_FILE=$(find "$RESULTS_DIR" -name "*.err" -o -name "slurm*.err" | head -1)

# Check if results exist
if [ ! -d "$RESULTS_DIR" ] || [ -z "$DUMP_FILE" ]; then
    echo "Error: No simulation results found in $RESULTS_DIR"
    echo "Expected: dump.*.lammpstrj file"
    exit 1
fi

# Check if already analyzed
if [ -f "$ANALYSIS_DIR/leaderboard_row.csv" ] && [ "$FORCE" = false ]; then
    echo "✓ Analysis already complete (leaderboard_row.csv exists)"
    echo "  Skipping re-analysis. Use --force to re-analyze."
    echo ""
    echo "Existing analysis:"
    echo "  $ANALYSIS_DIR/leaderboard_row.csv"
    exit 0
fi

if [ "$FORCE" = true ]; then
    echo "⚠ Force mode: Re-analyzing even if results exist"
fi

echo "=========================================="
echo "Analyzing Submission Results"
echo "=========================================="
echo "Problem:    $PROBLEM_ID"
echo "Username:  $USERNAME"
echo "Results:    $RESULTS_DIR"
echo ""

# Check for problem-specific analysis script
PROBLEM_DIR="problems/$PROBLEM_ID"
ANALYSIS_SCRIPT="$PROBLEM_DIR/analyze_submission.sh"

if [ ! -f "$ANALYSIS_SCRIPT" ]; then
    echo "Error: Problem-specific analysis script not found: $ANALYSIS_SCRIPT"
    echo "Each problem must have its own analyze_submission.sh script"
    exit 1
fi

# Create analysis directory
mkdir -p "$ANALYSIS_DIR"

# Run problem-specific analysis script
echo "Running problem-specific analysis script..."
echo "Script: $ANALYSIS_SCRIPT"
echo ""

# The analysis script should:
# 1. Generate timeseries CSV from trajectory
# 2. Generate leaderboard_row.csv with metrics
# 3. Output to $ANALYSIS_DIR/
#
# Expected outputs:
# - $ANALYSIS_DIR/timeseries.csv
# - $ANALYSIS_DIR/leaderboard_row.csv

bash "$ANALYSIS_SCRIPT" \
    --submission-dir "$SUBMISSION_DIR" \
    --results-dir "$RESULTS_DIR" \
    --analysis-dir "$ANALYSIS_DIR" \
    --dump "$DUMP_FILE" \
    --thermo "${THERMO_FILE:-$RESULTS_DIR/lammps_stdout.log}" \
    --events "${EVENTS_FILE:-$RESULTS_DIR/stderr.err}" \
    --username "$USERNAME" \
    --problem-id "$PROBLEM_ID"

# Verify outputs
TIMESERIES_CSV="$ANALYSIS_DIR/timeseries.csv"
LEADERBOARD_CSV="$ANALYSIS_DIR/leaderboard_row.csv"

if [ ! -f "$LEADERBOARD_CSV" ]; then
    echo "Error: Analysis script did not produce leaderboard_row.csv"
    exit 1
fi

# Add metadata to leaderboard row
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

# Try to get policy name from submission
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

# Write updated row
with open(leaderboard_path, "w", newline="") as f:
    fieldnames = ["problem_id", "username", "team_name", "policy_name", "submission_date", "policy_version"] + \
                 [k for k in row.keys() if k not in ["problem_id", "username", "team_name", "policy_name", "submission_date", "policy_version"]]
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerow(row)
EOF

echo ""
echo "=========================================="
echo "Analysis Complete!"
echo "=========================================="
echo "Results in: $ANALYSIS_DIR"
echo "  - timeseries.csv      # Full timeseries data"
echo "  - leaderboard_row.csv # Single row for leaderboard"
echo ""
echo "Next step: Aggregate into problem leaderboard"
echo "  tools/update_leaderboard.sh $PROBLEM_ID"
echo ""

