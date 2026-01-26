#!/bin/bash
# Aggregate all submission results for a problem into a leaderboard
#
# Usage:
#   tools/update_leaderboard.sh problem-001-dimer-ksat
#
# This collects all leaderboard_row.csv files from submissions and aggregates
# them into problems/{problem_id}/leaderboard.csv

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$SCRIPT_DIR"

PROBLEM_ID="${1}"

if [ -z "$PROBLEM_ID" ]; then
    echo "Usage: tools/update_leaderboard.sh problem-001-dimer-ksat"
    exit 1
fi

SUBMISSIONS_DIR="submissions/$PROBLEM_ID"
LEADERBOARD_DIR="problems/$PROBLEM_ID"
LEADERBOARD_CSV="$LEADERBOARD_DIR/leaderboard.csv"

if [ ! -d "$SUBMISSIONS_DIR" ]; then
    echo "Error: Submissions directory not found: $SUBMISSIONS_DIR"
    exit 1
fi

echo "=========================================="
echo "Updating Leaderboard"
echo "=========================================="
echo "Problem: $PROBLEM_ID"
echo "Output:  $LEADERBOARD_CSV"
echo ""

# Find all leaderboard row files
LEADERBOARD_FILES=$(find "$SUBMISSIONS_DIR" -name "leaderboard_row.csv" | sort)

if [ -z "$LEADERBOARD_FILES" ]; then
    echo "Warning: No leaderboard_row.csv files found in $SUBMISSIONS_DIR"
    echo "Make sure to run analyze_submission_results.sh for each submission first"
    exit 1
fi

COUNT=$(echo "$LEADERBOARD_FILES" | wc -l)
echo "Found $COUNT submission result(s)"
echo ""

# Create leaderboard directory
mkdir -p "$LEADERBOARD_DIR"

# Aggregate using the existing aggregate_leaderboard.py
TEMP_DIR=$(mktemp -d -t leaderboard_XXXXXX)
trap "rm -rf $TEMP_DIR" EXIT

# Copy all leaderboard rows to temp directory with unique names
COUNTER=1
for FILE in $LEADERBOARD_FILES; do
    USERNAME=$(echo "$FILE" | sed -n "s|.*submissions/$PROBLEM_ID/\([^/]*\)/.*|\1|p")
    cp "$FILE" "$TEMP_DIR/${USERNAME}_${COUNTER}.leaderboard.csv"
    ((COUNTER++))
done

# Aggregate
echo "Aggregating results..."
python3 core/benchmark/aggregate_leaderboard.py \
    --dir "$TEMP_DIR" \
    --out "$LEADERBOARD_CSV" \
    --sort-by "work_per_yield" \
    --desc

echo ""
echo "=========================================="
echo "Leaderboard Updated!"
echo "=========================================="
echo "Leaderboard: $LEADERBOARD_CSV"
echo ""
echo "Top entries:"
head -6 "$LEADERBOARD_CSV" | column -t -s, 2>/dev/null || head -6 "$LEADERBOARD_CSV"
echo ""

