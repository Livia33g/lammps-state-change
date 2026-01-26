#!/bin/bash
# Process all submissions for a given problem
#
# Usage:
#   tools/process_all_submissions.sh problem-001-dimer-ksat
#
# This processes all competitor submissions in submissions/{problem_id}/

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$SCRIPT_DIR"

PROBLEM_ID="${1}"

if [ -z "$PROBLEM_ID" ]; then
    echo "Usage: tools/process_all_submissions.sh problem-001-dimer-ksat"
    exit 1
fi

SUBMISSIONS_DIR="submissions/$PROBLEM_ID"

if [ ! -d "$SUBMISSIONS_DIR" ]; then
    echo "Error: Submissions directory not found: $SUBMISSIONS_DIR"
    exit 1
fi

echo "=========================================="
echo "Processing All Submissions"
echo "=========================================="
echo "Problem: $PROBLEM_ID"
echo "Directory: $SUBMISSIONS_DIR"
echo ""

# Find all competitor directories
COMPETITORS=$(find "$SUBMISSIONS_DIR" -mindepth 1 -maxdepth 1 -type d ! -name ".*" | sort)

if [ -z "$COMPETITORS" ]; then
    echo "No submissions found in $SUBMISSIONS_DIR"
    exit 0
fi

COMPETITOR_COUNT=$(echo "$COMPETITORS" | wc -l)
echo "Found $COMPETITOR_COUNT submission(s)"
echo ""

# Process each submission
SUCCESS_COUNT=0
FAIL_COUNT=0

for COMPETITOR_DIR in $COMPETITORS; do
    COMPETITOR_NAME=$(basename "$COMPETITOR_DIR")
    echo "----------------------------------------"
    echo "Processing: $COMPETITOR_NAME"
    echo "----------------------------------------"
    
    if tools/process_submission.sh "$COMPETITOR_DIR"; then
        ((SUCCESS_COUNT++))
        echo "✓ Successfully processed $COMPETITOR_NAME"
    else
        ((FAIL_COUNT++))
        echo "✗ Failed to process $COMPETITOR_NAME"
    fi
    echo ""
done

echo "=========================================="
echo "Processing Complete"
echo "=========================================="
echo "Successful: $SUCCESS_COUNT"
echo "Failed:     $FAIL_COUNT"
echo "Total:      $COMPETITOR_COUNT"
echo ""

