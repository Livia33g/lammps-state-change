#!/bin/bash
#
# Process all new submissions in submissions-inbox/
#
# This script automatically validates and submits all pending submissions
# to the cluster for evaluation.
#
# Usage:
#   ./tools/process_new_submissions.sh
#
# Directory structure expected:
#   submissions-inbox/
#     ├── username1_problem-001/
#     │   ├── policy.json
#     │   ├── params.json (optional)
#     │   └── submission.json
#     └── username2_problem-001/
#         └── ...
#

set -e

INBOX_DIR="submissions-inbox"
PRIVATE_DIR="submissions-private"
TOOLS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$TOOLS_DIR")"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "========================================="
echo "  Submission Batch Processor"
echo "========================================="
echo ""

# Check if inbox exists
if [ ! -d "$INBOX_DIR" ]; then
    echo -e "${YELLOW}Warning: $INBOX_DIR does not exist${NC}"
    echo "Creating inbox directory..."
    mkdir -p "$INBOX_DIR"
    echo "Place submissions in $INBOX_DIR/<username>_<problem-id>/"
    exit 0
fi

# Check if there are any submissions
submission_count=$(find "$INBOX_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)

if [ "$submission_count" -eq 0 ]; then
    echo "No submissions found in $INBOX_DIR"
    echo "Nothing to process."
    exit 0
fi

echo "Found $submission_count submission(s) to process"
echo ""

# Process each submission
processed=0
failed=0

for submission_dir in "$INBOX_DIR"/*/ ; do
    [ -e "$submission_dir" ] || continue

    dirname=$(basename "$submission_dir")
    echo "----------------------------------------"
    echo "Processing: $dirname"
    echo "----------------------------------------"

    # Extract username and problem from directory name
    # Expected format: username_problem-id
    if [[ "$dirname" =~ ^(.+)_(problem-.+)$ ]]; then
        username="${BASH_MATCH[1]}"
        problem_id="${BASH_MATCH[2]}"
    else
        echo -e "${RED}✗ Invalid directory name format: $dirname${NC}"
        echo "  Expected: username_problem-id"
        ((failed++))
        continue
    fi

    echo "  Username: $username"
    echo "  Problem:  $problem_id"

    # Check if submission.json exists
    if [ ! -f "$submission_dir/submission.json" ]; then
        echo -e "${RED}✗ Missing submission.json${NC}"
        ((failed++))
        continue
    fi

    # Extract problem_id from submission.json and verify it matches
    json_problem=$(python3 -c "import json; print(json.load(open('$submission_dir/submission.json'))['problem_id'])" 2>/dev/null || echo "")

    if [ "$json_problem" != "$problem_id" ]; then
        echo -e "${YELLOW}Warning: Problem ID mismatch${NC}"
        echo "  Directory name: $problem_id"
        echo "  submission.json: $json_problem"
        echo "  Using: $json_problem"
        problem_id="$json_problem"
    fi

    # Validate submission
    echo "  Validating..."
    if python3 "$TOOLS_DIR/validate_submission.py" "$submission_dir" > /dev/null 2>&1; then
        echo -e "  ${GREEN}✓ Validation passed${NC}"
    else
        echo -e "  ${RED}✗ Validation failed${NC}"
        echo ""
        echo "  Error details:"
        python3 "$TOOLS_DIR/validate_submission.py" "$submission_dir" 2>&1 | sed 's/^/    /'
        echo ""
        echo "  Action: Email participant with error message"
        ((failed++))
        continue
    fi

    # Create timestamp for this submission
    timestamp=$(date +%Y-%m-%d_%H%M%S)
    dest_dir="$PRIVATE_DIR/$problem_id/${username}_${timestamp}"

    # Create destination directory
    mkdir -p "$dest_dir"

    # Move submission to private repo
    echo "  Moving to: $dest_dir"
    cp -r "$submission_dir"/* "$dest_dir/"

    # Submit to cluster
    echo "  Submitting to cluster..."

    # Check if we're on a SLURM system
    if command -v sbatch &> /dev/null; then
        job_id=$(sbatch "$TOOLS_DIR/evaluate_submission.slurm" "$problem_id" "${username}_${timestamp}" | awk '{print $4}')
        echo -e "  ${GREEN}✓ Submitted (Job ID: $job_id)${NC}"
    else
        echo -e "  ${YELLOW}⚠ SLURM not available - manual evaluation required${NC}"
        echo "  Run: bash $TOOLS_DIR/evaluate_submission.slurm $problem_id ${username}_${timestamp}"
    fi

    # Archive the inbox submission
    archive_dir="$INBOX_DIR/.processed"
    mkdir -p "$archive_dir"
    mv "$submission_dir" "$archive_dir/${dirname}_${timestamp}"

    echo -e "  ${GREEN}✓ Processing complete${NC}"
    ((processed++))
    echo ""
done

echo "========================================="
echo "  Summary"
echo "========================================="
echo "  Processed: $processed"
echo "  Failed:    $failed"
echo ""

if [ "$processed" -gt 0 ]; then
    echo "Next steps:"
    echo "  1. Monitor cluster jobs (squeue)"
    echo "  2. Check results when jobs complete"
    echo "  3. Update leaderboard with new scores"
    echo "  4. Email participants with results"
fi

exit 0
