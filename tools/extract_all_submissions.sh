#!/bin/bash
#
# Batch extract all submissions from a directory
#
# Usage:
#   1. Save all email attachments to ~/Downloads/submissions/
#      Structure:
#        ~/Downloads/submissions/
#          ├── alice/
#          │   ├── policy.json
#          │   ├── submission.json
#          │   └── params.json
#          ├── diana/
#          │   ├── policy.json
#          │   └── submission.json
#          └── bob/
#
#   2. Run: ./tools/extract_all_submissions.sh ~/Downloads/submissions/ problem-001-dimer-ksat
#
#   3. All submissions extracted automatically!

set -e

SUBMISSIONS_DIR="$1"
PROBLEM_ID="$2"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Validate arguments
if [ -z "$SUBMISSIONS_DIR" ] || [ -z "$PROBLEM_ID" ]; then
    echo "Usage: $0 <submissions_directory> <problem_id>"
    echo ""
    echo "Example:"
    echo "  $0 ~/Downloads/submissions problem-001-dimer-ksat"
    echo ""
    echo "Expected directory structure:"
    echo "  submissions_directory/"
    echo "    ├── alice/            # Username from email/folder name"
    echo "    │   ├── policy.json"
    echo "    │   ├── submission.json"
    echo "    │   └── params.json (optional)"
    echo "    └── bob/"
    echo "        └── ..."
    exit 1
fi

if [ ! -d "$SUBMISSIONS_DIR" ]; then
    echo -e "${RED}Error: Directory not found: $SUBMISSIONS_DIR${NC}"
    exit 1
fi

# Configuration
PRIVATE_REPO="${SUBMISSIONS_PRIVATE_REPO:-$HOME/lammps-state-change-private}"
INBOX_DIR="$PRIVATE_REPO/submissions-inbox"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo -e "${BLUE}=========================================${NC}"
echo -e "${BLUE}  Batch Extract Submissions${NC}"
echo -e "${BLUE}=========================================${NC}"
echo "  Source dir: $SUBMISSIONS_DIR"
echo "  Problem ID: $PROBLEM_ID"
echo "  Target: $INBOX_DIR"
echo ""

# Find all subdirectories (each is a username)
submission_count=0
success_count=0
fail_count=0

for user_dir in "$SUBMISSIONS_DIR"/*/ ; do
    [ -e "$user_dir" ] || continue

    username=$(basename "$user_dir")

    # Skip hidden directories
    if [[ "$username" == .* ]]; then
        continue
    fi

    ((submission_count++))

    echo -e "${YELLOW}----------------------------------------${NC}"
    echo "Processing: $username"
    echo -e "${YELLOW}----------------------------------------${NC}"

    # Check for required files
    if [ ! -f "$user_dir/policy.json" ]; then
        echo -e "  ${RED}✗ Missing policy.json${NC}"
        ((fail_count++))
        continue
    fi

    if [ ! -f "$user_dir/submission.json" ]; then
        echo -e "  ${RED}✗ Missing submission.json${NC}"
        ((fail_count++))
        continue
    fi

    # Create destination directory
    dest_dir="$INBOX_DIR/${username}_${PROBLEM_ID}"
    mkdir -p "$dest_dir"

    # Copy files
    cp "$user_dir/policy.json" "$dest_dir/"
    cp "$user_dir/submission.json" "$dest_dir/"

    if [ -f "$user_dir/params.json" ]; then
        cp "$user_dir/params.json" "$dest_dir/"
        echo -e "  ${GREEN}✓ Extracted (with params.json)${NC}"
    else
        echo -e "  ${GREEN}✓ Extracted (no params)${NC}"
    fi

    ((success_count++))
done

echo ""
echo -e "${BLUE}=========================================${NC}"
echo -e "${BLUE}  Summary${NC}"
echo -e "${BLUE}=========================================${NC}"
echo "  Found:     $submission_count"
echo -e "  ${GREEN}Success:   $success_count${NC}"
echo -e "  ${RED}Failed:    $fail_count${NC}"
echo ""

if [ "$success_count" -gt 0 ]; then
    echo -e "${GREEN}✓ Extracted $success_count submissions to inbox${NC}"
    echo ""
    echo "Next steps:"
    echo "  1. Check status:     tools/workflow status $PROBLEM_ID"
    echo "  2. Process new only: tools/workflow process-new $PROBLEM_ID -j 4"
    echo "  3. Or run all:       tools/workflow run-all $PROBLEM_ID -j 4"
else
    echo -e "${YELLOW}No submissions extracted${NC}"
fi

exit 0
