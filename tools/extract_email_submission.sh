#!/bin/bash
#
# Extract submission files from email attachments
#
# Usage:
#   1. Save email attachments to a folder (e.g., ~/Downloads/alice_submission/)
#   2. Run: ./tools/extract_email_submission.sh ~/Downloads/alice_submission/ alice problem-001-dimer-ksat
#
# This will:
#   - Validate files are present
#   - Copy to private repo inbox
#   - Ready for batch processing

set -e

ATTACHMENT_DIR="$1"
USERNAME="$2"
PROBLEM_ID="$3"

# Configuration
PRIVATE_REPO="${SUBMISSIONS_PRIVATE_REPO:-$HOME/lammps-state-change-private}"
INBOX_DIR="$PRIVATE_REPO/submissions-inbox"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Validate arguments
if [ -z "$ATTACHMENT_DIR" ] || [ -z "$USERNAME" ] || [ -z "$PROBLEM_ID" ]; then
    echo "Usage: $0 <attachment_dir> <username> <problem_id>"
    echo ""
    echo "Example:"
    echo "  $0 ~/Downloads/alice_submission alice problem-001-dimer-ksat"
    exit 1
fi

if [ ! -d "$ATTACHMENT_DIR" ]; then
    echo -e "${RED}Error: Directory not found: $ATTACHMENT_DIR${NC}"
    exit 1
fi

echo "========================================="
echo "  Extract Email Submission"
echo "========================================="
echo "  Attachment dir: $ATTACHMENT_DIR"
echo "  Username: $USERNAME"
echo "  Problem: $PROBLEM_ID"
echo ""

# Check for required files
required_files=("policy.json" "submission.json")
missing_files=()

for file in "${required_files[@]}"; do
    if [ ! -f "$ATTACHMENT_DIR/$file" ]; then
        missing_files+=("$file")
    fi
done

if [ ${#missing_files[@]} -ne 0 ]; then
    echo -e "${RED}✗ Missing required files:${NC}"
    for file in "${missing_files[@]}"; do
        echo "  - $file"
    done
    echo ""
    echo "Required files: policy.json, submission.json"
    echo "Optional files: params.json"
    exit 1
fi

echo -e "${GREEN}✓ Found all required files${NC}"

# Create inbox directory
DEST_DIR="$INBOX_DIR/${USERNAME}_${PROBLEM_ID}"
mkdir -p "$DEST_DIR"

# Copy files
echo ""
echo "Copying files to inbox..."
cp "$ATTACHMENT_DIR/policy.json" "$DEST_DIR/"
cp "$ATTACHMENT_DIR/submission.json" "$DEST_DIR/"

if [ -f "$ATTACHMENT_DIR/params.json" ]; then
    cp "$ATTACHMENT_DIR/params.json" "$DEST_DIR/"
    echo -e "  ${GREEN}✓ policy.json${NC}"
    echo -e "  ${GREEN}✓ submission.json${NC}"
    echo -e "  ${GREEN}✓ params.json${NC}"
else
    echo -e "  ${GREEN}✓ policy.json${NC}"
    echo -e "  ${GREEN}✓ submission.json${NC}"
    echo -e "  ${YELLOW}○ params.json (optional - not provided)${NC}"
fi

echo ""
echo -e "${GREEN}✓ Submission extracted to inbox${NC}"
echo "  Location: $DEST_DIR"
echo ""
echo "Next steps:"
echo "  1. Review submission files (optional)"
echo "  2. Run: bash tools/process_new_submissions.sh"
echo "  3. Monitor cluster jobs: squeue -u \$USER"

exit 0
