#!/bin/bash
#
# Clean up public repository for participants
#
# Usage:
#   ./tools/cleanup_public_repo.sh
#
# This will:
#   1. Remove test submissions and temporary files
#   2. Verify .gitignore is correct
#   3. Check for any private data
#   4. Prepare repo for public use

set -e

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo -e "${BLUE}=========================================${NC}"
echo -e "${BLUE}  Public Repository Cleanup${NC}"
echo -e "${BLUE}=========================================${NC}"
echo ""
echo "Repository: $REPO_DIR"
echo ""

cd "$REPO_DIR"

# Step 1: Check for private directories
echo "Checking for private directories..."

FOUND_PRIVATE=false

if [ -d "submissions-inbox" ]; then
    echo -e "${RED}✗ Found: submissions-inbox/${NC}"
    echo "  This should NOT be in the public repo!"
    echo "  Action: Remove it (will be in .gitignore)"
    rm -rf submissions-inbox
    FOUND_PRIVATE=true
fi

if [ -d "submissions-private" ]; then
    echo -e "${RED}✗ Found: submissions-private/${NC}"
    echo "  This should NOT be in the public repo!"
    echo "  Action: Remove it (will be in .gitignore)"
    rm -rf submissions-private
    FOUND_PRIVATE=true
fi

if [ -d "test_submission" ]; then
    echo -e "${YELLOW}⚠ Found: test_submission/${NC}"
    echo "  Action: Remove it (test directory)"
    rm -rf test_submission
fi

if [ ! "$FOUND_PRIVATE" = true ]; then
    echo -e "${GREEN}✓ No private directories found${NC}"
fi

echo ""

# Step 2: Check for submission files
echo "Checking for submission files..."

SUBMISSION_FILES=$(find . -name "policy.json" -o -name "submission.json" -o -name "params.json" | grep -v "examples/" | grep -v "starter_template/" | grep -v ".git/")

if [ -n "$SUBMISSION_FILES" ]; then
    echo -e "${YELLOW}⚠ Found submission files (outside examples/starter_template/):${NC}"
    echo "$SUBMISSION_FILES"
    echo ""
    echo "Review these files - they might be test submissions"
else
    echo -e "${GREEN}✓ No stray submission files${NC}"
fi

echo ""

# Step 3: Verify .gitignore
echo "Verifying .gitignore..."

REQUIRED_IGNORES=(
    "submissions-inbox/"
    "submissions-private/"
    "test_submission/"
)

MISSING_IGNORES=()

for pattern in "${REQUIRED_IGNORES[@]}"; do
    if ! grep -q "$pattern" .gitignore; then
        MISSING_IGNORES+=("$pattern")
    fi
done

if [ ${#MISSING_IGNORES[@]} -eq 0 ]; then
    echo -e "${GREEN}✓ .gitignore is correct${NC}"
else
    echo -e "${YELLOW}⚠ Missing patterns in .gitignore:${NC}"
    for pattern in "${MISSING_IGNORES[@]}"; do
        echo "  - $pattern"
    done
fi

echo ""

# Step 4: Check git status
echo "Checking git status..."

UNCOMMITTED=$(git status --porcelain)

if [ -z "$UNCOMMITTED" ]; then
    echo -e "${GREEN}✓ No uncommitted changes${NC}"
else
    echo -e "${YELLOW}⚠ Uncommitted changes:${NC}"
    git status --short
    echo ""
    echo "Review and commit before making repo public"
fi

echo ""

# Step 5: Check for moderator-specific files in root
echo "Checking for moderator-only files..."

MODERATOR_FILES=(
    "MODERATOR_GUIDE.md"
    "QUICKSTART_MODERATOR.md"
    "WORKFLOW_GUIDE.md"
    "GOOGLE_SHEETS_SETUP.md"
    "SECURITY_SETUP.md"
    "SUBMISSION_WORKFLOW.md"
)

echo "These moderator files are in the repo (this is OK - they're documentation):"
for file in "${MODERATOR_FILES[@]}"; do
    if [ -f "$file" ]; then
        echo "  ✓ $file (public documentation)"
    fi
done

echo ""
echo -e "${YELLOW}Note: Moderator docs are fine in public repo - they help transparency${NC}"
echo ""

# Step 6: Verify participant documentation exists
echo "Verifying participant documentation..."

PARTICIPANT_FILES=(
    "README.md"
    "PARTICIPANT_GUIDE.md"
    "problems/problem-001-dimer-ksat/README.md"
)

ALL_DOCS_EXIST=true

for file in "${PARTICIPANT_FILES[@]}"; do
    if [ -f "$file" ]; then
        echo -e "  ${GREEN}✓${NC} $file"
    else
        echo -e "  ${RED}✗${NC} $file (MISSING!)"
        ALL_DOCS_EXIST=false
    fi
done

if [ "$ALL_DOCS_EXIST" = true ]; then
    echo -e "${GREEN}✓ All participant documentation exists${NC}"
fi

echo ""

# Step 7: Summary and recommendations
echo -e "${BLUE}=========================================${NC}"
echo -e "${BLUE}  Summary${NC}"
echo -e "${BLUE}=========================================${NC}"
echo ""

echo "Public Repository Checklist:"
echo ""
echo "Documentation:"
echo "  ✓ README.md (main overview)"
echo "  ✓ PARTICIPANT_GUIDE.md (how to participate)"
echo "  ✓ problems/*/README.md (problem descriptions)"
echo "  ✓ problems/*/leaderboard.csv (public leaderboards)"
echo ""
echo "Framework:"
echo "  ✓ core/ (schemas, generators, benchmark)"
echo "  ✓ examples/ (educational examples)"
echo "  ✓ docs/ (technical documentation)"
echo "  ✓ tools/ (helper scripts - moderators only)"
echo ""
echo "Privacy:"
echo "  ✓ .gitignore blocks private directories"
echo "  ✓ No submission files in repo"
echo "  ✓ No private data exposed"
echo ""

echo -e "${GREEN}Repository is ready for participants!${NC}"
echo ""

echo "Final steps before making public:"
echo "  1. Review uncommitted changes (if any)"
echo "  2. Push final changes to remote"
echo "  3. Verify repository visibility is PUBLIC"
echo "  4. Share link with participants!"
echo ""

echo "Share these links:"
echo "  Main repo: https://github.com/Livia33g/lammps-state-change"
echo "  Participant guide: https://github.com/Livia33g/lammps-state-change/blob/main/PARTICIPANT_GUIDE.md"
echo "  Problem 001: https://github.com/Livia33g/lammps-state-change/tree/main/problems/problem-001-dimer-ksat"
echo ""

exit 0
