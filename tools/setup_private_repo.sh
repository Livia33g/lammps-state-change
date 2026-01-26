#!/bin/bash
#
# Set up the private submissions repository
#
# Usage:
#   ./tools/setup_private_repo.sh
#
# This will:
#   1. Clone your private repo
#   2. Create directory structure
#   3. Add moderator documentation
#   4. Configure environment

set -e

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

PRIVATE_REPO_SSH="git@github.com:Livia33g/lammps-state-change-private.git"
PRIVATE_REPO_HTTPS="https://github.com/Livia33g/lammps-state-change-private.git"
PRIVATE_REPO_DIR="$HOME/lammps-state-change-private"

echo -e "${BLUE}=========================================${NC}"
echo -e "${BLUE}  Private Repository Setup${NC}"
echo -e "${BLUE}=========================================${NC}"
echo ""

# Step 1: Clone private repo
if [ -d "$PRIVATE_REPO_DIR" ]; then
    echo -e "${YELLOW}Private repo already exists at: $PRIVATE_REPO_DIR${NC}"
    echo "Updating instead of cloning..."
    cd "$PRIVATE_REPO_DIR"
    git pull
else
    echo "Cloning private repository..."
    echo ""
    echo "Choose cloning method:"
    echo "  1. SSH (recommended if you have SSH keys set up)"
    echo "  2. HTTPS (requires personal access token)"
    echo "  3. Manual (I'll set it up myself)"
    echo ""
    read -p "Enter choice [1-3]: " choice

    cd ~

    case $choice in
        1)
            echo "Using SSH..."
            git clone "$PRIVATE_REPO_SSH"
            if [ $? -ne 0 ]; then
                echo ""
                echo -e "${RED}SSH clone failed. You may need to set up SSH keys.${NC}"
                echo "See: https://docs.github.com/en/authentication/connecting-to-github-with-ssh"
                echo ""
                echo "Or run this script again and choose option 2 (HTTPS) or 3 (Manual)"
                exit 1
            fi
            ;;
        2)
            echo "Using HTTPS..."
            echo ""
            echo -e "${YELLOW}You'll need a GitHub Personal Access Token${NC}"
            echo "Create one at: https://github.com/settings/tokens"
            echo "Scopes needed: repo (full control of private repositories)"
            echo ""
            read -p "Press Enter when ready..."

            git clone "$PRIVATE_REPO_HTTPS"
            if [ $? -ne 0 ]; then
                echo ""
                echo -e "${RED}HTTPS clone failed.${NC}"
                echo "Make sure you entered your personal access token (not password)"
                echo ""
                echo "Or run this script again and choose option 3 (Manual)"
                exit 1
            fi
            ;;
        3)
            echo "Manual setup selected."
            echo ""
            echo "Please clone the repository manually:"
            echo ""
            echo "  cd ~"
            echo "  git clone git@github.com:Livia33g/lammps-state-change-private.git"
            echo ""
            echo "Then run this script again."
            exit 0
            ;;
        *)
            echo "Invalid choice. Exiting."
            exit 1
            ;;
    esac

    cd "$PRIVATE_REPO_DIR"
fi

echo -e "${GREEN}✓ Repository ready${NC}"
echo ""

# Step 2: Create directory structure
echo "Creating directory structure..."

mkdir -p submissions-inbox
mkdir -p submissions-archive
mkdir -p moderator-notes
mkdir -p .sheets-tracking

echo -e "${GREEN}✓ Directories created${NC}"
echo ""

# Step 3: Create README
echo "Creating README..."

cat > README.md << 'EOF'
# Competition Submissions (PRIVATE)

⚠️ **KEEP THIS REPOSITORY PRIVATE** ⚠️

This repository contains participant submissions and evaluation results.

## Directory Structure

```
lammps-state-change-private/
├── submissions-inbox/           # Place new submissions here
│   ├── alice_problem-001/
│   ├── bob_problem-001/
│   └── ...
├── submissions-archive/         # Evaluated submissions (auto-populated)
│   ├── problem-001-dimer-ksat/
│   │   ├── .workflow/
│   │   │   └── status.json     # Workflow tracking
│   │   ├── alice_2026-01-26/
│   │   │   ├── policy.json
│   │   │   └── generated/
│   │   │       └── analysis/
│   │   │           └── scores.json
│   │   └── bob_2026-01-27/
│   └── problem-002-octahedron/
├── moderator-notes/             # Private notes
└── .sheets-tracking/            # Google Sheets tracking
```

## Workflows

### Email Submissions

```bash
# 1. Save email attachments to ~/Downloads/username/
# 2. Extract to inbox
cd ~/lammps-state-change
tools/extract_email_submission.sh ~/Downloads/alice/ alice problem-001-dimer-ksat

# 3. Process
tools/workflow process-new problem-001-dimer-ksat -j 4
```

### Google Sheets Submissions

```bash
# 1. Export Google Sheet as CSV
# 2. Extract new submissions
cd ~/lammps-state-change
python3 tools/extract_from_sheets.py ~/Downloads/responses.csv problem-001-dimer-ksat

# 3. Process
tools/workflow process-new problem-001-dimer-ksat -j 4
```

### Batch Email Submissions

```bash
# 1. Save all submissions to ~/Downloads/submissions/alice/, bob/, etc.
# 2. Extract all
cd ~/lammps-state-change
tools/extract_all_submissions.sh ~/Downloads/submissions problem-001-dimer-ksat

# 3. Process
tools/workflow process-new problem-001-dimer-ksat -j 4
```

## Daily Routine

**Morning (5 minutes):**
```bash
# Check for new submissions (email/sheets)
# Extract to inbox (one of the methods above)
# Process new submissions
tools/workflow process-new problem-001-dimer-ksat -j 8
```

**Evening (2 minutes):**
```bash
# Check status
tools/workflow status problem-001-dimer-ksat

# Update leaderboard in public repo
cd ~/lammps-state-change
git add problems/*/leaderboard.csv
git commit -m "Update leaderboard"
git push
```

## Environment Setup

Make sure this is set in your `~/.bashrc`:

```bash
export SUBMISSIONS_PRIVATE_REPO="$HOME/lammps-state-change-private"
export LAMMPS_SRC="/path/to/lammps/src"
```

## Access Control

**Who should have access:**
- You (moderator)
- Other competition moderators
- Cluster (for evaluation)

**Who should NOT have access:**
- Participants
- Public

## Backups

Regularly backup this repository:

```bash
# Commit and push changes
cd ~/lammps-state-change-private
git add .
git commit -m "Archive submissions: $(date +%Y-%m-%d)"
git push

# Or create local backup
tar -czf submissions-backup-$(date +%Y%m%d).tar.gz ~/lammps-state-change-private/
```

## Documentation

See public repo for complete documentation:
- `MODERATOR_GUIDE.md` - Complete operations manual
- `QUICKSTART_MODERATOR.md` - Quick start guide
- `WORKFLOW_GUIDE.md` - Workflow system guide
- `GOOGLE_SHEETS_SETUP.md` - Google Sheets setup
- `SIMPLE_WORKFLOW.md` - Daily workflow examples

## Security

**NEVER:**
- Make this repository public
- Commit participant policies to public repo
- Share submission files with other participants
- Push to public repo by mistake

**Always:**
- Keep this repository private
- Add only moderators as collaborators
- Verify you're in the correct repo before pushing
- Use .gitignore in public repo to block private files
EOF

echo -e "${GREEN}✓ README created${NC}"
echo ""

# Step 4: Create .gitignore
echo "Creating .gitignore..."

cat > .gitignore << 'EOF'
# Python cache
__pycache__/
*.pyc
*.pyo

# OS files
.DS_Store
Thumbs.db

# Editor files
.vscode/
.idea/
*.swp
*.swo

# Temporary files
*.tmp
*.temp
.~*

# Don't ignore anything by default - we want to track all submissions
# But ignore common artifacts
*.log
*.err
*.out
slurm-*.out
EOF

echo -e "${GREEN}✓ .gitignore created${NC}"
echo ""

# Step 5: Create initial moderator notes
echo "Creating initial moderator notes..."

cat > moderator-notes/NOTES.md << 'EOF'
# Moderator Notes

## Competition Start Date
2026-01-26

## Active Problems
- problem-001-dimer-ksat

## Moderator Contact
[Your email]

## Important Decisions
[Document important decisions here]

## Participant Questions
[Track common questions and answers]

## Issues Encountered
[Document any issues and solutions]
EOF

echo -e "${GREEN}✓ Notes created${NC}"
echo ""

# Step 6: Commit everything
echo "Committing initial setup..."

git add .
git commit -m "Initial private repository setup

Directory structure:
- submissions-inbox/
- submissions-archive/
- moderator-notes/
- .sheets-tracking/

Documentation and .gitignore included." || echo "Nothing to commit (already set up)"

echo -e "${GREEN}✓ Changes committed${NC}"
echo ""

# Step 7: Push to remote
echo "Pushing to GitHub..."

git push || echo "Already up to date"

echo -e "${GREEN}✓ Pushed to remote${NC}"
echo ""

# Step 8: Configure environment
echo "Configuring environment..."

if ! grep -q "SUBMISSIONS_PRIVATE_REPO" ~/.bashrc; then
    echo "" >> ~/.bashrc
    echo "# Competition submission paths" >> ~/.bashrc
    echo "export SUBMISSIONS_PRIVATE_REPO=\"$PRIVATE_REPO_DIR\"" >> ~/.bashrc
    echo -e "${GREEN}✓ Added SUBMISSIONS_PRIVATE_REPO to ~/.bashrc${NC}"
else
    echo -e "${YELLOW}SUBMISSIONS_PRIVATE_REPO already in ~/.bashrc${NC}"
fi

echo ""
echo -e "${BLUE}=========================================${NC}"
echo -e "${BLUE}  Setup Complete!${NC}"
echo -e "${BLUE}=========================================${NC}"
echo ""
echo "Private repository: $PRIVATE_REPO_DIR"
echo ""
echo "Next steps:"
echo "  1. Reload environment: source ~/.bashrc"
echo "  2. Set LAMMPS_SRC if not already set:"
echo "     export LAMMPS_SRC=/path/to/lammps/src"
echo "     (Add to ~/.bashrc to make permanent)"
echo "  3. Test extraction:"
echo "     cd ~/lammps-state-change"
echo "     tools/workflow status problem-001-dimer-ksat"
echo ""
echo "Documentation: ~/lammps-state-change/MODERATOR_GUIDE.md"
echo ""

exit 0
