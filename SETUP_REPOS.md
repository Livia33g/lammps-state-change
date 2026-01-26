# Repository Setup Guide

**Setting up the two-repository system for your competition**

---

## 🎯 Overview

You'll manage **two repositories**:

### 1. **Public Repo** (lammps-state-change)
- Problems, framework, documentation
- Leaderboards (scores only)
- Visible to everyone
- **Already exists:** https://github.com/Livia33g/lammps-state-change

### 2. **Private Repo** (lammps-state-change-private)
- Actual submissions
- Evaluation results
- Moderator notes
- Visible only to moderators
- **Already created:** https://github.com/Livia33g/lammps-state-change-private

---

## 🚀 Setup (10 Minutes)

### Step 1: Set Up Private Repository

```bash
cd ~/lammps-state-change
tools/setup_private_repo.sh
```

**What this does:**
- ✅ Clones https://github.com/Livia33g/lammps-state-change-private
- ✅ Creates directory structure (submissions-inbox/, submissions-archive/, etc.)
- ✅ Adds README and documentation
- ✅ Configures environment variable `SUBMISSIONS_PRIVATE_REPO`
- ✅ Commits and pushes initial setup

**Output:**
```
=========================================
  Private Repository Setup
=========================================

Cloning private repository...
✓ Repository ready

Creating directory structure...
✓ Directories created

Creating README...
✓ README created

Committing initial setup...
✓ Changes committed

Pushing to GitHub...
✓ Pushed to remote

Configuring environment...
✓ Added SUBMISSIONS_PRIVATE_REPO to ~/.bashrc

=========================================
  Setup Complete!
=========================================

Private repository: /home/user/lammps-state-change-private

Next steps:
  1. Reload environment: source ~/.bashrc
  2. Set LAMMPS_SRC if not already set
  3. Test extraction
```

### Step 2: Reload Environment

```bash
source ~/.bashrc

# Verify it's set
echo $SUBMISSIONS_PRIVATE_REPO
# Output: /home/user/lammps-state-change-private
```

### Step 3: Set LAMMPS Path

```bash
# Add to ~/.bashrc
echo 'export LAMMPS_SRC="/path/to/your/lammps/src"' >> ~/.bashrc
source ~/.bashrc

# Verify
echo $LAMMPS_SRC
```

### Step 4: Clean Up Public Repository

```bash
cd ~/lammps-state-change
tools/cleanup_public_repo.sh
```

**What this does:**
- ✅ Removes any test submissions
- ✅ Verifies .gitignore is correct
- ✅ Checks for private data
- ✅ Ensures repo is ready for participants

**Output:**
```
=========================================
  Public Repository Cleanup
=========================================

Checking for private directories...
✓ No private directories found

Checking for submission files...
✓ No stray submission files

Verifying .gitignore...
✓ .gitignore is correct

Checking git status...
✓ No uncommitted changes

=========================================
  Summary
=========================================

Public Repository Checklist:
  ✓ Documentation exists
  ✓ Framework code ready
  ✓ .gitignore blocks private directories
  ✓ No private data exposed

✓ Repository is ready for participants!
```

### Step 5: Verify Setup

```bash
# Check that both repos exist
ls -ld ~/lammps-state-change
ls -ld ~/lammps-state-change-private

# Test workflow status
cd ~/lammps-state-change
tools/workflow status problem-001-dimer-ksat

# Output:
#   Total submissions: 0
#   (no submissions yet - this is correct!)
```

---

## 📂 Directory Structure (After Setup)

### Your Machine

```
~/
├── lammps-state-change/              # PUBLIC repo (clone of GitHub)
│   ├── problems/
│   │   └── problem-001-dimer-ksat/
│   │       ├── README.md
│   │       ├── leaderboard.csv       # Public leaderboard
│   │       └── starter_template/
│   ├── core/                         # Framework code
│   ├── docs/                         # Documentation
│   ├── tools/                        # Helper scripts
│   ├── README.md                     # Main overview
│   └── PARTICIPANT_GUIDE.md
│
└── lammps-state-change-private/      # PRIVATE repo (clone of GitHub)
    ├── submissions-inbox/            # NEW submissions go here
    ├── submissions-archive/          # EVALUATED submissions (auto-populated)
    │   └── problem-001-dimer-ksat/
    │       ├── .workflow/
    │       │   └── status.json       # Workflow tracking
    │       ├── alice_2026-01-26/
    │       └── bob_2026-01-27/
    ├── moderator-notes/              # Your private notes
    ├── .sheets-tracking/             # Google Sheets tracking (auto-created)
    └── README.md                     # Moderator instructions
```

### GitHub

**Public:** https://github.com/Livia33g/lammps-state-change
- Visible to everyone
- Players can clone and explore

**Private:** https://github.com/Livia33g/lammps-state-change-private
- Visible only to you and moderators
- Contains all submissions and results

---

## ✅ Quick Test

### Test Submission Extraction

```bash
# 1. Create a test submission
mkdir -p ~/test_submission
cat > ~/test_submission/policy.json << 'EOF'
{
  "policy_name": "test_policy",
  "check_frequency": 100,
  "state_transitions": []
}
EOF

cat > ~/test_submission/submission.json << 'EOF'
{
  "problem_id": "problem-001-dimer-ksat",
  "username": "test_user"
}
EOF

# 2. Extract to inbox
cd ~/lammps-state-change
tools/extract_email_submission.sh ~/test_submission test_user problem-001-dimer-ksat

# Output:
#   ✓ Submission extracted to inbox
#   Location: /home/user/lammps-state-change-private/submissions-inbox/test_user_problem-001

# 3. Verify it's in the private repo
ls ~/lammps-state-change-private/submissions-inbox/

# Output: test_user_problem-001-dimer-ksat

# 4. Check workflow status
tools/workflow status problem-001-dimer-ksat

# Output:
#   Total submissions: 1
#   validated: 0 complete, 1 pending

# 5. Clean up test
rm -rf ~/lammps-state-change-private/submissions-inbox/test_user_*
```

---

## 🎮 Daily Workflow (After Setup)

### Receive Submissions

**Email method:**
```bash
cd ~/lammps-state-change
tools/extract_email_submission.sh ~/Downloads/alice/ alice problem-001-dimer-ksat
```

**Google Sheets method:**
```bash
cd ~/lammps-state-change
python3 tools/extract_from_sheets.py ~/Downloads/responses.csv problem-001-dimer-ksat
```

**Batch method:**
```bash
cd ~/lammps-state-change
tools/extract_all_submissions.sh ~/Downloads/submissions problem-001-dimer-ksat
```

### Process Submissions

```bash
# Process only new submissions
cd ~/lammps-state-change
tools/workflow process-new problem-001-dimer-ksat -j 4

# Check status
tools/workflow status problem-001-dimer-ksat

# Clean up disk space
tools/workflow cleanup problem-001-dimer-ksat
```

### Update Leaderboard

```bash
# Leaderboard is auto-updated by workflow
# Just commit and push to public repo

cd ~/lammps-state-change
git add problems/*/leaderboard.csv
git commit -m "Update leaderboard: $(date +%Y-%m-%d)"
git push
```

---

## 🔐 Security Checklist

### Public Repo (lammps-state-change)

✅ **Should contain:**
- Problem definitions
- Framework code
- Documentation
- Leaderboards (scores only)
- Example policies

❌ **Should NOT contain:**
- Participant submissions
- Actual policy.json files from participants
- Private moderator notes
- Evaluation results (except scores)

### Private Repo (lammps-state-change-private)

✅ **Should contain:**
- All submissions
- Evaluation results
- Moderator notes
- Workflow tracking data

❌ **Should NOT be:**
- Made public
- Shared with participants
- Cloned to public machines

---

## 🛠️ Troubleshooting

### Error: "SUBMISSIONS_PRIVATE_REPO not set"

```bash
# Check if it's set
echo $SUBMISSIONS_PRIVATE_REPO

# If empty, add to ~/.bashrc
echo 'export SUBMISSIONS_PRIVATE_REPO="$HOME/lammps-state-change-private"' >> ~/.bashrc
source ~/.bashrc
```

### Error: Private repo not found

```bash
# Check if it exists
ls ~/lammps-state-change-private

# If not, run setup again
cd ~/lammps-state-change
tools/setup_private_repo.sh
```

### Error: Cannot push to private repo

```bash
# Check remote URL
cd ~/lammps-state-change-private
git remote -v

# Should show:
#   origin  https://github.com/Livia33g/lammps-state-change-private.git

# If wrong, update it
git remote set-url origin https://github.com/Livia33g/lammps-state-change-private.git
```

### Accidentally committed to wrong repo

```bash
# If you committed submission to PUBLIC repo:
# 1. Don't push!
# 2. Reset last commit
git reset HEAD~1

# 3. Verify you're in correct repo
pwd
# Should be ~/lammps-state-change-private for submissions
```

---

## 📚 Documentation

After setup, reference these guides:

**For daily work:**
- `SIMPLE_WORKFLOW.md` - Daily workflow (2 commands)
- `WORKFLOW_GUIDE.md` - Complete workflow system guide
- `QUICKSTART_MODERATOR.md` - 5-minute moderator guide

**For submission methods:**
- `GOOGLE_SHEETS_SETUP.md` - Google Sheets setup (recommended)
- `SUBMISSION_WORKFLOW.md` - All submission methods

**For operations:**
- `MODERATOR_GUIDE.md` - Complete operations manual
- `SECURITY_SETUP.md` - Two-repo security model
- `~/lammps-state-change-private/README.md` - Private repo guide

---

## 🎉 You're Ready!

**Setup complete when you can:**

```bash
# 1. Check both repos exist
ls ~/lammps-state-change
ls ~/lammps-state-change-private

# 2. Environment variables set
echo $SUBMISSIONS_PRIVATE_REPO
echo $LAMMPS_SRC

# 3. Workflow status works
cd ~/lammps-state-change
tools/workflow status problem-001-dimer-ksat

# 4. Test extraction works
tools/extract_email_submission.sh ~/test_submission test problem-001-dimer-ksat
```

**Next steps:**
1. ✅ Set up Google Form (optional) - see `GOOGLE_SHEETS_SETUP.md`
2. ✅ Announce competition to participants
3. ✅ Share public repo link
4. ✅ Start accepting submissions!

---

## 🔄 Backup Strategy

### Daily (Automatic)

```bash
# Commits happen automatically in private repo
cd ~/lammps-state-change-private
git add .
git commit -m "Daily backup: $(date +%Y-%m-%d)"
git push
```

### Weekly (Manual Backup)

```bash
# Create local archive
tar -czf ~/backups/submissions-$(date +%Y%m%d).tar.gz ~/lammps-state-change-private/
```

---

**Questions?** See `MODERATOR_GUIDE.md` or create an issue in the public repo!
