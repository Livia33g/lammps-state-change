# Moderator Privacy & Security Setup

## 🚨 Critical Privacy Issue Resolved

**Problem:** If `submissions-inbox/` or `submissions-private/` are in the public repo, participant submissions would be visible to everyone!

**Solution:** These directories are now in `.gitignore` and will NEVER be committed to the public repo.

---

## 🔒 Recommended Repository Structure

### Two-Repo Setup (Recommended)

**PUBLIC REPO** (`lammps-state-change`) - What players see:
```
lammps-state-change/                    # Public GitHub repo
├── problems/                           # Problem definitions
├── core/                              # Framework code
├── docs/                              # Documentation
├── examples/                          # Educational examples
├── tools/                             # Tools (scripts only)
├── README.md                          # Public documentation
├── PARTICIPANT_GUIDE.md
└── problems/*/leaderboard.csv         # Public leaderboards
```

**PRIVATE REPO** (`lammps-state-change-private`) - Moderators only:
```
lammps-state-change-private/           # Private GitHub repo
├── submissions-inbox/                 # Incoming submissions
│   ├── alice_problem-001/
│   ├── bob_problem-001/
│   └── ...
├── submissions-archive/               # Evaluated submissions
│   ├── problem-001/
│   │   ├── alice_2026-01-26_143022/
│   │   ├── bob_2026-01-27_091534/
│   │   └── ...
│   └── problem-002/
├── moderator-notes/                   # Private notes
├── evaluation-logs/                   # Processing logs
└── README.md                          # Moderator instructions
```

---

## 🛠️ Setup Instructions

### Step 1: Create Private Repository

```bash
# On GitHub: Create new PRIVATE repository
# Name: lammps-state-change-private
# Visibility: Private
# Add moderators as collaborators

# Clone private repo
cd ~/
git clone git@github.com:YourOrg/lammps-state-change-private.git
cd lammps-state-change-private

# Create directory structure
mkdir -p submissions-inbox
mkdir -p submissions-archive
mkdir -p moderator-notes
mkdir -p evaluation-logs

# Initialize
cat > README.md << 'EOF'
# Competition Submissions (PRIVATE)

⚠️ **KEEP THIS REPOSITORY PRIVATE** ⚠️

This repository contains participant submissions and evaluation results.

## Directory Structure

- `submissions-inbox/` - Place new submissions here for processing
- `submissions-archive/` - Evaluated submissions (organized by problem)
- `moderator-notes/` - Private moderator notes and decisions
- `evaluation-logs/` - Processing and evaluation logs

## Workflow

1. Receive submission via email
2. Extract to `submissions-inbox/username_problem-id/`
3. Process: `bash ~/lammps-state-change/tools/process_new_submissions.sh`
4. Monitor cluster jobs
5. Update leaderboard in PUBLIC repo
6. Archive moves submission to `submissions-archive/`

## Access Control

Only moderators should have access to this repository.
NEVER commit this directory to the public repo!
EOF

git add .
git commit -m "Initial moderator private repository structure"
git push
```

### Step 2: Configure Environment

```bash
# Add to your ~/.bashrc or ~/.bash_profile
export SUBMISSIONS_PRIVATE_REPO="$HOME/lammps-state-change-private"
export LAMMPS_SRC="/path/to/lammps/src"

# Reload
source ~/.bashrc
```

The `process_new_submissions.sh` script automatically uses this environment variable.

### Step 3: Verify .gitignore Protection

Already done! These directories are protected:
```gitignore
# CRITICAL: Moderator-only directories - NEVER commit these!
submissions-inbox/
submissions-private/
test_submission/
```

---

## 🔐 Security Checklist

### ✅ What's Protected

- [x] `submissions-inbox/` - Never committed (in .gitignore)
- [x] `submissions-private/` - Never committed (in .gitignore)
- [x] `test_submission/` - Never committed (in .gitignore)
- [x] Private repo separate from public
- [x] Only moderators have access to private repo
- [x] Participant policies never visible to other participants

### ❌ What GitHub Cannot Do

- [ ] **Branch privacy** - GitHub doesn't support private branches in public repos
  - Solution: Use separate private repository instead

---

## 🚨 Important Reminders

### For Moderators

1. **NEVER** run `git add submissions-inbox/` or `git add submissions-private/` in public repo
2. **ALWAYS** check `git status` before committing to public repo
3. **KEEP** private repo truly private (GitHub settings → Visibility → Private)
4. **ADD** collaborators carefully to private repo
5. **USE** email or private repo for submission intake

### Workflow Summary

**Public Repo (players see):**
- Problem definitions
- Framework code
- Documentation
- Leaderboards (scores only, no policies)

**Private Repo (moderators only):**
- Actual submission files
- Evaluation results
- Generated C++ fixes
- Trajectory files
- Moderator notes

---

## 📂 Complete Directory Layout

### Your Local Machine

```
~/lammps-state-change/                          # PUBLIC repo clone
├── problems/
├── core/
├── tools/
└── README.md

~/lammps-state-change-private/                  # PRIVATE repo clone
├── submissions-inbox/                          # Place new submissions here
│   ├── alice_problem-001-dimer-ksat/
│   │   ├── policy.json
│   │   ├── params.json
│   │   └── submission.json
│   └── bob_problem-001-dimer-ksat/
├── submissions-archive/                        # Auto-populated after evaluation
│   ├── problem-001-dimer-ksat/
│   │   ├── alice_2026-01-26_143022/
│   │   │   ├── policy.json
│   │   │   ├── submission.json
│   │   │   ├── generated/
│   │   │   │   ├── fix_state_change_*.cpp
│   │   │   │   ├── data.*.lammps
│   │   │   │   ├── dump.*.lammpstrj
│   │   │   │   └── analysis/scores.json
│   │   │   └── eval_12345.out
│   │   └── bob_2026-01-27_091534/
│   └── problem-002-octahedron/
├── moderator-notes/
│   ├── problem-001-notes.md
│   └── decisions.md
└── evaluation-logs/
    └── 2026-01-26-processing.log
```

---

## ⚙️ Updated Workflow

### Receive & Process Submission

```bash
# 1. Receive email with attachments
cd ~/Downloads

# 2. Move to PRIVATE repo inbox
mkdir -p ~/lammps-state-change-private/submissions-inbox/alice_problem-001-dimer-ksat/
cp policy.json params.json submission.json \
   ~/lammps-state-change-private/submissions-inbox/alice_problem-001-dimer-ksat/

# 3. Process from public repo tools
cd ~/lammps-state-change
bash tools/process_new_submissions.sh

# The script will:
#   - Find submissions in ~/lammps-state-change-private/submissions-inbox/
#   - Validate them
#   - Move to ~/lammps-state-change-private/submissions-archive/
#   - Submit to SLURM cluster
```

### Update Leaderboard (Public)

```bash
# 1. Update leaderboard from archived submission
python3 tools/update_leaderboard.py \
    --submission ~/lammps-state-change-private/submissions-archive/problem-001/alice_2026-01-26/

# 2. Commit ONLY leaderboard to public repo
cd ~/lammps-state-change
git add problems/problem-001-dimer-ksat/leaderboard.csv
git commit -m "Update leaderboard: alice submission"
git push

# 3. Commit submission to private repo
cd ~/lammps-state-change-private
git add submissions-archive/problem-001/alice_2026-01-26/
git commit -m "Archive: alice problem-001 submission"
git push
```

---

## 🎯 Why This Works

| Concern | Solution |
|---------|----------|
| Submissions visible to players | Private repo, never in public repo |
| Branch privacy | Separate repos (not branches) |
| Accidental commits | .gitignore blocks sensitive dirs |
| Moderator collaboration | Private repo with collaborator access |
| Public leaderboard | Only scores committed to public repo |
| Cluster file exposure | Use private directories, clean up after eval |

---

## 🔧 Alternative: Local-Only Directories

If you don't want a separate private repo, you can use local directories only:

```bash
# Create local-only directories (already in .gitignore)
mkdir -p ~/submissions-private-local/inbox
mkdir -p ~/submissions-private-local/archive

# Configure script to use local path
export SUBMISSIONS_PRIVATE_REPO="$HOME/submissions-private-local"

# Backup regularly (not via git)
tar -czf submissions-backup-$(date +%Y%m%d).tar.gz ~/submissions-private-local/
```

**Trade-offs:**
- ✅ Simpler setup
- ✅ No second repo to manage
- ❌ No version control for submissions
- ❌ Harder to collaborate with other moderators
- ❌ Manual backups required

---

## 📚 Related Documentation

- [MODERATOR_GUIDE.md](MODERATOR_GUIDE.md) - Complete operations manual
- [QUICKSTART_MODERATOR.md](QUICKSTART_MODERATOR.md) - 5-minute setup
- [tools/README.md](tools/README.md) - Tools reference

**All recommend the two-repo setup for maximum security!**
