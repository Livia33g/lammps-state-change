# Submission Workflow Guide

Complete guide for how participants submit and how moderators process submissions.

---

## 🎯 Overview

**The Challenge:** Participants need an easy way to submit, but submissions must remain private.

**The Solution:** Multiple intake methods → Private repository → Automated processing

---

## 📋 Submission Methods (Choose One)

### Option 1: Email Submissions (Recommended - Easiest for Moderators)

**How it works:**
1. Participants email JSON files as attachments
2. Moderator saves attachments to local folder
3. Run extraction script → Private repo inbox
4. Run batch processor → Automated evaluation

**Pros:**
- ✅ Simple for participants (everyone knows email)
- ✅ Completely private (no public visibility)
- ✅ Semi-automated with extraction script
- ✅ Email provides audit trail

**Cons:**
- ⚠️ Manual email checking
- ⚠️ File extraction step required

**Setup:** See [Email Submission Setup](#email-submission-setup) below

---

### Option 2: Google Forms (Easiest for Participants)

**How it works:**
1. Create Google Form with text fields for JSON
2. Participants paste policy.json content
3. Responses go to Google Sheets
4. Export/script extracts to private repo

**Pros:**
- ✅ Very user-friendly for participants
- ✅ Structured data collection
- ✅ Timestamps automatic
- ✅ Can require login (prevents spam)

**Cons:**
- ⚠️ Manual extraction from sheets
- ⚠️ Large JSON may exceed field limits
- ⚠️ Requires Google account setup

**Setup:** See [Google Forms Setup](#google-forms-setup) below

---

### Option 3: GitHub Issues with Templates (Good for Tech-Savvy Participants)

**How it works:**
1. Create issue template in public repo
2. Participants create issue with `[SUBMISSION]` tag
3. Paste JSON in issue body
4. Moderator extracts to private repo
5. Close/lock issue immediately

**Pros:**
- ✅ Native to GitHub workflow
- ✅ Issue templates enforce structure
- ✅ Participants already on GitHub

**Cons:**
- ⚠️ **Temporarily public** (visible until moderator closes/locks)
- ⚠️ Other participants could see submissions briefly
- ⚠️ Manual extraction required

**Setup:** See [GitHub Issues Setup](#github-issues-setup) below

---

### Option 4: Dedicated Web Portal (Most Automated - Advanced)

**How it works:**
1. Build simple web form (Flask/Django/Node.js)
2. Participants upload JSON files
3. Server validates and saves to private repo
4. Triggers automated processing

**Pros:**
- ✅ Fully automated
- ✅ Real-time validation feedback
- ✅ Professional appearance
- ✅ Can integrate with authentication

**Cons:**
- ⚠️ Requires web development
- ⚠️ Hosting/maintenance overhead
- ⚠️ Security considerations

**Setup:** Not covered in this guide (requires custom development)

---

## 📧 Email Submission Setup (Recommended)

### Step 1: Set Up Private Repository

```bash
# On GitHub: Create new repository
# Name: lammps-state-change-private
# Visibility: PRIVATE ⚠️
# Add only moderators as collaborators

# Clone locally
cd ~
git clone git@github.com:YourOrg/lammps-state-change-private.git
cd lammps-state-change-private

# Create directory structure
mkdir -p submissions-inbox
mkdir -p submissions-archive
mkdir -p moderator-notes

# Create README
cat > README.md << 'EOF'
# Competition Submissions (PRIVATE)

⚠️ **KEEP THIS REPOSITORY PRIVATE** ⚠️

This repository contains participant submissions and evaluation results.

## Directories

- `submissions-inbox/` - Incoming submissions (place new ones here)
- `submissions-archive/` - Evaluated submissions (auto-populated)
- `moderator-notes/` - Private notes and decisions

## Access

Only moderators should have access.
NEVER make this repository public!
EOF

git add .
git commit -m "Initial private submission repository"
git push
```

### Step 2: Create Submission Email Template

Create an email template to send to participants:

```
Subject: How to Submit Your Solution - Problem {PROBLEM_ID}

Hi {NAME},

Ready to submit your solution? Here's how:

📧 SUBMISSION INSTRUCTIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Email your submission to: competition@yourlab.edu

Subject line: [SUBMISSION] {PROBLEM_ID} - {YOUR_USERNAME}

Attach the following files:

1. policy.json (required)
   Your state-change policy design

2. submission.json (required)
   Metadata about your submission

3. params.json (optional)
   Custom parameter values (if allowed for this problem)

📝 FILE EXAMPLES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

See examples in the starter template:
https://github.com/YourOrg/lammps-state-change/tree/main/problems/{PROBLEM_ID}/starter_template/

⏱️ EVALUATION TIMELINE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

- Submissions processed within 24 hours
- Evaluation takes up to 24 hours on cluster
- You'll receive results via email
- Leaderboard updates posted to GitHub

🔒 PRIVACY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Your submission will NOT be shared with other participants.
Only your score appears on the public leaderboard.

❓ QUESTIONS?

Reply to this email or create an issue:
https://github.com/YourOrg/lammps-state-change/issues

Good luck!
Competition Team
```

### Step 3: Install Email Extraction Script

```bash
cd ~/lammps-state-change
# Script already created - make executable
chmod +x tools/extract_email_submission.sh
```

### Step 4: Process Email Submissions

**When you receive a submission email:**

```bash
# 1. Save attachments to Downloads (or any folder)
# Example: ~/Downloads/alice_submission_20260126/

# 2. Extract to private repo inbox
cd ~/lammps-state-change
bash tools/extract_email_submission.sh \
    ~/Downloads/alice_submission_20260126/ \
    alice \
    problem-001-dimer-ksat

# Output:
#   =========================================
#     Extract Email Submission
#   =========================================
#     Attachment dir: ~/Downloads/alice_submission_20260126/
#     Username: alice
#     Problem: problem-001-dimer-ksat
#
#   ✓ Found all required files
#
#   Copying files to inbox...
#     ✓ policy.json
#     ✓ submission.json
#     ✓ params.json
#
#   ✓ Submission extracted to inbox
#     Location: ~/lammps-state-change-private/submissions-inbox/alice_problem-001-dimer-ksat

# 3. Process all submissions in inbox
bash tools/process_new_submissions.sh

# 4. Monitor cluster jobs
squeue -u $USER

# 5. Update leaderboard when complete
python3 tools/update_leaderboard.py \
    --submission ~/lammps-state-change-private/submissions-archive/problem-001-dimer-ksat/alice_2026-01-26_143022/

# 6. Email results to participant
```

### Step 5: Automate Email Checking (Optional)

Set up email forwarding or filtering:

**Gmail:**
1. Create filter: subject contains `[SUBMISSION]`
2. Auto-forward to dedicated account
3. Auto-label as "Competition Submissions"

**Then check regularly:**
```bash
# Check for new emails
# Extract attachments
# Run extraction script
```

Or use email API (advanced):
- Gmail API
- IMAP with Python
- Auto-extract attachments to inbox folder

---

## 📝 Google Forms Setup

### Step 1: Create Google Form

1. Go to https://forms.google.com
2. Create new form: "Problem Submission - {PROBLEM_ID}"

**Form fields:**

```
Title: Submit Your Solution - Problem {PROBLEM_ID}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Field 1: Username *
  Type: Short answer
  Description: Your username for the leaderboard
  Validation: Required

Field 2: Problem ID *
  Type: Multiple choice
  Options:
    - problem-001-dimer-ksat
    - problem-002-octahedron
  Validation: Required

Field 3: policy.json *
  Type: Paragraph
  Description: Paste the ENTIRE contents of your policy.json file
  Validation: Required

Field 4: params.json (optional)
  Type: Paragraph
  Description: If this problem allows custom parameters, paste params.json contents here. Otherwise leave blank.
  Validation: Optional

Field 5: Email *
  Type: Short answer
  Description: Your email for notification of results
  Validation: Email address, Required

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Submission confirmation message:
"Thank you! Your submission has been received. You'll receive evaluation results within 48 hours."
```

### Step 2: Link to Google Sheets

1. Responses tab → "View in Sheets"
2. Creates linked spreadsheet

### Step 3: Extract from Sheets

**Manual extraction:**
```bash
# 1. Open Google Sheet
# 2. For each new row:
#    - Copy username (column A)
#    - Copy problem_id (column B)
#    - Copy policy.json (column C)
#    - Copy params.json (column D) if present

# 3. Create files locally
mkdir -p ~/temp_submission/alice_problem-001-dimer-ksat/

# 4. Paste policy.json content
cat > ~/temp_submission/alice_problem-001-dimer-ksat/policy.json << 'EOF'
{paste content here}
EOF

# 5. Create submission.json manually
cat > ~/temp_submission/alice_problem-001-dimer-ksat/submission.json << EOF
{
  "problem_id": "problem-001-dimer-ksat",
  "username": "alice",
  "timestamp": "$(date -Iseconds)",
  "submission_method": "google_forms"
}
EOF

# 6. Extract to inbox
bash tools/extract_email_submission.sh \
    ~/temp_submission/alice_problem-001-dimer-ksat/ \
    alice \
    problem-001-dimer-ksat
```

**Semi-automated with script:**

Create `tools/extract_from_sheets.py`:
```python
#!/usr/bin/env python3
"""Extract submissions from Google Sheets export"""
import json
import sys
from pathlib import Path
from datetime import datetime

# Usage: python3 extract_from_sheets.py responses.csv

# Read CSV, extract JSON fields, create submission directories
# See tools/extract_from_sheets.py for full implementation
```

---

## 🐙 GitHub Issues Setup

### Step 1: Create Issue Template

Create `.github/ISSUE_TEMPLATE/submission.md`:

```markdown
---
name: Submit Solution
about: Submit your solution for evaluation
title: '[SUBMISSION] problem-ID - username'
labels: submission
assignees: ''
---

## Submission Information

**Username:** your_username_here
**Problem ID:** problem-001-dimer-ksat

## policy.json

```json
{
  "policy_name": "your_policy_name",
  "check_frequency": 100,
  "state_transitions": [
    {
      "from_species": "A",
      "to_species": "C",
      "trigger": {
        "contact_required": {
          "species": "B",
          "cutoff": 0.7,
          "min_contacts": 1
        }
      },
      "flip_probability": 1.0,
      "hysteresis_checks": 5
    }
  ]
}
```

## params.json (optional)

```json
{
  "morse_depth_AB": 1.2
}
```

## Checklist

- [ ] I have tested my JSON syntax is valid
- [ ] I have read the problem description
- [ ] I understand my score will be public but my policy will remain private

---

**Note:** This issue will be closed immediately after processing to maintain privacy.
```

### Step 2: Process GitHub Issue Submissions

**When a submission issue is created:**

```bash
# 1. Copy policy.json content from issue

# 2. Create local files
mkdir -p ~/temp_submission/alice_problem-001/

# 3. Save policy.json (copy from issue)
nano ~/temp_submission/alice_problem-001/policy.json

# 4. Create submission.json
cat > ~/temp_submission/alice_problem-001/submission.json << EOF
{
  "problem_id": "problem-001-dimer-ksat",
  "username": "alice",
  "timestamp": "$(date -Iseconds)",
  "submission_method": "github_issue",
  "issue_number": 42
}
EOF

# 5. Extract to inbox
bash tools/extract_email_submission.sh \
    ~/temp_submission/alice_problem-001/ \
    alice \
    problem-001-dimer-ksat

# 6. IMMEDIATELY close and lock the issue
gh issue close 42 --comment "Thank you! Submission received and processing. You'll be notified when evaluation is complete."
gh issue lock 42
```

**⚠️ Important:** Close and lock issues immediately to prevent other participants from seeing submissions!

---

## 🔄 Complete Workflow Comparison

| Step | Email | Google Forms | GitHub Issues |
|------|-------|--------------|---------------|
| **Participant submits** | Attach files to email | Fill out web form | Create issue with JSON |
| **Moderator receives** | Check email inbox | Check form responses | Get GitHub notification |
| **Extraction** | Save attachments, run script | Copy from sheets, create files | Copy from issue, create files |
| **Privacy** | ✅ Completely private | ✅ Private (form owner only) | ⚠️ Temporarily public |
| **Automation potential** | ⚙️ Medium (email API) | ⚙️ Medium (Sheets API) | ⚙️ Low (manual copy) |
| **Ease for participants** | ⭐⭐⭐⭐ Easy | ⭐⭐⭐⭐⭐ Very easy | ⭐⭐⭐ Moderate |
| **Ease for moderators** | ⭐⭐⭐ Good with script | ⭐⭐ Manual extraction | ⭐⭐ Manual + must close fast |

**Recommendation:** Start with **Email** for complete privacy and semi-automation.

---

## 🛠️ Setting Up Private Repo (Required for All Methods)

No matter which submission method you choose, you need the private repository:

```bash
# 1. Create on GitHub
# Repository name: lammps-state-change-private
# Visibility: PRIVATE ⚠️
# Initialize with README

# 2. Clone
git clone git@github.com:YourOrg/lammps-state-change-private.git
cd lammps-state-change-private

# 3. Create structure
mkdir -p submissions-inbox
mkdir -p submissions-archive
mkdir -p moderator-notes
echo "# Private Submissions - Moderators Only" > README.md

# 4. Configure environment
echo 'export SUBMISSIONS_PRIVATE_REPO="$HOME/lammps-state-change-private"' >> ~/.bashrc
source ~/.bashrc

# 5. Commit
git add .
git commit -m "Initial structure"
git push
```

**✅ Private repo is now ready!**

---

## 📊 Daily Workflow Summary

**Morning:**
```bash
# Check for new submissions (email/forms/issues)
# Extract any new submissions to inbox
cd ~/lammps-state-change
bash tools/process_new_submissions.sh
```

**Afternoon:**
```bash
# Check cluster job status
squeue -u $USER

# Update leaderboards for completed jobs
for submission in ~/lammps-state-change-private/submissions-archive/problem-*/*/; do
    if [ -f "$submission/generated/analysis/scores.json" ]; then
        python3 tools/update_leaderboard.py --submission "$submission"
    fi
done

# Commit leaderboard updates to PUBLIC repo
cd ~/lammps-state-change
git add problems/*/leaderboard.csv
git commit -m "Update leaderboards: $(date +%Y-%m-%d)"
git push
```

**Evening:**
```bash
# Email participants with results
# Archive processed submissions in private repo
cd ~/lammps-state-change-private
git add submissions-archive/
git commit -m "Archive submissions: $(date +%Y-%m-%d)"
git push
```

---

## ❓ FAQ

**Q: Do participants need GitHub accounts?**
A: Only for GitHub Issues method. Email and Google Forms don't require GitHub.

**Q: How long does processing take?**
A: Validation is instant. SLURM evaluation takes up to 24 hours depending on cluster load.

**Q: Can participants resubmit?**
A: Yes! New submissions replace old ones on leaderboard. Username remains the same.

**Q: What if someone submits malicious code?**
A: The validation script checks for suspicious patterns. JSON doesn't execute code - it generates C++ that you compile.

**Q: How do I handle multiple problems?**
A: Same workflow. Problem ID is part of submission metadata and directory structure.

---

## 📚 Related Documentation

- [SECURITY_SETUP.md](SECURITY_SETUP.md) - Private repository setup
- [QUICKSTART_MODERATOR.md](QUICKSTART_MODERATOR.md) - Quick start guide
- [tools/README.md](tools/README.md) - Tools reference

---

**Ready to accept submissions!** Choose your method and set up the private repo. 🚀
