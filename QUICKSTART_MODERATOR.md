# Moderator Quick Start Guide

Get your competition running in **5 minutes**! ⚡

---

## Prerequisites

✅ LAMMPS installed with RIGID and MOLECULE packages
✅ Python 3.7+
✅ SLURM cluster access (optional - can run locally)
✅ Git repository access

---

## 🚀 Setup (One-Time)

### 1. Clone Repository
```bash
git clone https://github.com/Livia33g/lammps-state-change.git
cd lammps-state-change
```

### 2. Create Inbox Directory
```bash
mkdir -p submissions-inbox
echo "# Place submissions here: username_problem-id/" > submissions-inbox/README.txt
```

### 3. Create Private Submissions Directory (Optional)
```bash
# Option A: Separate private repo (recommended)
cd ..
git clone <your-private-repo-url> lammps-state-change-submissions
cd lammps-state-change-submissions
mkdir -p problem-001-dimer-ksat/
cd ../lammps-state-change

# Option B: Local directory only
mkdir -p submissions-private
```

### 4. Verify Tools
```bash
# Check all tools are executable
ls -l tools/*.py tools/*.sh

# If not executable:
chmod +x tools/*.py tools/*.sh
```

---

## 📧 Processing Submissions

### Receive Submission via Email

**Participant sends:**
- `policy.json`
- `params.json` (optional)
- `submission.json`

**You extract to inbox:**
```bash
# Extract email attachments to downloads folder
cd ~/Downloads

# Create inbox directory for this submission
mkdir -p ~/lammps-state-change/submissions-inbox/alice_problem-001-dimer-ksat/

# Copy files
cp policy.json params.json submission.json \
   ~/lammps-state-change/submissions-inbox/alice_problem-001-dimer-ksat/
```

### Validate Submission
```bash
cd ~/lammps-state-change

python3 tools/validate_submission.py \
    submissions-inbox/alice_problem-001-dimer-ksat/
```

**If validation passes:**
```
✓ policy.json: Valid JSON
✓ submission.json: Valid JSON
✓ params.json: Valid JSON
✅ VALIDATION PASSED
```

**If validation fails:**
```
✗ Error: Missing required field 'flip_probability'
❌ VALIDATION FAILED
```
→ Email participant with error details and ask them to resubmit.

### Process All Submissions (Batch Mode)
```bash
# Process all submissions in inbox at once
bash tools/process_new_submissions.sh
```

**Output:**
```
=========================================
  Submission Batch Processor
=========================================

Found 1 submission(s) to process

----------------------------------------
Processing: alice_problem-001-dimer-ksat
----------------------------------------
  Username: alice
  Problem:  problem-001-dimer-ksat
  Validating...
  ✓ Validation passed
  Moving to: submissions-private/problem-001-dimer-ksat/alice_2026-01-26_143022
  Submitting to cluster...
  ✓ Submitted (Job ID: 12345)
  ✓ Processing complete
```

---

## 🖥️ Monitor Evaluation

### Check SLURM Queue
```bash
squeue -u $USER
```

### View Job Output (Real-Time)
```bash
# Find the job output file
ls submissions-private/problem-001-dimer-ksat/alice_*/eval_*.out

# Watch output
tail -f submissions-private/problem-001-dimer-ksat/alice_2026-01-26_143022/eval_12345.out
```

### Check for Errors
```bash
# Check LAMMPS stderr
cat submissions-private/problem-001-dimer-ksat/alice_*/generated/lammps_stderr.log

# Look for failures
grep -i "error\|fail" submissions-private/problem-001-dimer-ksat/alice_*/eval_*.out
```

---

## 📊 Update Leaderboard

### When Job Completes

**Check scores exist:**
```bash
ls submissions-private/problem-001-dimer-ksat/alice_*/generated/analysis/scores.json
```

**Update leaderboard (auto-detect mode):**
```bash
python3 tools/update_leaderboard.py \
    --submission submissions-private/problem-001-dimer-ksat/alice_2026-01-26_143022/
```

**Output:**
```
Adding new entry for alice
✓ Updated leaderboard: problems/problem-001-dimer-ksat/leaderboard.csv
  Total entries: 2
  Primary metric: work_per_yield (ascending)

Top 5:
  1. alice: work_per_yield=7.2000
  2. baseline: work_per_yield=8.3000
```

### Publish to Public Repo
```bash
git add problems/problem-001-dimer-ksat/leaderboard.csv
git commit -m "Update leaderboard: alice submission"
git push
```

---

## 📧 Notify Participant

### Results Ready Email Template

```
Subject: Evaluation Complete - problem-001-dimer-ksat

Hi [Name],

Your submission for problem-001-dimer-ksat has been evaluated!

📊 Results:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Final Yield:       0.78
Work per Yield:    7.2 (lower is better)
Flip Count:        412
Cumulative Work:   5.62
Time to Threshold: 980,000 steps
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

🏆 Leaderboard Position: #1

Congratulations! Your policy beat the baseline!

See the full leaderboard:
https://github.com/Livia33g/lammps-state-change/blob/main/problems/problem-001-dimer-ksat/leaderboard.csv

Want to improve further? You can submit updated policies anytime!

Best,
Competition Team
```

### Validation Failed Email Template

```
Subject: Submission Validation Error - problem-001-dimer-ksat

Hi [Name],

We received your submission for problem-001-dimer-ksat, but it failed validation.

❌ Error Details:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Missing required field 'flip_probability' in policy.json
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Please fix this issue and resubmit your files.

Need help? Check the documentation:
- Participant Guide: [link]
- Policy Schema: [link]
- Example Policies: [link]

Or reply to this email with questions!

Best,
Competition Team
```

---

## ⚡ Automated Workflow (Advanced)

### Set Up Cron Job

Process submissions automatically every hour:

```bash
# Edit crontab
crontab -e

# Add this line:
0 * * * * cd /path/to/lammps-state-change && bash tools/process_new_submissions.sh >> logs/processing.log 2>&1
```

### Create Log Directory
```bash
mkdir -p logs
```

### View Processing Logs
```bash
tail -f logs/processing.log
```

---

## 🚨 Troubleshooting

### Issue: Validation always fails

**Check:**
```bash
# Test with example submission
cp -r problems/problem-001-dimer-ksat/starter_template/* test_submission/
python3 tools/validate_submission.py test_submission/
```

**If still fails:** Check Python version (`python3 --version` should be 3.7+)

---

### Issue: SLURM job fails immediately

**Check SLURM output:**
```bash
cat eval_*.out
```

**Common causes:**
- LAMMPS not in PATH → Add to `.bashrc`
- Wrong partition name in `evaluate_submission.slurm` → Edit `#SBATCH --partition`
- Insufficient resources → Increase `#SBATCH --mem` or `--time`

---

### Issue: Leaderboard not updating

**Check scores exist:**
```bash
find submissions-private -name "scores.json"
```

**If missing:** Evaluation didn't complete or analysis failed

**Re-run analysis manually:**
```bash
python3 core/benchmark/score_policy_from_timeseries.py \
    --dump submissions-private/.../dump.*.lammpstrj \
    --thermo submissions-private/.../lammps_stdout.log \
    --output submissions-private/.../generated/analysis/
```

---

## 📚 Reference Documentation

Quick links to complete documentation:

| Document | Purpose |
|----------|---------|
| [MODERATOR_GUIDE.md](MODERATOR_GUIDE.md) | Complete operations manual |
| [tools/README.md](tools/README.md) | Detailed tools reference |
| [AUTOMATION_COMPLETE.md](AUTOMATION_COMPLETE.md) | Automation overview |
| [PARTICIPANT_GUIDE.md](PARTICIPANT_GUIDE.md) | For participants |

---

## ✅ Daily Checklist

**Morning:**
- [ ] Check inbox for new submissions
- [ ] Run `bash tools/process_new_submissions.sh`
- [ ] Check SLURM queue status

**Afternoon:**
- [ ] Check completed jobs
- [ ] Update leaderboards
- [ ] Commit and push leaderboard updates
- [ ] Email participants with results

**Evening:**
- [ ] Respond to participant questions
- [ ] Archive processed submissions
- [ ] Backup private submissions repo

---

## 🎓 First Competition Tips

1. **Start small:** Invite 5-10 participants first to test the pipeline
2. **Set deadlines:** Weekly or bi-weekly submission windows
3. **Be responsive:** Quick responses keep engagement high
4. **Celebrate wins:** Announce new leaderboard leaders
5. **Document issues:** Keep notes on common problems and solutions

---

## 🚀 You're Ready!

**Next steps:**
1. Set up inbox directory ✅
2. Announce competition to participants 📢
3. Start accepting submissions! 🎉

**Questions?** Check the full documentation or open an issue on GitHub.

**Good luck with your competition! 🏆**
