# Moderator Tools

This directory contains automation scripts for processing submissions and maintaining the competition.

---

## 📋 Tool Overview

| Tool | Purpose | Usage |
|------|---------|-------|
| `validate_submission.py` | Validate submission format and security | Run before accepting submissions |
| `process_new_submissions.sh` | Batch process submissions from inbox | Run periodically to handle new submissions |
| `evaluate_submission.slurm` | Evaluate single submission on cluster | Called by process script or manually |
| `update_leaderboard.py` | Update leaderboard from scores | Run after evaluation completes |

---

## 🔧 Individual Tools

### 1. validate_submission.py

**Purpose:** Validate submission structure, JSON syntax, and security before processing.

**Usage:**
```bash
python3 tools/validate_submission.py submissions-inbox/alice_problem-001/
```

**Checks:**
- Required files present (`policy.json`, `submission.json`)
- Valid JSON syntax
- Schema compliance
- Parameter ranges
- Suspicious content detection
- File size limits

**Exit codes:**
- `0`: Validation passed
- `1`: Validation failed

**Example output:**
```
✓ Found policy.json
✓ Found submission.json
✓ Valid JSON syntax
✓ Policy schema valid
✓ Parameters within allowed ranges
✓ No security issues detected

Submission is valid!
```

---

### 2. process_new_submissions.sh

**Purpose:** Automatically process all submissions in the inbox directory.

**Usage:**
```bash
# Process all submissions in submissions-inbox/
bash tools/process_new_submissions.sh
```

**Directory structure:**
```
submissions-inbox/
├── alice_problem-001/
│   ├── policy.json
│   ├── params.json
│   └── submission.json
└── bob_problem-001/
    └── ...
```

**Workflow:**
1. Scans `submissions-inbox/` for new submissions
2. Validates each submission
3. Moves valid submissions to `submissions-private/`
4. Submits to cluster via SLURM
5. Archives processed submissions to `.processed/`

**Output:**
```
=========================================
  Submission Batch Processor
=========================================

Found 2 submission(s) to process

----------------------------------------
Processing: alice_problem-001
----------------------------------------
  Username: alice
  Problem:  problem-001
  Validating...
  ✓ Validation passed
  Moving to: submissions-private/problem-001/alice_2026-01-26_143022
  Submitting to cluster...
  ✓ Submitted (Job ID: 12345)
  ✓ Processing complete

=========================================
  Summary
=========================================
  Processed: 2
  Failed:    0
```

**Failed submissions:**
- Remain in inbox
- Error details printed to console
- Moderator should email participant with error message

---

### 3. evaluate_submission.slurm

**Purpose:** Complete evaluation pipeline for a single submission (runs on cluster).

**Usage:**
```bash
sbatch tools/evaluate_submission.slurm problem-001-dimer-ksat alice_2026-01-26
```

**Pipeline steps:**
1. Validate submission
2. Generate C++ fix from policy
3. Generate LAMMPS data/input files
4. Compile LAMMPS with custom fix
5. Run simulation (up to 24 hours)
6. Analyze trajectory and calculate scores
7. Generate scores.json

**SLURM configuration:**
```bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=4
#SBATCH --mem=8G
```

**Output files:**
```
submissions-private/problem-001/alice_2026-01-26/
├── policy.json                      # Original submission
├── params.json
├── submission.json
├── generated/
│   ├── fix_state_change_*.cpp       # Auto-generated C++ fix
│   ├── fix_state_change_*.h
│   ├── data.001_dimer_ksat          # LAMMPS data file
│   ├── in.001_dimer_ksat            # LAMMPS input script
│   ├── dump.001_dimer_ksat.lammpstrj # Trajectory
│   ├── lammps_stdout.log            # LAMMPS output
│   ├── lammps_stderr.log            # LAMMPS errors + events
│   └── analysis/
│       ├── scores.json              # Final scores
│       └── analysis.png             # Visualization
└── eval_JOBID.out                   # SLURM output
```

**Monitoring:**
```bash
# Check job status
squeue -u $USER

# View SLURM output (while running)
tail -f eval_JOBID.out

# View scores (when complete)
cat submissions-private/problem-001/alice/generated/analysis/scores.json
```

---

### 4. update_leaderboard.py

**Purpose:** Update leaderboard CSV from evaluation results.

**Usage:**

**Option A: Manual mode**
```bash
python3 tools/update_leaderboard.py \
    --problem problem-001-dimer-ksat \
    --username alice \
    --scores submissions-private/problem-001/alice_2026-01-26/generated/analysis/scores.json
```

**Option B: Auto-detect mode**
```bash
python3 tools/update_leaderboard.py \
    --submission submissions-private/problem-001/alice_2026-01-26/
```

**Option C: Anonymize leaderboard**
```bash
python3 tools/update_leaderboard.py \
    --submission submissions-private/problem-001/alice_2026-01-26/ \
    --anonymize
```

**Features:**
- Automatically calculates `work_per_yield` if not present
- Updates existing entries (if username exists)
- Sorts by primary metric (from `problem.json`)
- Preserves all existing entries
- Optional anonymization (removes usernames, shows ranks)

**Output:**
```
Updating existing entry for alice
✓ Updated leaderboard: problems/problem-001-dimer-ksat/leaderboard.csv
  Total entries: 3
  Primary metric: work_per_yield (ascending)

Top 5:
  1. alice: work_per_yield=7.2000
  2. baseline: work_per_yield=8.3000
  3. bob: work_per_yield=9.1000
```

**Leaderboard format:**
```csv
username,final_yield,work_per_yield,flip_count,cumulative_work,date
alice,0.7800,7.2000,412,5.6160,2026-01-26
baseline,0.7200,8.3000,450,5.9760,2026-01-20
bob,0.6900,9.1000,520,6.2790,2026-01-25
```

---

## 🔄 Complete Workflow

### Email-Based Submissions

**Step 1: Receive submission via email**
```bash
# Extract attachments to inbox
mkdir -p submissions-inbox/alice_problem-001/
cp ~/Downloads/{policy,params,submission}.json submissions-inbox/alice_problem-001/
```

**Step 2: Process submissions**
```bash
# Batch process all inbox submissions
bash tools/process_new_submissions.sh
```

**Step 3: Monitor cluster jobs**
```bash
# Check status
squeue -u $USER

# View output
tail -f submissions-private/problem-001/alice_*/eval_*.out
```

**Step 4: Update leaderboard when complete**
```bash
# Auto-detect mode
python3 tools/update_leaderboard.py \
    --submission submissions-private/problem-001/alice_2026-01-26_143022/
```

**Step 5: Commit leaderboard to public repo**
```bash
git add problems/problem-001-dimer-ksat/leaderboard.csv
git commit -m "Update leaderboard: alice submission"
git push
```

**Step 6: Email participant with results**
```
Subject: Evaluation Complete - problem-001

Hi Alice,

Your submission has been evaluated!

Results:
- Final Yield: 0.78
- Work per Yield: 7.2 (lower is better)
- Leaderboard Rank: #1 🎉

See the updated leaderboard:
https://github.com/Livia33g/lammps-state-change/blob/main/problems/problem-001-dimer-ksat/leaderboard.csv

Great job! You beat the baseline!
```

---

### GitHub Issue Submissions

**Step 1: Participant creates issue**
```
Title: [Submission] problem-001: Alice's Adaptive Policy
Body: <paste policy.json content>
```

**Step 2: Extract to inbox**
```bash
mkdir -p submissions-inbox/alice_problem-001/
# Copy from issue → policy.json
# Create submission.json manually
```

**Step 3: Continue with standard workflow** (same as email-based)

---

## 🔒 Security Considerations

### What validate_submission.py checks:

**File size limits:**
- `policy.json`: < 1 MB
- `params.json`: < 100 KB
- `submission.json`: < 100 KB

**Suspicious patterns:**
```python
FORBIDDEN_PATTERNS = [
    r'import\s+os',           # OS module
    r'import\s+subprocess',   # Subprocess
    r'eval\(',                # Code execution
    r'exec\(',                # Code execution
    r'__import__',            # Dynamic imports
    r'system\(',              # System calls
]
```

**Parameter ranges:**
- `flip_probability`: [0.0, 1.0]
- `hysteresis_checks`: [1, 100]
- `contact_cutoff`: [0.1, 5.0]
- `morse_depth_*`: [0.0, 100.0]
- `morse_alpha`: [0.1, 50.0]

### What to do if validation fails:

1. **JSON syntax error** → Email participant with specific line number
2. **Security violation** → Manual review required, contact participant
3. **Parameter out of range** → Email participant with allowed ranges
4. **Schema violation** → Email participant with schema documentation

---

## 📊 Leaderboard Management

### Standard leaderboard (with usernames)

```csv
username,final_yield,work_per_yield,flip_count,cumulative_work,date
alice,0.7800,7.2000,412,5.6160,2026-01-26
baseline,0.7200,8.3000,450,5.9760,2026-01-20
bob,0.6900,9.1000,520,6.2790,2026-01-25
```

### Anonymized leaderboard (ranks only)

```bash
python3 tools/update_leaderboard.py \
    --submission submissions-private/problem-001/alice/ \
    --anonymize
```

```csv
rank,final_yield,work_per_yield,flip_count,cumulative_work,date
1,0.7800,7.2000,412,5.6160,2026-01-26
2,0.7200,8.3000,450,5.9760,2026-01-20
3,0.6900,9.1000,520,6.2790,2026-01-25
```

### Re-running evaluations

If you need to re-evaluate (e.g., with different parameters):

```bash
# Re-submit to cluster
sbatch tools/evaluate_submission.slurm problem-001-dimer-ksat alice_2026-01-26

# Scores will be updated in the same directory
# Re-run leaderboard update to reflect new scores
python3 tools/update_leaderboard.py \
    --submission submissions-private/problem-001/alice_2026-01-26/
```

---

## 🚨 Troubleshooting

### Submission won't validate

**Error:** `Missing required field 'flip_probability'`

**Fix:**
```bash
# Show detailed validation errors
python3 tools/validate_submission.py submissions-inbox/alice_problem-001/

# Email participant with specific error
```

### SLURM job fails

**Error:** LAMMPS compilation fails

**Debug:**
```bash
# Check SLURM output
cat eval_JOBID.out

# Check stderr
cat submissions-private/problem-001/alice/generated/lammps_stderr.log

# Try manual compilation
cd submissions-private/problem-001/alice/generated/
# Follow compilation errors
```

### Leaderboard not updating

**Error:** `scores.json not found`

**Check:**
```bash
# Verify evaluation completed
ls submissions-private/problem-001/alice/generated/analysis/

# If scores.json missing, re-run analysis manually
python3 core/benchmark/score_policy_from_timeseries.py \
    --dump submissions-private/problem-001/alice/generated/dump.*.lammpstrj \
    --thermo submissions-private/problem-001/alice/generated/lammps_stdout.log \
    --output submissions-private/problem-001/alice/generated/analysis/
```

### Invalid scores

**Error:** `final_yield > 1.0`

**Possible causes:**
- Trajectory analysis bug
- Incorrect species type assignment
- LAMMPS simulation crashed

**Fix:**
```bash
# Re-run analysis with different parameters
# Check trajectory file exists and is complete
# Verify simulation completed successfully
```

---

## 💡 Best Practices

### 1. Batch processing schedule

Run `process_new_submissions.sh` on a schedule:

```bash
# Add to crontab (process inbox every hour)
0 * * * * cd /path/to/lammps-state-change && bash tools/process_new_submissions.sh >> logs/processing.log 2>&1
```

### 2. Email templates

Create templates for common scenarios:
- Submission received
- Validation failed
- Evaluation complete
- Error during evaluation

### 3. Backup private repo

Regularly backup `submissions-private/`:

```bash
# Create timestamped archive
tar -czf submissions-backup-$(date +%Y%m%d).tar.gz submissions-private/
```

### 4. Final evaluation with replicas

For final ranking, re-run top N submissions with multiple replicas:

```bash
# Run 5 replicas for top 3 submissions
for user in alice bob charlie; do
    for rep in {1..5}; do
        sbatch tools/evaluate_submission.slurm problem-001 ${user}_replica${rep}
    done
done

# Average scores across replicas
# Update leaderboard with averaged scores + error bars
```

---

## 📚 Related Documentation

- [MODERATOR_GUIDE.md](../MODERATOR_GUIDE.md) - Complete moderator manual
- [PARTICIPANT_GUIDE.md](../PARTICIPANT_GUIDE.md) - For participants
- [core/schemas/](../core/schemas/) - JSON schema specifications
- [problems/](../problems/) - Problem definitions

---

**Questions?** Update this README as you discover best practices!
