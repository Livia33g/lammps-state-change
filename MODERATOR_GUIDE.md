# Moderator Guide: Running the Competition

This guide is for competition organizers/moderators who process submissions and maintain the leaderboard.

---

## 📋 Overview

**Public Repo** (this repo): Problems, framework, documentation, leaderboard
**Private Repo** (recommended): Actual submissions (not public)
**Your Cluster**: Where evaluations run

---

## 🔄 Submission Workflow

### Option A: Email-Based (Simplest)

1. **Participant emails** their files to `competition@yourlab.edu`
2. **You receive:**
   - `policy.json`
   - `params.json` (optional)
   - `submission.json` (metadata)

3. **You process:**
```bash
# Create submission directory
mkdir -p submissions-private/problem-001/alice_2026-01-23/
cd submissions-private/problem-001/alice_2026-01-23/

# Add their files
cp ~/Downloads/policy.json .
cp ~/Downloads/params.json .
cp ~/Downloads/submission.json .

# Validate
python3 ../../tools/validate_submission.py .

# Submit to cluster
sbatch ../../tools/evaluate_submission.slurm problem-001 alice_2026-01-23
```

4. **Wait for results**, then update leaderboard

### Option B: GitHub Issues

Participants submit via GitHub Issue:
- Issue title: `[Submission] problem-001: Alice's Policy`
- Issue body: Paste `policy.json` content
- You copy to private repo and process

### Option C: Google Form

Create form with:
- Team name
- Problem ID (dropdown)
- Upload policy.json
- Upload params.json (optional)

Responses go to private Google Drive → you download and process

---

## 🛠️ Processing a Submission

### Step 1: Validate

```bash
python3 tools/validate_submission.py submissions-private/problem-001/username/
```

**Checks:**
- Required files present
- Valid JSON syntax
- Policy structure correct
- No suspicious content
- Parameters in reasonable ranges

**If validation fails:** Email participant with error message

### Step 2: Evaluate on Cluster

```bash
sbatch tools/evaluate_submission.slurm problem-001 username
```

**This automatically:**
1. Generates C++ fix from policy
2. Generates LAMMPS data/input files
3. Compiles LAMMPS with custom fix
4. Runs simulation (up to 24 hours)
5. Analyzes trajectory
6. Calculates scores

### Step 3: Check Results

```bash
# Check SLURM output
cat eval_JOBID.out

# Check scores
cat submissions-private/problem-001/username/generated/analysis/scores.json
```

Example `scores.json`:
```json
{
  "final_yield": 0.78,
  "work_per_yield": 7.2,
  "flip_count": 412,
  "time_to_threshold": 980000,
  "cumulative_work": 5.62
}
```

### Step 4: Update Leaderboard

**Manual (for now):**
```bash
# Edit problems/problem-001/leaderboard.csv
echo "alice,0.78,7.2,412,2026-01-23" >> problems/problem-001/leaderboard.csv

# Commit to public repo
cd /path/to/lammps-state-change
git add problems/problem-001/leaderboard.csv
git commit -m "Update leaderboard: alice submission"
git push
```

**Automated (future):**
```bash
python3 tools/update_leaderboard.py \
  --problem problem-001 \
  --username alice \
  --scores submissions-private/problem-001/alice/generated/analysis/scores.json
```

---

## 📊 Leaderboard Format

`problems/problem-001/leaderboard.csv`:
```csv
username,final_yield,work_per_yield,flip_count,date
baseline,0.72,8.3,450,2026-01-20
alice,0.78,7.2,412,2026-01-23
bob,0.81,6.8,380,2026-01-24
```

**Anonymized option:**
```csv
rank,final_yield,work_per_yield,flip_count,date
1,0.81,6.8,380,2026-01-24
2,0.78,7.2,412,2026-01-23
3,0.72,8.3,450,2026-01-20
```

---

## 🔒 Keeping Submissions Private

### **DO NOT:**
- ❌ Commit submissions to public repo
- ❌ Share policy.json files between participants
- ❌ Reveal strategies before competition ends

### **DO:**
- ✅ Use a separate private repo or local directory
- ✅ Add `submissions/` to `.gitignore` (already done)
- ✅ Only publish leaderboard scores (not policies)
- ✅ After competition: optionally make top solutions public

### Private Repo Structure

**Recommended: Create `lammps-state-change-submissions` (private repo)**

Access:
- ✅ You (organizer)
- ✅ Other moderators
- ❌ Participants (no access)

```
lammps-state-change-submissions/
├── problem-001/
│   ├── alice_2026-01-23/
│   ├── bob_2026-01-24/
│   └── charlie_2026-01-25/
├── problem-002/
└── archive/
    └── finished-competitions/
```

---

## 🚨 Troubleshooting

### Submission won't validate
- Check JSON syntax (use jsonlint.com)
- Verify all required fields present
- Email participant with specific error

### LAMMPS compilation fails
- Check C++ fix for syntax errors
- Verify LAMMPS version compatibility
- Try compiling manually to see full error

### Simulation crashes
- Check `stderr.log` for LAMMPS errors
- Common issues:
  - Atoms lost (bad initial config)
  - Energy NaN (unstable parameters)
  - Timeout (increase SLURM time limit)

### Scores seem wrong
- Verify dump file exists
- Check thermo output is complete
- Re-run analysis script manually:
  ```bash
  python3 core/benchmark/score_policy_from_timeseries.py \
    --dump dump.001_dimer_ksat.lammpstrj \
    --thermo lammps_stdout.log
  ```

---

## 📅 Competition Timeline

**Phase 1: Soft Launch (Week 1-2)**
- Announce competition
- Accept first submissions
- Debug any issues
- Build up baseline leaderboard

**Phase 2: Active Competition (Week 3-6)**
- Regular submissions
- Update leaderboard weekly
- Answer participant questions
- Monitor for issues

**Phase 3: Final Evaluation (Week 7)**
- Submission deadline
- Re-evaluate top N submissions with multiple replicas
- Announce winners
- Optionally publish top solutions

---

## 🔧 Automation Scripts

### Create these for full automation:

**`tools/process_new_submissions.sh`**
```bash
#!/bin/bash
# Process all new submissions in inbox/
for dir in submissions-inbox/*/; do
    username=$(basename $dir)
    problem=$(cat $dir/submission.json | jq -r .problem_id)

    echo "Processing: $username for $problem"

    # Validate
    python3 tools/validate_submission.py $dir
    if [ $? -eq 0 ]; then
        # Move to private repo
        mv $dir submissions-private/$problem/$username/

        # Submit to cluster
        sbatch tools/evaluate_submission.slurm $problem $username
    else
        echo "FAILED validation: $username"
        # Email participant
    fi
done
```

**`tools/update_leaderboard.py`**
```python
#!/usr/bin/env python3
import json, csv, sys
from pathlib import Path

def update_leaderboard(problem, username, scores_json):
    # Load scores
    with open(scores_json) as f:
        scores = json.load(f)

    # Update CSV
    leaderboard_path = Path(f"problems/{problem}/leaderboard.csv")
    # ... append new row ...
    # ... sort by primary metric ...
    # ... write back ...

    print(f"Updated leaderboard for {username}")

if __name__ == '__main__':
    update_leaderboard(sys.argv[1], sys.argv[2], sys.argv[3])
```

---

## 📧 Communication Templates

### Submission Received
```
Subject: Submission Received - problem-001

Hi [name],

We've received your submission for problem-001 and added it to the evaluation queue.

Estimated processing time: 6-24 hours
You'll receive an email when results are ready.

Thanks for participating!
```

### Validation Failed
```
Subject: Submission Validation Error - problem-001

Hi [name],

Your submission for problem-001 failed validation:

Errors:
- policy.json: Missing 'flip_probability' field
- params.json: Invalid JSON syntax on line 5

Please fix these issues and resubmit.

Need help? Reply to this email.
```

### Results Ready
```
Subject: Evaluation Complete - problem-001

Hi [name],

Your submission has been evaluated!

Results:
- Final Yield: 0.78
- Work per Yield: 7.2 (lower is better)
- Leaderboard Rank: #2

See the updated leaderboard: [link]

Want to improve? Submit a new version anytime!
```

---

## 💡 Tips for Moderators

1. **Set clear deadlines** - Weekly or monthly submission windows
2. **Communicate often** - Update participants on queue status
3. **Start small** - Test with 5-10 participants before scaling
4. **Archive everything** - Keep all submissions even after competition
5. **Be responsive** - Answer questions quickly to keep engagement high
6. **Celebrate wins** - Announce new leaderboard leaders
7. **Plan for replicas** - Re-run top submissions with N=5 for final ranking

---

## 🎓 Advanced: Multi-Problem Tournaments

For multiple problems, consider:
- **Aggregate score** across all problems
- **Weighted leaderboard** (harder problems worth more)
- **Pareto frontier** visualization (yield vs energy vs time)
- **Categories** (beginner, advanced, expert)

---

**Questions?** Update this guide as you learn! Make this your ops manual.
