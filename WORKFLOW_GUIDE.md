# Workflow Management System

**Signac-flow style workflow for competition submissions** - Process submissions efficiently with automatic tracking, concurrency, and cleanup.

---

## 🎯 Key Features

✅ **Action-based** - Separate commands for each step (validate, generate, compile, simulate, analyze, cleanup)
✅ **Idempotent** - Never repeat work already done
✅ **Concurrent** - Run multiple simulations in parallel
✅ **Progress tracking** - See what's pending/running/complete
✅ **Automatic cleanup** - Delete intermediate files after scoring
✅ **Smart leaderboard** - Keep only best submission per user
✅ **Process only new** - Skip already-processed submissions

---

## 🚀 Quick Start

### 1. Check Status

```bash
# See current state of all submissions
tools/workflow status problem-001-dimer-ksat
```

**Output:**
```
================================================================================
  Workflow Status: problem-001-dimer-ksat
================================================================================

Total submissions: 5

Action          Complete   Pending    Failed
--------------------------------------------------
validated       3          2          0
generated       2          1          0
compiled        2          0          0
simulated       1          1          0
analyzed        1          0          0

Next available actions:
--------------------------------------------------
  • validated: 2 submissions ready
  • generated: 1 submissions ready
  • simulated: 1 submissions ready
```

### 2. Process Everything

```bash
# Process all pending actions automatically
tools/workflow run-all problem-001-dimer-ksat
```

Or step-by-step:

```bash
tools/workflow validate problem-001-dimer-ksat
tools/workflow generate problem-001-dimer-ksat
tools/workflow compile problem-001-dimer-ksat
tools/workflow simulate problem-001-dimer-ksat -j 4  # 4 concurrent simulations
tools/workflow analyze problem-001-dimer-ksat
tools/workflow cleanup problem-001-dimer-ksat
tools/workflow update-leaderboard problem-001-dimer-ksat
```

### 3. Process Only New Submissions

```bash
# When you get new submissions, process only those (skip old ones)
tools/workflow process-new problem-001-dimer-ksat
```

---

## 📋 All Commands

| Command | Description | Example |
|---------|-------------|---------|
| `status` | Show current workflow state | `tools/workflow status problem-001` |
| `validate` | Validate pending submissions | `tools/workflow validate problem-001` |
| `generate` | Generate C++ fixes | `tools/workflow generate problem-001` |
| `compile` | Compile LAMMPS | `tools/workflow compile problem-001` |
| `simulate` | Run simulations | `tools/workflow simulate problem-001 -j 8` |
| `analyze` | Analyze results | `tools/workflow analyze problem-001` |
| `cleanup` | Remove intermediate files | `tools/workflow cleanup problem-001` |
| `update-leaderboard` | Update public leaderboard | `tools/workflow update-leaderboard problem-001` |
| `run-all` | Execute all pending actions | `tools/workflow run-all problem-001` |
| `process-new` | Process only new submissions | `tools/workflow process-new problem-001` |

---

## 🔄 Complete Workflow Example

### Receive New Submissions

```bash
# 1. Extract email submissions to inbox
tools/extract_email_submission.sh ~/Downloads/alice/ alice problem-001-dimer-ksat
tools/extract_email_submission.sh ~/Downloads/bob/ bob problem-001-dimer-ksat
tools/extract_email_submission.sh ~/Downloads/charlie/ charlie problem-001-dimer-ksat
```

### Check What Needs Processing

```bash
# 2. Check status
tools/workflow status problem-001-dimer-ksat

# Output:
#   Total submissions: 8
#   validated: 3 complete, 3 pending
#   generated: 3 complete, 0 pending
#   compiled: 3 complete, 0 pending
#   simulated: 2 complete, 1 pending
#   analyzed: 2 complete, 0 pending
#
#   Next available actions:
#     • validated: 3 submissions ready
#     • simulated: 1 submission ready
```

### Process New Submissions Only

```bash
# 3. Process only the 3 new submissions (alice, bob, charlie)
#    Skips the 5 already-processed submissions automatically
tools/workflow process-new problem-001-dimer-ksat -j 4
```

**This will:**
1. Validate the 3 new submissions
2. Generate C++ for validated ones
3. Compile LAMMPS (sequential to avoid conflicts)
4. Run simulations (4 concurrent jobs)
5. Analyze results
6. Clean up intermediate files
7. Update leaderboard with best scores

### Monitor Progress

```bash
# While simulations run, check status
tools/workflow status problem-001-dimer-ksat

# Output:
#   simulated: 2 complete, 3 running
```

### When Complete, Update Leaderboard

```bash
# Leaderboard is already updated by process-new
# But you can manually update if needed:
tools/workflow update-leaderboard problem-001-dimer-ksat

# Output:
#   Updating leaderboard...
#   Found 6 users with submissions
#   Updating: 100%|████████████████| 6/6
#   ✓ Leaderboard updated with best scores
```

---

## 🎛️ Advanced Usage

### Concurrent Simulations

```bash
# Run 8 simulations in parallel
tools/workflow simulate problem-001-dimer-ksat -j 8

# Default is 4 if not specified
tools/workflow simulate problem-001-dimer-ksat
```

### Reprocess Failed Submissions

If a submission failed during simulate, fix the issue and reprocess:

```bash
# The workflow system tracks status, so you can just run simulate again
# It will only process the failed ones (idempotent)
tools/workflow simulate problem-001-dimer-ksat
```

### Clean Up Disk Space

```bash
# Remove trajectory files, C++ sources, and other intermediate files
# Keeps: policy.json, scores.json, analysis results
tools/workflow cleanup problem-001-dimer-ksat

# This frees up ~90% of disk space per submission
```

### Process Everything at Once

```bash
# Process all pending actions in dependency order
tools/workflow run-all problem-001-dimer-ksat -j 8

# Equivalent to running each action individually:
#   validate → generate → compile → simulate → analyze → cleanup → update-leaderboard
```

---

## 📊 Status Tracking

The workflow system maintains state in:

```
~/lammps-state-change-private/
└── submissions-archive/
    └── problem-001-dimer-ksat/
        ├── .workflow/
        │   ├── status.json          # Current status of all submissions
        │   ├── all_scores.json      # All scores (not just best)
        │   └── slurm_jobs.json      # SLURM job IDs (if using SLURM)
        ├── alice_2026-01-26_143022/
        │   ├── policy.json
        │   ├── submission.json
        │   └── generated/
        │       ├── analysis/
        │       │   └── scores.json  # Alice's score
        │       └── lammps_stderr.log # Events log (kept for debugging)
        └── bob_2026-01-27_091534/
```

**status.json format:**
```json
{
  "alice_2026-01-26_143022": {
    "submission_id": "alice_2026-01-26_143022",
    "username": "alice",
    "timestamp": "2026-01-26_143022",
    "problem_id": "problem-001-dimer-ksat",
    "validated": "complete",
    "generated": "complete",
    "compiled": "complete",
    "simulated": "complete",
    "analyzed": "complete",
    "cleaned": "complete",
    "score": 7.2,
    "error_message": null,
    "created_at": "2026-01-26T14:30:22",
    "completed_at": "2026-01-26T18:45:10"
  },
  "bob_2026-01-27_091534": {
    "submission_id": "bob_2026-01-27_091534",
    "username": "bob",
    "timestamp": "2026-01-27_091534",
    "problem_id": "problem-001-dimer-ksat",
    "validated": "complete",
    "generated": "complete",
    "compiled": "complete",
    "simulated": "running",
    "analyzed": "pending",
    "cleaned": "pending",
    "score": null,
    "error_message": null,
    "created_at": "2026-01-27T09:15:34",
    "completed_at": null
  }
}
```

---

## 🧹 Cleanup Details

When you run `tools/workflow cleanup problem-001`, it removes:

**Removed files (can regenerate):**
- `dump.*.lammpstrj` - Trajectory files (~500 MB each)
- `fix_state_change_*.cpp` - C++ source (can regenerate from policy.json)
- `fix_state_change_*.h` - C++ header (can regenerate)
- `data.*.lammps` - Data file (can regenerate)
- `in.*.lammps` - Input script (can regenerate)
- `lammps_stdout.log` - LAMMPS stdout (not needed after analysis)

**Kept files:**
- `policy.json` - Original submission
- `submission.json` - Metadata
- `params.json` - Parameters
- `lammps_stderr.log` - Event log (useful for debugging)
- `analysis/scores.json` - Final scores
- `analysis/*.png` - Visualizations

**Disk savings:** ~90% reduction per submission

---

## 🏆 Leaderboard Management

### Multiple Submissions Per User

When a user submits multiple times, the workflow keeps track of all submissions:

```bash
# Alice submits 3 times
alice_2026-01-26_143022/  # score: 8.5
alice_2026-01-27_091534/  # score: 7.2  ← best
alice_2026-01-28_152341/  # score: 7.8
```

**Leaderboard shows only best:**

```bash
tools/workflow update-leaderboard problem-001-dimer-ksat

# Public leaderboard (problems/problem-001/leaderboard.csv) shows:
# username,score,...
# alice,7.2,...        ← best score only
```

**All scores tracked internally:**

The workflow maintains `all_scores.json` with complete history:

```json
{
  "alice": [
    {"score": 8.5, "submission_id": "alice_2026-01-26_143022", "timestamp": "2026-01-26_143022"},
    {"score": 7.2, "submission_id": "alice_2026-01-27_091534", "timestamp": "2026-01-27_091534"},
    {"score": 7.8, "submission_id": "alice_2026-01-28_152341", "timestamp": "2026-01-28_152341"}
  ]
}
```

### Force Leaderboard Update

```bash
# Manually trigger leaderboard update
tools/workflow update-leaderboard problem-001-dimer-ksat

# This:
# 1. Finds all analyzed submissions
# 2. Groups by username
# 3. Keeps best score per user
# 4. Updates public leaderboard in public repo
```

---

## 🔧 Integration with SLURM

The workflow automatically detects if SLURM is available:

**With SLURM:**
```bash
tools/workflow simulate problem-001-dimer-ksat

# Submits all pending simulations to SLURM
# Tracks job IDs in .workflow/slurm_jobs.json
# Returns immediately (jobs run on cluster)
```

**Without SLURM:**
```bash
tools/workflow simulate problem-001-dimer-ksat -j 4

# Runs locally with 4 concurrent processes
# Blocks until all simulations complete
```

**Check SLURM job status:**
```bash
squeue -u $USER

# Or check workflow status
tools/workflow status problem-001-dimer-ksat
```

---

## 💡 Best Practices

### Daily Workflow

**Morning: Check for new submissions**
```bash
# Extract email submissions
tools/extract_email_submission.sh ~/Downloads/alice/ alice problem-001

# Check status
tools/workflow status problem-001-dimer-ksat

# Process only new ones
tools/workflow process-new problem-001-dimer-ksat -j 8
```

**Evening: Clean up and update**
```bash
# Clean up completed submissions
tools/workflow cleanup problem-001-dimer-ksat

# Update leaderboard
tools/workflow update-leaderboard problem-001-dimer-ksat

# Commit leaderboard to public repo
cd ~/lammps-state-change
git add problems/*/leaderboard.csv
git commit -m "Update leaderboard: $(date +%Y-%m-%d)"
git push
```

### Batch Processing

```bash
# Process all problems at once
for problem in problem-001-dimer-ksat problem-002-octahedron problem-003-hamiltonian; do
    echo "Processing $problem..."
    tools/workflow process-new $problem -j 4
done
```

### Monitor Long-Running Simulations

```bash
# Check status periodically
watch -n 60 'tools/workflow status problem-001-dimer-ksat'

# Or with SLURM
watch -n 60 'squeue -u $USER'
```

---

## 🐛 Troubleshooting

### Submission Stuck in "running"

```bash
# Check if SLURM job is actually running
squeue -u $USER

# If job failed, re-run simulate
tools/workflow simulate problem-001-dimer-ksat

# Check error messages
cat ~/lammps-state-change-private/.workflow/status.json | grep error_message
```

### Failed Validation

```bash
# See which submissions failed
tools/workflow status problem-001-dimer-ksat

# Check error in status.json
cat .workflow/status.json

# Fix the submission and re-validate
tools/workflow validate problem-001-dimer-ksat
```

### Reprocess Everything

```bash
# Delete status.json to reset (WARNING: this will reprocess everything)
rm ~/lammps-state-change-private/submissions-archive/problem-001-dimer-ksat/.workflow/status.json

# Re-run workflow
tools/workflow run-all problem-001-dimer-ksat
```

---

## 📊 Comparison: Old vs New Workflow

| Task | Old Workflow | New Workflow |
|------|-------------|--------------|
| Process new submissions | Manual script per submission | `workflow process-new` |
| Check status | Check files manually | `workflow status` |
| Avoid reprocessing | Manual tracking | Automatic (idempotent) |
| Concurrent simulations | Manual SLURM submission | `-j 8` flag |
| Cleanup | Manual file deletion | `workflow cleanup` |
| Multiple submissions | Manual leaderboard logic | Automatic (keeps best) |
| Progress tracking | None | Progress bars + status |

---

## 🚀 Quick Reference

```bash
# Most common commands

# Check what needs to be done
workflow status problem-001-dimer-ksat

# Process everything
workflow run-all problem-001-dimer-ksat -j 8

# Process only new submissions
workflow process-new problem-001-dimer-ksat -j 8

# Clean up disk space
workflow cleanup problem-001-dimer-ksat

# Update leaderboard
workflow update-leaderboard problem-001-dimer-ksat
```

---

**That's it!** The workflow system handles all the complexity for you. 🎉
