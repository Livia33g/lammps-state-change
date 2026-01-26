# Simple Daily Workflow

**How to process submissions in 2 commands** 🚀

---

## 📧 When You Get Submissions

### Method 1: Batch Extract (Recommended)

```bash
# 1. Save ALL email attachments to one folder:
~/Downloads/submissions/
  ├── alice/
  │   ├── policy.json
  │   ├── submission.json
  │   └── params.json
  ├── diana/
  │   ├── policy.json
  │   └── submission.json
  └── bob/
      └── ...

# 2. Extract all at once
tools/extract_all_submissions.sh ~/Downloads/submissions problem-001-dimer-ksat

# Output:
#   =========================================
#     Batch Extract Submissions
#   =========================================
#   ----------------------------------------
#   Processing: alice
#   ----------------------------------------
#     ✓ Extracted (with params.json)
#   ----------------------------------------
#   Processing: diana
#   ----------------------------------------
#     ✓ Extracted (no params)
#
#   =========================================
#     Summary
#   =========================================
#     Found:     2
#     Success:   2
#     Failed:    0
#
#   ✓ Extracted 2 submissions to inbox

# 3. Process ONLY the new ones (automatically detects alice + diana)
tools/workflow process-new problem-001-dimer-ksat -j 4
```

**That's it!** 2 commands total.

---

### Method 2: Individual Extract

```bash
# 1. Extract each submission individually
tools/extract_email_submission.sh ~/Downloads/alice/ alice problem-001-dimer-ksat
tools/extract_email_submission.sh ~/Downloads/diana/ diana problem-001-dimer-ksat

# 2. Process ONLY the new ones (automatically detects alice + diana)
tools/workflow process-new problem-001-dimer-ksat -j 4
```

---

## 🔄 What Happens Automatically

When you run `workflow process-new`:

```bash
tools/workflow process-new problem-001-dimer-ksat -j 4
```

**The system automatically:**

1. ✅ **Detects** which submissions are new (haven't been validated yet)
2. ✅ **Skips** all old submissions that are already processed
3. ✅ **Processes** only the new ones through the full pipeline:
   - Validate
   - Generate C++
   - Compile LAMMPS
   - Simulate (4 concurrent jobs)
   - Analyze
   - Cleanup
   - Update leaderboard

**You DON'T need to:**
- ❌ Specify which submissions are new
- ❌ Tell it to skip old ones
- ❌ Manually run each step
- ❌ Track which ones are done

**It's fully automatic!**

---

## 📊 Example: Real Scenario

### Day 1: First Batch (5 submissions)

```bash
# Extract 5 submissions
tools/extract_all_submissions.sh ~/Downloads/week1/ problem-001-dimer-ksat

# Process all 5
tools/workflow process-new problem-001-dimer-ksat -j 4

# Output:
#   Found 5 new submissions
#   Processing new submissions only...
#   Validating: 100%|████| 5/5
#   Generating: 100%|████| 5/5
#   Compiling: 100%|█████| 5/5
#   Simulating: 100%|████| 5/5
#   Analyzing: 100%|█████| 5/5
#   ✓ All actions complete!
```

### Day 3: Second Batch (2 new + Alice resubmits)

```bash
# Extract new submissions (alice_v2, diana, bob)
tools/extract_all_submissions.sh ~/Downloads/week2/ problem-001-dimer-ksat

# Process ONLY the 3 new ones (automatically skips the 5 from Day 1)
tools/workflow process-new problem-001-dimer-ksat -j 4

# Output:
#   Found 3 new submissions
#   Processing new submissions only...
#   [processes only alice_v2, diana, bob]
#   [automatically skips the 5 old submissions]
#   ✓ All actions complete!

# Leaderboard automatically updated:
#   - If alice_v2 is better → shows new score
#   - If alice_v2 is worse → keeps old score
#   - Diana and Bob added
```

---

## 🎯 Complete Daily Routine

### Morning (5 minutes)

```bash
# 1. Check email for submissions
# 2. Save all attachments to ~/Downloads/today/

# 3. Extract all
tools/extract_all_submissions.sh ~/Downloads/today problem-001-dimer-ksat

# 4. Process new ones
tools/workflow process-new problem-001-dimer-ksat -j 8

# Done! Go get coffee while simulations run ☕
```

### Evening (2 minutes)

```bash
# 1. Check status
tools/workflow status problem-001-dimer-ksat

# 2. Commit leaderboard to public repo
cd ~/lammps-state-change
git add problems/*/leaderboard.csv
git commit -m "Update leaderboard: $(date +%Y-%m-%d)"
git push

# 3. Email participants their results (optional)
```

**Total time: 7 minutes per day**
**Everything else: Automated** ✅

---

## 🔍 Check What's New Before Processing

```bash
# See current status
tools/workflow status problem-001-dimer-ksat

# Output shows what's pending:
#   Total submissions: 8
#
#   Action          Complete   Pending    Failed
#   --------------------------------------------------
#   validated       5          3          0
#   generated       5          0          0
#   compiled        5          0          0
#   simulated       5          0          0
#   analyzed        5          0          0
#
#   Next available actions:
#   --------------------------------------------------
#     • validated: 3 submissions ready  ← These are the new ones!
```

The 3 pending validations = 3 new submissions that haven't been processed yet.

---

## 💡 Pro Tips

### Tip 1: Just run `process-new` whenever

```bash
# Run this anytime - it's idempotent
tools/workflow process-new problem-001-dimer-ksat -j 4

# If there are new submissions → processes them
# If there are no new submissions → does nothing
# Never repeats work!
```

### Tip 2: Process multiple problems at once

```bash
# Extract to all problems
tools/extract_all_submissions.sh ~/Downloads/today problem-001-dimer-ksat
tools/extract_all_submissions.sh ~/Downloads/today problem-002-octahedron

# Process both
tools/workflow process-new problem-001-dimer-ksat -j 4
tools/workflow process-new problem-002-octahedron -j 4

# Or use a loop
for problem in problem-001-dimer-ksat problem-002-octahedron; do
    tools/workflow process-new $problem -j 4
done
```

### Tip 3: Set up a cron job (fully automated)

```bash
# Add to crontab:
# Every day at 9 AM, process new submissions
0 9 * * * cd ~/lammps-state-change && tools/workflow process-new problem-001-dimer-ksat -j 8 >> ~/logs/workflow.log 2>&1
```

---

## ❓ FAQ

**Q: Do I need to tell the workflow which submissions are new?**
A: **No!** It auto-detects based on validation status.

**Q: What if I run `process-new` twice?**
A: **Nothing happens.** It's idempotent - won't repeat completed work.

**Q: What if Alice submits 5 times?**
A: **All 5 get processed.** Leaderboard shows her best score automatically.

**Q: Can I process just one specific submission?**
A: **Not directly.** The workflow processes all pending submissions. But you can delete unwanted ones from the inbox before running `process-new`.

**Q: What if a simulation fails?**
A: **Check status** with `workflow status`. The failed submission will show up in the "Failed" column. Fix the issue and run `workflow simulate` again - it will retry only the failed one.

**Q: How do I reprocess everything from scratch?**
A: **Delete the status file:**
```bash
rm ~/lammps-state-change-private/submissions-archive/problem-001-dimer-ksat/.workflow/status.json
tools/workflow run-all problem-001-dimer-ksat -j 4
```

---

## 🚀 Quick Reference

```bash
# Batch extract all submissions
tools/extract_all_submissions.sh ~/Downloads/submissions problem-001-dimer-ksat

# Process only new ones (most common)
tools/workflow process-new problem-001-dimer-ksat -j 4

# Check status
tools/workflow status problem-001-dimer-ksat

# That's all you need! 🎉
```

---

## Summary: You DON'T Need To Specify

✅ **Automatic detection** - The workflow finds new submissions
✅ **Automatic skipping** - Old submissions are never reprocessed
✅ **Automatic leaderboard** - Best scores per user
✅ **Automatic cleanup** - Intermediate files removed

**You only specify:**
- Which submissions to extract from email (batch or individual)
- Which problem to process
- How many concurrent jobs (-j flag)

**Everything else is automatic!** 🎉
