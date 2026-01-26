# Quick Start: Accepting Submissions

**Your questions answered:** How do people submit? Do you have to manually process emails?

---

## ✅ Quick Answers

**Q: Does the private repo exist yet?**
**A:** No, you need to create it (takes 2 minutes). See [Step 1](#step-1-create-private-repository-2-minutes) below.

**Q: How do people submit?**
**A:** Email you their JSON files as attachments (simplest method).

**Q: Do they submit to the private repo?**
**A:** No, they don't have access. They email YOU → you extract to private repo.

**Q: Would you have to manually go through each email?**
**A:** Semi-manual: You run a script that extracts attachments → then batch processing is fully automated.

---

## 🚀 5-Minute Setup

### Step 1: Create Private Repository (2 minutes)

```bash
# On GitHub web interface:
# 1. Click "New repository"
# 2. Name: lammps-state-change-private
# 3. Visibility: PRIVATE ⚠️ (this is critical!)
# 4. Initialize with README: Yes
# 5. Create repository

# Then on your machine:
cd ~
git clone git@github.com:YourUsername/lammps-state-change-private.git
cd lammps-state-change-private

# Create directory structure
mkdir -p submissions-inbox submissions-archive moderator-notes

# Add to your ~/.bashrc
echo 'export SUBMISSIONS_PRIVATE_REPO="$HOME/lammps-state-change-private"' >> ~/.bashrc
source ~/.bashrc

# Commit
git add .
git commit -m "Initial structure"
git push
```

**✅ Done! Private repo ready.**

### Step 2: Set Up Submission Email (3 minutes)

Create a template email to send to participants:

```
Subject: How to Submit - Problem 001

Hi!

To submit your solution:

1. Email to: your-email@university.edu
2. Subject: [SUBMISSION] problem-001 - YourUsername
3. Attach these files:
   - policy.json (required)
   - submission.json (required)
   - params.json (optional)

Example files: https://github.com/YourOrg/lammps-state-change/tree/main/problems/problem-001-dimer-ksat/starter_template/

You'll get results within 48 hours!
```

**✅ Done! Ready to accept submissions.**

---

## 📧 Processing Submissions (5 minutes per submission)

### When you receive a submission email:

```bash
# 1. Save email attachments to Downloads folder
#    (Most email clients: right-click attachment → Save All)
#    Example: saves to ~/Downloads/

# 2. Extract to private repo inbox
cd ~/lammps-state-change
bash tools/extract_email_submission.sh \
    ~/Downloads/ \
    alice \
    problem-001-dimer-ksat

# Output:
#   =========================================
#     Extract Email Submission
#   =========================================
#   ✓ Found all required files
#   ✓ Submission extracted to inbox
#
#   Next steps:
#     1. Review submission files (optional)
#     2. Run: bash tools/process_new_submissions.sh
#     3. Monitor cluster jobs

# 3. Process ALL submissions in inbox (automated from here)
bash tools/process_new_submissions.sh

# This automatically:
#   ✓ Validates all submissions
#   ✓ Moves to archive
#   ✓ Submits to SLURM cluster
#   ✓ Generates C++, compiles LAMMPS, runs simulation

# 4. Wait for cluster jobs to finish
squeue -u $USER

# 5. Update leaderboard when done
python3 tools/update_leaderboard.py \
    --submission ~/lammps-state-change-private/submissions-archive/problem-001-dimer-ksat/alice_2026-01-26_*/

# 6. Commit leaderboard to PUBLIC repo
cd ~/lammps-state-change
git add problems/problem-001-dimer-ksat/leaderboard.csv
git commit -m "Update leaderboard: alice"
git push

# 7. Email alice her results
```

**Time breakdown:**
- Save attachments: 30 seconds
- Run extraction script: 10 seconds
- Run batch processor: 10 seconds (then automated)
- Update leaderboard: 30 seconds
- Total manual time: ~2 minutes per submission

**Everything else is automated!**

---

## 🔄 Batch Processing Multiple Submissions

If you get 5 emails with submissions:

```bash
# 1. Save all attachments to separate folders:
#    ~/Downloads/alice_submission/
#    ~/Downloads/bob_submission/
#    ~/Downloads/charlie_submission/
#    ~/Downloads/diana_submission/
#    ~/Downloads/eve_submission/

# 2. Extract all to inbox
cd ~/lammps-state-change
bash tools/extract_email_submission.sh ~/Downloads/alice_submission alice problem-001-dimer-ksat
bash tools/extract_email_submission.sh ~/Downloads/bob_submission bob problem-001-dimer-ksat
bash tools/extract_email_submission.sh ~/Downloads/charlie_submission charlie problem-001-dimer-ksat
bash tools/extract_email_submission.sh ~/Downloads/diana_submission diana problem-001-dimer-ksat
bash tools/extract_email_submission.sh ~/Downloads/eve_submission eve problem-001-dimer-ksat

# 3. Process ALL at once
bash tools/process_new_submissions.sh

# Processes all 5 automatically!
# Output:
#   =========================================
#     Submission Batch Processor
#   =========================================
#
#   Found 5 submission(s) to process
#
#   ----------------------------------------
#   Processing: alice_problem-001-dimer-ksat
#   ----------------------------------------
#     ✓ Validation passed
#     ✓ Submitted (Job ID: 12345)
#
#   ----------------------------------------
#   Processing: bob_problem-001-dimer-ksat
#   ----------------------------------------
#     ✓ Validation passed
#     ✓ Submitted (Job ID: 12346)
#
#   ... (continues for all 5)
#
#   =========================================
#     Summary
#   =========================================
#     Processed: 5
#     Failed:    0

# Total manual time: ~10 minutes for 5 submissions
# Everything else: automated
```

---

## 🤖 Automation Level

| Task | Manual? | Time |
|------|---------|------|
| Check email for submissions | Manual | 1 min |
| Save attachments to folder | Manual | 30 sec |
| Run extraction script | Semi-auto | 10 sec |
| Validate submissions | **Automated** | instant |
| Generate C++ fix | **Automated** | instant |
| Compile LAMMPS | **Automated** | ~5 min |
| Run simulation | **Automated** | 1-24 hrs |
| Analyze results | **Automated** | ~1 min |
| Update leaderboard | Semi-auto | 30 sec |
| Email results | Manual | 2 min |

**Total manual time per submission: ~4 minutes**

**Everything in between: Fully automated!**

---

## 💡 Tips to Save Time

### Tip 1: Email Filters

Set up Gmail filter:
- If subject contains `[SUBMISSION]`
- Label as "Competition"
- Mark as important

Now you can quickly see new submissions.

### Tip 2: Batch Weekly

Instead of processing each submission immediately:
- Collect for 1 week
- Process all on Friday afternoon
- Update leaderboard once
- Send results Monday morning

Reduces context switching!

### Tip 3: Template Emails

Save result email templates:

**Success email:**
```
Subject: Results Ready - Problem 001

Hi {USERNAME},

Your submission has been evaluated!

Results:
- Final Yield: {YIELD}
- Work per Yield: {WORK} (lower is better)
- Rank: #{RANK}

Full leaderboard: {LEADERBOARD_URL}

Submit again anytime to improve!
```

**Failure email:**
```
Subject: Submission Error - Problem 001

Hi {USERNAME},

Your submission failed validation:

{ERROR_MESSAGE}

Please fix and resubmit!
See examples: {STARTER_TEMPLATE_URL}
```

Use find/replace to fill in values!

---

## 📊 Example Weekly Schedule

**Monday:**
- Announce new problem
- Email submission instructions to participants

**Tuesday-Thursday:**
- Check email once per day
- Extract new submissions (5 min)
- Run batch processor

**Friday:**
- Process all week's submissions
- Monitor cluster jobs
- Update leaderboards
- Prepare result emails

**Weekend:**
- Send result emails
- Archive in private repo

**Total moderator time: ~1 hour per week for 10-20 submissions**

---

## ❌ What You DON'T Have To Do

You **do NOT** need to:

- ❌ Write any C++ code
- ❌ Manually compile LAMMPS for each submission
- ❌ Manually run simulations
- ❌ Manually analyze trajectories
- ❌ Manually calculate scores
- ❌ Manually update leaderboard entries
- ❌ Give participants access to private repo
- ❌ Review participant code (JSON is validated automatically)

**All of this is automated!**

---

## ✅ What You DO Have To Do

You **do** need to:

- ✅ Check email for submissions
- ✅ Save attachments to folder (30 sec)
- ✅ Run extraction script (10 sec)
- ✅ Run batch processor (10 sec)
- ✅ Update leaderboard when jobs finish (30 sec)
- ✅ Email results to participants (2 min)

**Total: ~4 minutes per submission**

---

## 🎯 Summary

**Setup once:** Create private repo (2 minutes)

**Per submission:**
1. Save email attachments → folder
2. Run: `bash tools/extract_email_submission.sh ...`
3. Run: `bash tools/process_new_submissions.sh`
4. Wait for cluster (automated)
5. Run: `python3 tools/update_leaderboard.py ...`
6. Email results

**That's it!** The system handles everything else automatically.

---

## 📚 More Details

- **All submission methods:** [SUBMISSION_WORKFLOW.md](SUBMISSION_WORKFLOW.md)
- **Security setup:** [SECURITY_SETUP.md](SECURITY_SETUP.md)
- **Tools reference:** [tools/README.md](tools/README.md)

**Ready to start accepting submissions!** 🚀
