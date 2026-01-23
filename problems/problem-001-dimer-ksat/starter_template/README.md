# Starter Template for Problem 001

This directory contains scaffolding to help you get started with your submission.

## 📁 Files in This Template

- `policy_template.json` - Copy this and fill in your strategy
- `params_template.json` - Copy this and tune parameters
- `submission_template.json` - Metadata for your submission

## 🚀 Quick Start

### Step 1: Copy Templates

```bash
# Create your submission directory
mkdir -p ../../submissions/problem-001/your_username/

# Copy templates
cp policy_template.json ../../submissions/problem-001/your_username/policy.json
cp params_template.json ../../submissions/problem-001/your_username/params.json
cp submission_template.json ../../submissions/problem-001/your_username/submission.json
```

### Step 2: Design Your Policy

Edit `policy.json`:

```json
{
  "policy_name": "my_clever_strategy",
  "check_frequency": 100,
  "state_transitions": [
    {
      "from_species": "A",
      "to_species": "C",
      "trigger": {
        "contact_required": {
          "species": "B",
          "cutoff": 0.7,    // ← Tune this
          "min_contacts": 1
        }
      },
      "flip_probability": 0.9,  // ← Try <1.0 for exploration
      "hysteresis_checks": 5     // ← Adjust for stability
    }
  ]
}
```

### Step 3: Tune Parameters (Optional)

Edit `params.json`:

```json
{
  "morse_depth_AB": 1.2,  // Stronger A-B binding
  "morse_depth_BC": 0.9,  // Weaker C-B binding
  "morse_alpha": 5.0
}
```

### Step 4: Validate Your Submission

```bash
# From repo root
python tools/validate_submission.py submissions/problem-001/your_username/
```

This checks:
- ✅ JSON syntax is valid
- ✅ Policy conforms to schema
- ✅ Parameters are in allowed ranges
- ✅ All required files present

### Step 5: Generate C++ Code (Optional - we do this server-side)

```bash
python core/generators/generate_fix_from_policy.py \
  submissions/problem-001/your_username/policy.json \
  submissions/problem-001/your_username/generated/
```

This creates:
- `fix_state_change_my_clever_strategy.h`
- `fix_state_change_my_clever_strategy.cpp`

### Step 6: Submit

**Option A: Git commit + push**
```bash
git add submissions/problem-001/your_username/
git commit -m "Add submission for problem-001"
git push
```

**Option B: Share directory**
Zip and send `submissions/problem-001/your_username/` to us.

We'll run your simulation on our cluster and update the leaderboard!

---

## 💡 Strategy Ideas

### Beginner: Tune Baseline

Start with `baseline_policy.json` and just change:
- `flip_probability`: Try 0.8, 0.9, 0.95
- `hysteresis_checks`: Try 3, 5, 10
- `cutoff`: Try 0.6, 0.7, 0.8

Goal: Beat baseline `work_per_yield = 8.3`

### Intermediate: Add Reversibility

Allow C→A to escape kinetic traps:

```json
{
  "from_species": "C",
  "to_species": "A",
  "trigger": {
    "contact_required": {
      "species": "B",
      "cutoff": 0.7,
      "min_contacts": 0,
      "logic": "NOT"  // Flip when NOT touching B
    }
  },
  "flip_probability": 0.1,
  "hysteresis_checks": 10
}
```

### Advanced: Multi-Stage Policy

Different strategies early vs late:

**Early game:** Aggressive conversion (high probability)
**Late game:** Refinement (low probability, high hysteresis)

*Note: Current schema doesn't support time-dependent policies yet, but you can approximate with tuned parameters.*

---

## 🧪 Local Testing (If You Have LAMMPS)

### Build LAMMPS with Your Fix

```bash
# Copy generated C++ to LAMMPS src/
cp generated/fix_state_change_*.{h,cpp} /path/to/lammps/src/

# Recompile
cd /path/to/lammps/src/
make yes-RIGID
make yes-MOLECULE
make mpi -j8
```

### Generate System

```bash
cd problems/problem-001-dimer-ksat/
python ../../tools/generate_system.py \
  --problem problem.json \
  --policy ../../submissions/problem-001/your_username/policy.json \
  --params ../../submissions/problem-001/your_username/params.json \
  --output test_run/
```

### Run Simulation

```bash
cd test_run/
mpirun -np 4 lmp_mpi -in in.dimer_ksat_1core
```

### Analyze Results

```bash
python ../../analyze_trajectory_target_yield_and_work.py \
  --dump dump.dimer_ksat_1core.lammpstrj \
  --thermo lammps_stdout.log \
  --events stderr.log \
  --site-types 2 3 4 \
  --bond-cutoff 0.7 \
  --yield-mode species_fraction \
  --species-label 4 \
  --out analysis/my_test
```

Check `analysis/my_test.csv` for scores!

---

## ❓ FAQ

**Q: Can I use machine learning to find the best parameters?**
A: Absolutely! Train offline, then encode your learned policy in JSON.

**Q: What if I want to use different logic than contact-based?**
A: Contact us! We're building more trigger types (energy-based, composition-based, time-based).

**Q: How many iterations can I submit?**
A: Unlimited! We track all versions and show your best score on leaderboard.

**Q: Can I see others' policies?**
A: After the initial leaderboard is established (1 month), all policies become public for learning.

---

## 📚 Resources

- **Problem README**: `../README.md` - Full problem description
- **Participant Guide**: `../../PARTICIPANT_GUIDE.md` - Competition overview
- **State Change Explained**: `../../docs/state_change_explained.md` - Technical deep dive
- **Baseline Performance**: `../baseline_policy.json` - What to beat

**Need help?** Open an issue or join the discussion forum!

Good luck! 🚀
