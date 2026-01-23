# Problem 001: Catalytic Species Conversion (A→C via B)

**Difficulty:** Beginner | **Estimated Time:** 2-4 hours | **Category:** Species Conversion

---

## 🎯 Problem Statement

Design a **state-change policy** that efficiently converts A-labeled monomers to C-labeled monomers using B as a catalyst, while minimizing fuel consumption.

### Physical System
- **30 patchy monomers** (each: 1 core + 3 patches)
- **Initial composition:** 20 A-monomers, 10 B-monomers, 0 C-monomers
- **Allowed reaction:** A + B → C + B (B acts as catalyst)
- **Your goal:** Maximize C-yield while minimizing flip count

### Success Criteria
- **Primary metric:** `work_per_yield = cumulative_work / final_yield` (minimize)
- **Threshold:** Reach ≥60% C-yield within 2M timesteps
- **Scoring:** Lower work per yield is better (more efficient)

---

## 📖 Connection to Research Proposal

This problem implements the **3-SAT encoding** from the proposal (Section: Benchmark of Computational Challenges):

> *"Variables are encoded as particles with two internal states (''true'' and ''false''); clauses are encoded as particles with three binding sites, each corresponding to one variable. [...] An unsatisfied clause particle can consume fuel to flip the state of a bound variable from unsatisfying to satisfying."*

### Molecular Encoding

| Proposal Concept | Implementation |
|-----------------|----------------|
| **Variable (true/false)** | Not explicitly represented (simplified) |
| **Clause (satisfied/unsatisfied)** | A-patch (unsatisfied), C-patch (satisfied) |
| **Catalyst/Enzyme** | B-patch (flips A→C when in contact) |
| **Fuel consumption** | Each A→C flip event |
| **Target distribution** | High fraction of satisfied clauses (C) |

### Measurement → Perturbation → Relaxation

From Table 1 in the proposal:

| Component | Implementation |
|-----------|----------------|
| **Measurement 𝓜** | Detect when A-patch is within cutoff of B-patch |
| **Perturbation π** | Flip A→C with probability `pflip` (your policy) |
| **Relaxation ℛ** | LAMMPS rigid-body dynamics between checks |

Your policy defines **π(j|s,H)** - the decision rule for when and how to flip states.

---

## 🧪 Physical Details

### Geometry
Each monomer has:
- **1 core atom** (type 1, mass 0.6) - excluded volume repulsion
- **3 patch atoms** (type 2/3/4, mass 0.1 each) - binding sites

**Patch positions** (relative to core at origin):
```
Patch 1: (0.5, 0.000, 0.100)
Patch 2: (0.5, 0.087, -0.05)
Patch 3: (0.5, -0.087, -0.05)
```

### Interactions

**Only two attractive pairs:**
- **B-A** (type 3-2): Morse potential with depth `D0_AB` (catalyst binds unsatisfied)
- **B-C** (type 3-4): Morse potential with depth `D0_BC` (catalyst binds satisfied)

**Default:** `D0_AB = D0_BC = 1.0` (equal affinity)

**All other pairs:** Neutral (D0=0) except core-core repulsion (soft potential)

### Why This Design?

From the proposal:

> *"Since Hopfield's foundational work on kinetic proofreading, it has been understood that allosteric interactions enhance computational capacity beyond equilibrium bounds, provided an energy cost is paid."*

The **non-equilibrium driving** (A→C flips) allows the system to:
1. Escape local minima where A's are stuck
2. Bias toward higher C-yield than equilibrium allows
3. Trade fuel (flip events) for computational progress

---

## 🚀 Getting Started

### 1. Understand the Baseline Policy

See `baseline_policy.json`:

```json
{
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

**Logic:** When A is within 0.7 units of B for 5 consecutive checks, flip A→C with 100% probability.

**Performance:** `final_yield = 0.72`, `work_per_yield = 8.3`

### 2. Design Your Policy

Key parameters to tune:

**Trigger:**
- `cutoff` - How close must B be to trigger flip? (0.5-1.0)
  - Too small: Misses opportunities
  - Too large: False positives

**Flip behavior:**
- `flip_probability` - Chance of flip when triggered (0-1)
  - 1.0: Greedy (flip every time)
  - <1.0: Stochastic (explore alternative paths)

- `hysteresis_checks` - Consecutive checks required (1-20)
  - Low: Responsive but might oscillate
  - High: Stable but slow to react

### 3. Advanced Strategies

**Reversible flips** (optional):
```json
{
  "from_species": "C",
  "to_species": "A",
  "trigger": {
    "contact_required": {
      "species": "B",
      "cutoff": 0.7,
      "min_contacts": 0,
      "logic": "NOT"
    }
  },
  "flip_probability": 0.1,
  "hysteresis_checks": 10
}
```

**Logic:** Allow C→A if C is NOT touching B (escape kinetic traps)

**Parameter tuning:**
Modify `params.json`:
```json
{
  "morse_depth_AB": 1.5,  // Stronger A-B binding
  "morse_depth_BC": 0.8,  // Weaker C-B binding
  "morse_alpha": 7.0      // Sharper potential
}
```

**Effect:** If `D_AB > D_BC`, system favors A-B contacts → more flip opportunities.

---

## 📊 Evaluation

### What We Measure

**From trajectory analysis:**
- `final_yield` - Fraction of C at end
- `flip_count` - Total A→C events (from stderr logs)
- `cumulative_work` - Sum of ΔPE between thermo outputs
- `time_to_threshold` - Steps to reach 60% C

**Your score:**
```
work_per_yield = cumulative_work / final_yield
```

Lower is better! This captures the **energy-yield tradeoff** from the proposal's Pareto frontier.

### Baseline Comparison

| Policy | Final Yield | Work/Yield | Time to 60% |
|--------|-------------|------------|-------------|
| Random | 0.45 | 15.2 | Never |
| Greedy | 0.72 | 8.3 | 1.2M steps |
| **Your Goal** | >0.80 | <6.0 | <1M steps |

---

## 🧠 Design Principles

From the proposal (Section: Optimization Methods):

> *"The non-equilibrium policy can be optimized with policy gradients. [...] Good policies balance exploration (try new configurations) vs. exploitation (reinforce promising paths)."*

### Key Insights

1. **Hysteresis prevents thrashing**
   - Without it: A flips to C, immediately flips back
   - With it: Wait for stable contact before committing

2. **Probability enables exploration**
   - 100%: Greedy, finds local optimum fast
   - <100%: Stochastic, explores alternatives, may find global optimum

3. **Contact cutoff matters**
   - Must match interaction range (morse_rcut ≈ 7 × patch_radius = 0.7)
   - Too large: Spurious flips from distant B's
   - Too small: Miss valid contacts

4. **Reversibility helps escape traps**
   - If all A→C too quickly, system gets stuck
   - Small C→A probability allows rearrangement

---

## 📁 Submission Format

Create a directory `submissions/problem-001/your_username/`:

```
submissions/problem-001/your_username/
├── submission.json      # Metadata (name, date, description)
├── policy.json          # Your state transition rules
└── params.json          # Tunable parameters (optional)
```

### Example submission.json

```json
{
  "problem_id": "problem-001-dimer-ksat",
  "team_name": "YourUsername",
  "submission_date": "2026-01-23",
  "policy_version": "v1.0",
  "description": "Adaptive hysteresis with probabilistic exploration",
  "expected_performance": {
    "final_yield": 0.85,
    "work_per_yield": 5.5
  }
}
```

---

## 🔬 Understanding State Changes in LAMMPS

**Critical reading:** `docs/state_change_explained.md`

The **unfix → change → refix cycle** is necessary because:
1. LAMMPS rigid fixes cache atom-to-molecule mappings
2. Changing atom types breaks these mappings
3. Temporarily removing the fix, changing types, then recreating the fix solves this

**Your policy.json automatically generates C++ code that implements this correctly!**

---

## 💡 Tips & Tricks

### Common Mistakes

❌ **Cutoff too large** → Flips from non-bonded B's → wasted fuel
✅ Set `cutoff ≈ morse_rcut` for consistency

❌ **Hysteresis too low** → Oscillation (A↔C thrashing)
✅ Use hysteresis ≥ 5 checks

❌ **100% probability always** → Gets stuck in local minima
✅ Try 80-90% for better exploration

### Debugging

Check stderr for flip events:
```bash
grep "STATECHANGE" slurm_*.err | wc -l  # Count flips
```

Visualize trajectory:
```bash
python analyze_trajectory_target_yield_and_work.py \
  --dump dump.dimer_ksat_1core.lammpstrj \
  --thermo lammps_stdout.log \
  --events slurm_*.err \
  --site-types 2 3 4 \
  --bond-cutoff 0.7 \
  --yield-mode species_fraction \
  --species-label 4 \
  --out analysis/my_run
```

Check `analysis/my_run.png` for yield vs. time curve.

---

## 📚 Further Reading

**Proposal sections:**
- Section 3 (Physical Primitives) - Why non-equilibrium?
- Section 4 (Benchmark Challenges) - 3-SAT encoding
- Section 5 (Theory) - π(φ,θ) optimization framework
- Figure 3 - Policy comparison on 3-SAT instances

**Code references:**
- `dimer_ksat/variants/1core/generate.py` - System generator
- `dimer_ksat/variants/1core/fix_state_change_dimer_ksat.cpp` - Reference C++ fix
- `benchmark/score_policy_from_timeseries.py` - Scoring implementation

---

## 🏁 Ready to Compete!

1. **Read** the baseline policy (`baseline_policy.json`)
2. **Modify** flip probability, hysteresis, or cutoff
3. **Validate** your submission
4. **Run** evaluation (we'll do this server-side)
5. **Check** leaderboard for your ranking!

**Goal:** Beat the greedy baseline (`work_per_yield < 8.3`)

Good luck! 🚀
