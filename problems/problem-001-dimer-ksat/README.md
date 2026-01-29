# Problem 001: 2-SAT Dimer Problem (A,B,C and D,E,F channels)

**Difficulty:** Beginner | **Estimated Time:** 2-4 hours | **Category:** 2-SAT Dimer

---

## 🎯 Problem Statement

Design a **state-change policy** that efficiently maximizes the concentration of true particles (C and F) in a 2-SAT dimer system with two independent conditions, while minimizing energy consumption (state change events).

### Physical System
- **50 patchy monomers** (each: 1 core + 3 patches for A,C,E,F, or 1 core + 6 patches for B,D)
- **Two independent channels:**
  - **Channel 1 (A,B,C):** 10 A, 5 B, 10 C initially
  - **Channel 2 (D,E,F):** 10 E, 5 D, 10 F initially
- **Allowed reactions:**
  - Channel 1: A,A on both B faces → flip higher ID A→C
  - Channel 2: E on D (any face) → flip E→F
- **Your goal:** Maximize C+F concentration while minimizing flip count

### Success Criteria
- **Primary metric:** `concentration_of_true_particles = (n_C + n_F) / (n_A + n_C + n_E + n_F)` (maximize)
- **Threshold:** Reach ≥70% true particle concentration within 2M timesteps
- **Scoring:** Higher concentration is better

---

## 📖 Theoretical Background

This problem implements a **2-SAT formula** with two independent conditions:

### Condition B (TF - True OR False)
- **Monomer B**: Condition monomer with two patch faces
- **Monomer A**: False (switchable) - can flip to C
- **Monomer C**: True (non-switchable)
- **Rule**: Allows only 1 false (A) monomer to be attached at a time
  - If **A, A** are attached simultaneously to B (on both faces), flip the one with **higher molecule ID** to C (true)
  - If **C, C** or **A, C** are attached, nothing happens

### Condition D (TT - True AND True)
- **Monomer D**: Condition monomer with two patch faces
- **Monomer E**: False (switchable) - can flip to F
- **Monomer F**: True (non-switchable)
- **Rule**: Both need to be true (F)
  - If even **one false (E)** attaches to D, switch it to type F (true), regardless of what is attached to the other side

### Independence
The two sets **A,B,C** and **D,E,F** are **independent** and do **not interact** with each other.

### The Solution
**If we set B to True (ON), both rules are satisfied:**
- **Rule 1 (B=TF)**: (Anything ∨ True) = True ✓
- **Rule 2 (D=TT)**: (Anything ∨ True) = True ✓

---

## 🧪 Physical Details

### Geometry
Each monomer has:
- **1 core atom** (type 1 for A,B,C channel, type 12 for D,E,F channel, mass 0.6) - excluded volume repulsion
- **3 patch atoms** (A,C,E,F: types 2,4,8,10, mass 0.1 each) - binding sites
- **6 patch atoms** (B,D: types 3,5,9,11, mass 0.1 each) - two faces with 3 patches each

**Patch positions** (relative to core at origin):
```
Face +x: (0.5, 0.000, 0.100)
         (0.5, 0.087, -0.05)
         (0.5, -0.087, -0.05)

Face -x: (-0.5, 0.000, -0.100)
         (-0.5, 0.087, 0.05)
         (-0.5, -0.087, 0.05)
```

### Interactions

**Channel ABC attractions:**
- **A-B_face1** (type 2-3): Morse potential
- **A-B_face2** (type 2-5): Morse potential
- **C-B_face1** (type 4-3): Morse potential
- **C-B_face2** (type 4-5): Morse potential

**Channel EFD attractions:**
- **E-D_face1** (type 8-9): Morse potential
- **E-D_face2** (type 8-11): Morse potential
- **F-D_face1** (type 10-9): Morse potential
- **F-D_face2** (type 10-11): Morse potential

**Default:** All Morse depths = 1.0, alpha = 5.0, r0 = 0.0

**All other pairs:** Neutral (D0=0) except core-core repulsion (soft potential)

**Channels are independent:** No cross-channel interactions

---

## 🚀 Getting Started

### 1. Understand the State-Change Rules

**Channel 1 (A,B,C):**
- Trigger: A,A attached to B on both faces simultaneously
- Action: Flip the A molecule with **higher molecule ID** to C
- Hysteresis: Require N consecutive checks before flipping

**Channel 2 (D,E,F):**
- Trigger: E attached to D (any face, even one)
- Action: Flip E to F
- Hysteresis: Require N consecutive checks before flipping

### 2. Design Your Policy

See `baseline_policy.json` for a reference implementation. Key parameters to tune:

**Trigger:**
- `cutoff` - How close must patches be to trigger flip? (0.5-1.0)
  - Should match morse_rcut ≈ 0.7 for consistency

**Flip behavior:**
- `flip_probability` - Chance of flip when triggered (0-1)
  - 1.0: Deterministic (flip every time)
  - <1.0: Stochastic (explore alternative paths)

- `hysteresis_checks` - Consecutive checks required (1-20)
  - Low: Responsive but might oscillate
  - High: Stable but slow to react

**Channel-specific logic:**
- Channel 1: Must detect A on **both** B faces (face1 AND face2)
- Channel 2: Detects E on **any** D face (face1 OR face2)

### 3. Example Policy Structure

```json
{
  "state_transitions": [
    {
      "from_species": "A",
      "to_species": "C",
      "trigger": {
        "contact_required": {
          "species": "B_face1",
          "cutoff": 0.7,
          "min_contacts": 1,
          "logic": "AND",
          "also_requires": {
            "species": "B_face2",
            "cutoff": 0.7,
            "min_contacts": 1
          }
        }
      },
      "flip_probability": 1.0,
      "hysteresis_checks": 5,
      "flip_higher_id": true
    },
    {
      "from_species": "E",
      "to_species": "F",
      "trigger": {
        "contact_required": {
          "species": ["D_face1", "D_face2"],
          "cutoff": 0.7,
          "min_contacts": 1,
          "logic": "OR"
        }
      },
      "flip_probability": 1.0,
      "hysteresis_checks": 5
    }
  ]
}
```

---

## 📊 Evaluation

### What We Measure

**From trajectory analysis:**
- `concentration_of_true_particles` - (n_C + n_F) / (n_A + n_C + n_E + n_F) at end
- `total_state_changes` - Total A→C and E→F events (from stderr logs)
- `cumulative_work` - Sum of ΔPE between thermo outputs
- `time_to_solution` - Steps to reach 70% true particle concentration

**Your score:**
```
concentration_of_true_particles = (n_C + n_F) / (n_A + n_C + n_E + n_F)
```

Higher is better! The solution B=True should maximize this concentration.

### Baseline Comparison

| Policy | Concentration | State Changes | Time to 70% |
|--------|--------------|--------------|-------------|
| Random | 0.50 | 5 | Never |
| Greedy | 0.75 | 15 | 1.5M steps |
| **Your Goal** | >0.80 | <20 | <1M steps |

---

## 🧠 Design Principles

The non-equilibrium policy balances exploration vs. exploitation:

1. **Hysteresis prevents thrashing**
   - Without it: A flips to C, immediately flips back
   - With it: Wait for stable contact before committing

2. **Probability enables exploration**
   - 100%: Deterministic, finds local optimum fast
   - <100%: Stochastic, explores alternatives, may find global optimum

3. **Contact cutoff matters**
   - Must match interaction range (morse_rcut ≈ 7 × patch_radius = 0.7)
   - Too large: Spurious flips from distant patches
   - Too small: Miss valid contacts

4. **Channel independence**
   - Channels operate independently
   - Optimize each channel separately or together

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
  "submission_date": "2026-01-30",
  "policy_version": "v1.0",
  "description": "Optimized two-channel policy with adaptive hysteresis",
  "expected_performance": {
    "concentration": 0.85,
    "total_state_changes": 12
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

❌ **Cutoff too large** → Flips from non-bonded patches → wasted energy
✅ Set `cutoff ≈ morse_rcut` for consistency

❌ **Hysteresis too low** → Oscillation (A↔C thrashing)
✅ Use hysteresis ≥ 5 checks

❌ **Missing channel logic** → Channel 1 needs BOTH faces, Channel 2 needs ANY face
✅ Implement correct trigger logic for each channel

### Debugging

Check stderr for flip events:
```bash
grep "STATECHANGE" slurm_*.err | wc -l  # Count flips
```

Visualize trajectory:
```bash
python analyze_trajectory_target_yield_and_work.py \
  --dump dump.dimer_ksat_1core_twosideB_twins.lammpstrj \
  --thermo lammps_stdout.log \
  --events slurm_*.err \
  --site-types 2 3 4 5 8 9 10 11 \
  --bond-cutoff 0.7 \
  --yield-mode concentration_of_true \
  --species-labels 4 10 \
  --out analysis/my_run
```

---

## 🏁 Ready to Compete!

1. **Read** the baseline policy (`baseline_policy.json`)
2. **Modify** flip probability, hysteresis, or cutoff
3. **Validate** your submission
4. **Run** evaluation (we'll do this server-side)
5. **Check** leaderboard for your ranking!

**Goal:** Beat the greedy baseline (`concentration > 0.75`)

Good luck! 🚀
