# Participant Guide: Advanced Molecular Computing Competition

Welcome to the **advanced** molecular computing competition! This competition explores **what class of problems can molecular computing systems solve in a scalable way** as problem size increases.

## 🎯 Competition Vision

Our approach decomposes molecular computation into three essential stages:

1. **ENCODE** - Problem variables/constraints → molecular components (you design this!)
2. **DESIGN** - Protocol to generate solutions via state-change policy (you design this!)
3. **DECODE** - Nonequilibrium steady state → solution readout (you design this!)

**Your task** is to design the **complete solution**: encoding, state-change policy, and decoding strategy.

---

## 🧬 Physical Substrate: State-Changing Patchy Particles

### The Core Mechanism

We use **patchy particles** (molecules with directional binding sites) that can:
- **Bind** to each other based on patch types (you choose the types!)
- **Flip internal states** (A→B→C) by consuming fuel (you design the rules!)
- **Self-assemble** into clusters via diffusion + binding

The key innovation: particles can **change their state during simulation** based on local observations (who they're touching). This enables non-equilibrium computation beyond equilibrium limits.

### Why Non-Equilibrium?

- **Equilibrium self-assembly** fails for complex structures (yield degrades sharply)
- **Hamiltonian Path Problems** solvable only up to ~19 nodes at equilibrium
- **Non-equilibrium driving** (fuel-powered state changes) overcomes these limits

By carefully designing *when* and *how* particles flip states, you can:
- Escape local energy minima
- Suppress off-target assembly
- Accelerate convergence to solutions

---

## 📋 How to Participate

### Step 1: Choose a Problem

Browse `problems/` directory. Each problem has:
- `problem.json` - Formal problem definition (abstract problem, not encoding!)
- `README.md` - Human-readable description of the abstract problem

Example problems:
- **problem-001-ksat-advanced**: 2-SAT problem (full encode/design/decode via `submission.py`)
- More problems coming soon...

### Step 2: Design Your Solution

For the **advanced competition**, you submit a single Python file that implements three methods:

1. **`encode()`** - Design your own molecular encoding of the problem
   - Choose particle types, geometry, and interactions
   - Create LAMMPS data and input files
   - You have complete freedom in how you encode the abstract problem

2. **`design_policy()`** - Generate C++ code for state-change logic
   - Write C++ fix files that implement your state-change rules
   - Define when and how particles flip states
   - Your fix must compile with our LAMMPS build

3. **`decode()`** - Extract the solution from simulation results
   - Read LAMMPS trajectory/output files
   - Interpret the results according to your encoding
   - Return a solution that satisfies the problem

See `advance/submission_template.py` for the required class structure and `advance/example_submission_2sat.py` for a complete working example.

**Key points:**
- You design the entire encoding (not prescribed like in basic competition)
- Your encoding, policy, and decoding must be coherent (see `advance/COMPILATION_TIPS.md`)
- The example shows ONE possible approach - you're free to use a completely different strategy

### Step 3: Validate Your Submission

```bash
# Check your submission structure
python advance/check_submission.py \
  --submission your_submission.py \
  --problem problems/problem-001-ksat-advanced/problem.json
```

This checks:
- Correct API structure (class name, method signatures)
- `encode()` creates required LAMMPS files
- `design_policy()` generates valid C++ fix files
- `decode()` returns proper format

**Important**: Passing the checker is necessary but not sufficient. Your submission must also pass our sandbox (compilation + smoke simulation).

### Step 4: Submit Your Solution

**Submit via Google Form:**

1. **Open the submission form**: [Submit Your Solution](https://docs.google.com/forms/d/e/1FAIpQLSd6NBCwdrk_zeaSg8l27eCohOHZ_ZJpdmSYZyqaSlzj6QI8qg/viewform?usp=sharing)

2. **Fill out the form**:
   - **Username**: Your username for the leaderboard
   - **Email**: Your email for notifications
   - **Problem ID**: Select the problem you're submitting for
   - **Attempt #**: If this is a resubmission, increment the attempt number
   - **Submission Python**: Paste the **entire contents** of your `submission.py` file

3. **Submit**: Click submit and you'll receive a confirmation

**What happens next:**
- Your submission is automatically organized by problem ID in our tracking sheet
- We process submissions within 24-48 hours
- You'll receive an email when your results are ready (or if there are errors)
- Scores appear on the public leaderboard

**Before submitting, test locally:**
```bash
# Test your submission structure
python advance/check_submission.py \
  --submission your_submission.py \
  --problem problems/problem-001-ksat-advanced/problem.json
```

**Note**: Only your score appears on the leaderboard. Your submission code remains private.

---

## 🏆 Scoring: The Pareto Frontier

Within this model system, we aim to characterize the fundamental tradeoffs between **time, memory, and energy consumption** required to perform a given computational task.

You're optimized on multiple objectives:

### Primary Metrics
1. **Solution Correctness** - Does your decoded solution satisfy the problem?
   *For binary problems: Only solved solutions are ranked*
   *For continuous problems: Ranked by your player-defined score*

2. **Work** - Total fuel consumed (∑ ΔPE or flip count)
   *How energy-efficient is your solution?*

3. **Time to Solution** - Steps to reach satisfying state
   *How fast do you converge?*

### Leaderboard Ranking

**For binary problems** (e.g., 2-SAT):
- Only solutions that satisfy the problem are ranked
- Ranked by total work (lower is better)
- Unsolved submissions appear at the bottom (marked "Not Solved")

**For continuous problems**:
- Ranked by your player-defined score (higher is better)
- You define the scoring metric in your decoding policy

This captures the competition's core question: *Can molecular systems solve problems more efficiently than classical computers on the energy-time-memory Pareto frontier?*

---

## 🔬 Component Coherence: Critical for Success

Your submission has three components that **MUST be coherent**:

1. **`encode()`** - Creates the LAMMPS system (atom types, geometry, interactions)
2. **`design_policy()`** - Generates C++ fix that implements state-change logic
3. **`decode()`** - Reads simulation results to extract solution

### Why Coherence Matters

If these components are **not coherent**, your submission will produce **null/invalid results**:

- ❌ **encode()** creates type 2 for A, but **design_policy()** looks for type 3 → **No state changes occur**
- ❌ **encode()** creates geometry with 3 patches, but **design_policy()** expects 6 patches → **Crashes or wrong behavior**
- ❌ **design_policy()** flips type 2→4, but **decode()** counts type 5 → **Wrong solution measurement**

### How to Ensure Coherence

1. **Define atom type mapping once** and use it consistently across all three methods
2. **Pass metadata from encode() to design_policy()** via the return value
3. **Use same types in decode()** as you created in encode()

See `advance/example_submission_2sat.py` for a complete example demonstrating coherence.

---

## 📚 Resources

- **`advance/README.md`** - Complete guide to the advanced competition format
- **`advance/submission_template.py`** - Required class structure
- **`advance/example_submission_2sat.py`** - Complete working example for 2-SAT problem
- **`advance/COMPILATION_TIPS.md`** - Common sandbox failures and fixes
- **`advance/check_submission.py`** - Local validation tool

---

## 🚨 Common Pitfalls

1. **Missing mass definitions** - LAMMPS requires `mass` commands for all atom types
2. **Incorrect FixStyle registration** - Must use `#ifdef FIX_CLASS` pattern
3. **Component incoherence** - Types/geometry must match between encode/design/decode
4. **File naming errors** - Fix files must start with `fix_state_change_`
5. **Not testing locally** - Always run `check_submission.py` before submitting

See `advance/COMPILATION_TIPS.md` for detailed explanations and fixes.

---

## ❓ FAQ

**Q: Do I have to use the same encoding as the example?**
A: No! The example shows ONE possible approach. You're free to design a completely different encoding.

**Q: Can I use different particle geometries?**
A: Yes! You have complete freedom in your encoding design.

**Q: What if my solution doesn't satisfy the problem?**
A: For binary problems, unsolved solutions appear on the leaderboard but ranked last. You'll receive an email notification. You can resubmit improved versions.

**Q: How do I know if my submission passed?**
A: Check the leaderboard. If your submission appears, it passed all checks. If you receive an email about failures, fix the issues and resubmit.

**Q: Can I see other participants' solutions?**
A: No, submissions remain private. Only scores appear on the leaderboard.

---

Good luck! 🚀
