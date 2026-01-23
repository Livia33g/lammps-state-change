# Participant Guide: Molecular Computing Competition

Welcome to the molecular computing benchmark! This competition explores **what class of problems can molecular computing systems solve in a scalable way** as problem size increases.

## 🎯 Competition Vision

Our approach decomposes molecular computation into three essential stages:

1. **ENCODE** - Problem variables/constraints → molecular components
2. **DESIGN** - Protocol to generate solutions (your contribution!)
3. **DECODE** - Nonequilibrium steady state → solution readout

**Your task** is to design the **optimal state-change policy** that solves each problem efficiently.

---

## 🧬 Physical Substrate: State-Changing Patchy Particles

### The Core Mechanism

We use **patchy particles** (molecules with directional binding sites) that can:
- **Bind** to each other based on patch types (A, B, C, etc.)
- **Flip internal states** (A→B→C) by consuming fuel
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
- `problem.json` - Formal problem definition
- `README.md` - Human-readable description
- `starter_template/` - Scaffold to get started

Example problems:
- **problem-001-dimer-ksat**: Convert A→C via B catalyst (species conversion)
- **problem-002-octahedron**: Assemble exact octahedral structures
- **problem-003-hamiltonian-path**: Find valid paths in graphs
- **problem-004-frustrated-sampling**: Sample from frustrated systems

### Step 2: Design Your Policy

Create a **policy.json** file that specifies state transition rules:

```json
{
  "policy_name": "my_clever_strategy",
  "policy_version": "1.0.0",
  "check_frequency": 100,

  "state_transitions": [
    {
      "from_species": "A",
      "to_species": "C",
      "trigger": {
        "contact_required": {
          "species": "B",
          "cutoff": 1.6,
          "min_contacts": 1
        }
      },
      "flip_probability": 1.0,
      "hysteresis_checks": 5
    }
  ]
}
```

**Key parameters to tune:**
- `flip_probability` - How aggressively to flip when conditions met
- `hysteresis_checks` - How many consecutive checks before committing (prevents oscillation)
- `cutoff` - Distance threshold for detecting contacts
- Trigger logic - Which contacts trigger which flips

### Step 3: Validate Your Submission

```bash
# Check your policy is valid
python tools/validate_submission.py submissions/problem-001/your_username/

# Generate C++ code from your policy
python core/generators/generate_fix_from_policy.py \
  submissions/problem-001/your_username/policy.json \
  submissions/problem-001/your_username/generated/

# Compile LAMMPS with your fix (we do this server-side)
./tools/compile_fix.sh submissions/problem-001/your_username/

# Run evaluation (we do this server-side on cluster)
python tools/run_evaluation.py submissions/problem-001/your_username/ \
  --problem problem-001 --replicas 1
```

### Step 4: Submit

Create a directory:
```
submissions/problem-001/your_username/
├── policy.json          # Your state transition rules
├── submission.json      # Metadata
└── params.json          # Tunable parameters (optional)
```

We'll run your simulation on our cluster and update the leaderboard.

---

## 🏆 Scoring: The Pareto Frontier


Within this model system, we aim to characterize the fundamental tradeoffs between **time, memory, and energy consumption** required to perform a given computational task.

You're optimized on multiple objectives:

### Primary Metrics
1. **Yield** - Fraction of molecules in target configuration
   *How well did you solve the problem?*

2. **Work** - Total fuel consumed (∑ ΔPE or flip count)
   *How energy-efficient is your solution?*

3. **Time to Threshold** - Steps to reach target yield
   *How fast do you converge?*

### Leaderboard Ranking

Most problems use **work per yield** as primary metric:
```
Score = (Fuel consumed) / (Final yield)
Lower is better!
```

This captures the competition's core question: *Can molecular systems solve problems more efficiently than classical computers on the energy-time-memory Pareto frontier?*

---

## 🔬 The Theory Behind Your Policy

### Mathematical Framework

Your policy implements the **perturbation function π(j|s,H)** from stochastic thermodynamics:

```
𝓛(φ,θ) = p₀(x₀) ∏ₜ π_φ(jₜ|xₜ) · exp(-βH_θ(bₜ₊₁;sₜ₊₁))/Z_θ(sₜ₊₁)
         \_____/   \___________/   \________________________________/
          Initial      Policy           Equilibrium relaxation
```

Where:
- **π_φ** = your policy (which state to flip given local observations)
- **H_θ** = binding landscape (interaction energies)
- **Z_θ** = partition function (equilibrium physics)

The factorization separates:
- **Equilibrium** (what structures are thermodynamically favored)
- **Non-equilibrium** (how you actively drive the system)

Your policy.json defines π_φ - the decision rule for state changes.

### Policy Design Principles

From our preliminary work on 3-SAT (proposal Figure 3):

1. **Random policy** (flip any bound variable) - baseline performance
2. **Locally informed policy** (bias toward unsatisfied clauses) - 2-3× better
3. **Optimized policy** (gradient-based learning) - near-optimal (future work!)

Good policies balance:
- **Exploration** (try new configurations)
- **Exploitation** (reinforce promising paths)
- **Hysteresis** (don't thrash between states)

---

## 🛠️ Technical Details

### LAMMPS Implementation

Your policy.json auto-generates a C++ LAMMPS fix:

```
fix state_change patches state/change/your_policy nevery cutoff pflip ...
```

This fix:
1. **Measures** local environment every `nevery` steps
2. **Perturbs** according to your transition rules
3. **Relaxes** via LAMMPS molecular dynamics between checks

The generated code handles:
- Contact detection with periodic boundaries
- Hysteresis tracking per molecule
- Probabilistic flips with RNG
- MPI communication (ghost atoms)

### Validation Checks

We automatically verify:
- **Schema compliance** (policy.json format)
- **Energy conservation** (no unphysical behavior)
- **Molecule count** (atoms don't disappear)
- **Topology validity** (rigid bodies stay intact)

### Reproducibility

All evaluations use:
- Fixed random seeds
- Identical initial configurations
- Standardized analysis scripts
- LAMMPS version consistency

Your submission's generated C++ code is version-controlled for transparency.

---

## 📚 Resources

### Understanding State Changes

Read `docs/state_change_explained.md` for deep dive on:
- The unfix→change→refix cycle
- Why it's necessary for rigid bodies
- How coordination detection works
- Common pitfalls and solutions

### Example Policies

See `problems/*/baseline_solutions/` for reference implementations:
- `greedy/` - Simple greedy policies
- `random/` - Random baseline for comparison
- `optimized/` (future) - Gradient-optimized policies

### Benchmark Paper

Our proposal (included in repo) describes:
- Theoretical foundations
- Benchmark problem suite
- Scoring methodology
- Long-term research vision

---

## 🎓 Learning Path

### Beginner Track
1. Start with **problem-001** (simple species conversion)
2. Try the random baseline policy
3. Incrementally improve flip_probability and hysteresis
4. Compare your score to baseline

### Intermediate Track
1. Tackle **problem-002** (octahedron assembly)
2. Design multi-stage policies (exploration → exploitation)
3. Tune interaction parameters (morse_depth, etc.)
4. Optimize for work_per_yield

### Advanced Track
1. Attempt **problem-003** (Hamiltonian path)
2. Design compositional policies (different rules for different states)
3. Implement feedback loops (previous flip success → adjust probability)
4. Beat equilibrium baselines from our paper!

---

## ❓ FAQ

**Q: Can I write custom C++ instead of JSON?**
A: Advanced mode coming soon! For now, JSON keeps policies comparable.

**Q: How many replicas are averaged?**
A: Default is 1 run. Top leaderboard entries re-evaluated with N=5 for final ranking.

**Q: Can I use machine learning to optimize?**
A: Absolutely! Train offline, then encode your learned policy in JSON.

**Q: What if my simulation crashes?**
A: Check validation output. Common issues: cutoff too large, probability too high, hysteresis too low.

**Q: How do I debug my policy?**
A: Run locally with `--verbose` flag. Check `stderr` for `STATECHANGE` event logs.

**Q: Can I see others' policies?**
A: After submission deadline, all policies become public for learning.

---

## 🌟 Competition Goals

Our ambitious goal is to discover problems where we can clearly show through theory and computation that **molecular computation dominates standard computing protocols in the Pareto Frontier**, in the scaling sense as problem size grows.

By participating, you're helping answer:
1. Which problems are uniquely suited to molecular computing?
2. What design principles lead to scalable solutions?
3. How do energy-time-memory tradeoffs compare to silicon/quantum systems?

Your innovations could help establish molecular computing as a viable technology!

---

## 📬 Support

- **Issues**: Report bugs at `github.com/your-repo/issues`
- **Discussions**: Join `discussions/` for strategy sharing
- **Documentation**: See `docs/` for technical deep dives

**Good luck, and happy molecular designing!** 🧪
