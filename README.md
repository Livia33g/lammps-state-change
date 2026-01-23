# Molecular Computing Competition Framework

**A Kaggle-style benchmark for out-of-equilibrium molecular computation**

## 🎯 What Is This?

This repository implements a **competition framework** for designing non-equilibrium molecular computing systems using state-changing patchy particles in LAMMPS.

### The Core Question

**Which class of problems can molecular computing systems solve in a scalable way as problem size increases?**

Participants design **state-change policies** (when and how molecules flip between states) to solve computational challenges while optimizing the **energy-time-memory Pareto frontier**.

---

## 🚀 Quick Start

### For Participants

1. **Read the guide**: [PARTICIPANT_GUIDE.md](PARTICIPANT_GUIDE.md)
2. **Choose a problem**: Browse `problems/` directory
3. **Design your policy**: Fill in `policy.json` (declarative, no coding!)
4. **Submit**: We run it on our cluster, update leaderboard

### For Problem Creators

1. **Define encoding**: Create `problem.json` (species, interactions, scoring)
2. **Write description**: Problem README connecting to proposal framework
3. **Provide baseline**: Reference solution for comparison
4. **Set up analysis**: Configure benchmarking pipeline

---

## 📂 Repository Structure

```
lammps-state-change/
├── README.md                          # This file
├── PARTICIPANT_GUIDE.md               # Comprehensive participant guide
├── FRAMEWORK_SUMMARY.md               # Implementation details
│
├── problems/                          # Competition challenges
│   ├── problem-001-dimer-ksat/        # A→C catalytic conversion
│   │   ├── problem.json               # Problem definition
│   │   ├── README.md                  # Problem description
│   │   ├── baseline_policy.json       # Reference solution
│   │   ├── baseline_params.json       # Default parameters
│   │   └── starter_template/          # Scaffolding for participants
│   │       ├── README.md              # Quick start guide
│   │       ├── policy_template.json   # Fill this in!
│   │       ├── params_template.json   # Tune these values
│   │       └── submission_template.json
│   │
│   ├── problem-002-octahedron/        # (Future) Structure assembly
│   └── problem-003-hamiltonian-path/  # (Future) Graph search
│
├── submissions/                       # User submissions (gitignored)
│   └── problem-001/
│       └── username/
│           ├── policy.json            # State transition rules
│           ├── params.json            # Tunable parameters
│           └── submission.json        # Metadata
│
├── core/                              # Framework infrastructure
│   ├── schemas/                       # JSON schemas
│   │   ├── policy_schema.json         # Policy specification
│   │   ├── problem_schema.json        # Problem definition
│   │   └── example_*.json             # Examples
│   │
│   ├── generators/                    # Code generators
│   │   └── generate_fix_from_policy.py   # JSON → C++ converter
│   │
│   ├── benchmark/                     # Scoring system
│   │   ├── run_task.py
│   │   ├── score_policy_from_timeseries.py
│   │   └── aggregate_leaderboard.py
│   │
│   └── analysis/                      # Analysis tools
│       └── analyze_trajectory_target_yield_and_work.py
│
├── examples/                          # Educational examples
│   ├── dimer/                         # Simple 2-species system
│   ├── octahedron/                    # Rigid-body assembly
│   ├── dimer_ksat/                    # k-SAT variants
│   └── inverse-design/                # Equilibrium optimization tools
│
├── docs/                              # Documentation
│   ├── STATE_CHANGE_EXPLANATION.md    # How state changes work
│   ├── INSTRUCTIONS.md                # Rigid-body setup guide
│   ├── FIXES_APPLIED.md               # MPI safety notes
│   └── ...
│
└── tools/                             # Helper scripts
    ├── rebuild_manual.sh              # Build LAMMPS with fixes
    ├── add_new_fix.sh                 # Add custom fixes
    └── (validation scripts - future)
```

---

## 🧬 The Framework: Encode / Design / Decode

From the research proposal, molecular computation has three essential stages:

### 1. **ENCODE** (Problem Creator)
Map problem variables and constraints onto molecular components:
- **Species labels** (A, B, C) → LAMMPS atom types
- **Interactions** (binding affinities) → Morse/Soft potentials
- **Initial composition** → Number of each species

**Defined in:** `problem.json`

### 2. **DESIGN** (Participant)
Protocol to generate solutions via non-equilibrium driving:
- **State transitions** (A→B→C) based on local observations
- **Flip probabilities** (deterministic vs stochastic)
- **Hysteresis** (stability vs responsiveness)

**Defined in:** `policy.json` (your submission!)

### 3. **DECODE** (Automated)
Extract computational answer from molecular state:
- **Yield measurement** (fraction in target configuration)
- **Work calculation** (fuel consumed, ΔPE)
- **Scoring** (Pareto frontier: yield vs energy vs time)

**Defined in:** `problem.json` + `benchmark/score_*.py`

---

## 🏆 Competition Metrics

Following the proposal's focus on **fundamental tradeoffs between time, memory, and energy**, we score on multiple objectives:

### Primary Metric (Most Problems)
```
work_per_yield = cumulative_work / final_yield
```
Lower is better! Measures **energy efficiency** of your policy.

### Secondary Metrics
- **Final yield** - Did you solve the problem?
- **Time to threshold** - How fast did you converge?
- **Flip count** - Total fuel consumed (proxy for thermodynamic work)

### Leaderboard Ranking
Policies are ranked on the **Pareto frontier**:
- Best yield at given energy budget
- Fastest time to threshold
- Lowest energy for target yield

No single "winner" - multiple optimal strategies in different regimes!

---

## 📊 Current Problems

### Problem 001: Catalytic Species Conversion
**Difficulty:** Beginner | **Status:** Active

Convert A→C using B as catalyst. Simplified k-SAT encoding testing basic non-equilibrium fuel-driven computation.

- **Encoding:** 30 monomers (20 A, 10 B, 0 C initially)
- **Goal:** Maximize C-yield, minimize flips
- **Baseline:** Greedy policy achieves `yield=0.72`, `work/yield=8.3`
- **Your challenge:** Beat the baseline!

[→ Problem 001 Details](problems/problem-001-dimer-ksat/README.md)

### Problem 002: Octahedron Assembly
**Difficulty:** Intermediate | **Status:** Coming Soon

Exact structure assembly with proofreading. From proposal Figure 2.

### Problem 003: Hamiltonian Path
**Difficulty:** Advanced | **Status:** Coming Soon

Graph search encoding. From proposal Figure 2 (equilibrium limit: 19 nodes max).

---

## 🔬 Technical Innovation: Declarative Policies

**No C++ coding required!** Participants submit JSON policies that we auto-generate into production C++ LAMMPS fixes.

### Example Policy

```json
{
  "policy_name": "my_strategy",
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
      "flip_probability": 0.9,
      "hysteresis_checks": 5
    }
  ]
}
```

**Automatically generates:**
- MPI-safe C++ LAMMPS fix
- Contact detection with PBC
- Hysteresis tracking
- Event logging

---

## 🎓 Educational Resources

### Getting Started
1. [PARTICIPANT_GUIDE.md](PARTICIPANT_GUIDE.md) - Full walkthrough
2. [problems/problem-001-dimer-ksat/README.md](problems/problem-001-dimer-ksat/README.md) - First challenge
3. [docs/STATE_CHANGE_EXPLANATION.md](docs/STATE_CHANGE_EXPLANATION.md) - How it works under the hood

### Understanding State Changes
- **The unfix→change→refix cycle**: Why it's necessary for rigid bodies
- **Hysteresis**: Preventing oscillation between states
- **MPI safety**: Ghost atoms and local-only flips
- **Contact detection**: PBC-aware neighbor searches

### Theory Background
- **Research proposal** (included in repo) - Mathematical framework
- **π(φ,θ) optimization** - Policy gradients and landscape design
- **Stochastic thermodynamics** - Work, entropy production, Pareto frontiers

---

## 🛠️ For Developers

### Building LAMMPS with Custom Fixes

```bash
# Generate C++ from policy
python core/generators/generate_fix_from_policy.py \
  submissions/problem-001/username/policy.json \
  output_dir/

# Copy to LAMMPS source
cp output_dir/fix_state_change_*.{h,cpp} /path/to/lammps/src/

# Compile
cd /path/to/lammps/src/
make yes-RIGID yes-MOLECULE
make mpi -j8
```

### Adding New Problems

1. Create `problems/problem-NNN-name/`
2. Write `problem.json` following schema
3. Provide `baseline_policy.json` and `README.md`
4. Set up `starter_template/`
5. Configure `benchmark/` task

See [FRAMEWORK_SUMMARY.md](FRAMEWORK_SUMMARY.md) for details.

---

## 📈 Project Status

### ✅ Phase 1: Core Infrastructure (Complete)
- JSON schema system
- Policy → C++ code generator
- Problem-001 (dimer_ksat) fully specified
- Participant documentation
- Baseline solutions

### 🚧 Phase 2: Automation (In Progress)
- [ ] Validation pipeline (`tools/validate_submission.py`)
- [ ] System generator (`problem.json` + `params.json` → LAMMPS inputs)
- [ ] Automated evaluation (`tools/run_evaluation.py`)
- [ ] Leaderboard aggregation

### 🔮 Phase 3: Expansion (Planned)
- [ ] Problem-002 (octahedron assembly)
- [ ] Problem-003 (Hamiltonian path)
- [ ] Problem-004 (frustrated sampling)
- [ ] Policy gradient optimization tools
- [ ] Community submissions & leaderboard

---

## 🤝 Contributing

### Submit a Solution
1. Copy starter template
2. Design your policy
3. Validate with `tools/validate_submission.py` (coming soon)
4. Submit via pull request or email

### Propose a Problem
1. Create `problem.json` following schema
2. Write problem README
3. Provide baseline and analysis scripts
4. Open pull request to `problems/`

### Improve Infrastructure
- Add trigger types to code generator
- Enhance validation checks
- Optimize LAMMPS performance
- Extend benchmark metrics

---

## 📜 Citation

If you use this framework in your research, please cite:

```bibtex
@software{molecular_computing_competition,
  title = {Molecular Computing Competition Framework},
  author = {Brenner Group},
  year = {2026},
  url = {https://github.com/Livia33g/lammps-state-change},
  note = {Kaggle-style benchmark for non-equilibrium molecular computation}
}
```

And the research proposal:
```bibtex
@article{brenner2026molecular,
  title = {Towards Design Principles for Molecular Computation},
  author = {Brenner, Michael P.},
  year = {2026},
  journal = {In preparation}
}
```

---

## 📬 Contact & Community

- **Issues**: [GitHub Issues](https://github.com/Livia33g/lammps-state-change/issues)
- **Discussions**: [GitHub Discussions](https://github.com/Livia33g/lammps-state-change/discussions)
- **Email**: [your-email@example.com]

---

## ⚖️ License

MIT License - See LICENSE file for details.

---

**Ready to compete?** Start with [Problem 001](problems/problem-001-dimer-ksat/README.md)! 🚀
