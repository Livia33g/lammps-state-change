# Repository Reorganization Plan: Kaggle-Style Policy Competition Framework

## Current Problems

1. **No standardized policy input format** - Each fix is hand-coded C++
2. **Documentation chaos** - 20+ scattered markdown files, many outdated
3. **Mixed purposes** - Framework code mixed with problem instances and submissions
4. **Manual workflow** - No automated pipeline from policy → results
5. **Unclear how to add new problems** or submit new policies

## Proposed Architecture

### 1. Input Format: Policy Specification (YAML)

**Key Innovation**: Policies are specified declaratively, not as C++ code.

```yaml
# submissions/octahedron/avoid_same_v1/policy.yaml
metadata:
  name: "avoid_same_type_v1"
  author: "livia"
  description: "Higher-ID monomer changes to different type when touching same type"
  version: "1.0"

rules:
  # Rule 1: Unreacted monomers (type 1) touching type 1
  - name: "type1_collision"
    trigger:
      condition: "same_type_contact"
      my_types: [1]          # I am type 1
      neighbor_types: [1]     # Touching type 1
      cutoff: 2.5
      confidence: 1000        # Hysteresis: wait this many steps

    action:
      selector: "higher_molecule_id"    # Which monomer changes
      change_type: "random_from_list"
      options: [3, 4, 5]                # Pick one uniformly
      probability: 1.0                  # Always change (if selected)

    cooldown:
      base_steps: 10000
      jitter_range: [0.3, 2.0]         # Multiply by random [0.3, 2.0]

  # Rule 2: Reacted monomers (types 3/4/5) touching same type
  - name: "same_type_conflict"
    trigger:
      condition: "same_type_contact"
      my_types: [3, 4, 5]
      neighbor_types: [3, 4, 5]  # Each checks its own type
      cutoff: 2.5
      confidence: 1000

    action:
      selector: "higher_molecule_id"
      change_type: "other_types_only"  # CRITICAL: exclude current type
      available_types: [3, 4, 5]
      probability: 1.0

    cooldown:
      base_steps: 10000
      jitter_range: [0.3, 2.0]

parameters:
  check_every: 100
  max_changes_per_step: 1
  global_rate_limit: true
```

**Alternative for advanced users**: Allow direct C++ code injection:

```yaml
metadata:
  name: "custom_advanced_policy"
  type: "cpp_template"

cpp_code:
  includes: |
    #include <custom_lib.h>

  change_logic: |
    // Direct C++ code for expert users
    if (my_type == 3 && coord_3 > 2.0) {
      new_type = 4;
    }
```

### 2. Problem Definition Format

```yaml
# problems/octahedron/problem.yaml
problem:
  name: "octahedron_assembly"
  description: "Assemble 6-monomer octahedrons with all-different-type rule"

  system:
    geometry: "octahedron"           # References geometry.py
    atom_types:
      - id: 1
        name: "core"
        mass: 1.0
      - id: 2
        name: "patch_initial"        # Effective type 1 for policy
        mass: 1.0
      - id: 3
        name: "patch_state_A"
        mass: 1.0
      - id: 4
        name: "patch_state_B"
        mass: 1.0
      - id: 5
        name: "patch_state_C"
        mass: 1.0

    patch_group: [2, 3, 4, 5]        # Which types are "patches" for state changes

    simulation:
      n_monomers: 100
      concentration: 0.001           # monomers/sigma^3
      temperature: 0.5               # kT in LJ units
      timestep: 0.005
      steps: 50000000
      thermo_every: 1000
      dump_every: 1000

    potentials:
      - style: "morse"
        cutoff: 2.5
        alpha: 5.0
        pairs:
          - types: [1, 1]
            style: "soft"
            params: {A: 500.0, rc: 2.0}
          - types: [2, 3]
            D0: 6.0
            r0: 0.0
          - types: [3, 3]
            D0: 12.0
            r0: 0.0

      - style: "exclude"
        condition: "molecule/intra all"

  target:
    structure: "cluster"
    size: 6                          # 6-monomer clusters
    composition: "all_different"     # Each monomer different type

    yield_definition:
      mode: "fraction_molecules"     # % of molecules in target clusters
      label_mode: "majority_site_type"
      site_types: [2, 3, 4, 5]

  scoring:
    metrics:
      - name: "final_yield"
        weight: 0.5
        maximize: true
      - name: "work_per_yield"
        weight: 0.3
        maximize: false              # Lower is better
      - name: "time_to_50pct_yield"
        weight: 0.2
        maximize: false

    score_formula: "final_yield * 1000 - work_per_yield * 0.1 - time_to_50pct / 1000"
```

### 3. Directory Structure (Clean Separation)

```
lammps-state-change/
│
├── README.md                      # Main entry: "Welcome to Policy Competition"
├── GETTING_STARTED.md             # Quick start guide
├── .gitignore                     # Ignore generated/ and results/
│
├── framework/                     # CORE INFRASTRUCTURE (stable, rarely changes)
│   ├── __init__.py
│   ├── generators/
│   │   ├── policy_parser.py      # Parse policy.yaml
│   │   ├── fix_codegen.py        # Generate C++ fix from policy
│   │   ├── lammps_generator.py   # Generate data.* and in.* files
│   │   └── templates/             # Jinja2 templates
│   │       ├── fix_base.cpp.j2
│   │       ├── fix_base.h.j2
│   │       └── lammps_input.j2
│   │
│   ├── analysis/
│   │   ├── trajectory.py          # Parse dump files
│   │   ├── thermo.py              # Parse thermo output
│   │   ├── work_calculator.py     # Calculate work from dE
│   │   ├── yield_calculator.py    # Calculate yield from clusters
│   │   └── scorer.py              # Compute final score
│   │
│   ├── benchmark/
│   │   ├── task_runner.py         # Orchestrate full pipeline
│   │   ├── leaderboard.py         # Aggregate scores
│   │   └── schemas/
│   │       ├── policy_schema.json
│   │       └── problem_schema.json
│   │
│   ├── build/
│   │   ├── lammps_builder.py      # Auto-rebuild LAMMPS
│   │   └── fix_installer.py       # Copy fix to LAMMPS src
│   │
│   └── validation/
│       ├── policy_validator.py    # Check policy.yaml is valid
│       └── problem_validator.py   # Check problem.yaml is valid
│
├── problems/                      # PROBLEM DEFINITIONS (one per challenge)
│   ├── octahedron/
│   │   ├── problem.yaml           # Problem specification
│   │   ├── geometry.py            # Octahedron monomer geometry
│   │   ├── README.md              # Problem description for users
│   │   └── baseline/              # Reference baseline policy
│   │       └── policy.yaml
│   │
│   ├── dimer/
│   │   ├── problem.yaml
│   │   ├── geometry.py
│   │   ├── README.md
│   │   └── baseline/
│   │       └── policy.yaml
│   │
│   └── dimer_ksat/
│       └── variants/
│           ├── 1core/
│           │   ├── problem.yaml
│           │   └── geometry.py
│           └── base/
│               └── problem.yaml
│
├── submissions/                   # USER SUBMISSIONS (gitignored except structure)
│   ├── octahedron/
│   │   ├── .gitkeep
│   │   ├── avoid_same_v1/         # Submission 1
│   │   │   ├── policy.yaml        # ONLY file user writes
│   │   │   ├── generated/         # Auto-generated (gitignored)
│   │   │   │   ├── fix_*.cpp
│   │   │   │   ├── fix_*.h
│   │   │   │   ├── data.*
│   │   │   │   └── in.*
│   │   │   ├── results/           # Auto-generated (gitignored)
│   │   │   │   ├── dump.*
│   │   │   │   ├── log.lammps
│   │   │   │   ├── slurm-*.err
│   │   │   │   ├── timeseries.csv
│   │   │   │   └── score.json
│   │   │   └── metadata.json      # Timestamp, status, score
│   │   │
│   │   └── avoid_same_v2/         # Submission 2 (iterate)
│   │       └── policy.yaml
│   │
│   └── dimer/
│       └── simple_flip/
│           └── policy.yaml
│
├── docs/                          # CONSOLIDATED DOCUMENTATION
│   ├── user_guide/
│   │   ├── getting_started.md
│   │   ├── writing_policies.md
│   │   ├── policy_format.md
│   │   ├── problem_format.md
│   │   └── scoring.md
│   │
│   ├── developer_guide/
│   │   ├── adding_problems.md
│   │   ├── codegen_internals.md
│   │   └── extending_framework.md
│   │
│   └── reference/
│       ├── api.md
│       ├── yaml_schemas.md
│       └── examples/
│           ├── simple_policy.yaml
│           └── advanced_policy.yaml
│
├── legacy/                        # ARCHIVED OLD CODE (don't delete, just move)
│   ├── manual_fixes/              # Old hand-coded fixes
│   │   ├── fix_state_change/
│   │   ├── dimer/fix_state_change_dimer.*
│   │   └── octahedron/fix_state_change_octahedron.*
│   │
│   ├── old_docs/                  # Outdated documentation
│   │   ├── DIAGNOSTIC_ISSUES.md
│   │   ├── MASS_BEHAVIOR_*.md
│   │   ├── CONSISTENCY_SWEEP_*.md
│   │   └── NaN_DIAGNOSIS.md
│   │
│   └── old_generators/
│       ├── generate_change.py
│       └── generate_change_cpp.py
│
├── scripts/                       # UTILITY SCRIPTS
│   ├── submit_policy.py           # Main entry point for users
│   ├── validate_policy.py         # Validate policy.yaml
│   ├── show_leaderboard.py        # Display rankings
│   └── cleanup_results.py         # Clean generated files
│
├── tests/                         # UNIT TESTS
│   ├── test_policy_parser.py
│   ├── test_fix_codegen.py
│   └── test_scorer.py
│
└── inverse-design/                # KEEP AS-IS (separate tool)
    └── ...
```

### 4. User Workflow (Dead Simple)

```bash
# Step 1: Clone repo
git clone https://github.com/Livia33g/lammps-state-change.git
cd lammps-state-change

# Step 2: Read problem
cat problems/octahedron/README.md

# Step 3: Write policy (ONLY thing user does)
mkdir -p submissions/octahedron/my_policy
vim submissions/octahedron/my_policy/policy.yaml

# Step 4: Submit and run
python scripts/submit_policy.py \
    --problem octahedron \
    --policy my_policy \
    --steps 50000000 \
    --partition gpu

# Behind the scenes (automated):
# 1. Validate policy.yaml against schema
# 2. Generate C++ fix from policy + templates
# 3. Copy to LAMMPS src, rebuild
# 4. Generate LAMMPS input files from problem.yaml
# 5. Submit SLURM job
# 6. Monitor completion
# 7. Analyze results (yield, work)
# 8. Compute score
# 9. Update leaderboard

# Step 5: Check results
python scripts/show_leaderboard.py --problem octahedron

# Output:
# Rank | Policy          | Author | Final Yield | Work/Yield | Score    | Date
# -----+-----------------+--------+-------------+------------+----------+----------
# 1    | avoid_same_v2   | livia  | 0.87        | 12.3       | 845.2    | 2026-01-19
# 2    | avoid_same_v1   | livia  | 0.82        | 15.1       | 802.4    | 2026-01-19
# 3    | baseline        | system | 0.45        | 8.2        | 442.8    | 2026-01-15
```

### 5. Documentation Consolidation Plan

**KEEP (move to docs/):**
- ✅ GETTING_STARTED.md → docs/user_guide/getting_started.md
- ✅ WORK_ANALYSIS_IMPROVEMENTS.md → docs/reference/work_analysis.md
- ✅ MULTIPLE_FIXES_GUIDE.md → docs/developer_guide/multiple_fixes.md
- ✅ REBUILD_INSTRUCTIONS.md → docs/developer_guide/building_lammps.md

**ARCHIVE (move to legacy/old_docs/):**
- 📦 octahedron/DIAGNOSTIC_ISSUES.md (debugging notes, now fixed)
- 📦 octahedron/MASS_BEHAVIOR_*.md (investigation notes)
- 📦 octahedron/CONSISTENCY_SWEEP_*.md (MPI fix notes)
- 📦 octahedron/NaN_DIAGNOSIS.md (bug hunt notes)
- 📦 octahedron/LOGIC_EVOLUTION_EXPLANATION.md (development history)
- 📦 octahedron/FIXES_APPLIED.md (changelog, consolidated)
- 📦 STATE_CHANGE_EXPLANATION.md (old approach, superseded)
- 📦 NEXT_STEPS_AFTER_REBUILD.md (temporary instructions)

**CONSOLIDATE (merge into new docs):**
- ✅ octahedron/README.md + octahedron/INSTALL_FIX.md → problems/octahedron/README.md
- ✅ octahedron/FIX_SAME_TYPE_CONVERGENCE_BUG.md → docs/reference/bugfixes.md
- ✅ ANALYSIS_WORK_YIELD_FROM_TRAJECTORY.md → docs/reference/analysis_tools.md

**DELETE (outdated by framework):**
- ❌ example_work_analysis_workflow.sh (replaced by submit_policy.py)
- ❌ Individual submit*.slurm files (generated automatically)

### 6. Code Generation Strategy

**Policy → C++ Fix (Automated)**

The `framework/generators/fix_codegen.py` will:

1. Parse `policy.yaml`
2. Validate against schema
3. Fill in Jinja2 template `fix_base.cpp.j2`:

```cpp
// Auto-generated from policy: {{ policy.metadata.name }}
// Author: {{ policy.metadata.author }}
// Generated: {{ timestamp }}

#include "fix_state_change_{{ problem_name }}_{{ policy_name }}.h"
// ... standard includes ...

void FixStateChange{{ ProblemName }}{{ PolicyName }}::check_and_change() {
    {% for rule in policy.rules %}
    // Rule: {{ rule.name }}
    if (my_type in {{ rule.trigger.my_types }}) {
        // Check condition: {{ rule.trigger.condition }}
        {% if rule.trigger.condition == "same_type_contact" %}
        double coord = get_coordination(i, my_type);
        if (coord >= 1.0 && contact_timer[i] > {{ rule.trigger.confidence }}) {
            {% if rule.action.selector == "higher_molecule_id" %}
            if (my_mol_id > all_neighbor_mol_ids) {
                {% if rule.action.change_type == "other_types_only" %}
                // EXCLUDE current type
                new_type = pick_random_excluding({{ rule.action.available_types }}, my_type);
                {% elif rule.action.change_type == "random_from_list" %}
                new_type = pick_random_from({{ rule.action.options }});
                {% endif %}
            }
            {% endif %}
        }
        {% endif %}
    }
    {% endfor %}
}
```

This allows:
- **Non-programmers** to write policies declaratively
- **Experts** to inject custom C++ if needed
- **Rapid iteration** without manual C++ editing

### 7. Implementation Phases

**Phase 1: Core Framework (Week 1-2)**
- ✅ Define YAML schemas (policy, problem)
- ✅ Implement policy_parser.py
- ✅ Implement fix_codegen.py with templates
- ✅ Implement submit_policy.py orchestrator
- ✅ Test on octahedron (hand-migrate current fix to YAML)

**Phase 2: Reorganize Repo (Week 2-3)**
- ✅ Create new directory structure
- ✅ Move code to framework/
- ✅ Move problems to problems/
- ✅ Archive legacy/ docs
- ✅ Write consolidated docs/
- ✅ Update .gitignore

**Phase 3: Analysis Pipeline (Week 3-4)**
- ✅ Integrate analyze_work_statechange_frames.py
- ✅ Implement scorer.py
- ✅ Implement leaderboard.py
- ✅ Add visualization (plots)

**Phase 4: Polish & Testing (Week 4-5)**
- ✅ Unit tests
- ✅ End-to-end test (submit baseline policy)
- ✅ Documentation review
- ✅ Example policies for each problem

**Phase 5: Multi-Problem Support (Week 5+)**
- ✅ Migrate dimer, dimer_ksat to new format
- ✅ Add 2-3 new challenge problems
- ✅ Beta test with external users

### 8. Key Design Decisions

**Q: YAML vs JSON vs TOML?**
**A:** YAML - best for human readability, supports comments

**Q: How to handle C++ code injection for advanced users?**
**A:** Allow `cpp_code:` section in policy.yaml, bypass template

**Q: How to version policies?**
**A:** Metadata includes version, Git tracks policy.yaml changes

**Q: How to prevent cheating (hardcoding specific problem)?**
**A:**
- Hidden test problems (same format, different parameters)
- Require policies to generalize (e.g., work for different n_monomers)

**Q: How to handle MPI/parallel?**
**A:** Framework auto-detects #procs from SLURM env, scales accordingly

**Q: What if policy.yaml can't express complex logic?**
**A:** Allow fallback to manual C++ submission (advanced track)

### 9. Benefits of This Refactor

**For Users (Policy Submitters):**
- ✅ Write policies in 50 lines of YAML vs 1000 lines of C++
- ✅ No C++ knowledge required
- ✅ Instant validation & feedback
- ✅ Automatic analysis & scoring
- ✅ Compare policies easily (leaderboard)

**For Problem Designers:**
- ✅ Add new problems without touching framework
- ✅ Standardized problem.yaml format
- ✅ Automatic LAMMPS input generation

**For Maintainers:**
- ✅ Clean separation: framework vs problems vs submissions
- ✅ Consolidated documentation
- ✅ Easier to test & debug
- ✅ Version control friendly

**For Kaggle-Style Competition:**
- ✅ Easy to submit (just upload policy.yaml)
- ✅ Reproducible (problem + policy fully specified)
- ✅ Scalable (add problems without code changes)
- ✅ Fair (hidden test problems prevent overfitting)

---

## Next Steps

I can help you implement any of these phases. Which would you like to start with?

1. **Define YAML schemas** (policy.yaml & problem.yaml format)
2. **Create directory structure** (reorganize current repo)
3. **Implement code generator** (policy → C++ fix)
4. **Write submit_policy.py** (orchestration script)
5. **Consolidate documentation** (move/archive files)

Or would you like me to start with a **proof-of-concept** where I migrate the octahedron fix to a YAML policy and demonstrate the full pipeline?
