# Revised Architecture: Multi-Level Parameter Control

## The Challenge: What's Fixed vs What's Tunable?

For a fair competition, we need to distinguish:
- **Problem constraints** (FIXED for everyone - defines the challenge)
- **System design** (TUNABLE - part of the solution space)
- **Policy** (TUNABLE - the main optimization target)

## Proposed Solution: Two-File Submission

### Submission Structure

```
submissions/octahedron/my_submission_v1/
├── policy.yaml           # State change rules (always required)
├── system.yaml           # System design choices (optional, uses defaults if missing)
└── metadata.json         # Auto-generated: timestamp, score, status
```

**Users submit BOTH files** (or just policy.yaml if using default system)

---

## File 1: policy.yaml (State Change Rules)

This remains as before - defines the state change logic:

```yaml
# submissions/octahedron/my_submission_v1/policy.yaml
metadata:
  name: "avoid_same_v1"
  author: "livia"

parameters:
  check_every: 100        # How often to check for state changes
  hysteresis_threshold: 1000  # Consecutive steps before trigger
  max_changes_per_step: 1

rules:
  - name: "same_type_conflict"
    trigger:
      when: "same_type_contact"
      my_types: [3, 4, 5]
      cutoff: 2.5         # Can reference system.yaml value
    action:
      change_type: "exclude_current_type"
    cooldown:
      base_steps: 10000
      jitter: [0.3, 2.0]
```

---

## File 2: system.yaml (System Design Choices)

This is NEW - specifies all tunable system parameters:

```yaml
# submissions/octahedron/my_submission_v1/system.yaml
metadata:
  name: "high_temp_weak_attraction"
  description: "Uses higher temperature and weaker patch attractions"

# GEOMETRY (if problem allows customization)
geometry:
  use_default: true             # Use problem-defined geometry
  # OR specify custom (if problem allows):
  # monomer_type: "octahedron"
  # custom_geometry:
  #   core_positions: [[0, 0, 0]]
  #   patch_positions: [[1, 0, 0], [0, 1, 0], ...]

# PHYSICAL PARAMETERS (tunable within problem constraints)
physical:
  temperature: 0.7              # kT (LJ units) - TUNABLE
  concentration: 0.001          # monomers/sigma^3 - TUNABLE
  n_monomers: 100               # Total system size - TUNABLE (within limits)

  # Constraints from problem (shown for reference, enforced by framework)
  # min_temperature: 0.1
  # max_temperature: 2.0
  # min_concentration: 0.0001
  # max_concentration: 0.01

# INTERACTION POTENTIALS (tunable within problem constraints)
potentials:
  morse:
    alpha: 5.0                  # Morse width parameter - TUNABLE
    cutoff: 2.5                 # Morse cutoff - TUNABLE

    # Pair-specific well depths (TUNABLE)
    pairs:
      # Type 2-3 interaction (Red-Blue in dimer, Initial-A in octahedron)
      - types: [2, 3]
        D0: 8.0                 # Well depth (default: 6.0) - TUNED UP
        r0: 0.0                 # Equilibrium distance

      # Type 3-3 interaction (Blue-Blue, A-A)
      - types: [3, 3]
        D0: 15.0                # Well depth (default: 12.0) - TUNED UP
        r0: 0.0

      # Can add more pairs if needed
      - types: [4, 4]
        D0: 15.0
        r0: 0.0

      - types: [5, 5]
        D0: 15.0
        r0: 0.0

      # Cross-interactions between different types
      - types: [3, 4]
        D0: 10.0                # Inter-type attraction
        r0: 0.0

      - types: [3, 5]
        D0: 10.0
        r0: 0.0

      - types: [4, 5]
        D0: 10.0
        r0: 0.0

  soft:
    # Core-core repulsion
    - types: [1, 1]
      A: 500.0                  # Repulsion strength - TUNABLE
      rc: 2.0                   # Cutoff

# SIMULATION PARAMETERS
simulation:
  timestep: 0.005               # MD timestep (LJ units) - TUNABLE (with care!)
  run_steps: 50000000           # Total steps - TUNABLE

  # Integration method
  integrator:
    type: "rigid/nve"           # Fixed for rigid bodies (not tunable)
    molecule_based: true

  # Thermostat
  thermostat:
    type: "langevin"            # Type: langevin, nose-hoover, etc. - TUNABLE
    temperature: 0.7            # References physical.temperature
    damping: 0.5                # Damping parameter - TUNABLE
    seed: 12345                 # Random seed

  # Neighbor list settings
  neighbor:
    skin: 1.0                   # Neighbor list skin - TUNABLE
    every: 1                    # Rebuild frequency
    delay: 0
    check: true

# OUTPUT CONTROL (tunable but affects computational cost)
output:
  thermo_every: 1000            # Thermo output frequency - TUNABLE
  dump_every: 1000              # Dump trajectory frequency - TUNABLE

  # What to include in dump (affects file size)
  dump_columns: ["id", "mol", "type", "x", "y", "z"]
  # Could add: "vx", "vy", "vz", "fx", "fy", "fz" if needed

  # Thermo style
  thermo_style: ["step", "temp", "pe", "ke", "etotal"]

# ADVANCED: Resource constraints (auto-calculated if not specified)
resources:
  walltime: "24:00:00"          # SLURM walltime request
  partition: "gpu"              # SLURM partition
  nodes: 1
  ntasks_per_node: 4            # MPI tasks
```

---

## File 3: problem.yaml (Defines Challenge & Constraints)

This is set by the problem designer and defines what's FIXED vs TUNABLE:

```yaml
# problems/octahedron/problem.yaml
problem:
  name: "octahedron_assembly"
  description: "Assemble 6-monomer octahedrons with all-different-type rule"

# FIXED: Target structure (defines the challenge)
target:
  structure: "cluster"
  size: 6
  composition: "all_different"

  yield_definition:
    mode: "fraction_molecules"
    label_mode: "majority_site_type"
    site_types: [2, 3, 4, 5]

# FIXED: Geometry (competitors cannot change monomer shape)
geometry:
  type: "octahedron"
  customizable: false           # If true, users can modify geometry

  default:
    monomer_file: "octahedron_geometry.py"
    atoms_per_monomer: 30
    core_atoms: 1
    patch_atoms: 24
    atom_types:
      - id: 1
        name: "core"
        group: "cores"
      - id: 2
        name: "patch_initial"
        effective_type: 1       # Maps to policy's "type 1"
        group: "patches"
      - id: 3
        name: "patch_A"
        group: "patches"
      - id: 4
        name: "patch_B"
        group: "patches"
      - id: 5
        name: "patch_C"
        group: "patches"

# TUNABLE: Physical parameters (with constraints)
physical_constraints:
  temperature:
    tunable: true
    min: 0.1
    max: 2.0
    default: 0.5
    warning: "High temps (>1.0) may prevent assembly; low temps (<0.3) may cause freezing"

  concentration:
    tunable: true
    min: 0.0001
    max: 0.01
    default: 0.001
    warning: "High concentration increases collision rate but may cause jamming"

  n_monomers:
    tunable: true
    min: 20
    max: 500
    default: 100
    note: "Larger systems are more statistically robust but computationally expensive"

# TUNABLE: Interaction potentials (with constraints)
potential_constraints:
  morse:
    alpha:
      tunable: true
      min: 1.0
      max: 10.0
      default: 5.0
      warning: "Low alpha (<3) gives soft potentials; high alpha (>7) gives stiff potentials"

    cutoff:
      tunable: true
      min: 1.5
      max: 4.0
      default: 2.5
      warning: "Must be large enough to capture attractive well"

    D0_range:
      tunable: true
      min: 0.0
      max: 50.0
      default: 6.0
      note: "Well depth controls binding strength"

    r0:
      tunable: false            # Fixed at 0.0 for patch overlap
      value: 0.0

# TUNABLE: Simulation parameters (with stability constraints)
simulation_constraints:
  timestep:
    tunable: true
    min: 0.001
    max: 0.01
    default: 0.005
    warning: "Smaller timestep (<0.003) is more stable but slower; larger (>0.007) may be unstable"

  run_steps:
    tunable: true
    min: 1000000
    max: 100000000
    default: 50000000
    note: "Longer runs allow more equilibration but cost more compute time"

  thermostat_damping:
    tunable: true
    min: 0.1
    max: 2.0
    default: 0.5
    note: "Lower damping (<0.3) couples weakly to thermostat; higher (>1.0) may overdamp"

# TUNABLE: Output control
output_constraints:
  thermo_every:
    tunable: true
    min: 100
    max: 10000
    default: 1000
    note: "More frequent output helps analysis but increases file size"

  dump_every:
    tunable: true
    min: 100
    max: 10000
    default: 1000
    note: "Trajectory files can become very large with frequent dumps"

# SCORING: How submissions are evaluated
scoring:
  metrics:
    - name: "final_yield"
      weight: 0.5
      maximize: true
    - name: "work_per_yield"
      weight: 0.3
      maximize: false
    - name: "time_to_50pct_yield"
      weight: 0.2
      maximize: false

  score_formula: "final_yield * 1000 - work_per_yield * 0.1 - time_to_50pct / 1000"

  # Penalize violations (e.g., instability, NaNs)
  penalties:
    nan_energy: -1000
    lost_atoms: -1000
    incomplete_run: -500
```

---

## Validation Hierarchy

When a user submits `policy.yaml` + `system.yaml`:

```python
# Framework validation logic
def validate_submission(policy_yaml, system_yaml, problem_yaml):
    # Step 1: Validate policy.yaml against schema
    validate_policy_schema(policy_yaml)

    # Step 2: Validate system.yaml against schema
    validate_system_schema(system_yaml)

    # Step 3: Check system.yaml against problem constraints
    check_physical_constraints(system_yaml.physical, problem_yaml.physical_constraints)
    check_potential_constraints(system_yaml.potentials, problem_yaml.potential_constraints)
    check_simulation_constraints(system_yaml.simulation, problem_yaml.simulation_constraints)

    # Step 4: Check policy.yaml references valid system parameters
    check_policy_system_consistency(policy_yaml, system_yaml)

    # Step 5: Warn about risky choices
    issue_warnings(system_yaml, problem_yaml)
```

### Example Validation Output

```bash
$ python scripts/validate_submission.py \
    --problem octahedron \
    --policy submissions/octahedron/my_v1/policy.yaml \
    --system submissions/octahedron/my_v1/system.yaml

✓ Policy schema valid
✓ System schema valid

⚠ WARNING: system.yaml
  Line 12: physical.temperature = 1.2

  You set temperature to 1.2, which is above the recommended range (0.3-0.8).
  High temperatures may prevent cluster formation.

  Hint: Try temperature in [0.4, 0.7] for typical assembly problems.

⚠ WARNING: system.yaml
  Line 28: potentials.morse.pairs[0].D0 = 3.0

  Well depth D0=3.0 is quite weak (default: 6.0).
  Weak attractions may result in low yield.

✓ All constraints satisfied
✓ Policy-system consistency ok

VALIDATION PASSED (with 2 warnings)

Ready to submit? (y/n)
```

---

## Default Behavior: system.yaml is Optional

If user only submits `policy.yaml`, use problem defaults:

```bash
# Only policy.yaml provided
submissions/octahedron/simple_submission/
└── policy.yaml

# Framework auto-generates:
submissions/octahedron/simple_submission/
├── policy.yaml
└── system.yaml          # Auto-generated from problem.yaml defaults
```

---

## User Workflow Examples

### Beginner: Just Change Policy

```bash
# User only writes policy.yaml (uses all defaults)
mkdir -p submissions/octahedron/my_first_policy
vim submissions/octahedron/my_first_policy/policy.yaml

python scripts/submit_policy.py \
    --problem octahedron \
    --policy my_first_policy
# Uses default temperature, potentials, etc.
```

### Intermediate: Tune System Parameters

```bash
# User writes both files
mkdir -p submissions/octahedron/high_temp_policy
vim submissions/octahedron/high_temp_policy/policy.yaml
vim submissions/octahedron/high_temp_policy/system.yaml  # Override temperature, D0

python scripts/submit_policy.py \
    --problem octahedron \
    --policy high_temp_policy \
    --system high_temp_policy
# Uses custom temperature and potentials
```

### Advanced: Full Control

```bash
# User specifies everything
mkdir -p submissions/octahedron/optimized_full
vim submissions/octahedron/optimized_full/policy.yaml
vim submissions/octahedron/optimized_full/system.yaml

# system.yaml includes:
# - Custom temperature profile
# - Optimized D0 values for each pair
# - Fine-tuned timestep
# - Custom thermostat damping
# - Optimized output frequencies

python scripts/submit_policy.py \
    --problem octahedron \
    --policy optimized_full \
    --system optimized_full \
    --validate-strict  # Extra validation for advanced users
```

---

## Leaderboard Shows System Choices

```
Rank | Policy           | Temp | D0(2-3) | D0(3-3) | Final Yield | Work/Yield | Score
-----+------------------+------+---------+---------+-------------+------------+-------
1    | optimized_full   | 0.65 | 7.5     | 14.0    | 0.89        | 11.2       | 867.8
2    | high_temp_policy | 0.85 | 6.0     | 12.0    | 0.78        | 9.8        | 770.2
3    | my_first_policy  | 0.50 | 6.0     | 12.0    | 0.71        | 15.3       | 694.7
```

**Hovering over submission shows full system.yaml** for reproducibility.

---

## Hidden Test Problems: Generalization Check

To prevent overfitting to specific parameter choices:

```yaml
# problems/octahedron/problem.yaml (public)
scoring:
  public_test: true              # Users can see this problem

  hidden_tests:
    enabled: true
    description: "Your submission will also be tested on hidden variants"

    variants:
      - name: "octahedron_high_temp"
        override:
          physical.temperature: 0.9
          physical.concentration: 0.002

      - name: "octahedron_low_attraction"
        override:
          potentials.morse.D0_scale: 0.5    # All D0 values multiplied by 0.5

      - name: "octahedron_large_system"
        override:
          physical.n_monomers: 300

    scoring:
      public_weight: 0.6        # 60% of score from public test
      hidden_weight: 0.4        # 40% from average of hidden tests
```

This incentivizes policies that generalize well rather than overfit to one parameter set.

---

## Implementation: How Framework Handles This

```python
# framework/generators/lammps_generator.py
def generate_lammps_input(problem_yaml, policy_yaml, system_yaml):
    # Merge system.yaml with problem defaults
    system_params = merge_with_defaults(system_yaml, problem_yaml['default_system'])

    # Validate against constraints
    validate_constraints(system_params, problem_yaml['constraints'])

    # Generate LAMMPS input script
    input_script = f"""
# Auto-generated from policy: {policy_yaml['metadata']['name']}
# System: {system_yaml['metadata']['name']}

units           lj
atom_style      full
boundary        p p p

read_data       data.{problem_yaml['name']}

# Temperature: {system_params['physical']['temperature']} (default: {problem_yaml['physical_constraints']['temperature']['default']})
variable T equal {system_params['physical']['temperature']}

# Potentials (from system.yaml)
pair_style hybrid morse {system_params['potentials']['morse']['cutoff']} soft 2.0

# Pair coefficients (user-tuned)
"""

    for pair in system_params['potentials']['morse']['pairs']:
        types = pair['types']
        D0 = pair['D0']
        alpha = system_params['potentials']['morse']['alpha']
        r0 = pair['r0']
        input_script += f"pair_coeff {types[0]} {types[1]} morse {D0} {alpha} {r0}\n"

    input_script += f"""
# Integration (from system.yaml)
fix fx_nve all rigid/nve molecule
fix fx_langevin all langevin $T $T {system_params['simulation']['thermostat']['damping']} {system_params['simulation']['thermostat']['seed']}

# State change policy (from policy.yaml)
fix sc patches state/change/{problem_yaml['name']}/{policy_yaml['metadata']['name']} \\
    {policy_yaml['parameters']['check_every']} \\
    {policy_yaml['parameters']['hysteresis_threshold']} \\
    ...

# Output (from system.yaml)
thermo_style custom {' '.join(system_params['output']['thermo_style'])}
thermo {system_params['output']['thermo_every']}
dump d1 all custom {system_params['output']['dump_every']} dump.lammpstrj {' '.join(system_params['output']['dump_columns'])}

# Run (from system.yaml)
timestep {system_params['simulation']['timestep']}
run {system_params['simulation']['run_steps']}
"""

    return input_script
```

---

## Summary: What's Tunable at Each Level

| Parameter | Controlled By | Tunable? | Validated Against |
|-----------|---------------|----------|-------------------|
| Target structure | `problem.yaml` | ❌ FIXED | N/A (defines challenge) |
| Monomer geometry | `problem.yaml` | ✅/❌ (per problem) | `problem.yaml:geometry.customizable` |
| Temperature | `system.yaml` | ✅ YES | `problem.yaml:physical_constraints.temperature` |
| Concentration | `system.yaml` | ✅ YES | `problem.yaml:physical_constraints.concentration` |
| Potential D0 | `system.yaml` | ✅ YES | `problem.yaml:potential_constraints.D0_range` |
| Potential alpha | `system.yaml` | ✅ YES | `problem.yaml:potential_constraints.alpha` |
| Timestep | `system.yaml` | ✅ YES | `problem.yaml:simulation_constraints.timestep` |
| Check frequency | `policy.yaml` | ✅ YES | No constraint (part of policy) |
| Hysteresis | `policy.yaml` | ✅ YES | No constraint (part of policy) |
| Change rules | `policy.yaml` | ✅ YES | Schema validation only |
| Dump frequency | `system.yaml` | ✅ YES | `problem.yaml:output_constraints.dump_every` |

**Key Insight**: The problem defines the CHALLENGE (what to assemble), and competitors design both the SYSTEM (how monomers interact) and the POLICY (how they adapt).

---

This architecture gives users full control while maintaining fairness through:
1. **Constraints** prevent unrealistic choices
2. **Warnings** guide users toward reasonable ranges
3. **Hidden tests** prevent overfitting
4. **Leaderboard transparency** shows what system choices led to success
