# Complete Submission Example: All Three Files

This document shows a complete, concrete example of how problem.yaml, system.yaml, and policy.yaml work together for an octahedron assembly challenge.

---

## Example Scenario

**Problem**: Assemble octahedrons where each monomer connects only to different-type neighbors

**User**: "livia" wants to try a high-temperature, strong-attraction strategy with aggressive state-change policy

**Submission files**:
1. `problems/octahedron/problem.yaml` (provided by competition)
2. `submissions/octahedron/aggressive_high_temp_v1/system.yaml` (user's system design)
3. `submissions/octahedron/aggressive_high_temp_v1/policy.yaml` (user's state-change rules)

---

## File 1: problem.yaml (Provided by Competition)

```yaml
# problems/octahedron/problem.yaml
# This file is PROVIDED - users cannot modify it

problem:
  id: "octahedron_v1"
  name: "octahedron_assembly"
  description: |
    Assemble patchy monomers into 6-monomer octahedrons.
    Rule: Each monomer must be connected ONLY to different-type neighbors.
    Monomers can be type A, B, or C (3, 4, or 5 in LAMMPS).

  version: "1.0"
  created: "2026-01-15"

# ========== FIXED: Cannot be changed by competitors ==========

target:
  structure: "cluster"
  size: 6                       # 6-monomer octahedrons
  composition: "all_different"  # Each monomer must be different type from neighbors

  yield_definition:
    mode: "fraction_molecules"
    label_mode: "majority_site_type"
    site_types: [2, 3, 4, 5]    # Patches that determine molecule type
    bond_cutoff: 0.7            # Distance to consider molecules "bonded"

geometry:
  type: "octahedron"
  customizable: false           # ❌ Users CANNOT change monomer shape
  file: "octahedron_geometry.py"

  structure:
    atoms_per_monomer: 30
    composition:
      - type: 1                 # Core
        count: 1
        mass: 1.0
        group: "cores"
      - type: 2                 # Initial patch (maps to policy "type 1")
        count: 24
        mass: 1.0
        group: "patches"
        effective_type: 1       # For policy logic
      - type: 3                 # Patch state A
        count: 0                # Initially none
        mass: 1.0
        group: "patches"
      - type: 4                 # Patch state B
        count: 0
        mass: 1.0
        group: "patches"
      - type: 5                 # Patch state C
        count: 0
        mass: 1.0
        group: "patches"

scoring:
  metrics:
    - name: "final_yield"
      description: "Fraction of molecules in target octahedrons at end"
      weight: 0.5
      maximize: true

    - name: "work_per_yield"
      description: "Total work divided by final yield (efficiency)"
      weight: 0.3
      maximize: false

    - name: "time_to_50pct_yield"
      description: "Timesteps to reach 50% yield"
      weight: 0.2
      maximize: false

  formula: |
    score = final_yield * 1000
            - work_per_yield * 0.1
            - time_to_50pct / 1000

  penalties:
    nan_energy: -1000
    lost_atoms: -1000
    incomplete_run: -500

# ========== TUNABLE: Competitors can change within constraints ==========

constraints:
  physical:
    temperature:
      min: 0.1
      max: 2.0
      default: 0.5
      recommended: [0.4, 0.8]
      units: "kT in LJ"
      warning: >
        Temperature >1.0 may prevent assembly due to thermal disruption.
        Temperature <0.3 may cause kinetic trapping.

    concentration:
      min: 0.0001
      max: 0.01
      default: 0.001
      recommended: [0.0005, 0.002]
      units: "monomers per sigma^3"
      note: >
        Higher concentration increases collision rate but may cause jamming.
        Lower concentration slows assembly but improves quality.

    n_monomers:
      min: 20
      max: 500
      default: 100
      recommended: [50, 200]
      note: >
        Larger systems give better statistics but cost more compute time.
        Minimum 20 monomers to form 3-4 octahedrons.

  potentials:
    morse:
      alpha:
        min: 1.0
        max: 10.0
        default: 5.0
        recommended: [4.0, 6.0]
        note: "Controls steepness of Morse potential"

      cutoff:
        min: 1.5
        max: 4.0
        default: 2.5
        recommended: [2.0, 3.0]
        note: "Must be large enough to capture attractive well"

      D0_per_pair:
        min: 0.0
        max: 50.0
        default: 6.0
        recommended: [4.0, 15.0]
        note: >
          Well depth. Higher = stronger binding.
          Different pairs can have different values.

      r0:
        fixed: 0.0              # ❌ Cannot change (patch overlap geometry)

    soft:
      core_repulsion:
        min: 100.0
        max: 1000.0
        default: 500.0

  simulation:
    timestep:
      min: 0.001
      max: 0.01
      default: 0.005
      recommended: [0.003, 0.007]
      warning: >
        Timestep >0.007 may cause instability.
        Timestep <0.003 makes simulation very slow.

    run_steps:
      min: 1000000
      max: 100000000
      default: 50000000
      note: "Longer runs allow more equilibration"

    thermostat:
      type:
        options: ["langevin", "nose-hoover"]
        default: "langevin"

      damping:
        min: 0.1
        max: 2.0
        default: 0.5
        recommended: [0.3, 0.8]

  output:
    thermo_every:
      min: 100
      max: 10000
      default: 1000

    dump_every:
      min: 100
      max: 10000
      default: 1000

  policy:
    check_every:
      min: 10
      max: 10000
      default: 100
      note: "How often to check for state changes (timesteps)"

    hysteresis_threshold:
      min: 0
      max: 10000
      default: 1000
      note: "Consecutive steps of contact before triggering"

    cooldown_base:
      min: 0
      max: 100000
      default: 10000
      note: "Minimum steps before molecule can change again"

hidden_tests:
  enabled: true
  description: "Final score includes performance on hidden parameter variations"

  variants:
    - name: "high_temp_variant"
      description: "Tests robustness to higher temperature"
      override:
        physical.temperature: 0.9

    - name: "low_attraction_variant"
      description: "Tests with weaker attractions"
      override:
        potentials.morse.D0_scale: 0.6

    - name: "high_density_variant"
      description: "Tests with higher monomer concentration"
      override:
        physical.concentration: 0.005
        physical.n_monomers: 200

  scoring:
    public_weight: 0.6
    hidden_weight: 0.4          # Average of 3 hidden tests
```

---

## File 2: system.yaml (User's System Design)

```yaml
# submissions/octahedron/aggressive_high_temp_v1/system.yaml
# This is WRITTEN BY THE USER - their system design choices

metadata:
  name: "aggressive_high_temp_v1"
  author: "livia"
  description: |
    High-temperature strategy with very strong attractions.
    Hypothesis: Fast dynamics + strong binding = quick assembly.
  version: "1.0"
  created: "2026-01-19"

# USER'S CHOICES: Physical parameters
physical:
  temperature: 0.75             # ⬆️ HIGHER than default (0.5)
  concentration: 0.0015         # ⬆️ HIGHER than default (0.001)
  n_monomers: 150               # ⬆️ LARGER than default (100)

  reasoning: |
    Higher temperature increases collision rate and helps escape local minima.
    Higher concentration increases encounter frequency.
    Larger system improves statistics.

# USER'S CHOICES: Interaction potentials
potentials:
  morse:
    alpha: 5.5                  # Slightly steeper than default
    cutoff: 2.5                 # Keep default

    # Pair-specific well depths (key tuning!)
    pairs:
      # Type 2-3: Initial patches attracting to state A
      - types: [2, 3]
        D0: 8.0                 # ⬆️ STRONGER than default (6.0)
        r0: 0.0

      # Same-type attractions (for conflict detection)
      - types: [3, 3]           # A-A
        D0: 16.0                # ⬆️ VERY STRONG (default: 12.0)
        r0: 0.0

      - types: [4, 4]           # B-B
        D0: 16.0
        r0: 0.0

      - types: [5, 5]           # C-C
        D0: 16.0
        r0: 0.0

      # Cross-type attractions (for octahedron assembly)
      - types: [3, 4]           # A-B
        D0: 12.0                # ⬆️ STRONG (reward different-type contacts)
        r0: 0.0

      - types: [3, 5]           # A-C
        D0: 12.0
        r0: 0.0

      - types: [4, 5]           # B-C
        D0: 12.0
        r0: 0.0

    reasoning: |
      Strong same-type attraction (16.0) ensures conflicts are detected.
      Strong cross-type attraction (12.0) stabilizes octahedrons.
      Difference (16 vs 12) provides thermodynamic drive to avoid same-type.

  soft:
    - types: [1, 1]
      A: 600.0                  # ⬆️ STRONGER core repulsion (default: 500)
      rc: 2.0

    reasoning: |
      Stronger core repulsion prevents overlapping at high temperature.

# USER'S CHOICES: Simulation parameters
simulation:
  timestep: 0.004               # ⬇️ SMALLER (more stable at high T)
  run_steps: 75000000           # ⬆️ LONGER (default: 50M)

  integrator:
    type: "rigid/nve"
    molecule_based: true

  thermostat:
    type: "langevin"
    temperature: 0.75           # Matches physical.temperature
    damping: 0.4                # ⬇️ WEAKER damping (faster dynamics)
    seed: 42                    # Different seed for this run

  neighbor:
    skin: 1.2                   # ⬆️ LARGER (fewer rebuilds)
    every: 1
    delay: 0
    check: true

  reasoning: |
    Smaller timestep (0.004) maintains stability at high temperature.
    Longer run (75M steps) allows system to fully equilibrate.
    Weaker damping (0.4) couples less to thermostat, faster dynamics.

# USER'S CHOICES: Output control
output:
  thermo_every: 1000            # Keep default
  dump_every: 5000              # ⬆️ LESS FREQUENT (save disk space)

  dump_columns: ["id", "mol", "type", "x", "y", "z"]
  thermo_style: ["step", "temp", "pe", "ke", "etotal"]

  reasoning: |
    Less frequent dumps (5000) reduces file size for long run.
    Still have enough resolution for analysis.

# Resource allocation
resources:
  walltime: "48:00:00"          # 48 hours for 75M steps
  partition: "gpu"
  nodes: 1
  ntasks_per_node: 4
```

---

## File 3: policy.yaml (User's State-Change Rules)

```yaml
# submissions/octahedron/aggressive_high_temp_v1/policy.yaml
# This is WRITTEN BY THE USER - their state-change policy

metadata:
  name: "aggressive_flip_v1"
  author: "livia"
  description: |
    Aggressive state-change policy with short hysteresis.
    Rapid flipping to explore configuration space quickly.
  version: "1.0"
  created: "2026-01-19"

# Policy-level parameters
parameters:
  check_every: 50               # ⬇️ MORE FREQUENT checks (default: 100)
  max_changes_per_step: 2       # ⬆️ ALLOW MORE changes per step (default: 1)

  reasoning: |
    Frequent checks (every 50 steps) catch conflicts quickly.
    Allow 2 changes per step to speed up equilibration.

# State-change rules
rules:
  # Rule 1: Unreacted monomers (type 1) colliding
  - name: "type1_fast_reaction"
    description: "Type 1 quickly evolves to A/B/C"

    trigger:
      when: "same_type_contact"
      my_types: [1]             # Effective type 1 (LAMMPS type 2)
      neighbor_types: [1]
      cutoff: 2.5               # From system.yaml

      hysteresis:
        enabled: true
        threshold: 500          # ⬇️ SHORT hysteresis (default: 1000)
        mode: "consecutive"

    action:
      selector:
        method: "higher_molecule_id"

      change_type:
        method: "random_from_list"
        options: [3, 4, 5]      # Equal chance → A, B, or C

      probability: 1.0

    cooldown:
      base_steps: 5000          # ⬇️ SHORT cooldown (default: 10000)
      jitter:
        enabled: true
        range: [0.5, 1.5]       # ⬇️ LESS jitter (default: [0.3, 2.0])
        per_molecule: true

    reasoning: |
      Short hysteresis (500) and cooldown (5000) allow rapid exploration.
      At high temperature, monomers move fast, so we can be aggressive.

  # Rule 2: Reacted monomers (A/B/C) touching same type
  - name: "same_type_aggressive_flip"
    description: "Quickly resolve conflicts with short hysteresis"

    trigger:
      when: "same_type_contact"
      my_types: [3, 4, 5]
      neighbor_types: "same_as_mine"
      cutoff: 2.5

      hysteresis:
        enabled: true
        threshold: 300          # ⬇️⬇️ VERY SHORT (default: 1000)
        mode: "consecutive"

    action:
      selector:
        method: "higher_molecule_id"

      change_type:
        method: "exclude_current_type"  # MUST change to different type
        available_types: [3, 4, 5]

      probability: 1.0

    cooldown:
      base_steps: 3000          # ⬇️⬇️ VERY SHORT (default: 10000)
      jitter:
        enabled: true
        range: [0.5, 1.5]
        per_molecule: true

    reasoning: |
      Very aggressive: hysteresis=300, cooldown=3000.
      Strategy: At high T, conflicts are frequent, so resolve them quickly.
      Risk: May cause oscillations, but strong attractions should stabilize.

# Synchronization (MPI settings)
synchronization:
  forward_comm: true
  consistency_sweep:
    enabled: true
    timestamp_threshold: 100
    respect_cooldown: true

# Summary of strategy
strategy_summary: |
  This submission uses an "aggressive high-temperature" approach:

  SYSTEM DESIGN:
  - High temperature (0.75) for fast dynamics
  - Strong attractions (D0=16 same-type, D0=12 cross-type)
  - High concentration (0.0015) for frequent encounters
  - Small timestep (0.004) for stability

  POLICY:
  - Short hysteresis (300-500) catches conflicts quickly
  - Short cooldowns (3000-5000) allow rapid re-attempts
  - Frequent checks (every 50 steps)
  - Allow 2 changes per step

  HYPOTHESIS:
  High temperature overcomes barriers, strong attractions lock in correct
  structures, and aggressive policy explores space quickly. May have higher
  work_per_yield but faster time_to_yield.

  RISKS:
  - High temp may prevent formation (entropy wins over enthalpy)
  - Aggressive policy may cause thrashing
  - Short cooldowns may not allow structures to stabilize

  Expected tradeoff: Fast assembly, but possibly lower final yield.
```

---

## Validation Output

```bash
$ python scripts/validate_submission.py \
    --problem octahedron \
    --policy aggressive_high_temp_v1 \
    --system aggressive_high_temp_v1

Validating submission: aggressive_high_temp_v1
Problem: octahedron_v1

✓ Policy schema valid
✓ System schema valid

Checking constraints...

✓ physical.temperature = 0.75 [OK]
  (within range [0.1, 2.0], but outside recommended [0.4, 0.8])

⚠ WARNING: system.yaml:12
  physical.temperature = 0.75

  You set temperature to 0.75, which is above the recommended range (0.4-0.8).
  Temperature >1.0 may prevent assembly due to thermal disruption.

  Your reasoning: "Higher temperature increases collision rate..."
  This is a valid strategy, but be aware of the risk.

✓ physical.concentration = 0.0015 [OK]
✓ physical.n_monomers = 150 [OK]

✓ potentials.morse.alpha = 5.5 [OK]
✓ potentials.morse.cutoff = 2.5 [OK]
✓ potentials.morse.pairs[*].D0 = [8.0, 16.0, 16.0, ...] [OK]
  (all within range [0.0, 50.0])

✓ simulation.timestep = 0.004 [OK]
✓ simulation.run_steps = 75000000 [OK]
✓ simulation.thermostat.damping = 0.4 [OK]

✓ output.dump_every = 5000 [OK]

✓ policy.parameters.check_every = 50 [OK]

⚠ WARNING: policy.yaml:14
  parameters.max_changes_per_step = 2

  You allow 2 molecules to change per timestep (default: 1).
  This may speed up equilibration but increases risk of cascading changes.

✓ policy.rules[0].trigger.hysteresis.threshold = 500 [OK]
  (within range, but lower than default)

⚠ INFO: policy.yaml:35
  rules[0].cooldown.base_steps = 5000

  Cooldown is shorter than default (10000).
  This allows faster re-attempts but may cause oscillations.

⚠ INFO: policy.yaml:68
  rules[1].trigger.hysteresis.threshold = 300

  Very short hysteresis (300 vs default 1000).
  Aggressive strategy - ensure your system can handle rapid changes.

⚠ INFO: policy.yaml:80
  rules[1].cooldown.base_steps = 3000

  Very short cooldown (3000 vs default 10000).
  This is quite aggressive - may cause thrashing.

✓ All constraints satisfied
✓ Policy-system consistency OK

========================================
VALIDATION SUMMARY
========================================
Status: PASSED (with 5 warnings)

Your submission is valid and can be run.

Warnings indicate deviations from defaults/recommendations.
Review your reasoning for each choice.

Estimated cost:
  - Run time: ~48 hours (75M steps, timestep=0.004)
  - Output size: ~500 MB (dump every 5000)

Ready to submit? (y/n)
```

---

## Submission & Execution

```bash
$ python scripts/submit_policy.py \
    --problem octahedron \
    --policy aggressive_high_temp_v1 \
    --system aggressive_high_temp_v1 \
    --validate

Validating... ✓ PASSED (with 5 warnings)

Generating LAMMPS files...
  ✓ Generated: aggressive_high_temp_v1/generated/data.octahedron
  ✓ Generated: aggressive_high_temp_v1/generated/in.octahedron
  ✓ Generated: aggressive_high_temp_v1/generated/fix_state_change_octahedron_aggressive_flip_v1.cpp
  ✓ Generated: aggressive_high_temp_v1/generated/fix_state_change_octahedron_aggressive_flip_v1.h

Compiling LAMMPS...
  ✓ Copied fix to LAMMPS src
  ✓ Rebuilt LAMMPS with new fix
  ✓ Verified fix registration

Submitting SLURM job...
  ✓ Job submitted: 12345678
  ✓ Partition: gpu
  ✓ Walltime: 48:00:00
  ✓ Output: aggressive_high_temp_v1/results/

Monitor job with:
  squeue -j 12345678
  tail -f aggressive_high_temp_v1/results/slurm-12345678.out

View live analysis:
  python scripts/live_analysis.py --job 12345678
```

---

## Leaderboard After Completion

```
OCTAHEDRON ASSEMBLY - PUBLIC LEADERBOARD
Last updated: 2026-01-20 14:32:00

Rank | Submission             | Temp | D0_same | Final | Work/ | Time   | Score  | Status
     |                        |      |         | Yield | Yield | to 50% |        |
-----+------------------------+------+---------+-------+-------+--------+--------+---------
1    | optimized_balanced_v2  | 0.55 | 12.0    | 0.87  | 11.5  | 8.2M   | 854.3  | ✓ DONE
2    | aggressive_high_temp_v1| 0.75 | 16.0    | 0.81  | 15.2  | 3.5M   | 792.3  | ✓ DONE
3    | conservative_low_T_v1  | 0.40 | 10.0    | 0.78  | 9.8   | 15.1M  | 765.1  | ✓ DONE
4    | baseline               | 0.50 | 12.0    | 0.71  | 14.3  | 12.0M  | 695.7  | ✓ DONE

Hidden test results will be revealed after competition ends.
Current ranking is based on public test only (60% of final score).

Click submission name for full details (system.yaml, policy.yaml, analysis)
```

---

## Analysis Summary (Auto-generated)

```
========================================
SUBMISSION ANALYSIS
========================================
Submission: aggressive_high_temp_v1
Problem: octahedron_assembly_v1
Status: COMPLETED

System Design:
  Temperature: 0.75 (50% above default)
  Concentration: 0.0015 (50% above default)
  D0 (same-type): 16.0 (33% above default)
  D0 (cross-type): 12.0 (100% above default)

Policy:
  Hysteresis: 300-500 (50-70% below default)
  Cooldown: 3000-5000 (50-70% below default)
  Check frequency: Every 50 steps (2x default)

Results:
  Final yield: 0.81 (81% of monomers in target octahedrons)
  Total state changes: 1,247
  Work per yield: 15.2 kT
  Time to 50% yield: 3.5M timesteps
  Simulation time: 45.2 hours

Score breakdown:
  final_yield * 1000 = 810.0
  - work_per_yield * 0.1 = -1.52
  - time_to_50pct / 1000 = -3.5
  = 792.3 (RANK: #2)

Strategy assessment:
  ✓ Fast assembly (time_to_50% = 3.5M, best)
  ✓ Good final yield (0.81, 2nd best)
  ⚠ Higher work_per_yield (15.2, higher than optimal)

  Your aggressive policy achieved rapid assembly due to:
  - High temperature → fast dynamics
  - Short hysteresis → quick conflict detection
  - Short cooldowns → rapid re-attempts

  Tradeoff: Higher work_per_yield suggests more "trial-and-error"
  than balanced approaches, but time savings compensate.

  Suggestion: Try reducing max_changes_per_step to 1 to reduce
  work while maintaining fast assembly time.

Files:
  - Full trajectory: aggressive_high_temp_v1/results/dump.octahedron.lammpstrj
  - Work analysis: aggressive_high_temp_v1/results/work_timeseries.csv
  - Yield plot: aggressive_high_temp_v1/results/yield_vs_time.png
```

---

## Summary

This example shows how **three files work together**:

1. **problem.yaml** (fixed) - defines the challenge and constraints
2. **system.yaml** (user designs) - physical parameters and potentials
3. **policy.yaml** (user designs) - state-change rules

**Key benefits**:
- ✅ Users can tune everything that matters (potentials, T, policy)
- ✅ Competition is fair (everyone works within same constraints)
- ✅ Submissions are self-documenting (all choices explicit)
- ✅ Results are reproducible (just re-run the 3 files)
- ✅ Leaderboard shows design tradeoffs (fast vs efficient, hot vs cold, etc.)

**This architecture supports**:
- Beginners: just write policy.yaml (use defaults)
- Intermediate: tune system.yaml (optimize potentials, T)
- Advanced: full control over everything
