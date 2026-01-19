# Example Policy YAML Files

This document shows concrete examples of what policy.yaml files would look like for different problems.

## Example 1: Octahedron "Avoid Same Type" Policy

This is your current octahedron fix translated to YAML format:

```yaml
# submissions/octahedron/avoid_same_v1/policy.yaml

metadata:
  name: "avoid_same_type_v1"
  author: "livia"
  description: "Monomers change type when touching same-type neighbors. Higher-ID monomer changes to a different type."
  version: "1.0"
  created: "2026-01-19"

parameters:
  check_every: 100                    # Check for state changes every N steps
  max_changes_per_step: 1             # Global rate limiter
  patch_types: [2, 3, 4, 5]           # Which atom types are patches
  cutoff: 2.5                         # Distance cutoff for "contact"

rules:
  - name: "unreacted_type1_collision"
    description: "Type 1 monomers touching type 1 → higher ID changes to {3,4,5}"

    trigger:
      when: "same_type_contact"
      my_types: [1]                   # I am effective type 1 (LAMMPS type 2)
      neighbor_types: [1]             # Touching effective type 1
      cutoff: 2.5                     # Use parameter above
      hysteresis:
        enabled: true
        threshold: 1000               # Must be in contact for 1000 steps
        mode: "consecutive"           # Consecutive steps, not cumulative

    action:
      selector:
        method: "higher_molecule_id"  # Which monomer changes
        # Compares molecule IDs of all conflicting neighbors
        # Only the highest-ID monomer will change

      change_type:
        method: "random_from_list"
        options: [3, 4, 5]            # Pick uniformly from these
        # Note: Type 1 can evolve to any of {3, 4, 5}

      probability: 1.0                # Always change (if selected)

    cooldown:
      base_steps: 10000               # Base cooldown period
      jitter:
        enabled: true
        range: [0.3, 2.0]             # Multiply base by random [0.3, 2.0]
        per_molecule: true            # Same jitter for all atoms in molecule

  - name: "reacted_same_type_conflict"
    description: "Types {3,4,5} touching same type → higher ID changes to DIFFERENT type"

    trigger:
      when: "same_type_contact"
      my_types: [3, 4, 5]             # I am one of these
      neighbor_types: "same_as_mine"  # Touching MY specific type
      cutoff: 2.5
      hysteresis:
        enabled: true
        threshold: 1000
        mode: "consecutive"

    action:
      selector:
        method: "higher_molecule_id"

      change_type:
        method: "exclude_current_type"  # CRITICAL: must change to DIFFERENT type
        available_types: [3, 4, 5]      # Can become any of these...
        # ... but EXCLUDES my current type automatically
        # Example: If I'm type 3, picks from {4, 5} with 50% each

      probability: 1.0

    cooldown:
      base_steps: 10000
      jitter:
        enabled: true
        range: [0.3, 2.0]
        per_molecule: true

synchronization:
  forward_comm: true                  # Enable MPI ghost atom synchronization
  consistency_sweep:
    enabled: true                     # Fix "split brain" molecules across processors
    timestamp_threshold: 100          # Only update if timestamp differs by 100+ steps
    respect_cooldown: true            # Don't override molecules in cooldown
```

## Example 2: Dimer "Simple Flip" Policy

Your current dimer fix:

```yaml
# submissions/dimer/simple_flip_v1/policy.yaml

metadata:
  name: "simple_red_blue_flip"
  author: "livia"
  description: "Red patches flip to blue when touching blue patches"
  version: "1.0"

parameters:
  check_every: 100
  cutoff: 1.6
  patch_types: [2, 3]                 # Red (2) and Blue (3)

rules:
  - name: "red_touches_blue"
    description: "Red (type 2) patches flip to blue (type 3) when touching blue"

    trigger:
      when: "different_type_contact"  # Different from previous examples
      my_types: [2]                   # I am Red
      neighbor_types: [3]             # Touching Blue
      cutoff: 1.6
      hysteresis:
        enabled: true
        threshold: 5                  # Must be in contact for 5 checks (500 steps)
        mode: "consecutive"

    action:
      selector:
        method: "all"                 # All Red monomers touching Blue can change

      change_type:
        method: "fixed"
        new_type: 3                   # Always become Blue

      probability: 1.0                # 100% flip probability

    cooldown:
      base_steps: 0                   # No cooldown for simple dimer
      jitter:
        enabled: false
```

## Example 3: Dimer KSAT "Catalyzed Flip" Policy

A → C catalyzed by B contact:

```yaml
# submissions/dimer_ksat/catalyzed_flip_v1/policy.yaml

metadata:
  name: "A_to_C_catalyzed_by_B"
  author: "livia"
  description: "A molecules flip to C when in contact with B molecules"
  version: "1.0"

parameters:
  check_every: 100
  cutoff: 0.7
  patch_types: [2, 3, 4]              # A (2), B (3), C (4)

rules:
  - name: "A_catalyzed_by_B"
    description: "A becomes C when touching B (catalysis)"

    trigger:
      when: "catalyzed_contact"       # Special: requires different-type neighbor
      my_types: [2]                   # I am A
      catalyst_types: [3]             # Need to touch B (catalyst)
      cutoff: 0.7
      hysteresis:
        enabled: true
        threshold: 5
        mode: "consecutive"

    action:
      selector:
        method: "all"                 # All A's touching B can change

      change_type:
        method: "fixed"
        new_type: 4                   # A → C (type 4)

      probability: 1.0

    cooldown:
      base_steps: 0
      jitter:
        enabled: false

  # Could add reverse reaction C → A if desired:
  # - name: "C_reverse_to_A"
  #   trigger:
  #     when: "no_catalyst_contact"
  #     my_types: [4]
  #     catalyst_types: [3]
  #     mode: "absence"              # Change when NOT touching B
  #   ...
```

## Example 4: Advanced Custom Logic (Expert Mode)

For users who need logic that can't be expressed in YAML:

```yaml
# submissions/octahedron/custom_weighted_v1/policy.yaml

metadata:
  name: "custom_weighted_selection"
  author: "expert_user"
  description: "Uses neighborhood composition to weight type selection"
  version: "1.0"
  type: "cpp_template"                # Bypass YAML generator, use C++ directly

parameters:
  check_every: 100
  cutoff: 2.5
  patch_types: [2, 3, 4, 5]

# For advanced users: inject C++ code directly
cpp_code:
  includes: |
    #include <map>
    #include <algorithm>

  helper_functions: |
    // Custom helper function
    int compute_weighted_type(int my_mol, int my_type,
                              std::map<int, int>& neighbor_types) {
      // Count neighbor types
      int count_3 = 0, count_4 = 0, count_5 = 0;
      for (auto& pair : neighbor_types) {
        if (pair.second == 3) count_3++;
        if (pair.second == 4) count_4++;
        if (pair.second == 5) count_5++;
      }

      // Weight selection based on what's LEAST common in neighborhood
      // This promotes diversity
      if (count_3 <= count_4 && count_3 <= count_5) return 3;
      if (count_4 <= count_5) return 4;
      return 5;
    }

  change_logic: |
    // Custom state change logic
    if (my_eff_type >= 3 && my_eff_type <= 5) {
      // Build neighbor type map
      std::map<int, int> neighbor_types;
      for (int j = 0; j < nall; j++) {
        if (is_neighbor(i, j)) {
          neighbor_types[j] = effective_type[j];
        }
      }

      // Use custom function to pick type
      new_type = compute_weighted_type(molecule[i], my_eff_type, neighbor_types);

      // IMPORTANT: Still exclude current type
      if (new_type == my_eff_type) {
        // Pick randomly from other two
        std::vector<int> others;
        for (int t : {3, 4, 5}) {
          if (t != my_eff_type) others.push_back(t);
        }
        new_type = others[random.uniform() * others.size()];
      }
    }
```

## Example 5: Multi-Rule Policy (Complex)

A policy with multiple triggers for different scenarios:

```yaml
# submissions/octahedron/multi_rule_v1/policy.yaml

metadata:
  name: "multi_rule_adaptive"
  author: "livia"
  description: "Different rules based on local density and composition"
  version: "1.0"

parameters:
  check_every: 100
  cutoff: 2.5
  patch_types: [2, 3, 4, 5]

rules:
  # Rule 1: Low-density regions (few neighbors) → more aggressive flipping
  - name: "low_density_fast_flip"
    description: "In sparse regions, flip quickly to find partners"

    trigger:
      when: "same_type_contact"
      my_types: [3, 4, 5]
      neighbor_types: "same_as_mine"
      cutoff: 2.5

      # Additional condition: low total coordination
      coordination_filter:
        total_neighbors: "<= 2"       # Sparse neighborhood

      hysteresis:
        enabled: true
        threshold: 500                # Shorter wait in sparse regions
        mode: "consecutive"

    action:
      selector:
        method: "higher_molecule_id"
      change_type:
        method: "exclude_current_type"
        available_types: [3, 4, 5]
      probability: 1.0

    cooldown:
      base_steps: 5000                # Shorter cooldown in sparse regions

  # Rule 2: High-density regions (many neighbors) → slower, more careful
  - name: "high_density_careful_flip"
    description: "In dense regions, wait longer to avoid disrupting clusters"

    trigger:
      when: "same_type_contact"
      my_types: [3, 4, 5]
      neighbor_types: "same_as_mine"
      cutoff: 2.5

      coordination_filter:
        total_neighbors: "> 2"        # Dense neighborhood

      hysteresis:
        enabled: true
        threshold: 2000               # Longer wait in dense regions
        mode: "consecutive"

    action:
      selector:
        method: "higher_molecule_id"
      change_type:
        method: "exclude_current_type"
        available_types: [3, 4, 5]
      probability: 1.0

    cooldown:
      base_steps: 20000               # Longer cooldown in dense regions

  # Rule 3: Type 1 (unreacted) always uses same rule
  - name: "type1_standard"
    description: "Unreacted monomers follow standard rule regardless of density"

    trigger:
      when: "same_type_contact"
      my_types: [1]
      neighbor_types: [1]
      cutoff: 2.5
      hysteresis:
        enabled: true
        threshold: 1000
        mode: "consecutive"

    action:
      selector:
        method: "higher_molecule_id"
      change_type:
        method: "random_from_list"
        options: [3, 4, 5]
      probability: 1.0

    cooldown:
      base_steps: 10000
      jitter:
        enabled: true
        range: [0.3, 2.0]
```

## Policy Validation

The framework will validate policies against a JSON schema. Invalid policies are rejected with helpful error messages:

```bash
$ python scripts/validate_policy.py submissions/octahedron/my_policy/policy.yaml

✗ VALIDATION ERROR: policy.yaml
  Line 23: rules[0].action.change_type.method = "invalid_method"

  Expected one of: ["fixed", "random_from_list", "exclude_current_type", "weighted"]

  Did you mean "exclude_current_type"?

  See docs/reference/yaml_schemas.md for valid values.
```

## Benefits of YAML Format

1. **No C++ knowledge required** - declarative, readable
2. **Validation** - catch errors before compilation
3. **Versioning** - easy to track changes in git
4. **Sharing** - policies are self-documenting
5. **Comparison** - diff two policies easily
6. **Portability** - works across problems (if logic is general)

## Next: Generating C++ from YAML

The framework's `fix_codegen.py` will:

1. Parse YAML
2. Validate against schema
3. Fill Jinja2 template
4. Generate `.cpp` and `.h` files
5. Compile into LAMMPS

User never sees C++ code unless they want to inspect it!
