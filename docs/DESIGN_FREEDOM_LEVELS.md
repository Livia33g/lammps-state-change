# Design Freedom Levels in Molecular Computing Problems

Different problems offer different levels of design freedom to participants. This document explains what users can modify at each level.

---

## 📊 Four Levels of Design Freedom

### **Level 1: Policy Only**
**What you design:** State transition rules (π in the proposal framework)
**What's fixed:** Everything else (encoding, parameters, geometry)

**Example problem:**
```json
{
  "design_freedom": {
    "level": "policy_only",
    "allowed_modifications": ["state_transitions"],
    "description": "Design when and how molecules flip states. All physical parameters are fixed."
  }
}
```

**User submits:**
- `policy.json` - State transition rules only

**Cannot modify:**
- Interaction strengths (morse_depth, etc.)
- Initial composition (# of A vs B)
- Geometry (1core_3patch is fixed)
- Species definitions

**Use case:** Beginner problems focusing on policy logic without physics tuning.

---

### **Level 2: Policy and Parameters** ⭐ (Current: problem-001)
**What you design:** State transitions + interaction strengths
**What's fixed:** Encoding (species, geometry, composition)

**Example problem:**
```json
{
  "design_freedom": {
    "level": "policy_and_parameters",
    "allowed_modifications": [
      "state_transitions",
      "interaction_strengths",
      "contact_cutoffs"
    ],
    "description": "Design policy and tune interaction parameters. Encoding is fixed."
  }
}
```

**User submits:**
- `policy.json` - State transition rules
- `params.json` - Morse depths, cutoffs, etc.

**Can modify:**
- Flip probabilities, hysteresis
- `morse_depth_AB`, `morse_depth_BC`
- `morse_alpha` (interaction steepness)
- `contact_cutoff` (trigger distance)

**Cannot modify:**
- Species labels (A, B, C fixed)
- Geometry (1core_3patch fixed)
- Initial composition (20 A, 10 B fixed)
- Temperature, density

**Use case:** Intermediate problems with two-level optimization (policy + physics).

---

### **Level 3: Policy and Encoding**
**What you design:** State transitions + molecular components
**What's fixed:** Only the computational goal

**Example problem:**
```json
{
  "design_freedom": {
    "level": "policy_and_encoding",
    "allowed_modifications": [
      "state_transitions",
      "interaction_strengths",
      "initial_composition",
      "species_definition"
    ],
    "description": "Design the molecular encoding and policy. Goal: solve 3-SAT with N=50 variables."
  }
}
```

**User submits:**
- `policy.json` - State transition rules
- `encoding.json` - Species definitions, composition, interactions
- `params.json` - Tunable parameters

**Can modify:**
- Which species exist (A, B, C, D, ...)
- How many of each species initially
- Who interacts with whom (topology)
- All interaction strengths

**Cannot modify:**
- Computational goal (e.g., "solve 3-SAT instance X")
- Geometry options (must choose from allowed list)

**Use case:** Advanced problems where encoding is part of the challenge.

**Example use case:** "Design a molecular computer to solve this specific 3-SAT instance. You choose how to encode variables and clauses."

---

### **Level 4: Full System Design**
**What you design:** Everything (policy, encoding, even geometry)
**What's fixed:** Only the abstract computational task

**Example problem:**
```json
{
  "design_freedom": {
    "level": "full_system",
    "allowed_modifications": [
      "state_transitions",
      "interaction_strengths",
      "initial_composition",
      "geometry_choice",
      "species_definition",
      "interaction_topology",
      "particle_count"
    ],
    "description": "Open-ended design. Task: sample from Ising model on 10×10 grid. Design the full molecular system."
  }
}
```

**User submits:**
- `policy.json` - State transition rules
- `encoding.json` - Full molecular system definition
- `geometry.json` - Custom geometry (or choice from library)
- `params.json` - All tunable parameters

**Can modify:**
- Everything!
- Geometry (1core_3patch vs 1core_6patch vs custom)
- Particle count (10 vs 50 vs 100)
- Species topology (which types interact)

**Cannot modify:**
- The abstract task ("sample from this distribution")

**Use case:** Research-level open problems or competition grand challenges.

**Example use case:** "Design a molecular system that finds ground states of frustrated spin systems. Full design freedom."

---

## 🎯 How This Maps to the Proposal

From the proposal's **Encode / Design / Decode** framework:

| Level | ENCODE | DESIGN | DECODE |
|-------|--------|--------|--------|
| **1: Policy Only** | Fixed by problem | **User: π only** | Fixed by problem |
| **2: Policy + Params** | Fixed by problem | **User: π + θ** | Fixed by problem |
| **3: Policy + Encoding** | **User: species + H** | **User: π + θ** | Fixed by problem |
| **4: Full System** | **User: full system** | **User: π + θ + H** | **User: readout** |

Where:
- **π** = perturbation policy (state transitions)
- **θ** = interaction parameters (morse depths, etc.)
- **H** = Hamiltonian (binding landscape, species definitions)

---

## 📝 Example Problems by Level

### Level 1: Policy Only
**Problem:** "Convert A→C via B using fixed interactions"
```
Allowed: flip_probability, hysteresis_checks, trigger conditions
Fixed: morse_depth=1.0, composition=20A+10B, geometry=1core_3patch
```

### Level 2: Policy + Parameters (Current: problem-001)
**Problem:** "Optimize A→C conversion efficiency"
```
Allowed: policy + morse_depth_AB, morse_depth_BC, morse_alpha, contact_cutoff
Fixed: species (A,B,C), composition (20A+10B), geometry (1core_3patch)
```

### Level 3: Policy + Encoding
**Problem:** "Solve 3-SAT instance with 20 variables, 40 clauses"
```
Allowed: policy + species design + composition + interaction topology
Fixed: 3-SAT instance, available geometries (choose from 5 options)
Goal: Decode satisfying assignment from steady-state
```

### Level 4: Full System
**Problem:** "Sample Boltzmann distribution for 2D Ising model (10×10)"
```
Allowed: Full system design (geometry, species, interactions, policy, readout)
Fixed: Target distribution only
Goal: Achieve correct sampling statistics
```

---

## 🚀 Implementation Notes

### For Problem Creators

Specify `design_freedom` in `problem.json`:

```json
{
  "design_freedom": {
    "level": "policy_and_parameters",
    "allowed_modifications": ["state_transitions", "interaction_strengths"],
    "description": "..."
  },
  "constraints": {
    "fixed_parameters": ["geometry", "initial_composition"],
    "tunable_parameters": [
      {"name": "morse_depth_AB", "min": 0.1, "max": 5.0}
    ]
  }
}
```

### For Participants

Check `design_freedom.level` to understand what you can submit:

- **policy_only**: Submit `policy.json`
- **policy_and_parameters**: Submit `policy.json` + `params.json`
- **policy_and_encoding**: Submit `policy.json` + `encoding.json` + `params.json`
- **full_system**: Submit full system specification

### Validation

The validation pipeline checks:
1. User didn't modify fixed parameters
2. Tunable parameters are in allowed ranges
3. Submission includes all required files for that freedom level
4. Encoding follows allowed topology (if applicable)

---

## 🎓 Pedagogical Progression

Recommended learning path:

1. **Start:** Level 1 or 2 (policy focus)
2. **Intermediate:** Level 2 (policy + physics tuning)
3. **Advanced:** Level 3 (encoding design)
4. **Expert:** Level 4 (full system design)

This mirrors the proposal's research agenda:
- **Level 1-2:** Understanding non-equilibrium protocols (π optimization)
- **Level 3:** Co-designing encoding and protocols (H + π)
- **Level 4:** Discovering new computational primitives

---

## 🔬 Future Extensions

**Possible additions:**

- **Geometry library:** Users choose from {1core_3patch, 1core_6patch, dimer, octahedron, ...}
- **Interaction templates:** Users choose from {morse, LJ, soft, patchy_morse, ...}
- **Hybrid levels:** Some parameters fixed, others tunable (fine-grained control)
- **Time-varying freedom:** Early stage (fixed), late stage (tunable)

**Example hybrid:**
```json
{
  "design_freedom": {
    "level": "custom",
    "allowed_modifications": ["state_transitions", "morse_depth_AB"],
    "fixed": ["morse_depth_BC", "composition", "geometry"]
  }
}
```

---

**Current Status:**
- ✅ Level 2 implemented (problem-001)
- ⏳ Level 1, 3, 4 (future problems)
- 📝 Schema supports all levels

**This flexibility allows the competition to grow from beginner-friendly policy design to research-level open challenges!**
