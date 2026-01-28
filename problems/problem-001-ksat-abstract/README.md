# Problem 001: Abstract 3-SAT Instance

## Problem Statement

You are given a **3-SAT formula** (conjunctive normal form with 3 literals per clause). Your task is to:

1. **Encode** this formula into a patchy-particle LAMMPS system of your own design
2. **Design** a custom state-change policy (C++ fix) that drives the system toward a satisfying assignment
3. **Decode** a truth assignment from the simulation trajectory

The formula is:

```
φ = (x₁ ∨ ¬x₂ ∨ x₃) ∧ (¬x₁ ∨ x₄ ∨ x₅) ∧ (x₂ ∨ x₃ ∨ ¬x₆) ∧ (¬x₃ ∨ ¬x₄ ∨ x₆) ∧ (x₁ ∨ x₅ ∨ ¬x₆)
```

**In CNF notation:**
```
(1, -2, 3) ∧ (-1, 4, 5) ∧ (2, 3, -6) ∧ (-3, -4, 6) ∧ (1, 5, -6)
```

Where positive integers represent variables and negative integers represent negated variables.

## Problem Definition

- **Number of variables:** 6 (x₁ through x₆)
- **Number of clauses:** 5
- **Clause structure:** Each clause contains exactly 3 literals

The problem definition is provided in `problem.json` as:

```json
{
  "sat_instance": {
    "n_variables": 6,
    "clauses": [
      [1, -2, 3],
      [-1, 4, 5],
      [2, 3, -6],
      [-3, -4, 6],
      [1, 5, -6]
    ]
  }
}
```

## Your Challenge

**Design freedom:** You have complete freedom to choose:

- **Patchy particle geometry:** How many patches? What arrangement? What types?
- **Encoding strategy:** How do variables map to particles or states? How do clauses map to interactions or constraints?
- **State-change policy:** When and how do particles flip states? What drives the system toward satisfaction?
- **Decoding strategy:** How do you read a satisfying assignment from the trajectory?

**Constraints:**

- Maximum 500 particles
- Maximum 2,000,000 simulation steps
- Maximum 24 hours walltime
- Your `submission.py` must pass `advance/check_submission.py` **before** you submit
- If your submission fails to run on our machines, it is **automatically disqualified**

## Scoring

- **Primary metric:** Work per satisfied clause (lower is better)
  - Work = number of state-change events (or equivalent energy measure)
- **Tie-breakers:** Total work, time to first satisfying assignment

## Verification

We will independently verify that your decoded assignment satisfies all clauses. If it does not, your submission is disqualified regardless of efficiency.

## Example Solution

See `advance/example_submission.py` for a complete working example that demonstrates one possible encoding/design/decode strategy. **You are not required to follow this approach**—it is provided only as a reference.

## Submission Format

Submit a single Python file (`submission.py`) that implements the `StateChangeSolution` class from `advance/submission_template.py`. Copy your file content into the Google Form text field.

**Before submitting, run:**
```bash
python advance/check_submission.py --submission your_submission.py --problem problems/problem-001-ksat-abstract/problem.json
```

If the checker fails, your submission will not run on our side and will be automatically disqualified.

