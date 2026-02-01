# Problem 001: Abstract 2-SAT Dimer Instance

## Problem Statement

You are given a **2-SAT problem** with two independent logical conditions. Your task is to:

1. **Encode** this abstract 2-SAT problem into a patchy-particle LAMMPS system (you design the encoding)
2. **Design** a custom state-change policy (C++ fix) that implements the logical rules
3. **Decode** the solution from your simulation results

## Abstract Problem Definition

This is a **2-SAT problem** with two independent conditions that must be satisfied:

### Condition 1: TF (True OR False)
- **Logical meaning**: At least one of the two inputs must be True
- **Truth table**: 
  - (False, False) → **Unsatisfied** (both false)
  - (True, False) → Satisfied
  - (False, True) → Satisfied  
  - (True, True) → Satisfied

**State-change rule**: If both inputs are False simultaneously, flip one of them to True (you decide which one and how).

### Condition 2: TT (True AND True)
- **Logical meaning**: Both inputs must be True
- **Truth table**:
  - (False, False) → **Unsatisfied**
  - (True, False) → **Unsatisfied** (one false)
  - (False, True) → **Unsatisfied** (one false)
  - (True, True) → Satisfied

**State-change rule**: If either input is False, flip it to True (regardless of the other input's state).

### Independence
The two conditions are **independent** - they do not interact with each other. You can encode them as separate molecular channels or in any way you choose.

## The Solution

**If we set Condition 1 to True, both conditions are satisfied:**
- **Condition 1 (TF)**: (Anything ∨ True) = True ✓
- **Condition 2 (TT)**: (Anything ∨ True) = True ✓

(Note: This is one possible solution. Your encoding may have different solutions.)

## Your Challenge

**Design your own encoding and state-change policy** to solve this 2-SAT problem. You have complete freedom to:

- Choose particle types and geometry
- Design interaction potentials
- Implement state-change logic in your C++ fix
- Define how to decode the solution

### Example Encoding (For Reference Only)

One possible encoding (shown as an example - you don't have to use this):

- **Condition 1 (TF)**: 
  - Monomer types: A (false, switchable), B (condition with two faces), C (true, non-switchable)
  - Rule: If A,A attached to B on both faces → flip higher-ID A→C
  
- **Condition 2 (TT)**:
  - Monomer types: E (false, switchable), D (condition with two faces), F (true, non-switchable)
  - Rule: If E attaches to D → flip E→F

This is just **one example** - you are free to design a completely different encoding!

**Constraints:**
- Maximum 500 particles
- Maximum 2,000,000 simulation steps
- Maximum 24 hours walltime
- Your `submission.py` must pass `advance/check_submission.py` **before** you submit
- If your submission fails to run on our machines, it is **automatically disqualified**

## Scoring

- **Primary metric:** Solution correctness (does your decoded solution satisfy both conditions?)
- **For binary problems:** Only solved solutions are ranked, by total work (lower is better)
- **For continuous problems:** Ranked by your player-defined score (see decoding policy)

## Verification

We will independently verify that your decoded solution satisfies both conditions. If your solution does not satisfy both conditions, your submission is disqualified regardless of efficiency.

## Submission Format

Submit a single Python file (`submission.py`) that implements the `StateChangeSolution` class from `advance/submission_template.py`. Copy your file content into the Google Form text field.

### Submit (Google Form)

Submit your `submission.py` by pasting it into the form:

`https://docs.google.com/forms/d/e/1FAIpQLSd6NBCwdrk_zeaSg8l27eCohOHZ_ZJpdmSYZyqaSlzj6QI8qg/viewform?usp=sharing`

**Before submitting, run:**
```bash
python advance/check_submission.py --submission your_submission.py --problem problems/problem-001-ksat-advanced/problem.json
```

If the checker fails, your submission will not run on our side and will be automatically disqualified.

## Example Submission

See `advance/example_submission_2sat.py` for a complete working example that demonstrates one possible encoding strategy. **You are free to use a completely different approach!**
