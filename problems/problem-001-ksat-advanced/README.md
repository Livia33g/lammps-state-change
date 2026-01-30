# Problem 001: Abstract 2-SAT Dimer Instance

## Problem Statement

You are given a **2-SAT formula** encoded as a dimer patchy-particle system with two independent conditions. Your task is to:

1. **Encode** this 2-SAT problem into a patchy-particle LAMMPS system using a dimer structure
2. **Design** a custom state-change policy (C++ fix) that implements the logical rules for each condition
3. **Decode** the solution by measuring the concentration of true/false particles at the end of the simulation

## Problem Definition

This is a **2-SAT problem** with two independent conditions:

### Condition B (TF - True OR False)
- **Monomer B**: A condition monomer with two patch faces
- **Monomer A**: False (switchable) - can flip to C
- **Monomer C**: True (non-switchable)
- **Rule**: Allows only 1 false (A) monomer to be attached at a time
  - If **A, A** are attached simultaneously to B (on both faces), flip the one with **higher molecule ID** to C (true)
  - If **C, C** are attached, nothing happens
  - If **A, C** are attached, nothing happens

### Condition D (TT - True AND True)
- **Monomer D**: A condition monomer with two patch faces
- **Monomer E**: False (switchable) - can flip to F
- **Monomer F**: True (non-switchable)
- **Rule**: Both need to be true (F)
  - If even **one false (E)** attaches to D, switch it to type F (true), regardless of what is attached to the other side

### Independence
The two sets **A, B, C** and **D, E, F** are **independent** and do **not interact** with each other.

## The Solution

**If we set B to True (ON), both rules are satisfied:**

- **Rule 1 (B=TF)**: (Anything ∨ True) = True ✓
- **Rule 2 (D=TT)**: (Anything ∨ True) = True ✓

## Your Challenge

**Design a dimer system** following the structure similar to `state-change/dimer_ksat/twosideB_twin`:

- Use the same potentials and structure as the reference implementation
- Implement two independent channels: **A,B,C** and **D,E,F**
- Each condition monomer (B and D) has two patch faces
- Implement the state-change logic described above in your C++ fix

**Constraints:**

- Maximum 500 particles
- Maximum 2,000,000 simulation steps
- Maximum 24 hours walltime
- Your `submission.py` must pass `advance/check_submission.py` **before** you submit
- If your submission fails to run on our machines, it is **automatically disqualified**

## Scoring

- **Primary metric:** Concentration of true particles (C and F) at the end of simulation (higher is better)
  - The solution B=True should maximize this concentration
- **Tie-breakers:** Total state changes, time to solution

## Verification

We will independently verify that your decoded assignment satisfies both conditions. The solution is measured by the concentration of true/false particles at the end of the simulation. If your solution does not satisfy both conditions, your submission is disqualified regardless of efficiency.

## Reference Structure

See `state-change/dimer_ksat/twosideB_twin` for the general structure and how it works. The logic and structure should be similar—use the same potentials and structure, but implement the specific state-change rules described above.

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
