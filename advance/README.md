## Advanced open-design competition interface

This directory defines an **advanced** competition format where participants
submit **one Python file** that specifies:

- how to **encode** an abstract problem instance (e.g. a 3‑SAT formula) into a
  patchy-particle LAMMPS system,
- how to **design** a custom state-change **C++ fix** for that encoding, and
- how to **decode** a solution from the resulting trajectory.

Your single file will be copied directly from a Google Form text field and
saved as `submission.py`. Our tools will then:

1. Import your `StateChangeSolution` class.
2. Call `encode()` to generate the LAMMPS system under a sandboxed `work_dir`.
3. Call `design_policy()` to generate C++ fix files.
4. Compile LAMMPS with your fix (on our side).
5. Run the simulation using our standard workflow.
6. Call `decode()` to obtain your claimed solution, which we will verify.

If any of these steps fails or violates the rules, the submission is
**automatically disqualified**.

For local development, we provide:

- `submission_template.py` – a template you can copy and modify.
- `check_submission.py` – a structural checker you should run before
  submitting. If it fails, your submission will not run on our side.


