## Advanced Open-Design Competition Interface

### 🔴 FIRST STEP: Checkout the Correct Branch

**You MUST be on the `advance-interface` branch to access this competition format.**

If you cloned the repository, checkout this branch:
```bash
git checkout advance-interface
```

If you're starting fresh:
```bash
git clone https://github.com/Livia33g/lammps-state-change.git
cd lammps-state-change
git checkout advance-interface
```

**The `advance/` directory and `problems/problem-001-ksat-advanced/` only exist on this branch.**

---

### Competition Format

This directory defines an **advanced** competition format where participants
submit **one Python file** that specifies:

- how to **encode** an abstract problem instance (e.g. a 3‑SAT formula) into a
  patchy-particle LAMMPS system,
- how to **design** a custom state-change **C++ fix** for that encoding, and
- how to **decode** a solution from the resulting trajectory.

### ⚠️ CRITICAL: Strict Guidelines and Auto-Disqualification

**Your submission MUST:**

1. **Pass the local checker** (`check_submission.py`) before submitting
   - If the checker fails, your submission will **NOT run on our machines**
   - Failed submissions are **automatically disqualified**

2. **Follow the exact API** defined in `submission_template.py`:
   - `API_VERSION = "v1"` (do not change)
   - Class name: `StateChangeSolution` (exact spelling)
   - Method signatures: `encode()`, `design_policy(system_meta)`, `decode()` (exact names and arguments)
   - Return types: Must match the documented schemas

3. **Write files only in the provided `work_dir`**:
   - Create `work_dir/simulation/` with LAMMPS input files
   - Create `work_dir/generated/` with C++ fix files
   - **Never** read or write outside `work_dir`

4. **Follow C++ fix naming rules**:
   - Files must start with `fix_state_change_`
   - Must be valid C++ that compiles with our LAMMPS build
   - Must not conflict with base fixes (dimer, octahedron, ksat variants)

5. **Respect resource limits**:
   - Maximum particles, steps, and walltime are specified per problem
   - Exceeding limits results in disqualification

**If your submission fails to run on our machines for ANY reason (import error, missing files, compilation failure, runtime error, etc.), it is automatically disqualified.**

### Submission Process

Your single file will be copied directly from a Google Form text field and
saved as `submission.py`. Our tools will then:

1. Import your `StateChangeSolution` class (must succeed or disqualify)
2. Call `encode()` to generate the LAMMPS system under a sandboxed `work_dir` (must succeed or disqualify)
3. Call `design_policy()` to generate C++ fix files (must succeed or disqualify)
4. Compile LAMMPS with your fix (must succeed or disqualify)
5. Run the simulation using our standard workflow (must complete or disqualify)
6. Call `decode()` to obtain your claimed solution (must succeed or disqualify)
7. Verify your decoded solution satisfies the problem constraints (must pass or disqualify)

### Submit your solution (Google Form)

Use this form for the advanced competition submissions:

`https://docs.google.com/forms/d/e/1FAIpQLSd6NBCwdrk_zeaSg8l27eCohOHZ_ZJpdmSYZyqaSlzj6QI8qg/viewform?usp=sharing`

### Files Provided

- **`submission_template.py`** – Template with required class structure and method signatures
- **`check_submission.py`** – **RUN THIS BEFORE SUBMITTING** to verify your code structure
- **`example_submission.py`** – Complete working example (one possible strategy; you can use a different approach)

### Before You Submit

**ALWAYS run the checker:**
```bash
python advance/check_submission.py --submission your_submission.py --problem problems/problem-001-ksat-advanced/problem.json
```

If the checker reports any errors, **fix them before submitting**. Submissions that fail the checker will not run on our side and will be automatically disqualified.


