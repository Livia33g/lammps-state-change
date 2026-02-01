# Compilation Tips and Common Sandbox Failures

This document summarizes everything we've learned about what causes submissions to fail in the sandbox and what works correctly.

## ⚠️ Critical: Component Coherence

**Your submission has three components that MUST be coherent:**

1. **`encode()`** - Creates the LAMMPS system (atom types, geometry, interactions)
2. **`design_policy()`** - Generates C++ fix that implements state-change logic
3. **`decode()`** - Reads simulation results to extract solution

### Why Coherence Matters

If these components are **not coherent**, your submission will produce **null/invalid results**:

- ❌ **encode()** creates type 2 for A, but **design_policy()** looks for type 3 → **No state changes occur**
- ❌ **encode()** creates geometry with 3 patches, but **design_policy()** expects 6 patches → **Crashes or wrong behavior**
- ❌ **design_policy()** flips type 2→4, but **decode()** counts type 5 → **Wrong solution measurement**
- ❌ **encode()** uses channel ABC (types 1-5), but **design_policy()** checks channel EFD (types 8-12) → **No flips happen**

### ✅ How to Ensure Coherence

1. **Define atom type mapping once** and use it consistently:
   ```python
   # In __init__ or encode():
   self.type_map = {
       "A": 2, "B_face1": 3, "C": 4, "B_face2": 5,
       "E": 8, "D_face1": 9, "F": 10, "D_face2": 11
   }
   ```

2. **Pass metadata from encode() to design_policy()**:
   ```python
   def encode(self):
       # ... create system ...
       return {
           "atom_types": self.type_map,  # Pass to design_policy
           "geometry": "1core_twosideB_twins"
       }
   
   def design_policy(self, system_meta):
       type_map = system_meta["atom_types"]  # Use same mapping
       # Generate C++ with correct types
   ```

3. **Use same types in decode()**:
   ```python
   def decode(self):
       # Count types 4 (C) and 10 (F) - same as encode() created
       n_C = count_type(4)
       n_F = count_type(10)
   ```

## 🔴 Common Sandbox Failures

### 1. Missing Mass Definitions

**Error:**
```
ERROR: Not all per-type masses are set ...
```

**Cause:** LAMMPS requires `mass` commands for all atom types before `velocity` or dynamics.

**Fix:**
```lammps
# In your in.* file, BEFORE velocity or fix commands:
mass 1 0.6
mass 2 0.1
mass 3 0.1
# ... define mass for ALL types you use (1-12 for 2-SAT problem)
```

**What works:**
- ✅ Define `mass` for every atom type used in `data.*`
- ✅ Use `atom_style full` (recommended) - requires masses for all types
- ✅ Place `mass` commands after `read_data` but before `velocity` or `fix`

**What doesn't work:**
- ❌ Assuming LAMMPS will use default masses
- ❌ Defining masses only for some types
- ❌ Using `atom_style atomic` without proper mass definitions

### 2. Incorrect FixStyle Registration

**Error:**
```
ERROR: Unknown fix style 'state/change/yourname'
```

**Cause:** Fix header doesn't use the `#ifdef FIX_CLASS` pattern that our LAMMPS build requires.

**Fix:**
```cpp
// In fix_state_change_<name>.h:
#ifdef FIX_CLASS
// clang-format off
FixStyle(state/change/<name>,FixStateChangeName);
// clang-format on
#else

#ifndef LMP_FIX_STATE_CHANGE_NAME_H
#define LMP_FIX_STATE_CHANGE_NAME_H

#include "fix.h"
// ... rest of header ...

#endif
#endif
```

**What works:**
- ✅ Using `#ifdef FIX_CLASS` pattern in header
- ✅ Exact spelling: `FixStyle(state/change/<name>,FixStateChangeName);`
- ✅ Class name matches: `FixStateChangeName` in header and implementation

**What doesn't work:**
- ❌ Trying to `#include "fix_style.h"` (not in our LAMMPS tree)
- ❌ Registering fix in `.cpp` file instead of header
- ❌ Wrong FixStyle syntax or class name mismatch

### 3. MPI World Access Error

**Error:**
```
Segmentation fault or MPI error
```

**Cause:** Using `world` variable directly instead of `comm->world`.

**Fix:**
```cpp
// WRONG:
MPI_Allreduce(MPI_IN_PLACE, &nchanges, 1, MPI_INT, MPI_SUM, world);

// CORRECT:
MPI_Allreduce(MPI_IN_PLACE, &nchanges, 1, MPI_INT, MPI_SUM, comm->world);
```

**What works:**
- ✅ Using `comm->world` for MPI operations
- ✅ Using `comm->me` for rank ID
- ✅ Using `comm->nprocs` for number of processes

**What doesn't work:**
- ❌ Direct access to `world` variable
- ❌ Assuming MPI_COMM_WORLD is available

### 4. Periodic Boundary Condition (PBC) Errors

**Error:**
```
NaN values or division by zero
```

**Cause:** Not checking if `prd` (periodic box dimensions) is valid before division.

**Fix:**
```cpp
// WRONG:
dx -= prd[0] * std::round(dx / prd[0]);  // Can divide by zero!

// CORRECT:
if (prd && prd[0] > 0.0) {
    dx -= prd[0] * std::round(dx / prd[0]);
}
```

**What works:**
- ✅ Checking `prd != NULL` before use
- ✅ Checking `prd[i] > 0.0` before division
- ✅ Using `domain->prd` to get periodic dimensions

**What doesn't work:**
- ❌ Dividing by `prd[i]` without checking if it's zero
- ❌ Assuming periodic boundaries are always enabled
- ❌ Not handling non-periodic dimensions

### 5. Atom Type Mismatches

**Error:**
```
No state changes occur, or wrong particles flip
```

**Cause:** Atom types in C++ fix don't match types created in `encode()`.

**Fix:**
- Use constants in C++ that match your Python encoding:
  ```cpp
  // In fix_state_change_2sat.cpp:
  constexpr int TYPE_A = 2;      // Must match encode() output
  constexpr int TYPE_C = 4;
  constexpr int TYPE_E = 8;
  constexpr int TYPE_F = 10;
  ```

**What works:**
- ✅ Using named constants for atom types
- ✅ Documenting type mapping in comments
- ✅ Verifying types match between encode() and design_policy()

**What doesn't work:**
- ❌ Hardcoding type numbers without checking encode()
- ❌ Using different types in different components
- ❌ Assuming types are the same as examples

### 6. Missing Group Definitions

**Error:**
```
ERROR: Group 'patches' not found
```

**Cause:** Fix references a group that doesn't exist in LAMMPS input.

**Fix:**
```lammps
# In your in.* file, BEFORE the fix command:
group patches type 2 3 4 5 8 9 10 11
fix sc patches state/change/2sat 0.7 patches 5
```

**What works:**
- ✅ Defining groups before using them in fixes
- ✅ Using `group` command with correct type list
- ✅ Matching group name in fix command

**What doesn't work:**
- ❌ Referencing groups that don't exist
- ❌ Using wrong type numbers in group definition
- ❌ Defining group after fix command

### 7. File Naming Mismatches

**Error:**
```
ERROR: listed fix file does not exist
```

**Cause:** `design_policy()` returns filenames that don't match actual files created.

**Fix:**
```python
def design_policy(self, system_meta):
    # Create files
    cpp_file = self.gen_dir / 'fix_state_change_2sat.cpp'
    h_file = self.gen_dir / 'fix_state_change_2sat.h'
    # ... write files ...
    
    # Return EXACT filenames (relative to generated/)
    return {
        "fix_files": [
            "fix_state_change_2sat.cpp",  # Must match actual filename
            "fix_state_change_2sat.h"
        ]
    }
```

**What works:**
- ✅ Filenames in return dict match actual files
- ✅ Files start with `fix_state_change_`
- ✅ Both `.cpp` and `.h` files are listed

**What doesn't work:**
- ❌ Returning filenames that don't exist
- ❌ Typos in filenames (e.g., `fix_state_chage_` instead of `fix_state_change_`)
- ❌ Missing `.cpp` or `.h` file

### 8. LAMMPS Input File Errors

**Error:**
```
ERROR: Invalid command in input script
```

**Common causes:**
- Missing required commands (e.g., `read_data` before `fix`)
- Wrong command order
- Invalid syntax

**What works:**
- ✅ Standard LAMMPS command order:
  1. `units`, `atom_style`, `boundary`
  2. `read_data`
  3. `mass` (for all types)
  4. `group`
  5. `pair_style`, `pair_coeff`
  6. `neighbor`
  7. `fix` (dynamics and state-change)
  8. `thermo`, `dump`
  9. `run`

**What doesn't work:**
- ❌ Using commands before required setup
- ❌ Missing `read_data` before defining interactions
- ❌ Defining groups after fixes that use them

### 9. Smoke Simulation Failures

**Error:**
```
ERROR: Smoke simulation failed
```

**Cause:** Simulation crashes during the 1000-step smoke test (runtime error, not compilation).

**Common causes:**
- Invalid pair coefficients
- Missing required fixes
- Division by zero in fix code
- Array bounds errors

**What works:**
- ✅ Testing your fix locally before submission
- ✅ Checking for null pointers before dereferencing
- ✅ Validating array indices before access
- ✅ Using safe distance calculations with PBC checks

**What doesn't work:**
- ❌ Assuming compilation success means runtime success
- ❌ Not testing with a short simulation locally
- ❌ Accessing arrays without bounds checking

## ✅ Best Practices

1. **Test locally first:**
   ```bash
   python advance/check_submission.py --submission your_submission.py --problem problems/problem-001-ksat-advanced/problem.json
   ```

2. **Use the example as a template:**
   - See `example_submission_2sat.py` for a complete working example
   - Follow the same structure and patterns

3. **Verify coherence:**
   - Check that atom types match between encode/design/decode
   - Verify geometry matches between components
   - Ensure state-change logic matches problem requirements

4. **Test with minimal system:**
   - Start with small particle counts
   - Verify state changes occur
   - Check that decode() reads correct types

5. **Read error messages carefully:**
   - LAMMPS error messages are usually very specific
   - Check line numbers in error output
   - Verify file paths and names

## 📝 Example: Working 2-SAT Submission

See `example_submission_2sat.py` for a complete working example that:
- ✅ Properly encodes the 2-SAT problem
- ✅ Generates correct C++ fix with proper FixStyle registration
- ✅ Implements coherent state-change logic
- ✅ Decodes solution correctly
- ✅ Passes sandbox checks

Use this as a reference for your own submission!

