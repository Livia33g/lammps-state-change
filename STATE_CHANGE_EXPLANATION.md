# How State Changes Work in LAMMPS with Rigid Bodies

## The Core Challenge

The fundamental problem is that **LAMMPS rigid body fixes maintain internal state** about which atoms belong to which rigid body. When you change atom types directly while a rigid fix is active, LAMMPS gets confused because:
1. The rigid fix has cached which atoms are in which molecule
2. Changing types can break the molecule-to-atom mapping
3. This leads to "Lost atoms" errors or "atoms time integrated more than once" errors

## The Solution: Unfix → Change → Refix Cycle

The key insight is that we must **temporarily remove the rigid body constraint**, change the types, then **recreate the rigid body constraint**. This is done in a loop during the simulation.

## Detailed Mechanism Breakdown

### 1. **Coordination Detection** (Lines 288-289)

```lammps
compute cA all coord/atom cutoff 2.5 group patches_A
compute cB all coord/atom cutoff 2.5 group patches_B
```

**What this does:**
- `compute coord/atom` counts how many neighbors each atom has within a cutoff distance
- `cA` counts how many type-2 patches are near each atom
- `cB` counts how many type-3 patches are near each atom
- The cutoff (2.5) matches the Morse potential cutoff, so we detect when patches are close enough to interact

**Why "all" instead of the group:**
- Using `all` means the compute works for ALL atoms
- For atoms not in the group, it returns 0.0
- This allows us to use the same compute for all atoms in the variable formula

### 2. **Random Probability** (Lines 292-293)

```lammps
variable rand_seed equal 12345
variable rand_val atom random(0,1,${rand_seed})
```

**What this does:**
- Generates a random number between 0 and 1 for EACH atom
- The seed is incremented each iteration to get new random values
- This implements the probabilistic nature of state changes

### 3. **The State Change Formula** (Line 304)

```lammps
variable new_type atom "type + (type==2) * (c_cA>0.1) * (v_rand_val<0.7) - (type==3) * (c_cB>0.1) * (v_rand_val<0.7)"
```

**Breaking this down:**

The formula uses LAMMPS's arithmetic expression syntax:
- `type` - the current atom type (1, 2, or 3)
- `(type==2)` - equals 1 if type is 2, else 0
- `(c_cA>0.1)` - equals 1 if coordination > 0.1 (patch is near another type-2 patch), else 0
- `(v_rand_val<0.7)` - equals 1 if random < 0.7 (70% probability), else 0
- Multiplication `*` acts as logical AND
- Addition/subtraction modifies the type

**How it works for each atom type:**

**Type 1 (body atoms):**
- `type = 1`
- `(type==2) = 0` → first term is 0
- `(type==3) = 0` → second term is 0
- Result: `1 + 0 - 0 = 1` ✓ (stays type 1)

**Type 2 (patch A):**
- `type = 2`
- `(type==2) = 1`
- If coordinated AND random < 0.7: `1 * 1 * 1 = 1`
- Result: `2 + 1 - 0 = 3` ✓ (changes to type 3)
- If not coordinated OR random >= 0.7: `1 * 0 * 1 = 0` or `1 * 1 * 0 = 0`
- Result: `2 + 0 - 0 = 2` ✓ (stays type 2)

**Type 3 (patch B):**
- `type = 3`
- `(type==3) = 1`
- If coordinated AND random < 0.7: `1 * 1 * 1 = 1`
- Result: `3 + 0 - 1 = 2` ✓ (changes to type 2)
- If not coordinated OR random >= 0.7: `1 * 0 * 1 = 0` or `1 * 1 * 0 = 0`
- Result: `3 + 0 - 0 = 3` ✓ (stays type 3)

### 4. **The Critical Loop** (Lines 311-369)

The main simulation loop performs these steps:

#### Step A: Update Groups (Lines 325-326)
```lammps
group patches_A type 2
group patches_B type 3
```
- Redefines groups based on CURRENT types
- This is necessary because types changed in the previous iteration

#### Step B: Recreate Computes (Lines 321-330)
```lammps
uncompute cA
uncompute cB
compute cA all coord/atom cutoff 2.5 group patches_A
compute cB all coord/atom cutoff 2.5 group patches_B
```
- Must delete old computes before redefining (LAMMPS requirement)
- Recreate with updated groups so coordination is recalculated

#### Step C: Regenerate Random Values (Lines 333-334)
```lammps
variable rand_seed equal ${rand_seed}+1000
variable rand_val atom random(0,1,${rand_seed})
```
- Increment seed to get new random numbers
- Each atom gets a fresh random value for this iteration

#### Step D: Run Dynamics (Lines 341-342)
```lammps
run 5000
run 5000
```
- Split into two smaller runs for stability
- During this time, patches can move and coordinate

#### Step E: **CRITICAL - Unfix Rigid Bodies** (Line 345)
```lammps
unfix rigid_nvt
```
- **This is the key step!** Must remove rigid constraint before type changes
- Without this, LAMMPS will crash or lose atoms

#### Step F: Remap Atoms (Line 348)
```lammps
change_box all remap
```
- Safety measure: ensures all atoms are properly in the box
- Helps prevent image flag issues

#### Step G: Apply Type Changes (Line 355)
```lammps
set atom * type v_new_type
```
- Actually changes the atom types
- Uses the variable formula we defined earlier
- Works for all atoms simultaneously

#### Step H: Remap Again (Line 358)
```lammps
change_box all remap
```
- Another safety remap after type changes

#### Step I: **CRITICAL - Recreate Rigid Bodies** (Line 361)
```lammps
fix rigid_nvt all rigid/nvt molecule temp 1.0 1.0 0.01
```
- Recreates the rigid body constraint
- LAMMPS automatically detects molecules from the data file
- Each molecule (3 body atoms + 1 patch) becomes rigid again

#### Step J: Reinitialize Velocities (Line 364)
```lammps
velocity all scale 1.0
```
- Rescales velocities to maintain temperature
- Important because recreating the fix can affect velocity distribution

#### Step K: Run 0 Steps (Line 367)
```lammps
run 0
```
- Forces LAMMPS to recalculate forces and update internal state
- Ensures everything is synchronized after the fix recreation

## Why This Works

1. **Temporary Removal**: By unfixing before type changes, we avoid the rigid fix's internal state conflicts
2. **Molecule Preservation**: The `molecule` keyword in the fix means LAMMPS uses the molecule IDs from the data file, which don't change
3. **Atomic Formula**: The variable formula ensures integer results and handles all atom types correctly
4. **Coordination-Based**: Only patches that are actually close enough to interact can change state
5. **Probabilistic**: The random component makes state changes stochastic, not deterministic

## Evolution of the Solution

### Early Attempts (That Failed):
1. **Direct type changes**: Changed types while rigid fix was active → "Lost atoms" error
2. **Property/atom approach**: Tried storing state in properties instead of types → too complex
3. **Restart file approach**: Write restart, modify data, restart → works but inefficient

### Final Solution:
- **Unfix/Refix cycle**: Simple, reliable, works within a single simulation
- **Single variable formula**: One formula handles all atom types elegantly
- **Coordination detection**: Only changes state when patches are actually interacting

## Key Insights

1. **Rigid fixes are stateful**: They cache atom-to-molecule mappings, so they must be recreated after type changes
2. **Molecule IDs are persistent**: Even when types change, molecule IDs from the data file remain, so LAMMPS can reconstruct rigid bodies
3. **Variable formulas must be integer**: The `set atom` command requires integer types, so the formula must always evaluate to integers
4. **Coordination cutoff matters**: Must match the interaction cutoff to detect when patches can actually interact
5. **Random seed increment**: Incrementing the seed each iteration ensures we get new random values, not the same ones

## Performance Considerations

- The unfix/refix cycle adds overhead, but it's necessary for correctness
- Running in smaller chunks (e.g., 5000 steps) helps stability
- The `run 0` after refixing ensures proper initialization without wasting simulation time

This mechanism allows for **reversible, coordination-dependent, probabilistic state changes** in a rigid body system, which is exactly what we need for simulating patchy monomers with dynamic patch states!

