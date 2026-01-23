# Diagnostic Analysis: Mass Synchronized State Changes

## Problem Statement
Despite multiple layers of anti-synchronization mechanisms, monomers continue to mass-switch to the same color simultaneously. This document identifies potential root causes for external debugging assistance.

## Current Anti-Synchronization Mechanisms Implemented

1. **Hysteresis/Contact Timer**: Requires `hysteresis_threshold` (75 steps) of consecutive contact before triggering
2. **Symmetry Breaking**: Molecule ID-based tie-breaker (lower ID changes)
3. **Cooldown System**: Randomized cooldown durations (0.3x to 2.0x of base `cooldown_steps`)
4. **Probabilistic Activation**: Reduced probabilities (5-10% of base probability)
5. **Global Rate Limiter**: Maximum 1 molecule change per timestep globally (via `MPI_Allreduce`)
6. **Consistency Sweep**: Synchronizes `effective_type` across MPI boundaries using newest `last_change` timestamp

## Critical Issues Identified

### Issue 1: Consistency Sweep Timing and Cascade Risk
**Location**: `fix_state_change_octahedron.cpp`, lines 687-737

**Problem**: The consistency sweep runs AFTER state changes are applied and communicated. It forces all patches in a molecule to match the "newest" `effective_type` found across all processors (including ghosts).

**Cascade Scenario**:
1. Processor A: Molecule 10 changes to type 3 at timestep T
2. Processor B: Molecule 10 (ghost) still has type 4 with `last_change = T-1000`
3. Consistency sweep sees Processor A's change is "newer" (T > T-1000)
4. Processor B forces molecule 10 to type 3
5. If molecule 10 was attached to molecule 11 (also type 3), this creates a NEW conflict
6. Next timestep: Molecule 11 sees conflict and changes, triggering another cascade

**Why This Causes Mass Changes**:
- The consistency sweep can create NEW conflicts by forcing state changes
- These forced changes bypass all cooldown/probability checks
- The sweep doesn't respect the rate limiter (it happens after rate limiting)

**Evidence**: The consistency sweep modifies `effective_type` and `last_change` without incrementing `nchanges`, so it's invisible to the rate limiter.

---

### Issue 2: Rate Limiter Only Limits Count, Not Type Distribution
**Location**: `fix_state_change_octahedron.cpp`, lines 601-638

**Problem**: The rate limiter ensures ≤1 molecule changes per timestep globally, but it doesn't prevent multiple molecules from changing to the SAME type over consecutive timesteps.

**Cascade Scenario**:
- Timestep T: Molecule 5 changes to type 3 (allowed, rate limiter = 1)
- Timestep T+1: Molecule 6 changes to type 3 (allowed, rate limiter = 1)
- Timestep T+2: Molecule 7 changes to type 3 (allowed, rate limiter = 1)
- Result: Over 3 timesteps, 3 molecules all become type 3, creating a cluster

**Why This Causes Mass Changes**:
- The rate limiter prevents simultaneous changes but allows sequential cascades
- If many molecules are "ready" to change (past hysteresis threshold), they'll change one per timestep in sequence
- All changing to the same type creates visual "mass change" effect

**Evidence**: The rate limiter only checks `mol_changes.size()`, not the distribution of `new_type` values.

---

### Issue 3: Hysteresis Trigger Moment Synchronization
**Location**: `fix_state_change_octahedron.cpp`, lines 406-407, 419, 435

**Problem**: The trigger condition `contact_timer[i] == (hysteresis_threshold + 1)` means ALL molecules that reach the threshold at the same time will trigger on the SAME timestep.

**Cascade Scenario**:
1. Many molecules attach at timestep T-75
2. All `contact_timer` values increment together: T-74, T-73, ..., T-1, T
3. At timestep T, all molecules have `contact_timer = 75` (threshold = 75)
4. At timestep T+1, ALL molecules have `contact_timer = 76` (threshold + 1)
5. ALL molecules trigger simultaneously, rate limiter reduces to 1, but they all "want" to change

**Why This Causes Mass Changes**:
- Even with rate limiting, if 50 molecules all trigger at once, the random selection might pick molecules that all choose the same type
- The probabilistic activation (5-10%) doesn't help if they all pass the probability check

**Evidence**: The `trigger_moment` check creates a "synchronization point" where many molecules are ready simultaneously.

---

### Issue 4: Random Seed Determinism Across Processors
**Location**: `fix_state_change_octahedron.cpp`, line 361

**Problem**: The random seed is `seed + timestep + comm->me`. While `comm->me` differs per processor, if the logic is deterministic, processors might still make similar decisions.

**Cascade Scenario**:
- Processor 0: Molecule 10 (ID=10) vs Molecule 20 (ID=20) → 10 changes (lower ID)
- Processor 1: Molecule 11 (ID=11) vs Molecule 21 (ID=21) → 11 changes (lower ID)
- Both processors use the same random sequence (just offset by `comm->me`)
- If molecule IDs are clustered (low IDs on proc 0, high IDs on proc 1), all low-ID molecules change

**Why This Causes Mass Changes**:
- The ID-based tie-breaker is deterministic, not random
- If molecule IDs are non-uniformly distributed, many molecules might have "winning" IDs
- The probabilistic activation (5-10%) might not be enough to break the pattern

**Evidence**: The random seed doesn't include molecule ID, so molecules with similar IDs might get similar random values.

---

### Issue 5: Ghost Atom State Staleness
**Location**: `fix_state_change_octahedron.cpp`, lines 466-515 (neighbor checking)

**Problem**: When checking neighbors for conflicts, the code uses `effective_type[j]` from ghost atoms. If ghost atoms have stale `effective_type` values (not yet communicated), incorrect conflict detection occurs.

**Cascade Scenario**:
1. Processor A: Molecule 10 changes to type 3 at timestep T
2. Processor B: Molecule 10 (ghost) still has `effective_type = 4` (stale)
3. Processor B: Molecule 11 sees molecule 10 (ghost) as type 4, no conflict detected
4. Processor B: Molecule 11 changes to type 4
5. Next communication: Processor A sees molecule 11 as type 4, but molecule 10 is type 3 → conflict created

**Why This Causes Mass Changes**:
- Stale ghost data causes incorrect conflict detection
- Molecules change based on outdated information
- Creates cascading conflicts as communication catches up

**Evidence**: The neighbor check uses `nall = atom->nlocal + atom->nghost`, but ghost `effective_type` might be stale if `comm->forward_comm(this)` hasn't been called recently.

---

### Issue 6: Contact Timer Reset After State Change
**Location**: `fix_state_change_octahedron.cpp`, lines 430, 585, 661

**Problem**: After a state change, `contact_timer[i] = 0` is reset. If many molecules change simultaneously, they all reset their timers together, potentially synchronizing future triggers.

**Cascade Scenario**:
1. Timestep T: 10 molecules all change (rate limited to 1, but consistency sweep forces 9 more)
2. All 10 molecules reset `contact_timer = 0`
3. All 10 molecules are still attached (didn't move)
4. All 10 molecules increment `contact_timer` together: 1, 2, 3, ..., 75, 76
5. All 10 molecules trigger again at timestep T+76

**Why This Causes Mass Changes**:
- Resetting timers synchronizes future triggers
- Even with randomized cooldowns, the contact timers are synchronized
- Creates periodic "mass change" events every `hysteresis_threshold + 1` steps

**Evidence**: The contact timer reset happens immediately after state change, without jitter.

---

### Issue 7: Cooldown Duration Assignment Timing
**Location**: `fix_state_change_octahedron.cpp`, lines 643-651

**Problem**: Cooldown durations are randomized per molecule, but if many molecules change at once, they might get similar random values if the random number generator state is similar.

**Cascade Scenario**:
1. Timestep T: 5 molecules change (rate limited, but consistency sweep forces more)
2. All 5 molecules get cooldown durations from the same random sequence
3. If the random sequence is correlated, cooldown durations might be similar
4. All 5 molecules "wake up" from cooldown at similar times
5. All 5 molecules trigger together again

**Why This Causes Mass Changes**:
- The random seed is `seed + timestep + comm->me`, which is deterministic
- Molecules changing at the same timestep get similar random values
- Cooldown durations might cluster, causing synchronized "re-awakening"

**Evidence**: The cooldown randomization uses `random.uniform()` from the same `RanPark` instance, which might have correlated sequences.

---

### Issue 8: Consistency Sweep Doesn't Respect Cooldown
**Location**: `fix_state_change_octahedron.cpp`, lines 719-734

**Problem**: The consistency sweep forces `effective_type` changes without checking cooldown. A molecule that just changed (and is in cooldown) can be forced to change again by the consistency sweep.

**Cascade Scenario**:
1. Processor A: Molecule 10 changes to type 3 at timestep T (enters cooldown)
2. Processor B: Molecule 10 (ghost) has type 4 with `last_change = T-1000`
3. Consistency sweep: Processor A's change (T) is newer than Processor B's (T-1000)
4. Processor B forces molecule 10 to type 3, even though it's in cooldown
5. This creates a conflict with molecule 11 (type 3), triggering another change

**Why This Causes Mass Changes**:
- Consistency sweep bypasses cooldown checks
- Forces state changes that shouldn't happen
- Creates cascading conflicts

**Evidence**: The consistency sweep doesn't check `cooldown_duration[i]` before forcing changes.

---

## Recommended Fixes (Priority Order)

### Priority 1: Fix Consistency Sweep to Respect Cooldown
**Action**: Add cooldown check before forcing changes in consistency sweep.
**Location**: Line 728, before `effective_type[i] = target_type;`
**Code**:
```cpp
// Check if this molecule is in cooldown before forcing change
bool in_cooldown = false;
if (last_change[i] >= 0) {
    bigint steps_since_change = timestep - last_change[i];
    if (steps_since_change < cooldown_duration[i]) {
        in_cooldown = true;
    }
}
if (!in_cooldown && effective_type[i] != target_type) {
    effective_type[i] = target_type;
    last_change[i] = newest_time_per_mol[mol_id];
    consistency_fixes++;
}
```

### Priority 2: Add Jitter to Contact Timer Reset
**Action**: Add random jitter when resetting contact timer after state change.
**Location**: Lines 430, 585, 661
**Code**:
```cpp
// Reset timer with jitter to prevent synchronization
contact_timer[i] = static_cast<int>(random.uniform() * 10);  // 0-9 steps jitter
```

### Priority 3: Include Molecule ID in Random Seed
**Action**: Add molecule ID to random seed to break correlation.
**Location**: Line 361
**Code**:
```cpp
int my_seed = seed + static_cast<int>(timestep) + comm->me + molecule[i];
```

### Priority 4: Rate Limiter Should Track Type Distribution
**Action**: Limit not just count, but also prevent too many changes to the same type.
**Location**: Lines 601-638
**Code**: Add a check to ensure no more than 1 molecule changes to the same type per timestep.

### Priority 5: Ensure Ghost State is Fresh Before Conflict Detection
**Action**: Call `comm->forward_comm(this)` before neighbor checking.
**Location**: Before line 466
**Code**: Add `comm->forward_comm(this);` to ensure ghost `effective_type` is up-to-date.

---

## Testing Recommendations

1. **Add Debug Output**: Log every state change with molecule ID, new type, timestep, and processor ID to identify patterns.
2. **Visualize Type Distribution**: Plot histogram of types over time to see if mass changes are sequential or simultaneous.
3. **Check Consistency Sweep Impact**: Count how many molecules are forced to change by consistency sweep vs. natural triggers.
4. **Monitor Rate Limiter**: Log `total_changes_global` and `max_changes_per_step` to verify rate limiting is working.
5. **Track Contact Timer Distribution**: Log distribution of `contact_timer` values to see if they're synchronized.

---

## Files to Review

1. `/work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/octahedron/fix_state_change_octahedron.cpp`
   - Lines 353-596: Main state change logic
   - Lines 601-638: Global rate limiter
   - Lines 687-737: Consistency sweep
   - Lines 430, 585, 661: Contact timer resets

2. `/work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/octahedron/fix_state_change_octahedron.h`
   - Member variables: `contact_timer`, `cooldown_duration`, `max_changes_per_step`

3. `/work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/octahedron/octahedron_simulation_cpp/dump.octahedron_monomers.lammpstrj`
   - Visual inspection: Look for timesteps where many molecules change to the same type

4. `/work/nvme/bewl/lguttieres/sims/self_processors/sim_templates/state_change/octahedron/octahedron_simulation_cpp/lammps_stdout.log`
   - Debug output: Check for "DEBUG[step X]" messages showing trigger attempts and changes

---

## Questions for External Help

1. **MPI Communication Timing**: Is `comm->forward_comm(this)` called frequently enough to keep ghost `effective_type` values fresh?
2. **Random Number Generation**: Is `RanPark` with `seed + timestep + comm->me` sufficiently random, or should we use a better RNG?
3. **Consistency Sweep Design**: Should the consistency sweep run BEFORE state changes (to prevent conflicts) or AFTER (to fix splits)? Current design runs AFTER, which might cause cascades.
4. **Rate Limiter Algorithm**: Is proportional reduction the right approach, or should we use a different algorithm (e.g., priority queue based on molecule ID)?
5. **Hysteresis Threshold**: Is 75 steps too low? Should we increase it or add jitter to the threshold itself?

---

## Summary

The most likely root cause is **Issue 1 (Consistency Sweep Timing)** combined with **Issue 8 (Consistency Sweep Doesn't Respect Cooldown)**. The consistency sweep can force state changes that bypass all anti-synchronization mechanisms, creating cascading conflicts that lead to mass synchronized changes.

The second most likely cause is **Issue 3 (Hysteresis Trigger Moment Synchronization)**, where many molecules reach the threshold simultaneously and all trigger together, overwhelming the rate limiter.

Recommended immediate action: Fix the consistency sweep to respect cooldown periods and add jitter to contact timer resets.

