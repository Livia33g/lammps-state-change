# Consistency Sweep Diagnostics - Monitoring Guide

## Overview

Diagnostic counters have been added to monitor the consistency sweep behavior and verify it's not causing mass synchronized changes.

## Debug Output Format

Every 10,000 steps, the simulation prints debug information:

```
DEBUG[step X]: Confident=N, TriggerAttempts=M, CooldownBlocked=K, Changes=L | SWEEP: Attempts=A, Blocked(Cooldown=C,Timestamp=T), Applied=P
```

### Main Counters:
- **Confident**: Number of patches that reached hysteresis threshold (contact_timer > threshold)
- **TriggerAttempts**: Number of molecules that attempted to trigger a state change
- **CooldownBlocked**: Number of changes blocked by cooldown in main logic
- **Changes**: Number of molecules that actually changed state (via main logic)

### Sweep Counters:
- **Attempts (A)**: Number of molecules the sweep tried to update
- **Blocked Cooldown (C)**: Number blocked by Gate 1 (molecule in cooldown)
- **Blocked Timestamp (T)**: Number blocked by Gate 2 (timestamp difference too small or future timestamp)
- **Applied (P)**: Number actually updated by the sweep

## Interpreting the Output

### Healthy Behavior (No Mass Switching):

```
DEBUG[step 100000]: Confident=2, TriggerAttempts=1, CooldownBlocked=5, Changes=1 | SWEEP: Attempts=0, Blocked(Cooldown=0,Timestamp=0), Applied=0
```

**Interpretation**:
- Main logic: 1 change (rate limited correctly)
- Sweep: 0 attempts (no molecules need syncing)
- **This is GOOD** - no mass changes

### Problem Behavior (Mass Switching via Sweep):

```
DEBUG[step 100000]: Confident=1, TriggerAttempts=1, CooldownBlocked=0, Changes=1 | SWEEP: Attempts=50, Blocked(Cooldown=5,Timestamp=10), Applied=35
```

**Interpretation**:
- Main logic: 1 change (rate limited correctly)
- Sweep: Tried to update 50 molecules, applied 35 updates
- **This is BAD** - sweep is causing mass changes!

### Expected After Fix:

```
DEBUG[step 100000]: Confident=1, TriggerAttempts=1, CooldownBlocked=0, Changes=1 | SWEEP: Attempts=5, Blocked(Cooldown=3,Timestamp=2), Applied=0
```

**Interpretation**:
- Main logic: 1 change (rate limited correctly)
- Sweep: Tried to update 5 molecules, but all were blocked by gates
- **This is GOOD** - sweep is conservative, not causing cascades

## Diagnostic Counters Location

### Header File:
`fix_state_change_octahedron.h`, lines ~75-78:
```cpp
int nsweep_attempts;          // Count of molecules sweep tries to update
int nsweep_blocked_cooldown;  // Count blocked by Gate 1 (cooldown)
int nsweep_blocked_timestamp; // Count blocked by Gate 2 (timestamp difference too small)
int nsweep_applied;           // Count actually updated by sweep
```

### Implementation:
- Counters incremented in consistency sweep loop (lines ~767-815)
- Counters output in debug section (lines ~815-830)
- Counters reset after output (line ~833)

## What to Look For

### Red Flags (Mass Switching Indicators):

1. **High Sweep Applied Count**:
   - If `Applied > 10` per 10,000 steps, the sweep is forcing too many updates
   - This indicates the gates are not working correctly

2. **Sweep Applied > Main Changes**:
   - If sweep applies more changes than main logic, it's bypassing the rate limiter
   - This is the "backdoor" problem we're trying to fix

3. **Low Blocked Counts**:
   - If `Blocked(Cooldown=X,Timestamp=Y)` is very low but `Attempts` is high
   - This means the gates are not blocking effectively

### Green Flags (Healthy Behavior):

1. **Sweep Applied = 0 or Very Low**:
   - Most molecules are properly synced across processors
   - Only occasional updates needed

2. **High Blocked Counts**:
   - Gates are working correctly
   - Cooldown and timestamp checks are preventing cascades

3. **Sweep Applied ≤ Main Changes**:
   - Sweep is not bypassing the rate limiter
   - Main logic is controlling state changes

## Troubleshooting

### If Sweep Applied is High:

1. **Check Gate 1 (Cooldown)**:
   - If `Blocked(Cooldown)` is low, cooldown checks might not be working
   - Verify `cooldown_duration[i]` is being set correctly

2. **Check Gate 2 (Timestamp)**:
   - If `Blocked(Timestamp)` is low, timestamp threshold (100 steps) might be too small
   - Consider increasing the threshold to 500 or 1000 steps

3. **Check Gate 3 (Current Timestep)**:
   - Verify `last_change[i] = timestep` is being set (not target_timestamp)
   - This prevents backdating

### If Sweep Attempts is High:

- This might be normal if many molecules are split across processors
- Focus on the `Applied` count - that's what actually causes changes
- High `Attempts` with low `Applied` means gates are working

## Next Steps

1. Run a test simulation and monitor the debug output
2. Check if `Sweep Applied` is consistently 0 or very low
3. If mass switching occurs, check the sweep counters to see if they're high
4. Adjust gate thresholds if needed (cooldown duration, timestamp threshold)

## Summary

The diagnostic counters allow you to verify:
- ✅ Is the consistency sweep causing mass changes? (Check `Applied` count)
- ✅ Are the gates working? (Check `Blocked` counts)
- ✅ Is the sweep bypassing the rate limiter? (Compare `Applied` vs `Changes`)

**Goal**: `Sweep Applied` should be 0 or very low (< 2 per 10,000 steps), indicating the gates are preventing cascades.

