# Diagnostic Counters Added for Consistency Sweep

## Summary

Diagnostic counters have been successfully added to monitor the consistency sweep behavior and verify it's not causing mass synchronized changes.

## What Was Added

### 1. New Counter Variables (Header File)

Added to `fix_state_change_octahedron.h`:
```cpp
int nsweep_attempts;          // Count of molecules sweep tries to update
int nsweep_blocked_cooldown;  // Count blocked by Gate 1 (cooldown)
int nsweep_blocked_timestamp; // Count blocked by Gate 2 (timestamp difference too small)
int nsweep_applied;           // Count actually updated by sweep
```

### 2. Counter Incrementation (Implementation)

In the consistency sweep loop (`fix_state_change_octahedron.cpp`, lines ~760-815):
- `nsweep_attempts++` when sweep considers a molecule for update
- `nsweep_blocked_cooldown++` when Gate 1 (cooldown) blocks an update
- `nsweep_blocked_timestamp++` when Gate 2 (timestamp) blocks an update
- `nsweep_applied++` when an update is actually applied

### 3. Debug Output Enhancement

Debug output now includes sweep statistics:
```
DEBUG[step X]: Confident=N, TriggerAttempts=M, CooldownBlocked=K, Changes=L | SWEEP: Attempts=A, Blocked(Cooldown=C,Timestamp=T), Applied=P
```

### 4. Counter Initialization and Reset

- Counters initialized in constructor (line ~96-100)
- Counters reset after debug output (line ~833-836)
- Counters also initialized before consistency sweep (line ~759-762)

## How to Use

### Monitor During Simulation

The debug output appears every 10,000 steps in `stderr` (or `lammps_stdout.log`). Look for lines like:

```
DEBUG[step 100000]: Confident=2, TriggerAttempts=1, CooldownBlocked=5, Changes=1 | SWEEP: Attempts=0, Blocked(Cooldown=0,Timestamp=0), Applied=0
```

### Key Metrics to Watch

1. **Sweep Applied Count**: Should be 0 or very low (< 2 per 10,000 steps)
   - If high (> 10), the sweep is causing mass changes

2. **Sweep Blocked Counts**: Should be high relative to attempts
   - Indicates gates are working correctly

3. **Comparison**: Sweep Applied vs Main Changes
   - Sweep Applied should NOT exceed Main Changes
   - If it does, sweep is bypassing the rate limiter

## Status

✅ **Code compiled successfully**
✅ **All counters implemented**
✅ **Debug output enhanced**
✅ **Documentation created**: `CONSISTENCY_SWEEP_DIAGNOSTICS.md`

## Next Steps

1. Run a test simulation
2. Monitor the debug output for sweep statistics
3. Verify `Sweep Applied` is consistently low
4. If mass switching occurs, check if `Sweep Applied` is high - this would confirm the consistency sweep is the cause

## Files Modified

1. `fix_state_change_octahedron.h` - Added counter declarations
2. `fix_state_change_octahedron.cpp` - Added counter incrementation and output
3. `CONSISTENCY_SWEEP_DIAGNOSTICS.md` - Created monitoring guide

## Expected Behavior

After the fixes, you should see:
- **Sweep Applied = 0 or very low** (gates are blocking updates)
- **Sweep Blocked counts > 0** (gates are working)
- **Main Changes ≤ 1 per timestep** (rate limiter working)
- **No mass synchronized changes** (all mechanisms working together)

