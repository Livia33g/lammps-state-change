# Mass Behavior Investigation

## Status
✅ **You're NOT in trouble!** The login node policy just automatically kills long-running processes. It's a common mistake and completely normal. The system just enforces the policy automatically.

## Analysis Results

From quick analysis of the dump file:
- **Steps 500K-1500K**: 86% Type 3 (43/50 molecules) - ⚠️ **Mass behavior detected**
- **Step 2M**: 52% max distribution (T3:23, T5:26) - ✅ **Good distribution**

## Potential Root Causes

### 1. Trigger Moment Synchronization (Likely Issue)
**Location**: `fix_state_change_octahedron.cpp`, lines 421-422

The current trigger logic allows molecules to trigger every 50 steps after threshold:
```cpp
bool trigger_moment = (contact_timer[i] == (hysteresis_threshold + 1)) || 
                      (contact_timer[i] > hysteresis_threshold && (contact_timer[i] - hysteresis_threshold) % 50 == 0);
```

**Problem**: If many molecules reach the threshold around the same time (e.g., during initial cluster formation), they will all trigger at timesteps that are multiples of 50 (e.g., step 500030, 500080, 500130...). Even with the rate limiter (max 1 per timestep), if 50 molecules all want to change at step 500030, they'll queue up and change sequentially over the next 50 timesteps, creating a cascade.

**Why Type 3 Dominates**: Once type 3 becomes dominant (even slightly), the "smart selection" logic might favor type 3 because it's more likely to avoid conflicts with neighbors (who are also type 3). This creates a positive feedback loop.

### 2. Rate Limiter May Not Prevent Cascades
**Location**: `fix_state_change_octahedron.cpp`, lines 627-660

The rate limiter limits to **1 change per timestep globally**, but if 50 molecules all trigger at the same timestep, they'll be processed sequentially. Even though only 1 changes per timestep, over 50 timesteps, 50 molecules can change. If they all choose the same type (type 3) due to smart selection, this creates mass behavior.

### 3. Smart Selection May Create Positive Feedback
**Location**: `fix_state_change_octahedron.cpp`, lines 588-600

The "smart selection" logic picks a type that avoids conflicts with neighbors. If type 3 is already dominant:
- Most neighbors are type 3
- Smart selection avoids type 3 (to avoid conflicts)
- BUT: if ALL neighbors are type 3, then types 4 and 5 also avoid conflicts
- Random selection between 4 and 5 might still favor 3 in some edge cases
- OR: type 3 might be selected because it's "safe" if most neighbors are 3

## Recommended Fixes

### Fix 1: Make Trigger Moments More Random (Priority 1)
Instead of triggering every 50 steps (deterministic), add more randomness:

```cpp
// Current (allows synchronized triggering):
bool trigger_moment = (contact_timer[i] == (hysteresis_threshold + 1)) || 
                      (contact_timer[i] > hysteresis_threshold && (contact_timer[i] - hysteresis_threshold) % 50 == 0);

// Proposed (more random):
bool trigger_moment = (contact_timer[i] == (hysteresis_threshold + 1)) || 
                      (contact_timer[i] > hysteresis_threshold && 
                       random.uniform() < (1.0 / 50.0));  // 2% chance per step after threshold
```

This desynchronizes triggers: molecules won't all trigger at the same deterministic timesteps.

### Fix 2: Increase Rate Limiter Strictness
Instead of allowing 1 change per timestep, consider a sliding window:
- Max 5 changes per 1000 timesteps
- This prevents cascades even if many molecules want to change

### Fix 3: Add Anti-Dominance Bias
Modify smart selection to actively avoid types that are already dominant:
- Count how many neighbors are each type
- Bias selection away from the most common neighbor type
- This breaks positive feedback loops

## Next Steps

1. ✅ Created `analyze_dump.slurm` to analyze the full dump file on a compute node (submitted)
2. ⏳ Wait for analysis results to confirm the pattern
3. 🔧 Implement Fix 1 (randomize trigger moments) if analysis confirms synchronized triggering
4. 🧪 Test with new simulation run

## Current Status

The simulation completed successfully. Final state shows good distribution (52% max), suggesting the system eventually corrects itself, but the intermediate mass behavior (86% type 3) indicates the anti-synchronization mechanisms need strengthening.

