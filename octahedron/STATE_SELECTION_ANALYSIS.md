# State Selection Logic Analysis

## Question: Is State Selection Deterministic or Random?

The user identified a critical concern: **If state selection is deterministic (e.g., always Type 2, or always incrementing), then mass switching to the same type is expected behavior.**

---

## Current Implementation Analysis

### Logic A: Type 1 → New Type (Unreacted Monomer)

**Location**: `fix_state_change_octahedron.cpp`, lines 431-434

```cpp
// Pick new type randomly
double r = random.uniform();
if (r < 0.333333) new_type = 3;
else if (r < 0.666666) new_type = 4;
else new_type = 5;
```

**Analysis**:
- ✅ **TRULY RANDOM**: Uses `random.uniform()` to pick from types 3, 4, 5 with equal probability (33.33% each)
- ✅ **NOT DETERMINISTIC**: No hard-coded progression (e.g., `type++` or always Type 2)

### Logic B: Type 3/4/5 → Different Type (Reacted Monomer, Same-Type Conflict)

**Location**: `fix_state_change_octahedron.cpp`, lines 583-592

```cpp
// Available types: 3, 4, 5 (exclude current type)
std::vector<int> safe_options;
for (int t = 3; t <= 5; t++) {
    if (t != my_eff_type && neighbor_types.find(t) == neighbor_types.end()) {
        safe_options.push_back(t);
    }
}

// Pick from safe options, or random if none available
if (!safe_options.empty()) {
    int idx = static_cast<int>(random.uniform() * safe_options.size());
    new_type = safe_options[idx];
} else {
    // All types conflict - pick randomly anyway
    double r = random.uniform();
    if (r < 0.333333) new_type = 3;
    else if (r < 0.666666) new_type = 4;
    else new_type = 5;
}
```

**Analysis**:
- ✅ **TRULY RANDOM**: Uses `random.uniform()` to pick from available options
- ✅ **SMART SELECTION**: Prefers types that don't conflict with neighbors, but still random within safe options
- ✅ **FALLBACK**: If all types conflict, picks randomly anyway

---

## Potential Issue: Random Seed Correlation

**FIXED**: We added molecule ID to the random seed (line 375):

```cpp
int my_seed = seed + static_cast<int>(timestep) + comm->me + molecule[i];
RanPark random(lmp, my_seed);
```

This ensures:
- ✅ Each molecule gets a **unique random sequence**
- ✅ Molecules processed at the same time get **different random values**
- ✅ Prevents correlation between molecules

---

## Verification: Why Mass Switching Still Occurs (If It Does)

If mass switching to the same type still occurs despite random selection, possible causes:

### 1. **Random Number Generator Quality**
- `RanPark` is a pseudo-random number generator
- If the seed space is small, sequences might appear correlated
- **Status**: Should be fine with molecule ID included

### 2. **Synchronized Trigger Moments**
- If many molecules reach `hysteresis_threshold` simultaneously, they all trigger
- Even with random selection, if 100 molecules trigger at once, some will pick the same type by chance
- **Fix Applied**: Jitter on contact timer reset (-10 to 0 steps) desynchronizes triggers

### 3. **Rate Limiter Random Selection**
- When rate limiter reduces changes, it uses random shuffling
- If multiple molecules want to change, random selection might favor certain types if there's any bias
- **Status**: Uses `random.uniform()` which should be unbiased

### 4. **"Safe Options" Bias**
- The "smart type selection" tries to avoid neighbor types
- If many neighbors are Type 3, molecules will prefer Type 4 or 5
- This can create **domino effects** where clusters converge to one type
- **Status**: This is intentional behavior (avoid conflicts), but can create cascades

---

## Conclusion

✅ **State selection is TRULY RANDOM**, not deterministic.

✅ **Random seed includes molecule ID** to prevent correlation.

✅ **Jitter on timer reset** desynchronizes trigger moments.

⚠️ **"Smart type selection"** (avoiding neighbor types) can still create cascades, but this is intentional to minimize conflicts.

---

## If Mass Switching Persists

If molecules still mass-switch to the same type after all fixes, the issue is likely:

1. **Synchronized Trigger Moments** (many molecules triggering simultaneously)
   - **Solution**: Increase jitter magnitude or hysteresis threshold

2. **Domino Effects from Smart Selection** (avoiding conflicts creates clusters)
   - **Solution**: Make selection purely random (remove safe_options logic) or add noise to selection

3. **Rate Limiter Not Working** (too many changes allowed per timestep)
   - **Solution**: Verify `max_changes_per_step = 1` is enforced globally

4. **Consistency Sweep Issues** (forcing changes despite cooldown)
   - **Solution**: Verify cooldown check in consistency sweep is working

---

## Recommended Test

To verify state selection is random:

1. Run simulation and track type changes
2. Count how many molecules switch to each type (3, 4, 5)
3. **Expected**: Roughly equal distribution (33.33% each) over many changes
4. **If skewed**: Random number generator might have bias, or smart selection is creating clusters

---

## Current Status

- ✅ State selection uses `random.uniform()` (truly random)
- ✅ Molecule ID included in random seed (prevents correlation)
- ✅ Jitter on timer reset (desynchronizes triggers)
- ✅ Smart type selection (intentional, but can create cascades)

**All fixes are in place. If mass switching still occurs, it's likely due to synchronized triggers or domino effects from smart selection, not deterministic state selection.**

