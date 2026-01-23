# Evolution of State Change Logic: Why Random Changes Were Happening

## Your Goal (From the Start)
**"Molecules should ONLY change when patches form bonds with same-type patches (conflict resolution)"**

This is a clear, unambiguous requirement:
- ✅ Change when: Type 1 patch touches another Type 1 patch → conflict → change
- ✅ Change when: Type 3 patch touches another Type 3 patch → conflict → change
- ❌ Don't change when: Type 1 patch touches Type 3 patch → no conflict → stay unchanged

---

## What Was Actually Implemented (Before Today's Fix)

### **Logic B (Type 3/4/5) - CORRECT ✅**
```cpp
// Logic B: Reacted monomers
if ((my_eff_type == 3 || 4 || 5) && confident && trigger_moment) {
    bool touching_same = false;
    if (my_eff_type == 3 && coord_3 >= 1.0) touching_same = true;
    if (my_eff_type == 4 && coord_4 >= 1.0) touching_same = true;
    if (my_eff_type == 5 && coord_5 >= 1.0) touching_same = true;
    
    if (touching_same) {  // ✅ Only changes when touching same-type!
        // ... change logic ...
    }
}
```
**This was CORRECT from the start** - it only triggers when touching same-type patches.

### **Logic A (Type 1) - INCORRECT ❌ (Before Today's Fix)**
```cpp
// Logic A: Unreacted monomers (BEFORE FIX)
if (my_eff_type == 1 && confident && trigger_moment) {
    // ❌ MISSING: No check for same-type contact!
    // ❌ Would change when touching ANY type (Type 1, 3, 4, or 5)
    // ... change logic ...
}
```
**This was WRONG** - it was missing the same-type check that Logic B had!

---

## Why This Happened (The Confusion)

Looking at the code comments and evolution, there appears to have been a conceptual mismatch:

### **Original Intent (Likely)**
The comment said "Confident Contact" which probably meant:
- Type 1 molecules form bonds → become "confident" about their state → change to a specific type

But "confident contact" was interpreted as:
- "Has been attached to ANY patch for 75+ steps" (not "attached to same-type patch")

### **What Should Have Been**
"Confident Contact" should have meant:
- "Has been attached to a SAME-TYPE patch for 75+ steps" (conflict situation)

---

## The Two-Logic Problem

The code had **two separate logics** for state changes:

1. **Logic A**: Type 1 → Type 3/4/5 (initial reaction)
   - **Original**: Changed when touching ANY type
   - **Your Goal**: Should change ONLY when touching Type 1 (same-type conflict)

2. **Logic B**: Type 3/4/5 → Different type (conflict resolution)
   - **Always**: Changed ONLY when touching same-type (correct from start)

**Why Two Logics?**
- Logic A was meant to represent "initial activation" - Type 1 molecules becoming "active"
- Logic B was meant to represent "conflict resolution" - active molecules avoiding same-type conflicts

**The Misunderstanding:**
- Logic A was implemented as "activation upon contact with anything"
- But your goal was "activation upon conflict with same-type"

---

## Why It Seemed Like "Random" Changes

Even though Logic A required `confident` (75+ steps of contact), it was:
1. **Too permissive**: Type 1 touching Type 3 → changes (shouldn't!)
2. **Triggered deterministically**: Every 50 steps after threshold → synchronized
3. **No conflict requirement**: Changed without checking for same-type conflict

So from your perspective, molecules were changing "randomly" because:
- They changed when attached to DIFFERENT types (shouldn't happen)
- Many changed together (synchronized due to deterministic triggers)
- No clear reason for the change (no conflict to resolve)

---

## Today's Fix

### **Logic A (Type 1) - NOW CORRECT ✅**
```cpp
// Logic A: Unreacted monomers (AFTER FIX)
bool touching_same_type_1 = (my_eff_type == 1 && coord_1 >= 1.0);  // ✅ NEW: Check for same-type!
if (my_eff_type == 1 && confident && trigger_moment && touching_same_type_1) {
    // ✅ NOW: Only changes when touching Type 1 (same-type conflict)
    // ... change logic ...
}
```

Now **both Logic A and Logic B** follow your requirement:
- **Only change when patches form bonds with same-type patches**

---

## Summary

| Aspect | Logic A (Before) | Logic A (After) | Logic B (Always) |
|--------|------------------|-----------------|------------------|
| **Triggers on** | Any contact | Same-type only ✅ | Same-type only ✅ |
| **Matches your goal?** | ❌ No | ✅ Yes | ✅ Yes |
| **Conflict required?** | ❌ No | ✅ Yes | ✅ Yes |

**The Root Cause:**
- Logic A was missing the same-type check that Logic B always had
- This made Type 1 molecules change when they shouldn't (attached to different types)
- Combined with synchronized triggering, this created the appearance of "random" mass changes

**Now:**
- Both logics require same-type conflict
- Changes only happen to resolve conflicts
- Matches your original goal perfectly! ✅




