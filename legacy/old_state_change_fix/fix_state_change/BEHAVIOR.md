# State Change Fix Behavior - ONE-SHOT MECHANISM

## Overview
The `fix state/change` implements reversible state changes for patchy monomers with a **one-shot mechanism**. Patches get exactly ONE chance to switch when they attach. If neither switches, they stay fixed while attached.

## Key Behavior: One-Shot on Attachment

### When Patches Attach
1. **Two patches of the same type come close** (e.g., two type 2 patches)
2. **New attachment detected**: Coordination goes from < 0.1 to > 0.1
3. **ONE chance to switch**: Each patch gets a probability-based chance (default: 70%)
4. **If probability triggers**: Patch switches to opposite type
5. **If probability doesn't trigger**: Patch stays the same type
6. **Shot is consumed**: Whether they switch or not, they've had their shot
7. **No more checks while attached**: They won't check again until they detach

### While Attached
- **If one switched**: They become different types (2 and 3), coordination drops, they detach
- **If neither switched**: They stay the same type, stay attached, but **won't check again**
- **State is fixed**: Their types remain unchanged while attached

### When Patches Detach
- **Natural detachment**: Bond breaks due to other interactions (coordination drops to < 0.1)
- **Shot status resets**: `last_change[i] = -1` (ready for a shot)
- **Become viable again**: If they reattach later, they get another shot

### Reattachment
- **New attachment detected**: Coordination goes from < 0.1 to > 0.1 again
- **Get another shot**: Since `last_change[i] == -1`, they're ready
- **Same process**: One chance to switch based on probability

## Example Timeline

```
Timestep 0:    Two type 2 patches are separate (not attached)
               → last_change[i] = -1 (ready for shot)

Timestep 100:  Patches attach (coordination > 0.1)
               → New attachment detected!
               → Patch A: 70% chance → switches to type 3 ✓
               → Patch B: 70% chance → doesn't switch (30% chance)
               → Both mark last_change[i] = timestep (shot consumed)
               → Patch A is now type 3, Patch B is type 2
               → Different types → coordination drops → they detach

Timestep 200:  Patches detach (coordination < 0.1)
               → last_change[i] = -1 (reset, ready for next shot)

Timestep 300:  Patches reattach (both are type 2 again somehow, or new pair)
               → New attachment detected!
               → Get another shot (70% probability each)
               → Process repeats...
```

## Edge Case: Forever Attached Without Switching

**Scenario**: Two patches attach, neither switches, and they stay attached forever

- **Initial attachment**: Both get their shot, neither switches
- **Shot consumed**: `last_change[i] = timestep` (marked as had shot)
- **Stay attached**: They remain the same type (both type 2 or both type 3)
- **No more checks**: They won't check again while attached
- **State is fixed**: They stay attached with same types until something external breaks the bond

## Key Points

1. **One shot per attachment**: Each attachment gives exactly one chance
2. **Shot is consumed regardless**: Whether they switch or not, shot is used
3. **Detachment resets**: Only when they detach do they become viable again
4. **No time-based cooldown**: The "cooldown" parameter is now unused (kept for compatibility)
5. **Probability is per-shot**: Each shot is independent (70% chance)

## Parameters

- `check_every`: How often to check for new attachments (default: 100 steps)
- `cooldown_steps`: **UNUSED** (kept for compatibility, but one-shot mechanism doesn't use it)
- `probability`: Probability of switching when shot is given (default: 0.7 = 70%)
- `cutoff`: Coordination cutoff distance (default: 2.5)

## Summary

✅ **One shot per attachment** - Patches get exactly one chance when they attach  
✅ **Shot consumed regardless** - Whether they switch or not, shot is used  
✅ **Fixed while attached** - If neither switches, they stay fixed until detachment  
✅ **Detachment resets** - Only when they detach do they become viable again  
✅ **Reattachment gives new shot** - Each new attachment is a fresh opportunity
