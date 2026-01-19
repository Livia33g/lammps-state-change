# Instructions: Rigid-Body State-Change Simulations in LAMMPS (HOOMD-like Behavior)

This doc is a practical checklist for implementing **rigid-body** state-change simulations in LAMMPS that behave like a HOOMD-blue rigid body (e.g., `md.constrain.Rigid` / “glue” style): stable rotations, predictable thermostating, and safe mid-run “state changes”.

It is written for future humans/agents extending the `sim_templates/state_change/` workflows (octahedron / dimer / ksat).

---

## Core principles (read first)

- **Rigid bodies are not “automatic” in LAMMPS the way they often feel in HOOMD.**
  - LAMMPS rigid fixes maintain internal state and require correct integrator/thermostat choices.
  - If you change atom types or per-atom properties incorrectly, you can break rigid-body assumptions and get artifacts (spinning, “lost atoms”, or wrong dynamics).

- **Mass affects rotation; repulsion affects overlap.**  
  In LAMMPS, you can (and often should) give patch sites non-negligible mass to stabilize rotational inertia **without adding patch repulsion**.

- **State-change “logic correctness” and “MPI correctness” are separate problems.**
  - In serial you can get correct logic but still fail under MPI due to ghost atoms and stale state.

---

## Rigid body setup: avoid “crazy spin” and match HOOMD expectations

### Use the “two-fix” pattern (recommended)

Prefer:
- `fix rigid/nve molecule` **(integrator, no thermostat)**  
- `fix langevin ...` **on the same rigid group** **(thermostat/friction)**

This tends to be more stable and less “self-spinning” than `rigid/nvt` in many setups, and mirrors the “integrate rigid + thermostat” mental model from HOOMD.

Example pattern (conceptual):
- `group all_monomers molecule > 0`
- `fix fx_nve all_monomers rigid/nve molecule`
- `fix fx_langevin all_monomers langevin T T damp seed`

### Damping (`Tdamp`) matters a lot

- Smaller `damp` in `fix langevin` = **stronger friction**.
- If monomers spin too much, reduce `damp` (e.g. `1.0 → 0.5 → 0.2`).

### Mass distribution matters (moment of inertia)

To avoid near-singular inertia:
- Give patches real mass (not ~0).
- Keep total mass per rigid body consistent across species if possible.

Rule of thumb:
- **Cores**: majority mass
- **Patches**: enough mass so \(I = \sum m_i r_i^2\) is not tiny

### Time step

When using stiff short-range attractions (e.g., Morse with large `alpha`), use smaller dt.
If you see NaNs / instability:
- Reduce `timestep` (e.g. `0.005 → 0.002 → 0.001`)

---

## Patch attraction / overlap: “why aren’t dimers forming?”

### Ensure patches are allowed to overlap (no patch–patch repulsion unless intended)

If you want patch “overlap” assembly:
- Patch–patch repulsion should be **off** (or very weak)
- Core–core repulsion should be **on** (excluded volume)

### Verify you’re using Morse correctly

LAMMPS `pair_style morse` takes parameters:
- `D0 alpha r0` (and the style’s global cutoff)

Important:
- `r0` sets the **minimum location**.  
If you conceptually want “contact/overlap is best”, you likely want `r0 = 0.0`.

If `r0` is accidentally set to the cutoff distance, you can effectively remove attraction at short range.

### Relative strengths

If you want “ratcheting”:
- Ensure **sink** interaction is much stronger than **transition** interaction.
  - e.g. \(D0_{33} = 2 \times D0_{23}\) or more.

---

## State change in LAMMPS: safe patterns

There are two broad patterns:

### Pattern A — “Type flip only”

If species differ only by type labels (same geometry):
- The safest is to flip only patch atom types (e.g., `2 → 3`) for the *whole molecule*.

Keep in mind:
- Your dynamics must allow the RB dimer to form before flip if you want to visually see it.
  - Use hysteresis (see below).

### Pattern B — “Geometry changes” (avoid if possible)

If species differ by geometry:
- A type flip is not enough.
- You must either:
  - Represent both geometries as sites and activate/deactivate them (dummy types), **or**
  - Rebuild the rigid body definition mid-run (hard and error-prone).

---

## Kinetics: prevent flips on fleeting contacts (hysteresis)

If you check every `nevery` steps and `P_flip=1.0`, you will flip on the first momentary 2–3 touch.

To make flips correspond to real dimers:
- Require sustained contact for N consecutive checks:
  - `hysteresis_checks = 5` means `5 * nevery` steps of continuous contact.

Tuning:
- Smaller hysteresis → more flips, less “dimer-looking”
- Larger hysteresis → fewer flips, more stable RB dimers before flip

---

## MPI safety (especially for complex fixes like octahedron)

If your fix stores per-atom state (timestamps, effective types, counters), you must:

- **Forward-communicate** state to ghosts (`comm_forward` and pack/unpack) if ghost atoms participate in checks.
- Be careful with “consistency sweeps”:
  - They should fix split-brain states across MPI boundaries **without creating fake state changes**.
  - Never set “last_change = timestep” for a sweep correction unless you intend to trigger kinetics.

---

## Debugging checklist

When behavior looks wrong:

- **No assembly / no dimers**
  - Check `r0` in Morse
  - Check patch–patch repulsion is not accidentally on
  - Increase D0 and/or lower temperature
  - Reduce initial spacing (without changing concentration) to increase encounter rate

- **Crazy spinning**
  - Increase patch mass / improve mass distribution
  - Switch to `rigid/nve + langevin`
  - Reduce Langevin `damp`
  - Reduce timestep

- **State changes happen “without reason”**
  - Add hysteresis
  - Ensure neighbor/contact checks include ghost atoms if running MPI
  - Ensure any “consistency sweep” does not invent new timestamps

---

## Where things live in this repo

- Dimer generator: `sim_templates/state_change/generate_dimer_cpp.py`
- Dimer fix: `sim_templates/state_change/dimer/fix_state_change_dimer.*`
- Dimer slurm: `sim_templates/state_change/submit_dimer.slurm`
- Octahedron fix and docs: `sim_templates/state_change/octahedron/`


