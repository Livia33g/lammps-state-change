# legacy/

This folder contains **older / one-off / intermediate** workflows that we keep for reference but don’t want cluttering the main `state_change/` directory.

Nothing here is required for the **current active workflows**:
- Dimer: `fix_state_change_dimer.*`, `generate_dimer_cpp.py`, `submit_dimer.slurm`
- Octahedron: `octahedron/`
- Ksat: `ksat/`
- Inverse design: `inverse-design/`

## What’s inside

- `docs/`
  - Older notes like `README_CPP.md`, `NEXT_STEPS.md`
- `old_state_change_fix/`
  - The older `fix_state_change/` pipeline (pre dimer-specific fix)
- `old_scripts/`
  - One-off scripts (e.g., old GitHub push helper, old analysis scripts)
- `old_slurm/`
  - Older SLURM submit scripts used during debugging/optimization
- `test_generation/`
  - Small historical test inputs


