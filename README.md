# `state_change/` (template + benchmarked state-change simulations)

This repo contains **LAMMPS state-change fixes + simulation templates** where flips represent **fuel consumption** (out-of-equilibrium protocols).

## Directory convention (important for git hygiene)

### What stays in the repo root (shared across all systems)

- **Shared analysis tools**: `analyze_*.py`, `ANALYSIS_WORK_YIELD_FROM_TRAJECTORY.md`
- **Benchmark harness**: `benchmark/` (task runner + scoring + aggregation)
- **Build/rebuild helpers**: `rebuild_manual.sh`, `REBUILD_*.md`, `rebuild_all_fixes.slurm`, `add_new_fix.sh`
- **High-level docs**: `INSTRUCTIONS.md`, `STATE_CHANGE_EXPLANATION.md`, etc.

### What belongs in a system subdirectory

Each system folder should contain only **templates**, not outputs:

- **C++ fix source** (if system/variant has custom logic)
- **Generator(s)** (produce `simulation_cpp/` or another ignored output folder)
- **Template submit script(s)** (`submit.slurm`) to run the simulation
- **Local analysis wrappers**:
  - `analysis/tasks/*.json` (benchmark task configs)
  - `analysis/submit_analysis.slurm` (runs the shared benchmark runner)

Examples:
- `dimer/`
- `dimer_ksat/variants/<variant>/`
- `octahedron/`
- `ksat/`

### What must NOT be committed

The `.gitignore` is set up to prevent committing:
- **Trajectories/dumps**: `*.lammpstrj`
- **Runtime logs**: `*.out`, `*.err`, `*.log`, `log.lammps`
- **Generated run directories**: `**/simulation_cpp/` and other legacy generated folders
- **Analysis outputs**: `**/analysis/*.csv`, `**/analysis/*.png`, `**/analysis/*.leaderboard.csv`
- **Local scratch archives**: `back/` and `legacy/`

If you accidentally tracked an output in the past, you must untrack it once:

```bash
git rm -r --cached <path>
```

### About `back/` and `legacy/`

These directories are **local-only** archival space. They are ignored and are not meant to be pushed.



