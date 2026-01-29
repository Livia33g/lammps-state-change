# State-change policy benchmark (Kaggle-style)

This benchmark is meant to compare **state-change policies** for out-of-equilibrium LAMMPS simulations where flips represent **energy consumption** (directed energy input).

It follows the framing in [`file://Simons_Molecular_Computing.pdf`](file://Simons_Molecular_Computing.pdf):

- **Encode** a problem into interactions/components
- **Drive** a non-equilibrium protocol (your state-change policy + energy input)
- **Decode/read out** the produced “soup” into an answer (a target structure / yield)

The goal is to rank policies by **how much target progress they achieve per unit energy consumed**, while also penalizing slow protocols.

## The template workflow (minimal changes per new simulation)

1. **Run your simulation** (any statechange system).
2. Create a **task JSON** describing:
   - where the dump/thermo/events files are
   - what the target is (species yield or cluster yield)
3. Run the task runner:
   - `run_task.py` produces a standardized timeseries CSV/PNG + a 1-row leaderboard CSV
4. Repeat for many runs/policies and aggregate:
   - `aggregate_leaderboard.py` combines `*.leaderboard.csv` files into one table.

Task spec documentation:
- `task_schema.md`

Example tasks:
- `tasks/`

## Inputs (what a simulation must produce)

For each run directory, you should generate a **timeseries CSV** using:

- `sim_templates/state_change/analyze_trajectory_target_yield_and_work.py`

This CSV contains:
- `timestep`
- `yield` (target yield at that timestep)
- `dpe` (ΔPotEng between consecutive thermo outputs; a “kick/work proxy”)
- `n_statechanges_interval` (# flips in that thermo interval)

## Scoring outputs

Use:

- `score_policy_from_timeseries.py`

It produces a single “leaderboard row” per run:
- `final_yield`
- `t_reach_yield_<threshold>` (first timestep reaching a chosen yield threshold)
- `n_flips_total`
- `work_abs_total` = \(\sum |\Delta PE|\)
- `work_abs_flip_intervals` = \(\sum |\Delta PE|\) over intervals where flips occurred
- `work_per_yield` = `work_abs_flip_intervals / final_yield`
- `score` (one scalar; see script for exact definition)

## Why two “work” columns?

Your current LAMMPS outputs do not report the instantaneous per-flip energy change.
So:
- `work_abs_total` captures “energy activity” overall (includes thermal/relaxation)
- `work_abs_flip_intervals` focuses attention on intervals where flips happened, which is a better proxy for energy-driven cost.

If/when we modify the fix to print per-event `dE`, we can replace the proxy with a more faithful energy cost.

## Example

```bash
python3 sim_templates/state_change/benchmark/run_task.py \
  --task sim_templates/state_change/dimer_ksat/variants/1core/analysis/tasks/dimer_ksat_1core_C_yield_15466594.json
```

Aggregate many runs into one leaderboard:

```bash
python3 sim_templates/state_change/benchmark/aggregate_leaderboard.py \
  --dir sim_templates/state_change/analysis \
  --out sim_templates/state_change/analysis/leaderboard.csv \
  --sort-by score --desc
```

## Next iteration (recommended)

For serious comparisons we should extend the benchmark to:
- evaluate multiple seeds / replicates and average scores
- validate “correctness” (e.g. avoid undesired byproducts) via structure checks
- record energy cost as an explicit model parameter (e.g. Δμ per flip) and score in that unit


