# Benchmark task spec (JSON)

Each **benchmark task** describes one simulation run to analyze + score.

The goal is: **minimal changes across new statechange simulations** — you only write a new
`task.json`, then run the same benchmark runner.

## Required fields

- `name` *(string)*: unique identifier for this run (used for output filenames)
- `dump` *(string path)*: LAMMPS dump file (`.lammpstrj`) with `id mol type x y z`
- `thermo` *(string path)*: LAMMPS thermo output (`log.lammps` or `lammps_stdout.log`)
- `site_types` *(list[int])*: site atom types used to label molecules and/or form bonds
- `bond_cutoff` *(float)*: distance cutoff used for site–site bonds (required even for species yield; used for clustering modes)

## Optional fields

- `events` *(string path)*: file containing `STATECHANGE ... step <t>` lines (e.g. slurm `.err`)

### Yield definition

Choose one of:

#### (A) Species yield (common for “A→C” style runs)

- `yield_mode`: `"species_fraction"`
- `species_label`: *(int)* label to count as “target” (e.g. 4 for C in A(2)→C(4))

#### (B) Cluster yield

- `yield_mode`: `"fraction_molecules"` or `"n_clusters"`
- `target_size`: *(int)* cluster size in molecules
- `target_composition`: *(string, optional)* exact label composition like `"2:1,3:1"`

### Labeling

- `label_mode` *(string, default `"majority_site_type"`)*:
  - `"majority_site_type"`, `"min_site_type"`, `"max_site_type"`

### Sampling + performance

- `sample` *(string, default `"thermo"`)*: `"thermo"` or `"all"`
- `max_frames` *(int, optional)*: cap frames processed (debugging)

### Scoring

- `yield_threshold` *(float, default 0.6)*: threshold for “time-to-yield” metric

### Output control

- `out_dir` *(string path, default: `sim_templates/state_change/analysis/`)*: where to write outputs

## Example task

```json
{
  "name": "dimer_ksat_1core_C_yield_15466594",
  "dump": "../dimer_ksat_1core_simulation_cpp/dump.dimer_ksat_1core.lammpstrj",
  "thermo": "../dimer_ksat_1core_simulation_cpp/lammps_stdout.log",
  "events": "../slurm_dimer_ksat_1core-15466594.err",
  "site_types": [2, 3, 4],
  "bond_cutoff": 0.7,
  "yield_mode": "species_fraction",
  "species_label": 4,
  "yield_threshold": 0.6
}
```


