# Trajectory-based yield + “work” (ΔPE) analysis

This folder has state-change simulations that are out of equilibrium. A convenient proxy for the “kick/work” you described is:

- **Work(t)** = \( \Delta PE \) between two consecutive **thermo** outputs (e.g. every 1000 steps).

Separately, if you have a **target assembly** in mind, you usually want to track its **yield vs time** from the trajectory.

This repo now includes a single script that does both:

- `analyze_trajectory_target_yield_and_work.py`

## What the script needs

- **Dump trajectory**: a LAMMPS `.lammpstrj` produced by `dump ... custom ... id mol type x y z`
- **Thermo output**: `log.lammps` or `lammps_stdout.log` that contains a thermo table including `Step` and `PotEng`
- **Optional**: Slurm `.err` (or any file) containing `STATECHANGE ... step <t>` lines. If provided, the script counts how many statechanges happened in each thermo interval.

## How “yield” is defined (cluster-based)

The script detects clusters at a timestep by building a graph on **molecules**:

- **Nodes**: molecule ids (`mol` column in the dump)
- **Edges**: between molecule A and B if any pair of “site” atoms (selected by `--site-types`) are within `--bond-cutoff` (with PBC).
- **Clusters**: connected components of this graph

Then you define a “target structure” as:

- **Target size**: `--target-size K` (cluster of K molecules)
- **Optional** exact composition**: `--target-composition "label:count,..."` where each molecule’s label is derived from its site atom types (see `--label-mode`).

Default yield is the **fraction of molecules** that are in target clusters.

## Example: C-yield for dimer_ksat (A(2) → C(4))

If your goal is simply “fraction of molecules that are now C”, use `species_fraction`:

```bash
python3 analyze_trajectory_target_yield_and_work.py \
  --dump dimer_ksat_1core_simulation_cpp/dump.dimer_ksat_1core.lammpstrj \
  --thermo dimer_ksat_1core_simulation_cpp/lammps_stdout.log \
  --events slurm_dimer_ksat_1core-15466594.err \
  --site-types 2 3 4 \
  --bond-cutoff 0.7 \
  --yield-mode species_fraction \
  --species-label 4 \
  --out analysis/dimer_ksat_1core_C_yield
```

## Example: dimer target (cluster size 2)

From `state_change/`:

```bash
python3 analyze_trajectory_target_yield_and_work.py \
  --dump dimer_ksat_1core_simulation_cpp/dump.dimer_ksat_1core.lammpstrj \
  --thermo dimer_ksat_1core_simulation_cpp/lammps_stdout.log \
  --events slurm_dimer_ksat_1core-15466594.err \
  --site-types 2 3 4 \
  --bond-cutoff 0.7 \
  --target-size 2 \
  --out analysis/dimer_ksat_1core_target_dimer
```

This writes:
- `analysis/dimer_ksat_1core_target_dimer.csv`
- `analysis/dimer_ksat_1core_target_dimer.png`

## Example: require exact A:B composition (1 A + 1 B dimer)

If your labels are literally the patch type integers (common in your dimer/ksat generators), you can require composition.

Example “one type-2 + one type-3 molecule” in the cluster:

```bash
python3 analyze_trajectory_target_yield_and_work.py \
  --dump dimer_ksat_1core_simulation_cpp/dump.dimer_ksat_1core.lammpstrj \
  --thermo dimer_ksat_1core_simulation_cpp/lammps_stdout.log \
  --site-types 2 3 4 \
  --bond-cutoff 0.7 \
  --target-size 2 \
  --target-composition "2:1,3:1" \
  --out analysis/dimer_ksat_1core_AB_dimer
```

## Output columns

The CSV contains, at thermo timesteps:

- `timestep`
- `pe`: PotEng at this step
- `dpe`: \(PE(t_i) - PE(t_{i-1})\) (your “work/kick” proxy)
- `cum_work`: cumulative sum of `dpe`
- `yield`: target yield (by `--yield-mode`)
- `n_target_clusters`
- `n_molecules`
- `n_statechanges_interval`: number of `STATECHANGE` events in \((t_{i-1}, t_i]\) (0 if no `--events`)

## Notes / limitations

- A LAMMPS dump **does not contain energies** by default, so this script uses the thermo stream for PE.
- The current “target structure” detection is **cluster-based** (size + optional composition). If your target is *topologically* specific (e.g. a particular graph), we can extend this to motif matching once you define the target constraints.


