## Summary

- Purpose: run and analyze small three-particle patchy particle simulations with Morse potential patch interactions. Several variants of LAMMPS input files and recorded output trajectories are included, plus lightweight analysis CSVs and a Python counting script.

## Files

- `part1tripatch_mid.mol`
  - A `.mol` formatted file containing particle/patch specification. 

- `three_particles_morse_large_fixedT_s_shielding_test.lammps`
  - Primary LAMMPS input script. Contains system setup, potentials (Morse/patch), integrator settings, and dump/write commands. Use this as the canonical input to reproduce the large-scheme three-particle run.

- `count_abc_python.py`
  - Python script used to process LAMMPS trajectory/dump files and produce the `abc_counts_*.csv` outputs. Counts the number of A-B-C structures in three ways:
  nABC_pure (A-B-C, no branching at all), nABC_clean (B has *only* one A neighbour and *only* one C neighbour, but these can bond to other B's), nABC_all (B has *at least* one A neighbour and *at least* one B neighbour, branching allowed.)

- `three_particles_morse_large.lammps`, `three_particles_morse_large_fixedT_s.lammps`, `three_particles_morse_large_fixedT_s_shielding.lammps`
  - Variants of the main LAMMPS input that apply a fixed-temperature thermostat and different attempts at shielding/steric guards. 

- `abc_counts_26.csv`, `abc_counts_30.csv`, `abc_counts_30_long.csv`, `abc_counts_34.csv`
  - CSV result files produced with morse potential "steric guarding" $\sigma \in \{0.26, 0.30, 0.34\}$. "long" = longer simulation time.
    
## Dump directory

- `dump/`
  - Contains LAMMPS trajectory files in `.lammpstrj` format and a `test/` subdirectory with additional small test trajectories.

  - `dump/dump_tripatch_3distinct_morsechanges_new_large_longtime.lammpstrj`
    - A long-time trajectory from a run that explores Morse parameter changes. Use this for extended-time analysis or visualization.

  - `dump/test/dump_tripatch_3distinct_shielding_26.lammpstrj`
  - `dump/test/dump_tripatch_3distinct_small_test_shielding_26.lammpstrj`
  - `dump/test/dump_tripatch_3distinct_small_test_shielding_30_long.lammpstrj`
  - `dump/test/dump_tripatch_3distinct_small_test_shielding_30.lammpstrj`
  - `dump/test/dump_tripatch_3distinct_small_test_shielding_34.lammpstrj`
    - Small test trajectories with different parameter sets (numbers 26, 30, 34 in the names). 

## Quick tips
- Inspect LAMMPS input files
  - Open the `three_particles_*.lammps` files to confirm which pair/style potentials and patch interactions are used, and which dump commands are enabled. Look for `pair_style`, `pair_coeff`, `fix`, and `dump` statements.

- View trajectories
  - Visualize `.lammpstrj` files with OVITO (recommended) or convert them for other viewers. OVITO opens `.lammpstrj` directly and can animate, color by type, and compute cluster/motif analyses.

- Reproduce a run (example)

```bash
# Single-process LAMMPS (executable name may vary on your system)
lammps -in three_particles_morse_large.lammps

# Or with MPI (adjust -np and binary name as appropriate):
mpirun -np 4 lammps -in three_particles_morse_large.lammps
```

- Run the analysis script
  - Inspect `count_abc_python.py` to see expected input arguments. A typical pattern is to feed it a dump file and capture the CSV output. For example:

```bash
python3 count_abc_python.py dump/dump_tripatch_3distinct_small_test_shielding_30.lammpstrj > abc_counts_30.csv
```
