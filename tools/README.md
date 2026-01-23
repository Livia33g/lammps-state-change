# Regulator Tools

These tools are for **competition regulators only** (not participants). They help process and manage submissions locally.

## Tools Overview

### `import_submission.sh`
Import competitor submission files (zip/tar/directory) and organize by problem number.

**Usage:**
```bash
tools/import_submission.sh submission_001_username.zip
```

**What it does:**
- Extracts problem number from filename (e.g., "001")
- Extracts username from filename or submission.json
- Organizes into: `submissions/problem-001-dimer-ksat/username/`

### `test_pipeline.sh`
Test the pipeline with baseline policies before processing real submissions.

**Usage:**
```bash
tools/test_pipeline.sh [problem_id] [test_name]
tools/test_pipeline.sh problem-001-dimer-ksat baseline_test
```

**Output:** `tests/{problem_id}/{test_name}/`

### `process_submission.sh`
Process a single competitor submission: validate, generate C++ fix, generate LAMMPS files.

**Usage:**
```bash
tools/process_submission.sh submissions/problem-001-dimer-ksat/competitor_username/
```

**What it does:**
1. Validates submission (JSON format, schema compliance)
2. Generates C++ fix files → `submissions/.../generated/`
3. Generates LAMMPS input files → `submissions/.../simulation/`
4. Creates `results/` directory for outputs

### `process_all_submissions.sh`
Process all submissions for a given problem.

**Usage:**
```bash
tools/process_all_submissions.sh problem-001-dimer-ksat
```

**What it does:**
- Finds all competitor directories in `submissions/{problem_id}/`
- Processes each one using `process_submission.sh`
- Reports success/failure counts

### `analyze_submission_results.sh`
Analyze simulation results using problem-specific analysis script.

**Usage:**
```bash
tools/analyze_submission_results.sh submissions/problem-001-dimer-ksat/username/
```

**What it does:**
- Calls `problems/{problem_id}/analyze_submission.sh`
- Problem-specific script generates:
  - `results/analysis/timeseries.csv`
  - `results/analysis/leaderboard_row.csv`

### `update_leaderboard.sh`
Aggregate all submission results into problem leaderboard.

**Usage:**
```bash
tools/update_leaderboard.sh problem-001-dimer-ksat
```

**What it does:**
- Collects all `leaderboard_row.csv` files from submissions
- Aggregates into `problems/{problem_id}/leaderboard.csv`
- Sorts by primary metric

## Important Notes

- **Local Development Only**: These tools and their outputs are for local use
- **Not Pushed to Remote**: `submissions/` and `tests/` directories are gitignored
- **Results Later**: After processing and running simulations, results will be aggregated into leaderboards (separate process)

## Workflow

1. **Test pipeline** (before competition):
   ```bash
   tools/test_pipeline.sh problem-001-dimer-ksat baseline_test
   ```

2. **Receive submissions** (competitors create directories):
   - `submissions/problem-001-dimer-ksat/competitor_alice/`
   - `submissions/problem-001-dimer-ksat/competitor_bob/`

3. **Process submissions**:
   ```bash
   tools/process_all_submissions.sh problem-001-dimer-ksat
   ```

4. **Review generated files**:
   - Check `submissions/.../generated/` for C++ fixes
   - Check `submissions/.../simulation/` for LAMMPS inputs

5. **Run simulations** (separate process, not automated here)

6. **Aggregate results** (later, separate process):
   - Extract scores from simulation outputs
   - Update `problems/{problem_id}/leaderboard.csv`

## Directory Structure

After processing, each submission contains:

```
submissions/problem-001-dimer-ksat/competitor_username/
├── submission.json          # Original (from competitor)
├── policy.json              # Original (from competitor)
├── params.json              # Original (optional)
├── generated/               # Auto-generated C++ (gitignored)
├── simulation/              # Auto-generated LAMMPS (gitignored)
└── results/                 # Simulation outputs (gitignored)
```

All generated files are gitignored - only the original submission files would be tracked (if you choose to track them separately).

