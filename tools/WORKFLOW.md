# Competition Workflow for Regulators

Complete workflow from receiving submissions to updating leaderboards.

## Overview

```
Import → Process → Run Simulation → Analyze → Update Leaderboard
```

## Step-by-Step Workflow

### 1. Import Submission

Receive submission file (zip/tar/directory) with problem number in filename.

```bash
# Import automatically detects problem number from filename
tools/import_submission.sh submission_001_alice.zip
tools/import_submission.sh submission_001_bob.tar.gz
tools/import_submission.sh submission_directory/
```

**What it does:**
- Extracts problem number (001, 002, etc.) from filename
- Extracts username from filename or submission.json
- Organizes into: `submissions/problem-001-dimer-ksat/username/`

**Output structure:**
```
submissions/problem-001-dimer-ksat/username/
├── submission.json
├── policy.json
└── params.json (optional)
```

### 2. Process Submission

Generate C++ fix and LAMMPS simulation files.

```bash
tools/process_submission.sh submissions/problem-001-dimer-ksat/username/
```

**What it does:**
- Validates submission (JSON format, schema compliance)
- Generates C++ fix → `submissions/.../generated/`
- Generates LAMMPS files → `submissions/.../simulation/`
- Creates `results/` directory

**Output structure:**
```
submissions/problem-001-dimer-ksat/username/
├── submission.json
├── policy.json
├── params.json
├── generated/               # C++ fix files
│   ├── fix_state_change_*.h
│   └── fix_state_change_*.cpp
└── simulation/              # LAMMPS input files
    ├── data.*
    └── in.*
```

### 3. Run Simulation

Run LAMMPS simulation (manual or via SLURM).

```bash
cd submissions/problem-001-dimer-ksat/username/simulation/
lmp_mpi -in in.*
# Or submit via SLURM
```

**Output:**
- Trajectory file: `dump.*.lammpstrj`
- Thermodynamic output: `lammps_stdout.log`
- Event logs: `*.err` (contains state change events)

**Move outputs to results/:**
```bash
mv dump.*.lammpstrj ../results/
mv lammps_stdout.log ../results/
mv *.err ../results/
```

### 4. Analyze Results

Run problem-specific analysis script.

```bash
tools/analyze_submission_results.sh submissions/problem-001-dimer-ksat/username/
```

**What it does:**
- Calls `problems/problem-001-dimer-ksat/analyze_submission.sh`
- Problem-specific script:
  - Generates timeseries CSV from trajectory
  - Computes metrics (yield, work, score)
  - Creates leaderboard_row.csv with metadata

**Output structure:**
```
submissions/problem-001-dimer-ksat/username/
└── results/
    ├── dump.*.lammpstrj
    ├── lammps_stdout.log
    ├── *.err
    └── analysis/
        ├── timeseries.csv        # Full timeseries data
        ├── leaderboard_row.csv   # Single row for leaderboard
        └── analysis.log          # Analysis output log
```

### 5. Update Leaderboard

Aggregate all submission results into problem leaderboard.

```bash
tools/update_leaderboard.sh problem-001-dimer-ksat
```

**What it does:**
- Finds all `leaderboard_row.csv` files in submissions
- Aggregates into single CSV
- Sorts by primary metric (work_per_yield)
- Outputs to: `problems/problem-001-dimer-ksat/leaderboard.csv`

**Leaderboard structure:**
```
problems/problem-001-dimer-ksat/
└── leaderboard.csv  # Aggregated results, sorted by score
```

## Batch Processing

### Process All Submissions for a Problem

```bash
# Process all submissions (generate C++ and LAMMPS files)
tools/process_all_submissions.sh problem-001-dimer-ksat

# After running simulations, analyze all
for dir in submissions/problem-001-dimer-ksat/*/; do
    tools/analyze_submission_results.sh "$dir"
done

# Update leaderboard
tools/update_leaderboard.sh problem-001-dimer-ksat
```

## Problem-Specific Analysis

Each problem has its own analysis script:

```
problems/
├── problem-001-dimer-ksat/
│   ├── problem.json
│   ├── analyze_submission.sh    # Problem-specific analysis
│   └── leaderboard.csv          # Final leaderboard
│
├── problem-002-octahedron/
│   ├── problem.json
│   ├── analyze_submission.sh    # Different analysis for octahedron
│   └── leaderboard.csv
│
└── problem-003-hamiltonian-path/
    ├── problem.json
    ├── analyze_submission.sh    # Different analysis for graph problem
    └── leaderboard.csv
```

Each `analyze_submission.sh`:
- Reads problem-specific parameters from `problem.json`
- Uses appropriate analysis tools for that problem type
- Generates metrics relevant to that problem
- Outputs standardized `leaderboard_row.csv` format

## Leaderboard Format

The final `leaderboard.csv` contains:

| Column | Description |
|--------|-------------|
| `problem_id` | Problem identifier |
| `username` | Competitor username |
| `team_name` | Team/display name |
| `policy_name` | Policy name from policy.json |
| `submission_date` | When submitted |
| `policy_version` | Policy version |
| `final_yield` | Final target yield achieved |
| `work_per_yield` | Primary metric (lower is better) |
| `n_flips_total` | Total state changes |
| `time_to_threshold` | Steps to reach threshold |
| `score` | Overall score |

## Directory Organization Summary

```
submissions/                    # All competitor submissions (gitignored)
└── problem-001-dimer-ksat/
    └── username/
        ├── submission.json      # Original submission
        ├── policy.json
        ├── params.json
        ├── generated/           # C++ files (gitignored)
        ├── simulation/          # LAMMPS files (gitignored)
        └── results/            # Simulation outputs (gitignored)
            └── analysis/
                ├── timeseries.csv
                └── leaderboard_row.csv

problems/                       # Problem definitions (tracked)
└── problem-001-dimer-ksat/
    ├── problem.json            # Problem definition
    ├── analyze_submission.sh   # Problem-specific analysis
    └── leaderboard.csv         # Final leaderboard (tracked)
```

## Notes

- **Submissions are gitignored**: Only problem definitions and leaderboards are tracked
- **Analysis is problem-specific**: Each problem has its own analysis script
- **Results stay local**: All simulation outputs stay on your machine
- **Leaderboards are tracked**: Final aggregated leaderboards can be committed

