# Competition Organization Structure

This document describes the directory organization for the molecular computing competition framework.

## Directory Overview

```
lammps-state-change/
├── problems/              # Competition problem definitions
│   └── problem-001-dimer-ksat/
│       ├── problem.json              # Problem definition
│       ├── baseline_policy.json      # Reference solution
│       ├── baseline_params.json     # Default parameters
│       ├── analyze_submission.sh     # Problem-specific analysis script
│       ├── leaderboard.csv           # Final aggregated leaderboard
│       └── starter_template/         # Participant templates
│
├── submissions/          # Competitor submissions (gitignored)
│   └── problem-001-dimer-ksat/
│       └── competitor_username/
│           ├── submission.json       # Original (from competitor)
│           ├── policy.json           # Original (from competitor)
│           ├── params.json           # Original (optional)
│           ├── generated/            # Auto-generated C++ (gitignored)
│           ├── simulation/            # Auto-generated LAMMPS (gitignored)
│           └── results/               # Simulation outputs (gitignored)
│               ├── dump.*.lammpstrj
│               ├── lammps_stdout.log
│               └── analysis/
│                   ├── timeseries.csv
│                   └── leaderboard_row.csv
│
├── tests/                # Test outputs (gitignored)
│   └── problem-001-dimer-ksat/
│       └── baseline_test/
│           ├── generated/
│           └── simulation/
│
├── core/                 # Framework infrastructure
│   ├── generators/      # Code generators
│   ├── schemas/         # JSON schemas
│   ├── benchmark/        # Scoring system
│   └── analysis/        # Analysis tools (used by problem scripts)
│
└── tools/               # Regulator tools
    ├── import_submission.sh          # Import submission files
    ├── process_submission.sh         # Process single submission
    ├── process_all_submissions.sh    # Process all for a problem
    ├── analyze_submission_results.sh # Analyze simulation results
    ├── update_leaderboard.sh         # Aggregate into leaderboard
    └── test_pipeline.sh              # Test pipeline
```

## Parallel Structure: Problems ↔ Submissions

The `submissions/` directory mirrors `problems/` structure:

| Problems Directory | Submissions Directory |
|-------------------|----------------------|
| `problems/problem-001-dimer-ksat/` | `submissions/problem-001-dimer-ksat/competitor_username/` |
| `problems/problem-002-octahedron/` | `submissions/problem-002-octahedron/competitor_username/` |
| `problems/problem-003-hamiltonian-path/` | `submissions/problem-003-hamiltonian-path/competitor_username/` |

## Problem-Specific Analysis

**Key Design Decision**: Each problem has its own analysis script because:
- Different problems have different species/types
- Different yield calculations (species fraction vs cluster-based)
- Different metrics and scoring
- Different trajectory analysis requirements

**Location**: `problems/{problem_id}/analyze_submission.sh`

**Called by**: `tools/analyze_submission_results.sh`

**Output**: `submissions/.../results/analysis/leaderboard_row.csv`

## Complete Workflow

### 1. Import Submission
```bash
tools/import_submission.sh submission_001_username.zip
# → submissions/problem-001-dimer-ksat/username/
```

### 2. Process Submission
```bash
tools/process_submission.sh submissions/problem-001-dimer-ksat/username/
# → Generates C++ fix and LAMMPS files
```

### 3. Run Simulation
```bash
cd submissions/problem-001-dimer-ksat/username/simulation/
lmp_mpi -in in.*
# → Outputs to results/
```

### 4. Analyze Results
```bash
tools/analyze_submission_results.sh submissions/problem-001-dimer-ksat/username/
# → Calls problems/problem-001-dimer-ksat/analyze_submission.sh
# → Generates results/analysis/leaderboard_row.csv
```

### 5. Update Leaderboard
```bash
tools/update_leaderboard.sh problem-001-dimer-ksat
# → Aggregates all leaderboard_row.csv files
# → Outputs to problems/problem-001-dimer-ksat/leaderboard.csv
```

## Git Organization

### Tracked (Committed)
- `problems/` - Problem definitions and analysis scripts
- `problems/*/leaderboard.csv` - Final leaderboards
- `core/` - Framework code
- `tools/` - Processing scripts
- `submissions/README.md` - Documentation
- `submissions/.gitkeep` - Directory structure

### Ignored (Gitignored)
- `submissions/*/` - All competitor submissions and generated files
- `tests/` - Test outputs
- `**/generated/` - Auto-generated C++ files
- `**/simulation/` - LAMMPS input files (generated)
- `**/results/` - Simulation outputs

## Leaderboard Data Flow

```
submissions/problem-001-dimer-ksat/username/results/analysis/
└── leaderboard_row.csv  (individual result)
         ↓
tools/update_leaderboard.sh aggregates
         ↓
problems/problem-001-dimer-ksat/
└── leaderboard.csv  (aggregated, sorted, tracked)
```

## Key Files Per Submission

After complete processing:

```
submissions/problem-001-dimer-ksat/username/
├── submission.json          # Original (from competitor)
├── policy.json              # Original (from competitor)
├── params.json              # Original (optional)
│
├── generated/               # Auto-generated by process_submission.sh
│   ├── fix_state_change_*.h
│   └── fix_state_change_*.cpp
│
├── simulation/              # Auto-generated by process_submission.sh
│   ├── data.*
│   └── in.*
│
└── results/                 # Created by simulation + analysis
    ├── dump.*.lammpstrj     # Trajectory (from simulation)
    ├── lammps_stdout.log    # Thermodynamic output
    ├── *.err                # Event logs
    └── analysis/            # Created by analyze_submission_results.sh
        ├── timeseries.csv   # Full timeseries data
        └── leaderboard_row.csv  # Single row for leaderboard
```

## Benefits of This Structure

1. **Clear Separation**: Problems (definitions) vs Submissions (competitor work)
2. **Parallel Organization**: Easy to find submissions for each problem
3. **Problem-Specific Analysis**: Each problem can have custom analysis
4. **Automated Processing**: Tools can iterate over submissions systematically
5. **Scalable**: Works for any number of problems and competitors
6. **Clean Repository**: Only essential files are tracked
7. **Local Development**: All submission processing stays local (not pushed to remote)

## Adding New Problems

When adding a new problem:

1. Create `problems/problem-NNN-name/`
2. Create `problem.json` with problem definition
3. Create `analyze_submission.sh` with problem-specific analysis
4. Create `submissions/problem-NNN-name/` (empty, with .gitkeep)
5. Update tools to handle new problem ID format

## Example: Complete Workflow

```bash
# 1. Import submission
tools/import_submission.sh submission_001_alice.zip

# 2. Process submission
tools/process_submission.sh submissions/problem-001-dimer-ksat/alice/

# 3. Run simulation (manual)
cd submissions/problem-001-dimer-ksat/alice/simulation/
lmp_mpi -in in.*
mv dump.*.lammpstrj ../results/
mv lammps_stdout.log ../results/

# 4. Analyze results
tools/analyze_submission_results.sh submissions/problem-001-dimer-ksat/alice/

# 5. Update leaderboard (after processing all submissions)
tools/update_leaderboard.sh problem-001-dimer-ksat

# 6. View leaderboard
cat problems/problem-001-dimer-ksat/leaderboard.csv
```
