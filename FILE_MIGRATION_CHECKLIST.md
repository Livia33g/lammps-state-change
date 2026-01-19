# File Migration Checklist

This document lists every file that needs to be moved, archived, or deleted during the reorganization.

## Phase 1: Create New Directory Structure

```bash
# Create new directories
mkdir -p framework/{generators,analysis,benchmark,build,validation}
mkdir -p framework/generators/templates
mkdir -p framework/benchmark/schemas
mkdir -p problems/{octahedron,dimer,dimer_ksat/variants/{1core,base}}
mkdir -p submissions/{octahedron,dimer}/.gitkeep
mkdir -p docs/{user_guide,developer_guide,reference,reference/examples}
mkdir -p legacy/{manual_fixes,old_docs,old_generators}
mkdir -p scripts
mkdir -p tests
```

## Phase 2: Move Analysis Scripts → framework/analysis/

| Current Location | New Location | Action |
|-----------------|--------------|---------|
| `analyze_work_statechange_frames.py` | `framework/analysis/work_calculator.py` | **MOVE** |
| `analyze_trajectory_target_yield_and_work.py` | `framework/analysis/trajectory_analyzer.py` | **MOVE** |
| `analyze_work_from_statechanges.py` | `framework/analysis/work_from_events.py` | **MOVE** |
| `analyze_yield_from_statechanges.py` | `framework/analysis/yield_from_events.py` | **MOVE** |
| `analyze_yield_over_time.py` | `framework/analysis/yield_timecourse.py` | **MOVE** |
| `analyze_yield_timecourse.py` | `framework/analysis/` | **MOVE** |
| `analyze_yields_from_dump.py` | `framework/analysis/` | **MOVE** |
| `check_state_changes.py` | `framework/analysis/check_state_changes.py` | **MOVE** |
| `change_types.py` | `legacy/old_generators/` | **ARCHIVE** |
| `change_types_simple.py` | `legacy/old_generators/` | **ARCHIVE** |

## Phase 3: Move Benchmark Code → framework/benchmark/

| Current Location | New Location | Action |
|-----------------|--------------|---------|
| `benchmark/run_task.py` | `framework/benchmark/task_runner.py` | **MOVE** |
| `benchmark/score_policy_from_timeseries.py` | `framework/benchmark/scorer.py` | **MOVE** |
| `benchmark/aggregate_leaderboard.py` | `framework/benchmark/leaderboard.py` | **MOVE** |
| `benchmark/task_schema.md` | `framework/benchmark/schemas/task_schema.md` | **MOVE** |
| `benchmark/README.md` | `docs/reference/benchmark.md` | **MOVE** |

## Phase 4: Move Build Scripts → framework/build/

| Current Location | New Location | Action |
|-----------------|--------------|---------|
| `add_new_fix.sh` | `framework/build/add_fix.sh` | **MOVE** |
| `rebuild_all_fixes.slurm` | `framework/build/rebuild_all.slurm` | **MOVE** |
| `REBUILD_COMMANDS.sh` | `framework/build/rebuild_commands.sh` | **MOVE** |

## Phase 5: Archive Manual Fixes → legacy/manual_fixes/

| Current Location | New Location | Action |
|-----------------|--------------|---------|
| `dimer/fix_state_change_dimer.cpp` | `legacy/manual_fixes/dimer/` | **ARCHIVE** |
| `dimer/fix_state_change_dimer.h` | `legacy/manual_fixes/dimer/` | **ARCHIVE** |
| `octahedron/fix_state_change_octahedron.cpp` | `legacy/manual_fixes/octahedron/` | **ARCHIVE** |
| `octahedron/fix_state_change_octahedron.h` | `legacy/manual_fixes/octahedron/` | **ARCHIVE** |
| `dimer_ksat/variants/*/fix_*.cpp` | `legacy/manual_fixes/dimer_ksat/variants/` | **ARCHIVE** |
| `dimer_ksat/variants/*/fix_*.h` | `legacy/manual_fixes/dimer_ksat/variants/` | **ARCHIVE** |
| `fix_state_change/` | `legacy/manual_fixes/fix_state_change/` | **ARCHIVE** |

## Phase 6: Archive Old Generators → legacy/old_generators/

| Current Location | New Location | Action |
|-----------------|--------------|---------|
| `generate_change.py` | `legacy/old_generators/` | **ARCHIVE** |
| `generate_change_cpp.py` | `legacy/old_generators/` | **ARCHIVE** |
| `dimer/generate_dimer_cpp.py` | `legacy/old_generators/dimer/` | **ARCHIVE** |
| `octahedron/generate_octahedron_cpp.py` | `legacy/old_generators/octahedron/` | **ARCHIVE** |
| `dimer_ksat/variants/*/generate.py` | `legacy/old_generators/dimer_ksat/variants/` | **ARCHIVE** |

## Phase 7: Consolidate Documentation

### MOVE to docs/user_guide/

| Current Location | New Location | Action |
|-----------------|--------------|---------|
| `README.md` | Keep, but rewrite as landing page | **EDIT** |
| Create new | `docs/user_guide/getting_started.md` | **CREATE** |
| Create new | `docs/user_guide/writing_policies.md` | **CREATE** |
| Create new | `docs/user_guide/policy_format.md` | **CREATE** |
| Create new | `docs/user_guide/scoring.md` | **CREATE** |

### MOVE to docs/developer_guide/

| Current Location | New Location | Action |
|-----------------|--------------|---------|
| `MULTIPLE_FIXES_GUIDE.md` | `docs/developer_guide/multiple_fixes.md` | **MOVE** |
| `REBUILD_INSTRUCTIONS.md` | `docs/developer_guide/building_lammps.md` | **MOVE** |
| Create new | `docs/developer_guide/adding_problems.md` | **CREATE** |
| Create new | `docs/developer_guide/codegen_internals.md` | **CREATE** |

### MOVE to docs/reference/

| Current Location | New Location | Action |
|-----------------|--------------|---------|
| `WORK_ANALYSIS_IMPROVEMENTS.md` | `docs/reference/work_analysis.md` | **MOVE** |
| `ANALYSIS_WORK_YIELD_FROM_TRAJECTORY.md` | `docs/reference/trajectory_analysis.md` | **MOVE** |
| `octahedron/FIX_SAME_TYPE_CONVERGENCE_BUG.md` | `docs/reference/bugfixes/octahedron_same_type.md` | **MOVE** |
| Create new | `docs/reference/api.md` | **CREATE** |
| Create new | `docs/reference/yaml_schemas.md` | **CREATE** |

### ARCHIVE to legacy/old_docs/

| Current Location | New Location | Action |
|-----------------|--------------|---------|
| `STATE_CHANGE_EXPLANATION.md` | `legacy/old_docs/` | **ARCHIVE** |
| `NEXT_STEPS_AFTER_REBUILD.md` | `legacy/old_docs/` | **ARCHIVE** |
| `octahedron/DIAGNOSTIC_ISSUES.md` | `legacy/old_docs/octahedron/` | **ARCHIVE** |
| `octahedron/MASS_BEHAVIOR_INVESTIGATION.md` | `legacy/old_docs/octahedron/` | **ARCHIVE** |
| `octahedron/MASS_BEHAVIOR_ANALYSIS.md` | `legacy/old_docs/octahedron/` | **ARCHIVE** |
| `octahedron/CONSISTENCY_SWEEP_FIX.md` | `legacy/old_docs/octahedron/` | **ARCHIVE** |
| `octahedron/CONSISTENCY_SWEEP_DIAGNOSTICS.md` | `legacy/old_docs/octahedron/` | **ARCHIVE** |
| `octahedron/DIAGNOSTIC_COUNTERS_ADDED.md` | `legacy/old_docs/octahedron/` | **ARCHIVE** |
| `octahedron/NaN_DIAGNOSIS.md` | `legacy/old_docs/octahedron/` | **ARCHIVE** |
| `octahedron/LOGIC_EVOLUTION_EXPLANATION.md` | `legacy/old_docs/octahedron/` | **ARCHIVE** |
| `octahedron/FIXES_APPLIED.md` | `legacy/old_docs/octahedron/` | **ARCHIVE** |
| `octahedron/SIMULATION_ANALYSIS.md` | `legacy/old_docs/octahedron/` | **ARCHIVE** |
| `octahedron/STATE_SELECTION_ANALYSIS.md` | `legacy/old_docs/octahedron/` | **ARCHIVE** |
| `octahedron/INSTALL_FIX.md` | `legacy/old_docs/octahedron/` | **ARCHIVE** |

### DELETE (replaced by framework)

| Current Location | Reason | Action |
|-----------------|--------|---------|
| `example_work_analysis_workflow.sh` | Replaced by submit_policy.py | **DELETE** |
| `submit.slurm` | Auto-generated per submission | **DELETE** |
| `submit_continuous.slurm` | Auto-generated per submission | **DELETE** |
| `submit_cpp_fix.slurm` | Auto-generated per submission | **DELETE** |
| `dimer/submit_dimer.slurm` | Auto-generated per submission | **DELETE** |
| `octahedron/submit_octahedron.slurm` | Auto-generated per submission | **DELETE** |
| `dimer_ksat/variants/*/submit.slurm` | Auto-generated per submission | **DELETE** |

## Phase 8: Reorganize Problem Definitions

### Octahedron

| Current Location | New Location | Action |
|-----------------|--------------|---------|
| `octahedron/README.md` | `problems/octahedron/README.md` | **MOVE + EDIT** |
| Create new | `problems/octahedron/problem.yaml` | **CREATE** |
| Create new | `problems/octahedron/geometry.py` | **CREATE** |
| Create new | `problems/octahedron/baseline/policy.yaml` | **CREATE** |

### Dimer

| Current Location | New Location | Action |
|-----------------|--------------|---------|
| Create new | `problems/dimer/problem.yaml` | **CREATE** |
| Create new | `problems/dimer/geometry.py` | **CREATE** |
| Create new | `problems/dimer/README.md` | **CREATE** |
| Create new | `problems/dimer/baseline/policy.yaml` | **CREATE** |

### Dimer KSAT

| Current Location | New Location | Action |
|-----------------|--------------|---------|
| Create new | `problems/dimer_ksat/variants/1core/problem.yaml` | **CREATE** |
| Create new | `problems/dimer_ksat/variants/1core/geometry.py` | **CREATE** |
| Create new | `problems/dimer_ksat/variants/base/problem.yaml` | **CREATE** |
| Create new | `problems/dimer_ksat/variants/base/geometry.py` | **CREATE** |

## Phase 9: Update .gitignore

Add these patterns:

```gitignore
# Generated simulation files
submissions/**/generated/
submissions/**/results/
**/simulation_cpp/
**/simulation_continuous/
**/rigid_patchy_simulation*/

# Analysis outputs
*.leaderboard.csv
*_timeseries.csv
*_timeseries.png

# SLURM outputs
slurm-*.out
slurm-*.err

# Python cache
__pycache__/
*.pyc
.pytest_cache/

# Keep submission structure
!submissions/**/.gitkeep
!submissions/**/policy.yaml
```

## Automated Migration Script

```bash
#!/bin/bash
# migrate.sh - Automated file migration

set -e

echo "Starting repository reorganization..."

# Phase 1: Create directories
echo "Creating directory structure..."
mkdir -p framework/{generators,analysis,benchmark,build,validation}
mkdir -p framework/generators/templates
mkdir -p framework/benchmark/schemas
mkdir -p problems/{octahedron,dimer,dimer_ksat/variants/{1core,base}}
mkdir -p submissions/{octahedron,dimer}
mkdir -p docs/{user_guide,developer_guide,reference/examples,reference/bugfixes}
mkdir -p legacy/{manual_fixes,old_docs,old_generators}
mkdir -p scripts
mkdir -p tests

# Phase 2: Move analysis scripts
echo "Moving analysis scripts..."
git mv analyze_work_statechange_frames.py framework/analysis/work_calculator.py
git mv analyze_trajectory_target_yield_and_work.py framework/analysis/trajectory_analyzer.py
git mv analyze_work_from_statechanges.py framework/analysis/work_from_events.py
git mv analyze_yield_from_statechanges.py framework/analysis/yield_from_events.py
git mv analyze_yield_over_time.py framework/analysis/yield_timecourse.py
git mv analyze_yield_timecourse.py framework/analysis/
git mv analyze_yields_from_dump.py framework/analysis/
git mv check_state_changes.py framework/analysis/

# Phase 3: Move benchmark code
echo "Moving benchmark code..."
git mv benchmark/run_task.py framework/benchmark/task_runner.py
git mv benchmark/score_policy_from_timeseries.py framework/benchmark/scorer.py
git mv benchmark/aggregate_leaderboard.py framework/benchmark/leaderboard.py
git mv benchmark/task_schema.md framework/benchmark/schemas/
git mv benchmark/README.md docs/reference/benchmark.md
rmdir benchmark

# Phase 4: Move build scripts
echo "Moving build scripts..."
git mv add_new_fix.sh framework/build/add_fix.sh
git mv rebuild_all_fixes.slurm framework/build/
git mv REBUILD_COMMANDS.sh framework/build/rebuild_commands.sh

# Phase 5: Archive manual fixes
echo "Archiving manual fixes..."
mkdir -p legacy/manual_fixes/{dimer,octahedron,dimer_ksat/variants}
git mv dimer/fix_state_change_dimer.* legacy/manual_fixes/dimer/
git mv octahedron/fix_state_change_octahedron.* legacy/manual_fixes/octahedron/
# (repeat for all dimer_ksat variants)

# Phase 6: Archive old generators
echo "Archiving old generators..."
git mv generate_change.py legacy/old_generators/
git mv generate_change_cpp.py legacy/old_generators/
git mv change_types*.py legacy/old_generators/

# Phase 7: Move documentation
echo "Moving documentation..."
git mv WORK_ANALYSIS_IMPROVEMENTS.md docs/reference/work_analysis.md
git mv ANALYSIS_WORK_YIELD_FROM_TRAJECTORY.md docs/reference/trajectory_analysis.md
git mv MULTIPLE_FIXES_GUIDE.md docs/developer_guide/multiple_fixes.md
git mv REBUILD_INSTRUCTIONS.md docs/developer_guide/building_lammps.md

# Archive old docs
git mv STATE_CHANGE_EXPLANATION.md legacy/old_docs/
git mv NEXT_STEPS_AFTER_REBUILD.md legacy/old_docs/
mkdir -p legacy/old_docs/octahedron
git mv octahedron/DIAGNOSTIC_*.md legacy/old_docs/octahedron/
git mv octahedron/MASS_BEHAVIOR_*.md legacy/old_docs/octahedron/
git mv octahedron/CONSISTENCY_SWEEP_*.md legacy/old_docs/octahedron/
git mv octahedron/NaN_DIAGNOSIS.md legacy/old_docs/octahedron/
git mv octahedron/LOGIC_EVOLUTION_EXPLANATION.md legacy/old_docs/octahedron/
git mv octahedron/FIXES_APPLIED.md legacy/old_docs/octahedron/

# Keep FIX_SAME_TYPE_CONVERGENCE_BUG.md (recent, important)
git mv octahedron/FIX_SAME_TYPE_CONVERGENCE_BUG.md docs/reference/bugfixes/octahedron_same_type.md

# Phase 8: Delete obsolete files
echo "Deleting obsolete files..."
git rm example_work_analysis_workflow.sh
git rm submit*.slurm
git rm dimer/submit_dimer.slurm
git rm octahedron/submit_octahedron.slurm

# Phase 9: Create .gitkeep files
touch submissions/octahedron/.gitkeep
touch submissions/dimer/.gitkeep

# Phase 10: Update .gitignore
cat >> .gitignore << 'EOF'

# Generated simulation files
submissions/**/generated/
submissions/**/results/
**/simulation_cpp/
**/simulation_continuous/
**/rigid_patchy_simulation*/

# Analysis outputs
*.leaderboard.csv
*_timeseries.csv
*_timeseries.png

# SLURM outputs
slurm-*.out
slurm-*.err

# Python cache
__pycache__/
*.pyc
.pytest_cache/

# Keep submission structure
!submissions/**/.gitkeep
!submissions/**/policy.yaml
EOF

echo "Migration complete! Please review changes with 'git status'"
echo "Next steps:"
echo "1. Create new framework code (generators, validators)"
echo "2. Create problem.yaml files"
echo "3. Write new documentation"
echo "4. Test end-to-end workflow"
```

## Manual Review Required

After running the migration script, manually review:

1. **README.md** - Rewrite as landing page for competition
2. **docs/user_guide/** - Write new user-facing documentation
3. **problems/**/problem.yaml** - Create problem specifications
4. **problems/**/baseline/policy.yaml** - Migrate current fixes to YAML format

## Rollback Plan

If migration causes issues:

```bash
git reset --hard HEAD~1  # Undo last commit
git clean -fd             # Remove untracked files
```

Or create a backup branch first:

```bash
git checkout -b pre-reorganization-backup
git checkout testing
# ... do migration ...
```
