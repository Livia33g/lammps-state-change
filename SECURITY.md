# Security and Privacy: Submissions

## Submission Privacy

**All competitor submissions are kept private and local-only.**

### What is NOT pushed to remote:

- ✅ **All submission files** (`submissions/*/`) - Completely gitignored
- ✅ **All generated files** (`submissions/*/generated/`, `simulations/`, `results/`)
- ✅ **All test outputs** (`tests/`)
- ✅ **All competitor policy.json, params.json, submission.json files**

### What IS pushed to remote:

- ✅ **Problem definitions** (`problems/*/problem.json`)
- ✅ **Baseline policies** (`problems/*/baseline_policy.json`) - Reference solutions only
- ✅ **Analysis scripts** (`problems/*/analyze_submission.sh`)
- ✅ **Final leaderboards** (`problems/*/leaderboard.csv`) - Aggregated scores only
- ✅ **Framework code** (`core/`, `tools/`)
- ✅ **Documentation** (`README.md`, `PARTICIPANT_GUIDE.md`, etc.)

## How It Works

### For Competitors (Remote Users)

When competitors clone the repository, they see:
- Problem definitions and requirements
- Baseline policies (reference solutions)
- Framework documentation
- Starter templates

They **DO NOT** see:
- Other competitors' submissions
- Other competitors' policies
- Other competitors' results
- Generated C++ code from other competitors

### For Regulator (You)

You have access to:
- All submissions locally (in `submissions/` directory)
- All generated files
- All simulation results
- All analysis outputs

You can:
- Process submissions locally
- Run simulations
- Generate leaderboards
- Push only aggregated leaderboards (scores only, no policies)

## Gitignore Protection

The `.gitignore` file ensures:

```gitignore
# All submissions are ignored
submissions/
!submissions/.gitkeep          # Only directory structure tracked
!submissions/README.md          # Only documentation tracked
!submissions/**/EXAMPLE_STRUCTURE.md  # Only example structure tracked

# All generated files are ignored
submissions/**/generated/
submissions/**/simulation/
submissions/**/results/
```

## Verification

To verify submissions are not tracked:

```bash
# Check what git sees
git status

# Submissions should NOT appear (except .gitkeep and README.md)
# Only these should be visible:
# - submissions/.gitkeep
# - submissions/README.md
# - submissions/problem-XXX/EXAMPLE_STRUCTURE.md
```

## Leaderboard Privacy

The final `leaderboard.csv` contains:
- ✅ Username/team name
- ✅ Scores and metrics
- ✅ Submission dates
- ❌ **NO policy details** (no policy.json contents)
- ❌ **NO parameter values** (no params.json contents)

Only aggregated results are shared, not the actual submission strategies.

## Best Practices

1. **Never commit submission files**: Always verify with `git status` before committing
2. **Only push leaderboards**: After aggregating results, only push the CSV files
3. **Keep submissions local**: All competitor work stays on your machine
4. **Review before pushing**: Check `git diff` to ensure no submission data is included

## If You Need to Share Submissions

If you need to share specific submissions (e.g., for research):

1. **Create a separate branch**: `git checkout -b shared-submissions`
2. **Selectively add**: Only add the specific submission you want to share
3. **Get permission**: Ensure competitor has given permission
4. **Document**: Add a note explaining why this submission is public

**Default behavior**: All submissions remain private and local-only.

