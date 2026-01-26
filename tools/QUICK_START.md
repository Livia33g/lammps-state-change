# Quick Start: Efficient Submission Processing

## The Smart Way (Recommended)

**Process all submissions efficiently with automatic policy reuse:**

```bash
# 1. Set your LAMMPS directory (once)
export LAMMPS_DIR=/path/to/lammps

# 2. Process all submissions for a problem
tools/compile_and_run_batch.sh problem-001-dimer-ksat --compile-only

# This automatically:
# - Generates C++ files for all submissions
# - Groups by policy hash (same policy = one compilation)
# - Compiles in isolated directories
# - Cleans up after (saves disk space)
```

## How It Works

### Policy Hash Reuse
- **Same policy.json** → **Same hash** → **One compilation**
- If 10 submissions use identical policies, you compile once, not 10 times!

### Isolated Builds
- Each unique policy gets its own build directory
- No conflicts between different submissions
- Builds are in `.builds/` (automatically cleaned up)

### Example Scenario

```
Submissions:
- alice: policy_A.json
- bob:   policy_A.json  (same as alice)
- charlie: policy_B.json (different)

Result:
- Compile policy_A once → used by alice AND bob
- Compile policy_B once → used by charlie
- Total: 2 compilations (not 3!)
```

## Processing Options

### Option 1: Process All at Once (Best for Many Submissions)

```bash
tools/compile_and_run_batch.sh problem-001-dimer-ksat --compile-only
```

**When**: You have many submissions to process  
**Benefit**: Maximum efficiency, automatic cleanup

### Option 2: Process Specific Submissions

```bash
tools/compile_and_run_batch.sh problem-001-dimer-ksat \
    --username alice \
    --username bob \
    --compile-only
```

**When**: Only some submissions need processing  
**Benefit**: Faster, targeted processing

### Option 3: Process As They Arrive

```bash
# When new submission comes in
tools/import_submission.sh new_submission.zip
tools/process_submission.sh submissions/problem-001-dimer-ksat/newuser/
tools/compile_and_run_batch.sh problem-001-dimer-ksat \
    --username newuser \
    --compile-only
```

**When**: Submissions arrive incrementally  
**Benefit**: Will reuse build if policy matches existing

## What Gets Created

### Per Submission:
```
submissions/problem-001-dimer-ksat/username/
├── generated/          # C++ fix files (small, ~10KB)
└── simulation/         # LAMMPS input files (small, ~1KB)
```

### Per Unique Policy:
```
.builds/problem-001-dimer-ksat/
└── policy_abc12345_build/
    ├── lammps_build/   # Isolated LAMMPS copy (~500MB-2GB)
    └── lmp_mpi        # Compiled binary (~50MB)
```

**Note**: `.builds/` is automatically cleaned up (unless you use `--keep-builds`)

## Disk Space

- **One LAMMPS build**: ~500MB - 2GB
- **With reuse**: 10 same policies = 1 build (not 10!)
- **Without reuse**: 10 different policies = 10 builds

**With policy hash reuse, you save 90% disk space!**

## Complete Workflow

```bash
# 1. Import submissions
for file in new_submissions/*.zip; do
    tools/import_submission.sh "$file"
done

# 2. Process and compile (efficient batch)
export LAMMPS_DIR=/path/to/lammps
tools/compile_and_run_batch.sh problem-001-dimer-ksat --compile-only

# 3. Run simulations (manual or SLURM)
# (use the compiled binaries from .builds/)

# 4. Analyze results
for dir in submissions/problem-001-dimer-ksat/*/; do
    tools/analyze_submission_results.sh "$dir"
done

# 5. Update leaderboard
tools/update_leaderboard.sh problem-001-dimer-ksat
```

## Key Benefits

✅ **No conflicts**: Each policy gets isolated build  
✅ **Efficient**: Reuses compilations for same policies  
✅ **Clean**: Automatic cleanup saves disk space  
✅ **Flexible**: Process all at once or incrementally  
✅ **Smart**: Policy hash matching handles reuse automatically  

## Troubleshooting

**"LAMMPS_DIR not set"**
```bash
export LAMMPS_DIR=/path/to/your/lammps/source
```

**Want to keep builds for debugging?**
```bash
tools/compile_and_run_batch.sh problem-001-dimer-ksat --keep-builds
```

**Out of disk space?**
- Use `--cleanup` (default) to remove builds after
- Process in smaller batches
- Check `.builds/` directory size

## Summary

**Recommended approach**: Use `compile_and_run_batch.sh` with default options.

- ✅ Processes all submissions efficiently
- ✅ Reuses builds for same policies (saves time & space)
- ✅ Isolated builds (no conflicts)
- ✅ Automatic cleanup (no clutter)
- ✅ Minimal resource usage

**No virtual environments needed** - isolated build directories handle everything!

