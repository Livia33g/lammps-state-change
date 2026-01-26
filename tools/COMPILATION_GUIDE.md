# Efficient Compilation and Processing Guide

## Overview

This guide explains how to efficiently process multiple submissions while minimizing resource usage and avoiding conflicts.

## Key Concepts

### 1. Policy Hash Reuse
**Problem**: Multiple submissions might use the same policy (or very similar ones).  
**Solution**: Calculate MD5 hash of `policy.json` and reuse compiled LAMMPS builds for identical policies.

**Benefit**: If 10 submissions use the same policy, you only compile once!

### 2. Isolated Build Directories
**Problem**: Compiling different fixes in the same LAMMPS directory causes conflicts.  
**Solution**: Each unique policy gets its own isolated build directory.

**Structure**:
```
.builds/problem-001-dimer-ksat/
├── policy_abc12345_build/    # Build for policy hash abc12345
│   ├── lammps_build/         # Isolated LAMMPS copy
│   ├── lmp_mpi              # Compiled binary
│   └── *.h, *.cpp           # Fix files
└── policy_def67890_build/    # Different policy
    └── ...
```

### 3. Automatic Cleanup
**Problem**: Build directories accumulate and waste disk space.  
**Solution**: Automatic cleanup after processing (optional).

## Processing Strategies

### Strategy 1: One-by-One (Minimal Resources)

**Best for**: Limited disk space, processing as submissions arrive

```bash
# Process one submission at a time
for dir in submissions/problem-001-dimer-ksat/*/; do
    tools/process_submission.sh "$dir"
    # Compile and run
    # Clean up immediately
done
```

**Pros**:
- Minimal disk usage
- No conflicts possible
- Can process as submissions arrive

**Cons**:
- Slower (no build reuse)
- More manual work

### Strategy 2: Batch with Shared Builds (Recommended)

**Best for**: Processing many submissions efficiently

```bash
# Process all submissions, reuse builds for same policies
tools/compile_and_run_batch.sh problem-001-dimer-ksat --compile-only
```

**How it works**:
1. Process all submissions (generate C++ files)
2. Calculate policy hash for each
3. Group by hash (same policy = same build)
4. Compile once per unique policy
5. Clean up build directories

**Pros**:
- Efficient (reuses compilations)
- Automatic cleanup
- Handles many submissions

**Cons**:
- Requires more disk space temporarily

### Strategy 3: Batch with Isolated Builds

**Best for**: Debugging, when you need to keep builds

```bash
# Process all, keep builds separate
tools/compile_and_run_batch.sh problem-001-dimer-ksat --compile-only --keep-builds
```

**Pros**:
- Safe (no conflicts)
- Can debug builds later

**Cons**:
- Uses more disk space
- Slower (no reuse)

## Usage Examples

### Example 1: Process All Submissions Efficiently

```bash
# Set LAMMPS directory
export LAMMPS_DIR=/path/to/lammps

# Process all submissions for a problem
tools/compile_and_run_batch.sh problem-001-dimer-ksat --compile-only

# This will:
# 1. Generate C++ files for all submissions
# 2. Group by policy hash
# 3. Compile once per unique policy
# 4. Clean up automatically
```

### Example 2: Process Specific Submissions

```bash
# Process only specific competitors
tools/compile_and_run_batch.sh problem-001-dimer-ksat \
    --username alice \
    --username bob \
    --compile-only
```

### Example 3: Process and Run Simulations

```bash
# Compile and prepare for running
tools/compile_and_run_batch.sh problem-001-dimer-ksat --run

# Then run simulations (manual or via SLURM)
# The script will show you the commands to run
```

### Example 4: Process As Submissions Arrive

```bash
# When you receive a new submission
tools/import_submission.sh submission_001_newuser.zip
tools/process_submission.sh submissions/problem-001-dimer-ksat/newuser/

# Compile just this one
tools/compile_and_run_batch.sh problem-001-dimer-ksat \
    --username newuser \
    --compile-only
```

## Build Directory Management

### Location
Builds are stored in: `.builds/{problem_id}/`

### Cleanup Options

**Automatic cleanup (default)**:
```bash
tools/compile_and_run_batch.sh problem-001-dimer-ksat --compile-only
# Builds are removed after processing
```

**Keep builds (for debugging)**:
```bash
tools/compile_and_run_batch.sh problem-001-dimer-ksat --compile-only --keep-builds
# Builds remain in .builds/
```

**Manual cleanup**:
```bash
rm -rf .builds/
```

### Disk Space Considerations

- **One LAMMPS build**: ~500MB - 2GB (depending on compilation)
- **With reuse**: 10 submissions with same policy = 1 build (not 10)
- **Without reuse**: 10 submissions = 10 builds = 5-20GB

**Recommendation**: Use batch-shared strategy to minimize disk usage.

## Policy Hash Reuse Details

The system calculates MD5 hash of `policy.json`:

```bash
# Same policy = same hash = same build
md5sum submissions/problem-001-dimer-ksat/alice/policy.json
# abc12345...

md5sum submissions/problem-001-dimer-ksat/bob/policy.json  
# abc12345...  # Same hash = reuse alice's build!
```

**When reuse happens**:
- ✅ Identical `policy.json` files
- ✅ Same policy name and parameters

**When separate builds are created**:
- ❌ Different `policy.json` files
- ❌ Different policy names
- ❌ Different parameters

## Workflow Recommendations

### For Regular Processing (Recommended)

```bash
# 1. Import all new submissions
for file in new_submissions/*.zip; do
    tools/import_submission.sh "$file"
done

# 2. Process all at once (efficient)
tools/compile_and_run_batch.sh problem-001-dimer-ksat --compile-only

# 3. Run simulations (via SLURM or manual)
# (simulations use the compiled binaries)

# 4. Analyze results
for dir in submissions/problem-001-dimer-ksat/*/; do
    tools/analyze_submission_results.sh "$dir"
done

# 5. Update leaderboard
tools/update_leaderboard.sh problem-001-dimer-ksat
```

### For Incremental Processing

```bash
# When new submission arrives
tools/import_submission.sh new_submission.zip
tools/process_submission.sh submissions/problem-001-dimer-ksat/newuser/

# Check if policy matches existing (will reuse if same)
tools/compile_and_run_batch.sh problem-001-dimer-ksat \
    --username newuser \
    --compile-only
```

## Troubleshooting

### Issue: "LAMMPS_DIR not set"
```bash
export LAMMPS_DIR=/path/to/your/lammps/source
```

### Issue: Build conflicts
- Use `--keep-builds` to debug
- Check that each policy gets its own build directory
- Verify policy hashes are different for different policies

### Issue: Out of disk space
- Use `--cleanup` (default) to remove builds after
- Process in smaller batches
- Use one-by-one strategy

### Issue: Want to reuse existing builds
- Builds are stored in `.builds/{problem_id}/`
- If you keep builds (`--keep-builds`), they'll be reused automatically
- Policy hash matching handles reuse

## Best Practices

1. **Use batch-shared strategy** for efficiency
2. **Set LAMMPS_DIR** environment variable once
3. **Process in batches** rather than one-by-one when possible
4. **Clean up regularly** (automatic with `--cleanup`)
5. **Check disk space** before large batches
6. **Verify policy hashes** if you suspect issues

## Summary

- ✅ **Policy hash reuse**: Same policies = one compilation
- ✅ **Isolated builds**: No conflicts between submissions
- ✅ **Automatic cleanup**: Saves disk space
- ✅ **Flexible strategies**: Choose based on your needs
- ✅ **Batch processing**: Efficient for many submissions

**Recommended workflow**: Use `compile_and_run_batch.sh` with default options for maximum efficiency.

