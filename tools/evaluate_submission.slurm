#!/bin/bash
#SBATCH --job-name=eval_submission
#SBATCH --output=eval_%j.out
#SBATCH --error=eval_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=4
#SBATCH --mem=8G

# Automated submission evaluation script
# Usage: sbatch tools/evaluate_submission.slurm problem-001 username

PROBLEM=$1
USERNAME=$2
SUBMISSION_DIR="submissions/${PROBLEM}/${USERNAME}"

if [ ! -d "$SUBMISSION_DIR" ]; then
    echo "Error: Submission directory not found: $SUBMISSION_DIR"
    exit 1
fi

echo "=== Evaluating submission: ${USERNAME} for ${PROBLEM} ==="

# Step 1: Validate submission
echo "Step 1: Validating submission..."
python3 tools/validate_submission.py "$SUBMISSION_DIR"
if [ $? -ne 0 ]; then
    echo "FAILED: Validation errors"
    exit 1
fi

# Step 2: Generate C++ fix
echo "Step 2: Generating C++ fix..."
python3 core/generators/generate_fix_from_policy.py \
    "${SUBMISSION_DIR}/policy.json" \
    "${SUBMISSION_DIR}/generated/"

# Step 3: Generate LAMMPS files
echo "Step 3: Generating LAMMPS simulation files..."
python3 core/generators/generate_system_from_problem.py \
    --problem "problems/${PROBLEM}/problem.json" \
    --policy "${SUBMISSION_DIR}/policy.json" \
    --params "${SUBMISSION_DIR}/params.json" \
    --output "${SUBMISSION_DIR}/generated/"

# Step 4: Compile LAMMPS with custom fix
echo "Step 4: Compiling LAMMPS..."
FIX_NAME=$(ls ${SUBMISSION_DIR}/generated/fix_state_change_*.cpp | xargs basename | sed 's/.cpp$//')
cp ${SUBMISSION_DIR}/generated/fix_state_change_*.{h,cpp} $LAMMPS_SRC/
cd $LAMMPS_SRC
make mpi -j8
if [ $? -ne 0 ]; then
    echo "FAILED: LAMMPS compilation"
    exit 1
fi

# Step 5: Run simulation
echo "Step 5: Running LAMMPS simulation..."
cd ${SUBMISSION_DIR}/generated/
mpirun -np 4 $LAMMPS_SRC/../src/lmp_mpi \
    -in in.${PROBLEM/problem-/}.lammpstrj \
    > lammps_stdout.log 2> stderr.log

# Step 6: Analyze results
echo "Step 6: Analyzing results..."
python3 ../../core/benchmark/score_policy_from_timeseries.py \
    --dump dump.*.lammpstrj \
    --thermo lammps_stdout.log \
    --events stderr.log \
    --output analysis/

# Step 7: Update leaderboard
echo "Step 7: Updating leaderboard..."
python3 ../../core/benchmark/update_leaderboard.py \
    --problem ${PROBLEM} \
    --username ${USERNAME} \
    --results analysis/scores.json

echo "=== Evaluation complete! ==="
echo "Results: ${SUBMISSION_DIR}/generated/analysis/"
