#!/bin/bash
#SBATCH --job-name=analyze_dimer
#SBATCH --account=bewl-delta-gpu
#SBATCH --partition=gpuA100x4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --output=slurm_analyze-%j.out
#SBATCH --error=slurm_analyze-%j.err

set -euo pipefail

echo "=== Running analysis tasks for dimer/ ==="
cd "$SLURM_SUBMIT_DIR"

set +u
source /u/lguttieres/miniforge3/bin/activate
conda activate jax
set -u

TASK_DIR="tasks"
RUNNER="../../../benchmark/run_task.py"

shopt -s nullglob
for task in "$TASK_DIR"/*.json; do
  echo "--- Task: $task ---"
  python3 "$RUNNER" --task "$task"
done

echo "✅ Analysis complete"


