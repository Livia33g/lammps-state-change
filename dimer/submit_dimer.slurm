#!/bin/bash
#SBATCH --job-name=lammps_dimer
#SBATCH --account=bewl-delta-gpu
#SBATCH --partition=gpuA100x4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --output=slurm_dimer-%j.out
#SBATCH --error=slurm_dimer-%j.err

set -euo pipefail

echo "=== Starting DIMER state-change run ==="
cd "$SLURM_SUBMIT_DIR"
echo "Working directory: $(pwd)"

module purge
module load gcc/11.4.0 2>/dev/null || true
module load openmpi/4.1.6 2>/dev/null || true
module load libfabric/1.22.0 2>/dev/null || true

set +u
source /u/lguttieres/miniforge3/bin/activate
conda activate py-env
set -u

export PATH=/work/nvme/bewl/lguttieres/lammps_build/lammps/src:$PATH
LAMMPS_EXE=/work/nvme/bewl/lguttieres/lammps_build/lammps/src/lmp_mpi
echo "Using LAMMPS: $LAMMPS_EXE"

SIM_DIR="dimer_simulation_cpp"
FORCE_REGENERATE=${FORCE_REGENERATE:-1}

if [[ "$FORCE_REGENERATE" == "1" || ! -d "$SIM_DIR" ]]; then
  rm -rf "$SIM_DIR"
  echo "--- Generating dimer inputs ---"
  # Parameters (override at submit time, e.g. RUN_STEPS=50000000 sbatch submit_dimer.slurm)
  KT=${KT:-0.5}
  C_TOTAL=${C_TOTAL:-0.001}
  N_TOTAL=${N_TOTAL:-16}          # even for 50/50
  RED_FRAC=${RED_FRAC:-0.5}
  RUN_STEPS=${RUN_STEPS:-20000000}
  DUMP_EVERY=${DUMP_EVERY:-1000}
  THERMO_EVERY=${THERMO_EVERY:-1000}
  D23=${D23:-6.0}
  D33=${D33:-12.0}
  SEED=${SEED:-1234}

  python3 generate_dimer_cpp.py \
    --out_dir "$SIM_DIR" \
    --kt "$KT" \
    --c_total "$C_TOTAL" \
    --n_total "$N_TOTAL" \
    --red_fraction "$RED_FRAC" \
    --run_steps "$RUN_STEPS" \
    --dump_every "$DUMP_EVERY" \
    --thermo_every "$THERMO_EVERY" \
    --D23 "$D23" \
    --D33 "$D33" \
    --seed "$SEED"
else
  echo "Reusing existing $SIM_DIR (set FORCE_REGENERATE=1 to rebuild)"
fi

cd "$SIM_DIR"
if [[ ! -f in.dimer || ! -f data.dimer ]]; then
  echo "Missing input files in $SIM_DIR"
  exit 1
fi

# Run
unset SLURM_MEM_PER_CPU SLURM_MEM_PER_GPU SLURM_MEM_PER_NODE 2>/dev/null || true
srun "$LAMMPS_EXE" -in in.dimer 2>&1 | tee lammps_stdout.log

echo "✅ Dimer job complete. Outputs in $SIM_DIR"





