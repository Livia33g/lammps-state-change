#!/bin/bash
#SBATCH --job-name=lammps_dimer_ksat_twins
#SBATCH --account=bewl-delta-gpu
#SBATCH --partition=gpuA100x4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --output=slurm_dimer_ksat_1core_twosideB_twins-%j.out
#SBATCH --error=slurm_dimer_ksat_1core_twosideB_twins-%j.err

set -euo pipefail

echo "=== Starting DIMER_KSAT twins (ABC + EFD) run ==="
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

SIM_DIR="simulation_cpp"
FORCE_REGENERATE=${FORCE_REGENERATE:-1}

if [[ "$FORCE_REGENERATE" == "1" || ! -d "$SIM_DIR" ]]; then
  rm -rf "$SIM_DIR"
  echo "--- Generating inputs ---"
  python3 generate.py
else
  echo "Reusing existing $SIM_DIR (set FORCE_REGENERATE=1 to rebuild)"
fi

cd "$SIM_DIR"
if [[ ! -f in.dimer_ksat_1core_twosideB_twins || ! -f data.dimer_ksat_1core_twosideB_twins ]]; then
  echo "Missing input files in $SIM_DIR"
  exit 1
fi

unset SLURM_MEM_PER_CPU SLURM_MEM_PER_GPU SLURM_MEM_PER_NODE 2>/dev/null || true
srun "$LAMMPS_EXE" -in in.dimer_ksat_1core_twosideB_twins 2>&1 | tee lammps_stdout.log

echo "✅ Dimer_ksat twins job complete. Outputs in $SIM_DIR"



