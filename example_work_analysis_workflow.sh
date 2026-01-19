#!/bin/bash
# Example workflow for improved work analysis
# This script shows how to use the new instantaneous dE functionality

set -e  # Exit on error

echo "===================================================="
echo "Work Analysis Workflow Example"
echo "===================================================="
echo ""

# Step 1: Check if LAMMPS has the updated fix
echo "Step 1: Verifying LAMMPS fix..."
LAMMPS_BIN="/work/nvme/bewl/lguttieres/lammps_build/lammps/src/lmp_mpi"

if [ ! -f "$LAMMPS_BIN" ]; then
    echo "  Warning: LAMMPS binary not found at $LAMMPS_BIN"
    echo "  You'll need to rebuild LAMMPS with the modified fix"
else
    echo "  ✓ LAMMPS binary found"
fi
echo ""

# Step 2: Generate simulation (example for dimer)
echo "Step 2: Generating simulation input files..."
cd dimer
python generate_dimer_cpp.py \
    --n_total 16 \
    --c_total 0.001 \
    --kt 0.5 \
    --dump_every 1000 \
    --thermo_every 1000 \
    --run_steps 5000000 \
    --out_dir dimer_simulation_cpp
echo "  ✓ Created dimer_simulation_cpp/"
cd ..
echo ""

# Step 3: Run simulation (example with SLURM)
echo "Step 3: Submitting simulation..."
echo "  You would run:"
echo "    cd dimer && sbatch submit_dimer.slurm"
echo ""
echo "  This will create:"
echo "    - dimer_simulation_cpp/dump.dimer.lammpstrj (trajectory)"
echo "    - dimer_simulation_cpp/lammps_stdout.log (thermo output)"
echo "    - slurm-JOBID.err (STATECHANGE events with dE)"
echo ""

# Step 4: Analyze work (after simulation completes)
echo "Step 4: Analyzing work from state changes..."
echo "  After your simulation completes, run:"
echo ""
echo "    python analyze_work_statechange_frames.py \\"
echo "        --thermo dimer_simulation_cpp/lammps_stdout.log \\"
echo "        --events slurm-JOBID.err \\"
echo "        --out analysis/work_dimer"
echo ""
echo "  Output files:"
echo "    - analysis/work_dimer.csv (per-event work values)"
echo "    - analysis/work_dimer.png (work plots)"
echo ""

# Step 5: Inspect results
echo "Step 5: Inspecting results..."
echo "  Check if dE is being calculated:"
echo ""
echo "    grep 'STATECHANGE' slurm-JOBID.err | head -5"
echo ""
echo "  Expected format (with dE):"
echo "    STATECHANGE dimer: step 12500 mol 3 flipped 2->3 dE -2.134567890123456"
echo ""
echo "  Old format (without dE):"
echo "    STATECHANGE dimer: step 12500 mol 3 flipped 2->3"
echo ""

# Step 6: Compare with old analysis
echo "Step 6: Comparing with old analysis method..."
echo "  Run the old script for comparison:"
echo ""
echo "    python analyze_work_from_statechanges.py \\"
echo "        --thermo dimer_simulation_cpp/lammps_stdout.log \\"
echo "        --events slurm-JOBID.err \\"
echo "        --out analysis/work_dimer_old"
echo ""
echo "  The new script will automatically use instantaneous dE if available,"
echo "  while the old script may use coarser frame differences."
echo ""

echo "===================================================="
echo "Workflow Summary"
echo "===================================================="
echo ""
echo "Key Improvements:"
echo "  ✓ Instantaneous dE calculation (zero time window)"
echo "  ✓ No thermal fluctuations included"
echo "  ✓ Works with any thermo output frequency"
echo "  ✓ 15-decimal precision for accuracy"
echo ""
echo "Next Steps:"
echo "  1. Rebuild LAMMPS with modified fix"
echo "  2. Run a test simulation"
echo "  3. Use analyze_work_statechange_frames.py"
echo "  4. Compare results with old method"
echo ""
echo "Documentation:"
echo "  See WORK_ANALYSIS_IMPROVEMENTS.md for details"
echo ""
