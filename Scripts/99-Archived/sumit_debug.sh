#!/bin/bash
#SBATCH --job-name=rdr_debug            # A name for this test job
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G                      # Keep a high memory request for the test
#SBATCH --time=00:15:00                 # 15 minutes is plenty for this test
#SBATCH --output=slurm_logs/debug_run_%A.out
#SBATCH --error=slurm_logs/debug_run_%A.err
#SBATCH --qos=preemptable               # Use a high-priority queue if available for faster feedback

# This command ensures that the script will exit immediately if any command fails
set -e

# --- Your Job's Commands ---
mkdir -p slurm_logs
echo "--- SLURM SCRIPT STARTED ---"

echo "Step 1: Loading Anaconda module..."
source /curc/sw/anaconda3/latest
echo "Step 2: Anaconda module loaded."

echo "Step 3: Activating Conda environment..."
conda activate /projects/xuly4739/general_env
echo "Step 4: Conda environment activated."

echo "Step 5: Checking Python location and version..."
which python
python --version

echo "Step 6: Forcing libraries to be single-threaded..."
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
echo "Step 7: Environment variables set."

echo "Step 8: LAUNCHING PYTHON SCRIPT..."
# Use python -u for unbuffered output, so we see prints immediately
python -u run_rdr_analysis_DEBUG.py
echo "Step 9: PYTHON SCRIPT FINISHED."