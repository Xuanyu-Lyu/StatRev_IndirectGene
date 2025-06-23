#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --job-name=aggregate_results
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1                 # This is a sequential script, only needs 1 CPU
#SBATCH --mem=8G                          # Should be plenty of memory for this task
#SBATCH --time=0-01:00:00                 # 1 hour should be more than enough
#SBATCH --chdir /projects/xuly4739/Py_Projects/StatRev_IndirectGene/Scripts/01-Simulation
#SBATCH --output=slurm_logs/aggregation_%A.out
#SBATCH --error=slurm_logs/aggregation_%A.err

# --- Your Job's Commands ---
set -e
mkdir -p slurm_logs

echo "------------------------------------------------"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Running on host: $(hostname)"
echo "------------------------------------------------"

# Load your environment
module purge
source /curc/sw/anaconda3/latest
conda activate /projects/xuly4739/general_env

# Run the Python aggregation script
echo "Starting aggregation script..."
python -u extract_and_average_results.py

echo "Aggregation job finished."