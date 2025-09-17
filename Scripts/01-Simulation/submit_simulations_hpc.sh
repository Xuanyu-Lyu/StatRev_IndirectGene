#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --job-name=large_sim_batch_unitest      # A name for your large-scale job
#SBATCH --nodes=1                       # Each task requires 1 node
#SBATCH --ntasks=1                      # Each task is 1 main Python process
#SBATCH --cpus-per-task=4               # Match REPLICATIONS_PER_SLURM_TASK
#SBATCH --mem=300G                      # Large memory for a large population size
#SBATCH --time=6:00:00               # Request 1 full day (D-HH:MM:SS) for safety
#SBATCH --chdir /projects/xuly4739/Py_Projects/StatRev_IndirectGene/Scripts/01-Simulation
#SBATCH --exclude bmem-rico1
#SBATCH --output=slurm_logs/sim_run_%A_%a.out
#SBATCH --error=slurm_logs/sim_run_%A_%a.err

# --- Define the total number of tasks for the array ---
# CORRECTED: (4 conditions * 1000 reps) / 4 reps_per_task = 1000 tasks
#SBATCH --array=1-1000%25

# --- Your Job's Commands ---

mkdir -p slurm_logs

echo "------------------------------------------------"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Slurm Array Job ID: $SLURM_ARRAY_JOB_ID"
echo "Slurm Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "CPUs allocated to this task: $SLURM_CPUS_PER_TASK"
echo "Running on host: $(hostname)"
echo "------------------------------------------------"

# Load necessary modules
module purge
source /curc/sw/anaconda3/latest
conda activate /projects/xuly4739/general_env


# Run your simulation runner script
# The script now handles generating unique summary filenames internally
echo "Starting Python runner script..."
python -u run_simulations.py

echo "Slurm Array Task $SLURM_ARRAY_TASK_ID finished."