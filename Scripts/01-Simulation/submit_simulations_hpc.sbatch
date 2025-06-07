#!/bin/bash
#SBATCH --job-name=sim_batch_hpc        # A name for your job
#SBATCH --output=slurm_logs/sim_run_%A_%a.out  # Path to write stdout, %A is job ID, %a is array task ID
#SBATCH --error=slurm_logs/sim_run_%A_%a.err   # Path to write stderr
#SBATCH --nodes=1                       # Each task requires 1 node
#SBATCH --ntasks=1                      # Each task is 1 main Python process
#SBATCH --cpus-per-task=10              # <<< Request 10 CPUs for each task
#SBATCH --mem=16G                       # <<< Increased memory for running 10 sims at once
#SBATCH --time=0-08:00:00               # <<< Increased time limit for running 10 sims

# --- Define the total number of tasks for the array ---
# This is now (total replications) / (replications per task)
# (2 conditions * 100 reps/condition) / 10 reps/task = 20 tasks
# So the array will be indexed 1-20.
#SBATCH --array=1-20

# --- Your Job's Commands ---

# Create the log directory if it doesn't exist
mkdir -p slurm_logs

# Print some useful information to the output file
echo "------------------------------------------------"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Slurm Array Job ID: $SLURM_ARRAY_JOB_ID"
echo "Slurm Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "CPUs allocated to this task: $SLURM_CPUS_PER_TASK"
echo "Running on host: $(hostname)"
echo "------------------------------------------------"

# Load necessary modules for your environment (e.g., Python)
# module load python/3.9

# Activate your Python virtual environment if you have one
# source /path/to/your/virtual/environment/bin/activate

# Run your Python script. It will read the Slurm environment variables
# to determine which batch of 10 replications to run.
python run_slurm_batch_job.py

echo "Slurm Array Task $SLURM_ARRAY_TASK_ID finished."