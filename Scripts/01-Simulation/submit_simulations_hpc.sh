#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --job-name=sim_batch_hpc        # A name for your job
#SBATCH --nodes=1                       # Each task requires 1 node
#SBATCH --ntasks=1                      # Each task is 1 main Python process
#SBATCH --cpus-per-task=10              # <<< Request 10 CPUs for each task
#SBATCH --mem=32G                       # <<< Increased memory for running 10 sims at once
#SBATCH --time=08:00:00               # <<< Increased time limit for running 10 sims
#SBATCH --chdir /projects/xuly4739/Py_Projects/StatRev_IndirectGene/Scripts/01-Simulation
#SBATCH --exclude bmem-rico1
#SBATCH -o %x.out%A
#SBATCH -e %x.err%A



# --- Define the total number of tasks for the array ---
# This is now (total replications) / (replications per task)
# (2 conditions * 100 reps/condition) / 10 reps/task = 20 tasks
# So the array will be indexed 1-20.
#SBATCH --array=1-2

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
source /curc/sw/anaconda3/latest
conda activate /projects/xuly4739/general_env

# Activate your Python virtual environment if you have one
# source /path/to/your/virtual/environment/bin/activate

# Run your Python script. It will read the Slurm environment variables
# to determine which batch of 10 replications to run.
python run_slurm_batch_job.py

echo "Slurm Array Task $SLURM_ARRAY_TASK_ID finished."