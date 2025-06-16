#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --job-name=process_sims         # A name for your job
#SBATCH --nodes=1                       # Request 1 node
#SBATCH --ntasks=1                      # This script is 1 main process
#SBATCH --cpus-per-task=10              # Number of CPUs allocated to this task
#SBATCH --mem=200G                       # Memory for loading and processing data
#SBATCH --time=0-12:00:00               # Max job time: D-HH:MM:SS (e.g., 2 hours)
#SBATCH --chdir /projects/xuly4739/Py_Projects/StatRev_IndirectGene/Scripts/01-Simulation
#SBATCH --exclude bmem-rico1
#SBATCH --output=slurm_logs/processing_%A.out  # Path to write stdout (%A is job ID)
#SBATCH --error=slurm_logs/processing_%A.err   # Path to write stderr

# --- Your Job's Commands ---

# Create the log directory if it doesn't exist
mkdir -p slurm_logs

# Print job information
echo "------------------------------------------------"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Running on host: $(hostname)"
echo "CPUs available to this job: $SLURM_CPUS_PER_TASK"
echo "------------------------------------------------"

# Load the same Conda environment you used for the simulations
source /curc/sw/anaconda3/latest
conda activate /projects/xuly4739/general_env

# Run your new Python processing script
python process_simulation_results.py

echo "Processing job finished."