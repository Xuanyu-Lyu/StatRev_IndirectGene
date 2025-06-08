#!/bin/bash
#SBATCH --job-name=rdr_analysis_bivar   # A name for your job
#SBATCH --nodes=1                         # Request 1 node
#SBATCH --ntasks=1                        # This script is 1 main process
#SBATCH --cpus-per-task=5                # Request 20 CPUs for multiprocessing
#SBATCH --mem=100G                         # Memory for loading and processing data
#SBATCH --time=0-04:00:00                 # *** MODIFIED: Increased time for running both traits ***
#SBATCH --output=slurm_logs/rdr_processing_%A.out  # Path to write stdout
#SBATCH --error=slurm_logs/rdr_processing_%A.err   # Path to write stderr

# --- Your Job's Commands ---

mkdir -p slurm_logs

echo "------------------------------------------------"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Running on host: $(hostname)"
echo "CPUs available to this job: $SLURM_CPUS_PER_TASK"
echo "------------------------------------------------"

# Load the same Conda environment you used for the simulations
source /curc/sw/anaconda3/latest
conda activate /projects/xuly4739/general_env

# Run your new Python processing script
python run_rdr_analysis.py

echo "RDR processing job finished."