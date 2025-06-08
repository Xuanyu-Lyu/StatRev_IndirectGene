#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --job-name=rdr_analysis_bivar   # A name for your job
#SBATCH --nodes=1
#SBATCH --ntasks=1                        # This script is 1 main process
#SBATCH --cpus-per-task=1  
#SBATCH --mem=200G                         # Memory for loading and processing data
#SBATCH --time=0-06:00:00                 # *** MODIFIED: Increased time for running both traits ***
#SBATCH --chdir /projects/xuly4739/Py_Projects/StatRev_IndirectGene/Scripts/01-Simulation
#SBATCH --exclude bmem-rico1
#SBATCH --output=slurm_logs/rdr_processing_%A.out  # Path to write stdout
#SBATCH --error=slurm_logs/rdr_processing_%A.err   # Path to write stderr
module purge

# Load necessary modules for your environment (e.g., Python)
source /curc/sw/anaconda3/latest
conda activate /projects/xuly4739/general_env
# --- Your Job's Commands ---

mkdir -p slurm_logs

echo "------------------------------------------------"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Running on host: $(hostname)"
echo "CPUs available to this job: $SLURM_CPUS_PER_TASK"
echo "------------------------------------------------"

echo "Forcing libraries to be single-threaded..."
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Run your new Python processing script
python run_rdr_analysis.py

echo "RDR processing job finished."