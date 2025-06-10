#!/bin/bash
# --- Cluster Specific Settings ---
#SBATCH --qos=preemptable                             
#SBATCH --chdir /projects/xuly4739/Py_Projects/StatRev_IndirectGene/Scripts/Debug 
#SBATCH --exclude bmem-rico1                         

# --- Job Resource Settings ---
# These are tailored for a benchmark run.
# The job name is changed to distinguish it from your analysis runs.
# Time and memory should be sufficient for your LARGEST test case.
#SBATCH --job-name=rdr_benchmark        # A specific name for this benchmark job
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1                 # Correct for the sequential Python benchmark script
#SBATCH --mem=128G                        # Memory for the largest N x M test
#SBATCH --time=0-02:00:00                 # 2-hour time limit for the benchmark

# --- Log File Settings ---
#SBATCH --output=slurm_logs/benchmark_%A.out      # Specific log file for benchmark output
#SBATCH --error=slurm_logs/benchmark_%A.err       # Specific log file for benchmark errors

# --- Start of Job Commands ---

# Ensure a clean environment
module purge

# Load the same Conda environment you used for your other jobs
source /curc/sw/anaconda3/latest
conda activate /projects/xuly4739/general_env

# Create the log directory if it doesn't exist
mkdir -p slurm_logs

# Print job information to the log for easy debugging
echo "------------------------------------------------"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Running on host: $(hostname)"
echo "Working Directory: $(pwd)"
echo "------------------------------------------------"

# Force scientific libraries to be single-threaded to prevent hangs
echo "Forcing libraries to be single-threaded..."
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Run your Python benchmarking script
# Using python -u for unbuffered output to see prints immediately
echo "Starting benchmark script..."
python -u benchmark_rdr.py

echo "Benchmark job finished."