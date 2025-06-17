#!/bin/bash
# --- Cluster Specific Settings (Preserved from your original script) ---
#SBATCH --qos=preemptable
#SBATCH --chdir /projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/RDR
#SBATCH --exclude bmem-rico1

# --- Job Resource Settings ---
# Time and memory are reduced as this aggregation script is much less demanding
# than the main analysis jobs.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2      # Reduced CPUs, as aggregation is not highly parallel
#SBATCH --mem=8G               # Reduced memory, 8G is sufficient for reading text files
#SBATCH --time=0-01:00:00      # 1-hour time limit is ample for this task

# --- Job Naming and Log Files ---
# A specific name and log location for the aggregation job.
# Note: The %A_%a format is for job arrays; for a single job, %j is simpler.
# Using %x for job name, %j for job ID.
#SBATCH --job-name=aggregate_rdr
#SBATCH --output=slurm_logs/%x_%j.out
#SBATCH --error=slurm_logs/%x_%j.err

#
# --- Start of Job Commands ---
set -e
mkdir -p slurm_logs

# Load the same conda environment
echo "--- Loading Conda Environment ---"
module purge
source /curc/sw/anaconda3/latest
conda activate /projects/xuly4739/general_env
echo "--- Environment Loaded ---"

# --- Execute the aggregation script ---
# The python script uses absolute paths, so it can be run from here.
echo "--- Starting GCTA results aggregation ---"
python aggregate_gcta_results.py
echo "--- Aggregation script finished successfully. ---"