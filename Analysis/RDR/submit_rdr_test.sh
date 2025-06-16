#!/bin/bash
# --- Cluster Specific Settings ---
#SBATCH --qos=preemptable
#SBATCH --chdir /projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/RDR
#SBATCH --exclude bmem-rico1

# --- Job Resource Settings ---
# The job name is changed to distinguish it from your main analysis runs.
# Time and memory are reduced for a quicker test.
#SBATCH --job-name=gcta_rdr      # A specific name for this test job
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G                      # Keep memory high as each task is still large
#SBATCH --time=0-01:00:00               # 1-hour time limit is plenty for a 2-task test

# --- Log File Settings ---
#SBATCH --output=slurm_logs/gcta_rdr_%A_%a.out
#SBATCH --error=slurm_logs/gcta_rdr_%A_%a.err

# --- Define the array size for the test ---
# *** MODIFIED: Only run the first 2 tasks of the array ***
#SBATCH --array=1-2

# --- Your Job's Commands ---
set -e
mkdir -p slurm_logs

# Load your environment
module purge
source /curc/sw/anaconda3/latest
conda activate /projects/xuly4739/general_env

# Make the worker script executable
chmod +x run_gcta_for_samplesizes.sh

# --- Map Slurm Task ID to an input folder ---
# This finds all 'run_*' directories and uses the SLURM_ARRAY_TASK_ID to pick one
# *** IMPORTANT: Make sure this path points to ONE of your condition folders ***
CONDITION_DIR="/scratch/alpine/xuly4739/StatRev_IndirectGene/Data/ASHG_Final/phenoVT_phenoAM"
RUN_FOLDER=$(find ${CONDITION_DIR} -mindepth 1 -maxdepth 1 -type d | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Define where all results from this job array will go
BASE_OUTPUT_DIR="/projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/RDR_Results/Tests/${SLURM_ARRAY_JOB_ID}_${CONDITION_DIR##*/}"

if [ -z "${RUN_FOLDER}" ]; then
    echo "Error: Could not find a run folder for task ID ${SLURM_ARRAY_TASK_ID}. Exiting."
    exit 1
fi

# --- Execute the worker script ---
./run_gcta_for_samplesizes.sh "${RUN_FOLDER}" "${BASE_OUTPUT_DIR}"

echo "Slurm Array Task ${SLURM_ARRAY_TASK_ID} finished processing all sample sizes."