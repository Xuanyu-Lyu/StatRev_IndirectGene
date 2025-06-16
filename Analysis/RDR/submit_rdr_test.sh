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

# --- Start of Job Commands ---
set -e
mkdir -p slurm_logs

# Load environment
module purge
source /curc/sw/anaconda3/latest
conda activate /projects/xuly4739/general_env

# Make scripts executable
chmod +x run_gcta_rdr.sh
chmod +x prepare_combined_plink.py
chmod +x partition_grm.py

# --- Map Slurm Task ID to an input folder ---
CONDITION_NAME="phenoVT_phenoAM"
BASE_SIM_DIR="/scratch/alpine/xuly4739/StatRev_IndirectGene/Data/ASHG_Final/${CONDITION_NAME}"
RUN_FOLDER=$(find ${BASE_SIM_DIR} -mindepth 1 -maxdepth 1 -type d | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Define the final, permanent directory for the test results
FINAL_RESULTS_DIR="/projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/RDR_Results/Tests/${CONDITION_NAME}"

if [ -z "${RUN_FOLDER}" ]; then
    echo "Error: Could not find a run folder for task ID ${SLURM_ARRAY_TASK_ID}. Exiting."
    exit 1
fi

# --- Execute the worker script ---
echo "--- Running test for Task ID ${SLURM_ARRAY_TASK_ID} on folder ${RUN_FOLDER} ---"
./run_gcta_rdr.sh "${RUN_FOLDER}" "${FINAL_RESULTS_DIR}"

echo "--- Slurm Array Task ${SLURM_ARRAY_TASK_ID} finished. ---"