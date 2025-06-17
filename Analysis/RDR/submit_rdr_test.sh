#!/bin/bash
# --- Cluster Specific Settings ---
#SBATCH --qos=preemptable
#SBATCH --chdir /projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/RDR
#SBATCH --exclude bmem-rico1

# --- Job Resource Settings ---
# The job name is changed to distinguish it from your main analysis runs.
# Time and memory are reduced for a quicker test.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=100G                      # Keep memory high as each task is still large
#SBATCH --time=0-20:00:00               # 20-hour time limit is plenty for a 2-task test

# --- Define the array size for the test ---
# *** MODIFIED: Only run the first 2 tasks of the array ***
#SBATCH --array=1-1000%25

# --- Dynamic Job Name and Log File Settings ---
# The job name and log files are set dynamically using the CONDITION_NAME variable.
# This prevents outputs from different conditions from overwriting each other.
#
# *** MODIFIED: Using ${CONDITION_NAME} variable for unique naming ***
#SBATCH --job-name=gcta_rdr_${CONDITION_NAME}
#SBATCH --output=slurm_logs/gcta_rdr_${CONDITION_NAME}_%A_%a.out
#SBATCH --error=slurm_logs/gcta_rdr_${CONDITION_NAME}_%A_%a.err

# --- How to Submit the Job ---
#
# You must still provide the condition name when you submit the job.
#
# In your terminal, run this command:
#   sbatch --export=CONDITION_NAME="your_condition_name_here" your_script_name.sh
#
# --- Start of Job Commands ---
set -e
mkdir -p slurm_logs

# --- Validate Input: Check if CONDITION_NAME is set ---
if [ -z "${CONDITION_NAME}" ]; then
    echo "Error: The CONDITION_NAME environment variable is not set."
    echo "Please set it when submitting the job with sbatch."
    echo 'Example: sbatch --export=CONDITION_NAME="phenoVT_phenoAM"' "$0"
    exit 1
fi



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