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
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G                      # Keep memory high as each task is still large
#SBATCH --time=0-24:00:00               # 20-hour time limit is plenty for a 2-task test

# --- Define the array size for the test ---
#SBATCH --array=1-1000%50

#SBATCH --job-name=gcta_rdr_four_condition  # A specific name for this test job
#SBATCH --output=slurm_logs/gcta_rdr_four_%A_%a.out
#SBATCH --error=slurm_logs/gcta_rdr_four_%A_%a.err

#
# --- Start of Job Commands ---
set -e
mkdir -p slurm_logs




# Load environment
module purge
source /curc/sw/anaconda3/latest
conda activate /projects/xuly4739/general_env

# Make scripts executable
chmod +x run_gcta_rdr_v2.sh
chmod +x prepare_grm_noGCTA.py
#chmod +x prepare_combined_plink.py
#chmod +x partition_grm.py

# --- Map Slurm Task ID to an input folder ---
CONDITION_NAME="phenoVT_phenoAM"
BASE_SIM_DIR="/scratch/alpine/xuly4739/StatRev_IndirectGene/Data/ASHG_Final/${CONDITION_NAME}"
RUN_FOLDER=$(find ${BASE_SIM_DIR} -mindepth 1 -maxdepth 1 -type d | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Define the final, permanent directory for the test results
FINAL_RESULTS_DIR="/projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/RDR_Results/${CONDITION_NAME}"

if [ -z "${RUN_FOLDER}" ]; then
    echo "Error: Could not find a run folder for task ID ${SLURM_ARRAY_TASK_ID}. Exiting."
    exit 1
fi

# --- Execute the worker script ---
echo "--- Running test for Task ID ${SLURM_ARRAY_TASK_ID} on folder ${RUN_FOLDER} ---"
./run_gcta_rdr_v2.sh "${RUN_FOLDER}" "${FINAL_RESULTS_DIR}"

# run another condition
CONDITION_NAME="socialVT_phenoAM"

BASE_SIM_DIR="/scratch/alpine/xuly4739/StatRev_IndirectGene/Data/ASHG_Final/${CONDITION_NAME}"
RUN_FOLDER=$(find ${BASE_SIM_DIR} -mindepth 1 -maxdepth 1 -type d | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ -z "${RUN_FOLDER}" ]; then
    echo "Error: Could not find a run folder for task ID ${SLURM_ARRAY_TASK_ID}. Exiting."
    exit 1
fi
# Execute the worker script for the second condition
echo "--- Running test for Task ID ${SLURM_ARRAY_TASK_ID} on folder ${RUN_FOLDER} ---"
./run_gcta_rdr_v2.sh "${RUN_FOLDER}" "${FINAL_RESULTS_DIR}"
echo "--- Slurm Array Task ${SLURM_ARRAY_TASK_ID} finished for second condition. ---"
# --- End of Job Commands ---

# phenoVT_socialAM condition
CONDITION_NAME="phenoVT_socialAM"
BASE_SIM_DIR="/scratch/alpine/xuly4739/StatRev_IndirectGene/Data/ASHG_Final/${CONDITION_NAME}"
RUN_FOLDER=$(find ${BASE_SIM_DIR} -mindepth 1 -maxdepth 1 -type d | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ -z "${RUN_FOLDER}" ]; then
    echo "Error: Could not find a run folder for task ID ${SLURM_ARRAY_TASK_ID}. Exiting."
    exit 1 
fi

# Execute the worker script for the third condition
echo "--- Running test for Task ID ${SLURM_ARRAY_TASK_ID} on folder ${RUN_FOLDER} ---"
./run_gcta_rdr_v2.sh "${RUN_FOLDER}" "${FINAL_RESULTS_DIR}"
# --- End of Job Commands ---

# phenoVT_geneticAM condition
CONDITION_NAME="phenoVT_geneticAM"
BASE_SIM_DIR="/scratch/alpine/xuly4739/StatRev_IndirectGene/Data/ASHG_Final/${CONDITION_NAME}"
RUN_FOLDER=$(find ${BASE_SIM_DIR} -mindepth 1 -maxdepth 1 -type d | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ -z "${RUN_FOLDER}" ]; then
    echo "Error: Could not find a run folder for task ID ${SLURM_ARRAY_TASK_ID}. Exiting."
    exit 1
fi  
# Execute the worker script for the fourth condition
echo "--- Running test for Task ID ${SLURM_ARRAY_TASK_ID} on folder ${RUN_FOLDER} ---"
./run_gcta_rdr_v2.sh "${RUN_FOLDER}" "${FINAL_RESULTS_DIR}" 


# Print completion message
echo "All tasks for Slurm Array Task ID ${SLURM_ARRAY_TASK_ID} completed successfully."
echo "Job finished at: $(date)"
echo "Total duration: $(( $(date +%s) - ${SCRIPT_START_TIME} )) seconds"
echo "------------------------------------------------"