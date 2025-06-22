#!/bin/bash
# --- Cluster Specific Settings ---
#SBATCH --qos=preemptable
#SBATCH --chdir /projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/RDR
#SBATCH --exclude bmem-rico1

# --- Job Resource Settings ---
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=140G
#SBATCH --time=0-24:00:00

# --- Define the array size ---
#SBATCH --array=1-3

#SBATCH --job-name=gcta_rdr_check_variance
#SBATCH --output=slurm_logs/gcta_rdr_check_variance_%A_%a.out
#SBATCH --error=slurm_logs/gcta_rdr_check_variance_%A_%a.err

# --- Start of Job Commands ---
# The "set -e" command has been removed to ensure the script continues even if one condition fails.
mkdir -p slurm_logs

# Load environment
module purge
source /curc/sw/anaconda3/latest
conda activate /projects/xuly4739/general_env

# Make worker scripts executable
chmod +x run_gcta_rdr_v2.sh
chmod +x prepare_grm_noGCTA.py

# --- Define the list of conditions to run ---
CONDITIONS=(
    "phenoVT_phenoAM"
    #"socialVT_phenoAM"
    #"phenoVT_socialAM"
    #"phenoVT_geneticAM"
)

# --- Loop through each condition ---
for CONDITION_NAME in "${CONDITIONS[@]}"; do
    echo "================================================="
    echo "--- Starting condition: ${CONDITION_NAME} for Slurm Task ID ${SLURM_ARRAY_TASK_ID} ---"
    echo "================================================="

    BASE_SIM_DIR="/scratch/alpine/xuly4739/StatRev_IndirectGene/Data/ASHG_Final/Tests/${CONDITION_NAME}"
    RUN_FOLDER=$(find "${BASE_SIM_DIR}" -mindepth 1 -maxdepth 1 -type d | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
    
    # Define the final, permanent directory for results
    FINAL_RESULTS_DIR="/projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/RDR_Results/${CONDITION_NAME}"

    if [ -z "${RUN_FOLDER}" ]; then
        echo "Error: Could not find a run folder for task ID ${SLURM_ARRAY_TASK_ID} in condition ${CONDITION_NAME}. Skipping to next condition."
        continue # This command skips to the next item in the loop
    fi

    # --- Execute the worker script for the current condition ---
    ./run_gcta_rdr_v2.sh "${RUN_FOLDER}" "${FINAL_RESULTS_DIR}"

    # Check the exit code of the worker script to log success or failure
    if [ $? -ne 0 ]; then
        echo "WARNING: The worker script failed for condition: ${CONDITION_NAME}. Continuing to the next condition."
    else
        echo "--- Successfully finished condition: ${CONDITION_NAME} ---"
    fi
    echo
done

# --- Final completion message ---
echo "================================================="
echo "All conditions for Slurm Array Task ID ${SLURM_ARRAY_TASK_ID} have been attempted."
echo "Job finished at: $(date)"
echo "------------------------------------------------"