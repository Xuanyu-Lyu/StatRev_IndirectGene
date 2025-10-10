#!/bin/bash
# SLURM job submission script for RDR regression analysis - TEST VERSION
# Based on submit_rdr_test_he.sh but using Python regression instead of GCTA
# This version has most conditions commented out for testing

# --- Cluster Specific Settings ---
#SBATCH --qos=preemptable
#SBATCH --chdir /projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/RDR
#SBATCH --exclude bmem-rico1

# --- Job Resource Settings ---
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4      # Reduced CPUs since Python regression is less parallel than GCTA
#SBATCH --mem=32G              # Reduced memory requirement for Python regression
#SBATCH --time=0-1:00:00      # Reduced time since Python regression is typically faster

# --- Define the array size ---
#SBATCH --array=1-4%2        # Test with 50 runs and max 10 concurrent jobs

#SBATCH --job-name=rdr_regression_test
#SBATCH --output=slurm_logs/rdr_regression_test_%A_%a.out
#SBATCH --error=slurm_logs/rdr_regression_test_%A_%a.err

# --- Start of Job Commands ---
# The "set -e" command has been removed to ensure the script continues even if one condition fails.
mkdir -p slurm_logs

# Load environment
module purge
source /curc/sw/anaconda3/latest
conda activate /projects/xuly4739/general_env

# Make worker scripts executable
chmod +x run_rdr_regression.sh
chmod +x run_rdr_regression.py

# --- Define the list of conditions to run ---
CONDITIONS=(
    #"phenoVT_phenoAM"
    #"socialVT_phenoAM"
    #"phenoVT_socialAM"
    #"phenoVT_geneticAM"
    #"t1pheVT_t2socVT_uniphenoAM"
    #"01_t1pheVTnoAM_t2socVTnoAM"
    #"02_t1noVTpheAM_t2noVTnoAM"
    #"03_t1noVTsocAM_t2noVTnoAM"
    #"04_t1noVTgenAM_t2noVTnoAM"
    "05_t1pheVTnoAM_t2socVTnoAM_PGSall"
    "06_t1noVTpheAM_t2pheVTpheAM_PGSall"
    "07_t1noVTsocAM_t2pheVTsocAM_PGSall"
    "08_t1noVTgenAM_t2pheVTgenAM_PGSall"
)

# --- Loop through each condition ---
for CONDITION_NAME in "${CONDITIONS[@]}"; do
    echo "================================================="
    echo "--- Starting condition: ${CONDITION_NAME} for Slurm Task ID ${SLURM_ARRAY_TASK_ID} ---"
    echo "================================================="

    BASE_SIM_DIR="/scratch/alpine/xuly4739/StatRev_IndirectGene/Data/Paper/${CONDITION_NAME}"
    RUN_FOLDER=$(find "${BASE_SIM_DIR}" -mindepth 1 -maxdepth 1 -type d | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
    
    # Define the final, permanent directory for results  
    FINAL_RESULTS_DIR="/projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/RDR_Regression_Results/${CONDITION_NAME}"

    if [ -z "${RUN_FOLDER}" ]; then
        echo "Error: Could not find a run folder for task ID ${SLURM_ARRAY_TASK_ID} in condition ${CONDITION_NAME}. Skipping to next condition."
        continue # This command skips to the next item in the loop
    fi

    echo "Run folder: ${RUN_FOLDER}"
    echo "Results directory: ${FINAL_RESULTS_DIR}"

    # --- Execute the worker script for the current condition ---
    ./run_rdr_regression.sh "${RUN_FOLDER}" "${FINAL_RESULTS_DIR}"

    # Check the exit code of the worker script to log success or failure
    if [ $? -ne 0 ]; then
        echo "WARNING: The RDR regression script failed for condition: ${CONDITION_NAME}. Continuing to the next condition."
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