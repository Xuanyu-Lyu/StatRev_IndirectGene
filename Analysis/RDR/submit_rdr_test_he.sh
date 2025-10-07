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
#SBATCH --time=0-12:00:00

# --- Define the array size ---
#SBATCH --array=1-30%2

#SBATCH --job-name=gcta_rdr_hereg_check_variance
#SBATCH --output=slurm_logs/gcta_rdr_hereg_check_variance_%A_%a.out
#SBATCH --error=slurm_logs/gcta_rdr_hereg_check_variance_%A_%a.err

# --- Start of Job Commands ---
# The "set -e" command has been removed to ensure the script continues even if one condition fails.
mkdir -p slurm_logs

# Load environment
module purge
source /curc/sw/anaconda3/latest
conda activate /projects/xuly4739/general_env

# Make worker scripts executable
chmod +x run_GCTA_rdr_HEreg.sh
chmod +x prepare_grm_noGCTA.py

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
    "06_t1noVTpheAM_t2noVTnoAM_PGSall"
    "07_t1noVTsocAM_t2noVTnoAM_PGSall"
    "08_t1noVTgenAM_t2noVTnoAM_PGSall"
)

# --- Loop through each condition ---
for CONDITION_NAME in "${CONDITIONS[@]}"; do
    echo "================================================="
    echo "--- Starting condition: ${CONDITION_NAME} for Slurm Task ID ${SLURM_ARRAY_TASK_ID} ---"
    echo "================================================="

    BASE_SIM_DIR="/scratch/alpine/xuly4739/StatRev_IndirectGene/Data/Paper/${CONDITION_NAME}"
    RUN_FOLDER=$(find "${BASE_SIM_DIR}" -mindepth 1 -maxdepth 1 -type d | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
    
    # Define the final, permanent directory for results
    FINAL_RESULTS_DIR="/projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/RDR_HEreg_Results/${CONDITION_NAME}"

    if [ -z "${RUN_FOLDER}" ]; then
        echo "Error: Could not find a run folder for task ID ${SLURM_ARRAY_TASK_ID} in condition ${CONDITION_NAME}. Skipping to next condition."
        continue # This command skips to the next item in the loop
    fi

    # --- Execute the worker script for the current condition ---
    ./run_GCTA_rdr_HEreg.sh "${RUN_FOLDER}" "${FINAL_RESULTS_DIR}"

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