#!/bin/bash
"""
Shell script to run RDR regression analysis on a single simulation run
This script is similar to run_gcta_rdr_he.sh but uses Python regression instead of GCTA
"""

set -e 

if [[ -z "$1" || -z "$2" ]]; then
    echo "Usage: ./run_rdr_regression.sh <path_to_run_folder_on_scratch> <path_to_final_results_dir_on_projects>"
    exit 1
fi

RUN_FOLDER=$1
FINAL_RESULTS_DIR=$2
RUN_ID=$(basename ${RUN_FOLDER}) 

TARGET_SAMPLE_SIZES=(8000)

echo "--- Starting RDR regression analysis for ${RUN_ID} using Python ---"
echo "Analysis started at: $(date)"
SCRIPT_START_TIME=$(date +%s)

# Loop through each target sample size
for N in "${TARGET_SAMPLE_SIZES[@]}"; do
    echo -e "\n================== Processing for N=${N} =================="
    echo "Processing N=${N} started at: $(date)"
    N_START_TIME=$(date +%s)
    
    OUTPUT_DIR="${FINAL_RESULTS_DIR}/${RUN_ID}"
    
    mkdir -p ${OUTPUT_DIR}
    mkdir -p ${FINAL_RESULTS_DIR}

    # Run RDR regression analysis using Python
    echo "Running RDR regression analysis for N=${N}..."
    echo "RDR regression started at: $(date)"
    RDR_START_TIME=$(date +%s)
    
    python run_rdr_regression.py ${RUN_FOLDER} ${OUTPUT_DIR} ${N}
    
    RDR_END_TIME=$(date +%s)
    RDR_DURATION=$((RDR_END_TIME - RDR_START_TIME))
    echo "RDR regression completed at: $(date) (Duration: ${RDR_DURATION}s)"

    N_END_TIME=$(date +%s)
    N_DURATION=$((N_END_TIME - N_START_TIME))
    echo "--- Finished RDR regression analysis for N=${N} at: $(date) (Total duration: ${N_DURATION}s) ---"
done

SCRIPT_END_TIME=$(date +%s)
SCRIPT_DURATION=$((SCRIPT_END_TIME - SCRIPT_START_TIME))
echo "--- All sample sizes for ${RUN_ID} complete at: $(date) ---"
echo "Total script duration: ${SCRIPT_DURATION}s"