#!/bin/bash
set -e 

if [[ -z "$1" || -z "$2" ]]; then
    echo "Usage: ./run_gcta_rdr.sh <path_to_run_folder_on_scratch> <path_to_final_results_dir_on_projects>"
    exit 1
fi

RUN_FOLDER=$1
FINAL_RESULTS_DIR=$2
RUN_ID=$(basename ${RUN_FOLDER}) 

# Make the new script executable
chmod +x prepare_grm_noGCTA.py

TARGET_SAMPLE_SIZES=(2000 4000 8000 16000 32000)

echo "--- Starting RDR analysis for ${RUN_ID} using custom GRM calculation ---"
echo "Analysis started at: $(date)"
SCRIPT_START_TIME=$(date +%s)

# Loop through each target sample size
for N in "${TARGET_SAMPLE_SIZES[@]}"; do
    echo -e "\n================== Processing for N=${N} =================="
    echo "Processing N=${N} started at: $(date)"
    N_START_TIME=$(date +%s)
    
    WORK_DIR="${RUN_FOLDER}/gcta_analysis/N_${N}"
    OUTPUT_PREFIX="${FINAL_RESULTS_DIR}/${RUN_ID}_N${N}"
    
    mkdir -p ${WORK_DIR}
    mkdir -p ${FINAL_RESULTS_DIR}

    # Step 1: Prepare inputs AND calculate all three GRMs using the correct z-score method
    echo "Step 1: Calculating custom RDR GRMs for N=${N}..."
    echo "Step 1 started at: $(date)"
    STEP1_START_TIME=$(date +%s)
    
    python prepare_grm_noGCTA.py ${RUN_FOLDER} ${WORK_DIR} ${N}
    
    STEP1_END_TIME=$(date +%s)
    STEP1_DURATION=$((STEP1_END_TIME - STEP1_START_TIME))
    echo "Step 1 completed at: $(date) (Duration: ${STEP1_DURATION}s)"

    # Step 2 : Create the multi-GRM input file (mgrm.txt)
    echo "Step 2: Creating multi-GRM file for N=${N}..."
    echo "Step 2 started at: $(date)"
    STEP2_START_TIME=$(date +%s)
    
    printf "%s\n" \
    "${WORK_DIR}/rdr_grm_Ro_offspring" \
    "${WORK_DIR}/rdr_grm_Rp_parental" \
    "${WORK_DIR}/rdr_grm_Rop_cross" \
    > "${WORK_DIR}/mgrm.txt"
    
    STEP2_END_TIME=$(date +%s)
    STEP2_DURATION=$((STEP2_END_TIME - STEP2_START_TIME))
    echo "Step 2 completed at: $(date) (Duration: ${STEP2_DURATION}s)"

    # Step 3: Run Univariate GREML Analysis with the correctly calculated GRMs
    echo "Step 3: Running RDR GREML analysis for Trait 1 (Y1) with N=${N}..."
    echo "Step 3a (Y1) started at: $(date)"
    STEP3A_START_TIME=$(date +%s)
    
    gcta64 --reml --mgrm "${WORK_DIR}/mgrm.txt" \
           --reml-no-lrt \
           --reml-no-constrain \
           --pheno "${WORK_DIR}/offspring.phen" --mpheno 1 \
           --out "${OUTPUT_PREFIX}_Y1" \
           --reml-maxit 100 \
           --thread-num 5 || echo "WARNING: GCTA for Y1 failed at N=${N}."
    
    STEP3A_END_TIME=$(date +%s)
    STEP3A_DURATION=$((STEP3A_END_TIME - STEP3A_START_TIME))
    echo "Step 3a (Y1) completed at: $(date) (Duration: ${STEP3A_DURATION}s)"

    echo "Step 3b (Y2) started at: $(date)"
    STEP3B_START_TIME=$(date +%s)
    
    gcta64 --reml --mgrm "${WORK_DIR}/mgrm.txt" \
           --reml-no-lrt \
           --reml-no-constrain \
           --pheno "${WORK_DIR}/offspring.phen" --mpheno 2 \
           --out "${OUTPUT_PREFIX}_Y2" \
           --reml-maxit 100 \
           --thread-num 5 || echo "WARNING: GCTA for Y2 failed at N=${N}."

    STEP3B_END_TIME=$(date +%s)
    STEP3B_DURATION=$((STEP3B_END_TIME - STEP3B_START_TIME))
    echo "Step 3b (Y2) completed at: $(date) (Duration: ${STEP3B_DURATION}s)"

    N_END_TIME=$(date +%s)
    N_DURATION=$((N_END_TIME - N_START_TIME))
    echo "--- Finished analysis for N=${N} at: $(date) (Total duration: ${N_DURATION}s) ---"
    echo "Timing breakdown for N=${N}:"
    echo "  - Step 1 (GRM calculation): ${STEP1_DURATION}s"
    echo "  - Step 2 (Multi-GRM file): ${STEP2_DURATION}s"
    echo "  - Step 3a (Y1 GREML): ${STEP3A_DURATION}s"
    echo "  - Step 3b (Y2 GREML): ${STEP3B_DURATION}s"
    
    # rm -rf ${WORK_DIR} # Keep cleanup commented for debugging
done

SCRIPT_END_TIME=$(date +%s)
SCRIPT_DURATION=$((SCRIPT_END_TIME - SCRIPT_START_TIME))
echo "--- All sample sizes for ${RUN_ID} complete at: $(date) ---"
echo "Total script duration: ${SCRIPT_DURATION}s"