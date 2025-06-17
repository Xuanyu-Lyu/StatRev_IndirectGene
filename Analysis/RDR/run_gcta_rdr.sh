#!/bin/bash
set -e 

if [[ -z "$1" || -z "$2" ]]; then
    echo "Usage: ./run_gcta_rdr.sh <path_to_run_folder_on_scratch> <path_to_final_results_dir_on_projects>"
    exit 1
fi

RUN_FOLDER=$1
FINAL_RESULTS_DIR=$2
RUN_ID=$(basename ${RUN_FOLDER}) 

# Define the sample sizes to loop through
TARGET_SAMPLE_SIZES=(2000 4000 8000 16000 32000)

echo "--- Starting RDR analysis for ${RUN_ID} across multiple sample sizes ---"

# Loop through each target sample size
for N in "${TARGET_SAMPLE_SIZES[@]}"; do
    echo -e "\n================== Processing for N=${N} =================="
    
    # This WORK_DIR is now correctly created inside the SCRATCH RUN_FOLDER
    WORK_DIR="${RUN_FOLDER}/gcta_analysis/N_${N}"
    # This OUTPUT_PREFIX correctly points to the permanent PROJECTS directory
    OUTPUT_PREFIX="${FINAL_RESULTS_DIR}/${RUN_ID}_N${N}"
    
    mkdir -p ${WORK_DIR}
    mkdir -p ${FINAL_RESULTS_DIR}

    # Step 1: Prepare inputs for this specific subsample, passing N as an argument
    echo "Step 1: Preparing GCTA inputs for N=${N}..."
    python prepare_combined_plink.py ${RUN_FOLDER} ${WORK_DIR} ${N}
    
    # Step 2: Create PLINK binary file from the subsampled data
    plink --file "${WORK_DIR}/combined_genos" --make-bed --out "${WORK_DIR}/combined_plink" --noweb

    # Step 3: Calculate the large combined GRM for the subsample
    echo "Step 3: Calculating combined GRM for N=${N}..."
    gcta64 --bfile "${WORK_DIR}/combined_plink" --make-grm --out "${WORK_DIR}/grm_combined" --thread-num 2

    # Step 4: Partition the GRM into the three RDR components
    echo "Step 4: Partitioning GRM for N=${N}..."
    python partition_grm.py "${WORK_DIR}/grm_combined"

    # Step 5: Create the multi-GRM input file (mgrm.txt) - WITH TYPO FIXED
    echo "Step 5: Creating multi-GRM file for N=${N}..."
    echo "grm_O ${WORK_DIR}/grm_combined_Ro_offspring" > "${WORK_DIR}/mgrm.txt"
    echo "grm_P ${WORK_DIR}/grm_combined_Rp_parental" >> "${WORK_DIR}/mgrm.txt"
    echo "grm_OP ${WORK_DIR}/grm_combined_Rop_cross" >> "${WORK_DIR}/mgrm.txt" 

    # Step 6: Run Univariate GREML Analysis with 3 GRMs for each Trait
    echo "Step 6: Running RDR GREML analysis for Trait 1 (Y1) with N=${N}..."
    gcta64 --reml --mgrm "${WORK_DIR}/mgrm.txt" \
           --pheno "${WORK_DIR}/offspring.phen" --mpheno 1 \
           --out "${OUTPUT_PREFIX}_Y1" \
           --reml-maxit 100 \
           --thread-num 2

    # Step 7: Running RDR GREML analysis for Trait 2 (Y2) with N=${N}..."
    gcta64 --reml --mgrm "${WORK_DIR}/mgrm.txt" \
           --pheno "${WORK_DIR}/offspring.phen" --mpheno 2 \
           --out "${OUTPUT_PREFIX}_Y2" \
           --reml-maxit 100 \
           --thread-num 2
           
    echo "--- Finished analysis for N=${N}. Cleaning up intermediate files. ---"
    rm -rf ${WORK_DIR}
           
done

echo "--- All sample sizes for ${RUN_ID} complete. ---"