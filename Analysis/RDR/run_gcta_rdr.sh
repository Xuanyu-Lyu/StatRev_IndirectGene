#!/bin/bash
set -e 

if [[ -z "$1" || -z "$2" ]]; then
    echo "Usage: ./run_gcta_for_samplesizes.sh <path_to_run_folder> <path_to_base_output_dir>"
    exit 1
fi

RUN_FOLDER=$1
BASE_OUTPUT_DIR=$2
RUN_ID=$(basename ${RUN_FOLDER}) 

# Define the sample sizes to loop through
TARGET_SAMPLE_SIZES=(2000 4000 8000 16000 32000)

echo "--- Starting RDR analysis for ${RUN_ID} across multiple sample sizes ---"

# Loop through each target sample size
for N in "${TARGET_SAMPLE_SIZES[@]}"; do
    echo -e "\n================== Processing for N=${N} =================="
    
    # Create a unique working directory for this specific run AND sample size
    WORK_DIR="${BASE_OUTPUT_DIR}/${RUN_ID}_work/N_${N}"
    # Create a unique output prefix for the final GCTA results
    OUTPUT_PREFIX="${BASE_OUTPUT_DIR}/${RUN_ID}_N${N}"
    
    mkdir -p ${WORK_DIR}

    # Step 1: Prepare inputs for this specific subsample
    echo "Step 1: Preparing GCTA inputs for N=${N}..."
    python prepare_gcta_subsample.py ${RUN_FOLDER} ${WORK_DIR} ${N}
    
    # Step 2: Create PLINK binary file from the subsampled data
    plink --file "${WORK_DIR}/combined_genos" --make-bed --out "${WORK_DIR}/combined_plink" --noweb

    # Step 3: Calculate the large combined GRM for the subsample
    echo "Step 3: Calculating combined GRM for N=${N}..."
    gcta64 --bfile "${WORK_DIR}/combined_plink" --make-grm --out "${WORK_DIR}/grm_combined" --thread-num 2

    # Step 4: Partition the GRM into the three RDR components
    echo "Step 4: Partitioning GRM for N=${N}..."
    python partition_grm.py "${WORK_DIR}/grm_combined"

    # Step 5: Create the multi-GRM input file (mgrm.txt)
    echo "Step 5: Creating multi-GRM file for N=${N}..."
    echo -e "grm_O\n${WORK_DIR}/grm_combined_Ro_offspring" > "${WORK_DIR}/mgrm.txt"
    echo -e "grm_P\n${WORK_DIR}/grm_combined_Rp_parental" >> "${WORK_DIR}/mgrm.txt"
    echo -e "grm_OP\n${WORK_DIR}/grm_combined_Rop_cross" >> "${WORK_DIR}/mgrm.txt"

    # Step 6: Run Univariate GREML Analysis with 3 GRMs for each Trait
    echo "Step 6: Running RDR GREML analysis for Trait 1 (Y1) with N=${N}..."
    gcta64 --reml --mgrm "${WORK_DIR}/mgrm.txt" \
           --pheno "${WORK_DIR}/offspring.phen" --mpheno 1 \
           --out "${OUTPUT_PREFIX}_Y1_results" \
           --reml-maxit 100 \
           --thread-num 2

    echo "Step 7: Running RDR GREML analysis for Trait 2 (Y2) with N=${N}..."
    gcta64 --reml --mgrm "${WORK_DIR}/mgrm.txt" \
           --pheno "${WORK_DIR}/offspring.phen" --mpheno 2 \
           --out "${OUTPUT_PREFIX}_Y2_results" \
           --reml-maxit 100 \
           --thread-num 2
           
    echo "--- Finished analysis for N=${N} ---"
done

echo "--- All sample sizes for ${RUN_ID} complete. ---"