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
    echo ">>> WORK_DIR is now: ${WORK_DIR}"
    ls -ld "${RUN_FOLDER}/gcta_analysis"

    mkdir -p ${FINAL_RESULTS_DIR}

    # Step 1: Prepare inputs for this specific subsample, passing N as an argument
    echo "Step 1: Preparing GCTA inputs for N=${N}..."
    python prepare_combined_plink.py ${RUN_FOLDER} ${WORK_DIR} ${N}
    
    # Step 2: Create PLINK binary file from the subsampled data
    plink --file "${WORK_DIR}/combined_genos" --make-bed --out "${WORK_DIR}/combined_plink" --noweb --allow-no-sex

    # Step 3: Calculate the large combined GRM for the subsample
    echo "Step 3: Calculating combined GRM for N=${N}..."
    gcta64 --bfile "${WORK_DIR}/combined_plink" --make-grm --out "${WORK_DIR}/grm_combined" --thread-num 2

    # Step 4: Partition the GRM into the three RDR components
    echo "Step 4: Partitioning GRM for N=${N}..."
    python partition_grm.py "${WORK_DIR}/grm_combined"
    # echo ">>> Checking for partitioned GRMs:"
    # for comp in Ro_offspring Rp_parental Rop_cross; do
    # gzfile="${WORK_DIR}/grm_combined_${comp}.grm.gz"
    # if [[ -f "$gzfile" ]]; then
    #     echo "    found: $gzfile"
    # else
    #     echo "ERROR: missing $gzfile — partition_grm.py didn’t write it!"
    #     exit 1
    # fi
    # done

    # echo ">>> Decompressing text GRMs for GCTA multi‐GRM:"
    # for comp in Ro_offspring Rp_parental Rop_cross; do
    # gzfile="${WORK_DIR}/grm_combined_${comp}.grm.gz"
    # grmfile="${WORK_DIR}/grm_combined_${comp}.grm"
    # gzip -dc "$gzfile" > "$grmfile"
    # echo "    decompressed $gzfile → $grmfile"
    # done

    # echo ">>> Converting each gz GRM into binary form:"
    # for comp in Ro_offspring Rp_parental Rop_cross; do
    # prefix="${WORK_DIR}/grm_combined_${comp}"
    # echo "    converting ${prefix}.grm → binary"
    # gcta64 \
    #     --grm-gz "${prefix}.grm.gz" \
    #     --make-grm \
    #     --out "${prefix}"
    # echo "    produced: ${prefix}.grm.bin + ${prefix}.grm.N.bin + ${prefix}.grm.id"
    # done


    # # Step 5: Create the multi-GRM input file (mgrm.txt)
    # echo "Step 5: Creating multi-GRM file for N=${N}..."
    # printf "%s %s\n" \
    # "grm_O" "${WORK_DIR}/grm_combined_Ro_offspring.grm.gz" \
    # "grm_P" "${WORK_DIR}/grm_combined_Rp_parental.grm.gz" \
    # "grm_OP" "${WORK_DIR}/grm_combined_Rop_cross.grm.gz" \
    # > "${WORK_DIR}/mgrm.txt"

    # Step 5: Create the multi-GRM input file (mgrm.txt)
    echo "Step 5: Creating multi-GRM file for N=${N}..."
    printf "%s \n" \
    "${WORK_DIR}/grm_combined_Ro_offspring"  \
    "${WORK_DIR}/grm_combined_Rp_parental"  \
    "${WORK_DIR}/grm_combined_Rop_cross" \
    > "${WORK_DIR}/mgrm.txt"

    # Step 6: Run Univariate GREML Analysis with 3 GRMs for each Trait
    echo ">>> Contents of ${WORK_DIR}/mgrm.txt:"
    cat "${WORK_DIR}/mgrm.txt"
    echo "Step 6: Running RDR GREML analysis for Trait 1 (Y1) with N=${N}..."
    gcta64 --reml --mgrm "${WORK_DIR}/mgrm.txt" \
           --reml-no-lrt \
           --pheno "${WORK_DIR}/offspring.phen" --mpheno 1 \
           --out "${OUTPUT_PREFIX}_Y1" \
           --reml-maxit 100 \
           --reml-priors 0.68 0.1 0.2 0.3 \
           --thread-num 5 || echo "WARNING: GCTA for Y1 failed at N=${N}. Continuing to next step."

    # Step 7: Running RDR GREML analysis for Trait 2 (Y2) with N=${N}..."
    gcta64 --reml --mgrm "${WORK_DIR}/mgrm.txt" \
           --reml-no-lrt \
           --pheno "${WORK_DIR}/offspring.phen" --mpheno 2 \
           --out "${OUTPUT_PREFIX}_Y2" \
           --reml-maxit 100 \
           --reml-priors 0.68 0.1 0.2 0.3 \
           --thread-num 5 || echo "WARNING: GCTA for Y2 failed at N=${N}. Continuing to next sample size."

    echo "--- Finished analysis for N=${N}. Cleaning up intermediate files. ---"
    # rm -rf ${WORK_DIR}
           
done

echo "--- All sample sizes for ${RUN_ID} complete. ---"