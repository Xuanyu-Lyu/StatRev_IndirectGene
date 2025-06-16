#!/bin/bash
set -e 

if [[ -z "$1" || -z "$2" ]]; then
    echo "Usage: ./run_gcta_rdr.sh <path_to_run_folder> <path_to_output_dir>"
    exit 1
fi

RUN_FOLDER=$1
OUTPUT_DIR=$2
RUN_ID=$(basename ${RUN_FOLDER}) 
WORK_DIR="${OUTPUT_DIR}/${RUN_ID}_work"

mkdir -p ${WORK_DIR}
echo "--- Starting RDR analysis for ${RUN_ID} (New Workflow) ---"
echo "Working Directory: ${WORK_DIR}"

# --- Step 1: Prepare combined PLINK file ---
echo "Step 1: Preparing combined PLINK input..."
python prepare_combined_plink.py ${RUN_FOLDER} ${WORK_DIR}
plink --file "${WORK_DIR}/combined_genos" --make-bed --out "${WORK_DIR}/combined_plink" --noweb

# --- Step 2: Calculate ONE large GRM from the combined data ---
echo "Step 2: Calculating the large combined GRM..."
gcta64 --bfile "${WORK_DIR}/combined_plink" --make-grm --out "${WORK_DIR}/grm_combined" --thread-num 2

# --- Step 3: Partition the large GRM into three final matrices using Python ---
echo "Step 3: Partitioning the combined GRM..."
python partition_grm.py "${WORK_DIR}/grm_combined"

# --- Step 4: Create the multi-GRM input file (mgrm.txt) listing the three new GRMs ---
echo "Step 4: Creating multi-GRM file..."
echo -e "grm_O\n${WORK_DIR}/grm_combined_Ro_offspring" > "${WORK_DIR}/mgrm.txt"
echo -e "grm_P\n${WORK_DIR}/grm_combined_Rp_parental" >> "${WORK_DIR}/mgrm.txt"
echo -e "grm_OP\n${WORK_DIR}/grm_combined_Rop_cross" >> "${WORK_DIR}/mgrm.txt"

# --- Step 5: Run Univariate GREML Analysis with 3 GRMs for each Trait ---
echo "Step 5: Running RDR GREML analysis for Trait 1 (Y1)..."
gcta64 --reml --mgrm "${WORK_DIR}/mgrm.txt" \
       --pheno "${WORK_DIR}/offspring.phen" --mpheno 1 \
       --out "${OUTPUT_DIR}/${RUN_ID}_Y1_results" \
       --reml-maxit 100 \
       --thread-num 2

echo "Step 6: Running RDR GREML analysis for Trait 2 (Y2)..."
gcta64 --reml --mgrm "${WORK_DIR}/mgrm.txt" \
       --pheno "${WORK_DIR}/offspring.phen" --mpheno 2 \
       --out "${OUTPUT_DIR}/${RUN_ID}_Y2_results" \
       --reml-maxit 100 \
       --thread-num 2

echo "--- GCTA analysis for ${RUN_ID} complete. ---"