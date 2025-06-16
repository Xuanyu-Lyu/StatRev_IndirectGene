#!/bin/bash
# This script correctly implements the 3-GRM RDR analysis for a single replication.
set -e 

# Check for required arguments
if [[ -z "$1" || -z "$2" ]]; then
    echo "Usage: ./run_gcta_rdr.sh <path_to_run_folder> <path_to_final_results_dir>"
    exit 1
fi

# --- 1. Define Paths ---
RUN_FOLDER="$1"
FINAL_RESULTS_DIR="$2"
RUN_ID=$(basename "${RUN_FOLDER}") 

# Intermediate files are stored on scratch for performance
WORK_DIR="${RUN_FOLDER}/gcta_analysis_3grm" 

mkdir -p "${WORK_DIR}"
mkdir -p "${FINAL_RESULTS_DIR}"

echo "--- Starting 3-GRM RDR analysis for ${RUN_ID} ---"
echo "Working Directory: ${WORK_DIR}"

# --- Step 2: Prepare combined PLINK file using your Python script ---
echo "Step 2: Preparing combined PLINK input..."
python prepare_combined_plink.py "${RUN_FOLDER}" "${WORK_DIR}"

# --- Step 3: Convert to PLINK binary format ---
echo "Step 3: Creating PLINK binary file..."
plink --file "${WORK_DIR}/combined_genos" --make-bed --out "${WORK_DIR}/combined_plink" --noweb

# --- Step 4: Calculate ONE large GRM from the combined data ---
echo "Step 4: Calculating the large combined GRM..."
gcta64 --bfile "${WORK_DIR}/combined_plink" --make-grm --out "${WORK_DIR}/grm_combined" --thread-num 2

# --- Step 5: Partition the large GRM into three final matrices using Python ---
echo "Step 5: Partitioning the combined GRM..."
python partition_grm.py "${WORK_DIR}/grm_combined"

# --- Step 6: Create the multi-GRM input file (mgrm.txt) listing the three GRMs ---
# This format with labels and paths is correct for GCTA.
echo "Step 6: Creating multi-GRM file..."
echo -e "Offspring_GRM\n${WORK_DIR}/grm_combined_Ro_offspring" > "${WORK_DIR}/mgrm.txt"
echo -e "Parental_GRM\n${WORK_DIR}/grm_combined_Rp_parental" >> "${WORK_DIR}/mgrm.txt"
echo -e "Offspring-Parent_Cross_GRM\n${WORK_DIR}/grm_combined_Rop_cross" >> "${WORK_DIR}/mgrm.txt"

# --- Step 7: Run Univariate GREML Analysis with 3 GRMs for each Trait ---
OUTPUT_PREFIX="${FINAL_RESULTS_DIR}/${RUN_ID}"

echo "Step 7: Running RDR GREML analysis for Trait 1 (Y1)..."
gcta64 --reml --mgrm "${WORK_DIR}/mgrm.txt" \
       --pheno "${WORK_DIR}/offspring.phen" --mpheno 1 \
       --out "${OUTPUT_PREFIX}_Y1" \
       --reml-maxit 100 \
       --thread-num 2

echo "Step 8: Running RDR GREML analysis for Trait 2 (Y2)..."
gcta64 --reml --mgrm "${WORK_DIR}/mgrm.txt" \
       --pheno "${WORK_DIR}/offspring.phen" --mpheno 2 \
       --out "${OUTPUT_PREFIX}_Y2" \
       --reml-maxit 100 \
       --thread-num 2

echo "--- GCTA analysis for ${RUN_ID} complete. Final results are in ${FINAL_RESULTS_DIR}/ ---"