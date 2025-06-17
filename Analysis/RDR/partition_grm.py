# File: partition_grm.py

import numpy as np
import pandas as pd
import sys

def save_grm_binary(grm_matrix, id_df, output_prefix, num_snps):
    """Saves a GRM in GCTA’s binary format."""
    N = grm_matrix.shape[0]
    
    # 1. Write the .grm.id file
    id_df.to_csv(f"{output_prefix}.grm.id", sep='\t', header=False, index=False)

    # 2. Write the .grm.bin file (lower triangle elements as float32)
    tril_indices = np.tril_indices(N)
    grm_values = grm_matrix[tril_indices]
    grm_values.astype(np.float32).tofile(f"{output_prefix}.grm.bin")

    # 3. Write the .grm.N.bin file (number of SNPs for each element)
    num_elements = len(grm_values)
    n_snps_array = np.full(num_elements, num_snps, dtype=np.float32)
    n_snps_array.tofile(f"{output_prefix}.grm.N.bin")

    print(f"Saved BINARY GRM: {output_prefix}.grm.bin, .grm.N.bin, and .grm.id")

def main(gcta_grm_prefix):
    """Reads a large GCTA GRM and partitions it into RDR components."""
    print(f"--- Partitioning combined GRM from prefix: {gcta_grm_prefix} ---")

    # 1. Load the GCTA-generated GRM files
    grm_id = pd.read_csv(f"{gcta_grm_prefix}.grm.id", sep='\t', header=None, names=['FID', 'IID'])
    grm_bin = np.fromfile(f"{gcta_grm_prefix}.grm.bin", dtype=np.float32)

    # 2. Reconstruct the full 2N x 2N matrix
    n_total = len(grm_id)
    tril_indices = np.tril_indices(n_total)
    full_grm = np.zeros((n_total, n_total), dtype=np.float32)
    
    if len(grm_bin) != len(tril_indices[0]):
        print(f"ERROR: Size mismatch. Expected {len(tril_indices[0])} elements in .grm.bin but found {len(grm_bin)}.")
        sys.exit(1)
        
    full_grm[tril_indices] = grm_bin
    full_grm = full_grm + full_grm.T - np.diag(np.diag(full_grm))
    
    # 3. Split the matrix into four blocks
    n_offspring = n_total // 2
    
    grm_oo = full_grm[0:n_offspring, 0:n_offspring]
    grm_op = full_grm[0:n_offspring, n_offspring:]
    grm_po = full_grm[n_offspring:, 0:n_offspring]
    grm_pp = full_grm[n_offspring:, n_offspring:]

    # 4. Symmetrize the cross-term
    grm_op_sym = (grm_op + grm_po) / 2.0
    
    # 5. Save the three final N x N matrices in GCTA binary format
    offspring_ids = grm_id.iloc[0:n_offspring]
    
    # FIX: The log indicates 600 SNPs were used.
    num_snps = 600

    save_grm_binary(grm_oo, offspring_ids, f"{gcta_grm_prefix}_Ro_offspring", num_snps)
    save_grm_binary(grm_pp, offspring_ids, f"{gcta_grm_prefix}_Rp_parental", num_snps)
    save_grm_binary(grm_op_sym, offspring_ids, f"{grm_grm_prefix}_Rop_cross", num_snps)
    
    print("--- Finished partitioning GRM into BINARY format ---")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python partition_grm.py <gcta_grm_prefix>")
        sys.exit(1)
    main(sys.argv[1])