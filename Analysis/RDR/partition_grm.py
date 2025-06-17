# File: partition_grm.py

import numpy as np
import pandas as pd
import sys
import gzip


def save_grm_text(grm_matrix, id_df, output_prefix, num_snps):
    """Saves a GRM in GCTA’s zipped text format.

    grm_matrix : numpy array, shape (N,N)
    id_df      : pandas DataFrame with your sample IDs (two columns: FID IID)
    output_prefix : path/to/output/prefix (no extension)
    num_snps   : int, the number of SNPs used to build the GRM
    """
    N = grm_matrix.shape[0]
    # write the .grm.id
    id_df.to_csv(f"{output_prefix}.grm.id", sep='\t', header=False, index=False)

    # write the zipped text GRM
    with gzip.open(f"{output_prefix}.grm.gz", "wt") as f:
        # np.tril_indices includes the diagonal (i>=j)
        rows, cols = np.tril_indices(N)
        for i, j in zip(rows, cols):
            # +1 to convert 0-based python → 1-based GCTA
            f.write(f"{i+1}\t{j+1}\t{num_snps}\t{grm_matrix[i, j]}\n")

    print(f"Saved text GRM: {output_prefix}.grm.gz and .grm.id")

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
    full_grm[tril_indices] = grm_bin
    full_grm = full_grm + full_grm.T - np.diag(np.diag(full_grm))
    
    # 3. Split the matrix into four blocks
    # Assumes the first N individuals are offspring and the next N are parents
    n_offspring = n_total // 2
    
    grm_oo = full_grm[0:n_offspring, 0:n_offspring]
    grm_op = full_grm[0:n_offspring, n_offspring:]
    grm_po = full_grm[n_offspring:, 0:n_offspring] # This is grm_op.T
    grm_pp = full_grm[n_offspring:, n_offspring:]

    # 4. Symmetrize the cross-term
    grm_op_sym = (grm_op + grm_po) / 2.0
    
    # 5. Save the three final N x N matrices in GCTA text format
    offspring_ids = grm_id.iloc[0:n_offspring]
    global_M = 300

    save_grm_text(grm_oo, offspring_ids, f"{gcta_grm_prefix}_Ro_offspring", global_M)
    save_grm_text(grm_pp, offspring_ids, f"{gcta_grm_prefix}_Rp_parental", global_M)
    save_grm_text(grm_op_sym, offspring_ids, f"{gcta_grm_prefix}_Rop_cross", global_M)
    
    print("--- Finished partitioning GRM ---")

if __name__ == '__main__':
    main(sys.argv[1])