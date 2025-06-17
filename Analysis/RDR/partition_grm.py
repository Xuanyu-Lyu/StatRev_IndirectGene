# File: partition_grm.py

import numpy as np
import pandas as pd
import sys
import gzip

def save_grm_text(grm_matrix, id_df, output_prefix):
    """Saves a GRM in GCTA's zipped text format."""
    n = grm_matrix.shape[0]
    id_df.to_csv(f"{output_prefix}.grm.id", sep='\t', header=False, index=False)
    
    grm_text_data = []
    indices = np.tril_indices(n)
    for i, j in zip(indices[0], indices[1]):
        grm_text_data.append(f"{i+1}\t{j+1}\t{n}\t{grm_matrix[i, j]}")
        
    with gzip.open(f"{output_prefix}.grm.gz", "wt") as f:
        f.write("\n".join(grm_text_data))
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
    
    save_grm_text(grm_oo, offspring_ids, f"{gcta_grm_prefix}_Ro_offspring")
    save_grm_text(grm_pp, offspring_ids, f"{gcta_grm_prefix}_Rp_parental")
    save_grm_text(grm_op_sym, offspring_ids, f"{gcta_grm_prefix}_Rop_cross")
    
    print("--- Finished partitioning GRM ---")

if __name__ == '__main__':
    main(sys.argv[1])