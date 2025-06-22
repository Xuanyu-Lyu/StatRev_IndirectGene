import numpy as np
import pandas as pd
import sys
import os
import glob

def save_grm_binary(grm_matrix, id_df, output_prefix, num_snps):
    """Saves a GRM in GCTA’s binary format."""
    N = grm_matrix.shape[0]
    id_df.to_csv(f"{output_prefix}.grm.id", sep='\t', header=False, index=False)
    tril_indices = np.tril_indices(N)
    grm_values = grm_matrix[tril_indices]
    grm_values.astype(np.float32).tofile(f"{output_prefix}.grm.bin")
    n_snps_array = np.full(len(grm_values), num_snps, dtype=np.float32)
    n_snps_array.tofile(f"{output_prefix}.grm.N.bin")
    print(f"Saved BINARY GRM: {output_prefix}.grm.bin, .grm.N.bin, and .grm.id")

def main(run_folder_path, work_dir, sample_size):
    """
    Calculates RDR GRMs using the z-score method from raw data and saves them
    in GCTA binary format.
    """
    sample_size = int(sample_size)
    print(f"--- Calculating RDR GRMs for N={sample_size} using z-score method ---")

    # Step 1: Load and sample the raw data (logic from prepare_combined_plink.py)
    final_gen_num = 20
    phen_filepath = glob.glob(os.path.join(run_folder_path, f"*_phen_gen{final_gen_num}.tsv"))[0]
    xo_filepath = phen_filepath.replace('_phen_', '_xo_')
    xl_filepath = phen_filepath.replace('_phen_', '_xl_')
    parent_phen_filepath = glob.glob(os.path.join(run_folder_path, f"*_phen_gen{final_gen_num-1}.tsv"))[0]
    parent_xo_filepath = parent_phen_filepath.replace('_phen_', '_xo_')
    parent_xl_filepath = parent_phen_filepath.replace('_phen_', '_xl_')

    df_phen_offspring_full = pd.read_csv(phen_filepath, sep='\t')
    df_gene_o_full = pd.concat([pd.read_csv(xo_filepath, sep='\t', header=None), pd.read_csv(xl_filepath, sep='\t', header=None)], axis=1)
    df_phen_parents_full = pd.read_csv(parent_phen_filepath, sep='\t')
    df_gene_parents_full = pd.concat([pd.read_csv(parent_xo_filepath, sep='\t', header=None), pd.read_csv(parent_xl_filepath, sep='\t', header=None)], axis=1)
    
    parent_id_to_idx = {id_val: i for i, id_val in enumerate(df_phen_parents_full['ID'])}
    df_phen_offspring_full['father_idx'] = df_phen_offspring_full['Father.ID'].map(parent_id_to_idx)
    df_phen_offspring_full['mother_idx'] = df_phen_offspring_full['Mother.ID'].map(parent_id_to_idx)
    valid_trios_df = df_phen_offspring_full.dropna(subset=['father_idx', 'mother_idx'])
    independent_offspring_df = valid_trios_df.drop_duplicates(subset=['Father.ID'], keep='first')
    
    if len(independent_offspring_df) < sample_size:
        sample_size = len(independent_offspring_df)
    
    sampled_trios_df = independent_offspring_df.sample(n=sample_size, random_state=42)
    
    offspring_indices = sampled_trios_df.index
    father_indices = sampled_trios_df['father_idx'].astype(int).values
    mother_indices = sampled_trios_df['mother_idx'].astype(int).values

    df_phen_offspring = df_phen_offspring_full.loc[offspring_indices].reset_index(drop=True)
    
    # print the phenotypic variance of the offspring
    phen_variance = df_phen_offspring['Y1'].var()
    print(f"Phenotypic variance of the offspring Y1: {phen_variance}")
    # print the phenotypic variance of the offspring Y2
    phen_variance_y2 = df_phen_offspring['Y2'].var()
    print(f"Phenotypic variance of the offspring Y2: {phen_variance_y2}")
    
    df_gene_o = df_gene_o_full.loc[offspring_indices].reset_index(drop=True)
    df_gene_f = df_gene_parents_full.iloc[father_indices].reset_index(drop=True)
    df_gene_m = df_gene_parents_full.iloc[mother_indices].reset_index(drop=True)

    # Step 2: Calculate GRMs using the method from RDR.py
    num_snps = df_gene_o.shape[1]

    # Offspring GRM
    df_gene_o_std = (df_gene_o - df_gene_o.mean()) / df_gene_o.std()
    grm_oo = np.dot(df_gene_o_std, df_gene_o_std.T) / num_snps
    
    # Parental GRM (using the summed method from RDR.py)
    df_gene_p_sum = df_gene_f.values + df_gene_m.values
    df_gene_p_sum = pd.DataFrame(df_gene_p_sum, columns=df_gene_f.columns)
    df_gene_p_sum_std = (df_gene_p_sum - df_gene_p_sum.mean()) / (df_gene_p_sum.std()/np.sqrt(2)) # Original RDR.py std calculation may need checking
    grm_pp = np.dot(df_gene_p_sum_std, df_gene_p_sum_std.T) / (num_snps*2)

    # Cross GRM
    grm_op = (np.dot(df_gene_o_std, df_gene_p_sum_std.T) + np.dot(df_gene_p_sum_std, df_gene_o_std.T)) / (num_snps*2)

    # Step 3: Save the three matrices in GCTA binary format
    output_prefix = os.path.join(work_dir, "rdr_grm")
    offspring_ids = df_phen_offspring[['Father.ID', 'ID']].copy()
    offspring_ids.columns = ['FID', 'IID']
    
    save_grm_binary(grm_oo, offspring_ids, f"{output_prefix}_Ro_offspring", num_snps)
    save_grm_binary(grm_pp, offspring_ids, f"{output_prefix}_Rp_parental", num_snps)
    save_grm_binary(grm_op, offspring_ids, f"{output_prefix}_Rop_cross", num_snps)

    # Step 4: Also need to save the phenotype file for GCTA
    pheno_gcta_df = df_phen_offspring[['Father.ID', 'ID', 'Y1', 'Y2']].copy()
    pheno_gcta_df.to_csv(os.path.join(work_dir, "offspring.phen"), sep='\t', index=False, header=False)
    
    print("--- Finished calculating and saving custom RDR GRMs ---")

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python calculate_rdr_grms.py <path_to_run_folder> <output_work_dir> <sample_size>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])