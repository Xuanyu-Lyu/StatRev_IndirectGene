# File: prepare_combined_plink.py

import pandas as pd
import numpy as np
import sys
import os
import glob

def create_plink_files(phen_df, geno_df, output_prefix):
    """Helper function to create .ped and .map files."""
    num_snps = geno_df.shape[1]
    map_df = pd.DataFrame({'chr': 1, 'snp_id': [f'snp_{i+1}' for i in range(num_snps)], 'morgan': 0, 'bp_pos': list(range(1, num_snps + 1))})
    map_df.to_csv(f"{output_prefix}.map", sep='\t', index=False, header=False)
    
    genotype_alleles = []
    for _, row in geno_df.iterrows():
        alleles = ['A', 'A'] * row.shape[0]
        a_indices = np.where(row.values == 1)[0]
        b_indices = np.where(row.values == 2)[0]
        for idx in a_indices: alleles[2*idx+1] = 'B'
        for idx in b_indices: alleles[2*idx], alleles[2*idx+1] = 'B', 'B'
        genotype_alleles.append(alleles)
    
    alleles_df = pd.DataFrame(genotype_alleles)
    ped_df = pd.DataFrame({'FID': phen_df['FID'], 'IID': phen_df['IID'], 'PatID': phen_df['PatID'], 'MatID': phen_df['MatID'], 'Sex': phen_df['Sex'], 'Phenotype': phen_df['Phenotype']})
    full_ped_df = pd.concat([ped_df.reset_index(drop=True), alleles_df.reset_index(drop=True)], axis=1)
    full_ped_df.to_csv(f"{output_prefix}.ped", sep=' ', index=False, header=False, na_rep='0')

def main(run_folder_path, work_dir, sample_size):
    """Prepares GCTA inputs for a specific subsample size."""
    os.makedirs(work_dir, exist_ok=True)
    sample_size = int(sample_size)
    print(f"--- Preparing GCTA inputs for {os.path.basename(run_folder_path)} with N={sample_size} ---")
    
    final_gen_num = 20
    phen_filepath = glob.glob(os.path.join(run_folder_path, f"*_phen_gen{final_gen_num}.tsv"))[0]
    xo_filepath = phen_filepath.replace('_phen_', '_xo_')
    parent_phen_filepath = glob.glob(os.path.join(run_folder_path, f"*_phen_gen{final_gen_num-1}.tsv"))[0]
    parent_xo_filepath = parent_phen_filepath.replace('_phen_', '_xo_')

    df_phen_offspring_full = pd.read_csv(phen_filepath, sep='\t')
    df_gene_o_full = pd.read_csv(xo_filepath, sep='\t', header=None)
    df_phen_parents_full = pd.read_csv(parent_phen_filepath, sep='\t')
    df_gene_parents_full = pd.read_csv(parent_xo_filepath, sep='\t', header=None)

    # *** NEW: Subsample Data to Target Size (N) ***
    if len(df_phen_offspring_full) < sample_size:
        print(f"    Warning: Available individuals ({len(df_phen_offspring_full)}) is smaller than target size ({sample_size}). Using all available.")
        sample_size = len(df_phen_offspring_full)
    
    # Take a random sample of offspring indices from the full dataset
    sampled_indices = np.random.choice(df_phen_offspring_full.index, size=sample_size, replace=False)
    
    # Filter all dataframes based on the sampled offspring
    df_phen_offspring = df_phen_offspring_full.loc[sampled_indices].reset_index(drop=True)
    df_gene_o = df_gene_o_full.loc[sampled_indices].reset_index(drop=True)
    
    # Prepare parent data corresponding to the sampled offspring
    parent_id_to_idx = {id_val: i for i, id_val in enumerate(df_phen_parents_full['ID'])}
    father_indices = [parent_id_to_idx.get(fid) for fid in df_phen_offspring['Father.ID']]
    mother_indices = [parent_id_to_idx.get(mid) for mid in df_phen_offspring['Mother.ID']]
    valid_mask = np.array([(f is not None) and (m is not None) for f, m in zip(father_indices, mother_indices)])
    
    if not np.all(valid_mask):
        raise RuntimeError("Could not find all parents for the subsampled offspring.")
        
    df_gene_f = df_gene_parents_full.iloc[father_indices]
    df_gene_m = df_gene_parents_full.iloc[mother_indices]
    df_gene_p_midpoint = (df_gene_f.values + df_gene_m.values) / 2.0

    # Create combined pedigree and genotype data FOR THE SUBSAMPLE
    ped_offspring = pd.DataFrame({'FID': df_phen_offspring['Father.ID'], 'IID': df_phen_offspring['ID'], 'PatID': df_phen_offspring['Father.ID'], 'MatID': df_phen_offspring['Mother.ID'], 'Sex': df_phen_offspring['Sex'], 'Phenotype': df_phen_offspring['Y1']})
    ped_parents = pd.DataFrame({'FID': df_phen_offspring['Father.ID'], 'IID': [f"P_{iid}" for iid in df_phen_offspring['ID']], 'PatID': 0, 'MatID': 0, 'Sex': 0, 'Phenotype': -9})
    ped_combined = pd.concat([ped_offspring, ped_parents], ignore_index=True)
    geno_combined = pd.concat([df_gene_o, pd.DataFrame(df_gene_p_midpoint)], ignore_index=True)

    # Create combined PLINK files
    create_plink_files(ped_combined, geno_combined, os.path.join(work_dir, "combined_genos"))
    print("Saved combined .ped and .map files for subsample.")

    # Save phenotype file for GCTA
    pheno_gcta_df = df_phen_offspring[['Father.ID', 'ID', 'Y1', 'Y2']].copy()
    pheno_gcta_df.columns = ['FID', 'IID', 'Trait1', 'Trait2']
    pheno_gcta_df.to_csv(os.path.join(work_dir, "offspring.phen"), sep='\t', index=False, header=False)
    print("Saved GCTA phenotype file for subsample.")

    print(f"--- Finished data prep for N={sample_size} ---")

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])