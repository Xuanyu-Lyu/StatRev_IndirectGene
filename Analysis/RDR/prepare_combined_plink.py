# File: prepare_combined_plink.py

import pandas as pd
import numpy as np
import sys
import os
import glob

def main(run_folder_path, work_dir):
    """Prepares a single combined PLINK file for offspring and parental midpoints."""
    os.makedirs(work_dir, exist_ok=True)
    print(f"--- Preparing combined PLINK input for {run_folder_path} ---")

    # 1. Load data
    final_gen_num = 20
    phen_filepath = glob.glob(os.path.join(run_folder_path, f"*_phen_gen{final_gen_num}.tsv"))[0]
    xo_filepath = phen_filepath.replace('_phen_', '_xo_')
    parent_phen_filepath = glob.glob(os.path.join(run_folder_path, f"*_phen_gen{final_gen_num-1}.tsv"))[0]
    parent_xo_filepath = parent_phen_filepath.replace('_phen_', '_xo_')

    df_phen_offspring = pd.read_csv(phen_filepath, sep='\t')
    df_gene_o = pd.read_csv(xo_filepath, sep='\t', header=None)
    df_phen_parents = pd.read_csv(parent_phen_filepath, sep='\t')
    df_gene_parents_full = pd.read_csv(parent_xo_filepath, sep='\t', header=None)

    # 2. Filter for offspring with parents and get parental genotypes
    parent_id_to_idx = {id_val: i for i, id_val in enumerate(df_phen_parents['ID'])}
    father_indices = [parent_id_to_idx.get(fid) for fid in df_phen_offspring['Father.ID']]
    mother_indices = [parent_id_to_idx.get(mid) for mid in df_phen_offspring['Mother.ID']]
    valid_mask = np.array([(f is not None) and (m is not None) for f, m in zip(father_indices, mother_indices)])

    df_phen_offspring_final = df_phen_offspring[valid_mask].reset_index(drop=True)
    df_gene_o_final = df_gene_o[valid_mask].reset_index(drop=True)
    final_father_indices = [idx for i, idx in enumerate(father_indices) if valid_mask[i]]
    final_mother_indices = [idx for i, idx in enumerate(mother_indices) if valid_mask[i]]
    df_gene_f = df_gene_parents_full.iloc[final_father_indices]
    df_gene_m = df_gene_parents_full.iloc[final_mother_indices]
    df_gene_p_midpoint = (df_gene_f.values + df_gene_m.values) / 2.0

    # 3. Create combined pedigree and genotype data
    # Offspring PED info
    ped_offspring = pd.DataFrame({
        'FID': df_phen_offspring_final['Father.ID'], 'IID': df_phen_offspring_final['ID'],
        'PatID': df_phen_offspring_final['Father.ID'], 'MatID': df_phen_offspring_final['Mother.ID'],
        'Sex': df_phen_offspring_final['Sex'], 'Phenotype': df_phen_offspring_final['Y1']
    })
    
    # Parental Midpoint PED info
    ped_parents = pd.DataFrame({
        'FID': df_phen_offspring_final['Father.ID'], 'IID': [f"P_{iid}" for iid in df_phen_offspring_final['ID']],
        'PatID': 0, 'MatID': 0, 'Sex': 0, 'Phenotype': -9 # Phenotype is irrelevant here
    })

    # Combine them
    ped_combined = pd.concat([ped_offspring, ped_parents], ignore_index=True)
    geno_combined = pd.concat([df_gene_o_final, pd.DataFrame(df_gene_p_midpoint)], ignore_index=True)

    # 4. Create combined PLINK .ped and .map files
    num_snps = geno_combined.shape[1]
    map_df = pd.DataFrame({'chr': 1, 'snp_id': [f'snp_{i+1}' for i in range(num_snps)], 'morgan': 0, 'bp_pos': list(range(1, num_snps + 1))})
    map_df.to_csv(os.path.join(work_dir, "combined_genos.map"), sep='\t', index=False, header=False)
    
    genotype_alleles = []
    for _, row in geno_combined.iterrows():
        alleles = ['A', 'A'] * row.shape[0] # Faster pre-allocation
        a_indices = np.where(row.values == 1)[0]
        b_indices = np.where(row.values == 2)[0]
        for idx in a_indices: alleles[2*idx+1] = 'B'
        for idx in b_indices: alleles[2*idx], alleles[2*idx+1] = 'B', 'B'
        genotype_alleles.append(alleles)

    alleles_df = pd.DataFrame(genotype_alleles)
    full_ped_df = pd.concat([ped_combined, alleles_df], axis=1)
    full_ped_df.to_csv(os.path.join(work_dir, "combined_genos.ped"), sep=' ', index=False, header=False, na_rep='0')
    print("Saved combined .ped and .map files.")

    # 5. Save phenotype file for GCTA
    pheno_gcta_df = df_phen_offspring_final[['Father.ID', 'ID', 'Y1', 'Y2']].copy()
    pheno_gcta_df.columns = ['FID', 'IID', 'Trait1', 'Trait2']
    pheno_gcta_df.to_csv(os.path.join(work_dir, "offspring.phen"), sep='\t', index=False, header=False)
    print("Saved GCTA phenotype file.")

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])