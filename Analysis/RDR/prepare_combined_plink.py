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
    """Prepares GCTA inputs for a specific subsample of N independent trios."""
    os.makedirs(work_dir, exist_ok=True)
    sample_size = int(sample_size)
    print(f"--- Preparing GCTA inputs for {os.path.basename(run_folder_path)} with N={sample_size} independent trios ---")
    
    # 1. Load all necessary data files
    final_gen_num = 20
    phen_filepath = glob.glob(os.path.join(run_folder_path, f"*_phen_gen{final_gen_num}.tsv"))[0]
    xo_filepath = phen_filepath.replace('_phen_', '_xo_')
    xl_filepath = phen_filepath.replace('_phen_', '_xl_')
    parent_phen_filepath = glob.glob(os.path.join(run_folder_path, f"*_phen_gen{final_gen_num-1}.tsv"))[0]
    parent_xo_filepath = parent_phen_filepath.replace('_phen_', '_xo_')
    parent_xl_filepath = parent_phen_filepath.replace('_phen_', '_xl_')

    df_phen_offspring_full = pd.read_csv(phen_filepath, sep='\t')
    df_gene_o_full = pd.concat([
        pd.read_csv(xo_filepath, sep='\t', header=None),
        pd.read_csv(xl_filepath, sep='\t', header=None)
    ], axis=1, ignore_index=True)
    df_phen_parents_full = pd.read_csv(parent_phen_filepath, sep='\t')
    df_gene_parents_full = pd.concat([
        pd.read_csv(parent_xo_filepath, sep='\t', header=None),
        pd.read_csv(parent_xl_filepath, sep='\t', header=None)
    ], axis=1, ignore_index=True)

    # 2. First, identify all offspring who form complete trios
    parent_id_to_idx = {id_val: i for i, id_val in enumerate(df_phen_parents_full['ID'])}
    df_phen_offspring_full['father_idx'] = df_phen_offspring_full['Father.ID'].map(parent_id_to_idx)
    df_phen_offspring_full['mother_idx'] = df_phen_offspring_full['Mother.ID'].map(parent_id_to_idx)
    valid_trios_df = df_phen_offspring_full.dropna(subset=['father_idx', 'mother_idx']).copy()

    # 3. From the valid trios, get one offspring per family to ensure independence
    # Using Father.ID as the family identifier, similar to your R script
    independent_offspring_df = valid_trios_df.drop_duplicates(subset=['Father.ID'], keep='first')
    print(f"    Found {len(independent_offspring_df)} total independent trios in the raw data.")

    # 4. Now, sample N trios from this set of independent families
    if len(independent_offspring_df) < sample_size:
        print(f"    Warning: Available independent trios ({len(independent_offspring_df)}) is smaller than target N ({sample_size}). Using all available.")
        sample_size = len(independent_offspring_df)
    
    # Take a random sample of N trios from the independent set
    sampled_trios_df = independent_offspring_df.sample(n=sample_size, random_state=42)
    
    # 5. Get the final data based on the sampled trios
    offspring_indices = sampled_trios_df.index
    father_indices = sampled_trios_df['father_idx'].astype(int).values
    mother_indices = sampled_trios_df['mother_idx'].astype(int).values

    df_phen_offspring = df_phen_offspring_full.loc[offspring_indices].reset_index(drop=True)
    df_gene_o = df_gene_o_full.loc[offspring_indices].reset_index(drop=True)
    print(f"    Loaded and combined offspring & parent genotypes. Total SNPs: {df_gene_o_full.shape[1]}")

    # Get the corresponding parental genotypes
    df_gene_f = df_gene_parents_full.iloc[father_indices]
    df_gene_m = df_gene_parents_full.iloc[mother_indices]
    df_gene_p_midpoint = (df_gene_f.values + df_gene_m.values) / 2.0

    # 6. Create PLINK and phenotype files for the final N x M dataset
    phen_offspring_plink = pd.DataFrame({'FID': df_phen_offspring['Father.ID'], 'IID': df_phen_offspring['ID'], 'PatID': df_phen_offspring['Father.ID'], 'MatID': df_phen_offspring['Mother.ID'], 'Sex': df_phen_offspring['Sex'], 'Phenotype': df_phen_offspring['Y1']})
    create_plink_files(phen_offspring_plink, df_gene_o, os.path.join(work_dir, "offspring_genos"))
    
    phen_parental_plink = phen_offspring_plink.copy()
    create_plink_files(phen_parental_plink, pd.DataFrame(df_gene_p_midpoint), os.path.join(work_dir, "parental_midpoint_genos"))
    
    pheno_gcta_df = df_phen_offspring[['Father.ID', 'ID', 'Y1', 'Y2']].copy()
    pheno_gcta_df.columns = ['FID', 'IID', 'Trait1', 'Trait2']
    pheno_gcta_df.to_csv(os.path.join(work_dir, "offspring.phen"), sep='\t', index=False, header=False)
    
    print(f"--- Finished data prep for N={sample_size} independent trios ---")

if __name__ == '__main__':
    # Usage: python prepare_gcta_input.py <path_to_run_folder> <output_work_dir> <sample_size>
    if len(sys.argv) != 4:
        print("Usage: python prepare_gcta_input.py <path_to_run_folder> <output_work_dir> <sample_size>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])