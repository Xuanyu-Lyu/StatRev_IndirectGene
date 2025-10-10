#!/usr/bin/env python3
"""
RDR (Relatedness Disequilibrium Regression) Analysis using Python regression
Based on the method implemented in PrelimResults.ipynb

This script runs RDR analysis on a single simulation run, for both Y1 and Y2 traits.
Results are saved as JSON files for later aggregation.

Author: Xuanyu Lyu
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
import sys
import os
import glob
from datetime import datetime

def load_data_from_run_folder(run_folder_path, sample_size):
    """
    Load and prepare data from a simulation run folder, similar to prepare_grm_noGCTA.py
    """
    print(f"Loading data from {run_folder_path} with sample size {sample_size}")
    
    # Find the final generation files (Generation 20)
    final_gen_num = 20
    phen_filepath = glob.glob(os.path.join(run_folder_path, f"*_phen_gen{final_gen_num}.tsv"))[0]
    xo_filepath = phen_filepath.replace('_phen_', '_xo_')
    xl_filepath = phen_filepath.replace('_phen_', '_xl_')
    parent_phen_filepath = glob.glob(os.path.join(run_folder_path, f"*_phen_gen{final_gen_num-1}.tsv"))[0]
    parent_xo_filepath = parent_phen_filepath.replace('_phen_', '_xo_')
    parent_xl_filepath = parent_phen_filepath.replace('_phen_', '_xl_')

    # Load offspring data
    df_phen_offspring_full = pd.read_csv(phen_filepath, sep='\t')
    df_gene_o_full = pd.concat([
        pd.read_csv(xo_filepath, sep='\t', header=None), 
        pd.read_csv(xl_filepath, sep='\t', header=None)
    ], axis=1)
    
    # Load parent data
    df_phen_parents_full = pd.read_csv(parent_phen_filepath, sep='\t')
    df_gene_parents_full = pd.concat([
        pd.read_csv(parent_xo_filepath, sep='\t', header=None), 
        pd.read_csv(parent_xl_filepath, sep='\t', header=None)
    ], axis=1)
    
    # Create parent ID mapping
    parent_id_to_idx = {id_val: i for i, id_val in enumerate(df_phen_parents_full['ID'])}
    df_phen_offspring_full['father_idx'] = df_phen_offspring_full['Father.ID'].map(parent_id_to_idx)
    df_phen_offspring_full['mother_idx'] = df_phen_offspring_full['Mother.ID'].map(parent_id_to_idx)
    
    # Filter valid trios and ensure independent families
    valid_trios_df = df_phen_offspring_full.dropna(subset=['father_idx', 'mother_idx'])
    independent_offspring_df = valid_trios_df.drop_duplicates(subset=['Father.ID'], keep='first')
    
    # Sample data
    if len(independent_offspring_df) < sample_size:
        print(f"Warning: Only {len(independent_offspring_df)} independent offspring available, using all")
        sample_size = len(independent_offspring_df)
    
    sampled_trios_df = independent_offspring_df.sample(n=sample_size, random_state=42)
    
    # Extract indices
    offspring_indices = sampled_trios_df.index
    father_indices = sampled_trios_df['father_idx'].astype(int).values
    mother_indices = sampled_trios_df['mother_idx'].astype(int).values

    # Extract sampled data
    df_phen_offspring = df_phen_offspring_full.loc[offspring_indices].reset_index(drop=True)
    df_gene_o = df_gene_o_full.loc[offspring_indices].reset_index(drop=True)
    df_gene_f = df_gene_parents_full.iloc[father_indices].reset_index(drop=True)
    df_gene_m = df_gene_parents_full.iloc[mother_indices].reset_index(drop=True)

    return df_phen_offspring, df_gene_o, df_gene_f, df_gene_m

def calculate_grms(df_gene_o, df_gene_f, df_gene_m):
    """
    Calculate the three relatedness matrices following the RDR method
    """
    print("Calculating relatedness matrices...")
    
    num_snps = df_gene_o.shape[1]
    
    # Standardize offspring genotypes (z-score method)
    df_gene_o_std = (df_gene_o - df_gene_o.mean()) / df_gene_o.std()
    
    # Sum parental genotypes and standardize
    df_gene_p_sum = df_gene_f.values + df_gene_m.values
    df_gene_p_sum = pd.DataFrame(df_gene_p_sum, columns=df_gene_f.columns)
    df_gene_p_sum_std = (df_gene_p_sum - df_gene_p_sum.mean()) / (df_gene_p_sum.std()/np.sqrt(2))

    # Calculate the three relatedness matrices
    R_SNP_o = np.dot(df_gene_o_std, df_gene_o_std.T) / num_snps
    R_SNP_p = np.dot(df_gene_p_sum_std, df_gene_p_sum_std.T) / (num_snps * 2)
    R_SNP_op = (np.dot(df_gene_o_std, df_gene_p_sum_std.T) + np.dot(df_gene_p_sum_std, df_gene_o_std.T)) / (num_snps * 2)

    return R_SNP_o, R_SNP_p, R_SNP_op

def run_rdr_regression(Y_offspring, R_SNP_o, R_SNP_p, R_SNP_op, trait_name):
    """
    Run the RDR regression analysis for a single trait
    """
    print(f"Running RDR regression for {trait_name}...")
    
    # Create phenotypic covariance matrix
    Y_scale = Y_offspring - Y_offspring.mean()
    Y_scale = pd.DataFrame(Y_scale)
    COV_Y = Y_scale @ Y_scale.T

    # Extract lower triangular elements (off-diagonal)
    v_Y = COV_Y.values[np.tril_indices(COV_Y.shape[0], -1)]
    v_R_SNP_o = R_SNP_o[np.tril_indices(R_SNP_o.shape[0], -1)]
    v_R_SNP_p = R_SNP_p[np.tril_indices(R_SNP_p.shape[0], -1)]
    v_R_SNP_op = R_SNP_op[np.tril_indices(R_SNP_op.shape[0], -1)]

    # Run regression
    X = np.column_stack((v_R_SNP_o, v_R_SNP_p, v_R_SNP_op))
    X = sm.add_constant(X)
    model = sm.OLS(v_Y, X)
    results = model.fit()

    # Extract results
    regression_results = {
        'trait': trait_name,
        'phenotypic_variance': float(Y_offspring.var()),
        'sample_size': len(Y_offspring),
        'n_pairs': len(v_Y),
        'intercept_coef': float(results.params[0]),
        'intercept_se': float(results.bse[0]),
        'intercept_pvalue': float(results.pvalues[0]),
        'R_offspring_coef': float(results.params[1]),
        'R_offspring_se': float(results.bse[1]),
        'R_offspring_pvalue': float(results.pvalues[1]),
        'R_parental_coef': float(results.params[2]),
        'R_parental_se': float(results.bse[2]),
        'R_parental_pvalue': float(results.pvalues[2]),
        'R_cross_coef': float(results.params[3]),
        'R_cross_se': float(results.bse[3]),
        'R_cross_pvalue': float(results.pvalues[3]),
        'r_squared': float(results.rsquared),
        'adj_r_squared': float(results.rsquared_adj),
        'f_statistic': float(results.fvalue) if results.fvalue is not None else np.nan,
        'f_pvalue': float(results.f_pvalue) if results.f_pvalue is not None else np.nan,
        'aic': float(results.aic),
        'bic': float(results.bic),
        'log_likelihood': float(results.llf)
    }
    
    return regression_results

def save_results_to_rdreg_file(results, output_path):
    """
    Save results in a custom text format similar to GCTA .HEreg files
    """
    with open(output_path, 'w') as f:
        f.write("RDR-REG\n")
        f.write(f"Analysis timestamp: {results.get('analysis_timestamp', 'N/A')}\n")
        f.write(f"Run ID: {results.get('run_id', 'N/A')}\n")
        f.write(f"Trait: {results['trait']}\n")
        f.write(f"Sample size: {results['sample_size']}\n")
        f.write(f"Number of pairs: {results['n_pairs']}\n")
        f.write(f"Phenotypic variance: {results['phenotypic_variance']:.6f}\n")
        f.write("\n")
        
        f.write("Regression Results:\n")
        f.write("Component\tEstimate\tSE\tP-value\n")
        f.write(f"Intercept\t{results['intercept_coef']:.6f}\t{results['intercept_se']:.6f}\t{results['intercept_pvalue']:.6e}\n")
        f.write(f"R_offspring\t{results['R_offspring_coef']:.6f}\t{results['R_offspring_se']:.6f}\t{results['R_offspring_pvalue']:.6e}\n")
        f.write(f"R_parental\t{results['R_parental_coef']:.6f}\t{results['R_parental_se']:.6f}\t{results['R_parental_pvalue']:.6e}\n")
        f.write(f"R_cross\t{results['R_cross_coef']:.6f}\t{results['R_cross_se']:.6f}\t{results['R_cross_pvalue']:.6e}\n")
        f.write("\n")
        
        f.write("Model Statistics:\n")
        f.write(f"R-squared: {results['r_squared']:.6f}\n")
        f.write(f"Adjusted R-squared: {results['adj_r_squared']:.6f}\n")
        f.write(f"F-statistic: {results['f_statistic']:.6f}\n")
        f.write(f"F p-value: {results['f_pvalue']:.6e}\n")
        f.write(f"AIC: {results['aic']:.6f}\n")
        f.write(f"BIC: {results['bic']:.6f}\n")
        f.write(f"Log-likelihood: {results['log_likelihood']:.6f}\n")

def main(run_folder_path, output_dir, sample_size):
    """
    Main function to run RDR analysis on a single simulation run
    """
    sample_size = int(sample_size)
    run_id = os.path.basename(run_folder_path)
    
    print(f"=== Starting RDR regression analysis for {run_id} ===")
    print(f"Run folder: {run_folder_path}")
    print(f"Output directory: {output_dir}")
    print(f"Sample size: {sample_size}")
    print(f"Analysis started at: {datetime.now()}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Load data
        df_phen_offspring, df_gene_o, df_gene_f, df_gene_m = load_data_from_run_folder(run_folder_path, sample_size)
        
        # Calculate relatedness matrices
        R_SNP_o, R_SNP_p, R_SNP_op = calculate_grms(df_gene_o, df_gene_f, df_gene_m)
        
        # Run regression for both traits
        all_results = []
        
        for trait in ['Y1', 'Y2']:
            Y_offspring = df_phen_offspring[trait].values
            result = run_rdr_regression(Y_offspring, R_SNP_o, R_SNP_p, R_SNP_op, trait)
            
            # Add metadata
            result.update({
                'run_id': run_id,
                'analysis_timestamp': datetime.now().isoformat(),
                'actual_sample_size': sample_size
            })
            
            all_results.append(result)
            
            # Save individual result as text file (similar to .HEreg format)
            output_file = f"{run_id}_RDR_regression_N{sample_size}_{trait}.RDRreg"
            output_path = os.path.join(output_dir, output_file)
            
            save_results_to_rdreg_file(result, output_path)
            
            print(f"Saved results for {trait} to {output_path}")
        
        print(f"=== RDR regression analysis for {run_id} completed successfully ===")
        print(f"Analysis finished at: {datetime.now()}")
        
        return True
        
    except Exception as e:
        print(f"ERROR: RDR regression analysis failed for {run_id}: {str(e)}")
        return False

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python run_rdr_regression.py <path_to_run_folder> <output_dir> <sample_size>")
        sys.exit(1)
    
    success = main(sys.argv[1], sys.argv[2], sys.argv[3])
    sys.exit(0 if success else 1)