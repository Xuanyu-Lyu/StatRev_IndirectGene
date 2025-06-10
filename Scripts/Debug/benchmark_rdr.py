# File: benchmark_rdr.py

import os
import sys
import glob
import time
import resource # Used for checking memory usage on Linux/Unix systems
import pandas as pd
import numpy as np
from scipy.optimize import minimize

def run_rdr_on_subset(full_data, sample_size_n, snp_size_m):
    """
    Runs the RDR analysis on a subsampled set of data (N individuals, M SNPs)
    and returns performance metrics.
    
    Args:
        full_data (dict): A dictionary containing the full loaded DataFrames.
        sample_size_n (int): The number of individuals (N) to subsample to.
        snp_size_m (int): The number of SNPs (M) to subsample to.
        
    Returns:
        dict: A dictionary containing the analysis results and performance metrics.
    """
    df_phen_offspring = full_data['phen_offspring']
    df_gene_o = full_data['gene_o']
    df_phen_parents = full_data['phen_parents']
    df_gene_parents_full = full_data['gene_parents']
    
    # --- 1. Subsample Individuals (Rows) to size N ---
    if len(df_phen_offspring) < sample_size_n:
        print(f"    Warning: Available individuals ({len(df_phen_offspring)}) is smaller than target N ({sample_size_n}). Using all available.")
        sample_size_n = len(df_phen_offspring)
    
    sampled_indices = np.random.choice(df_phen_offspring.index, size=sample_size_n, replace=False)
    
    df_phen_offspring_sub_n = df_phen_offspring.loc[sampled_indices].reset_index(drop=True)
    df_gene_o_sub_n = df_gene_o.loc[sampled_indices].reset_index(drop=True)
    
    Y_offspring = df_phen_offspring_sub_n["Y1"].values # Benchmark on Y1

    # --- 2. Subsample SNPs (Columns) to size M ---
    total_snps_available = df_gene_o.shape[1]
    if snp_size_m > total_snps_available:
        print(f"    Warning: Target SNP size M ({snp_size_m}) is larger than available ({total_snps_available}). Using all SNPs.")
        snp_size_m = total_snps_available

    snp_indices_to_keep = np.random.choice(total_snps_available, size=snp_size_m, replace=False)
    
    # Apply SNP subsampling to the already row-subsampled genotype data
    df_gene_o_sub_nm = df_gene_o_sub_n.iloc[:, snp_indices_to_keep]
    
    # Also subsample the parent genotypes consistently
    parent_id_to_idx = {id_val: i for i, id_val in enumerate(df_phen_parents['ID'])}
    father_indices = [parent_id_to_idx.get(fid) for fid in df_phen_offspring_sub_n['Father.ID']]
    mother_indices = [parent_id_to_idx.get(mid) for mid in df_phen_offspring_sub_n['Mother.ID']]

    valid_parent_mask = np.array([(f is not None) and (m is not None) for f, m in zip(father_indices, mother_indices)])
    if not np.all(valid_parent_mask):
        return {'sample_size_n': sample_size_n, 'snp_size_m': snp_size_m, 'status': 'error', 'error': 'Could not find all parents for the subsample'}

    df_gene_f_sub_m = df_gene_parents_full.iloc[father_indices, snp_indices_to_keep]
    df_gene_m_sub_m = df_gene_parents_full.iloc[mother_indices, snp_indices_to_keep]
    
    # --- 3. Run the RDR analysis logic on the N x M subsampled data ---
    df_gene_o_std = (df_gene_o_sub_nm - df_gene_o_sub_nm.mean()) / df_gene_o_sub_nm.std()
    df_gene_p = df_gene_m_sub_m.values + df_gene_f_sub_m.values
    df_gene_p_std = (df_gene_p - df_gene_p.mean(axis=0)) / (df_gene_p.std(axis=0) / np.sqrt(2))
    
    R_SNP_o = np.dot(df_gene_o_std, df_gene_o_std.T) / snp_size_m
    R_SNP_p = np.dot(df_gene_p_std, df_gene_p_std.T) / (snp_size_m * 2)
    R_SNP_op = (np.dot(df_gene_o_std, df_gene_p_std.T) + np.dot(df_gene_p_std, df_gene_o_std.T)) / (snp_size_m * 2)
    
    y = Y_offspring - Y_offspring.mean()
    
    init_params = np.array([0, np.log(y.var()*0.5), np.log(y.var()*0.1), np.log(0.01), np.log(y.var()*0.4)])
    res = minimize(fun=neg_log_lik, x0=init_params, args=(y, R_SNP_o, R_SNP_p, R_SNP_op), method='L-BFGS-B')
    
    success = 1 if res.success else 0
    message = res.message if not res.success else "Success"
    
    return {'sample_size_n': sample_size_n, 'snp_size_m': snp_size_m, 'status': success, 'message': message}

# --- MLE Helper Functions (same as before) ---
def build_Sigma(params, R_snp, R_par, R_op):
    mu, alpha_g, alpha_e_g, alpha_g_e, alpha_sig = params; v_g, v_e_g, c_g_e, sig2 = np.exp(alpha_g), np.exp(alpha_e_g), np.exp(alpha_g_e), np.exp(alpha_sig)
    Sigma = v_g * R_snp + v_e_g * R_par + c_g_e * R_op + sig2 * np.eye(R_snp.shape[0]); return Sigma, mu
def neg_log_lik(params, y, R_snp, R_par, R_op):
    n = len(y); Sigma, mu = build_Sigma(params, R_snp, R_par, R_op); ym = y - mu
    try:
        sign, logdet = np.linalg.slogdet(Sigma);
        if sign <= 0: return 1e12
        invSy = np.linalg.solve(Sigma, ym); quadform = ym.dot(invSy); nll = 0.5 * (n * np.log(2 * np.pi) + logdet + quadform)
        return nll if np.isfinite(nll) else 1e12
    except np.linalg.LinAlgError: return 1e12


def main():
    # --- 1. Configuration for Benchmarking ---
    REPLICATION_FOLDER_TO_BENCHMARK = "/scratch/alpine/xuly4739/StatRev_IndirectGene/Data/ASHG_Preliminary/phenotypic_transmission/run_001"
    
    # *** MODIFIED: Define ranges for both N (individuals) and M (SNPs) ***
    BENCHMARK_N_VALUES = [500, 1000, 2000, 10000] 
    BENCHMARK_M_VALUES = [100, 200, 300] 

    # --- 2. Load the full dataset ONCE ---
    print(f"--- Loading full dataset from: {REPLICATION_FOLDER_TO_BENCHMARK} ---")
    try:
        phen_files = glob.glob(os.path.join(REPLICATION_FOLDER_TO_BENCHMARK, "*_phen_gen*.tsv"))
        if not phen_files: raise FileNotFoundError("No phenotype files found.")
        final_gen_num = max([int(f.split('_gen')[1].split('.tsv')[0]) for f in phen_files])
        phen_filepath = [f for f in phen_files if f"gen{final_gen_num}.tsv" in f][0]
        xo_filepath = phen_filepath.replace('_phen_', '_xo_')
        parent_phen_filepath = [f for f in phen_files if f"gen{final_gen_num-1}.tsv" in f][0]
        parent_xo_filepath = parent_phen_filepath.replace('_phen_', '_xo_')

        full_data = {
            'phen_offspring': pd.read_csv(phen_filepath, sep='\t'),
            'gene_o': pd.read_csv(xo_filepath, sep='\t', header=None),
            'phen_parents': pd.read_csv(parent_phen_filepath, sep='\t'),
            'gene_parents': pd.read_csv(parent_xo_filepath, sep='\t', header=None)
        }
        print(f"Successfully loaded data. Offspring: {len(full_data['phen_offspring'])}, SNPs: {full_data['gene_o'].shape[1]}")
    except Exception as e:
        print(f"Failed to load data for benchmarking: {e}"); sys.exit(1)
        
    # --- 3. Run benchmarks sequentially for each (N, M) combination ---
    print("\n--- Starting Benchmark Runs ---")
    benchmark_results = []
    
    # *** MODIFIED: Nested loop for N and M ***
    for n_val in BENCHMARK_N_VALUES:
        for m_val in BENCHMARK_M_VALUES:
            print(f"\n--- Benchmarking with N = {n_val}, M = {m_val} ---")
            
            start_time = time.time()
            result = run_rdr_on_subset(full_data, sample_size_n=n_val, snp_size_m=m_val)
            end_time = time.time()
            
            peak_memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            peak_memory_mb = peak_memory_kb / 1024.0
            
            elapsed_time = end_time - start_time
            
            result['time_seconds'] = elapsed_time
            result['peak_memory_mb'] = peak_memory_mb
            
            benchmark_results.append(result)
            
            print(f"  -> Time taken: {elapsed_time:.2f} seconds")
            print(f"  -> Peak Memory: {peak_memory_mb:.2f} MB")
            print(f"  -> Optimizer Success: {bool(result['status'])}")

    # --- 4. Display Final Report ---
    print("\n\n--- Benchmark Report ---")
    report_df = pd.DataFrame(benchmark_results)
    # Reorder columns for clarity in the final report
    report_cols = ['sample_size_n', 'snp_size_m', 'time_seconds', 'peak_memory_mb', 'status', 'message']
    print(report_df[report_cols].to_string())


if __name__ == '__main__':
    main()