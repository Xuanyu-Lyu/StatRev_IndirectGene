# File: run_rdr_analysis.py

import os
import sys
import glob
import multiprocessing
import pandas as pd
import numpy as np
from scipy.optimize import minimize

def _ensure_folder_exists(folder_path):
    """Helper function to create a folder if it doesn't exist."""
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

# In run_rdr_analysis.py

# ... (imports and other functions remain the same) ...

def run_rdr_on_replication(task_params):
    """
    Worker function to run RDR analysis on a single simulation replication for a single trait.
    """
    run_folder_path = task_params['run_folder_path']
    trait_to_analyze = task_params['trait']
    condition_name = task_params['condition_name']
    replication_id = int(os.path.basename(run_folder_path).split('_')[-1])
    
    base_result = {'replication': replication_id, 'trait': trait_to_analyze, 'condition_name': condition_name}

    try:
        # --- 1. Find and Load Data ---
        phen_files = glob.glob(os.path.join(run_folder_path, "*_phen_gen*.tsv"))
        if not phen_files:
            return {**base_result, 'status': 'failed', 'error': 'No phenotype files found'}
        
        final_gen_num = max([int(f.split('_gen')[1].split('.tsv')[0]) for f in phen_files])
        phen_filepath = [f for f in phen_files if f"gen{final_gen_num}.tsv" in f][0]
        
        xo_filepath = phen_filepath.replace('_phen_', '_xo_')
        if not os.path.exists(xo_filepath):
             return {**base_result, 'status': 'failed', 'error': f'XO file not found for {phen_filepath}'}

        # Load phenotype data (this has a header)
        df_phen_offspring = pd.read_csv(phen_filepath, sep='\t')
        
        # *** FIX IS HERE: Load genotype data with header=None ***
        df_gene_o = pd.read_csv(xo_filepath, sep='\t', header=None)
        
        # --- 2. Prepare Data and Matrices ---
        Y_offspring = df_phen_offspring[trait_to_analyze].values
        
        parent_gen_num = final_gen_num - 1
        parent_phen_filepath_list = [f for f in glob.glob(os.path.join(run_folder_path, "*_phen_gen*.tsv")) if f"gen{parent_gen_num}.tsv" in f]
        if not parent_phen_filepath_list:
            return {**base_result, 'status': 'failed', 'error': f'Parent generation file gen{parent_gen_num} not found'}
        parent_phen_filepath = parent_phen_filepath_list[0]
        parent_xo_filepath = parent_phen_filepath.replace('_phen_', '_xo_')

        df_phen_parents = pd.read_csv(parent_phen_filepath, sep='\t')
        
        # *** FIX IS HERE: Load parent genotype data with header=None ***
        df_gene_parents_full = pd.read_csv(parent_xo_filepath, sep='\t', header=None)
        
        parent_id_to_idx = {id_val: i for i, id_val in enumerate(df_phen_parents['ID'])}
        father_indices = [parent_id_to_idx.get(fid) for fid in df_phen_offspring['Father.ID']]
        mother_indices = [parent_id_to_idx.get(mid) for mid in df_phen_offspring['Mother.ID']]

        valid_parent_mask = np.array([(f is not None) and (m is not None) for f, m in zip(father_indices, mother_indices)])
        
        df_gene_o = df_gene_o.loc[valid_parent_mask].reset_index(drop=True)
        Y_offspring = Y_offspring[valid_parent_mask]
        
        final_father_indices = [idx for i, idx in enumerate(father_indices) if valid_parent_mask[i]]
        final_mother_indices = [idx for i, idx in enumerate(mother_indices) if valid_parent_mask[i]]

        df_gene_f = df_gene_parents_full.iloc[final_father_indices]
        df_gene_m = df_gene_parents_full.iloc[final_mother_indices]
        
        df_gene_o_std = (df_gene_o - df_gene_o.mean()) / df_gene_o.std()
        df_gene_p = df_gene_m.values + df_gene_f.values
        df_gene_p_std = (df_gene_p - df_gene_p.mean(axis=0)) / (df_gene_p.std(axis=0) / np.sqrt(2))
        
        n_snps = df_gene_o_std.shape[1]
        R_SNP_o = np.dot(df_gene_o_std, df_gene_o_std.T) / n_snps
        R_SNP_p = np.dot(df_gene_p_std, df_gene_p_std.T) / (n_snps * 2)
        R_SNP_op = (np.dot(df_gene_o_std, df_gene_p_std.T) + np.dot(df_gene_p_std, df_gene_o_std.T)) / (n_snps * 2)

        y = Y_offspring - Y_offspring.mean()

        # --- 3. Run Maximum Likelihood Estimation ---
        init_params = np.array([0, np.log(y.var()*0.5), np.log(y.var()*0.1), np.log(0.01), np.log(y.var()*0.4)])
        
        res = minimize(
            fun=neg_log_lik, x0=init_params, args=(y, R_SNP_o, R_SNP_p, R_SNP_op),
            method='L-BFGS-B',
        )
        
        # --- 4. Process and Return Results ---
        if res.success:
            p_hat = res.x
            estimates = {
                'mu_hat': p_hat[0], 'v_g_hat': np.exp(p_hat[1]), 'v_e_g_hat': np.exp(p_hat[2]),
                'c_g_e_hat': np.exp(p_hat[3]), 'sigma2_hat': np.exp(p_hat[4])
            }
            inv_hessian = res.hess_inv.todense() if hasattr(res.hess_inv, 'todense') else res.hess_inv
            se_log_params = np.sqrt(np.diag(inv_hessian))
            ses = {
                'se_mu': se_log_params[0], 'se_v_g': np.exp(p_hat[1]) * se_log_params[1],
                'se_v_e_g': np.exp(p_hat[2]) * se_log_params[2], 'se_c_g_e': np.exp(p_hat[3]) * se_log_params[3],
                'se_sigma2': np.exp(p_hat[4]) * se_log_params[4]
            }
            return {**base_result, 'status': 'success', **estimates, **ses}
        else:
            return {**base_result, 'status': 'failed', 'error': res.message}

    except Exception as e:
        return {**base_result, 'status': 'error', 'error': str(e)}


def build_Sigma(params, R_snp, R_par, R_op):
    mu, alpha_g, alpha_e_g, alpha_g_e, alpha_sig = params
    v_g, v_e_g, c_g_e, sig2 = np.exp(alpha_g), np.exp(alpha_e_g), np.exp(alpha_g_e), np.exp(alpha_sig)
    Sigma = v_g * R_snp + v_e_g * R_par + c_g_e * R_op + sig2 * np.eye(R_snp.shape[0])
    return Sigma, mu

def neg_log_lik(params, y, R_snp, R_par, R_op):
    n = len(y)
    Sigma, mu = build_Sigma(params, R_snp, R_par, R_op)
    ym = y - mu
    try:
        sign, logdet = np.linalg.slogdet(Sigma)
        if sign <= 0: return 1e12
        invSy = np.linalg.solve(Sigma, ym)
        quadform = ym.dot(invSy)
        nll = 0.5 * (n * np.log(2 * np.pi) + logdet + quadform)
        return nll if np.isfinite(nll) else 1e12
    except np.linalg.LinAlgError:
        return 1e12

def main():
    # --- 1. Configuration ---
    SOURCE_DATA_DIR = "/scratch/alpine/xuly4739/StatRev_IndirectGene/Data/ASHG_Preliminary" # Make sure this points to the right batch folder
    DESTINATION_DIR = "/projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/RDR_Results"
    CONDITIONS_TO_PROCESS = ["phenotypic_transmission", "social_transmission"]
    TRAITS_TO_ANALYZE = ["Y1", "Y2"] # *** NEW: Specify traits to analyze ***
    NUM_PROCESSES = int(os.environ.get('SLURM_CPUS_PER_TASK', 10))

    # --- 2. Generate Task List for All Replications and Traits ---
    tasks = []
    print("--- Preparing RDR Analysis Tasks ---")
    for condition in CONDITIONS_TO_PROCESS:
        condition_source_path = os.path.join(SOURCE_DATA_DIR, condition)
        if not os.path.isdir(condition_source_path):
            print(f"Warning: Source directory for condition '{condition}' not found. Skipping.")
            continue
            
        run_folders = glob.glob(os.path.join(condition_source_path, "run_*"))
        
        for run_folder in run_folders:
            # *** MODIFICATION: Create a task for each trait ***
            for trait in TRAITS_TO_ANALYZE:
                tasks.append({
                    'run_folder_path': run_folder,
                    'condition_name': condition,
                    'trait': trait
                })

    if not tasks:
        print("No tasks to run. Check source directory."); return

    print(f"Total analysis tasks to run: {len(tasks)} ({len(tasks)//len(TRAITS_TO_ANALYZE)} replications x {len(TRAITS_TO_ANALYZE)} traits)")
    print(f"Using {NUM_PROCESSES} parallel processes.")

    # --- 3. Run All Tasks in Parallel ---
    with multiprocessing.Pool(processes=NUM_PROCESSES) as pool:
        results_list = pool.map(run_rdr_on_replication, tasks)

    # --- 4. Aggregate and Save Results ---
    print("\n--- Aggregating and Saving Results ---")
    full_results_df = pd.DataFrame(results_list)
    
    _ensure_folder_exists(DESTINATION_DIR)
    
    # The output file will now contain results for both traits,
    # grouped by condition and sorted by replication and trait.
    for condition_name, group_df in full_results_df.groupby('condition_name'):
        output_filename = os.path.join(DESTINATION_DIR, f"RDR_results_{condition_name}_both_traits.txt")
        group_df.sort_values(by=['replication', 'trait'], inplace=True)
        group_df.to_csv(output_filename, sep='\t', index=False)
        print(f"Saved results for condition '{condition_name}' to: {output_filename}")

    print("\n--- All Processing Complete ---")

if __name__ == '__main__':
    main()