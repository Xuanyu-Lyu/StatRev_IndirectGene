# In run_rdr_analysis.py

# ... (imports and other functions remain the same) ...

def run_rdr_on_replication(task_params):
    """
    Worker function to run RDR analysis on a single simulation replication for a single trait.
    This encapsulates the logic from your RDR script.

    Args:
        task_params (dict): Contains the path to the run folder and the trait to analyze.

    Returns:
        dict: A dictionary containing the estimated parameters and standard errors.
    """
    # --- MODIFICATION START ---
    # Unpack ALL necessary parameters at the top
    run_folder_path = task_params['run_folder_path']
    trait_to_analyze = task_params['trait']
    condition_name = task_params['condition_name'] # <-- This was the missing line
    
    replication_id = int(os.path.basename(run_folder_path).split('_')[-1])
    
    # Create a base result dictionary that now INCLUDES the condition_name
    base_result = {
        'replication': replication_id, 
        'trait': trait_to_analyze,
        'condition_name': condition_name # <-- Add it to the base result for all return paths
    }
    # --- MODIFICATION END ---

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

        df_phen_offspring = pd.read_csv(phen_filepath, sep='\t')
        df_gene_o = pd.read_csv(xo_filepath, sep='\t')
        
        # --- 2. Prepare Data and Matrices ---
        Y_offspring = df_phen_offspring[trait_to_analyze].values 
        
        parent_gen_num = final_gen_num - 1
        parent_phen_filepath = [f for f in glob.glob(os.path.join(run_folder_path, "*_phen_gen*.tsv")) if f"gen{parent_gen_num}.tsv" in f]
        if not parent_phen_filepath:
            return {**base_result, 'status': 'failed', 'error': f'Parent generation file gen{parent_gen_num} not found'}
        parent_xo_filepath = parent_phen_filepath[0].replace('_phen_', '_xo_')

        df_phen_parents = pd.read_csv(parent_phen_filepath[0], sep='\t')
        df_gene_parents_full = pd.read_csv(parent_xo_filepath, sep='\t')
        
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
            # All return paths now correctly include the base_result dictionary
            return {**base_result, 'status': 'success', **estimates, **ses}
        else:
            return {**base_result, 'status': 'failed', 'error': res.message}

    except Exception as e:
        return {**base_result, 'status': 'error', 'error': str(e)}