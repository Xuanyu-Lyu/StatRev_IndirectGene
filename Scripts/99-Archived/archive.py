import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.linalg import svd as scipy_svd # Or numpy.linalg.svd
# from numpy.linalg import cholesky # For checking positive definiteness

# Utility functions (can be outside the class or static methods)
def is_even(x):
    return x % 2 == 0

def is_odd(x):
    return x % 2 != 0

def cor2cov(correlation_matrix, var1, var2):
    # Ensure correlation_matrix is a NumPy array
    correlation_matrix = np.array(correlation_matrix)
    sd_mat = np.array([[np.sqrt(var1), 0], [0, np.sqrt(var2)]])
    return sd_mat @ correlation_matrix @ sd_mat

class AssortativeMatingSimulation:
    def __init__(self, cv_info, num_generations, pop_size, avoid_inbreeding,
                 save_each_gen, save_covs, seed,
                 cove_mat, f_mat, a_mat, d_mat, am_list, covy_mat, k2_matrix):

        self.cv_info = pd.DataFrame(cv_info) # Ensure it's a DataFrame
        self.num_generations = int(num_generations)
        self.initial_pop_size = int(pop_size) # Assuming pop_size is scalar for initialization
        
        # pop_vector can be a list/array if pop_size changes per generation,
        # or a fixed array if pop_size is constant.
        if isinstance(pop_size, (list, np.ndarray)):
            self.pop_vector = np.array(pop_size, dtype=int)
            if len(self.pop_vector) != self.num_generations:
                raise ValueError("pop_size list/array length must match num_generations")
        else: # If scalar, assume constant population size
            self.pop_vector = np.full(self.num_generations, int(pop_size), dtype=int)
            
        self.avoid_inbreeding = avoid_inbreeding
        self.save_each_gen = save_each_gen
        self.save_covs = save_covs
        
        if seed != 0:
            np.random.seed(int(seed))

        # Store matrices as numpy arrays
        self.cove_mat = np.array(cove_mat)
        self.f_mat = np.array(f_mat)
        self.a_mat = np.array(a_mat)
        self.d_mat = np.array(d_mat)
        self.am_list = [np.array(m) for m in am_list] # List of assortative mating correlation matrices (mu)
        self.covy_mat = np.array(covy_mat)
        self.k2_matrix = np.array(k2_matrix) # Expected variance of BVo and BVl components

        self.num_cvs = len(self.cv_info)

        # Initialize state variables
        self.phen_df = None  # Pandas DataFrame for PHEN
        self.xo = None       # NumPy array for XO (Observed Genotypes)
        self.xl = None       # NumPy array for XL (Latent Genotypes)
        
        self.summary_results = []
        self.history = {'MATES': [], 'PHEN': [], 'XO': [], 'XL': []} if self.save_each_gen else None
        self.covariances_log = [] if self.save_covs else None
        
        # Define column names for PHEN dataframe consistently
        self.phen_column_names = [
            'ID','Father.ID','Mother.ID','Fathers.Father.ID','Fathers.Mother.ID',
            'Mothers.Father.ID','Mothers.Mother.ID','MALE',
            'AO1','AO2','AL1','AL2','F1','F2','E1','E2', # Genetic/Env components
            'AOy1','AOy2','ALy1','ALy2','Fy1','Fy2','Ey1','Ey2', # Y-scaled components
            'Y1','Y2', # Final Phenotypes
            'BV.NT.O1','BV.NT.O2', # Non-transmitted observed breeding values
            'TPO1','TPO2','TMO1','TMO2', # Transmitted Paternal/Maternal Observed
            'NTPO1','NTPO2','NTMO1','NTMO2', # Non-Transmitted Paternal/Maternal Observed
            'BV.NT.L1','BV.NT.L2', # Non-transmitted latent breeding values
            'TPL1','TPL2','TML1','TML2', # Transmitted Paternal/Maternal Latent
            'NTPL1','NTPL2','NTML1','NTML2', # Non-Transmitted Paternal/Maternal Latent
            'Y1P','Y2P','F1P','F2P', # Paternal Phenotypes Y and F components
            'Y1M','Y2M','F1M','F2M'  # Maternal Phenotypes Y and F components
        ]
        
        self._initialize_generation_zero()

    def _is_positive_definite(self, matrix):
        try:
            np.linalg.cholesky(matrix)
            return True
        except np.linalg.LinAlgError:
            return False

    def _find_unique_closest_indices(self, dist_matrix_input, max_ord=100):
        dist_matrix = dist_matrix_input.copy() 
        n_to_match_from = dist_matrix.shape[0] 
        n_targets = dist_matrix.shape[1]     

        if n_targets == 0 or n_to_match_from == 0:
            return np.array([], dtype=int)

        rand_index_cols = np.random.permutation(n_targets)
        dist_matrix_shuffled_cols = dist_matrix[:, rand_index_cols]
        
        dist_ord_mat = np.argsort(dist_matrix_shuffled_cols, axis=0)

        new_closest_indices_for_shuffled_targets = np.full(n_targets, -1, dtype=int) 
        
        for i in range(n_targets): 
            my_min_idx_val = -1 
            k = 0 
            
            potential_match_found = False
            while k < max_ord and k < n_to_match_from:
                current_preferred_individual_idx = dist_ord_mat[k, i]
                
                # Check if this individual has already been picked by a *previous* target
                # (targets with smaller indices in the shuffled order get priority)
                # Note: `new_closest_indices_for_shuffled_targets[:i]` checks against choices made by targets 0 to i-1
                if current_preferred_individual_idx not in new_closest_indices_for_shuffled_targets[new_closest_indices_for_shuffled_targets != -1][:i]: # Only check against valid previous assignments
                    my_min_idx_val = current_preferred_individual_idx
                    potential_match_found = True
                    break 
                k += 1
            
            if potential_match_found:
                new_closest_indices_for_shuffled_targets[i] = my_min_idx_val
            else:
                new_closest_indices_for_shuffled_targets[i] = -1 # R uses 1e9 + i like placeholder

        final_ordered_indices = np.full(n_targets, -1, dtype=int)
        original_indices_of_shuffled_cols = np.argsort(rand_index_cols) 
        final_ordered_indices = new_closest_indices_for_shuffled_targets[original_indices_of_shuffled_cols]
        
        # R code uses `new.closest.M[new.closest.M > 1e9] <- NA`. Here, -1 is the NA equivalent.
        return final_ordered_indices # This array may contain -1 for unmatchable targets

    def _assort_mate(self, phendata_df_current_gen, pheno_mate_corr_target_xsim, pop_size_target_offspring):
        num_males_total = int(phendata_df_current_gen['MALE'].sum())
        num_females_total = len(phendata_df_current_gen) - num_males_total

        males_phendata = phendata_df_current_gen[phendata_df_current_gen['MALE'] == 1].copy()
        females_phendata = phendata_df_current_gen[phendata_df_current_gen['MALE'] == 0].copy()

        if num_males_total > num_females_total:
            drop_indices = np.random.choice(males_phendata.index, num_males_total - num_females_total, replace=False)
            males_phendata.drop(drop_indices, inplace=True)
        elif num_females_total > num_males_total:
            drop_indices = np.random.choice(females_phendata.index, num_females_total - num_males_total, replace=False)
            females_phendata.drop(drop_indices, inplace=True)
        
        n_potential_pairs = len(males_phendata)
        if n_potential_pairs == 0:
             return {'males.PHENDATA': pd.DataFrame(columns=males_phendata.columns), 
                     'females.PHENDATA': pd.DataFrame(columns=females_phendata.columns)}

        matcor_xsim = np.array(pheno_mate_corr_target_xsim)
        if not self._is_positive_definite(matcor_xsim):
            print("MATCOR for Xsim is not positive definite. Attempting SVD correction.")
            eigenvalues, eigenvectors = np.linalg.eigh(matcor_xsim)
            eigenvalues_pd = np.maximum(eigenvalues, 1e-12) # Ensure positive
            matcor_xsim = eigenvectors @ np.diag(eigenvalues_pd) @ eigenvectors.T

        x_sim = np.random.multivariate_normal(np.zeros(4), matcor_xsim, size=n_potential_pairs, check_valid='warn')
        # To replicate R's empirical=TRUE, Xsim would need to be transformed so its sample cov is exactly matcor_xsim.
        # For simplicity here, we rely on large n_potential_pairs for approximation.
        
        x_sim_m = x_sim[:, 0:2] 
        x_sim_f = x_sim[:, 2:4] 

        # Process Males
        if not males_phendata.empty:
            m_y_values = males_phendata[['Y1', 'Y2']].values.astype(float)
            if m_y_values.shape[0] > 1 and np.all(np.std(m_y_values, axis=0, ddof=1) > 1e-9): # Check for zero std
                 m_y_scaled = (m_y_values - np.mean(m_y_values, axis=0)) / np.std(m_y_values, axis=0, ddof=1)
            else: 
                 m_y_scaled = m_y_values - np.mean(m_y_values, axis=0)
            dm_dist = cdist(m_y_scaled, x_sim_m) 
            chosen_male_indices_raw = self._find_unique_closest_indices(dm_dist)
            valid_chosen_male_indices = chosen_male_indices_raw[chosen_male_indices_raw != -1]
            # These indices are into males_phendata (which is already a slice of the original phendata)
            xm_ord2 = males_phendata.iloc[valid_chosen_male_indices]
        else:
            xm_ord2 = pd.DataFrame(columns=males_phendata.columns)

        # Process Females
        if not females_phendata.empty:
            f_y_values = females_phendata[['Y1', 'Y2']].values.astype(float)
            if f_y_values.shape[0] > 1 and np.all(np.std(f_y_values, axis=0, ddof=1) > 1e-9):
                 f_y_scaled = (f_y_values - np.mean(f_y_values, axis=0)) / np.std(f_y_values, axis=0, ddof=1)
            else:
                 f_y_scaled = f_y_values - np.mean(f_y_values, axis=0)
            df_dist = cdist(f_y_scaled, x_sim_f) 
            chosen_female_indices_raw = self._find_unique_closest_indices(df_dist)
            valid_chosen_female_indices = chosen_female_indices_raw[chosen_female_indices_raw != -1]
            xf_ord2 = females_phendata.iloc[valid_chosen_female_indices]
        else:
            xf_ord2 = pd.DataFrame(columns=females_phendata.columns)

        min_paired_len = min(len(xm_ord2), len(xf_ord2))
        males_paired = xm_ord2.iloc[:min_paired_len].reset_index(drop=True)
        females_paired = xf_ord2.iloc[:min_paired_len].reset_index(drop=True)

        if males_paired.empty:
             return {'males.PHENDATA': males_paired, 'females.PHENDATA': females_paired}

        if self.avoid_inbreeding:
            # Fill NaNs in ID columns with a sentinel value unlikely to match actual IDs
            id_cols_to_check = ['Father.ID', 'Mother.ID', 'Fathers.Father.ID', 'Fathers.Mother.ID', 'Mothers.Father.ID', 'Mothers.Mother.ID']
            sentinel_na = -999999 
            
            m_ids = males_paired[id_cols_to_check].fillna(sentinel_na).astype(int).values
            f_ids = females_paired[id_cols_to_check].fillna(sentinel_na).astype(int).values

            # no.sib.inbreeding <- males.PHENDATA[,"Father.ID"]!=females.PHENDATA[,"Father.ID"]
            # R's Father.ID is at index 1 (0-indexed) if using standard parent/grandparent order from phen_column_names
            no_sib_inbreeding = m_ids[:, 0] != f_ids[:, 0] # Comparing Father.IDs

            # no.cousin.inbreeding:
            # males.PHENDATA[,"Fathers.Father.ID"]!=females.PHENDATA[,"Fathers.Father.ID"] &
            # males.PHENDATA[,"Fathers.Father.ID"]!=females.PHENDATA[,"Mothers.Father.ID"] &
            # ... and so on for all 8 combinations of grandparental pairings for cousins.
            # Example: m_FF (idx 2), m_FM (idx 3), m_MF (idx 4), m_MM (idx 5)
            #          f_FF (idx 2), f_FM (idx 3), f_MF (idx 4), f_MM (idx 5)
            
            # This needs all combinations. E.g. male's Paternal Grandfather vs female's Paternal Grandfather
            no_cousin_mask = (m_ids[:, 2] != f_ids[:, 2]) & \
                             (m_ids[:, 2] != f_ids[:, 4]) & \
                             (m_ids[:, 4] != f_ids[:, 2]) & \
                             (m_ids[:, 4] != f_ids[:, 4]) & \
                             (m_ids[:, 3] != f_ids[:, 3]) & \
                             (m_ids[:, 3] != f_ids[:, 5]) & \
                             (m_ids[:, 5] != f_ids[:, 3]) & \
                             (m_ids[:, 5] != f_ids[:, 5])
            
            # Also ensure that an individual's ID isn't NA itself if that's possible.
            # The R code also checks that Father.ID is not NA for no_sib_inbreeding.
            # This is implicitly handled if NAs were replaced by sentinel_na, assuming sentinel_na != actual IDs.
            # A more robust check would be `(males_paired["Father.ID"].notna() & females_paired["Father.ID"].notna() & (m_ids[:, 0] != f_ids[:, 0]))`

            no_inbreeding_mask = no_sib_inbreeding & no_cousin_mask
            
            males_paired = males_paired[no_inbreeding_mask].reset_index(drop=True)
            females_paired = females_paired[no_inbreeding_mask].reset_index(drop=True)

        if males_paired.empty:
             return {'males.PHENDATA': males_paired, 'females.PHENDATA': females_paired}

        males_paired['Spouse.ID'] = females_paired['ID'].values
        females_paired['Spouse.ID'] = males_paired['ID'].values

        num_actual_mating_pairs = len(males_paired)
        if num_actual_mating_pairs == 0:
            num_offspring_per_pair = np.array([])
        else:
            lambda_offspring = pop_size_target_offspring / num_actual_mating_pairs
            num_offspring_per_pair = np.random.poisson(lambda_offspring, num_actual_mating_pairs)

            current_total_offspring = np.sum(num_offspring_per_pair)
            diff_offspring = pop_size_target_offspring - current_total_offspring
            
            if diff_offspring > 0: 
                add_indices = np.random.choice(num_actual_mating_pairs, diff_offspring, replace=True)
                for idx in add_indices:
                    num_offspring_per_pair[idx] += 1
            elif diff_offspring < 0: 
                eligible_indices_to_subtract = np.where(num_offspring_per_pair > 0)[0]
                if len(eligible_indices_to_subtract) > 0:
                    subtract_indices = np.random.choice(eligible_indices_to_subtract, abs(diff_offspring), replace=True)
                    for idx in subtract_indices:
                        if num_offspring_per_pair[idx] > 0:
                            num_offspring_per_pair[idx] -= 1
        
        males_paired['num.offspring'] = num_offspring_per_pair
        females_paired['num.offspring'] = num_offspring_per_pair
        
        return {'males.PHENDATA': males_paired, 'females.PHENDATA': females_paired}

    def _reproduce(self, mates_dict, xo_parent_gen_full, xl_parent_gen_full, phendata_parent_gen_df_full):
        males_phen_producing = mates_dict['males.PHENDATA'][mates_dict['males.PHENDATA']['num.offspring'] > 0]
        females_phen_producing = mates_dict['females.PHENDATA'][mates_dict['females.PHENDATA']['num.offspring'] > 0]

        if males_phen_producing.empty or len(males_phen_producing) != len(females_phen_producing):
            return {'PHEN': pd.DataFrame(columns=self.phen_column_names), 'XO': np.array([]), 'XL': np.array([])}

        # Create a mapping from ID to original row index in phendata_parent_gen_df_full
        id_to_original_idx_map = {id_val: i for i, id_val in enumerate(phendata_parent_gen_df_full['ID'])}
        
        male_original_indices = [id_to_original_idx_map[id_val] for id_val in males_phen_producing['ID'] if id_val in id_to_original_idx_map]
        female_original_indices = [id_to_original_idx_map[id_val] for id_val in females_phen_producing['ID'] if id_val in id_to_original_idx_map]

        # Get parental genotypes for those producing offspring
        males_gentp_obs_producers = xo_parent_gen_full[male_original_indices, :]
        males_gentp_lat_producers = xl_parent_gen_full[male_original_indices, :]
        females_gentp_obs_producers = xo_parent_gen_full[female_original_indices, :]
        females_gentp_lat_producers = xl_parent_gen_full[female_original_indices, :]

        num_offspring_counts = males_phen_producing['num.offspring'].astype(int).values
        total_offspring_generated = np.sum(num_offspring_counts)

        if total_offspring_generated == 0:
            return {'PHEN': pd.DataFrame(columns=self.phen_column_names), 'XO': np.array([]), 'XL': np.array([])}

        # Repeat parent data for each offspring
        male_rep_indices_in_producers = np.repeat(np.arange(len(males_phen_producing)), num_offspring_counts)
        female_rep_indices_in_producers = np.repeat(np.arange(len(females_phen_producing)), num_offspring_counts) # Should be same as male if pairs aligned

        fathers_phen_for_offspring = males_phen_producing.iloc[male_rep_indices_in_producers].reset_index(drop=True)
        mothers_phen_for_offspring = females_phen_producing.iloc[female_rep_indices_in_producers].reset_index(drop=True)
        
        males_gentp_obs_rep = males_gentp_obs_producers[male_rep_indices_in_producers, :]
        males_gentp_lat_rep = males_gentp_lat_producers[male_rep_indices_in_producers, :]
        females_gentp_obs_rep = females_gentp_obs_producers[female_rep_indices_in_producers, :]
        females_gentp_lat_rep = females_gentp_lat_producers[female_rep_indices_in_producers, :]

        # Haplotype Generation (Observed - XO)
        male_adder_obs = np.random.randint(0, 2, size=males_gentp_obs_rep.shape)
        xm1_obs = (males_gentp_obs_rep == 1).astype(int)
        xm2_obs = (males_gentp_obs_rep == 2).astype(int)
        males_haps_obs = (male_adder_obs * xm1_obs) + xm2_obs

        female_adder_obs = np.random.randint(0, 2, size=females_gentp_obs_rep.shape)
        xf1_obs = (females_gentp_obs_rep == 1).astype(int)
        xf2_obs = (females_gentp_obs_rep == 2).astype(int)
        females_haps_obs = (female_adder_obs * xf1_obs) + xf2_obs
        xo_new = males_haps_obs + females_haps_obs

        # Haplotype Generation (Latent - XL)
        male_adder_lat = np.random.randint(0, 2, size=males_gentp_lat_rep.shape)
        xl1_lat = (males_gentp_lat_rep == 1).astype(int)
        xl2_lat = (males_gentp_lat_rep == 2).astype(int)
        males_haps_lat = (male_adder_lat * xl1_lat) + xl2_lat

        female_adder_lat = np.random.randint(0, 2, size=females_gentp_lat_rep.shape)
        fl1_lat = (females_gentp_lat_rep == 1).astype(int)
        fl2_lat = (females_gentp_lat_rep == 2).astype(int)
        females_haps_lat = (female_adder_lat * fl1_lat) + fl2_lat
        xl_new = males_haps_lat + females_haps_lat

        # Non-Transmitted Haplotypes
        male_nt_adder_obs = 1 - male_adder_obs
        males_nt_haps_obs = (male_nt_adder_obs * xm1_obs) + xm2_obs
        female_nt_adder_obs = 1 - female_adder_obs
        females_nt_haps_obs = (female_nt_adder_obs * xf1_obs) + xf2_obs

        male_nt_adder_lat = 1 - male_adder_lat
        males_nt_haps_lat = (male_nt_adder_lat * xl1_lat) + xl2_lat
        female_nt_adder_lat = 1 - female_adder_lat
        females_nt_haps_lat = (female_nt_adder_lat * fl1_lat) + fl2_lat

        # Phenotype Components
        alphas = self.cv_info[['alpha1', 'alpha2']].values

        ao_new = xo_new @ alphas
        al_new = xl_new @ alphas
        
        tpo_new = males_haps_obs @ alphas
        tmo_new = females_haps_obs @ alphas
        ntpo_new = males_nt_haps_obs @ alphas
        ntmo_new = females_nt_haps_obs @ alphas
        bv_nt_o_new = ntpo_new + ntmo_new

        tpl_new = males_haps_lat @ alphas
        tml_new = females_haps_lat @ alphas
        ntpl_new = males_nt_haps_lat @ alphas
        ntml_new = females_nt_haps_lat @ alphas
        bv_nt_l_new = ntpl_new + ntmo_new
        
        f_new_father_contrib = fathers_phen_for_offspring[['Y1', 'Y2']].values @ self.f_mat.T
        f_new_mother_contrib = mothers_phen_for_offspring[['Y1', 'Y2']].values @ self.f_mat.T
        f_new = f_new_father_contrib + f_new_mother_contrib # This is F1, F2 for offspring

        e_new = np.random.multivariate_normal(np.zeros(2), self.cove_mat, size=total_offspring_generated)

        aoy_new = ao_new @ self.d_mat.T
        aly_new = al_new @ self.a_mat.T
        fy_new = f_new  # F component is already Y-scaled by f_mat
        ey_new = e_new  # E component is already Y-scaled (cove_mat is for Y-scaled E)
        y_new = aoy_new + aly_new + fy_new + ey_new

        # Assemble New PHEN DataFrame
        new_phen_df = pd.DataFrame(index=np.arange(total_offspring_generated), columns=self.phen_column_names)

        new_phen_df['ID'] = np.random.randint(1_000_000, 10_000_000, size=total_offspring_generated) # TODO: Ensure unique IDs if running many gens / large pops
        new_phen_df['Father.ID'] = fathers_phen_for_offspring['ID'].values
        new_phen_df['Mother.ID'] = mothers_phen_for_offspring['ID'].values
        new_phen_df['Fathers.Father.ID'] = fathers_phen_for_offspring['Father.ID'].values
        new_phen_df['Fathers.Mother.ID'] = fathers_phen_for_offspring['Mother.ID'].values
        new_phen_df['Mothers.Father.ID'] = mothers_phen_for_offspring['Father.ID'].values
        new_phen_df['Mothers.Mother.ID'] = mothers_phen_for_offspring['Mother.ID'].values

        sex_vec_offspring = np.zeros(total_offspring_generated, dtype=int)
        half_pop = total_offspring_generated // 2
        sex_vec_offspring[half_pop:] = 1 # Assign 1 to the second half for males
        if total_offspring_generated % 2 != 0: # If odd, assign last one randomly or fixed
            sex_vec_offspring[-1] = np.random.randint(0,2)
        new_phen_df['MALE'] = np.random.permutation(sex_vec_offspring)

        # Assign calculated values
        new_phen_df['AO1'], new_phen_df['AO2'] = ao_new[:,0], ao_new[:,1]
        new_phen_df['AL1'], new_phen_df['AL2'] = al_new[:,0], al_new[:,1]
        new_phen_df['F1'], new_phen_df['F2'] = f_new[:,0], f_new[:,1] # These are direct F contributions
        new_phen_df['E1'], new_phen_df['E2'] = e_new[:,0], e_new[:,1] # These are direct E contributions

        new_phen_df['AOy1'], new_phen_df['AOy2'] = aoy_new[:,0], aoy_new[:,1]
        new_phen_df['ALy1'], new_phen_df['ALy2'] = aly_new[:,0], aly_new[:,1]
        new_phen_df['Fy1'], new_phen_df['Fy2'] = fy_new[:,0], fy_new[:,1]
        new_phen_df['Ey1'], new_phen_df['Ey2'] = ey_new[:,0], ey_new[:,1]
        new_phen_df['Y1'], new_phen_df['Y2'] = y_new[:,0], y_new[:,1]

        new_phen_df['BV.NT.O1'], new_phen_df['BV.NT.O2'] = bv_nt_o_new[:,0], bv_nt_o_new[:,1]
        new_phen_df['TPO1'], new_phen_df['TPO2'] = tpo_new[:,0], tpo_new[:,1]
        new_phen_df['TMO1'], new_phen_df['TMO2'] = tmo_new[:,0], tmo_new[:,1]
        new_phen_df['NTPO1'], new_phen_df['NTPO2'] = ntpo_new[:,0], ntpo_new[:,1]
        new_phen_df['NTMO1'], new_phen_df['NTMO2'] = ntmo_new[:,0], ntmo_new[:,1]
        
        new_phen_df['BV.NT.L1'], new_phen_df['BV.NT.L2'] = bv_nt_l_new[:,0], bv_nt_l_new[:,1]
        new_phen_df['TPL1'], new_phen_df['TPL2'] = tpl_new[:,0], tpl_new[:,1]
        new_phen_df['TML1'], new_phen_df['TML2'] = tml_new[:,0], tml_new[:,1]
        new_phen_df['NTPL1'], new_phen_df['NTPL2'] = ntpl_new[:,0], ntpl_new[:,1]
        new_phen_df['NTML1'], new_phen_df['NTML2'] = ntml_new[:,0], ntml_new[:,1]

        new_phen_df['Y1P'] = fathers_phen_for_offspring['Y1'].values
        new_phen_df['Y2P'] = fathers_phen_for_offspring['Y2'].values
        new_phen_df['F1P'] = fathers_phen_for_offspring['F1'].values # Father's F component
        new_phen_df['F2P'] = fathers_phen_for_offspring['F2'].values # Father's F component
        new_phen_df['Y1M'] = mothers_phen_for_offspring['Y1'].values
        new_phen_df['Y2M'] = mothers_phen_for_offspring['Y2'].values
        new_phen_df['F1M'] = mothers_phen_for_offspring['F1'].values # Mother's F component
        new_phen_df['F2M'] = mothers_phen_for_offspring['F2'].values # Mother's F component
            
        return {'PHEN': new_phen_df, 'XO': xo_new, 'XL': xl_new}

    def _initialize_generation_zero(self):
        pop_size_gen0 = self.pop_vector[0] if len(self.pop_vector) > 0 else self.initial_pop_size
        
        maf_probs = self.cv_info['maf'].values.flatten()
        self.xo = np.random.binomial(2, maf_probs, size=(pop_size_gen0, self.num_cvs))
        self.xl = np.random.binomial(2, maf_probs, size=(pop_size_gen0, self.num_cvs))

        alphas = self.cv_info[['alpha1', 'alpha2']].values
        ao_gen0 = self.xo @ alphas
        al_gen0 = self.xl @ alphas

        f_parent1_contrib_y_scaled = np.random.multivariate_normal(np.zeros(2), self.covy_mat, size=pop_size_gen0)
        f_gen0_p1 = f_parent1_contrib_y_scaled @ self.f_mat.T 
        f_parent2_contrib_y_scaled = np.random.multivariate_normal(np.zeros(2), self.covy_mat, size=pop_size_gen0)
        f_gen0_p2 = f_parent2_contrib_y_scaled @ self.f_mat.T
        f_gen0_contrib = f_gen0_p1 + f_gen0_p2 # This is F1, F2 contributions

        e_gen0_contrib = np.random.multivariate_normal(np.zeros(2), self.cove_mat, size=pop_size_gen0)

        aoy_gen0 = ao_gen0 @ self.d_mat.T
        aly_gen0 = al_gen0 @ self.a_mat.T
        fy_gen0 = f_gen0_contrib # Already Y-scaled
        ey_gen0 = e_gen0_contrib # Already Y-scaled
        y_gen0 = aoy_gen0 + aly_gen0 + fy_gen0 + ey_gen0
        
        self.phen_df = pd.DataFrame(index=np.arange(pop_size_gen0), columns=self.phen_column_names)

        # Fictional IDs for Gen0 - R uses sample(10M:99M, size=N*7)
        # For simplicity, create unique IDs for individuals and make others distinct placeholders or NaNs
        self.phen_df['ID'] = np.arange(1, pop_size_gen0 + 1) + 10_000_000 
        for col in ['Father.ID', 'Mother.ID', 'Fathers.Father.ID', 'Fathers.Mother.ID', 
                    'Mothers.Father.ID', 'Mothers.Mother.ID']:
            self.phen_df[col] = -(np.arange(1, pop_size_gen0 + 1)) # Placeholder negative IDs or np.nan

        sex_vec_gen0 = np.zeros(pop_size_gen0, dtype=int)
        half_pop_g0 = pop_size_gen0 // 2
        sex_vec_gen0[half_pop_g0:] = 1
        if pop_size_gen0 % 2 != 0: sex_vec_gen0[-1] = np.random.randint(0,2)
        self.phen_df['MALE'] = np.random.permutation(sex_vec_gen0)

        self.phen_df['AO1'], self.phen_df['AO2'] = ao_gen0[:,0], ao_gen0[:,1]
        self.phen_df['AL1'], self.phen_df['AL2'] = al_gen0[:,0], al_gen0[:,1]
        self.phen_df['F1'], self.phen_df['F2'] = f_gen0_contrib[:,0], f_gen0_contrib[:,1]
        self.phen_df['E1'], self.phen_df['E2'] = e_gen0_contrib[:,0], e_gen0_contrib[:,1]
        self.phen_df['AOy1'], self.phen_df['AOy2'] = aoy_gen0[:,0], aoy_gen0[:,1]
        self.phen_df['ALy1'], self.phen_df['ALy2'] = aly_gen0[:,0], aly_gen0[:,1]
        self.phen_df['Fy1'], self.phen_df['Fy2'] = fy_gen0[:,0], fy_gen0[:,1]
        self.phen_df['Ey1'], self.phen_df['Ey2'] = ey_gen0[:,0], ey_gen0[:,1]
        self.phen_df['Y1'], self.phen_df['Y2'] = y_gen0[:,0], y_gen0[:,1]
        
        # Columns not present in Gen0 are filled with NaN
        for col in ['BV.NT.O1','BV.NT.O2','TPO1','TPO2','TMO1','TMO2','NTPO1','NTPO2','NTMO1','NTMO2',
                    'BV.NT.L1','BV.NT.L2','TPL1','TPL2','TML1','TML2','NTPL1','NTPL2','NTML1','NTML2',
                    'Y1P','Y2P','F1P','F2P','Y1M','Y2M','F1M','F2M']:
            self.phen_df[col] = np.nan

        if self.save_each_gen:
            self.history['MATES'].append(None) 
            self.history['PHEN'].append(self.phen_df.copy())
            self.history['XO'].append(self.xo.copy())
            self.history['XL'].append(self.xl.copy())
        
        if self.save_covs:
            self.covariances_log.append(None)

        summary_gen0 = {'GEN': 0, 'NUM.CVs': self.num_cvs, 
                        'MATE.COR': self.am_list[0].tolist() if len(self.am_list) > 0 else None,
                        'POPSIZE': pop_size_gen0}
        
        for comp in ['AOy', 'ALy', 'F', 'E', 'Y']: # Note: R uses F1,F2 and E1,E2 for VF, VE which are Fy, Ey
            var_comp_name = 'V' + comp.replace('y','').upper() # VAO, VAL, VF, VE, VP
            if comp == 'F' or comp == 'E': # R script uses 'F1'/'F2' for VF, which are Fy1/Fy2
                 cols_to_cov = [comp+'y1', comp+'y2']
            else:
                 cols_to_cov = [comp+'1', comp+'2'] if comp != 'Y' else ['Y1', 'Y2']


            if all(c in self.phen_df.columns for c in cols_to_cov):
                 df_subset = self.phen_df[cols_to_cov].dropna()
                 summary_gen0[var_comp_name] = np.cov(df_subset, rowvar=False).tolist() if len(df_subset) >=2 else np.array([[np.nan,np.nan],[np.nan,np.nan]]).tolist()
            else:
                 summary_gen0[var_comp_name] = np.array([[np.nan,np.nan],[np.nan,np.nan]]).tolist()
        
        # h2 calculations
        vp_diag = np.diag(np.array(summary_gen0.get('VP', [[np.nan,np.nan]])))
        vao_diag = np.diag(np.array(summary_gen0.get('VAO', [[np.nan,np.nan]])))
        val_diag = np.diag(np.array(summary_gen0.get('VAL', [[np.nan,np.nan]])))

        with np.errstate(divide='ignore', invalid='ignore'): # Suppress division by zero warnings
            summary_gen0['h2'] = ((vao_diag + val_diag) / vp_diag).tolist()
            summary_gen0['h2.obs'] = (vao_diag / vp_diag).tolist()
            summary_gen0['h2.lat'] = (val_diag / vp_diag).tolist()

        for key in ['covY', 'covG', 'covH', 'covI', 'w', 'q', 'covF', 'covE', 'hapsO.covs', 'hapsL.covs', 'omega', 'gamma', 'thetaNT', 'thetaT']:
            summary_gen0[key] = np.nan # These are not calculated for Gen0 in R
        self.summary_results.append(summary_gen0)

    def run_simulation(self):
        print(f"Starting simulation for {self.num_generations} generations.")
        if self.phen_df is None:
            print("Error: Generation 0 not initialized.")
            return None

        for cur_gen_py_idx in range(self.num_generations): 
            r_cur_gen_num = cur_gen_py_idx + 1 
            print(f"--- Beginning Generation {r_cur_gen_num} ---")

            pop_size_target_for_offspring = self.pop_vector[cur_gen_py_idx]
            
            # AM correlation matrix for parents of this generation's offspring
            # self.phen_df are the parents, cur_gen_py_idx is their generation index (0 to N-1)
            mate_cor_mu_for_current_parents = self.am_list[cur_gen_py_idx]

            pheno_mate_target_xsim = np.eye(4)
            y_parent_corr_df = self.phen_df[['Y1', 'Y2']].corr()
            y_parent_corr_val = y_parent_corr_df.iloc[0,1] if not y_parent_corr_df.empty else 0.0
            
            pheno_mate_target_xsim[0,1] = pheno_mate_target_xsim[1,0] = y_parent_corr_val
            pheno_mate_target_xsim[2,3] = pheno_mate_target_xsim[3,2] = y_parent_corr_val 
            pheno_mate_target_xsim[0:2, 2:4] = mate_cor_mu_for_current_parents 
            pheno_mate_target_xsim[2:4, 0:2] = mate_cor_mu_for_current_parents.T
            
            print(f"Generation {r_cur_gen_num}: Assortative mating...")
            mates = self._assort_mate(self.phen_df, pheno_mate_target_xsim, pop_size_target_for_offspring)
            
            num_m_mated = len(mates['males.PHENDATA'])
            sum_offspring = mates['males.PHENDATA']['num.offspring'].sum() if num_m_mated > 0 else 0

            if num_m_mated == 0 or sum_offspring == 0:
                print(f"Generation {r_cur_gen_num}: No reproducing pairs or no offspring produced. Simulation halting.")
                break 
            print(f"Generation {r_cur_gen_num}: Assortative mating done. {num_m_mated} pairs formed, targeting {sum_offspring} offspring.")
            
            print(f"Generation {r_cur_gen_num}: Reproduction...")
            offspring_data = self._reproduce(mates, self.xo, self.xl, self.phen_df)
            
            if offspring_data['PHEN'].empty:
                 print(f"Generation {r_cur_gen_num}: No offspring resulted from reproduction. Simulation halting.")
                 break
            print(f"Generation {r_cur_gen_num}: Reproduction done. {len(offspring_data['PHEN'])} offspring created.")

            self.xo = offspring_data['XO']
            self.xl = offspring_data['XL']
            self.phen_df = offspring_data['PHEN'] # self.phen_df is now the new offspring generation

            print(f"Generation {r_cur_gen_num}: Calculating summary statistics for new generation...")
            
            cols_for_full_cov_calc = [ # Columns used in R for `covs`
                'TPO1','TMO1','NTPO1','NTMO1','TPL1','TML1','NTPL1','NTML1',
                'TPO2','TMO2','NTPO2','NTMO2','TPL2','TML2','NTPL2','NTML2',
                'AO1','AO2','AL1','AL2','F1','F2','E1','E2', 
                'BV.NT.O1','BV.NT.O2','Y1','Y2',
                'Y1P','Y2P','Y1M','Y2M','F1P','F2P','F1M','F2M']
            
            valid_cols_for_cov_df = [col for col in cols_for_full_cov_calc if col in self.phen_df.columns and pd.api.types.is_numeric_dtype(self.phen_df[col])]
            full_cov_df = self.phen_df[valid_cols_for_cov_df].cov() # This is the `covs` matrix from R

            # --- Detailed Covariance Calculations (covG, covH, covI, w, q, omega, gamma, theta) ---
            # These are direct translations of the R script's logic for extracting these components.
            # They require careful indexing into `full_cov_df` (or sub-covariance matrices).
            
            # hapsO.covs & covG
            hapsO_cols = ['TPO1','TMO1','NTPO1','NTMO1','TPO2','TMO2','NTPO2','NTMO2']
            hapsO_covs_df = self.phen_df[hapsO_cols].cov() if all(c in self.phen_df.columns for c in hapsO_cols) else pd.DataFrame()
            
            if not hapsO_covs_df.empty:
                hapsO_covs1_vals = hapsO_covs_df.iloc[0:4, 0:4].values
                g11pre = hapsO_covs1_vals - np.eye(4) * (self.k2_matrix[0,0] / 2.0)
                g11 = g11pre[np.tril_indices_from(g11pre)]
                
                hapsO_covs2_vals = hapsO_covs_df.iloc[4:8, 4:8].values
                g22pre = hapsO_covs2_vals - np.eye(4) * (self.k2_matrix[1,1] / 2.0)
                g22 = g22pre[np.tril_indices_from(g22pre)]
                
                hapsO_covs12_vals = hapsO_covs_df.iloc[4:8, 0:4].values # cov(trait2_haps, trait1_haps)
                g12pre = hapsO_covs12_vals - np.eye(4) * (self.k2_matrix[0,1] / 2.0) # k2_matrix[0,1] is cov_bvo1_bvo2
                g12 = g12pre.flatten() # R uses lower.tri for symmetric, but this is full 4x4 block. R script seems to take mean of all elements.
                
                covG_val = np.array([[np.mean(g11), np.mean(g12)],
                                     [np.mean(g12), np.mean(g22)]])
            else:
                covG_val = np.full((2,2), np.nan)

            # hapsL.covs & covH
            hapsL_cols = ['TPL1','TML1','NTPL1','NTML1','TPL2','TML2','NTPL2','NTML2']
            hapsL_covs_df = self.phen_df[hapsL_cols].cov() if all(c in self.phen_df.columns for c in hapsL_cols) else pd.DataFrame()
            if not hapsL_covs_df.empty:
                hapsL_covs1_vals = hapsL_covs_df.iloc[0:4, 0:4].values
                h11pre = hapsL_covs1_vals - np.eye(4) * (self.k2_matrix[0,0] / 2.0) # Assuming k2_matrix is also for latent variance contribution from each allele
                h11 = h11pre[np.tril_indices_from(h11pre)]
                
                hapsL_covs2_vals = hapsL_covs_df.iloc[4:8, 4:8].values
                h22pre = hapsL_covs2_vals - np.eye(4) * (self.k2_matrix[1,1] / 2.0)
                h22 = h22pre[np.tril_indices_from(h22pre)]
                
                hapsL_covs12_vals = hapsL_covs_df.iloc[4:8, 0:4].values
                h12pre = hapsL_covs12_vals - np.eye(4) * (self.k2_matrix[0,1] / 2.0)
                h12 = h12pre.flatten()
                
                covH_val = np.array([[np.mean(h11), np.mean(h12)],
                                     [np.mean(h12), np.mean(h22)]])
            else:
                covH_val = np.full((2,2), np.nan)

            # covI (Covariance between Observed and Latent haplotypes)
            # R: (hapsOL1.covs <- cov(PHEN[,c('TPO1','TMO1','NTPO1','NTMO1','TPL1','TML1','NTPL1','NTML1')]))
            # R: (i11 <- as.vector(hapsOL1.covs[5:8,1:4])) # Cov(Latent1, Observed1)
            # This requires cov matrices like cov(O1_haps, L1_haps), cov(O1_haps, L2_haps) etc.
            # Simplified approach: use the full_cov_df
            o1_haps_cols = ['TPO1','TMO1','NTPO1','NTMO1']
            l1_haps_cols = ['TPL1','TML1','NTPL1','NTML1']
            o2_haps_cols = ['TPO2','TMO2','NTPO2','NTMO2']
            l2_haps_cols = ['TPL2','TML2','NTPL2','NTML2']
            
            if not full_cov_df.empty:
                i11_block = full_cov_df.loc[l1_haps_cols, o1_haps_cols].values.flatten() if all(c in full_cov_df.columns for c in o1_haps_cols+l1_haps_cols) else [np.nan]
                i22_block = full_cov_df.loc[l2_haps_cols, o2_haps_cols].values.flatten() if all(c in full_cov_df.columns for c in o2_haps_cols+l2_haps_cols) else [np.nan]
                # R: i12 is cov(O1,L2) -> L2_haps_cols vs O1_haps_cols
                i12_block = full_cov_df.loc[l2_haps_cols, o1_haps_cols].values.flatten() if all(c in full_cov_df.columns for c in o1_haps_cols+l2_haps_cols) else [np.nan]
                # R: i21 is cov(O2,L1) -> L1_haps_cols vs O2_haps_cols
                i21_block = full_cov_df.loc[l1_haps_cols, o2_haps_cols].values.flatten() if all(c in full_cov_df.columns for c in o2_haps_cols+l1_haps_cols) else [np.nan]
                
                covI_val = np.array([[np.mean(i11_block), np.mean(i12_block)],
                                     [np.mean(i21_block), np.mean(i22_block)]])
            else:
                covI_val = np.full((2,2), np.nan)

            # w, q, covF_calc, covE_calc
            # R: (covAF <- cov(PHEN[,c("F1","F2","AO1","AO2","AL1","AL2")]))
            # Note: F1,F2 here are Y-scaled environmental components from parents, NOT final Fy1,Fy2 of offspring.
            # AO1,AO2 are offspring's observed genetic values.
            cov_F_AO_AL_cols = ['F1','F2','AO1','AO2','AL1','AL2'] # Using F1,F2 from offspring which represent parental env.
            if all(c in self.phen_df.columns for c in cov_F_AO_AL_cols) and not self.phen_df[cov_F_AO_AL_cols].isnull().all().all():
                cov_F_AO_AL_df = self.phen_df[cov_F_AO_AL_cols].cov()
                w_val = cov_F_AO_AL_df.loc[['F1','F2'], ['AO1','AO2']].values # cov(F_parental_env, AO_offspring)
                q_val = cov_F_AO_AL_df.loc[['F1','F2'], ['AL1','AL2']].values # cov(F_parental_env, AL_offspring)
                covF_calc_val = cov_F_AO_AL_df.loc[['F1','F2'], ['F1','F2']].values
            else:
                w_val, q_val, covF_calc_val = np.full((2,2), np.nan), np.full((2,2), np.nan), np.full((2,2), np.nan)
            
            # covE_calc is Cov(E1, E2) of offspring, which are their unique environmental effects
            covE_cols = ['E1','E2']
            covE_calc_val = self.phen_df[covE_cols].cov().values if all(c in self.phen_df.columns for c in covE_cols) and not self.phen_df[covE_cols].isnull().all().all() else np.full((2,2), np.nan)


            # omega, gamma, thetaNT, thetaT
            # These are covariances between parental Y and offspring haplotypic values, or offspring Y and offspring haplotypic values
            # Example for omega: Cov(Parental Y, Offspring Transmitted Observed Haplotypes)
            # R: (omega.p.T <- cov(PHEN[,c("Y1P","Y2P","TPO1","TPO2")]))
            # This needs careful selection of columns from phen_df
            if not full_cov_df.empty:
                omega_p_T_val = full_cov_df.loc[['Y1P','Y2P'], ['TPO1','TPO2']].values if all(c in full_cov_df.columns for c in ['Y1P','Y2P','TPO1','TPO2']) else np.full((2,2),np.nan)
                omega_m_T_val = full_cov_df.loc[['Y1M','Y2M'], ['TMO1','TMO2']].values if all(c in full_cov_df.columns for c in ['Y1M','Y2M','TMO1','TMO2']) else np.full((2,2),np.nan)
                omega_T_val = (omega_p_T_val + omega_m_T_val) * 0.5
                
                omega_p_NT_val = full_cov_df.loc[['Y1P','Y2P'], ['NTPO1','NTPO2']].values if all(c in full_cov_df.columns for c in ['Y1P','Y2P','NTPO1','NTPO2']) else np.full((2,2),np.nan)
                omega_m_NT_val = full_cov_df.loc[['Y1M','Y2M'], ['NTMO1','NTMO2']].values if all(c in full_cov_df.columns for c in ['Y1M','Y2M','NTMO1','NTMO2']) else np.full((2,2),np.nan)
                omega_NT_val = (omega_p_NT_val + omega_m_NT_val) * 0.5
                omega_val = (omega_T_val + omega_NT_val) * 0.5 # R code was (omega.T + omega.NT)*.5

                # Gamma (Parental Y vs Offspring Latent Haps)
                gamma_p_T_val = full_cov_df.loc[['Y1P','Y2P'],['TPL1','TPL2']].values if all(c in full_cov_df.columns for c in ['Y1P','Y2P','TPL1','TPL2']) else np.full((2,2),np.nan)
                gamma_m_T_val = full_cov_df.loc[['Y1M','Y2M'],['TML1','TML2']].values if all(c in full_cov_df.columns for c in ['Y1M','Y2M','TML1','TML2']) else np.full((2,2),np.nan)
                gamma_T_val = (gamma_p_T_val + gamma_m_T_val) * 0.5
                
                gamma_p_NT_val = full_cov_df.loc[['Y1P','Y2P'],['NTPL1','NTPL2']].values if all(c in full_cov_df.columns for c in ['Y1P','Y2P','NTPL1','NTPL2']) else np.full((2,2),np.nan)
                gamma_m_NT_val = full_cov_df.loc[['Y1M','Y2M'],['NTML1','NTML2']].values if all(c in full_cov_df.columns for c in ['Y1M','Y2M','NTML1','NTML2']) else np.full((2,2),np.nan)
                gamma_NT_val = (gamma_p_NT_val + gamma_m_NT_val) * 0.5
                gamma_val = (gamma_T_val + gamma_NT_val) * 0.5

                # Theta (Offspring Y vs Offspring Non-Transmitted/Transmitted Observed Haps)
                theta_NTp_val = full_cov_df.loc[['Y1','Y2'],['NTPO1','NTPO2']].values if all(c in full_cov_df.columns for c in ['Y1','Y2','NTPO1','NTPO2']) else np.full((2,2),np.nan)
                theta_NTm_val = full_cov_df.loc[['Y1','Y2'],['NTMO1','NTMO2']].values if all(c in full_cov_df.columns for c in ['Y1','Y2','NTMO1','NTMO2']) else np.full((2,2),np.nan)
                thetaNT_val = theta_NTp_val + theta_NTm_val

                theta_Tp_val = full_cov_df.loc[['Y1','Y2'],['TPO1','TPO2']].values if all(c in full_cov_df.columns for c in ['Y1','Y2','TPO1','TPO2']) else np.full((2,2),np.nan)
                theta_Tm_val = full_cov_df.loc[['Y1','Y2'],['TMO1','TMO2']].values if all(c in full_cov_df.columns for c in ['Y1','Y2','TMO1','TMO2']) else np.full((2,2),np.nan)
                thetaT_val = theta_Tp_val + theta_Tm_val
            else: # full_cov_df is empty
                omega_val, gamma_val, thetaNT_val, thetaT_val = np.full((2,2),np.nan), np.full((2,2),np.nan), np.full((2,2),np.nan), np.full((2,2),np.nan)


            summary_this_gen = {'GEN': r_cur_gen_num, 'NUM.CVs': self.num_cvs,
                                'MATE.COR': self.am_list[r_cur_gen_num].tolist() if r_cur_gen_num < len(self.am_list) else None,
                                'POPSIZE': len(self.phen_df)}
            
            for comp_s in ['AOy', 'ALy', 'Fy', 'Ey', 'Y']: # Use Y-scaled components for variances
                 var_comp_name_s = 'V' + comp_s.replace('y','').upper() # VAO, VAL, VF, VE, VP
                 cols_to_cov_s = [comp_s+'1', comp_s+'2'] if comp_s != 'Y' else ['Y1', 'Y2']
                 
                 if all(c in self.phen_df.columns for c in cols_to_cov_s):
                     df_subset_s = self.phen_df[cols_to_cov_s].dropna()
                     summary_this_gen[var_comp_name_s] = np.cov(df_subset_s, rowvar=False).tolist() if len(df_subset_s) >=2 else np.array([[np.nan,np.nan],[np.nan,np.nan]]).tolist()
                 else:
                     summary_this_gen[var_comp_name_s] = np.array([[np.nan,np.nan],[np.nan,np.nan]]).tolist()
            
            vp_diag_s = np.diag(np.array(summary_this_gen.get('VP', [[np.nan,np.nan]])))
            vao_diag_s = np.diag(np.array(summary_this_gen.get('VAO', [[np.nan,np.nan]])))
            val_diag_s = np.diag(np.array(summary_this_gen.get('VAL', [[np.nan,np.nan]])))
            with np.errstate(divide='ignore', invalid='ignore'):
                summary_this_gen['h2'] = ((vao_diag_s + val_diag_s) / vp_diag_s).tolist()
                summary_this_gen['h2.obs'] = (vao_diag_s / vp_diag_s).tolist()
                summary_this_gen['h2.lat'] = (val_diag_s / vp_diag_s).tolist()

            summary_this_gen.update({
                'covY': full_cov_df.loc[['Y1P','Y2P','Y1M','Y2M','Y1','Y2'], ['Y1P','Y2P','Y1M','Y2M','Y1','Y2']].values.tolist() if 'Y1P' in full_cov_df and 'Y1' in full_cov_df else np.nan,
                'covG': covG_val.tolist(), 'covH': covH_val.tolist(), 'covI': covI_val.tolist(),
                'w': w_val.tolist(), 'q': q_val.tolist(),
                'covF': covF_calc_val.tolist(), 'covE': covE_calc_val.tolist(),
                'hapsO.covs': hapsO_covs_df.to_dict() if not hapsO_covs_df.empty else np.nan,
                'hapsL.covs': hapsL_covs_df.to_dict() if not hapsL_covs_df.empty else np.nan,
                'omega': omega_val.tolist(), 'gamma': gamma_val.tolist(),
                'thetaNT': thetaNT_val.tolist(), 'thetaT': thetaT_val.tolist()
            })
            self.summary_results.append(summary_this_gen)

            if self.save_each_gen:
                self.history['MATES'].append(mates)
                self.history['PHEN'].append(self.phen_df.copy())
                self.history['XO'].append(self.xo.copy())
                self.history['XL'].append(self.xl.copy())
            
            if self.save_covs:
                self.covariances_log.append(full_cov_df.round(3).to_dict() if not full_cov_df.empty else None)

            print(f"--- Generation {r_cur_gen_num} Processing Done ---")
            
        return {
            'SUMMARY.RES': self.summary_results,
            'XO': self.xo, 'XL': self.xl, 'PHEN': self.phen_df,
            'HISTORY': self.history,
            'COVARIANCES': self.covariances_log
        }