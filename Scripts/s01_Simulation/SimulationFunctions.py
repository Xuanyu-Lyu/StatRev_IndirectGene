# This is a Python script that contains functions to perform multivariate gene-evolve simulations. 
# This is an adapted version of the original R script, rewritten in Python.

# The following functions are Python equivalents of the provided R functions.
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.linalg import svd as scipy_svd # Or numpy.linalg.svd


# R: is.even <- function(x) {x%%2 == 0}
def is_even(x: int) -> bool:
    """Checks if an integer is even."""
    return x % 2 == 0

# R: is.odd <- function(x) {x%%2 != 0}
def is_odd(x: int) -> bool:
    """Checks if an integer is odd."""
    return x % 2 != 0

def next_smallest(x_list: list, value_thresh: float, return_index: bool = True):
    """
    Finds the largest element in x_list that is smaller than value_thresh.
    This is equivalent to finding the element in the sorted version of x_list
    at the highest possible index i such that x_sorted[i] < value_thresh.

    Args:
        x_list: A list of numbers.
        value_thresh: The threshold value.
        return_index: If True, returns the 0-based index of this element
                      in the sorted version of x_list.
                      If False, returns the value of this element.

    Returns:
        The 0-based index or the value, or None if no element in x_list is smaller than value_thresh.
    """
    # Sort a copy of the input list.
    x2_sorted = sorted(list(x_list))

    # Find 0-based indices of elements in x2_sorted that are less than value_thresh
    potential_indices = [i for i, val in enumerate(x2_sorted) if val < value_thresh]

    if not potential_indices:
        # R's `max(which(condition_is_all_false))` results in -Inf.
        # Accessing an array with -Inf index would be an error in R.
        # Pythonic behavior is to return None if no such element is found.
        return None

    # The R equivalent `max(which(x2<value))` gives the largest 1-based index.
    # We use the largest 0-based index here.
    actual_index_in_x2_sorted = max(potential_indices)

    if return_index:
        return actual_index_in_x2_sorted
    else:
        return x2_sorted[actual_index_in_x2_sorted]

def cor2cov(X, var1: float, var2: float):
    """
    Converts a 2x2 correlation matrix X to a covariance matrix using specified variances.

    Args:
        X: A 2x2 matrix-like object (e.g., list of lists or numpy array)
           representing the correlation matrix.
        var1: Variance of the first variable. Must be non-negative.
        var2: Variance of the second variable. Must be non-negative.

    Returns:
        A 2x2 numpy array representing the covariance matrix.

    Raises:
        ValueError: If X is not a 2x2 matrix or if variances are negative.
    """
    X_np = np.asarray(X)
    if X_np.shape != (2, 2):
        raise ValueError("Input matrix X must be a 2x2 matrix.")
    if var1 < 0 or var2 < 0:
        raise ValueError("Variances cannot be negative.")

    # sd_mat is the diagonal matrix of standard deviations
    sd_mat = np.array([
        [np.sqrt(var1), 0],
        [0,             np.sqrt(var2)]
    ])

    # Covariance matrix = sd_mat * Correlation_matrix * sd_mat
    cov_mat = sd_mat @ X_np @ sd_mat
    return cov_mat

## the following functions are the main functions to perform the simulation

def prepare_CV(n_CV, rg, maf_min, maf_max, maf_dist = "uniform"):
    """
    Prepares a list of causal variants (CV) for the simulation.
    
    Args:
        n_CV: Number of causal variants.
        rg: Genetic correlation between two traits. 
        maf_min: Minimum minor allele frequency.
        maf_max: Maximum minor allele frequency.
    Returns:
        A DataFrame containing the minor allele frequencies (MAF) and effect sizes for two traits.
    Check: It is working properly.
    """
    
    # Generate minor allele frequencies (MAF) based on the specified distribution
    if maf_dist == "uniform":
        maf = np.random.uniform(maf_min, maf_max, n_CV)
    elif maf_dist == "normal":
        maf = np.random.normal((maf_min + maf_max) / 2, (maf_max - maf_min) / 6, n_CV)
        maf = np.clip(maf, maf_min, maf_max)  # Ensure MAF is within bounds
    else:
        raise ValueError("Invalid MAF distribution specified.")

    # Calculate the genetic correlation matrix
    cor_mat = np.array([[1, rg], [rg, 1]])
    
    # get the genotypic variance for each allele
    # The genotypic variance is calculated as 2 * maf * (1 - maf)
    var_gene = 2 * maf * (1 - maf)
    # get the effect sizes for each allele using multivariate normal distribution based on the correlation matrix
    alpha_pre = np.random.multivariate_normal(mean=np.zeros(2), cov=cor_mat, size=n_CV)
    scaler = np.sqrt(1/(var_gene * n_CV))
    alpha_final = alpha_pre * scaler[:, np.newaxis]  # Scale the effect sizes
    
    # return the MAF and two effect sizes columns
    return pd.DataFrame({
        "maf": maf,
        "alpha1": alpha_final[:, 0],
        "alpha2": alpha_final[:, 1]
    })
        


def GE_simulate(CV, n_gen, n_pop, delta_mat, a_mat, fm_mat, fp_mat, VE_mat, VY_mat, k2_mat, matecor_list, avoid_inbreeding = True, save_intermediate = False, save_covs, seed=62):
    """
    Simulates the evolution of a population over multiple generations.

    Args:
        CV: Coefficient of variation.
        n_gen: Number of generations to simulate.
        n_pop: Number of individuals in the population.
        delta_mat: Matrix of genetic variances.
        a_mat: Matrix of additive genetic variances.
        fm_mat: Matrix of maternal vertical transmission.
        fp_mat: Matrix of paternal vertical transmission.
        VE_mat: Matrix of (residual/unique) environmental variances.
        VY_mat: Matrix of phenotypic variances.
        k2_mat: Matrix of standardized genetic covariance.
        matecor_list: List of mating correlations.
        avoid_inbreeding: Boolean indicating whether to avoid inbreeding.
        save_intermediate: Boolean indicating whether to save intermediate results.
        save_covs: Boolean indicating whether to save covariance matrices.
        seed: Random seed for reproducibility.

    Returns:
        A DataFrame containing the simulation results.

    Raises:
        ValueError: If any input matrix is not square or has invalid dimensions.
    """
    # Set random seed for reproducibility
    np.random.seed(seed)
    
    # create a vector for population sizes for all generations
    n_pop_vec = np.repeat(n_pop, n_gen)
    # get the number of CVs from the input CV dataframe
    n_CV = CV.shape[0]
    
    ##############################
    # Step 1: create observed and latent genotypes of the founder population (Generation 0)
    ##############################
    n_total_cv = n_CV * n_pop
    # create a matrix of genotypes for the founder population
    genotype_obs = np.random.binomial(2, CV["maf"].values, size=(n_pop, n_CV))
    # create a matrix of latent genotypes for the founder population
    genotype_latent = np.random.binomial(2, CV["maf"].values, size=(n_pop, n_CV))
    
    ##############################
    # Step 2: generate the phenotypes for the founder population (Generation 0), prepare output files for the founder population
    ##############################
    # create the standardized genotypic values for the two traits; note these have not been scaled by the path coefficients yet.
    AO_std = genotype_obs @ CV[["alpha1", "alpha2"]]
    AL_std = genotype_latent @ CV[["alpha1", "alpha2"]]
    # create empirical k*2 and j*2 matrices, may not be needed
    #k2_mat_emp = AO.cov()
    #j2_mat_emp = AL.cov()
    
    # Create components of the phenotypic variance for the founder population
    # F component: environmental variance from maternal and paternal vertical transmission
    F = np.random.multivariate_normal(mean=np.zeros(2), cov=VY_mat, size=n_pop) @ fm_mat.T + np.random.multivariate_normal(mean=np.zeros(2), cov=VY_mat, size=n_pop) @ fp_mat.T
    # E component: environmental variance from residual/unique environmental effects
    E = np.random.multivariate_normal(mean=np.zeros(2), cov=VE_mat, size=n_pop) 
    # Genetic components
    AO = AO_std @ delta_mat.T
    AL = AL_std @ delta_mat.T
    Y = AO + AL + F + E
    
    # generate ID for the founder population
    ID_vec = np.random.choice(np.arange(10000000, 100000000), size=n_pop * 7, replace=False).reshape(n_pop, 7)
    
    # Create SEXVEC (gender vector for Generation 0)
    if is_even(n_pop):
        # Ensure integer division for indices
        sex_vec = np.concatenate([np.zeros(n_pop // 2, dtype=int), np.ones(n_pop // 2, dtype=int)])
    else:
        sex_vec = np.concatenate([np.zeros(n_pop // 2, dtype=int), np.ones(n_pop // 2, dtype=int), np.array([0], dtype=int)])
    
    # combine all the info into a dataframe
    Phen = pd.DataFrame({
        "ID": ID_vec[:,0],
        "Father.ID": ID_vec[:,1],
        "Mother.ID": ID_vec[:,2],
        "Fathers.Father.ID": ID_vec[:,3],
        "Fathers.Mother.ID": ID_vec[:,4],
        "Mothers.Father.ID": ID_vec[:,5],
        "Mothers.Mother.ID": ID_vec[:,6],
        "Sex":sex_vec, # 0 for female, 1 for male
        "AO_std1": AO_std[:,0],
        "AO_std2": AO_std[:,1],
        "AL_std1": AL_std[:,0],
        "AL_std2": AL_std[:,1],
        "AO1": AO[:,0],
        "AO2": AO[:,1],
        "AL1": AL[:,0],
        "AL2": AL[:,1],
        "F1": F[:,0],
        "F2": F[:,1],
        "E1": E[:,0],
        "E2": E[:,1],
        "Y1": Y[:,0],
        "Y2": Y[:,1],
        "Y1P": np.full(n_pop, np.nan),
        "Y2P": np.full(n_pop, np.nan),
        "Y1M": np.full(n_pop, np.nan),
        "Y2M": np.full(n_pop, np.nan),
        "F1P": np.full(n_pop, np.nan), # I'm not sure if we need F1P to F2M
        "F2P": np.full(n_pop, np.nan),
        "F1M": np.full(n_pop, np.nan),
        "F2M": np.full(n_pop, np.nan),
        "TPO1": np.full(n_pop, np.nan),
        "TPO2": np.full(n_pop, np.nan),
        "TMO1": np.full(n_pop, np.nan),
        "TMO2": np.full(n_pop, np.nan),
        "NTPO1": np.full(n_pop, np.nan),
        "NTPO2": np.full(n_pop, np.nan),
        "NTMO1": np.full(n_pop, np.nan),
        "NTMO2": np.full(n_pop, np.nan),
        "TPL1": np.full(n_pop, np.nan),
        "TPL2": np.full(n_pop, np.nan),
        "TML1": np.full(n_pop, np.nan),
        "TML2": np.full(n_pop, np.nan),
        "NTPL1": np.full(n_pop, np.nan),
        "NTPL2": np.full(n_pop, np.nan),
        "NTML1": np.full(n_pop, np.nan),
        "NTML2": np.full(n_pop, np.nan)
    })
    # save the phenotypic data and covariance matrices for founder population if save_intermediate is True
    # also create a dictionary to match the index with the saved contents
    

    l_Mates = []
    l_Phen = []
    l_genotype_obs = []
    l_genotype_latent = []
    l_cov = []
    if save_intermediate:
        l_Mates.append(np.full((n_pop, 2), np.nan))
        l_Phen.append(Phen)
        l_genotype_obs.append(genotype_obs)
        l_genotype_latent.append(genotype_latent)
        l_cov.append(np.full((2, 2), np.nan))
        l_History = [l_Mates, l_Phen, l_genotype_obs, l_genotype_latent, l_cov]
        # dictionary to match the index with the saved contents
        d_History = {
            "Mates": 0,
            "Phen": 1,
            "genotype_obs": 2,
            "genotype_latent": 3,
            "Cov": 4
        }

    # initialize a list to save important summary statistics for each generation
    dict_summary = {
        "Generation": 0,
        "n_CV": 1,
        "n_pop": 2,
        "matecor": 3,
        "VAO": 4,
        "VAL": 5,
        "VF": 6,
        "VE": 7,
        "VY": 8,
        "heritability": 9,
        "h2_observed": 10,
        "h2_latent": 11,
        "covY": 12,
        "corY": 13,
        "covG": 14,
        "covH": 15,
        "covI": 16,
        "w": 17,
        "v": 18,
        "omega": 19,
        "gamma": 20,
        "thetaNT": 21,
        "thetaT": 22
    }
        
    l_summary_current = [0, n_CV, n_pop, matecor_list[0], 
                         Phen.loc[:,["AO1", "AO2"]].cov(), Phen.loc[:,["AL1", "AL2"]].cov(), Phen.loc[:,["F1", "F2"]].cov(), Phen.loc[:,["E1", "E2"]].cov(), Phen.loc[:,["Y1", "Y2"]].cov(),
                         (Phen.loc[:,["AO1", "AO2"]].cov() + Phen.loc[:,["AL1", "AL2"]].cov())/Phen.loc[:,["Y1", "Y2"]].cov(),
                         Phen.loc[:,["AO1", "AO2"]].cov()/Phen.loc[:,["Y1", "Y2"]].cov(),
                         Phen.loc[:,["AL1", "AL2"]].cov()/Phen.loc[:,["Y1", "Y2"]].cov(),
                         Phen.loc[:,["Y1","Y2", "Y1P", "Y2P", "Y1M", "Y2M"]].cov(),
                         Phen.loc[:,["Y1","Y2", "Y1P", "Y2P", "Y1M", "Y2M"]].corr(),
                         np.full((2, 2), np.nan),
                         np.full((2, 2), np.nan),
                         np.full((2, 2), np.nan),
                         np.full((2, 2), np.nan),
                         np.full((2, 2), np.nan),
                         np.full((2, 2), np.nan),
                         np.full((2, 2), np.nan),
                         np.full((2, 2), np.nan),
                         np.full((2, 2), np.nan)
                        ]
    l_summary = [l_summary_current]
    
    #################################
    # Step 3: Loop through generations
    #################################
    
                                 
    
def assortative_mating(current_phen, matecor, avoid_inbreeding, type="phenotypic"):
    """
    Performs assortative mating based on the current phenotypic data and mating correlation.

    Args:
        current_phen: DataFrame containing the current phenotypic data.
        matecor: Mating correlation.
        avoid_inbreeding: Boolean indicating whether to avoid inbreeding.
        type: Type of mating correlation to use. Default is "phenotypic". Other options are "genotypic" and "environmental".

    Returns:
        A DataFrame containing the updated mating pairs.
    """
    #### first step: specify mating pairs based on the mating correlation
    n = current_phen.shape[0]
    n_males = current_phen['Sex'].sum()
    n_females = n - n_males
    
    # spilit the phen dataframe to two dataframes for males and females where only have the ID columns and the columns used for mating
    if type == "phenotypic":
        columns_for_mating = ["ID", "Y1", "Y2"]
        males_condition = current_phen["Sex"] == 1
        females_condition = current_phen["Sex"] == 0
        males_phen = current_phen.loc[males_condition, columns_for_mating]
        females_phen = current_phen.loc[females_condition, columns_for_mating]
        # change the column names to be consistent across three ways
        males_phen = males_phen.rename(columns={"Y1": "mating1", "Y2": "mating2"})
        females_phen = females_phen.rename(columns={"Y1": "mating1", "Y2": "mating2"})
    elif type == "social":
        columns_for_mating = ["ID", "F1", "F2", "E1", "E2"]
        males_condition = current_phen["Sex"] == 1
        females_condition =current_phen["Sex"] == 0
        males_phen = current_phen.loc[males_condition, columns_for_mating]
        females_phen = current_phen.loc[females_condition, columns_for_mating]
        # add the F and E components to be the mating columns
        males_phen["mating1"] = males_phen["F1"] + males_phen["E1"]
        males_phen["mating2"] = males_phen["F2"] + males_phen["E2"]
        females_phen["mating1"] = females_phen["F1"] + females_phen["E1"]
        females_phen["mating2"] = females_phen["F2"] + females_phen["E2"]
        # drop the F and E columns
        males_phen = males_phen.drop(columns=["F1", "F2", "E1", "E2"])
        females_phen = females_phen.drop(columns=["F1", "F2", "E1", "E2"])
    elif type == "genotypic":
        columns_for_mating = ["ID", "AO1", "AO2", "AL1", "AL2"]
        males_condition = current_phen["Sex"] == 1
        females_condition =current_phen["Sex"] == 0
        males_phen = current_phen.loc[males_condition, columns_for_mating]
        females_phen = current_phen.loc[females_condition, columns_for_mating]
        # add the AO and AL components to be the mating columns
        males_phen["mating1"] = males_phen["AO1"] + males_phen["AL1"]
        males_phen["mating2"] = males_phen["AO2"] + males_phen["AL2"]
        females_phen["mating1"] = females_phen["AO1"] + females_phen["AL1"]
        females_phen["mating2"] = females_phen["AO2"] + females_phen["AL2"]
        # drop the AO and AL columns
        males_phen = males_phen.drop(columns=["AO1", "AO2", "AL1", "AL2"])
        females_phen = females_phen.drop(columns=["AO1", "AO2", "AL1", "AL2"])
    
    # make equal numbers of males and females
    if n_males > n_females:
        drop_idx = np.random.choice(males_phen.index, size=n_males - n_females)
        males_phen = males_phen.drop(drop_idx)
    elif n_females > n_males:
        drop_idx = np.random.choice(females_phen.index, size=n_females-n_males)
        females_phen = females_phen.drop(drop_idx)
    


    

        

        
        
    
    
    # This function should return a DataFrame with updated mating pairs
    pass
    
    
    
    
 
    
