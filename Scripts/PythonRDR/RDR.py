# This is a script to conduct RDR model on simulated data from Gene-Evolve.
# author: Xuanyu Lyu

import numpy as np
import pandas as pd
import statsmodels.api as sm
import numpy as np
from scipy.optimize import minimize

# read the offspring phenotype and genotype data
df_offspring = pd.read_table("Her.5-AM.4-Lat0-VF0.15-offspring-nosib.txt", sep="\t")

# extract the phenotype data
Y_offspring = df_offspring["Y"].values
# standardize the phenotype data
#Y_offspring = (Y_offspring - Y_offspring.mean()) / Y_offspring.std()
df_gene_o = df_offspring.iloc[:, 4:]

# read the mother and father genotype data
df_gene_m = pd.read_table("Her.5-AM.4-Lat0-VF0.15-mother-nosib.txt", sep="\t")
df_gene_f = pd.read_table("Her.5-AM.4-Lat0-VF0.15-father-nosib.txt", sep="\t")

# standardize all columns in the genotype data to have mean 0 and variance 1
df_gene_o_std = (df_gene_o - df_gene_o.mean()) / df_gene_o.std()


#  add the mother and father genotype data into parental genotype
# check how to deal with duos
df_gene_p = df_gene_m + df_gene_f
# standardize them
df_gene_p_std = (df_gene_p - df_gene_p.mean()) / (df_gene_p.std()/np.sqrt(2))

# calculate the relatedness matrix
R_SNP_o = np.dot(df_gene_o_std, df_gene_o_std.T)/df_gene_o_std.shape[1]
R_SNP_p = np.dot(df_gene_p_std, df_gene_p_std.T)/(df_gene_p_std.shape[1]*2)
R_SNP_op = (np.dot(df_gene_o_std, df_gene_p_std.T) + np.dot(df_gene_p_std, df_gene_o_std.T))/(df_gene_p_std.shape[1]*2)

# get the phenotypic covariance matrix
Y_scale = Y_offspring - Y_offspring.mean()
Y_scale = pd.DataFrame(Y_scale)
COV_Y = Y_scale @ Y_scale.T

# regress the off-diagenal phenotypic covariance matrix on the relatedness matrices
# extract the lower triangle of Y, excluding the diagonal
v_Y = COV_Y.values[np.tril_indices(COV_Y.shape[0], -1)]
# extract the lower triangle of R, excluding the diagonal
v_R_SNP_o = R_SNP_o[np.tril_indices(R_SNP_o.shape[0], -1)]
v_R_SNP_p = R_SNP_p[np.tril_indices(R_SNP_p.shape[0], -1)]
v_R_SNP_op = R_SNP_op[np.tril_indices(R_SNP_op.shape[0], -1)]
# get a

# run a maximum likelihood model that Y should follow a multivariate normal distribution with mean 0 and covariance
# matrix that is a linear combination of the relatedness matrices and the error matrix
# Y ~ R_SNP_o + R_SNP_p + R_SNP_op + I

def build_Sigma(params, R_snp, R_par, R_op):
    """
    params: [mu, alpha_g, alpha_e_g, alpha_g_e, alpha_sigma]
    R_snp:        (n x n)
    R_par:        (n x n)
    R_op:         (n x n)
    return:       covariance matrix Sigma
    """
    mu        = params[0]
    alpha_g   = params[1]
    alpha_e_g = params[2]
    alpha_g_e = params[3]
    alpha_sig = params[4]

    v_g     = np.exp(alpha_g)
    v_e_g   = np.exp(alpha_e_g)
    c_g_e   = np.exp(alpha_g_e)
    sig2    = np.exp(alpha_sig)

    # build covariance
    Sigma = v_g   * R_snp \
          + v_e_g * R_par \
          + c_g_e * R_op   \
          + sig2  * np.eye(R_snp.shape[0])
    return Sigma, mu

def neg_log_lik(params, y, R_snp, R_par, R_op):
    """
    Negative log-likelihood for
      y ~ N(mu, v_g R_snp + v_e_g R_par + c_g_e R_op + sigma^2 I).
    """
    n = len(y)
    Sigma, mu = build_Sigma(params, R_snp, R_par, R_op)
    
    # Compute y - mu
    ym = y - mu
    
    # slogdet for log|Sigma|
    sign, logdet = np.linalg.slogdet(Sigma)
    if sign <= 0:
        return 1e15  # penalize invalid Sigma
    
    # Solve for quadratic form
    invSy = np.linalg.solve(Sigma, ym)
    quadform = ym.dot(invSy)
    
    # - log L
    nll = 0.5 * (n * np.log(2 * np.pi) + logdet + quadform)
    return nll

# Example usage:
# y is your phenotype vector
# R_snp, R_par, R_op are (n x n) arrays
# Choose initial guesses:
mu0      = np.mean(v_Y)
alpha_g0 = np.log(0.1)
alpha_e0 = np.log(0.1)
alpha_c0 = np.log(0.1)
alpha_s0 = np.log(0.1)

init_params = np.array([mu0, alpha_g0, alpha_e0, alpha_c0, alpha_s0])
Y_scale_array = Y_scale.values.flatten()

# Define bounds for the parameters to prevent overflow
bounds = [(-100, 100)] * len(init_params)  # Example bounds, adjust as needed

res = minimize(
    fun=neg_log_lik,
    x0=init_params,
    args=(Y_scale_array, R_SNP_o , R_SNP_p, R_SNP_op),
    method='SLSQP',
    jac=None,  # You can provide a gradient if available
    bounds=bounds,
    options={'disp': True, 'maxls': 20, 'gtol': 1e-5}  # Increase maxls, adjust gtol
)

if res.success:
    p_hat = res.x
    mu_hat      = p_hat[0]
    v_g_hat     = np.exp(p_hat[1])
    v_e_g_hat   = np.exp(p_hat[2])
    c_g_e_hat   = np.exp(p_hat[3])
    sigma2_hat  = np.exp(p_hat[4])
    
    print("MLE results:")
    print("mu        =", mu_hat)
    print("v_g       =", v_g_hat)
    print("v_e_g     =", v_e_g_hat)
    print("c_g_e     =", c_g_e_hat)
    print("sigma^2   =", sigma2_hat)
else:
    print("Optimization failed:", res.message)
    
# run the regression
X = np.column_stack((v_R_SNP_o, v_R_SNP_p, v_R_SNP_op))
X = sm.add_constant(X)
model = sm.OLS(v_Y, X)
results = model.fit()
print("Her.5-AM0.4-Lat0-VF0.15-nosib\n")
print(results.summary())

