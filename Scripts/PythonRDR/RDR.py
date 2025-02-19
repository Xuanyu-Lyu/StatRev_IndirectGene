# This is a script to conduct RDR model on simulated data from Gene-Evolve.
# author: Xuanyu Lyu

import numpy as np
import pandas as pd
import statsmodels.api as sm

# read the offspring phenotype and genotype data
df_offspring = pd.read_table("Her.5-AM.4-Lat.5-VF0-offspring.txt", sep="\t")

# extract the phenotype data
Y_offspring = df_offspring["Y"].values
df_gene_o = df_offspring.iloc[:, 4:]

# read the mother and father genotype data
df_gene_m = pd.read_table("Her.5-AM.4-Lat.5-VF0-mother.txt", sep="\t")
df_gene_f = pd.read_table("Her.5-AM.4-Lat.5-VF0-father.txt", sep="\t")

# standardize all columns in the genotype data to have mean 0 and variance 1
df_gene_o_std = (df_gene_o - df_gene_o.mean()) / df_gene_o.std()


#  add the mother and father genotype data into parental genotype
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

# run the regression
X = np.column_stack((v_R_SNP_o, v_R_SNP_p, v_R_SNP_op))
X = sm.add_constant(X)
model = sm.OLS(v_Y, X)
results = model.fit()
print(results.summary())

