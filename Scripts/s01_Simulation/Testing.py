import sys
import os

# Define the absolute path to your project's root directory
# This is the directory that CONTAINS the 'Scripts' folder
project_root = '/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/coding/PyProject/StatRev_IndirectGene'

# Add the project root to sys.path if it's not already there
if project_root not in sys.path:
    sys.path.append(project_root)

# This is a temperary file for testing purposes.
from Scripts.s01_Simulation.SimulationFunctions import *

# test prepareCV function
n_gen = 40000
n_CV = 2000
rg = 0.5
maf_min = .1
maf_max = 0.5

CV_df = prepare_CV(n_CV, rg, maf_min, maf_max)

# test the binomial generation function
# n_total_cv = 1000
geno_df_test = np.random.binomial(2, CV_df["maf"].values, size=(n_gen, n_CV))

AO = geno_df_test @ CV_df[["alpha1", "alpha2"]]
AL = geno_df_test @ CV_df[["alpha1", "alpha2"]]

AO.cov()
