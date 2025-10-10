# RDR Regression Pipeline

This directory contains the complete pipeline for running Relatedness Disequilibrium Regression (RDR) analysis using Python regression instead of GCTA HE regression.

## Pipeline Overview

The RDR regression pipeline consists of four main components:

1. **Individual Run Analysis**: Analyzes a single simulation run
2. **Batch Job Submission**: Submits jobs for all conditions to the cluster
3. **Results Aggregation**: Combines results from all runs
4. **Aggregation Job Submission**: Submits the aggregation job to the cluster

## Files in the Pipeline

### Core Analysis Scripts
- `run_rdr_regression.py` - Main Python script that runs RDR regression on a single simulation run
- `run_rdr_regression.sh` - Shell wrapper script for running the Python analysis

### Job Submission Scripts
- `submit_rdr_regression.sh` - Full SLURM job array script for all conditions (500 runs)
- `submit_rdr_regression_test.sh` - Test version with limited conditions (50 runs)

### Aggregation Scripts
- `aggregate_rdr_regression_results.py` - Python script to aggregate all JSON results
- `aggregate_rdr_regression.sh` - SLURM job script for running aggregation

## Usage Instructions

### Step 1: Submit Analysis Jobs

For testing (limited conditions):
```bash
sbatch submit_rdr_regression_test.sh
```

For full analysis (all conditions, 1000 runs):
```bash
sbatch submit_rdr_regression.sh
```

### Step 2: Monitor Jobs
```bash
squeue -u $USER
```

### Step 3: Aggregate Results (after jobs complete)
```bash
sbatch aggregate_rdr_regression.sh
```

## Method Description

The RDR regression method follows the approach from `PrelimResults.ipynb`:

1. **Data Loading**: Loads offspring and parent genotype/phenotype data from generation 20
2. **Sampling**: Creates independent trios (one offspring per father) with specified sample size
3. **GRM Calculation**: Computes three relatedness matrices:
   - `R_offspring`: Offspring-offspring relatedness
   - `R_parental`: Parental-parental relatedness (sum of father + mother genotypes)
   - `R_cross`: Cross relatedness between offspring and parents
4. **Regression**: Regresses off-diagonal phenotypic covariances on the three relatedness matrices

The regression model is:
```
Cov(Y_i, Y_j) = β₀ + β₁×R_offspring(i,j) + β₂×R_parental(i,j) + β₃×R_cross(i,j) + ε
```

## Output Format

### Individual Results
Each run produces text files with detailed regression results (similar to GCTA .HEreg format):
- `{run_id}_RDR_regression_N{sample_size}_Y1.RDRreg`
- `{run_id}_RDR_regression_N{sample_size}_Y2.RDRreg`

### Aggregated Results
- `aggregated_rdr_regression_results_detailed.tsv` - All individual results
- `aggregated_rdr_regression_results_summary.tsv` - Summary statistics by condition

## Key Differences from GCTA HE Regression

1. **Method**: Uses OLS regression instead of Haseman-Elston regression
2. **Speed**: Generally faster than GCTA (no GRM file I/O)
3. **Memory**: Lower memory requirements
4. **Output**: More detailed statistical output (R², AIC, BIC, etc.)
5. **Flexibility**: Easier to modify and extend

## Resource Requirements

- **CPU**: 4 cores per job (reduced from GCTA's 8)
- **Memory**: 32GB per job (reduced from GCTA's 140GB)  
- **Time**: 6 hours per job (reduced from GCTA's 12 hours)
- **Storage**: Results stored as .RDRreg text files in `/projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/RDR_Regression_Results/`

## Troubleshooting

1. **Permission Issues**: Ensure scripts are executable:
   ```bash
   chmod +x run_rdr_regression.sh
   chmod +x run_rdr_regression.py
   ```

2. **Environment Issues**: Ensure conda environment is activated:
   ```bash
   conda activate /projects/xuly4739/general_env
   ```

3. **Path Issues**: Check that data paths in scripts match your cluster setup

4. **Memory Issues**: If jobs fail with OOM errors, increase `--mem` in submission scripts