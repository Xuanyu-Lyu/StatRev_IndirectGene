#!/usr/bin/env python3
"""
Aggregate RDR regression results from multiple simulation runs
This script processes JSON result files and creates summary statistics

Author: Xuanyu Lyu
"""

import os
import glob
import pandas as pd
import numpy as np
from datetime import datetime

def load_result_from_rdreg_file(filepath):
    """
    Load a single RDR regression result from .RDRreg file
    """
    try:
        results = {}
        with open(filepath, 'r') as f:
            content = f.read()
            
        lines = content.strip().split('\n')
        
        # Parse metadata
        for line in lines:
            if line.startswith('Run ID:'):
                results['run_id'] = line.split(':', 1)[1].strip()
            elif line.startswith('Trait:'):
                results['trait'] = line.split(':', 1)[1].strip()
            elif line.startswith('Sample size:'):
                results['sample_size'] = int(line.split(':', 1)[1].strip())
            elif line.startswith('Number of pairs:'):
                results['n_pairs'] = int(line.split(':', 1)[1].strip())
            elif line.startswith('Phenotypic variance:'):
                results['phenotypic_variance'] = float(line.split(':', 1)[1].strip())
            elif line.startswith('Analysis timestamp:'):
                results['analysis_timestamp'] = line.split(':', 1)[1].strip()
        
        # Parse regression results
        in_regression_section = False
        for line in lines:
            if line == "Regression Results:":
                in_regression_section = True
                continue
            elif line.startswith("Model Statistics:"):
                in_regression_section = False
                continue
            elif in_regression_section and '\t' in line and not line.startswith('Component'):
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    component = parts[0]
                    estimate = float(parts[1])
                    se = float(parts[2])
                    pvalue = float(parts[3])
                    
                    if component == 'Intercept':
                        results['intercept_coef'] = estimate
                        results['intercept_se'] = se
                        results['intercept_pvalue'] = pvalue
                    elif component == 'R_offspring':
                        results['R_offspring_coef'] = estimate
                        results['R_offspring_se'] = se
                        results['R_offspring_pvalue'] = pvalue
                    elif component == 'R_parental':
                        results['R_parental_coef'] = estimate
                        results['R_parental_se'] = se
                        results['R_parental_pvalue'] = pvalue
                    elif component == 'R_cross':
                        results['R_cross_coef'] = estimate
                        results['R_cross_se'] = se
                        results['R_cross_pvalue'] = pvalue
        
        # Parse model statistics
        for line in lines:
            if line.startswith('R-squared:'):
                results['r_squared'] = float(line.split(':', 1)[1].strip())
            elif line.startswith('Adjusted R-squared:'):
                results['adj_r_squared'] = float(line.split(':', 1)[1].strip())
            elif line.startswith('F-statistic:'):
                results['f_statistic'] = float(line.split(':', 1)[1].strip())
            elif line.startswith('F p-value:'):
                results['f_pvalue'] = float(line.split(':', 1)[1].strip())
            elif line.startswith('AIC:'):
                results['aic'] = float(line.split(':', 1)[1].strip())
            elif line.startswith('BIC:'):
                results['bic'] = float(line.split(':', 1)[1].strip())
            elif line.startswith('Log-likelihood:'):
                results['log_likelihood'] = float(line.split(':', 1)[1].strip())
                
        return results
    except Exception as e:
        print(f"Could not read or process file {filepath}: {e}")
        return None

def parse_filename(filepath):
    """
    Parse the filename to extract metadata
    Expected format: run_XXX_RDR_regression_N8000_Y1.RDRreg
    """
    basename = os.path.basename(filepath)
    
    try:
        # Remove .RDRreg extension
        name_parts = basename.replace('.RDRreg', '').split('_')
        
        # Extract information
        run_num = int(name_parts[1])  # run_XXX
        sample_size = int(name_parts[4].replace('N', ''))  # N8000
        trait = name_parts[5]  # Y1 or Y2
        
        return run_num, sample_size, trait
    except (IndexError, ValueError) as e:
        print(f"Warning: Could not parse filename '{basename}': {e}")
        return None, None, None

def main():
    """
    Main aggregation function
    """
    RESULTS_DIR = "/projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/RDR_Regression_Results"
    
    print(f"=== Starting RDR regression results aggregation ===")
    print(f"Results directory: {RESULTS_DIR}")
    print(f"Aggregation started at: {datetime.now()}")
    
    # Find all RDRreg result files
    rdreg_files = glob.glob(os.path.join(RESULTS_DIR, "**", "*_RDR_regression_*.RDRreg"), recursive=True)
    
    if not rdreg_files:
        print(f"No RDR regression .RDRreg files found in {RESULTS_DIR}. Exiting.")
        return

    print(f"Found {len(rdreg_files)} RDR regression result files to aggregate.")

    all_results = []
    
    for filepath in rdreg_files:
        # Parse filename
        run_num, sample_size, trait = parse_filename(filepath)
        if run_num is None:
            continue
            
        # Load results
        result_data = load_result_from_rdreg_file(filepath)
        if result_data is None:
            continue
            
        # Extract condition from directory structure
        dirname = os.path.basename(os.path.dirname(filepath))
        # The condition should be extracted from the run_id or directory structure
        # Assuming the directory name contains the condition info
        condition_parts = dirname.split('_')
        if len(condition_parts) >= 4:
            condition = '_'.join(condition_parts[:-4])  # Remove run_XXX_nfam8000 part
        else:
            condition = dirname
        
        # Add metadata to results
        result_data.update({
            'condition': condition,
            'replication': run_num,
            'sample_size': sample_size,
            'trait': trait,
            'file_path': filepath
        })
        
        all_results.append(result_data)

    # Check if results were found
    if not all_results:
        print("\nAggregation complete, but no valid RDR regression .RDRreg files could be parsed.")
        print("Please check file contents and naming conventions.")
        return

    # Create results DataFrame
    detailed_df = pd.DataFrame(all_results)
    detailed_df.sort_values(by=['condition', 'trait', 'sample_size', 'replication'], inplace=True)

    # Define columns to keep in the detailed output
    detail_columns = [
        'condition', 'replication', 'sample_size', 'trait',
        'phenotypic_variance', 'n_pairs',
        'intercept_coef', 'intercept_se', 'intercept_pvalue',
        'R_offspring_coef', 'R_offspring_se', 'R_offspring_pvalue',
        'R_parental_coef', 'R_parental_se', 'R_parental_pvalue', 
        'R_cross_coef', 'R_cross_se', 'R_cross_pvalue',
        'r_squared', 'adj_r_squared', 'f_statistic', 'f_pvalue',
        'aic', 'bic', 'log_likelihood', 'analysis_timestamp'
    ]
    
    # Keep only available columns
    available_columns = [col for col in detail_columns if col in detailed_df.columns]
    detailed_output_df = detailed_df[available_columns].copy()

    # Save detailed results
    detailed_output_path = os.path.join(RESULTS_DIR, "aggregated_rdr_regression_results_detailed.tsv")
    detailed_output_df.to_csv(detailed_output_path, sep='\t', index=False, float_format='%.6f')
    
    print(f"\nDetailed results saved to: {detailed_output_path}")

    # Create summary statistics by condition and trait
    summary_stats = []
    
    for (condition, trait), group in detailed_df.groupby(['condition', 'trait']):
        n_reps = len(group)
        
        summary = {
            'condition': condition,
            'trait': trait,
            'n_replications': n_reps,
            'sample_size': group['sample_size'].iloc[0],  # Should be constant
            # Intercept statistics
            'intercept_mean': group['intercept_coef'].mean(),
            'intercept_se_mean': group['intercept_se'].mean(),
            'intercept_sd': group['intercept_coef'].std(),
            # R_offspring statistics  
            'R_offspring_mean': group['R_offspring_coef'].mean(),
            'R_offspring_se_mean': group['R_offspring_se'].mean(),
            'R_offspring_sd': group['R_offspring_coef'].std(),
            # R_parental statistics
            'R_parental_mean': group['R_parental_coef'].mean(),
            'R_parental_se_mean': group['R_parental_se'].mean(),
            'R_parental_sd': group['R_parental_coef'].std(),
            # R_cross statistics
            'R_cross_mean': group['R_cross_coef'].mean(),
            'R_cross_se_mean': group['R_cross_se'].mean(),
            'R_cross_sd': group['R_cross_coef'].std(),
            # Model fit statistics
            'r_squared_mean': group['r_squared'].mean(),
            'r_squared_sd': group['r_squared'].std(),
            'phenotypic_variance_mean': group['phenotypic_variance'].mean(),
            'phenotypic_variance_sd': group['phenotypic_variance'].std()
        }
        
        summary_stats.append(summary)

    # Create summary DataFrame and save
    summary_df = pd.DataFrame(summary_stats)
    summary_df.sort_values(by=['condition', 'trait'], inplace=True)
    
    summary_output_path = os.path.join(RESULTS_DIR, "aggregated_rdr_regression_results_summary.tsv")
    summary_df.to_csv(summary_output_path, sep='\t', index=False, float_format='%.6f')
    
    print(f"Summary results saved to: {summary_output_path}")

    # Print final summary
    print(f"\n=== Aggregation complete ===")
    print(f"Total results aggregated: {len(detailed_df)}")
    print(f"Conditions: {sorted(detailed_df['condition'].unique())}")
    print(f"Traits: {sorted(detailed_df['trait'].unique())}")
    print(f"Sample sizes: {sorted(detailed_df['sample_size'].unique())}")
    print(f"Replications per condition: {detailed_df.groupby('condition')['replication'].nunique().to_dict()}")
    print(f"Aggregation finished at: {datetime.now()}")

if __name__ == '__main__':
    main()