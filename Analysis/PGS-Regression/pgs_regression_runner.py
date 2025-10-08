#!/usr/bin/env python3
"""
PGS Regression Analysis Runner

This script runs specific PGS regression analyses on conditions 5-8:
- Kong's approach: NTm only, then NTm + NTp 
- Full PGS approach: PGSo only, then PGSo + PGSm, then PGSo + PGSm + PGSp
"""

import os
import pandas as pd
from pgs_regression_functions import run_pgs_analysis, create_analysis_summary


def run_kong_analysis(condition, trait):
    """
    Run Kong's approach analyses:
    1. NTm only
    2. NTm + NTp
    """
    print(f"  Running Kong's analyses for {trait} on {condition}")
    
    results = {}
    
    # 1. NTm only (maternal non-transmitted)
    print(f"    - NTm only analysis...")
    _, ntm_results = run_pgs_analysis(condition, trait, 'kong_maternal', save_results=True)
    results['ntm_only'] = ntm_results
    
    # 2. NTm + NTp (both non-transmitted)
    print(f"    - NTm + NTp analysis...")
    _, kong_full_results = run_pgs_analysis(condition, trait, 'kong', save_results=True)
    results['kong_full'] = kong_full_results
    
    return results


def run_full_pgs_analysis(condition, trait):
    """
    Run Full PGS approach with incremental R²:
    1. PGSo only
    2. PGSo + PGSm (incremental R² for maternal)
    3. PGSo + PGSm + PGSp (incremental R² for all three)
    """
    print(f"  Running Full PGS analyses for {trait} on {condition}")
    
    results = {}
    
    # For full PGS incremental analysis, we need the standard approach
    # which will give us incremental R² for each predictor
    print(f"    - Full PGS incremental analysis...")
    _, full_pgs_results = run_pgs_analysis(condition, trait, 'full_pgs', save_results=True)
    results['full_pgs_incremental'] = full_pgs_results
    
    return results


def extract_key_metrics(results_dict, condition, trait, analysis_type):
    """Extract key R² metrics from results"""
    
    extracted_results = []
    
    if analysis_type == 'kong':
        # Extract R² for NTm only
        if 'ntm_only' in results_dict:
            ntm_df = results_dict['ntm_only']
            trait_num = '1' if trait == 'trait1' else '2'
            r2_col = f'r2_score'  # For single predictor, should be r2_score
            
            if not ntm_df.empty:
                for _, row in ntm_df.iterrows():
                    extracted_results.append({
                        'condition': condition,
                        'trait': trait,
                        'analysis_type': 'kong_maternal_only',
                        'predictor': f'NTm{trait_num}',
                        'r2_value': row.get(r2_col, row.get(f'total_r2_NTm{trait_num}', None)),
                        'n_samples': row.get('n_samples', None)
                    })
        
        # Extract R² for NTm + NTp
        if 'kong_full' in results_dict:
            kong_df = results_dict['kong_full']
            trait_num = '1' if trait == 'trait1' else '2'
            
            if not kong_df.empty:
                # Look for incremental R² columns
                ntm_r2_col = f'total_r2_NTm{trait_num}'
                ntp_r2_col = f'total_r2_NTp{trait_num}'
                ntm_incr_col = f'incremental_r2_NTm{trait_num}'
                ntp_incr_col = f'incremental_r2_NTp{trait_num}'
                
                for _, row in kong_df.iterrows():
                    # Total R² for NTm (first predictor)
                    if ntm_r2_col in row and pd.notna(row[ntm_r2_col]):
                        extracted_results.append({
                            'condition': condition,
                            'trait': trait,
                            'analysis_type': 'kong_ntm_total',
                            'predictor': f'NTm{trait_num}',
                            'r2_value': row[ntm_r2_col],
                            'n_samples': row.get('n_samples', None)
                        })
                    
                    # Total R² for NTp (second predictor)
                    if ntp_r2_col in row and pd.notna(row[ntp_r2_col]):
                        extracted_results.append({
                            'condition': condition,
                            'trait': trait,
                            'analysis_type': 'kong_ntp_total',
                            'predictor': f'NTp{trait_num}',
                            'r2_value': row[ntp_r2_col],
                            'n_samples': row.get('n_samples', None)
                        })
                    
                    # Incremental R² for NTp
                    if ntp_incr_col in row and pd.notna(row[ntp_incr_col]):
                        extracted_results.append({
                            'condition': condition,
                            'trait': trait,
                            'analysis_type': 'kong_ntp_incremental',
                            'predictor': f'NTp{trait_num}',
                            'r2_value': row[ntp_incr_col],
                            'n_samples': row.get('n_samples', None)
                        })
    
    elif analysis_type == 'full_pgs':
        # Extract incremental R² for full PGS
        if 'full_pgs_incremental' in results_dict:
            pgs_df = results_dict['full_pgs_incremental']
            trait_num = '1' if trait == 'trait1' else '2'
            
            if not pgs_df.empty:
                # Look for incremental R² columns
                pgso_r2_col = f'total_r2_PGSo{trait_num}'
                pgsm_r2_col = f'total_r2_PGSm{trait_num}'
                pgsp_r2_col = f'total_r2_PGSp{trait_num}'
                pgso_incr_col = f'incremental_r2_PGSo{trait_num}'
                pgsm_incr_col = f'incremental_r2_PGSm{trait_num}'
                pgsp_incr_col = f'incremental_r2_PGSp{trait_num}'
                
                for _, row in pgs_df.iterrows():
                    # Total R² for PGSo (offspring only)
                    if pgso_r2_col in row and pd.notna(row[pgso_r2_col]):
                        extracted_results.append({
                            'condition': condition,
                            'trait': trait,
                            'analysis_type': 'full_pgs_offspring_total',
                            'predictor': f'PGSo{trait_num}',
                            'r2_value': row[pgso_r2_col],
                            'n_samples': row.get('n_samples', None)
                        })
                    
                    # Incremental R² for PGSm (mother)
                    if pgsm_incr_col in row and pd.notna(row[pgsm_incr_col]):
                        extracted_results.append({
                            'condition': condition,
                            'trait': trait,
                            'analysis_type': 'full_pgs_mother_incremental',
                            'predictor': f'PGSm{trait_num}',
                            'r2_value': row[pgsm_incr_col],
                            'n_samples': row.get('n_samples', None)
                        })
                    
                    # Incremental R² for PGSp (father)
                    if pgsp_incr_col in row and pd.notna(row[pgsp_incr_col]):
                        extracted_results.append({
                            'condition': condition,
                            'trait': trait,
                            'analysis_type': 'full_pgs_father_incremental',
                            'predictor': f'PGSp{trait_num}', 
                            'r2_value': row[pgsp_incr_col],
                            'n_samples': row.get('n_samples', None)
                        })
    
    return extracted_results


def main():
    """Main execution function"""
    
    print("="*80)
    print("PGS REGRESSION ANALYSIS")
    print("Running conditions 5-8 with Kong's and Full PGS approaches")
    print("="*80)
    
    # Define conditions 5-8
    conditions = [
        '05_t1pheVTnoAM_t2socVTnoAM_PGSall',
        '06_t1noVTpheAM_t2noVTnoAM_PGSall', 
        '07_t1noVTsocAM_t2noVTnoAM_PGSall',
        '08_t1noVTgenAM_t2noVTnoAM_PGSall'
    ]
    
    # Define traits to analyze
    traits = ['trait1', 'trait2']
    
    # Store all extracted results
    all_extracted_results = []
    
    # Process each condition
    for condition in conditions:
        print(f"\n{'='*60}")
        print(f"PROCESSING CONDITION: {condition}")
        print(f"{'='*60}")
        
        for trait in traits:
            print(f"\nAnalyzing {trait}...")
            
            # Run Kong's analyses
            kong_results = run_kong_analysis(condition, trait)
            kong_extracted = extract_key_metrics(kong_results, condition, trait, 'kong')
            all_extracted_results.extend(kong_extracted)
            
            # Run Full PGS analyses
            pgs_results = run_full_pgs_analysis(condition, trait)
            pgs_extracted = extract_key_metrics(pgs_results, condition, trait, 'full_pgs')
            all_extracted_results.extend(pgs_extracted)
    
    # Create summary DataFrame
    print(f"\n{'='*60}")
    print("CREATING SUMMARY RESULTS")
    print(f"{'='*60}")
    
    if all_extracted_results:
        summary_df = pd.DataFrame(all_extracted_results)
        
        # Create output directory
        output_dir = '/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/coding/PyProject/StatRev_IndirectGene/Analysis/PGS-Regression/results'
        os.makedirs(output_dir, exist_ok=True)
        
        # Save extracted results
        output_file = os.path.join(output_dir, 'conditions_5_8_extracted_r2_results.csv')
        summary_df.to_csv(output_file, index=False)
        print(f"Extracted results saved to: {output_file}")
        
        # Print summary statistics
        print("\nSUMMARY OF EXTRACTED RESULTS:")
        print(f"Total records: {len(summary_df)}")
        print("\nBreakdown by analysis type:")
        print(summary_df['analysis_type'].value_counts())
        
        print("\nSample of results:")
        print(summary_df.head(10))
        
        # Create summary statistics by condition and analysis type
        summary_stats = summary_df.groupby(['condition', 'trait', 'analysis_type'])['r2_value'].agg([
            'count', 'mean', 'std', 'min', 'max', 'median'
        ]).round(4)
        
        stats_file = os.path.join(output_dir, 'conditions_5_8_r2_summary_statistics.csv')
        summary_stats.to_csv(stats_file)
        print(f"\nSummary statistics saved to: {stats_file}")
        
    else:
        print("No results were extracted!")
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE!")
    print(f"{'='*80}")


if __name__ == '__main__':
    main()