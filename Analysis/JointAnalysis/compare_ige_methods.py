#!/usr/bin/env python3
"""
IGE Methods Comparison Script

This script compares four methods for estimating indirect genetic effects (IGE):
1. RDR (Relatedness Disequilibrium Regression)
2. SEM-PGS (Structural Equation Modeling with PGS)
3. Full PGS Regression
4. Kong's Haplotypic PGS Regression

Usage:
    python compare_ige_methods.py --condition 05_t1pheVTnoAM_t2socVTnoAM_PGSall --trait 1
    python compare_ige_methods.py --condition 06_t1noVTpheAM_t2noVTnoAM_PGSall --trait 2
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys
from pathlib import Path


def load_rdr_data(condition):
    """Load RDR results for the specified condition."""
    try:
        # Load RDR aggregated results from the correct path
        rdr_file = "/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/coding/PyProject/StatRev_IndirectGene/Analysis/RDR/Results/aggregated_rdr_hereg_results_detailed.tsv"
        df_rdr = pd.read_csv(rdr_file, sep='\t')
        
        # Filter for the specific condition
        condition_data = df_rdr[df_rdr['condition'] == condition]
        
        if condition_data.empty:
            print(f"Warning: No RDR data found for condition {condition}")
            print(f"Available conditions: {df_rdr['condition'].unique()}")
            return None
        
        return condition_data
    except Exception as e:
        print(f"Error loading RDR data: {e}")
        return None


def load_sempgs_data(condition):
    """Load SEM-PGS method results for the specified condition."""
    try:
        # First try to load individual condition file
        individual_file = f'../BiSEMPGS/{condition}_parameters_with_effects.csv'
        try:
            df = pd.read_csv(individual_file)
            print(f"SEM-PGS: Loaded {len(df)} rows with {len(df.columns)} columns from individual file")
            return df
        except:
            # Fall back to combined file and filter by condition
            combined_file = '../BiSEMPGS/all_conditions_parameters.csv'
            df_all = pd.read_csv(combined_file)
            df = df_all[df_all['condition'] == condition] if 'condition' in df_all.columns else df_all
            print(f"SEM-PGS: Loaded {len(df)} rows with {len(df.columns)} columns from combined file")
            
            if df.empty:
                print(f"Warning: No SEM-PGS data found for condition {condition}")
                print(f"Available conditions: {df_all['condition'].unique() if 'condition' in df_all.columns else 'No condition column'}")
                return None
            
            return df
    except Exception as e:
        print(f"Error loading SEM-PGS data: {e}")
        return None


def load_pgs_regression_data(condition, trait, analysis_type):
    """Load PGS regression results for the specified condition and trait."""
    try:
        # Load PGS regression results with correct file naming
        pgs_file = f"/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/coding/PyProject/StatRev_IndirectGene/Analysis/PGS-Regression/results/{condition}_{analysis_type}_trait{trait}_r2_summary.csv"
        df_pgs = pd.read_csv(pgs_file)
        
        if df_pgs.empty:
            print(f"Warning: No {analysis_type} PGS data found for condition {condition}, trait {trait}")
            return None
        
        return df_pgs
    except Exception as e:
        print(f"Error loading {analysis_type} PGS data: {e}")
        return None


def calculate_rdr_ige_proportions(rdr_data, trait):
    """Calculate IGE proportions from RDR data - returns both VG2 and VG3 parameters."""
    try:
        # Filter for the specific trait
        trait_data = rdr_data[rdr_data['trait'] == f'Y{trait}']
        
        if trait_data.empty:
            return None, None
        
        # Return both VG2_Vp_est (second parameter) and VG3_Vp_est (third parameter)
        vg2_proportions = trait_data['VG2_Vp_est'].values
        vg3_proportions = trait_data['VG3_Vp_est'].values
        
        return vg2_proportions, vg3_proportions
    except Exception as e:
        print(f"Error calculating RDR IGE proportions: {e}")
        return None, None


def calculate_sempgs_ige_proportions(sempgs_data, trait):
    """Calculate direct IGE proportions from SEM-PGS data using f parameters."""
    try:
        trait_suffix = f'{trait}{trait}'  # e.g., '11' for trait 1, '22' for trait 2
        
        # Get genetic nurture effect (f parameter) and total phenotypic variance
        f_col = f'phi{trait_suffix}'  # f11 for trait 1, f22 for trait 2
        vy_col = f'VY{trait_suffix}'  # VY11 for trait 1, VY22 for trait 2
        
        if not all(col in sempgs_data.columns for col in [f_col, vy_col]):
            print(f"Warning: Missing required SEM-PGS columns for trait {trait}")
            print(f"Looking for: {f_col}, {vy_col}")
            print(f"Available columns: {[col for col in sempgs_data.columns if 'f' in col.lower() or 'vy' in col.lower()]}")
            return None
        
        # Calculate IGE proportions, filtering out invalid values
        f_vals = sempgs_data[f_col].values
        vy_vals = sempgs_data[vy_col].values
        
        # Filter out rows with NaN or zero VY values
        valid_mask = ~(np.isnan(f_vals) | np.isnan(vy_vals) | (vy_vals == 0))
        
        if not np.any(valid_mask):
            print(f"Warning: No valid SEM-PGS data for trait {trait}")
            return None
        
        f_vals = f_vals[valid_mask]
        vy_vals = vy_vals[valid_mask]
        
        # Direct IGE proportion: genetic nurture effect / total phenotypic variance
        # Note: f parameter represents the genetic nurture pathway strength
        direct_ige = np.abs(f_vals) / vy_vals  # Take absolute value for proportion
        
        return direct_ige
    except Exception as e:
        print(f"Error calculating SEM-PGS IGE proportions: {e}")
        return None


def calculate_fullpgs_ige_proportions(fullpgs_data, trait):
    """Calculate IGE proportions from Full PGS regression data."""
    try:
        # The IGE effect is the incremental R² from adding parental PGS (mother + father)
        # We have incremental R² for mother and father separately, so sum them
        
        mother_incr_col = f'incremental_r2_PGSm{trait}'
        father_incr_col = f'incremental_r2_PGSp{trait}'
        
        print(f"Looking for columns: {mother_incr_col}, {father_incr_col}")
        print(f"Available columns: {fullpgs_data.columns.tolist()}")
        
        if mother_incr_col not in fullpgs_data.columns or father_incr_col not in fullpgs_data.columns:
            print(f"Warning: Missing required incremental R² columns for Full PGS IGE calculation")
            return None
        
        # Calculate total IGE effect as sum of maternal and paternal incremental R²
        mother_incr = fullpgs_data[mother_incr_col].values
        father_incr = fullpgs_data[father_incr_col].values
        
        # Filter out invalid values
        valid_mask = ~(np.isnan(mother_incr) | np.isnan(father_incr))
        
        if not np.any(valid_mask):
            print(f"Warning: No valid Full PGS data for trait {trait}")
            return None
        
        mother_incr = mother_incr[valid_mask]
        father_incr = father_incr[valid_mask]
        
        # IGE proportion is the sum of incremental R² from both parents
        ige_proportions = mother_incr + father_incr
        
        return ige_proportions
    except Exception as e:
        print(f"Error calculating Full PGS IGE proportions: {e}")
        return None


def calculate_kong_ige_proportions(kong_data, trait):
    """Calculate IGE proportions from Kong's haplotypic PGS regression data."""
    try:
        # Look for incremental R² from both maternal and paternal non-transmitted alleles
        # These represent the IGE effects
        ntm_col = f'incremental_r2_NTm{trait}'
        ntp_col = f'incremental_r2_NTp{trait}'
        
        if ntm_col not in kong_data.columns or ntp_col not in kong_data.columns:
            print(f"Warning: Missing NT R² columns in Kong data")
            print(f"Available columns: {kong_data.columns.tolist()}")
            return None
        
        # The IGE effect is the sum of incremental R² from both parental non-transmitted alleles
        ntm_r2 = kong_data[ntm_col].values
        ntp_r2 = kong_data[ntp_col].values
        
        # Filter out invalid values
        valid_mask = ~(np.isnan(ntm_r2) | np.isnan(ntp_r2))
        
        if not np.any(valid_mask):
            print(f"Warning: No valid Kong PGS data for trait {trait}")
            return None
        
        ntm_r2 = ntm_r2[valid_mask]
        ntp_r2 = ntp_r2[valid_mask]
        
        # Total IGE proportion is the sum of maternal and paternal non-transmitted effects
        ige_proportions = ntm_r2 + ntp_r2
        
        return ige_proportions
    except Exception as e:
        print(f"Error calculating Kong IGE proportions: {e}")
        return None


def create_violin_plot(data_dict, condition, trait, save_path=None):
    """Create violin plot comparing the four methods."""
    
    # Prepare data for plotting
    plot_data = []
    for method, values in data_dict.items():
        if values is not None and len(values) > 0:
            for val in values:
                if not np.isnan(val):
                    plot_data.append({'Method': method, 'IGE_Proportion': val})
    
    if not plot_data:
        print("No valid data to plot")
        return
    
    df_plot = pd.DataFrame(plot_data)
    
    # Create the violin plot
    plt.figure(figsize=(12, 8))
    
    # Define colors for each method
    method_colors = {
        'RDR VG2': '#d62728',
        'RDR VG3': '#b91c1c',
        'SEM-PGS Direct': '#9467bd',
        'Full PGS': '#ff7f0e',
        'Kong PGS': '#2ca02c'
    }
    
    # Create violin plot
    ax = sns.violinplot(data=df_plot, x='Method', y='IGE_Proportion', 
                       palette=[method_colors.get(method, '#666666') for method in df_plot['Method'].unique()],
                       inner='box')
    
    # Customize the plot
    plt.title(f'Comparison of IGE Estimation Methods\nCondition: {condition}, Trait {trait}', 
              fontsize=14, fontweight='bold')
    plt.ylabel('Proportion of Variance Explained by IGE', fontsize=12)
    plt.xlabel('Method', fontsize=12)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    
    # Add grid
    plt.grid(True, alpha=0.3, axis='y')
    
    # Add summary statistics as text (using median instead of mean)
    summary_text = []
    for method in df_plot['Method'].unique():
        method_data = df_plot[df_plot['Method'] == method]['IGE_Proportion']
        median_val = method_data.median()
        std_val = method_data.std()
        n_val = len(method_data)
        summary_text.append(f'{method}: med={median_val:.4f}, σ={std_val:.4f}, n={n_val}')
    
    # Add summary text box
    textstr = '\n'.join(summary_text)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', bbox=props)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save plot if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {save_path}")
    else:
        # Show plot in notebook/interactive environment
        plt.show()


def main():
    """Main function to run the comparison."""
    parser = argparse.ArgumentParser(description='Compare IGE estimation methods')
    parser.add_argument('--condition', required=True, 
                       help='Condition to analyze (e.g., 05_t1pheVTnoAM_t2socVTnoAM_PGSall)')
    parser.add_argument('--trait', type=int, required=True, choices=[1, 2],
                       help='Trait number to analyze (1 or 2)')
    parser.add_argument('--save', 
                       help='Path to save the plot (optional)')
    
    args = parser.parse_args()
    
    condition = args.condition
    trait = args.trait
    
    print(f"Analyzing condition: {condition}, trait: {trait}")
    print("="*60)
    
    # Load data from all methods
    print("Loading data...")
    
    # 1. Load RDR data
    rdr_data = load_rdr_data(condition)
    if rdr_data is not None:
        rdr_vg2, rdr_vg3 = calculate_rdr_ige_proportions(rdr_data, trait)
    else:
        rdr_vg2, rdr_vg3 = None, None
    
    # 2. Load SEM-PGS data (direct IGE only)
    sempgs_data = load_sempgs_data(condition)
    if sempgs_data is not None:
        sempgs_direct = calculate_sempgs_ige_proportions(sempgs_data, trait)
    else:
        sempgs_direct = None
    
    # 3. Load Full PGS data
    fullpgs_data = load_pgs_regression_data(condition, trait, 'full_pgs')
    fullpgs_ige = calculate_fullpgs_ige_proportions(fullpgs_data, trait) if fullpgs_data is not None else None
    
    # 4. Load Kong PGS data
    kong_data = load_pgs_regression_data(condition, trait, 'kong')
    kong_ige = calculate_kong_ige_proportions(kong_data, trait) if kong_data is not None else None
    
    # Prepare data dictionary
    data_dict = {}
    
    if rdr_vg2 is not None:
        data_dict['RDR VG2'] = rdr_vg2
        print(f"RDR VG2: {len(rdr_vg2)} observations, median={np.median(rdr_vg2):.4f}")
    
    if rdr_vg3 is not None:
        data_dict['RDR VG3'] = rdr_vg3
        print(f"RDR VG3: {len(rdr_vg3)} observations, median={np.median(rdr_vg3):.4f}")
    
    if sempgs_direct is not None:
        data_dict['SEM-PGS Direct'] = sempgs_direct
        print(f"SEM-PGS Direct: {len(sempgs_direct)} observations, median={np.median(sempgs_direct):.4f}")
    
    if fullpgs_ige is not None:
        data_dict['Full PGS'] = fullpgs_ige
        print(f"Full PGS: {len(fullpgs_ige)} observations, median={np.median(fullpgs_ige):.4f}")
    
    if kong_ige is not None:
        data_dict['Kong PGS'] = kong_ige
        print(f"Kong PGS: {len(kong_ige)} observations, median={np.median(kong_ige):.4f}")
    
    # Create violin plot
    print("\nCreating violin plot...")
    create_violin_plot(data_dict, condition, trait, args.save)
    
    print("Analysis complete!")


if __name__ == '__main__':
    main()