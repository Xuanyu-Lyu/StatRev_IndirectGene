"""
PGS Regression Analysis Functions

This module contains functions for running PGS regression analyses on simulated data,
supporting both full PGS and haplotypic PGS approaches (Kong's method).
"""

import pandas as pd
import numpy as np
import os
import glob
import re
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import warnings

warnings.filterwarnings('ignore')


def read_trait_data(file_path, focal_trait='both', combine_pgs=False):
    """
    Read trait data from a file and extract specific traits.
    
    Parameters:
    file_path (str): Path to the data file
    focal_trait (str or int): Which trait(s) to extract
                              - 'trait1' or 1: Extract only trait 1 columns (ending with '1')
                              - 'trait2' or 2: Extract only trait 2 columns (ending with '2')  
                              - 'both' or 'all': Extract both trait 1 and trait 2 columns
    combine_pgs (bool): Whether to combine haplotypic PGS scores into full PGS scores
                        - True: Combine NTp+Tp->PGSp, NTm+Tm->PGSm, Tp+Tm->PGSo
                        - False: Keep original columns separate
    
    Returns:
    pandas.DataFrame: DataFrame containing the selected trait columns (and combined PGS if requested)
    """
    
    # Read the data file
    try:
        df = pd.read_csv(file_path, sep='\t')
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None
    except Exception as e:
        print(f"Error reading file: {e}")
        return None
    
    # Get all column names
    all_columns = df.columns.tolist()
    
    # Convert focal_trait to string for consistent handling
    if focal_trait == 1:
        focal_trait = 'trait1'
    elif focal_trait == 2:
        focal_trait = 'trait2'
    
    # Select columns based on focal_trait
    if focal_trait.lower() in ['trait1', '1']:
        # Select columns ending with '1'
        selected_columns = [col for col in all_columns if col.endswith('1')]
        
    elif focal_trait.lower() in ['trait2', '2']:
        # Select columns ending with '2'
        selected_columns = [col for col in all_columns if col.endswith('2')]
        
    elif focal_trait.lower() in ['both', 'all']:
        # Select all columns (both trait 1 and trait 2)
        selected_columns = all_columns
        
    else:
        print(f"Invalid focal_trait: {focal_trait}. Use 'trait1', 'trait2', or 'both'")
        return None
    
    # Get the DataFrame with selected columns
    result_df = df[selected_columns].copy()
    
    # Combine haplotypic PGS scores if requested
    if combine_pgs:
        # Helper function to combine PGS for a specific trait
        def combine_trait_pgs(df, trait_suffix):
            # Check if the required columns exist for this trait
            nt_col = f'NT{trait_suffix}'  # Non-transmitted PGS
            t_col = f'T{trait_suffix}'    # Transmitted PGS
            
            if nt_col in df.columns and t_col in df.columns:
                # Combine NTp + Tp -> PGSp (or NTm + Tm -> PGSm)
                pgs_col = f'PGS{trait_suffix}'
                df[pgs_col] = df[nt_col] + df[t_col]
        
        # Combine PGS for paternal (p) and maternal (m) scores
        if focal_trait.lower() in ['trait1', '1', 'both', 'all']:
            combine_trait_pgs(result_df, 'p1')
            combine_trait_pgs(result_df, 'm1')
            # For offspring PGS (Tp + Tm -> PGSo)
            if 'Tp1' in result_df.columns and 'Tm1' in result_df.columns:
                result_df['PGSo1'] = result_df['Tp1'] + result_df['Tm1']
                
        if focal_trait.lower() in ['trait2', '2', 'both', 'all']:
            combine_trait_pgs(result_df, 'p2')
            combine_trait_pgs(result_df, 'm2')
            # For offspring PGS (Tp + Tm -> PGSo)
            if 'Tp2' in result_df.columns and 'Tm2' in result_df.columns:
                result_df['PGSo2'] = result_df['Tp2'] + result_df['Tm2']
    
    return result_df


def custom_regression_analysis(data, filename=None, predictors=None, outcomes=None, 
                             multiple_regression=True, incremental_r2=False, **kwargs):
    """
    Generic regression analysis function that supports both simple and multiple regression
    
    Parameters:
    data (pd.DataFrame): The data to analyze
    filename (str): Optional filename for reference  
    predictors (list): List of predictor column names
    outcomes (list): List of outcome column names
    multiple_regression (bool): If True and len(predictors)>1, run multiple regression
                               If False, run simple regression for each predictor-outcome pair
    incremental_r2 (bool): If True and multiple_regression=True, calculate incremental R²
                          Shows R² contribution of each predictor when added sequentially
    **kwargs: Additional parameters for regression
    
    Returns:
    pd.DataFrame: Results from regression analyses
    """
    
    results_list = []
    
    # Default predictors and outcomes if not specified
    if predictors is None:
        predictors = [col for col in data.columns if 'PGS' in col or 'T' in col]
    if outcomes is None:
        outcomes = [col for col in data.columns if col.startswith('Y')]
    
    # Ensure predictors and outcomes are lists
    if isinstance(predictors, str):
        predictors = [predictors]
    if isinstance(outcomes, str):
        outcomes = [outcomes]
    
    # Check if we should run multiple regression
    if multiple_regression and len(predictors) > 1:
        # Run multiple regression for each outcome
        for outcome in outcomes:
            if outcome in data.columns:
                try:
                    # Check which predictors are available
                    available_predictors = [p for p in predictors if p in data.columns]
                    
                    if len(available_predictors) == 0:
                        continue
                    
                    # Prepare data for multiple regression
                    X = data[available_predictors].values
                    y = data[outcome].values
                    
                    # Remove any NaN values
                    mask = ~(pd.isna(X).any(axis=1) | pd.isna(y))
                    X_clean = X[mask]
                    y_clean = y[mask]
                    
                    if len(X_clean) > len(available_predictors) + 5:  # Need more samples than predictors
                        
                        if incremental_r2:
                            # Calculate incremental R² by adding predictors sequentially
                            previous_r2 = 0
                            
                            for i in range(1, len(available_predictors) + 1):
                                # Use first i predictors
                                current_predictors = available_predictors[:i]
                                X_current = X_clean[:, :i]
                                
                                # Fit model with current predictors
                                model_current = LinearRegression()
                                model_current.fit(X_current, y_clean)
                                y_pred_current = model_current.predict(X_current)
                                current_r2 = r2_score(y_clean, y_pred_current)
                                
                                # Calculate incremental R²
                                incremental_r2_value = current_r2 - previous_r2
                                
                                # Create incremental result
                                incremental_result = {
                                    'filename': filename,
                                    'predictors_up_to': ', '.join(current_predictors),
                                    'added_predictor': current_predictors[-1],
                                    'outcome': outcome,
                                    'n_samples': len(X_clean),
                                    'n_predictors': i,
                                    'intercept': model_current.intercept_,
                                    'total_r2': current_r2,
                                    'incremental_r2': incremental_r2_value,
                                    'r2_change': incremental_r2_value,
                                    'regression_type': 'incremental_multiple'
                                }
                                
                                # Add coefficients for current model
                                for j, pred in enumerate(current_predictors):
                                    incremental_result[f'coef_{pred}'] = model_current.coef_[j]
                                
                                results_list.append(incremental_result)
                                previous_r2 = current_r2
                            
                        else:
                            # Standard multiple regression (all predictors at once)
                            model = LinearRegression()
                            model.fit(X_clean, y_clean)
                            y_pred = model.predict(X_clean)
                            
                            # Create result for multiple regression
                            result = {
                                'filename': filename,
                                'predictors': ', '.join(available_predictors),
                                'outcome': outcome,
                                'n_samples': len(X_clean),
                                'n_predictors': len(available_predictors),
                                'intercept': model.intercept_,
                                'r2_score': r2_score(y_clean, y_pred),
                                'regression_type': 'multiple'
                            }
                            
                            # Add individual coefficients
                            for i, predictor in enumerate(available_predictors):
                                result[f'coef_{predictor}'] = model.coef_[i]
                            
                            results_list.append(result)
                        
                except Exception as e:
                    print(f"  Error in multiple regression {available_predictors} -> {outcome}: {str(e)}")
    
    else:
        # Run simple regression for each predictor-outcome pair
        for predictor in predictors:
            for outcome in outcomes:
                if predictor in data.columns and outcome in data.columns:
                    try:
                        X = data[[predictor]].values
                        y = data[outcome].values
                        
                        # Remove any NaN values
                        mask = ~(pd.isna(X).any(axis=1) | pd.isna(y))
                        X_clean = X[mask]
                        y_clean = y[mask]
                        
                        if len(X_clean) > 10:  # Minimum sample size
                            model = LinearRegression()
                            model.fit(X_clean, y_clean)
                            y_pred = model.predict(X_clean)
                            
                            result = {
                                'filename': filename,
                                'predictors': predictor,
                                'outcome': outcome,
                                'n_samples': len(X_clean),
                                'n_predictors': 1,
                                'coefficient': model.coef_[0],
                                'intercept': model.intercept_,
                                'r2_score': r2_score(y_clean, y_pred),
                                'regression_type': 'simple'
                            }
                            results_list.append(result)
                            
                    except Exception as e:
                        print(f"  Error in simple regression {predictor} -> {outcome}: {str(e)}")
    
    return pd.DataFrame(results_list) if results_list else None


def run_analysis_on_directory(directory_path, analysis_function, 
                             focal_trait='both', combine_pgs=False, 
                             file_pattern='*.txt', save_results=True, 
                             output_dir=None, **kwargs):
    """
    Generic function to loop through files in a directory and run any analysis on each file.
    
    Parameters:
    directory_path (str): Path to the directory containing data files
    analysis_function (callable): Function to run analysis on each DataFrame
    focal_trait (str): Which trait(s) to extract ('trait1', 'trait2', or 'both')
    combine_pgs (bool): Whether to combine haplotypic PGS scores
    file_pattern (str): Pattern to match files (default: '*.txt')
    save_results (bool): Whether to save results to files
    output_dir (str): Directory to save results (if None, uses directory_path)
    **kwargs: Additional keyword arguments to pass to the analysis_function
    
    Returns:
    dict: Dictionary with filename as key and analysis results as value
    """
    
    # Get all files matching the pattern
    file_pattern_full = os.path.join(directory_path, file_pattern)
    files = glob.glob(file_pattern_full)
    
    if not files:
        print(f"No files found matching pattern: {file_pattern_full}")
        return {}
    
    print(f"Found {len(files)} files to process")
    
    # Set output directory
    if output_dir is None:
        output_dir = directory_path
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    results = {}
    failed_files = []
    
    for i, file_path in enumerate(files, 1):
        filename = os.path.basename(file_path)
        print(f"\nProcessing file {i}/{len(files)}: {filename}")
        
        try:
            # Load data using our read_trait_data function
            data = read_trait_data(file_path, focal_trait=focal_trait, combine_pgs=combine_pgs)
            
            if data is None:
                print(f"  Failed to load data from {filename}")
                failed_files.append(filename)
                continue
            
            # Run analysis
            print(f"  Running analysis on {data.shape[0]} rows, {data.shape[1]} columns")
            analysis_result = analysis_function(data, filename=filename, **kwargs)
            
            # Store results
            results[filename] = analysis_result
            
            # Save results if requested
            if save_results and analysis_result is not None:
                output_file = os.path.join(output_dir, f"analysis_results_{filename.replace('.txt', '.csv')}")
                
                # Handle different types of analysis results
                if isinstance(analysis_result, pd.DataFrame):
                    analysis_result.to_csv(output_file, index=False)
                    print(f"  Saved results to: {output_file}")
                elif isinstance(analysis_result, dict):
                    # Convert dict to DataFrame for saving
                    pd.DataFrame([analysis_result]).to_csv(output_file, index=False)
                    print(f"  Saved results to: {output_file}")
            
        except Exception as e:
            print(f"  Error processing {filename}: {str(e)}")
            failed_files.append(filename)
    
    # Summary
    print(f"\n=== Processing Summary ===")
    print(f"Total files processed: {len(files)}")
    print(f"Successfully processed: {len(results)}")
    print(f"Failed: {len(failed_files)}")
    
    if failed_files:
        print(f"Failed files: {failed_files}")
    
    return results


def extract_key_results(analysis_results, key_columns=None, filename_pattern=None, 
                       pivot_columns=None, aggregate_func='first', 
                       include_metadata=True, sort_by=None):
    """
    Extract specified key results from run_analysis_on_directory output and consolidate into a DataFrame.
    
    Parameters:
    analysis_results (dict): Output from run_analysis_on_directory function
    key_columns (list): List of column names to extract from each analysis result
    filename_pattern (str): Optional regex pattern to filter filenames
    pivot_columns (list): Optional list of columns to pivot on
    aggregate_func (str or callable): How to aggregate multiple rows per file
    include_metadata (bool): Whether to include metadata columns like filename, n_samples
    sort_by (str or list): Column(s) to sort the final DataFrame by
    
    Returns:
    pd.DataFrame: Consolidated DataFrame with rows as data files and columns as extracted metrics
    """
    
    if not analysis_results:
        print("No analysis results provided")
        return pd.DataFrame()
    
    consolidated_data = []
    
    for filename, result in analysis_results.items():
        try:
            # Handle different types of results
            if result is None:
                continue
                
            # Convert to DataFrame if it's not already
            if isinstance(result, dict):
                result_df = pd.DataFrame([result])
            elif isinstance(result, pd.DataFrame):
                result_df = result.copy()
            else:
                print(f"Warning: Unsupported result type {type(result)} for {filename}")
                continue
            
            if result_df.empty:
                continue
            
            # Extract filename information
            base_info = {'filename': filename}
            
            # Extract run number or other info from filename using pattern
            if filename_pattern:
                match = re.search(filename_pattern, filename)
                if match:
                    if match.groups():
                        base_info['run_number'] = match.group(1)
                    else:
                        base_info['match'] = match.group(0)
            
            # Determine which columns to extract
            if key_columns is None:
                # Extract all numeric columns
                numeric_cols = result_df.select_dtypes(include=[np.number]).columns.tolist()
                extract_cols = numeric_cols
            else:
                # Use specified columns that exist in the result
                extract_cols = [col for col in key_columns if col in result_df.columns]
            
            # Include metadata columns if requested
            metadata_cols = []
            if include_metadata:
                possible_metadata = ['n_samples', 'n_predictors', 'outcome', 'predictors', 
                                   'regression_type', 'method', 'added_predictor']
                metadata_cols = [col for col in possible_metadata if col in result_df.columns]
            
            all_extract_cols = extract_cols + metadata_cols
            
            if not all_extract_cols:
                print(f"Warning: No columns to extract from {filename}")
                continue
            
            # Handle pivoting if requested
            if pivot_columns and any(col in result_df.columns for col in pivot_columns):
                existing_pivot_cols = [col for col in pivot_columns if col in result_df.columns]
                
                if len(existing_pivot_cols) == 1:
                    pivot_col = existing_pivot_cols[0]
                    
                    pivoted_data = base_info.copy()
                    
                    for extract_col in extract_cols:
                        if extract_col in result_df.columns:
                            for pivot_value in result_df[pivot_col].unique():
                                mask = result_df[pivot_col] == pivot_value
                                subset = result_df[mask]
                                
                                if not subset.empty:
                                    # Aggregate if multiple rows
                                    if len(subset) > 1:
                                        if aggregate_func == 'first':
                                            value = subset[extract_col].iloc[0]
                                        elif aggregate_func == 'mean':
                                            value = subset[extract_col].mean()
                                        else:
                                            value = subset[extract_col].iloc[0]
                                    else:
                                        value = subset[extract_col].iloc[0]
                                    
                                    # Create column name
                                    col_name = f"{extract_col}_{pivot_value}"
                                    pivoted_data[col_name] = value
                    
                    # Add metadata
                    for meta_col in metadata_cols:
                        if meta_col not in pivot_columns and meta_col in result_df.columns:
                            pivoted_data[meta_col] = result_df[meta_col].iloc[0]
                    
                    consolidated_data.append(pivoted_data)
            
            else:
                # No pivoting - handle multiple rows by aggregating
                if len(result_df) > 1:
                    row_data = base_info.copy()
                    
                    # Aggregate numeric columns
                    for col in extract_cols:
                        if col in result_df.columns:
                            if aggregate_func == 'first':
                                row_data[col] = result_df[col].iloc[0]
                            elif aggregate_func == 'mean':
                                row_data[col] = result_df[col].mean()
                            else:
                                row_data[col] = result_df[col].iloc[0]
                    
                    # Handle metadata columns
                    for col in metadata_cols:
                        if col in result_df.columns:
                            row_data[col] = result_df[col].iloc[0]
                    
                    consolidated_data.append(row_data)
                
                else:
                    # Single row - just extract the values
                    row_data = base_info.copy()
                    
                    for col in all_extract_cols:
                        if col in result_df.columns:
                            row_data[col] = result_df[col].iloc[0]
                    
                    consolidated_data.append(row_data)
        
        except Exception as e:
            print(f"Error processing {filename}: {str(e)}")
            continue
    
    # Create final DataFrame
    if not consolidated_data:
        print("No data was successfully extracted")
        return pd.DataFrame()
    
    final_df = pd.DataFrame(consolidated_data)
    
    # Sort if requested
    if sort_by and sort_by in final_df.columns:
        final_df = final_df.sort_values(sort_by).reset_index(drop=True)
    elif isinstance(sort_by, list):
        available_sort_cols = [col for col in sort_by if col in final_df.columns]
        if available_sort_cols:
            final_df = final_df.sort_values(available_sort_cols).reset_index(drop=True)
    
    return final_df


def extract_r2_results(analysis_results, r2_types=None, filename_pattern=r'run_(\d+)', 
                      include_predictors=True, sort_by_run=True):
    """
    Specialized function to extract R² values from regression analysis results.
    
    Parameters:
    analysis_results (dict): Output from run_analysis_on_directory function
    r2_types (list): List of R² column names to extract
    filename_pattern (str): Regex pattern to extract run numbers from filenames
    include_predictors (bool): Whether to include predictor information
    sort_by_run (bool): Whether to sort by run number
    
    Returns:
    pd.DataFrame: DataFrame with runs as rows and R² values as columns
    """
    
    # Auto-detect R² column types if not specified
    if r2_types is None:
        r2_types = []
        # Check a sample of results to find R² columns
        for result in analysis_results.values():
            if result is not None:
                if isinstance(result, dict):
                    sample_df = pd.DataFrame([result])
                elif isinstance(result, pd.DataFrame):
                    sample_df = result
                else:
                    continue
                
                # Find columns that likely contain R² values
                potential_r2_cols = [col for col in sample_df.columns 
                                   if 'r2' in col.lower() or 'r_squared' in col.lower()]
                r2_types.extend(potential_r2_cols)
        
        # Remove duplicates and sort
        r2_types = sorted(list(set(r2_types)))
        print(f"Auto-detected R² columns: {r2_types}")
    
    # Determine what to include based on the type of analysis
    key_columns = r2_types.copy()
    
    # Add predictor information if requested
    if include_predictors:
        key_columns.extend(['predictors', 'added_predictor', 'outcome'])
    
    # Check if we need to pivot (for incremental R² results)
    needs_pivot = any('incremental' in col.lower() for col in r2_types)
    pivot_columns = ['added_predictor'] if needs_pivot else None
    
    # Extract the results
    r2_df = extract_key_results(
        analysis_results=analysis_results,
        key_columns=key_columns,
        filename_pattern=filename_pattern,
        pivot_columns=pivot_columns,
        aggregate_func='first',
        include_metadata=True,
        sort_by='run_number' if sort_by_run else None
    )
    
    return r2_df


def run_pgs_analysis(condition, trait, analysis_type='full_pgs', data_base_dir=None,
                    output_suffix='', incremental_r2=True, save_results=True):
    """
    Unified function to run PGS analysis for both full PGS and Kong's haplotypic approaches.
    
    Parameters:
    condition (str): Data condition (e.g., 't1pheVT_t2socVT_uniphenoAM')
    trait (str): Trait to analyze ('trait1' or 'trait2') 
    analysis_type (str): Type of analysis
                        - 'full_pgs': Full PGS approach (PGSo, PGSp, PGSm)
                        - 'kong': Kong's approach (NTm, NTp)
                        - 'kong_maternal': Kong's approach with maternal only (NTm)
                        - 'kong_paternal': Kong's approach with paternal only (NTp)
    data_base_dir (str): Base directory containing data (if None, uses default)
    output_suffix (str): Suffix to add to output filenames
    incremental_r2 (bool): Whether to calculate incremental R²
    save_results (bool): Whether to save results to CSV
    
    Returns:
    tuple: (raw_results_dict, summary_dataframe)
    """
    
    print(f"=== Running {analysis_type} analysis for {trait} on {condition} ===")
    
    # Set up data directory
    if data_base_dir is None:
        data_base_dir = "/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/coding/PyProject/StatRev_IndirectGene/Data"
    
    data_directory = f"{data_base_dir}/{condition}/nfam8000"
    
    # Define predictors and outcomes based on analysis type and trait
    trait_num = '1' if trait == 'trait1' else '2'
    outcome = f'Yo{trait_num}'
    
    if analysis_type == 'full_pgs':
        predictors = [f'PGSo{trait_num}', f'PGSp{trait_num}', f'PGSm{trait_num}']
        combine_pgs = True
        print(f"  Using full PGS predictors: {predictors}")
    elif analysis_type == 'kong':
        predictors = [f'NTm{trait_num}', f'NTp{trait_num}']
        combine_pgs = False
        print(f"  Using Kong's approach predictors: {predictors}")
    elif analysis_type == 'kong_maternal':
        predictors = [f'NTm{trait_num}']
        combine_pgs = False
        print(f"  Using Kong's maternal predictor: {predictors}")
    elif analysis_type == 'kong_paternal':
        predictors = [f'NTp{trait_num}']
        combine_pgs = False
        print(f"  Using Kong's paternal predictor: {predictors}")
    else:
        raise ValueError(f"Unknown analysis_type: {analysis_type}")
    
    # Run analysis on all files
    files = glob.glob(os.path.join(data_directory, "*.txt"))
    print(f"  Processing {len(files)} files...")
    
    demo_results = {}
    for file_path in files:
        filename = os.path.basename(file_path)
        
        try:
            data = read_trait_data(file_path, focal_trait='both', combine_pgs=combine_pgs)
            if data is not None:
                result = custom_regression_analysis(
                    data, 
                    filename=filename,
                    predictors=predictors,
                    outcomes=[outcome],
                    multiple_regression=True,
                    incremental_r2=incremental_r2
                )
                demo_results[filename] = result
        except Exception as e:
            print(f"    Error processing {filename}: {str(e)}")
    
    print(f"  Successfully processed {len(demo_results)} files")
    
    # Extract R² results
    print("  Extracting R² values...")
    r2_summary = extract_r2_results(
        demo_results,
        r2_types=['r2_score', 'total_r2', 'incremental_r2'],
        filename_pattern=r'run_(\d+)',
        include_predictors=True,
        sort_by_run=True
    )
    
    # Save results if requested
    if save_results:
        output_dir = os.path.join(os.path.dirname(data_base_dir), 'Analysis', 'PGS-Regression', 'results')
        os.makedirs(output_dir, exist_ok=True)
        output_filename = f"{condition}_{analysis_type}_{trait}_r2_summary{output_suffix}.csv"
        output_path = os.path.join(output_dir, output_filename)
        r2_summary.to_csv(output_path, index=False)
        print(f"  Saved results to: {output_filename}")
    
    return demo_results, r2_summary


def create_analysis_summary(results_dict, condition_list, analysis_type, trait, 
                           output_filename, target_r2_columns=None):
    """
    Create summary statistics from multiple analysis results.
    
    Parameters:
    results_dict (dict): Dictionary with condition as key and summary_df as value
    condition_list (list): List of condition names
    analysis_type (str): Type of analysis for labeling
    trait (str): Trait analyzed
    output_filename (str): Name of output TSV file
    target_r2_columns (list): Specific R² columns to summarize (if None, auto-detect)
    
    Returns:
    pd.DataFrame: Summary statistics
    """
    
    # Combine all results
    combined_dfs = []
    for condition in condition_list:
        if condition in results_dict:
            df = results_dict[condition].copy()
            df['condition'] = condition
            combined_dfs.append(df)
    
    if not combined_dfs:
        print("No results to summarize")
        return pd.DataFrame()
    
    df_combined = pd.concat(combined_dfs, ignore_index=True)
    
    # Auto-detect R² columns if not specified
    if target_r2_columns is None:
        target_r2_columns = [col for col in df_combined.columns 
                           if any(r2_name in col.lower() for r2_name in ['r2_score', 'total_r2', 'incremental_r2'])]
    
    # Find existing columns
    existing_columns = [col for col in target_r2_columns if col in df_combined.columns]
    print(f"Creating summary for columns: {existing_columns}")
    
    # Create summary statistics
    summary_list = []
    for condition in df_combined['condition'].unique():
        condition_data = df_combined[df_combined['condition'] == condition]
        for col in existing_columns:
            values = condition_data[col].dropna()
            if len(values) > 0:
                summary_list.append({
                    'condition': condition,
                    'analysis_type': analysis_type,
                    'trait': trait,
                    'variable': col,
                    'count': len(values),
                    'mean': values.mean(),
                    'std': values.std(),
                    'min': values.min(),
                    'max': values.max(),
                    'median': values.median(),
                    'mad': np.median(np.abs(values - values.median()))
                })
    
    summary_stats = pd.DataFrame(summary_list).round(4)
    
    # Save summary statistics
    output_path = f'/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/coding/PyProject/StatRev_IndirectGene/Analysis/PGS-Regression/results/{output_filename}'
    summary_stats.to_csv(output_path, sep='\t', index=False)
    print(f"Saved summary statistics to: {output_filename}")
    
    return summary_stats