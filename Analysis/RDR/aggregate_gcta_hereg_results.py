import os
import glob
import pandas as pd
import re

def parse_hereg_file(filepath):
    """
    Parses a GCTA .HEreg file to extract variance estimates from HE-CP section only.
    Returns only estimates (no SEs or p-values).
    """
    results = {}
    
    try:
        with open(filepath, 'r') as f:
            content = f.read()
            
            # Split content into sections
            sections = content.split('\n\n')
            
            for section in sections:
                lines = section.strip().split('\n')
                if len(lines) < 2:
                    continue
                    
                # Only process HE-CP section
                if lines[0].strip() != 'HE-CP':
                    continue
                
                # Skip header line
                data_lines = lines[2:]  # Skip section name and column headers
                
                for line in data_lines:
                    line = line.strip()
                    if not line:
                        continue
                        
                    # Split by whitespace and extract values
                    parts = line.split()
                    if len(parts) >= 2:  # At least coefficient name and estimate
                        coeff_name = parts[0]
                        estimate = float(parts[1])
                        
                        # Clean coefficient name for column naming
                        clean_name = coeff_name.replace('/', '_').replace('(', '').replace(')', '').replace(' ', '_')
                        
                        # Store only estimates
                        results[f"{clean_name}_est"] = estimate
                        
    except Exception as e:
        print(f"Could not read or process file {filepath}: {e}")
    
    return results

def main():
    RESULTS_DIR = "/projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/RDR_HEreg_Results"
    all_results = []
    
    hereg_files = glob.glob(os.path.join(RESULTS_DIR, "**", "*.HEreg"), recursive=True)
    
    if not hereg_files:
        print(f"No .HEreg files found in {RESULTS_DIR}. Exiting.")
        return

    print(f"Found {len(hereg_files)} HEreg result files to aggregate.")

    for f in hereg_files:
        try:
            basename = os.path.basename(f)
            dirname = os.path.basename(os.path.dirname(f))
            
            # Parse filename: run_001_HEreg_N8000_Y1.HEreg
            parts = basename.replace('.HEreg', '').split('_')
            
            condition = dirname
            run_id = int(parts[1])
            # Skip 'HEreg' part at index 2
            sample_size = int(parts[3].replace('N',''))
            trait = parts[4]
            
            parsed_data = parse_hereg_file(f)
            
            if parsed_data:
                # Create base info
                base_info = {
                    'condition': condition,
                    'replication': run_id,
                    'sample_size': sample_size,
                    'trait': trait
                }
                
                # Add to results
                all_results.append({**base_info, **parsed_data})
                
        except (IndexError, ValueError) as e:
            print(f"Warning: Could not parse filename metadata for '{f}'. Error: {e}. Skipping.")
            continue
            
    # Check if results were found
    if not all_results:
        print("\nAggregation complete, but no valid GCTA HEreg result files could be parsed.")
        print("Please check file contents and naming conventions.")
        return

    # Create results DataFrame
    detailed_df = pd.DataFrame(all_results)
    detailed_df.sort_values(by=['condition', 'trait', 'sample_size', 'replication'], inplace=True)

    # Save results
    detailed_output_path = os.path.join(RESULTS_DIR, "aggregated_rdr_hereg_results_detailed.tsv")
    
    detailed_df.to_csv(detailed_output_path, sep='\t', index=False, float_format='%.6f')
    
    print(f"\nAggregation complete. {len(detailed_df)} results aggregated.")
    print(f"Results saved to:\n{detailed_output_path}")
    
    # Print summary of what was extracted
    print(f"\nSummary:")
    print(f"- Conditions: {sorted(detailed_df['condition'].unique())}")
    print(f"- Traits: {sorted(detailed_df['trait'].unique())}")
    print(f"- Sample sizes: {sorted(detailed_df['sample_size'].unique())}")
    print(f"- Replications per condition: {detailed_df.groupby('condition')['replication'].nunique().to_dict()}")

if __name__ == '__main__':
    main()