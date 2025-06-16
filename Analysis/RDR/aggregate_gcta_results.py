# File: aggregate_gcta_results.py
import os
import glob
import pandas as pd

def parse_hsq_file(filepath):
    """Parses a GCTA .hsq file to extract variance estimates and SEs."""
    results = {}
    try:
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 3:
                    source, estimate, se = parts
                    # Standardize keys for easy access, e.g., 'V(G1)/Vp' -> 'h2_G1'
                    clean_key = source.replace('/', '_').replace('(', '').replace(')', '')
                    results[f"{clean_key}_est"] = float(estimate)
                    results[f"{clean_key}_se"] = float(se)
    except Exception as e:
        print(f"Could not parse file {filepath}: {e}")
    return results

def main():
    # Directory where your final .hsq files are saved
    RESULTS_DIR = "/projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/RDR_Results/"
    
    all_results = []
    
    # Find all .hsq files recursively
    hsq_files = glob.glob(os.path.join(RESULTS_DIR, "**", "*.hsq"), recursive=True)
    
    if not hsq_files:
        print(f"No .hsq files found in {RESULTS_DIR}. Exiting.")
        return

    print(f"Found {len(hsq_files)} result files to aggregate.")

    for f in hsq_files:
        basename = os.path.basename(f)
        parts = basename.split('_')
        
        # Extract info from filename, e.g., "run_001_N2000_Y1_results.hsq"
        condition = os.path.basename(os.path.dirname(f))
        run_id = int(parts[1])
        sample_size = int(parts[2].replace('N',''))
        trait = parts[3]
        
        parsed_data = parse_hsq_file(f)
        
        if parsed_data:
            base_info = {
                'condition': condition,
                'replication': run_id,
                'sample_size': sample_size,
                'trait': trait
            }
            all_results.append({**base_info, **parsed_data})
            
    # Create a final DataFrame and save it
    final_df = pd.DataFrame(all_results)
    final_df.sort_values(by=['condition', 'trait', 'sample_size', 'replication'], inplace=True)
    
    output_path = os.path.join(RESULTS_DIR, "aggregated_rdr_gcta_results.tsv")
    final_df.to_csv(output_path, sep='\t', index=False)
    
    print(f"\nAggregation complete. Final results table saved to:\n{output_path}")

if __name__ == '__main__':
    main()