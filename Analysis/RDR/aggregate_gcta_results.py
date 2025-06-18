import os
import glob
import pandas as pd
import re

def parse_hsq_file(filepath):
    """
    Parses a GCTA .hsq file to extract variance estimates and SEs.
    This version robustly handles sources with spaces and special 2-column lines.
    """
    results = {}
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                # Use rsplit to handle sources that contain spaces, like "Sum of V(G)/Vp"
                # It splits from the right at most 2 times.
                parts = line.rsplit(None, 2)
                
                try:
                    # Handle standard 3-column variance estimate lines
                    if len(parts) == 3:
                        source, estimate, se = parts
                        est_val = float(estimate)
                        se_val = float(se)
                        
                        clean_key = source.replace('/', '_').replace('(', '').replace(')', '').replace(' ', '_')
                        results[f"{clean_key}_est"] = est_val
                        results[f"{clean_key}_se"] = se_val

                    # Handle special 2-column lines like logL and n
                    elif len(parts) == 2:
                        key, value = parts
                        val = float(value)
                        
                        if key.lower() == 'logl':
                            results['logL'] = val
                        elif key.lower() == 'n':
                            results['n_from_file'] = val

                except (ValueError, IndexError):
                    # This will skip the header or any other malformed lines
                    continue
    except Exception as e:
        print(f"Could not read or process file {filepath}: {e}")
    return results

def main():
    RESULTS_DIR = "/projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/RDR_Results/Tests"
    all_results = []
    
    hsq_files = glob.glob(os.path.join(RESULTS_DIR, "**", "*.hsq"), recursive=True)
    
    if not hsq_files:
        print(f"No .hsq files found in {RESULTS_DIR}. Exiting.")
        return

    print(f"Found {len(hsq_files)} result files to aggregate.")

    for f in hsq_files:
        try:
            basename = os.path.basename(f)
            dirname = os.path.basename(os.path.dirname(f))
            
            parts = basename.replace('.hsq', '').split('_')
            
            condition = dirname
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
        except (IndexError, ValueError) as e:
            print(f"Warning: Could not parse filename metadata for '{f}'. Error: {e}. Skipping.")
            continue
            
    # CRITICAL CHECK: Ensure results were actually found before creating a DataFrame
    if not all_results:
        print("\nAggregation complete, but no valid GCTA result files could be parsed.")
        print("Please check file contents and naming conventions.")
        return

    final_df = pd.DataFrame(all_results)
    final_df.sort_values(by=['condition', 'trait', 'sample_size', 'replication'], inplace=True)

    output_path = os.path.join(RESULTS_DIR, "aggregated_rdr_gcta_results.tsv")
    final_df.to_csv(output_path, sep='\t', index=False, float_format='%.6f')
    
    print(f"\nAggregation complete. {len(final_df)} results aggregated.")
    print(f"Final results table saved to:\n{output_path}")

if __name__ == '__main__':
    main()