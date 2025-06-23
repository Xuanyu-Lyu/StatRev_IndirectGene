# File: extract_and_average_results.py

import os
import re
import glob
import numpy as np
import pandas as pd

def parse_matrix_from_summary(file_path, generation, matrix_name):
    """
    Parses a specific 2x2 matrix from a given generation block in a summary file.
    
    Args:
        file_path (str): The path to the summary.txt file.
        generation (int): The generation number to look for (e.g., 20).
        matrix_name (str): The name of the matrix to extract (e.g., "VP", "w").

    Returns:
        A 2x2 numpy array, or None if the matrix is not found.
    """
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        return None

    in_correct_generation_block = False
    matrix_lines = []
    
    for i, line in enumerate(lines):
        # Check if we are entering the correct generation block
        if f"--- Summary for Final Generation (Generation {generation}) ---" in line or \
           f"--- Generation {generation} Summary ---" in line:
            in_correct_generation_block = True
            continue

        # If we are in the correct block, look for the matrix name
        if in_correct_generation_block:
            # Check if we've entered a new block, if so, stop
            if "---" in line and not (f"Generation {generation}" in line):
                in_correct_generation_block = False
                continue

            if line.strip().startswith(matrix_name + ":"):
                # The next two lines should contain the matrix
                try:
                    line1 = lines[i + 1].strip()
                    line2 = lines[i + 2].strip()
                    
                    # Use regex to find all floating point numbers in the lines
                    vals1 = re.findall(r"[-+]?\d*\.\d+|\d+", line1)
                    vals2 = re.findall(r"[-+]?\d*\.\d+|\d+", line2)
                    
                    if len(vals1) == 2 and len(vals2) == 2:
                        matrix = np.array([
                            [float(vals1[0]), float(vals1[1])],
                            [float(vals2[0]), float(vals2[1])]
                        ])
                        return matrix
                except (IndexError, ValueError):
                    # Could not parse the matrix lines, continue searching
                    continue
    return None


def main():
    """
    Main function to find summary files, extract matrices, average them,
    and save the final results.
    """
    # --- 1. CONFIGURATION ---
    # Top-level directory where your simulation condition folders are
    SOURCE_DIR = "/scratch/alpine/xuly4739/StatRev_IndirectGene/Data/ASHG_Final"
    
    # Where to save the final averaged output
    DESTINATION_DIR = "/projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/Averaged_Results"
    os.makedirs(DESTINATION_DIR, exist_ok=True)
    
    # List the conditions you want to process
    CONDITIONS_TO_PROCESS = ["phenoVT_phenoAM", "socialVT_phenoAM", "phenoVT_socialAM", "phenoVT_geneticAM"]
    
    # Configuration for the extraction
    GENERATION_TO_EXTRACT = 20
    MATRICES_TO_EXTRACT = ["VP", "w"]
    
    print("--- Starting Matrix Extraction and Averaging ---")

    # --- 2. Loop through each condition ---
    for condition in CONDITIONS_TO_PROCESS:
        print(f"\nProcessing condition: {condition}...")
        
        condition_path = os.path.join(SOURCE_DIR, condition)
        if not os.path.isdir(condition_path):
            print(f"  Warning: Directory not found for condition '{condition}'. Skipping.")
            continue
            
        # Find all summary files for this condition
        summary_files = glob.glob(os.path.join(condition_path, "run_*", "*_summary.txt"))
        
        if not summary_files:
            print(f"  Warning: No summary files found for condition '{condition}'. Skipping.")
            continue
        
        print(f"  Found {len(summary_files)} replication summaries to process.")
        
        # Initialize a dictionary to hold lists of matrices
        matrix_collections = {key: [] for key in MATRICES_TO_EXTRACT}
        
        # --- 3. Parse each file and collect matrices ---
        for f_path in summary_files:
            for matrix_name in MATRICES_TO_EXTRACT:
                matrix = parse_matrix_from_summary(f_path, GENERATION_TO_EXTRACT, matrix_name)
                if matrix is not None:
                    matrix_collections[matrix_name].append(matrix)

        # --- 4. Calculate averages and save to a single file ---
        output_filename = os.path.join(DESTINATION_DIR, f"Averaged_Matrices_{condition}_Gen{GENERATION_TO_EXTRACT}.txt")
        
        with open(output_filename, 'w') as f:
            f.write(f"--- Averaged Results for Condition: {condition} ---\n")
            f.write(f"Generation Analyzed: {GENERATION_TO_EXTRACT}\n")
            f.write(f"Total Replications Found: {len(summary_files)}\n")
            f.write("-" * 50 + "\n\n")

            for matrix_name, matrix_list in matrix_collections.items():
                if not matrix_list:
                    f.write(f"Matrix: {matrix_name}\n")
                    f.write("  -> No data found or could not be parsed.\n\n")
                    continue
                
                # Calculate the element-wise average of all matrices in the list
                avg_matrix = np.mean(matrix_list, axis=0)
                
                f.write(f"Average '{matrix_name}' Matrix (from {len(matrix_list)} replications):\n")
                f.write(np.array2string(avg_matrix, precision=4, separator=', ', floatmode='fixed'))
                f.write("\n\n")
        
        print(f"  -> Finished. Averaged results saved to: {output_filename}")

    print("\n--- All conditions processed. ---")


if __name__ == '__main__':
    main()