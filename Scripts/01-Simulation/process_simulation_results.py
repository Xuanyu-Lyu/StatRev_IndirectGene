# File: process_simulation_results.py

import os
import sys
import glob
import multiprocessing
import pandas as pd
import numpy as np
import re

def _ensure_folder_exists(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

# ... (other helper functions like _save_dataframe_to_tsv etc. can remain the same) ...


def process_single_run(task_params):
    """
    Worker function to process one simulation run folder for one target sample size.
    Uses regular expressions for precise matching of phenotype files.
    """
    condition_name = task_params['condition_name']
    run_folder_path = task_params['run_folder_path']
    dest_folder = task_params['dest_folder']
    target_n = task_params['target_n']
    columns_to_extract = task_params['columns_to_extract']
    
    run_id = os.path.basename(run_folder_path)
    output_filename = f"{condition_name}_{run_id}_nfam{int(target_n)}.txt"
    output_filepath = os.path.join(dest_folder, output_filename)

    try:
        if os.path.exists(output_filepath):
            return f"Skipped: {output_filename} (already exists)"

        print(f"  -> Worker starting on: {run_id} for n={target_n}")

        # Use a broad glob to get potential files first
        all_possible_files = glob.glob(os.path.join(run_folder_path, "*_phen_gen*.tsv"))
        if not all_possible_files:
            return f"Failed: No potentially matching phenotype files in {run_folder_path}"

        # --- START: ROBUST PARSING WITH REGEX ---
        valid_files = {} # Dictionary to map generation_number -> file_path
        # This pattern specifically looks for "_gen", followed by one or more digits,
        # right before the ".tsv" at the end of the filename.
        pattern = re.compile(r'_gen(\d+)\.tsv$')

        for f in all_possible_files:
            match = pattern.search(f)
            # If the filename matches the pattern...
            if match:
                # ...extract the number (the part in parentheses) and store it.
                gen_num = int(match.group(1))
                valid_files[gen_num] = f
        
        if not valid_files:
            return (f"Failed: No files matching the required '_gen<NUMBER>.tsv' format "
                    f"found in {run_folder_path}")
        
        # Find the highest generation number and its corresponding file path
        final_gen_num = max(valid_files.keys())
        final_gen_filepath = valid_files[final_gen_num]
        # --- END: ROBUST PARSING WITH REGEX ---

        print(f"     Reading file: {os.path.basename(final_gen_filepath)}")
        phen_df = pd.read_csv(final_gen_filepath, sep='\t')
        
        print(f"     Read {len(phen_df)} rows. Processing...")
        
        cols_to_load = ['ID', 'Father.ID', 'Mother.ID'] + list(columns_to_extract.keys())
        phen_df = phen_df[cols_to_load]
        phen_df.drop_duplicates(subset=['Father.ID'], inplace=True)
        
        num_families = len(phen_df)
        if num_families < target_n:
            print(f"     Warning: Only {num_families} families available, less than target_n={target_n}.")
            sampled_df = phen_df.copy()
        else:
            sampled_df = phen_df.sample(n=int(target_n), random_state=42)

        sampled_df.rename(columns=columns_to_extract, inplace=True)
        final_df = sampled_df[list(columns_to_extract.values())]
        
        _ensure_folder_exists(dest_folder)
        
        print(f"     Processing done. Writing to: {output_filepath}")
        final_df.to_csv(output_filepath, sep='\t', index=False, header=True)
        
        return f"Success: {output_filename}"

    except Exception as e:
        print(f"!!! ERROR processing {run_folder_path} for n={target_n}: {e} !!!")
        import traceback
        traceback.print_exc()
        return f"Failed: {run_folder_path}"


# --- MODIFIED MAIN FUNCTION ---
def main():
    """
    Main function to define conditions and orchestrate parallel processing
    using a more robust method to catch errors from worker processes.
    """
    SOURCE_DATA_DIR = "/scratch/alpine/xuly4739/StatRev_IndirectGene/Data/Paper" # Make sure this points to the right batch folder
    DESTINATION_DIR = "/projects/xuly4739/Py_Projects/StatRev_IndirectGene/Data/Paper"
    #CONDITIONS_TO_PROCESS = [#"phenoVT_phenoAM",
                             #"socialVT_phenoAM",
                             #"phenoVT_socialAM", 
                             #"phenoVT_geneticAM", 
                             #"socialphenoVT_phenoAM",
                             #"t1pheVT_t2socVT_uniphenoAM",
                             #"01_t1pheVTnoAM_t2socVTnoAM"] # Add or remove conditions as needed
    # CONDITIONS_TO_PROCESS = ["01_t1pheVTnoAM_t2socVTnoAM", "02_t1noVTpheAM_t2noVTnoAM"
    #                          #, "03_t1noVTsocAM_t2noVTnoAM", "04_t1noVTgenAM_t2noVTnoAM"
    #                          ]  
    CONDITIONS_TO_PROCESS = [#"05_t1pheVTnoAM_t2socVTnoAM_PGSall", 
                             #"06_t1noVTpheAM_t2pheVTpheAM_PGSall",
                             "07_t1noVTsocAM_t2pheVTsocAM_PGSall", 
                             "08_t1noVTgenAM_t2pheVTgenAM_PGSall"]     
    TARGET_SAMPLE_SIZES = [2000, 4000, 8000, 16000, 32000]
    NUM_PROCESSES = int(os.environ.get('SLURM_CPUS_PER_TASK', 10))
    DEFAULT_COLUMNS_TO_EXTRACT = {
        'Y1P': 'Yp1', 'Y2P': 'Yp2', 'Y1M': 'Ym1', 'Y2M': 'Ym2', 'Y1': 'Yo1', 'Y2': 'Yo2',
        'TPO1':'Tp1', 'TPO2':'Tp2', 'NTPO1':'NTp1','NTPO2':'NTp2',
        'TMO1':'Tm1', 'TMO2':'Tm2', 'NTMO1':'NTm1','NTMO2':'NTm2'
    }

    # --- 3. Generate Task List ---
    tasks = []
    print("--- Preparing Processing Tasks ---")
    for condition in CONDITIONS_TO_PROCESS:
        condition_source_path = os.path.join(SOURCE_DATA_DIR, condition)
        if not os.path.isdir(condition_source_path):
            print(f"Warning: Source directory for condition '{condition}' not found. Skipping."); continue
            
        run_folders = glob.glob(os.path.join(condition_source_path, "run_*"))
        
        for run_folder in run_folders:
            for target_n in TARGET_SAMPLE_SIZES:
                dest_folder_for_task = os.path.join(DESTINATION_DIR, condition, f"nfam{int(target_n)}")
                
                # Check if output exists BEFORE adding to task list
                output_filename = f"{condition}_{os.path.basename(run_folder)}_nfam{int(target_n)}.txt"
                if os.path.exists(os.path.join(dest_folder_for_task, output_filename)):
                    continue # Skip this task entirely if output already exists

                tasks.append({
                    'condition_name': condition, 'run_folder_path': run_folder,
                    'dest_folder': dest_folder_for_task, 'target_n': target_n,
                    'columns_to_extract': DEFAULT_COLUMNS_TO_EXTRACT
                })

    if not tasks:
        print("No new tasks to run. All output files may already exist.")
        return

    print(f"Total new files to generate: {len(tasks)}")
    print(f"Using {NUM_PROCESSES} parallel processes.")

    # --- 4. Run All Processing Tasks in Parallel (MODIFIED TO CATCH ERRORS) ---
    print("\n--- Starting Parallel Processing ---")
    
    results = []
    # Use a Pool of processes
    with multiprocessing.Pool(processes=NUM_PROCESSES) as pool:
        # Submit all tasks without blocking using apply_async
        async_results = [pool.apply_async(process_single_run, args=(task,)) for task in tasks]
        
        # Iterate through the results, .get() will block until that specific task is done
        # If a worker process crashed, .get() will raise the exception here.
        for i, res in enumerate(async_results):
            try:
                # Use a timeout to detect tasks that are truly hanging
                result_message = res.get(timeout=1800) # Wait up to 30 minutes for a single file
                results.append(result_message)
                print(f"({i+1}/{len(tasks)}) Task completed: {result_message}")
            except multiprocessing.TimeoutError:
                print(f"!!! TIMEOUT ERROR: Task {i+1} ({tasks[i]['run_folder_path']}) took too long and was aborted.")
                results.append("Failed: Timeout")
            except Exception as e:
                print(f"!!! WORKER ERROR: Task {i+1} ({tasks[i]['run_folder_path']}) failed with an exception.")
                results.append(f"Failed: {e}")

    print("\n--- All Processing Tasks Finished or Timed Out ---")
    
    success_count = sum(1 for r in results if str(r).startswith("Success"))
    skipped_count = sum(1 for r in results if str(r).startswith("Skipped"))
    failure_count = len(results) - success_count - skipped_count
    print(f"\nBatch Summary: {success_count} successful, {skipped_count} skipped, {failure_count} failed.")

if __name__ == '__main__':
    main()