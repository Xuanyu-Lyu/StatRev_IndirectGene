# File: run_simulations.py

import os
import sys
import multiprocessing
import datetime
import numpy as np
import pandas as pd

# Import the necessary class and functions from your other files
try:
    from SimulationFunctions import AssortativeMatingSimulation
    from save_simulation_data import save_simulation_results
except ImportError as e:
    print(f"Error: Could not import necessary components. {e}")
    print("Please ensure SimulationFunctions.py and save_simulation_data.py are in the same directory.")
    sys.exit(1)

def run_single_replication(task_params):
    """
    A worker function for the multiprocessing Pool.
    It runs ONE full simulation replication from start to finish and saves the results.
    """
    # Define a default name here in case an error happens before it's properly assigned
    run_folder_name_for_logging = task_params.get('condition_name', 'unknown') + '_run_' + str(task_params.get('replication_id', 'unknown'))
    
    try:
        # Unpack parameters for this specific replication
        replication_id = task_params['replication_id']
        condition_name = task_params['condition_name']
        base_output_folder = task_params['base_output_folder']
        
        # --- PATH AND FILENAME SETUP ---
        condition_folder = os.path.join(base_output_folder, condition_name)
        run_subfolder_name = f"run_{replication_id:03d}"
        run_output_folder = os.path.join(condition_folder, run_subfolder_name)
        
        # *** FIX #1: CREATE THE OUTPUT DIRECTORY BEFORE ANYTHING ELSE ***
        # This prevents the race condition where the script tries to write a file
        # before its directory has been created by another function.
        os.makedirs(run_output_folder, exist_ok=True)
        
        file_prefix = f"{condition_name}_run_{replication_id:03d}"
        summary_txt_filename = os.path.join(run_output_folder, f"{file_prefix}_summary.txt")

        # Create a unique, reproducible seed for each replication
        run_seed = task_params['simulation_params']['seed'] + replication_id
        
        # Update simulation parameters with the unique seed and the specific summary filename
        sim_params = task_params['simulation_params'].copy()
        sim_params['seed'] = run_seed
        sim_params['output_summary_filename'] = summary_txt_filename

        print(f"  -> Starting replication: {file_prefix} with seed {run_seed}")
        
        # --- Instantiate and Run Simulation ---
        sim_instance = AssortativeMatingSimulation(**sim_params)
        results = sim_instance.run_simulation()
        
        # --- Save All Results ---
        if results:
            save_simulation_results(
                results=results, 
                output_folder=run_output_folder, 
                file_prefix=file_prefix,
                scope=[1,19,20]
            )
        
        print(f"  -> Finished replication: {file_prefix}")
        return f"Success: {file_prefix}"

    except Exception as e:
        error_message = f"Failed: {run_folder_name_for_logging} with error: {e}"
        print(f"!!! ERROR in {run_folder_name_for_logging} !!!")
        import traceback
        traceback.print_exc()
        return error_message


def main():
    """
    Main function to define conditions and map a Slurm task ID to a BATCH of replications.
    """
    # --- 1. Main Configuration ---
    main_output_directory = "/scratch/alpine/xuly4739/StatRev_IndirectGene/Data/ASHG_Final"
    REPLICATIONS_PER_CONDITION = 1000 
    REPLICATIONS_PER_SLURM_TASK = 4

    # --- 2. Define Simulation Conditions ---
    # (This section is unchanged and correctly configured)
    base_params = {
        "num_generations": 20, "pop_size": 8e4, "n_CV": 300, "rg_effects": 0.1,
        "maf_min": 0.25, "maf_max": 0.45, "avoid_inbreeding": True,
        "save_each_gen": True, "save_covs": False, "summary_file_scope": "all",
        "seed": 202506
    }
    k2_val = np.array([[1.0, base_params["rg_effects"]], [base_params["rg_effects"], 1.0]]); d_mat_val = np.diag([np.sqrt(0.3), np.sqrt(0.2)]); a_mat_val = np.diag([np.sqrt(0.5), np.sqrt(0.6)]); fmat_val = np.array([[0,0],[0,0]]); s_mat_val = np.array([[0,0],[0,0]]); cove_val = np.array([[0.2, 0.05], [0.05, 0.2]]); covy_val = np.array([[1.0, 0.25], [0.25, 1.0]]); am_list_val = [np.array([[0.3, 0.05], [0.05, 0.3]])] * int(base_params["num_generations"])
    base_params.update({"k2_matrix": k2_val, "d_mat": d_mat_val, "a_mat": a_mat_val, "f_mat": fmat_val, "s_mat": s_mat_val, "cove_mat": cove_val, "covy_mat": covy_val, "am_list": am_list_val})
    f_mat_condition_A = np.array([[.10,.15],[.05,.15]]); s_mat_condition_A = np.array([[0,0],[0,0]])
    f_mat_condition_B = np.array([[0,0],[0,0]]); s_mat_condition_B = np.array([[.10,.15],[.05,.15]])
    simulation_conditions = [
        #{"condition_name": "phenoVT_phenoAM", "simulation_params": {**base_params, "mating_type": "phenotypic", "f_mat": f_mat_condition_A, "s_mat": s_mat_condition_A}},
        #{"condition_name": "socialVT_phenoAM", "simulation_params": {**base_params, "mating_type": "phenotypic", "f_mat": f_mat_condition_B, "s_mat": s_mat_condition_B}},
        {"condition_name": "phenoVT_socialAM", "simulation_params": {**base_params, "mating_type": "social", "f_mat": f_mat_condition_A, "s_mat": s_mat_condition_A}},
        {"condition_name": "phenoVT_geneticAM", "simulation_params": {**base_params, "mating_type": "genetic", "f_mat": f_mat_condition_A, "s_mat": s_mat_condition_A}}
    ]

    # --- 3. Determine Which Batch This Slurm Task Will Run (unchanged) ---
    try:
        task_id = int(os.environ.get('SLURM_ARRAY_TASK_ID', 1))
        num_cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', 1))
    except (ValueError, TypeError) as e:
        print(f"Could not read Slurm variables. Defaulting to task_id=1, cpus=1. Error: {e}")
        task_id = 1
        num_cpus = 1
        
    num_tasks_per_condition = REPLICATIONS_PER_CONDITION // REPLICATIONS_PER_SLURM_TASK
    condition_index = (task_id - 1) // num_tasks_per_condition
    block_index = (task_id - 1) % num_tasks_per_condition
    start_replication_id = block_index * REPLICATIONS_PER_SLURM_TASK + 1
    end_replication_id = start_replication_id + REPLICATIONS_PER_SLURM_TASK - 1
    
    if condition_index >= len(simulation_conditions):
        print(f"Error: Task ID {task_id} is out of bounds. Exiting."); sys.exit(1)

    current_condition = simulation_conditions[condition_index]

    # --- 4. Generate the List of Tasks, Skipping Completed Runs ---
    tasks_for_this_job = []
    for i in range(start_replication_id, end_replication_id + 1):
        
        # *** FIX #2: CHECK IF THE RUN DIRECTORY ALREADY EXISTS ***
        condition_name = current_condition["condition_name"]
        run_output_folder = os.path.join(main_output_directory, condition_name, f"run_{i:03d}")
        
        if os.path.isdir(run_output_folder):
            # A simple check is just for the folder, a more robust check could be
            # for a specific file, e.g., the summary.txt file.
            summary_file_path = os.path.join(run_output_folder, f"{condition_name}_run_{i:03d}_summary.txt")
            if os.path.exists(summary_file_path):
                 print(f"Skipping replication {i} for condition '{condition_name}': Output directory and summary file already exist.")
                 continue # Go to the next replication

        task = {
            "replication_id": i,
            "condition_name": condition_name,
            "base_output_folder": main_output_directory,
            "simulation_params": current_condition["simulation_params"]
        }
        tasks_for_this_job.append(task)
    
    # Check if there are any tasks left to run for this job
    if not tasks_for_this_job:
        print(f"--- Slurm Task {task_id}: All assigned replications already complete. Exiting. ---")
        return # Exit cleanly if all work is done

    print(f"--- Slurm Task {task_id} Configuration ---")
    print(f"Main output directory: {main_output_directory}")
    print(f"Condition: {current_condition['condition_name']}")
    print(f"Will run {len(tasks_for_this_job)} new replications from range: {start_replication_id} to {end_replication_id}")
    print(f"Using {num_cpus} CPUs for parallel processing within this task.")
    
    # --- 5. Run the Batch of Simulations in Parallel ---
    with multiprocessing.Pool(processes=num_cpus) as pool:
        pool.map(run_single_replication, tasks_for_this_job)

    print(f"--- Slurm Task {task_id} Finished ---")


if __name__ == '__main__':
    main()