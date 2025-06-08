# File: run_slurm_batch_job.py

import os
import sys
import multiprocessing
import datetime
import numpy as np
import pandas as pd

# Import the necessary class and functions from your other files
try:
    from SimulationFunctions import AssortativeMatingSimulation
    from save_simulation_data import *
except ImportError as e:
    print(f"Error: Could not import necessary components. {e}")
    print("Please ensure SimulationFunctions.py and save_simulation_data.py are in the same directory.")
    sys.exit(1)

def run_single_replication(task_params):
    """
    A worker function for the multiprocessing Pool.
    It runs ONE full simulation replication from start to finish and saves the results
    into a structured, condition-specific directory.

    Args:
        task_params (dict): A dictionary containing all parameters for a single replication.
    """
    # Define a default name here in case an error happens before it's properly assigned
    run_folder_name_for_logging = task_params.get('condition_name', 'unknown') + '_run_' + str(task_params.get('replication_id', 'unknown'))
    
    try:
        # Unpack parameters for this specific replication
        replication_id = task_params['replication_id']
        condition_name = task_params['condition_name']
        base_output_folder = task_params['base_output_folder']
        
        # --- CORRECTION IS HERE ---
        # 1. Define the path for the condition-specific folder
        condition_folder = os.path.join(base_output_folder, condition_name)
        
        # 2. Define the path for this specific run's subfolder
        #    Let's call the directory itself something simple like "run_001"
        run_subfolder_name = f"run_{replication_id:03d}"
        run_output_folder = os.path.join(condition_folder, run_subfolder_name)
        
        # 3. Define the descriptive name for logging and file prefixes
        #    This is the variable that was missing.
        run_folder_name = f"{condition_name}_run_{replication_id:03d}"

        # 4. Define the path for the summary text file for this run
        summary_txt_filename = os.path.join(run_output_folder, "run_summary.txt")

        # Create a unique, reproducible seed for each replication
        run_seed = task_params['simulation_params']['seed'] + replication_id
        
        # Update simulation parameters with the unique seed and summary filename
        sim_params = task_params['simulation_params'].copy()
        sim_params['seed'] = run_seed
        sim_params['output_summary_filename'] = summary_txt_filename

        # This print statement now correctly uses the defined 'run_folder_name' variable
        print(f"  -> Starting replication: {run_folder_name} (in folder {condition_name}) with seed {run_seed}")
        
        # --- Instantiate and Run Simulation ---
        sim_instance = AssortativeMatingSimulation(**sim_params)
        results = sim_instance.run_simulation()
        
        # --- Save All Results ---
        if results:
            save_simulation_results(
                results=results, 
                output_folder=run_output_folder, 
                file_prefix=run_folder_name, # Use the descriptive name for file prefixes
                scope="all"
            )
        
        print(f"  -> Finished replication: {run_folder_name} (in folder {condition_name})")
        return f"Success: {run_folder_name}"

    except Exception as e:
        # The logging variable is now defined at the top, so it will exist even if an early error occurs
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
    
    # Define a base output folder for this entire batch of simulations
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    main_output_directory = f"/scratch/alpine/xuly4739/StatRev_IndirectGene/Data/TestRc2"
    
    REPLICATIONS_PER_CONDITION = 100 
    REPLICATIONS_PER_SLURM_TASK = 10

    # --- 2. Define Simulation Conditions ---
    # Each dictionary in this list is a separate experimental condition.
    base_params = {
        "num_generations": 15, "pop_size": 4e4, "n_CV": 300, "rg_effects": 0.1,
        "maf_min": 0.25, "maf_max": 0.45, "avoid_inbreeding": True,
        "save_each_gen": True, "save_covs": False, "summary_file_scope": "all",
        "seed": 202506, "mating_type": "phenotypic"
    }
    k2_val = np.array([[1.0, base_params["rg_effects"]], [base_params["rg_effects"], 1.0]])
    d_mat_val = np.diag([np.sqrt(0.3), np.sqrt(0.2)]); a_mat_val = np.diag([np.sqrt(0.5), np.sqrt(0.6)])
    fmat_val = np.array([[0,0],[0,0]]); s_mat_val = np.array([[0,0],[0,0]])
    cove_val = np.array([[0.2, 0.05], [0.05, 0.2]]); covy_val = np.array([[1.0, 0.25], [0.25, 1.0]])
    am_list_val = [np.array([[0.3, 0.05], [0.05, 0.3]])] * base_params["num_generations"]
    base_params.update({
        "k2_matrix": k2_val, "d_mat": d_mat_val, "a_mat": a_mat_val, "f_mat": fmat_val, 
        "s_mat": s_mat_val, "cove_mat": cove_val, "covy_mat": covy_val, "am_list": am_list_val
    })
    
    # *** DEFINE DIFFERENT f_mat AND s_mat FOR EACH CONDITION ***
    f_mat_condition_A = np.array([[.10,.15],[.05,.15]]) 
    s_mat_condition_A = np.array([[0,0],[0,0]]) 

    f_mat_condition_B = np.array([[0,0],[0,0]])
    s_mat_condition_B = np.array([[.10,.15],[.05,.15]])

    simulation_conditions = [
        {
            "condition_name": "phenotypic_transmission",
            "simulation_params": {
                **base_params, 
                "mating_type": "phenotypic",
                "f_mat": f_mat_condition_A, # Use the first set of matrices
                "s_mat": s_mat_condition_A
            }
        },
        {
            "condition_name": "social_transmission",
            "simulation_params": {
                **base_params, 
                "mating_type": "phenotypic",
                "f_mat": f_mat_condition_B, # Use the second set of matrices
                "s_mat": s_mat_condition_B
            }
        }
    ]

    # --- 3. Determine Which Batch This Slurm Task Will Run ---
    try:
        # Slurm provides a 1-based task ID.
        task_id = int(os.environ.get('SLURM_ARRAY_TASK_ID', 1))
        # Determine number of CPUs allocated to this task by Slurm
        num_cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', 1))
    except (ValueError, TypeError) as e:
        print(f"Could not read Slurm environment variables. Defaulting to task_id=1, cpus=1. Error: {e}")
        task_id = 1
        num_cpus = 1
        
    # Total number of tasks per condition
    num_tasks_per_condition = REPLICATIONS_PER_CONDITION // REPLICATIONS_PER_SLURM_TASK

    # Determine which condition and which block of replications this task handles
    condition_index = (task_id - 1) // num_tasks_per_condition
    block_index = (task_id - 1) % num_tasks_per_condition
    
    start_replication_id = block_index * REPLICATIONS_PER_SLURM_TASK + 1
    end_replication_id = start_replication_id + REPLICATIONS_PER_SLURM_TASK - 1
    
    if condition_index >= len(simulation_conditions):
        print(f"Error: Task ID {task_id} is out of bounds for the defined conditions. Exiting.")
        sys.exit(1)

    current_condition = simulation_conditions[condition_index]

    # --- 4. Generate the List of 10 Tasks for This Slurm Job ---
    tasks_for_this_job = []
    for i in range(start_replication_id, end_replication_id + 1):
        task = {
            "replication_id": i,
            "condition_name": current_condition["condition_name"],
            "base_output_folder": main_output_directory,
            "simulation_params": current_condition["simulation_params"]
        }
        tasks_for_this_job.append(task)
    
    print(f"--- Slurm Task {task_id} Configuration ---")
    print(f"Main output directory: {main_output_directory}")
    print(f"Condition: {current_condition['condition_name']}")
    print(f"Running replications: {start_replication_id} through {end_replication_id}")
    print(f"Using {num_cpus} CPUs for parallel processing within this task.")
    
    # Create the main output directory if it doesn't exist
    # All tasks will try to create it, which is fine.
    os.makedirs(main_output_directory, exist_ok=True)

    # --- 5. Run the Batch of 10 Simulations in Parallel using multiprocessing ---
    with multiprocessing.Pool(processes=num_cpus) as pool:
        pool.map(run_single_replication, tasks_for_this_job)

    print(f"--- Slurm Task {task_id} Finished ---")


if __name__ == '__main__':
    main()