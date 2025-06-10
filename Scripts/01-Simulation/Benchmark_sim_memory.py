# File: benchmark_memory.py

import os
import sys
import multiprocessing
import time
import numpy as np
import pandas as pd

try:
    from SimulationFunctions import AssortativeMatingSimulation
    # We don't need the saving functions for this memory test
except ImportError as e:
    print(f"Error: Could not import AssortativeMatingSimulation. {e}")
    sys.exit(1)

def run_one_test_simulation(task_params):
    """
    A worker function that runs a single, large-scale simulation.
    """
    run_id = task_params['run_id']
    sim_params = task_params['simulation_params']
    
    # Use a unique seed for each parallel run
    sim_params['seed'] = sim_params['seed'] + run_id
    
    print(f"  -> Starting test replication #{run_id} with pop_size={sim_params['pop_size']}...")
    start_time = time.time()
    
    try:
        sim_instance = AssortativeMatingSimulation(**sim_params)
        results = sim_instance.run_simulation()
        end_time = time.time()
        
        elapsed = end_time - start_time
        print(f"  -> Finished test replication #{run_id}. Time taken: {elapsed:.2f} seconds.")
        return f"Success: Run {run_id}"
    except Exception as e:
        print(f"!!! ERROR in test replication #{run_id}: {e} !!!")
        import traceback
        traceback.print_exc()
        return f"Failure: Run {run_id}"


def main():
    """
    Sets up and runs 2 parallel simulations for benchmarking.
    """
    # --- 1. Configuration for the Benchmark ---
    
    # Define the parameters for the large-scale run you want to test
    # This uses the same parameters from your run_simulations.py script
    benchmark_params = {
        "num_generations": 20,
        "pop_size": 8e4, # <-- The large population size you want to test
        "n_CV": 300, 
        "rg_effects": 0.1,
        "maf_min": 0.25, 
        "maf_max": 0.45, 
        "avoid_inbreeding": True,
        "mating_type": "phenotypic",
        "seed": 202508,
        # Disable file I/O for a pure memory/CPU test
        "save_each_gen": False, 
        "save_covs": False, 
        "output_summary_filename": None,
        "summary_file_scope": "final"
    }

    k2_val = np.array([[1.0, benchmark_params["rg_effects"]], [benchmark_params["rg_effects"], 1.0]])
    d_mat_val = np.diag([np.sqrt(0.3), np.sqrt(0.2)])
    a_mat_val = np.diag([np.sqrt(0.5), np.sqrt(0.6)])
    fmat_val = np.array([[0.10, 0.15], [0.05, 0.15]]) # Using one of your example matrices
    s_mat_val = np.array([[0,0],[0,0]])
    cove_val = np.array([[0.2, 0.05], [0.05, 0.2]])
    covy_val = np.array([[1.0, 0.25], [0.25, 1.0]])
    am_list_val = [np.array([[0.3, 0.05], [0.05, 0.3]])] * benchmark_params["num_generations"]
    
    benchmark_params.update({
        "k2_matrix": k2_val, "d_mat": d_mat_val, "a_mat": a_mat_val, "f_mat": fmat_val, 
        "s_mat": s_mat_val, "cove_mat": cove_val, "covy_mat": covy_val, "am_list": am_list_val
    })
    
    # --- 2. Set up the 2 Parallel Tasks ---
    num_parallel_runs = 2
    tasks = []
    for i in range(num_parallel_runs):
        tasks.append({
            "run_id": i + 1,
            "simulation_params": benchmark_params
        })
        
    # Get number of CPUs allocated by Slurm for the pool
    num_cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', num_parallel_runs))
    
    print(f"--- Starting Memory Benchmark ---")
    print(f"Running {num_parallel_runs} simulations in parallel with pop_size={int(benchmark_params['pop_size'])}.")
    print(f"Using a pool of {num_cpus} worker processes.")

    # --- 3. Run the simulations in parallel ---
    with multiprocessing.Pool(processes=num_cpus) as pool:
        pool.map(run_one_test_simulation, tasks)

    print("\n--- Benchmark Run Finished ---")
    print("Check job completion status and use 'sacct' to find memory usage.")


if __name__ == '__main__':
    main()