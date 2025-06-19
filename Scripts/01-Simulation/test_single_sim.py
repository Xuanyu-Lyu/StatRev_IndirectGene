#!/usr/bin/env python3

import numpy as np
import os
from SimulationFunctions import AssortativeMatingSimulation

def main():
    """
    Run a single simulation with specified parameters:
    - s matrix with non-zero values but all less than 0.4
    - f matrix set to zero
    - Only save summary file, not all generations' data
    """
    
    # Base simulation parameters
    base_params = {
        "num_generations": 10,     # Number of generations to simulate
        "pop_size": 10000,          # Population size (smaller for testing)
        "n_CV": 100,               # Number of causal variants
        "rg_effects": 0,         # Genetic correlation between effects
        "maf_min": 0.25,           # Minimum minor allele frequency
        "maf_max": 0.45,           # Maximum minor allele frequency
        "avoid_inbreeding": True,
        "save_each_gen": False,    # Don't save all generations' data
        "save_covs": False,        # Don't save covariances
        "summary_file_scope": "final",  # Only save final summary
        "seed": 12345
    }
    
    # Define the matrices as per your requirements
    
    # k2_matrix: genetic correlation matrix
    k2_val = np.array([[1.0, base_params["rg_effects"]], 
                       [base_params["rg_effects"], 1.0]])
    
    # d_mat: environmental variance matrix (diagonal)
    d_mat_val = np.diag([np.sqrt(0), np.sqrt(0)])
    
    # a_mat: additive genetic variance matrix (diagonal)
    a_mat_val = np.diag([np.sqrt(0), np.sqrt(0)])
    
    # f_mat: ZERO matrix as requested
    f_mat_val = np.array([[0.0, 0.0], 
                          [0.0, 0.0]])
    
    # s_mat: NON-ZERO values but all less than 0.4 as requested
    s_mat_val = np.array([[0.3, 0.3], 
                          [0.3, 0.3]])  # All values < 0.4
    
    # cove_mat: environmental covariance matrix
    cove_val = np.array([[0.2, 0.05], 
                         [0.05, 0.2]])
    
    # covy_mat: phenotypic covariance matrix
    covy_val = np.array([[1.0, 0.25], 
                         [0.25, 1.0]])
    
    # am_list: assortative mating correlation matrices for each generation
    am_correlation = np.array([[0, 0], 
                               [0, 0]])
    am_list_val = [am_correlation] * base_params["num_generations"]
    
    # Update base parameters with matrices
    base_params.update({
        "k2_matrix": k2_val,
        "d_mat": d_mat_val,
        "a_mat": a_mat_val,
        "f_mat": f_mat_val,
        "s_mat": s_mat_val,
        "cove_mat": cove_val,
        "covy_mat": covy_val,
        "am_list": am_list_val,
        "mating_type": "phenotypic"  # Type of assortative mating
    })
    
    # Set up output file
    output_dir = "/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/coding/PyProject/StatRev_IndirectGene/Scripts/01-Simulation"
    summary_filename = os.path.join(output_dir, "test_simulation_summary.txt")
    base_params["output_summary_filename"] = summary_filename
    
    print("Starting simulation with parameters:")
    print(f"- Generations: {base_params['num_generations']}")
    print(f"- Population size: {base_params['pop_size']}")
    print(f"- Number of causal variants: {base_params['n_CV']}")
    print(f"- Mating type: {base_params['mating_type']}")
    print(f"- Save each generation: {base_params['save_each_gen']}")
    print(f"- Summary scope: {base_params['summary_file_scope']}")
    print(f"- F matrix (should be zeros): \\n{f_mat_val}")
    print(f"- S matrix (non-zero, all < 0.4): \\n{s_mat_val}")
    print()
    
    try:
        # Create and run the simulation
        sim = AssortativeMatingSimulation(**base_params)
        
        print("Simulation initialized successfully!")
        print("Running simulation...")
        
        # Run the simulation
        results = sim.run_simulation()
        
        if results is not None:
            print("\\nSimulation completed successfully!")
            print(f"Summary saved to: {summary_filename}")
            print(f"Final generation population size: {len(results['PHEN'])}")
            print(f"Number of summary results: {len(results['SUMMARY.RES'])}")
            
            # Print some basic statistics from the final generation
            final_summary = results['SUMMARY.RES'][-1]
            print(f"\\nFinal generation statistics:")
            print(f"- Generation: {final_summary['GEN']}")
            print(f"- Population size: {final_summary['POPSIZE']}")
            print(f"- Heritability (h2): {final_summary.get('h2', 'N/A')}")
            
        else:
            print("Simulation failed to complete.")
            
    except Exception as e:
        print(f"Error running simulation: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
