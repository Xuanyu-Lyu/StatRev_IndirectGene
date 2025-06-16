#!/bin/bash
# --- Cluster Specific Settings (from your script) ---
#SBATCH --qos=preemptable
#SBATCH --chdir /projects/xuly4739/Py_Projects/StatRev_IndirectGene/Scripts/01-Simulation
#SBATCH --exclude bmem-rico1

# --- Job Resource Settings ---
# These are set for a single test run on one replication.
# The job name is changed to distinguish it from your main analysis runs.
# Time and memory should be generous enough for one large run to complete.
#SBATCH --job-name=gcta_rdr_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2               # Appropriate for GCTA's --thread-num 2 option
#SBATCH --mem=200G                      # A safe, large memory allocation for testing
#SBATCH --time=0-04:00:00               # 4-hour time limit for the test run

# --- Log File Settings ---
#SBATCH --output=slurm_logs/gcta_test_%A.out      # Specific log file for this test
#SBATCH --error=slurm_logs/gcta_test_%A.err       # Specific log file for errors

# --- Start of Job Commands ---

# Ensure a clean environment before loading your own
module purge

# Load the same Conda environment you use for your other jobs
source /curc/sw/anaconda3/latest
conda activate /projects/xuly4739/general_env

# Create the log directory if it doesn't exist
mkdir -p slurm_logs

# Print job information to the log for easy debugging
echo "------------------------------------------------"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Running on host: $(hostname)"
echo "Working Directory: $(pwd)"
echo "------------------------------------------------"

# Make the runner script executable, just in case
chmod +x run_gcta_rdr.sh

# --- Define Paths and Run the Test ---
# Define the base directory of your simulation output
BASE_SIM_DIR="/scratch/alpine/xuly4739/StatRev_IndirectGene/Data/ASHG_Final/phenoVT_geneticAM"
# Define a dedicated output directory for this test's results
BASE_OUTPUT_DIR="/projects/xuly4739/Py_Projects/StatRev_IndirectGene/Analysis/RDR_Results_Test"

echo "--- Running test for run_001 ---"
# Call your worker script with the path to the simulation data and where to save results
./run_gcta_rdr.sh "${BASE_SIM_DIR}/run_001" ${BASE_OUTPUT_DIR}

echo "--- Running test for run_002 ---"
# You can uncomment this line to run a second test case
./run_gcta_rdr.sh "${BASE_SIM_DIR}/run_002" ${BASE_OUTPUT_DIR}

echo "--- Test job finished ---"