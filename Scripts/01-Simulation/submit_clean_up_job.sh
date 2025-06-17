#!/bin/bash
# --- Cluster Specific Settings ---
#SBATCH --qos=preemptable
#SBATCH --chdir /projects/xuly4739/Py_Projects/StatRev_IndirectGene/Scripts/01-Simulation
#SBATCH --exclude bmem-rico1

#!/bin/bash
#SBATCH --job-name=cleanup_gcta         # A name for your cleanup job
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00                 # 30 minutes should be more than enough
#SBATCH --mem=2G                        # Memory requirement is minimal
#SBATCH --output=slurm_logs/cleanup_gcta_%A.out
#SBATCH --error=slurm_logs/cleanup_gcta_%A.err

# This command ensures that the script will exit immediately if any command fails
set -e

# --- CONFIGURATION ---
# !!! CRITICAL: Make sure this path is exactly correct before running !!!
# This is the top-level directory where the script will start searching.
BASE_DATA_DIR="/scratch/alpine/xuly4739/StatRev_IndirectGene/Data/ASHG_Final"


# --- SCRIPT LOGIC ---

# Create the log directory if it doesn't exist
mkdir -p slurm_logs

echo "--- GCTA Analysis Cleanup Job ---"
echo "Job ID: $SLURM_JOB_ID"
echo "Searching in: ${BASE_DATA_DIR}"
echo ""

# Step 1: Find all target directories and print them to the log file.
# This creates a record of what was targeted for deletion.
echo "--- Finding all directories named 'gcta_analysis' to delete... ---"
find "${BASE_DATA_DIR}" -type d -name "gcta_analysis" -print
echo "----------------------------------------------------------------"

# Step 2: Find the directories again and execute the deletion.
# The script will proceed automatically without asking for confirmation.
echo "Now deleting the directories listed above..."

# The -exec flag runs the 'rm -rf' command on each found item.
# The '+' at the end groups paths together for efficiency.
find "${BASE_DATA_DIR}" -type d -name "gcta_analysis" -exec rm -rf {} +

echo "Cleanup complete."