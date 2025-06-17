#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --job-name=cleanup_loop
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --chdir /projects/xuly4739/Py_Projects/StatRev_IndirectGene/Scripts/01-Simulation
#SBATCH --output=slurm_logs/cleanup_all_runs_%A.out
#SBATCH --error=slurm_logs/cleanup_all_runs_%A.err

# --- Job Start ---
echo "------------------------------------------------"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Running on host: $(hostname)"
echo "Job started at: $(date)"
echo "------------------------------------------------"

# --- Configuration ---
# The base directory where all condition folders are located.
BASE_DATA_DIR="/scratch/alpine/xuly4739/StatRev_IndirectGene/Data/ASHG_Final"
# The total number of generations in the simulation (e.g., 20 means gens 0-19).
TOTAL_GENERATIONS=20

# --- Calculate Generation Indices ---
# Generations are 0-indexed. If TOTAL_GENERATIONS is 20, indices are 0 to 19.
MAX_GEN_INDEX=$((TOTAL_GENERATIONS - 1))
PENULTIMATE_GEN_INDEX=$((TOTAL_GENERATIONS - 2))

# Define the range of generations to delete.
# We will delete from generation 1 up to (but not including) the penultimate one.
START_DELETE_GEN=1
END_DELETE_GEN=$((PENULTIMATE_GEN_INDEX - 1))

echo "Cleanup Configuration:"
echo "  - Base Directory: $BASE_DATA_DIR"
echo "  - Keeping generation: 0"
echo "  - Keeping generations: $PENULTIMATE_GEN_INDEX and $MAX_GEN_INDEX"
echo "  - Deleting generations from $START_DELETE_GEN to $END_DELETE_GEN"

# --- Find all target directories ---
echo
echo "Finding all 'run_*' directories..."
# The 'find' command populates an array with all the directory paths.
readarray -t ALL_RUN_DIRS < <(find "$BASE_DATA_DIR" -mindepth 2 -maxdepth 2 -type d -name "run_*" | sort)

NUM_DIRS=${#ALL_RUN_DIRS[@]}
if [ "$NUM_DIRS" -eq 0 ]; then
    echo "No 'run_*' directories found. Exiting."
    exit 0
fi
echo "Found $NUM_DIRS directories to process."


# --- Main Cleanup Loop ---
PROCESSED_COUNT=0
for dir in "${ALL_RUN_DIRS[@]}"; do
    PROCESSED_COUNT=$((PROCESSED_COUNT + 1))
    echo "----------------------------------------"
    echo "Processing directory ($PROCESSED_COUNT / $NUM_DIRS): $dir"

    # Loop through the intermediate generation numbers to delete
    for i in $(seq $START_DELETE_GEN $END_DELETE_GEN); do
        # Use find to locate the files for a specific generation within the directory
        # The \( ... \) groups the -o conditions.
        files_to_delete=$(find "$dir" -maxdepth 1 -type f \( -name "*_gen${i}.tsv" -o -name "*_gen${i}_*.tsv" \))

        # If the 'find' command returned any file paths...
        if [ -n "$files_to_delete" ]; then
            echo "  Deleting files for generation $i:"
            # Use xargs to safely remove the found files.
            # -t echoes the command for the log, -r prevents running if input is empty.
            echo "$files_to_delete" | xargs -r -t rm
        fi
    done
    echo "Finished processing $dir"
done

echo "========================================"
echo "Cleanup script finished."
echo "Job ended at: $(date)"