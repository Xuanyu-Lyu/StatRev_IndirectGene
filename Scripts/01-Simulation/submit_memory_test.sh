#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --job-name=memory_benchmark     # A specific name for this test job
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2               # <<< CRITICAL: Request 4 CPUs for the 4 parallel runs
#SBATCH --mem=256G                      # <<< CRITICAL: Request a large amount of memory to avoid OOM
#SBATCH --time=0-04:00:00               # 4-hour time limit (adjust if one run takes longer)
#SBATCH --chdir /projects/xuly4739/Py_Projects/StatRev_IndirectGene/Scripts/01-Simulation
#SBATCH --exclude bmem-rico1
#SBATCH --output=slurm_logs/memory_test_%A.out
#SBATCH --error=slurm_logs/memory_test_%A.err

module purge
source /curc/sw/anaconda3/latest
conda activate /projects/xuly4739/general_env

echo "------------------------------------------------"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Running memory benchmark..."
echo "------------------------------------------------"

python -u Benchmark_sim_memory.py

echo "Benchmark script finished."