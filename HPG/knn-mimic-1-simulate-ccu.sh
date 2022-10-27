#!/bin/sh
#SBATCH --job-name=knn_mimic_ccu    # Job name
#SBATCH --mail-type=ALL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jason.brunson@medicine.ufl.edu
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --cpus-per-task=1           # Use 1 core
#SBATCH --mem=600mb                # Memory limit
#SBATCH --time=12:00:00             # Time limit hrs:min:sec
#SBATCH --output=knn-mimic-1-simulate-ccu.out # Std output and error log

pwd; hostname; date

module load R

echo "Running `knn_mimic_ccu` on one core..."

Rscript ~/lastfirst/scripts/knn-mimic-1-simulate-ccu.r

date
