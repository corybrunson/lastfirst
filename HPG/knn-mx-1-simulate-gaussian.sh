#!/bin/sh
#SBATCH --job-name=knn_mx_sim_gauss # Job name
#SBATCH --mail-type=ALL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jason.brunson@medicine.ufl.edu
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --cpus-per-task=1           # Use 1 core
#SBATCH --mem=3200mb                # Memory limit
#SBATCH --time=12:00:00             # Time limit hrs:min:sec
#SBATCH --output=knn-mx-1-simulate-gaussian.out # Std output and error log

pwd; hostname; date

module load R

echo "Running `knn-mx-1-simulate-gaussian.r` on one core..."

Rscript ~/lastfirst/scripts/knn-mx-1-simulate-gaussian.r

date
