#!/bin/sh
#SBATCH --job-name=lmk_knn_covid19mx_1_simulate_cohort  # Job name
#SBATCH --mail-type=ALL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jason.brunson@medicine.ufl.edu
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --cpus-per-task=1           # Use 1 core
#SBATCH --mem=1200mb                 # Memory limit
#SBATCH --time=01:00:00             # Time limit hrs:min:sec
#SBATCH --output=lmk_knn_covid19mx_1_simulate_cohort.out # Std output and error log

pwd; hostname; date

module load R

echo "Running `lmk_knn_covid19mx_1_simulate_cohort` on one core..."

Rscript ~/lastfirst/scripts/lmk-knn-covid19mx-1-simulate-cohort.r

date
