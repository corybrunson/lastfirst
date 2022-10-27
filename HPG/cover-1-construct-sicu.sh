#!/bin/sh
#SBATCH --job-name=sicu_covers      # Job name
#SBATCH --mail-type=ALL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jason.brunson@medicine.ufl.edu
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --cpus-per-task=1           # Use 1 core
#SBATCH --mem=1200mb                # Memory limit
#SBATCH --time=12:00:00             # Time limit hrs:min:sec
#SBATCH --output=cover-1-construct-sicu.out # Std output and error log

pwd; hostname; date

module load R

echo "Running `cover-1-construct-sicu.r` on one core..."

Rscript ~/lastfirst/scripts/cover-1-construct-sicu.r

date
