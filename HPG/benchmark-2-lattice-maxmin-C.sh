#!/bin/sh
#SBATCH --job-name=bm_lat_mm_C      # Job name
#SBATCH --mail-type=ALL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jason.brunson@medicine.ufl.edu
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --cpus-per-task=1           # Use 1 core
#SBATCH --mem=1200mb                # Memory limit
#SBATCH --time=01:00:00             # Time limit hrs:min:sec
#SBATCH --output=benchmark-2-lattice-maxmin-C.out # Std output and error log

pwd; hostname; date

module load R

echo "Running `bm_lat_mm_C` on one core..."

Rscript ~/lastfirst/scripts/benchmark-2-lattice-maxmin-C.r

date
