#!/bin/bash

###############################################################################
# Script: relatedness.R.sh
# Description: Submits an R script for pairwise relatedness analysis on a SLURM cluster.
# Usage: Update the R script path as needed.
# Requirements: R, SLURM
###############################################################################

#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=120:00:00
#SBATCH --mem=96GB
#SBATCH --job-name=Rscript
#SBATCH --mail-type=END,FAIL
#SBATCH --output=log/rscript%j.log
#SBATCH --partition=cpu_medium,cpu_long

echo "[INFO] Job started at $(date)"

module load r/3.6.1-x11

# Run the relatedness R script (update script path as needed)
Rscript relatedness.R

echo "[STATUS] Rscript exit code: $?"
echo "[INFO] Job completed at $(date)"
