#!/bin/bash

###############################################################################
# Script: logit_firth.sh
# Description: Submits an R script for Firth logistic regression analysis.
# Usage: Submit as a SLURM job. Update R script variable as needed.
# Requirements: R, SLURM
###############################################################################

#SBATCH --job-name=rscript
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96gb
#SBATCH --time=72:00:00
#SBATCH --output=log/rscript_%j.log
#SBATCH --partition=a100_short,a100_long

echo "[INFO] Job started at $(date)"
module load r/4.1.2

# Specify the R script to run (update as needed)
rscript=logit_firth.R 

Rscript $rscript

echo "[STATUS] Rscript exit code: $?"
echo "[INFO] Job completed at $(date)"



echo _ESTATUS_ [ rscript ]: $?
echo _END_  $(date)




