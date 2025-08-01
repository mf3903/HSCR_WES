#!/bin/bash

###############################################################################
# Script: haplocheck.sh
# Description: Runs haplocheck for contamination detection on BAM files.
# Usage: Submit as a SLURM job. Update input/output/reference paths as needed.
# Requirements: haplocheck, SLURM
###############################################################################

#SBATCH --job-name=haplocheck
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=11:00:00
#SBATCH --output=logs-contamination/haplocheck_%j.log
#SBATCH --partition=cpu_short,cpu_medium,cpu_long

# Reference: https://mitoverse.readthedocs.io/haplocheck/haplocheck/

echo "[INFO] Job started at $(date)"

# Directory containing sorted BAM and corresponding BAI files (update path as needed)
bam_dir=/path/to/bam_sorted_dir

# Ensure BAM and BAI files are together in the input directory
# cp /path/to/source/*sorted* $bam_dir

# Change to directory where haplocheck is installed (update as needed)
cd /path/to/haplocheck_install_dir

# Run haplocheck on all BAM files in the directory
./cloudgene run haplocheck@1.3.2 --files $bam_dir --format bam --output results --threads 4

echo "[STATUS] haplocheck exit code: $?"

echo "[INFO] Job completed at $(date)"
