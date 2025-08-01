#!/bin/bash

###############################################################################
# Script: sex_check_ratio.sh
# Description: Estimates sex by calculating X and Y chromosome coverage ratios.
#              Designed for SLURM array jobs. Each task processes one BAM sample.
# Usage: Submit as a SLURM array job. Reads BAM file paths from a text file.
# Requirements: Samtools, R, SLURM, bc
###############################################################################

#SBATCH --job-name=sex_check
#SBATCH --mem=10gb
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=sex/sex_ratio%j.log
#SBATCH --array=1-N  # Set N to the number of samples
#SBATCH --partition=cpu_dev,cpu_short,cpu_medium,cpu_long,fn_short,fn_medium,fn_long

echo "[INFO] Job started at $(date)"

# Load required modules
module load samtools/1.9
module load r/4.0.3

# Input BAM list (one BAM path per line)
table=samples.bam-dedup-bqsr-mapped.txt
line="$(head -n $SLURM_ARRAY_TASK_ID $table | tail -n 1)"
sample_bam="$(basename $line)"
sample_id=$(basename $sample_bam .bam)

echo "[INFO] Processing array index: $SLURM_ARRAY_TASK_ID sample: $sample_id"

# Calculate X and Y chromosome coverage ratios
x_cov=$(echo "scale=4; $(samtools idxstats $line | awk '$1==\"chrX\" {print $3}')/$(samtools idxstats $line | awk '$1==\"chrX\" {print $2}')" | bc)
y_cov=$(echo "scale=4; $(samtools idxstats $line | awk '$1==\"chrY\" {print $3}')/$(samtools idxstats $line | awk '$1==\"chrY\" {print $2}')" | bc)

echo "[INFO] X chromosome coverage: $x_cov"
echo "[INFO] Y chromosome coverage: $y_cov"

xy_ratio=$(echo "scale=4; ${x_cov}/${y_cov}" | bc)
echo "[RESULT] X:Y coverage ratio for $sample_id: $xy_ratio"

echo "[INFO] Job completed at $(date)"




