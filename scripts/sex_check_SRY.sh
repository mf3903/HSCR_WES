#!/bin/bash

###############################################################################
# Script: sex_check.sh
# Description: Estimates sex by calculating read depth over the SRY region on chrY.
#              Designed for SLURM array jobs. Each task processes one BAM sample.
# Usage: Submit as a SLURM array job. Reads BAM file paths from a text file.
# Requirements: Samtools, R, SLURM
###############################################################################

#SBATCH --job-name=sex_check
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --time=4:00:00
#SBATCH --output=logs-sex/sex_check%j.log
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

# Reference genome (update path as needed)
ref_genome=/path/to/reference_genome.fasta

# Output directories
out_dir=sex
mkdir -p $out_dir

# Calculate read depth over SRY region on chrY
samtools depth -r chrY:2786840-2787751 --reference $ref_genome $line -a -q 0 -Q 0 > $out_dir/${sample_id}.SRY.readdepth.txt

# Summarize read depth using R
cut -f3 $out_dir/${sample_id}.SRY.readdepth.txt | Rscript -e 'print(summary(scan("stdin")));' > $out_dir/${sample_id}.SRY.readdepth.summary.txt

echo "[STATUS] samtools exit code: $?"

echo "[INFO] Job completed at $(date)"
