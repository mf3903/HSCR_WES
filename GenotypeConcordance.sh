#!/bin/bash

###############################################################################
# Script: GenotypeConcordance.sh
# Description: Runs Picard GenotypeConcordance to compare test and truth VCFs.
# Usage: Submit as a SLURM job. Update input/output/reference paths as needed.
# Requirements: Picard, GATK, R, SLURM
###############################################################################

#SBATCH --job-name=genotypeconcordance
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60gb
#SBATCH --time=5-00:00:00
#SBATCH --output=log/genotypeconcordance_%j.log
#SBATCH --partition=cpu_medium

echo "[INFO] Job started at $(date)"

module purge
module load gatk/4.2.1.0
module load r/3.5.1

# Input directories and files (update paths as needed)
sample_dir=/path/to/selectVariant_dir
test_sample_vcf=$sample_dir/test_sample.vcf.gz
truth_sample_vcf=$sample_dir/truth_sample.vcf.gz

# Output file
output_vcf=genotype_concordance_output.vcf

# Run Picard GenotypeConcordance
java -Xmx48g -Xms48g -jar /path/to/picard.jar GenotypeConcordance \
    CALL_VCF=$test_sample_vcf \
    O=$output_vcf \
    TRUTH_VCF=$truth_sample_vcf

echo "[STATUS] GenotypeConcordance exit code: $?"

echo "[INFO] Job completed at $(date)"
