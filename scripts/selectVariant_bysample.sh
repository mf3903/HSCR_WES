#!/bin/bash

###############################################################################
# Script: selectVariant_bysample.sh
# Description: Selects variants for a specified sample set from a VCF using GATK.
# Usage: Submit as a SLURM job. Update input/output/reference/sample list paths as needed.
# Requirements: GATK 4.x, R, SLURM
###############################################################################

#SBATCH --job-name=selectvariant
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200gb
#SBATCH --time=24:00:00
#SBATCH --output=log/selectvariant_%j.log
#SBATCH --partition=cpu_medium

echo "[INFO] Job started at $(date)"

module purge
module load gatk/4.1.2.0
module load r/3.5.1

# Reference genome (update path as needed)
ref_genome=/path/to/reference_genome.fasta

# Input VCF and output directory (update paths as needed)
vcf_file=/path/to/input_variants.vcf.gz
outdir=/path/to/output_selectVariant_dir

# Sample list file (one sample name per line)
sample_list=/path/to/sample_list.txt

mkdir -p $outdir

gatk --java-options '-Xmx180g -Xms160g' \
    SelectVariants \
    -R $ref_genome \
    -V $vcf_file \
    -O ${outdir}/selected_variants.vcf.gz \
    -sn $sample_list \
    --exclude-filtered true \
    --exclude-non-variants true \
    --remove-unused-alternates true

echo "[STATUS] SelectVariants exit code: $?"

echo "[INFO] Job completed at $(date)"
