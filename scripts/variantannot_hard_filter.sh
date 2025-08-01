#!/bin/bash

###############################################################################
# Script: variantannot_hard_filter.sh
# Description: Annotates variants in a VCF file using GATK VariantAnnotator.
#              Designed for SLURM jobs. Processes a single VCF file.
# Usage: Submit as a SLURM job. Update input/output/reference paths as needed.
# Requirements: GATK 3.8, Java, SLURM
###############################################################################

#SBATCH --job-name=variantannot
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200gb
#SBATCH --time=5-00:00:00
#SBATCH --output=log/variantannot_%j.log
#SBATCH --partition=cpu_medium

module purge
module load gatk/3.8.0

# Reference genome (update path as needed)
ref_genome=/path/to/reference_genome.fasta

# Input VCF and directories (update paths as needed)
in_dir=/path/to/input_vcf_dir
in_vcf=$in_dir/input_variants.vcf.gz

out_dir=/path/to/output_annot_dir

# GATK resource bundle (update path as needed)
gatk_bundle_dir=/path/to/gatk_bundle
dbsnp=$gatk_bundle_dir/dbsnp.vcf

# Run GATK VariantAnnotator
java -Xmx180g -Xms160g -jar /path/to/GenomeAnalysisTK.jar \
    -T VariantAnnotator \
    -R $ref_genome \
    --variant $in_vcf \
    -L $in_vcf \
    --out $out_dir/variants.extraannot.vcf.gz \
    -A GCContent -A HomopolymerRun -A HardyWeinberg \
    --dbsnp $dbsnp \
    --reference_window_stop 200

echo "[STATUS] VariantAnnotator exit code: $?"

echo "[INFO] Job completed at $(date)"

# Example of original usage:
# java -Xmx150g -jar GenomeAnalysisTK.jar -R hs37d5.fa -T VariantAnnotator --out hscrPlusControlsWithExtraAnnotation.vcf -A GCContent -A HomopolymerRun -A HardyWeinberg --variant hscrPlusControls.vcf -L hscrPlusControls.vcf --dbsnp dbsnp_135.b37.vcf
