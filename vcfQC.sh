#!/bin/bash

###############################################################################
# Script: vcfQC.sh
# Description: Runs GATK CollectVariantCallingMetrics on VCF files for QC.
# Usage: Submit as a SLURM job. Update input/output/reference paths as needed.
# Requirements: GATK 4.x, R, SLURM
###############################################################################

#SBATCH --job-name=CollectVariantCallingMetrics
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200gb
#SBATCH --time=12:00:00
#SBATCH --output=log/CollectVariantCallingMetrics_%j.log
#SBATCH --partition=cpu_short

echo "[INFO] Job started at $(date)"

module purge
module load gatk/4.1.2.0
module load r/3.5.1

# Reference genome and dictionary (update paths as needed)
ref_genome=/path/to/reference_genome.fasta
ref_dict=/path/to/reference.dict

# Input VCF directory and files (update as needed)
vcf_dir=/path/to/vcf_dir
vcf_file=input_variants.vcf.gz

# Input VCF file path
vcf_filein=$vcf_dir/$vcf_file

# GATK resource bundle (update path as needed)
gatk_bundle_dir=/path/to/gatk_bundle
dbsnp=$gatk_bundle_dir/dbsnp.vcf

# Output directory and file
outdir=/path/to/output_qc_dir
outfile=$outdir/variant_calling_metrics.QC

mkdir -p $outdir

# Run GATK CollectVariantCallingMetrics
gatk --java-options '-Xmx180g -Xms160g' \
    CollectVariantCallingMetrics \
    -I $vcf_filein \
    --DBSNP $dbsnp \
    -O $outfile

echo "[STATUS] CollectVariantCallingMetrics exit code: $?"

echo "[INFO] Job completed at $(date)"
