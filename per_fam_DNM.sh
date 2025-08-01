#!/bin/bash

###############################################################################
# Script: per_fam_DNM.sh
# Description: Performs per-family variant selection, genotype refinement, and
#              de novo mutation (DNM) annotation using GATK.
# Usage: Submit as a SLURM job. Update input/output/reference/sample list paths as needed.
# Requirements: GATK 4.x, SLURM
###############################################################################

#SBATCH --job-name=selectvariant
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48gb
#SBATCH --time=12:00:00
#SBATCH --output=log/per_fam_DNM_entire_flow_%j.log
#SBATCH --partition=gpu4_short,gpu4_medium,gpu4_long,gpu8_short,gpu8_medium,gpu8_long

echo "[INFO] Job started at $(date)"

module purge
module load gatk/4.1.2.0

# Define family/sample input (update as needed)
famID=FAMILY_ID
per_fam_info_dir=/path/to/per_fam_sampleID
per_fam_wesID_list=${per_fam_info_dir}/${famID}.list

echo "[INFO] Processing array index: $SLURM_ARRAY_TASK_ID famID: $famID"

# Reference genome (update path as needed)
ref_genome=/path/to/reference_genome.fasta

# Input VCF and output directories (update as needed)
vcf_dir=/path/to/input_vcf_dir
vcf_file=$vcf_dir/input_variants.vcf.gz

output=/path/to/per_fam_vcf
mkdir -p $output/$famID
outdir=$output/$famID
my_sample=$per_fam_wesID_list

# Select variants for family samples
gatk --java-options '-Xmx42g -Xms36g' \
    SelectVariants \
    -R $ref_genome \
    -V $vcf_file \
    -O ${outdir}/selected_${famID}.vcf.gz \
    -sn $my_sample \
    --restrict-alleles-to BIALLELIC \
    --exclude-filtered true \
    --exclude-non-variants true \
    --remove-unused-alternates true

echo "[STATUS] SelectVariants exit code: $?"

# Genotype refinement using family pedigree
per_fam_vcf_dir=/path/to/per_fam_vcf
per_fam_vcf=${per_fam_vcf_dir}/${famID}/selected_${famID}.vcf.gz
per_fam_ped_dir=/path/to/per_fam_ped
per_fam_ped=${per_fam_ped_dir}/${famID}.ped

ref_genome=/path/to/reference_genome.fasta

outdir=/path/to/per_fam_posterior_GT
mkdir -p ${outdir}/${famID}
outdir_per_fam=${outdir}/${famID}

gatk --java-options '-Xmx40g -Xms32g -XX:ParallelGCThreads=2' \
    CalculateGenotypePosteriors \
    -V $per_fam_vcf \
    -ped $per_fam_ped \
    -XL chrX \
    -XL chrY \
    -O ${outdir_per_fam}/posterior_GT_${famID}.vcf.gz \
    --skip-population-priors

echo "[STATUS] CalculateGenotypePosteriors exit code: $?"

# Tag variants with low GQ
per_fam_vcf_tagged=${outdir_per_fam}/posterior_GT_${famID}.tagged.vcf.gz

gatk --java-options '-Xmx40g -Xms32g -XX:ParallelGCThreads=2' \
    VariantFiltration \
    -R $ref_genome \
    -V ${outdir_per_fam}/posterior_GT_${famID}.vcf.gz \
    -O $per_fam_vcf_tagged \
    --genotype-filter-expression "GQ < 20.0" \
    --genotype-filter-name "low_GQ_20"

echo "[STATUS] VariantFiltration exit code: $?"

# Annotate possible de novo mutations
outdir=/path/to/per_fam_DNM
mkdir -p ${outdir}/${famID}
outdir_per_fam=${outdir}/${famID}

gatk --java-options '-Xmx40g -Xms32g -XX:ParallelGCThreads=2' \
    VariantAnnotator \
    -R $ref_genome \
    -V $per_fam_vcf_tagged \
    -ped $per_fam_ped \
    -A PossibleDeNovo \
    -O ${outdir_per_fam}/possibleDNM_${famID}.vcf.gz

echo "[STATUS] VariantAnnotator exit code: $?"

# Convert annotated VCF to table
per_fam_vcf=${outdir_per_fam}/possibleDNM_${famID}.vcf.gz

outdir=/path/to/per_fam_DNM_filtered_table
mkdir -p ${outdir}/${famID}
outdir_per_fam=${outdir}/${famID}

gatk --java-options '-Xmx40g -Xms32g -XX:ParallelGCThreads=2' \
    VariantsToTable \
    -V $per_fam_vcf \
    -O ${outdir_per_fam}/possibleDNM_${famID}.table.txt \
    -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F TYPE -F hiConfDeNovo -F loConfDeNovo \
    -GF GT -GF DP -GF GQ

echo "[STATUS] VariantsToTable exit code: $?"

echo "[INFO] Job completed at $(date)"
