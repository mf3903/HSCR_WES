#!/bin/bash

###############################################################################
# Script: gnomadv4_AF_by_population.sh
# Description: Extracts allele frequency and annotation data by population from
#              gnomAD v4 VCF files for specified chromosomes/regions.
# Usage: Submit as a SLURM array job. Update input/output/reference paths as needed.
# Requirements: bcftools, R, SLURM
###############################################################################

#SBATCH --job-name=gnomadv4
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=72gb
#SBATCH --time=12:00:00
#SBATCH --output=log/dwld_subextract_%j.log
#SBATCH --partition=a100_short,a100_long
#SBATCH --array=1-N  # Set N to the number of regions/chromosomes

echo "[INFO] Job started at $(date)"

module load bz2/1.0.6
module load gsl/2.5
module load bcftools/1.15.1

# Input chromosome/region list (one per line)
table=chrom_list.txt
line="$(head -n $SLURM_ARRAY_TASK_ID $table | tail -n 1)"
region_name=$line

echo "[INFO] Processing array index: $SLURM_ARRAY_TASK_ID for $region_name"

# Set input/output directories (update as needed)
filedir=/path/to/gnomad_v4_vcf
file=$filedir/gnomad.genomes.v4.0.sites.${region_name}.vcf.bgz

outdir=/path/to/gnomad_AF_by_population_output

# Extract desired fields using bcftools
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/allele_type\t%INFO/phylop\t%INFO/revel_max\t%INFO/spliceai_ds_max\t%INFO/AF_joint\t%INFO/AF_joint_nfe\t%INFO/AF_joint_eas\t%INFO/AF_joint_afr\t%INFO/AF_joint_amr\t%INFO/vep\n' "${file}" > ${region_name}.annot.AF.population.extract.txt

mv ${region_name}.annot.AF.population.extract.txt $outdir/${region_name}.annot.AF.population.extract.txt

echo "[STATUS] bcftools exit code: $?"

# Filter for coding and synonymous variants
var_dat=$outdir/${region_name}.annot.AF.population.extract.txt
grep -E 'frameshift_variant|inframe_deletion|inframe_insertion|missense_variant|splice_acceptor_variant|splice_donor_variant|stop_gained|synonymous_variant' $var_dat > $outdir/${region_name}.annot.AF.population.codingvar.txt

# Filter for high-quality PASS variants
grep 'PASS' $outdir/${region_name}.annot.AF.population.codingvar.txt > $outdir/${region_name}.annot.AF.population.codingvar.PASS.txt

echo "[STATUS] grep exit code: $?"



