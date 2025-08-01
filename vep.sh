#!/bin/bash

###############################################################################
# Script: vep.sh
# Description: Runs Ensembl VEP for functional annotation of variants in a VCF.
# Usage: Submit as a SLURM job. Update input/output/reference/plugin paths as needed.
# Requirements: VEP, Samtools, Tabix, SLURM
###############################################################################

#SBATCH --job-name=vep
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=150gb
#SBATCH --time=5-00:00:00
#SBATCH --output=log/vep_%j.log
#SBATCH --partition=fn_medium,fn_long

echo "[INFO] Job started at $(date)"

# Load required modules
module purge
module load vep/101
module load samtools
module load tabix

# Reference genome (update path as needed)
ref_genome=/path/to/reference_genome.fasta

# Input VCF file (update path as needed)
in_vcf=/path/to/input_variants.vcf.gz

# Plugin and cache directories (update paths as needed)
plugin_dir=/path/to/vep_plugins
plugin_file_dir=/path/to/vep_plugin_files
cache_dir=/path/to/vep_cache

# Custom annotation resources (update paths as needed)
phylop_100way=/path/to/phylop100way.bw
clinvar=/path/to/clinvar.vcf.gz
CADD1=$plugin_file_dir/whole_genome_SNVs.tsv.gz
CADD2=$plugin_file_dir/gnomad.genomes.indel.tsv.gz

# Output directory and file
out_dir=/path/to/vep_output
out_name=${out_dir}/vep_annotated_output.vcf

mkdir -p $out_dir

# Run VEP annotation
vep -i $in_vcf \
    --dir_plugins $plugin_dir \
    --dir_cache $cache_dir \
    --offline \
    --force_overwrite \
    --vcf \
    --refseq \
    --fasta $ref_genome \
    --use_transcript_ref \
    --output_file $out_name \
    --plugin AncestralAllele,$plugin_file_dir/homo_sapiens_ancestor_GRCh38.fa \
    --plugin CADD,$CADD1,$CADD2 \
    --sift b \
    --polyphen p \
    --canonical \
    --symbol \
    --custom $phylop_100way,phylop_100way,bigwig,exact,0 \
    --af_1kg \
    --af_gnomad \
    --distance 1000 

echo _ESTATUS_ [ vep for all_chr ]: $?
echo _END_  $(date)
