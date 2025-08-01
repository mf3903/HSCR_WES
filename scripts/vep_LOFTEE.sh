#!/bin/bash

###############################################################################
# Script: vep_LOFTEE.sh
# Description: Runs Ensembl VEP with the LOFTEE plugin for loss-of-function annotation.
# Usage: Submit as a SLURM job. Update input/output/reference/plugin paths as needed.
# Requirements: VEP, LOFTEE, Samtools, Tabix, SLURM
###############################################################################

#SBATCH --job-name=vep_loftee
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96gb
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
in_vcf=/path/to/input_variants.vcf

# LOFTEE plugin and cache directories (update paths as needed)
loftee_dir=/path/to/loftee
cache_dir=/path/to/vep_cache

# Output directory and file
out_dir=/path/to/vep_loftee_output
out_name=${out_dir}/vep_annotated_LOFTEE.vcf

mkdir -p $out_dir

# Run VEP with LOFTEE plugin
vep -i $in_vcf \
    --dir_plugins $loftee_dir \
    --dir_cache $cache_dir \
    --offline \
    --force_overwrite \
    --vcf \
    --refseq \
    --fasta $ref_genome \
    --use_transcript_ref \
    --output_file $out_name \
    --plugin LoF,loftee_path:$loftee_dir,gerp_bigwig:$loftee_dir/gerp_conservation_scores.homo_sapiens.GRCh38.bw,human_ancestor_fa:$loftee_dir/human_ancestor.fa.gz,conservation_file:$loftee_dir/loftee.sql \
    --canonical \
    --symbol \
    --distance 1000

echo "[STATUS] VEP LOFTEE exit code: $?"
echo "[INFO] Job completed at $(date)"

