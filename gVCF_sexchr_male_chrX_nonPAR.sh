#!/bin/bash

###############################################################################
# Script: gVCF_sexchr_male_chrX_nonPAR.sh
# Description: Runs GATK HaplotypeCaller on male samples for chrX non-PAR regions.
#              Designed for SLURM array jobs. Each task processes one sample.
# Usage: Submit as a SLURM array job. Reads sample info from a text file.
# Requirements: GATK, SLURM
###############################################################################

#SBATCH --job-name=gVCF
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=12:00:00
#SBATCH --output=log/gVCF_%j.log
#SBATCH --array=1-N  # Set N to the number of samples
#SBATCH --partition=cpu_short,cpu_medium,cpu_long

echo "[INFO] Job started at $(date)"

module purge
module load gatk/4.1.2.0

# Input sample list (tab-delimited: sample_id, ..., bam_path)
sample_table=/path/to/male_samples.txt
sample_line="$(head -n $SLURM_ARRAY_TASK_ID $sample_table | tail -n 1)"
sample_id="$(printf "%s" "${sample_line}" | cut -f1)"
bam_path="$(printf "%s" "${sample_line}" | cut -f5)"

echo "[INFO] Processing array index: $SLURM_ARRAY_TASK_ID sample: $sample_id"

# Reference genome and target intervals (update paths as needed)
ref_genome=/path/to/reference_genome.fasta
interval_bed=/path/to/chrX_nonPAR.interval.bed

# Output directories
out_dir=/path/to/gVCF_output
out_chr_dir=$out_dir/chrX_nonPAR
bamout_dir=$out_dir/bamout

mkdir -p $out_chr_dir $bamout_dir

gatk --java-options '-Xmx15g -Xms12g' HaplotypeCaller \
    --genotyping-mode DISCOVERY \
    -R $ref_genome \
    -I $bam_path \
    -O $out_chr_dir/${sample_id}.chrX_nonPAR.g.vcf \
    -ERC GVCF \
    --intervals $interval_bed \
    --interval-padding 10 \
    --bam-output $bamout_dir/${sample_id}.chrX_nonPAR.gvcf.out.bam \
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
    -G StandardHCAnnotation \
    --sample-ploidy 1

echo "[STATUS] GATK-HaplotypeCaller exit code: $?"

echo "[INFO] Job completed at $(date)"
    -O $out_chr/${samplename}.${union_chr}.g.vcf \
    -ERC GVCF \
    --intervals $union_region \
    --interval-padding 10 \
    --bam-output $out_dir/bamout/${samplename}.${union_chr}.gvcf.out.bam \
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
    -G StandardHCAnnotation \
    --sample-ploidy 1


echo _ESTATUS_ [ GATK-HaplotypeCaller for $samplename for $union_chr ]: $?


echo _END_  $(date)



