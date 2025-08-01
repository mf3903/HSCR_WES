#!/bin/bash

###############################################################################
# Script: gVCF.sh
# Description: Runs GATK HaplotypeCaller to generate gVCF files from BAM inputs.
#              Designed for SLURM array jobs. Each task processes one sample.
# Usage: Submit as a SLURM array job. Reads BAM file paths from a text file.
# Requirements: GATK, SLURM
###############################################################################

#SBATCH --job-name=gVCF
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=99:00:00
#SBATCH --output=logs-gVCF/gVCF_%j.log
#SBATCH --array=1-N  # Set N to the number of samples
#SBATCH --partition=cpu_short,cpu_medium,cpu_long

echo "[INFO] Job started at $(date)"

module purge
module load gatk/4.1.2.0

# Generate sample BAM list if not already present
ls /path/to/bam_dir/*mapped.sorted.bam > samples.bam-dedup-bqsr-mapped.txt

table=samples.bam-dedup-bqsr-mapped.txt
line="$(head -n $SLURM_ARRAY_TASK_ID $table | tail -n 1)"
sample_bam="$(basename $line)"
sample_id=$(basename $sample_bam .mapped.sorted.bam)

echo "[INFO] Processing array index: $SLURM_ARRAY_TASK_ID sample: $sample_id"

# Reference genome and target regions (update paths as needed)
ref_genome=/path/to/reference_genome.fasta
capture_region=/path/to/target_regions.list

out_dir=gVCF
mkdir -p $out_dir

gatk --java-options '-Xmx28g' HaplotypeCaller \
    --genotyping-mode DISCOVERY \
    -R $ref_genome \
    -I $line \
    -O $out_dir/${sample_id}.g.vcf \
    -ERC GVCF \
    --intervals $capture_region \
    --interval-padding 10 \
    --bam-output $out_dir/${sample_id}.realigned.g.vcf.out.bam \
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
    -G StandardHCAnnotation

echo "[STATUS] GATK-HaplotypeCaller exit code: $?"

echo "[INFO] Job completed at $(date)"
    --interval-padding 10 \
    --bam-output $out_dir/${samplename}.realigned.g.vcf.out.bam \
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
    -G StandardHCAnnotation 


echo _ESTATUS_ [ GATK-HaplotypeCaller for $samplename ]: $?


echo _END_  $(date)


