#!/bin/bash

###############################################################################
# Script: verifybam.sh
# Description: Runs VerifyBamID2 to estimate contamination in BAM files.
# Usage: Submit as a SLURM array job. Update input/output/reference paths as needed.
# Requirements: VerifyBamID2, SLURM
###############################################################################

#SBATCH --job-name=freemix
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=4:00:00
#SBATCH --output=logs-contamination/freemix_%j.log
#SBATCH --array=1-N  # Set N to the number of samples
#SBATCH --partition=cpu_dev,cpu_short,cpu_medium,cpu_long

echo "[INFO] Job started at $(date)"

module load verifybamid/2.0.1

home_dir=$PWD

# Input BAM list (one BAM path per line)
table=samples.bam-dedup-bqsr.txt
line="$(head -n $SLURM_ARRAY_TASK_ID $table | tail -n 1)"
sample_bam="$(basename $line)"
sample_id=$(basename $sample_bam .bam)

echo "[INFO] Processing array index: $SLURM_ARRAY_TASK_ID sample: $sample_id"

# Create output directory for this sample
out_dir=contamination/$sample_id
mkdir -p $out_dir
cd $out_dir

# Set VerifyBamID2 resource and reference paths (update as needed)
VERIFY_BAM_ID_HOME=/path/to/VerifyBamID
ref_genome=/path/to/reference_genome.fasta

${VERIFY_BAM_ID_HOME}/bin/VerifyBamID \
    --UDPath ${VERIFY_BAM_ID_HOME}/resource/1000g.phase3.100k.b38.vcf.gz.dat.UD \
    --BedPath ${VERIFY_BAM_ID_HOME}/resource/1000g.phase3.100k.b38.vcf.gz.dat.bed \
    --MeanPath ${VERIFY_BAM_ID_HOME}/resource/1000g.phase3.100k.b38.vcf.gz.dat.mu \
    --Reference $ref_genome \
    --BamFile ${home_dir}/$line

echo "[STATUS] VerifyBamID exit code: $?"

echo "[INFO] Job completed at $(date)"

