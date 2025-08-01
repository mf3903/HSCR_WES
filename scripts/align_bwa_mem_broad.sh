#!/bin/bash

###############################################################################
# Script: align_bwa_mem.sh
# Description: Batch alignment of paired-end FASTQ files using BWA-MEM.
#              Converts SAM to sorted, indexed BAM and generates coverage bedgraph.
# Usage: Designed for SLURM array jobs. Reads sample info from a CSV table.
# Requirements: BWA, Samtools, Bedtools, SLURM
###############################################################################

#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=30:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=bwamem_array
#SBATCH --mail-type=END,FAIL
#SBATCH --output=logs-align-bwa-mem/bwa_mem_%j.log
#SBATCH --array=1-N  # Set N to the number of samples

# Load required modules
module add bwa/0.7.17
module add samtools
module load bedtools/2.27.1

# Input sample table (CSV: sample_id,fastq_R1,fastq_R2)
table=samples.fastq-trim-trimmomatic.csv

# Parse sample information for this array task
line="$(head -n $SLURM_ARRAY_TASK_ID $table | tail -n 1)"
sample_id="$(printf "%s" "${line}" | cut -f1 -d ",")"
fastq_R1="$(printf "%s" "${line}" | cut -f2 -d ",")"
fastq_R2="$(printf "%s" "${line}" | cut -f3 -d ",")"

echo "[INFO] Processing sample: $sample_id"

# Reference genome (update path as needed)
ref_bwa=/path/to/reference_genome.fasta

# Output directories
out_dir=BAM-BWA
out_bam_dir=$out_dir/sorted_indexed_bam
out_bedgraph_dir=$out_dir/bedgraph

mkdir -p $out_dir $out_bam_dir $out_bedgraph_dir

# Construct read group
readGroup="@RG\tID:${sample_id}\tSM:${sample_id}\tLB:${sample_id}\tPL:ILLUMINA"
echo "[INFO] Read group: $readGroup"

# Alignment with BWA-MEM
bwa mem -K 100000000 -M -Y -v 2 -t $SLURM_CPUS_PER_TASK -R "${readGroup}" $ref_bwa $fastq_R1 $fastq_R2 > $out_dir/${sample_id}_aligned_reads.sam
echo "[STATUS] bwa mem exit code: $?"

# Convert SAM to BAM, sort, index, and generate bedgraph
samtools view -bS -h -o $out_dir/${sample_id}_aligned_reads.bam $out_dir/${sample_id}_aligned_reads.sam
samtools sort -o $out_bam_dir/${sample_id}_aligned_reads_sorted.bam $out_dir/${sample_id}_aligned_reads.bam
samtools index $out_bam_dir/${sample_id}_aligned_reads_sorted.bam
samtools view -b $out_bam_dir/${sample_id}_aligned_reads_sorted.bam | genomeCoverageBed -ibam stdin -bg > $out_bedgraph_dir/${sample_id}_aligned_reads.bedgraph

echo "[STATUS] samtools pipeline exit code: $?"

echo "[INFO] Alignment pipeline completed at $(date)"

