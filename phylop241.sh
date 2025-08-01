#!/bin/bash

###############################################################################
# Script: phylop241.sh
# Description: Extracts phyloP scores from a bigWig file for specified regions.
# Usage: Submit as a SLURM array job. Update input/output/reference paths as needed.
# Requirements: UCSC utils (bigWigToBedGraph), SLURM
###############################################################################

#SBATCH --job-name=phylop
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24gb
#SBATCH --time=4:00:00
#SBATCH --output=log/phylop_%j.log
#SBATCH --partition=a100_dev,a100_short,a100_long
#SBATCH --array=1-N  # Set N to the number of regions

echo "[INFO] Job started at $(date)"

# Set table name and path (update as needed)
tablename=input_table
table=/path/to/region_sublist/${tablename}

# Parse region information for this array task
line="$(head -n $SLURM_ARRAY_TASK_ID $table | tail -n 1)"
region_id="$(printf "%s" "${line}" | cut -f1)"
chr="$(printf "%s" "${line}" | cut -f2)"
startpos="$(printf "%s" "${line}" | cut -f3)"
endpos="$(printf "%s" "${line}" | cut -f4)"

echo "[INFO] Processing array index: $SLURM_ARRAY_TASK_ID region: $region_id"

# Load UCSC utilities
module load ucscutils/374

# Reference bigWig file (update path as needed)
refbw=/path/to/phyloP_scores.bw

# Output directory and file
outdir=/path/to/phylop_output
outname=$outdir/${region_id}.bedGraph

mkdir -p $outdir

# Extract phyloP scores for the specified region
bigWigToBedGraph $refbw $outname -chrom=$chr -start=$startpos -end=$endpos

echo "[STATUS] bigWigToBedGraph exit code: $?"
echo "[INFO] Job completed at $(date)"

