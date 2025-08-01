#!/bin/bash

###############################################################################
# Script: plink_admixture.sh
# Description: Runs ADMIXTURE for a range of K values and collects cross-validation errors.
#              Also provides commands for extracting sample IDs and plotting ancestry.
# Usage: Update input/output paths as needed.
# Requirements: ADMIXTURE, awk, R, SLURM (if running on a cluster)
###############################################################################

# Reference: https://speciationgenomics.github.io/ADMIXTURE/

# Run ADMIXTURE for K=2 to K=5 and collect cross-validation errors
for K in {2..5}; do
    admixture --cv /path/to/input.bed $K > log${K}.out
done

# Collect CV error for each K
grep "CV" log*.out | awk '{print $3, $4}' | cut -c 4,7-20 > cv_error_summary.txt

# Extract filename/sampleID from .nosex files
awk '{split($1, name, "."); print $1, name[2]}' /path/to/*.nosex > sample_id_list.txt

# Note:
# - The .nosex file contains sample names corresponding to the *Q and *P files.
# - Manually annotate population groups as needed for downstream analysis.

# Plot individual ancestry using R (autoplot)
Rscript plotADMIXTURE.r -p <prefix> -i <sample_list.txt> -k <max_K> -l <population_labels>

# Example:
# Rscript plotADMIXTURE.r -p probands_postprune -i probands.nosex.list.cleaned.txt -k 5 -l Chaklab,AFR,AMR,EAS,EUR
