#!/bin/bash

###############################################################################
# Script: missingness_filter.sh
# Description: Filters variants with more than 10% missingness using bcftools.
#              Optionally verifies missingness with vcftools.
# Usage: Update input/output VCF paths as needed.
# Requirements: bcftools, vcftools
###############################################################################

# Filter variants with less than 10% missingness
bcftools view -i 'F_MISSING < 0.1' input_data.vcf.gz -Oz -o filtered_output.vcf.gz

# (Optional) Double-check missingness using vcftools
# Replace the following variables with actual file paths as needed:
#   SUBSET_VCF: path to filtered VCF
#   OUT: prefix for output files

vcftools --gzvcf $SUBSET_VCF --missing-site --out $OUT

# Check that the 'F_MISS' column (column 6) in the output is less than 0.1 for all sites
# This confirms that the data is correctly filtered for missingness
