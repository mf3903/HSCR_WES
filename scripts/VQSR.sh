#!/bin/bash

###############################################################################
# Script: VQSR.sh
# Description: Performs Variant Quality Score Recalibration (VQSR) on SNPs using GATK.
#              Applies recalibration to input VCF and outputs recalibrated VCF.
# Usage: Submit as a SLURM job. Update input/output/reference paths as needed.
# Requirements: GATK 4.x, R, SLURM
###############################################################################

#SBATCH --job-name=vqsr_snp
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=92gb
#SBATCH --time=5-00:00:00
#SBATCH --output=log/vqsr_snp_%j.log
#SBATCH --partition=cpu_medium

echo "[INFO] Job started at $(date)"

module purge
module load gatk/4.1.2.0
module load r/3.5.1

# Reference genome and input VCF (update paths as needed)
ref_genome=/path/to/reference_genome.fasta
indir=/path/to/input_vcf_dir

# GATK resource bundle (update paths as needed)
gatk_bundle_dir=/path/to/gatk_bundle
hapmap=$gatk_bundle_dir/hapmap.vcf.gz
kg_omni=$gatk_bundle_dir/omni.vcf.gz
kg_snps=$gatk_bundle_dir/1000G_snps.vcf.gz
dbsnp=$gatk_bundle_dir/dbsnp.vcf

# Output directory
outdir=/path/to/output_vqsr_dir

# Run VariantRecalibrator for SNPs
gatk --java-options '-Xmx72g -Xms72g -XX:ParallelGCThreads=4' \
    VariantRecalibrator \
    -R $ref_genome \
    -V $indir/input_jointcalls.vcf.gz \
    -mode SNP \
    -O $outdir/output.vqsr.snp.recal \
    --tranches-file $outdir/output.vqsr.snp.tranches \
    --rscript-file $outdir/output.vqsr.snp.plots.R \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 $kg_omni \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 $kg_snps \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -tranche 100.0 \
    -tranche 99.95 \
    -tranche 99.9 \
    -tranche 99.8 \
    -tranche 99.5 \
    -tranche 99.0 \
    -tranche 95.0 \
    -tranche 91.0 \
    -tranche 90.0

echo "[STATUS] VariantRecalibrator exit code: $?"

# Apply VQSR to SNPs
gatk --java-options '-Xmx240g -Xms200g -XX:ParallelGCThreads=4' \
    ApplyVQSR \
    -R $ref_genome \
    -V $indir/input_jointcalls.vcf.gz \
    -mode SNP \
    --truth-sensitivity-filter-level 99.8 \
    --recal-file $outdir/output.vqsr.snp.recal \
    --tranches-file $outdir/output.vqsr.snp.tranches \
    -O $outdir/output.snp.recalibrated_99.8.vcf.gz \
    --create-output-variant-index true

echo "[STATUS] ApplyVQSR exit code: $?"

echo "[INFO] Job completed at $(date)"
    ApplyVQSR \
    -R $ref_genome \
    -V $indir/AUTOSOMAL.jointcalls.vcf.gz \
    -mode SNP \
    --truth-sensitivity-filter-level 99.8 \
    --recal-file $outdir/AUTOSOMAL.vqsr.snp.recal \
    --tranches-file $outdir/AUTOSOMAL.vqsr.snp.tranches \
    -O $outdir/SNP.recalibrated_99.8.vcf.gz \
    --create-output-variant-index true 


echo _ESTATUS_ [ vqsr_p2 for $chr]: $?


echo _END_  $(date)
