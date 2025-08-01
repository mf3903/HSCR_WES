#!/bin/sh

###############################################################################
# Script: metaRNN.sh
# Description: Runs MetaRNN for variant annotation using TensorFlow.
# Usage: Submit as a SLURM job. Update input/output paths as needed.
# Requirements: MetaRNN, TensorFlow, Conda, SLURM
###############################################################################

#SBATCH --job-name=metaRNN
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=72gb
#SBATCH --time=72:00:00
#SBATCH --output=log/metaRNN_%j.log
#SBATCH --partition=gpu4_medium,gpu4_long,gpu8_medium,gpu8_long

echo "[INFO] Job started at $(date)"

module load condaenvs/gpu/tensorflow2.9
conda activate MetaRNN

# Run MetaRNN annotation (update input VCF path as needed)
python ./MetaRNN.py hg38 INDELs_to_annot.vcf

echo "[STATUS] MetaRNN exit code: $?"
echo "[INFO] Job completed at $(date)"

