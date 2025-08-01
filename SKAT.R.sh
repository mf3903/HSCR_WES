#!/bin/bash 
#SBATCH --job-name=skat # Job name 
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL) 
#SBATCH --mail-user=
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Run on a single node 
#SBATCH --cpus-per-task=4
#SBATCH --mem=48gb
#SBATCH --time=4:00:00
#SBATCH --output=log/skat_skato_%j.log
#SBATCH --partition=a100_dev,a100_short,a100_long
#SBATCH --array=1-xxx

module load r/4.1.2

echo _START_ $(date)

# ---- Set Input/Output Directories ----
INPUT_LIST_DIR="input_list_dir" # Set to your input list directory
TABLE_FILE="${INPUT_LIST_DIR}/perGene.skat.final.segmentaf"
OUTDIR="skat_out_dir" # Set to your output directory

# ---- Get Sample Name ----
line="$(head -n $SLURM_ARRAY_TASK_ID $TABLE_FILE | tail -n 1)"
samplename="$(basename "${line}" | cut -d'.' -f1)"

echo "Processing array index: $SLURM_ARRAY_TASK_ID for $samplename"

outfile="${OUTDIR}/${samplename}.skatbinary.skatO.pval.txt"

# ---- Run R Script ----
RSCRIPT_PATH="SKAT.R" # Set to your R script path
Rscript "$RSCRIPT_PATH" "${line}" "${outfile}"

echo _ESTATUS_ [ rscript ]: $?
echo _END_  $(date)




