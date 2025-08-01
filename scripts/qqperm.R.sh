#!/bin/bash 
#SBATCH --job-name=qqperm_job # Job name 
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL) 
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Run on a single node 
#SBATCH --cpus-per-task=2 #this needs to align with the t (thread) needed by the software
#SBATCH --mem=72gb # Job memory request 
#SBATCH --time=64:00:00 # Time limit hrs:min:sec 
#SBATCH --output=log/QQperm_%j.log
#SBATCH --array=1
#SBATCH --partition=gpu4_medium,gpu4_long,gpu8_medium,gpu8_long


# ---- Load Modules ----
module purge  
module load r/4.0.3


# ---- Input Table (set this variable as needed) ----
INPUT_TABLE="input.syn.var1pct_final.table"


# ---- Parse Input Line ----
line="$(head -n $SLURM_ARRAY_TASK_ID $INPUT_TABLE | tail -n 1)"
filedat="$(printf "%s" "${line}" | cut -f2)"
filestatus="$(printf "%s" "${line}" | cut -f3)"
patype="$(printf "%s" "${line}" | cut -f1)"
samplename="SAMPLE_NAME" # Set this variable as needed


# ---- R Script to Run ----
R_SCRIPT="qqperm.R" # Set to your R script filename


# ---- Run R Script ----
Rscript "$R_SCRIPT" "$filedat" "$filestatus" "$patype" "$samplename"
 
echo "_ESTATUS_ [ rscript ]: $?"
echo "_END_  $(date)"

Rscript $rcode $filedat $filestatus $patype $samplename
 
echo _ESTATUS_ [ rscript ]: $?
echo _END_  $(date)


