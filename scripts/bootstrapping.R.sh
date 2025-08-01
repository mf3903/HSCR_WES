#!/bin/bash 
#SBATCH --job-name=bootstrapping_job # Job name 
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL) 
#SBATCH --mail-user=mf3903@nyu.edu # Where to send mail 
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Run on a single node 
#SBATCH --cpus-per-task=4 #this needs to align with the t (thread) needed by the software
#SBATCH --mem=96gb # Job memory request 
#SBATCH --time=99:00:00 # Time limit hrs:min:sec 
#SBATCH --output=log/bootstrapping_%j.log
#SBATCH --partition=cpu_medium,cpu_long
 
# ---- Load Modules ----
module purge
module load r/4.0.3  

# ---- Set Script and Input Variables ----
R_SCRIPT="bootstrapping.R" # Set to your R script filename
CTRL_ALL_FILE="ctrl_all_file.csv"
CASE_ALL_FILE="case_all_file.csv"
CTRL_CTRL_POOL_FILE="ctrl_ctrl_pool_file.csv"
OUT_DIR="iteration_out"

# ---- Run Bootstrapping ----
Rscript "$R_SCRIPT" "$CTRL_ALL_FILE" "$CASE_ALL_FILE" "$CTRL_CTRL_POOL_FILE" "$OUT_DIR"
 
echo _ESTATUS_ [ rscript ]: $?
echo _END_  $(date)
