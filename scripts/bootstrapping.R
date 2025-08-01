# Bootstrapping script for per-gene variant counts

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 4) stop("Usage: Rscript bootstrapping.R <ctrl_all_file> <case_all_file> <ctrl_ctrl_pool_file> <out_dir>")

ctrl_all_file <- args[1]
case_all_file <- args[2]
ctrl_ctrl_pool_file <- args[3]
out_dir <- args[4]

library(stringr)
library(tidyverse)

# ---- Read Input Files ----
ctrl_all <- read.csv(ctrl_all_file, stringsAsFactors = FALSE)
case_all <- read.csv(case_all_file, stringsAsFactors = FALSE)
ctrl_ctrl_pool <- read.csv(ctrl_ctrl_pool_file, stringsAsFactors = FALSE)

# ---- Clean Control Sample IDs ----
ctrl_ids <- unique(str_trim(ctrl_ctrl_pool$ctrl_ID))
ctrl_ids <- ctrl_ids[ctrl_ids != "" & !is.na(ctrl_ids)]

# ---- Clean All Sample IDs ----
all_samples <- unique(str_trim(unlist(strsplit(toString(ctrl_all$allsample), ","))))
all_samples <- all_samples[all_samples != "" & !is.na(all_samples)]

# ---- Filter Control IDs to Only Those Present in All Samples ----
ctrl_ppl_uniq_clean <- intersect(ctrl_ids, all_samples)

# ---- Prepare Long Format Data ----
ctrl_all_long_clean <- ctrl_all %>%
  separate_rows(allsample, sep = ",", convert = FALSE) %>%
  filter(allsample %in% ctrl_ppl_uniq_clean) %>%
  distinct()

case_all_long <- case_all %>%
  separate_rows(allsample, sep = ",", convert = FALSE) %>%
  distinct()

n_sample_tot <- length(unique(case_all_long$allsample))

# ---- Per-Gene Counts for Cases ----
case_ct <- case_all_long %>%
  group_by(SYMBOL, var_type) %>%
  summarise(a = n_distinct(allsample), .groups = "drop") %>%
  mutate(c = n_sample_tot - a) %>%
  distinct()

# ---- Bootstrapping ----
n_iter <- 10000 # Number of bootstrap iterations

message("Bootstrapping iterations start")
for (i in seq_len(n_iter)) {
  # Sample controls
  sample_ctrls <- sample(ctrl_ppl_uniq_clean, n_sample_tot, replace = FALSE)
  ctrl_301 <- ctrl_all_long_clean %>%
    filter(allsample %in% sample_ctrls) %>%
    distinct()
  
  # Per-gene counts for controls
  ctrl_301_ct <- ctrl_301 %>%
    group_by(SYMBOL, var_type) %>%
    summarise(b = n_distinct(allsample), .groups = "drop") %>%
    mutate(d = n_sample_tot - b) %>%
    distinct()
  
  # Merge case and control counts
  case_vs_ctrl301 <- inner_join(
    case_ct,
    ctrl_301_ct,
    by = c('SYMBOL', 'var_type')
  ) %>%
    distinct()
  
  case_vs_ctrl301$iteration <- i
  
  # Ensure output directory exists
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  fname <- file.path(out_dir, paste0('a_b_c_d_perGene_', i, '.csv'))
  write.csv(case_vs_ctrl301, fname, row.names = FALSE)
  
  if(i %% 100 == 0) message("Iteration: ", i)
}
message("Bootstrapping iterations complete")
    ) %>%
    distinct()
  
  #put the two datsets together by gene & var (only think of the genes with >=1 ct in case & ctrl database)
  case_vs_ctrl = inner_join(
    case_ct,
    ctrl_301_ct, 
    by = c('SYMBOL', 'var_type')
  ) %>%
    distinct()

  case_vs_ctrl$iteration = i
  
  outdir=$out_dir
  
  fname=paste0(outdir, 'a_b_c_d_perGene_', i, '.csv')

  write.csv(case_vs_ctrl, fname, row.names = F)

  print(i)
}

 
print('iter done')
 
