args <- commandArgs()
print(args)

# Import required libraries
library(logistf)
library(dplyr)
library(tidyverse)

# Load input data (update file paths as needed)
input_data = read.csv('/path/to/input_data.csv', stringsAsFactors = FALSE)
covar_data = read.csv('/path/to/covariates.csv', stringsAsFactors = FALSE)
print('import done')

input_data_sub = input_data
gene_list = unique(input_data_sub$SYMBOL)
gene_coef = data.frame(genename=gene_list, ln_coef=NA, ln_lower0.95=NA, ln_upper0.95=NA, p=NA)

for(i in 1:length(gene_list)) {
  genename = gene_list[i]
  gene = input_data_sub %>%
    filter(SYMBOL == genename) %>%
    distinct()
  gene_covar = left_join(
    covar_data,
    gene[, c('IID', "n_PA")],
    by = 'IID'
  ) %>%
    rowwise() %>%
    mutate(
      PA_clean = ifelse(is.na(n_PA), 0, n_PA)
    ) %>%
    select(!c(n_PA)) %>%
    distinct()
  fit_gene = logistf(
    data = gene_covar,
    aff2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + sex + PA_clean,
    firth = TRUE,
    pl = TRUE,
    maxit = 1000,
    control = logistf.control(maxit = 1000)
  )
  result_log = cbind(lnOR = coef(fit_gene), confint(fit_gene), pval = fit_gene$prob)
  result_log_PA = result_log[rownames(result_log) == 'PA_clean']
  gene_coef$ln_coef[gene_coef$genename == genename] = result_log_PA[1]
  gene_coef$ln_lower0.95[gene_coef$genename == genename] = result_log_PA[2]
  gene_coef$ln_upper0.95[gene_coef$genename == genename] = result_log_PA[3]
  gene_coef$p[gene_coef$genename == genename] = result_log_PA[4]
  print(i)
}

print('logit_firth done')

setwd('/path/to/output_dir/')
write.csv(gene_coef, 'gene_coef2_PA.csv')
print('output done')


