# QQperm Lambda Calculation and QQ Plot Script

# ---- Libraries ----
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyverse)
  library(QQperm)
  library(ggrepel)
  library(gridExtra)
})

# ---- Argument Parsing ----
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 4) stop("Usage: Rscript qqperm.R <input_varinfo> <input_casestatus> <vartype> <samplename>")

input_varinfo <- args[1]
input_casestatus <- args[2]
vartype <- args[3]
samplename <- args[4]

# ---- Configurable Paths ----
input_dir <- "INPUT_DIRECTORY" # Set to your input directory
output_dir <- "OUTPUT_DIRECTORY" # Set to your output directory
lambda_dir <- file.path(output_dir, "lambda")
qqplot_dir <- file.path(output_dir, "qqplot")
postgc_dir <- file.path(output_dir, "postGCdat")

# ---- File Names ----
qqperm_out <- file.path(output_dir, paste0(vartype, "_", samplename, "_qqperm_allgene.txt"))

# ---- Data Import ----
setwd(input_dir)
case_status <- read.table(input_casestatus, stringsAsFactors = FALSE, header = TRUE, sep = '\t')
case_status_mat <- as.matrix(case_status)

per_gene_sample <- read.table(input_varinfo, stringsAsFactors = FALSE, sep = '\t', header = TRUE)
per_gene_sample_mat <- as.matrix(t(per_gene_sample))

# ---- QQperm Results ----
setwd(output_dir)
qqperm_results <- read.table(qqperm_out, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

# ---- Data Consistency Checks ----
stopifnot(all(rownames(case_status_mat) == colnames(per_gene_sample_mat)))
stopifnot(nrow(qqperm_results) == nrow(per_gene_sample_mat))

# ---- Annotate Genes ----
qqperm_results$gene <- rownames(per_gene_sample_mat)

# ---- Case/Control Split ----
case_samples <- case_status_mat[case_status_mat[,1] == TRUE, , drop = FALSE]
ctrl_samples <- case_status_mat[case_status_mat[,1] == FALSE, , drop = FALSE]

# ---- Per-Gene Summaries ----
sum_gene_counts <- function(samples, per_gene_sample) {
  per_gene_sample %>%
    as.data.frame() %>%
    dplyr::filter(rownames(.) %in% rownames(samples)) %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "genes") %>%
    gather(key = "genes", value = "n_uniq_ppl", -genes) %>%
    select(genes, n_uniq_ppl) %>%
    arrange(desc(n_uniq_ppl))
}

gene_case <- sum_gene_counts(case_samples, per_gene_sample)
gene_ctrl <- sum_gene_counts(ctrl_samples, per_gene_sample)

n_total_case = nrow(case_samples)
n_total_ctrl = nrow(ctrl_samples)

gene_case_ctrl <- full_join(
  gene_case, gene_ctrl, by = 'genes', suffix = c('.case', '.ctrl')
) %>%
  mutate(
    a = ifelse(is.na(n_uniq_ppl.case), 0, n_uniq_ppl.case),
    b = ifelse(is.na(n_uniq_ppl.ctrl), 0, n_uniq_ppl.ctrl),
    c = n_total_case - a,    
    d = n_total_ctrl - b
  )

gene_case_ctrl$fisherp <- apply(gene_case_ctrl[,c('a','c','b','d')], 1, function(x) fisher.test(matrix(x, nrow = 2))$p.value)
gene_case_ctrl$ORfisher <- apply(gene_case_ctrl[,c('a','c','b','d')], 1, function(x) fisher.test(matrix(x, nrow = 2))$estimate)

gene_case_ctrl_annot <- gene_case_ctrl %>%
  left_join(qqperm_results, by = c('genes' = 'gene'))

gene_case_ctrl_annot$chi.p.o <- qchisq(gene_case_ctrl_annot$observed, 1, lower.tail = FALSE)
gene_case_ctrl_annot$chi.p.e <- qchisq(gene_case_ctrl_annot$perm, 1, lower.tail = FALSE)
gene_case_ctrl_annot$chi.fisherp <- qchisq(gene_case_ctrl_annot$fisherp, 1, lower.tail = FALSE)
gene_case_ctrl_annot$chi.est <- apply(gene_case_ctrl_annot[,c('a','c','b','d')], 1, function(x) chisq.test(matrix(x, nrow = 2))$statistic)

# ---- Output Annotated Data ----
out_name <- paste0(vartype, "_", samplename, "_qqperm_annot.csv")
setwd(lambda_dir)
write.csv(gene_case_ctrl_annot, out_name, row.names = FALSE)

# ---- QQ Plot Function ----
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(aes(x = expected, ymin = clower, ymax = cupper), alpha = 0.1) +
    geom_point(aes(expected, observed), shape = 1, size = 2, color = 'black') +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5, color = 'red') +
    geom_line(aes(expected, cupper), linetype = 2, size = 0.5, color = 'blue') +
    geom_line(aes(expected, clower), linetype = 2, size = 0.5, color = 'blue') +
    xlab(log10Pe) +
    ylab(log10Po)
}

# ---- Lambda Estimation ----
lambda_out <- data.frame(allgene = NA, geneset1 = NA, geneset2 = NA)
lamda_all <- estlambda2(p.o = gene_case_ctrl_annot$observed, p.e = gene_case_ctrl_annot$perm)$estimate
lambda_out$allgene <- lamda_all

gene_case_ctrl_filter1 <- gene_case_ctrl_annot %>%
  filter(ORfisher > 1, ORfisher != Inf)
lamda_filter1 <- estlambda2(p.o = gene_case_ctrl_filter1$observed, p.e = gene_case_ctrl_filter1$perm)$estimate
lambda_out$geneset1 <- lamda_filter1

gene_case_ctrl_filter2 <- gene_case_ctrl_annot %>%
  filter(ORfisher > 1, ORfisher != Inf, a > 1, b >= 1)
lamda_filter2 <- estlambda2(p.o = gene_case_ctrl_filter2$observed, p.e = gene_case_ctrl_filter2$perm)$estimate
lambda_out$geneset2 <- lamda_filter2
lambda_out$method <- 'estlambda2'

# ---- Alternative Lambda Estimation ----
lambda2_out <- data.frame(allgene = NA, geneset1 = NA, geneset2 = NA)
gene_case_ctrl_annot_clean <- gene_case_ctrl_annot %>%
  filter(chi.p.o != 0 & chi.p.e != 0, !is.na(chi.p.o) & !is.na(chi.p.e))
lambda2_out$allgene <- summary(lm(sort(gene_case_ctrl_annot_clean$chi.p.o) ~ 0 + sort(gene_case_ctrl_annot_clean$chi.p.e)))$coeff[1]

gene_case_ctrl_annot_clean2 <- gene_case_ctrl_annot_clean %>%
  filter(ORfisher > 1, ORfisher != Inf)
lambda2_out$geneset1 <- summary(lm(sort(gene_case_ctrl_annot_clean2$chi.p.o) ~ 0 + sort(gene_case_ctrl_annot_clean2$chi.p.e)))$coeff[1]

gene_case_ctrl_annot_clean3 <- gene_case_ctrl_annot_clean %>%
  filter(ORfisher > 1, ORfisher != Inf, a > 1, b >= 1)
lambda2_out$geneset2 <- summary(lm(sort(gene_case_ctrl_annot_clean3$chi.p.o) ~ 0 + sort(gene_case_ctrl_annot_clean3$chi.p.e)))$coeff[1]
lambda2_out$method <- 'estlamda2-unfilter'

lambda_final <- rbind(lambda_out, lambda2_out)
lambda_final$dataset <- paste0(vartype, "_", samplename, "_lambada.csv")

setwd(lambda_dir)
write.csv(lambda_final, lambda_final$dataset, row.names = FALSE)

# ---- Post-GC Correction and Output ----
GC_lam_allgene <- lambda_final[1, 'allgene']
GC_lam_filter1 <- lambda_final[1, 'geneset1']

gene_case_ctrl_annot <- gene_case_ctrl_annot %>%
  rowwise() %>%
  mutate(
    fisherp.GC = pchisq((qchisq(fisherp, df = 1, lower.tail = FALSE) / GC_lam_allgene), df = 1, lower.tail = FALSE),
    fisherp.GC.filter1 = pchisq((qchisq(fisherp, df = 1, lower.tail = FALSE) / GC_lam_filter1), df = 1, lower.tail = FALSE)
  )

gene_case_ctrl_filter1 <- gene_case_ctrl_filter1 %>%
  rowwise() %>%
  mutate(
    fisherp.GC = pchisq((qchisq(fisherp, df = 1, lower.tail = FALSE) / GC_lam_filter1), df = 1, lower.tail = FALSE)
  )

allvar_GC_PA <- gene_case_ctrl_annot %>%
  mutate(pfinal.GC = fisherp.GC.filter1, SYMBOL = genes, gene = genes)
allinone_PA_filter <- gene_case_ctrl_filter1 %>%
  mutate(pfinal.GC = fisherp.GC, SYMBOL = genes, gene = genes)

sig_cutoff_PA <- 0.05 / length(unique(allinone_PA_filter$gene))

# ---- Output Post-GC Data ----
allvar_GC_PA_update <- allvar_GC_PA %>%
  mutate(pfinal.GC = ifelse(pfinal.GC == 0, 1E-200, pfinal.GC)) %>%
  filter(!(ORfisher <= 1 | ORfisher == Inf & pfinal.GC < sig_cutoff_PA)) %>%
  arrange(pfinal.GC)

datGC_out_name <- paste0(vartype, "_", samplename, "_post_GCgenesetctrl.csv")
setwd(postgc_dir)
write.csv(allvar_GC_PA_update, datGC_out_name, row.names = FALSE)

# ---- QQ Plot ----
pts_pfinal.GC_PA <- gg_qqplot(allvar_GC_PA_update$pfinal.GC)

# ---- Save Plot ----
plot_out_name <- paste0(vartype, "_", samplename, "_qqplot.pdf")
setwd(qqplot_dir)
pdf(plot_out_name, width = 12, height = 6)
print(pts_pfinal.GC_PA)
dev.off()

# ---- End of Script ----
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 2, color='black') +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5, color='red') +
    geom_line(aes(expected, cupper), linetype = 2, size = 0.5, color='blue') +
    geom_line(aes(expected, clower), linetype = 2, size = 0.5, color='blue') +
    xlab(log10Pe) +
    ylab(log10Po)
}


print('plot function obtained')

hscr_knowngene = c()

GRN_gene = c()

#check GRN/known genes

t = PA_gene_case_ctrl_allgene %>%
  dplyr::filter(genes %in% c(GRN_gene, hscr_knowngene))


GC_lam_allgene=lambda_final[1,'allgene']
GC_lam_filter1=lambda_final[1,'geneset1']

PA_gene_case_ctrl_allgene = PA_gene_case_ctrl_allgene %>%
  rowwise() %>%
  dplyr::mutate(
    fisherp.GC = pchisq((qchisq(fisherp, df=1, lower.tail = F)/GC_lam_allgene), df=1, lower.tail = F),
    fisherp.GC.filter1  = pchisq((qchisq(fisherp, df=1, lower.tail = F)/GC_lam_filter1), df=1, lower.tail = F)
  ) 


PA_gene_case_ctrl_filter1 = PA_gene_case_ctrl_filter1 %>%
  rowwise() %>%
  dplyr::mutate(
    fisherp.GC = pchisq((qchisq(fisherp, df=1, lower.tail = F)/GC_lam_filter1), df=1, lower.tail = F)
  ) 



allvar_GC_PA = PA_gene_case_ctrl_allgene %>% dplyr::mutate(pfinal.GC = fisherp.GC.filter1, SYMBOL = genes, gene = genes)
allinone_PA_filter = PA_gene_case_ctrl_filter1 %>% dplyr::mutate(pfinal.GC = fisherp.GC, SYMBOL = genes, gene = genes)


#use the geneset cutoff
sig_cutoff_PA = 0.05/( length(unique(allinone_PA_filter$gene)))

#get list of sig genes
allvar_GC_PA %>% filter(pfinal.GC <= sig_cutoff_PA) %>% dplyr::select(SYMBOL, pfinal.GC, a,b,c,d, ORfisher)

min(allvar_GC_PA$pfinal.GC)


#make change for plotting: if p.logitfirth.GC==0, use  pfinal.GC
allvar_GC_PA_update = allvar_GC_PA %>%
  dplyr::mutate(
    pfinal.GC = ifelse(pfinal.GC==0, 1E-200, pfinal.GC)
  )

allvar_GC_PA_update %>% filter(pfinal.GC <= sig_cutoff_PA) %>% dplyr::select(SYMBOL, pfinal.GC, a,b,c,d, ORfisher)


#plot after GC & annotate important genes 
library(ggrepel)
library(gridExtra)
#annotation:
#GRN & known genes - blue dots

#remove the top genes as ORfisher <1 & sig
dat_PA_to_rm = allvar_GC_PA_update %>%
  dplyr::filter(ORfisher <= 1 | ORfisher == Inf) %>%
  dplyr::filter(pfinal.GC < sig_cutoff_PA)


allvar_GC_PA_update_filter = allvar_GC_PA_update %>%
  dplyr::filter(!SYMBOL %in% dat_PA_to_rm$SYMBOL) %>%
  arrange(pfinal.GC)

datGC_out_name = paste0(vartype,"_", samplename, '_', 'post_GCgenesetctrl.csv')

setwd('/postGCdat/')
write.csv(allvar_GC_PA_update_filter, datGC_out_name, row.names = F)


allvar_GC_PA_update = allvar_GC_PA_update_filter

print('final dataset to plot obtained')

#use pfinal.GC
allvar_GC_PA_update$pfinal.GC[allvar_GC_PA_update$fisherp==1]=1
pts_pfinal.GC_PA = gg_qqplot(allvar_GC_PA_update$pfinal.GC)

print('basic plot done')

GRN_hscrknown = unique(c(hscr_knowngene, GRN_gene) )

dat_PA = allvar_GC_PA_update %>% arrange(desc(-log10(pfinal.GC)))
dat_PA$x = pts_pfinal.GC_PA$data$expected
dat_PA$y = pts_pfinal.GC_PA$data$observed


GRN_hscrknown_dat_PA = dat_PA %>% filter(SYMBOL %in% GRN_hscrknown)

sig_dat_PA = dat_PA %>% filter(pfinal.GC <= sig_cutoff_PA)

library(stringr)
p_PA = pts_pfinal.GC_PA +
  geom_point(data = GRN_hscrknown_dat_PA, aes(x=x, y=y), color = 'lightblue') +
  geom_point(data=sig_dat_PA, aes(x=x, y=y), color='red') +
  geom_text_repel(data=sig_dat_PA, aes(label=SYMBOL, x=x, y=y), 
                  size=3, box.padding = 1.3, point.padding = 0.3, segment.color='grey50', 
                  hjust='left', vjust='top', max.overlaps = 100, show.legend = F) +
  ggtitle(bquote("QQplot of Pathogenic variants, Post GC, " ~ lambda ~ "=" ~ .(  str_sub(as.character(GC_lam_filter1), 1, 4) ) ~ " ")) +
  #add signficant cutoff for the var & exome wide sig line
  geom_hline(yintercept = -log10(sig_cutoff_PA), color='blue', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05/(18000)), color='darkgray', linetype = 'dashed') +
  theme_bw()

plot_out_name = paste0(vartype,"_", samplename, '_', 'qqplot.pdf')

setwd('qqplot/')
pdf(plot_out_name, width = 12, height = 6)
print(p_PA)
dev.off()






