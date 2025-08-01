###Single-cell ENS human embryonic gut annotation ###

# Load required libraries
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(dplyr)
library(patchwork)
library(ggplot2)
library(AnnotationHub)
library(umap)

# Define directories
input_dir <- "/path/to/input"
output_dir <- "/path/to/output"

# Load data
data <- readRDS(file.path(input_dir, "human_ENS_scRNA.rds"))
meta_data <- data@meta.data
unique(meta_data$annotation)

# Extract normalized and scaled data
normalized_data <- data@assays$RNA@data
scaled_data <- data@assays$RNA@scale.data

# Scree plot
ElbowPlot(data)

# Match clusters to annotations
levels(data)
head(Idents(data), 1)
table(Idents(data))

# Rename clusters with neural annotations
new_cluster_ids <- c(
  'ENCC/glia Progenitor', 'cycling ENCC/glia', 'Neuroblast', 'cycling neuroblast',
  'Undefined', 'Branch A1 (iMN)', 'Branch A2 (IPAN/IN)', 'Branch A3 (IPAN/IN)',
  'Branch A4 (IN)', 'Branch B1 (eMN)', 'Branch B2 (eMN)', 'Branch B3 (IPAN)',
  'Differentiating glia', 'Undefined', 'Glia 1 (DHH+)', 'Glia 2 (ELN+)', 'Glia 3 (BCAN+)'
)
names(new_cluster_ids) <- levels(data)
data <- RenameIdents(data, new_cluster_ids)

# Update annotations
data$annotation <- Idents(data)

# Plot UMAP
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5)

# Fetch expression data for genes of interest
genes_of_interest = NA
gene_data <- FetchData(data, vars = genes_of_interest)

# Calculate percentage of cells expressing each gene in each cluster
results <- lapply(colnames(gene_data), function(gene) {
  gene_data$expressed <- gene_data[, gene] > 0
  tapply(gene_data$expressed, data$seurat_clusters, mean) * 100
})

# Convert results to a data frame
results_df <- as.data.frame(results)
results_df$cluster_annotation <- new_cluster_ids
results_df <- results_df %>%
  dplyr::select(cluster_annotation, everything())

# Transpose and filter results
results_df_t <- t(results_df)
colnames(results_df_t) <- results_df_t[1, ]
results_df_t <- data.frame(results_df_t[-1, ]) %>%
  dplyr::mutate_all(as.numeric)

results_df_t_filtered <- results_df_t %>%
  tibble::rownames_to_column(var = "gene_name") %>%
  rowwise() %>%
  mutate(max_pct = max(c_across(-gene_name), na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(gene_name, max_pct, everything())

# Save filtered results
write.csv(results_df_t_filtered, file.path(output_dir, 'gene_expression_by_cluster.csv'), row.names = FALSE)


