# ============================================================
# Spatial transcriptomics analysis using Seurat
# ------------------------------------------------------------
# Purpose:
#   - Load spatial transcriptomics data into a Seurat object
#   - Perform QC, normalization, clustering, and downstream analysis
#
# Reproducibility notes:
#   - Uses relative file paths
#   - All tunable parameters are defined below
#   - Raw data are not tracked in this repository
#
# Notes:
#   - Helper functions are defined in R/build_seurat_object.R
#   - Raw data are not included in this repository
#
# Author: Jessica Nicholson
# ============================================================


# ============================
# CONFIGURATION
# ============================

# Directories (relative to repo root)
DATA_DIR    <- "data"
RESULTS_DIR <- "results"
dir.create(file.path(DIR_FIGURES), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(DIR_TABLES), showWarnings = FALSE, recursive = TRUE)

# Analysis parameters
N_PCS      <- 20
RESOLUTION <- 0.6


# ============================
# LIBRARIES
# ============================
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(RImageJROI)
  library(dplyr)
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
  library(tidyverse)
  library(Matrix)
  })

# =========================
# HELPER SOURCES
# =========================
source("R/build_seurat_object.R")

# =========================
# =========================
# Define your samples
# =========================
SAMPLES <- tribble(
  ~sample_id, ~condition, ~sample_code, ~counts_path,                             ~roi_zip,
  "sample_untreated",   "untreated",      1,            "sample_untreated.csv",                         "sample_untreated_rois.zip",
  "sample_treated",   "treated",      2,            "sample_treated.csv",                         "sample_treated_rois.zip",
 )

# If some files are actually tab-separated, replace read_counts_csv inside build_sample
# with a small switch that uses fread(..., sep="\t") for those rows. For now, we assume comma CSV.

# =========================
# Build Seurat objects per sample
# =========================
objs <- purrr::pmap(SAMPLES, function(sample_id, condition, sample_code, counts_path, roi_zip) {
  message("Building ", sample_id, " …")
  build_sample(counts_path, roi_zip, sample_id, condition, sample_code)
})

names(objs) <- SAMPLES$sample_id

# =========================
# Merge and continue with pipeline
# =========================
seurat_combined <- Reduce(function(a, b) merge(a, y = b), objs)

# Define groups
seurat_combined$condition <- ifelse(
  grepl("untreated", seurat_combined$sample_id),
  "untreated",
  "treated"
)

# =========================
# Basic QC: filter cells based on feature and count thresholds
# =========================
# Ensure RNA is the default assay
DefaultAssay(seurat_combined) <- "RNA"

# ---- cells: keep cells with ≥10 detected genes ----
seurat_combined <- subset(seurat_combined, subset = nFeature_RNA >= 10)

# ---- genes: keep genes expressed in ≥10 cells (v5/v4 compatible) ----

# ---- QC BEFORE normalization ----
# 1) filter cells by detected genes
seurat_combined <- subset(seurat_combined, subset = nFeature_RNA >= 10)

# 2) filter genes expressed in ≥10 cells (use raw counts)
mat_counts <- get_raw_counts(seurat_combined, assay = "RNA")
gene_counts <- Matrix::rowSums(mat_counts > 0)
genes_to_keep <- names(gene_counts[gene_counts >= 10])
seurat_combined <- subset(seurat_combined, features = genes_to_keep)

# =========================
# Normalization
# =========================
seurat_combined <- NormalizeData(object = seurat_combined,
normalization.method = "LogNormalize",
scale.factor = 10000)

# Scale data
seurat_combined <- ScaleData(seurat_combined, features = rownames(seurat_combined))

# Run PCA
seurat_combined <- RunPCA(seurat_combined, features = rownames(seurat_combined))

# Run UMAP
seurat_combined <- FindNeighbors(seurat_combined, dims = 1:N_PCS)
seurat_combined <- FindClusters(seurat_combined, resolution = RESOLUTION)
seurat_combined <- RunUMAP(seurat_combined, dims = 1:N_PCS)

# =========================
# Separate by treatment
# =========================
seurat_combined$condition <- ifelse(
grepl("untreated", seurat_combined$sample_id),
"untreated",
"treated")
DimPlot(seurat_combined, group.by = "seurat_clusters", split.by = "condition", label = TRUE)
DimPlot(seurat_combined, group.by = "sample_id", split.by = "condition")

# =========================
# User-defined cluster labels
# =========================
combined_cluster_labels <- c("0" = "Neuro",
                             "1" = "Olig",
                             "2" = "Astro",
                             "3" = "GABA",
)
seurat_combined$cluster_named <- plyr::mapvalues(
  x = as.character(seurat_combined$seurat_clusters),
  from = names(combined_cluster_labels),
  to = combined_cluster_labels
)

# =========================
# Run DEG within each cluster for cluster gene expression
# =========================
Idents(seurat_combined) <- "seurat_clusters"
DefaultAssay(seurat_combined) <- "RNA"

# (Optional but safe for Seurat v5 objects)
for (a in Assays(seurat_combined)) {
  seurat_combined <- SeuratObject::JoinLayers(seurat_combined, assay = a)
}

clusters <- levels(Idents(seurat_combined))

deg_list <- map(clusters, function(clu) {
  obj_clu <- subset(seurat_combined, idents = clu)
  
  # If a cluster doesn't contain both conditions, skip it gracefully
  conds <- unique(as.character(obj_clu$condition))
  if (!all(c("treated", "untreated") %in% conds)) return(NULL)
  
  res <- FindMarkers(
    object          = obj_clu,
    ident.1         = "treated",
    ident.2         = "untreated",
    group.by        = "condition",
    test.use        = "wilcox",
    min.pct         = 0,
    logfc.threshold = 0
  )
  
  res %>%
    tibble::rownames_to_column("gene") %>%
    transmute(
      cluster = clu,
      gene,
      p_val,
      p_val_adj,
      avg_log2FC
    )
})

deg_all <- list_rbind(deg_list)

# Write single CSV
readr::write_csv(
  deg_all,
  file.path(DIR_TABLES, "DEG_treated_vs_untreated_all_clusters.csv")
)

# =========================
# Average expression per cluster
# =========================
avg_exp <- AverageExpression(seurat_combined, group.by = "seurat_clusters")$RNA 
readr::write_csv(
  avg_exp,
  file.path(DIR_TABLES, "Avg_Gene_Expression_Per_Cluster.csv") # Average expression
)

# =========================
# Cell counts per cluster
# =========================
# Subset to just the cells in the current cluster
seurat_subset <- subset(seurat_combined, idents = clust)
  
# Set Idents to 'condition' (treated vs untreated)
Idents(seurat_combined) <- "condition"

cell_counts <- table(
  factor(seurat_combined$seurat_clusters, levels = levels(Idents(seurat_combined))),
  seurat_combined$condition,
  seurat_combined$sample_id
)
cell_counts_df <- as.data.frame(cell_counts)
colnames(cell_counts_df) <- c("cluster", "condition", "sample_id", "cell_count")

readr::write_csv(
  cell_counts_df,
  file.path(DIR_TABLES, "Cell_Counts_Per_Cluster_By_Condition.csv", row.names = FALSE) # Cell counts
)

# =========================
# OPTIONAL: Cluster untreated separately for marker gene identification
# =========================
# Untreated
# seurat_untreated <- subset(seurat_combined, subset = condition == "untreated")
# seurat_untreated <- seurat_untreated |>
#  NormalizeData() |>
#  FindVariableFeatures() |>
#  ScaleData() |>
#  RunPCA() |>
#  FindNeighbors(dims = 1:N_PCS) |>
#  FindClusters(resolution = RESOLUTION) |>   # Adjust as needed
#  RunUMAP(dims = 1:N_PCS)

# Scale all genes
#   seurat_untreated <- ScaleData(seurat_untreated, features = rownames(seurat_untreated))
# PCA using all genes
#   seurat_untreated <- RunPCA(seurat_untreated, features = rownames(seurat_untreated))
# Find neighbors & clusters
#   seurat_untreated <- FindNeighbors(seurat_untreated, dims = 1:N_PCS)
#   seurat_untreated <- FindClusters(seurat_untreated, resolution = RESOLUTION)
# UMAP for visualization
#   seurat_untreated <- RunUMAP(seurat_untreated, dims = N_PCS)
# Optional: visualize UMAP
#   DimPlot(seurat_untreated, reduction = "umap", label = TRUE, group.by = "seurat_clusters") +
#   ggtitle("Untreated Sample Clusters")
# Identify top genes for untreated cluster identification
#   untreated_markers <- FindAllMarkers(
#   seurat_untreated,
#    only.pos = TRUE,
#   min.pct = 0.2,
#   logfc.threshold = 0.25
# )

#readr::write_csv(
#  untreated_markers,
#  file.path(DIR_TABLES, "Untreated_Cluster_Markers.csv", row.names = FALSE)
# )

# Treated
# Cluster treated & untreated separately
# Separate treated samples for marker identification
# seurat_treated <- subset(seurat_combined, subset = condition == "treated")
# Normalize (if not done yet)
# seurat_treated <- NormalizeData(seurat_treated)
# Skip variable feature selection
# Scale all genes
# seurat_treated <- ScaleData(seurat_treated, features = rownames(seurat_treated))
# PCA using all genes
# seurat_treated <- RunPCA(seurat_treated, features = rownames(seurat_treated))
# Find neighbors & clusters
# seurat_treated <- FindNeighbors(seurat_treated, dims = N_PCS)
# seurat_treated <- FindClusters(seurat_treated, resolution = RESOLUTION)
# UMAP for visualization
# seurat_treated <- RunUMAP(seurat_treated, dims = N_PCS)
# Optional: visualize UMAP
# DimPlot(seurat_treated, reduction = "umap", label = TRUE, group.by = "seurat_clusters") +
#  ggtitle("Treated Sample Clusters")
# Identify top genes for treated cluster identification
# treated_markers <- FindAllMarkers(
#  seurat_treated,
#  only.pos = TRUE,
#  min.pct = 0.2,
#  logfc.threshold = 0.25
#)

#readr::write_csv(
#  treated_markers,
#  file.path(DIR_TABLES, "Treated_Cluster_Markers.csv", row.names = FALSE)
#)

# =========================
# Anchoring for untreated label transfer to treated 
# - ONLY IF UNTREATED AND TREATED CLUSTERED SEPARATELY
# =========================
anchors <- FindTransferAnchors(
  reference = seurat_untreated,
  query = seurat_treated,
  dims = N_PCS,
  features = rownames(seurat_untreated)
)
# Transfer labels
treated_predictions <- TransferData(
  anchorset = anchors,
  refdata = seurat_untreated$seurat_clusters,
  dims = N_PCS
)

# =========================
# Average expression for certain gene, 'GENE'
# =========================
gene_name <- "GENE"   # <-- replace with the gene you want
gene_expr <- FetchData(seurat_combined, vars = gene_name)

# Get metadata
meta <- seurat_combined@meta.data

# Build dataframe
df_gene <- data.frame(
  cluster   = meta$cluster_named, # cluster_named for user-defined labels or else seurat_clusters
  sample    = meta$sample_id,
  condition = meta$condition,
  expr      = gene_expr[[gene_name]]
)

# Average per cluster per sample
avg_expr_gene <- df_gene %>%
  group_by(cluster, sample) %>%
  summarize(avg_expr = mean(expr, na.rm = TRUE), .groups = "drop")

# Pivot wide for easier export/heatmaps
avg_expr_wide <- avg_expr_gene %>%
  pivot_wider(names_from = sample, values_from = avg_expr)

# Save to CSV
out_file <- file.path(
  DIR_TABLES,
  paste0(gene_name, "_avg_expr_by_cluster_sample.csv")
)

write.csv(avg_expr_wide, out_file, row.names = FALSE)

# =========================
# Subset out clusters
# =========================
CELL_TYPE <- subset(seurat_combined, subset = seurat_clusters == "CLUSTER_NUMBER") # or cluster_named and user defined name
CELL_TYPE <- NormalizeData(CELL_TYPE) # Data normalization and scaling
CELL_TYPE <- FindVariableFeatures(CELL_TYPE)
CELL_TYPE <- ScaleData(CELL_TYPE)
CELL_TYPE <- RunPCA(CELL_TYPE)
CELL_TYPE <- FindNeighbors(CELL_TYPE, dims = N_PCS) # adjust dims as needed
CELL_TYPE <- FindClusters(CELL_TYPE, resolution = 0.4) # adjust resoltuion as needed
CELL_TYPE <- RunUMAP(CELL_TYPE, dims = N_PCS) # match dims to FindNeighbors

DimPlot(CELL_TYPE, reduction = "umap", label = TRUE, split.by = "condition") +
  ggtitle("Gene by cluster")

CELL_TYPE$seurat_clusters <- factor(CELL_TYPE$seurat_clusters)

# =========================
# DEG within clusters
# =========================
# Set identities to condition
Idents(CELL_TYPE) <- "condition"

# Initialize list to store results
deg_list_cell_type <- list()

# Loop through each cluster
clusters_cell_type <- levels(CELL_TYPE$cluster_named)
for (clust in clusters_cell_type) {
  message("Running DEG for cluster: ", clust)
  
  # Subset to just the current cluster
  cluster_subset_cell_type <- subset(CELL_TYPE, cluster_named == clust)
  
  # Set Idents again in the subset to condition
  Idents(cluster_subset_cell_type) <- "condition"
  
  # Count cells per condition
  condition_counts <- table(cluster_subset_cell_type$condition)
  
  # Skip if either condition has fewer than 3 cells
  if (length(condition_counts) < 2 || any(condition_counts < 3)) {
    message("Skipping cluster ", clust, ": too few cells in one or both conditions")
    next
  }
  
  # Run DEG
  deg_result_cell_type <- FindMarkers(
    cluster_subset_cell_type,
    ident.1 = "treated",
    ident.2 = "untreated",
    group.by = "condition",
    min.pct = 0.1,
    logfc.threshold = 0.25
  )
  
  # Annotate results
  deg_result_cell_type$cluster <- clust
  deg_result_cell_type$gene <- rownames(deg_result_cell_type)
  deg_list_astro[[clust]] <- deg_result_cell_type
}

# Combine results
cell_type_deg_combined <- bind_rows(deg_list_cell_type)
readr::write_csv(
  cell_type_deg_combined,
  file.path(DIR_TABLES, "CELL_TYPE_cells_DEGs_by_condition.csv", row.names = FALSE)
)

# =========================
# Heat map CSV file of avg gene expression
# =========================
avg_exp <- make_cluster_condition_csv_z(
  obj = seurat_combined,
  assay = "RNA",
  slot  = "data",
  top_n_by_var = NULL,   # or e.g., 500 for a compact heat map
  cap_z = 3, # good for visualization; set NULL to disable
  outdir = DIR_TABLES,
  prefix = NULL
)