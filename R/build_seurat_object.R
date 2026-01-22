# ============================
# Object build helper functions
# ============================


# ============================
# Build
# ============================
# Read a CSV (comma-separated) counts file like your A1 screenshot:
# first column = Gene-ID (rownames), remaining columns = cells
read_counts_csv <- function(path_csv) {
  stopifnot(file.exists(path_csv))
  dt <- fread(path_csv, sep = ",", header = TRUE, data.table = FALSE, check.names = FALSE)
  rn <- dt[[1]]; dt[[1]] <- NULL
  M <- as.matrix(dt)
  rownames(M) <- make.names(rn, unique = TRUE)
  # Drop completely empty columns if any
  keep <- colSums(!is.na(M)) > 0
  if (!all(keep)) M <- M[, keep, drop = FALSE]
  M
}

# Extract numeric ROI indices from an ROI zip; keep only ROIs that look like "Cell"
roi_indices_from_zip <- function(roi_zip) {
  stopifnot(file.exists(roi_zip))
  tmp <- tempfile(); dir.create(tmp)
  utils::unzip(roi_zip, exdir = tmp)
  files <- list.files(tmp, pattern = "\\.roi$", full.names = TRUE, recursive = TRUE)
  stopifnot(length(files) > 0)
  base  <- basename(files)
  keep  <- grepl("cell", base, ignore.case = TRUE)  # ignore "Annotation", etc.
  base  <- base[keep]
  idx   <- suppressWarnings(as.integer(stringr::str_extract(base, "\\d+")))
  sort(idx[!is.na(idx)])
}

# Build one Seurat object, forcing column names to ROI order and dropping any "Annotation" column
build_sample <- function(counts_path, roi_zip, sample_id, condition, sample_code) {
  # 1) counts
  M <- read_counts_csv(counts_path)
  
  # Drop any "Annotation" column if it sneaked into the matrix
  ann_cols <- which(grepl("annot", colnames(M), ignore.case = TRUE))
  if (length(ann_cols)) M <- M[, -ann_cols, drop = FALSE]
  
  # 2) ROI indices in order
  idx <- roi_indices_from_zip(roi_zip)
  
  # 3) Reconcile sizes (common cases)
  if (ncol(M) == length(idx)) {
    idx_use <- idx
  } else if (ncol(M) == length(idx) - 1L) {
    # counts already excluded the big outer ROI, ROI zip still has it → drop the biggest index
    idx_use <- idx[seq_len(ncol(M))]
  } else {
    stop("For ", sample_id, ": ncol(counts)=", ncol(M), " but ROI cells=", length(idx),
         ". Make sure the ROI zip matches this counts file and that the big 'Annotation' ROI ",
         "is removed from BOTH the zip and the CSV.")
  }
  
  # 4) Force cell names to match ROI index + sample code (unique across samples)
  colnames(M) <- paste0("Cell", idx_use, "_", sample_code)
  
  # 5) Seurat object + metadata
  sobj <- CreateSeuratObject(counts = M, project = sample_id, min.features = 0)
  sobj$sample_id  <- sample_id
  sobj$condition  <- condition
  sobj$roi_idx    <- idx_use
  sobj$roi_key    <- paste0(sample_id, "_", idx_use)
  sobj
}

# ============================
# QC helper
# ============================
# ---- get raw counts robustly (v5: layers; v4: slots) ----
get_raw_counts <- function(obj, assay = DefaultAssay(obj)) {
  a <- obj[[assay]]
  # v5 path: use LayerData
  if ("layers" %in% slotNames(a)) {
    lyr_names <- SeuratObject::Layers(a)
    lyr <- if ("counts" %in% lyr_names) "counts" else lyr_names[1]
    m <- SeuratObject::LayerData(a, layer = lyr)
    return(m)
  } else {
    # v4 fallback
    return(GetAssayData(obj, assay = assay, slot = "counts"))
  }
}

# ============================
# Heatmap/avg expression csv helper
# ============================
# Call JoinLayers only when available & needed (v5 layered assays)
maybe_join_layers <- function(obj, assay) {
  a <- obj[[assay]]
  has_layers <- "layers" %in% slotNames(a)
  can_join   <- "JoinLayers" %in% getNamespaceExports("SeuratObject")
  if (has_layers && can_join) obj <- SeuratObject::JoinLayers(obj, assay = assay)
  obj
}

make_cluster_condition_csv_z <- function(
    obj,
    outdir_tables,
    outdir_figures = NULL,
    prefix = NULL,
    assay = "RNA",
    slot  = "data",
    top_n_by_var = NULL,
    cap_z = 3,
    out_csv_means = "cluster_condition_means_hclust.csv",
    out_csv_z     = "cluster_condition_means_Z_hclust.csv"
) {
  dir.create(outdir_tables, showWarnings = FALSE, recursive = TRUE)
  
  if (!is.null(outdir_figures)) {
    dir.create(outdir_figures, showWarnings = FALSE, recursive = TRUE)
  }
  out_csv_means <- file.path(outdir_tables, out_csv_means)
  out_csv_z     <- file.path(outdir_tables, out_csv_z)
  
  stopifnot(all(c("seurat_clusters","condition","sample_id") %in% colnames(obj@meta.data)))
  
  obj <- maybe_join_layers(obj, assay = assay)
  DefaultAssay(obj) <- assay
  
  # Group by (sample | condition | cluster) to average per sample first
  obj$.__grp <- with(obj@meta.data, paste(sample_id, condition, seurat_clusters, sep = "|"))
  Idents(obj) <- ".__grp"
  
  avg <- AverageExpression(obj, assays = assay, slot = slot, verbose = FALSE)[[assay]]  # genes x groups
  
  # Parse group labels
  grp <- str_match(colnames(avg), "^(.+)\\|(.+)\\|(.+)$") |> as.data.frame()
  colnames(grp) <- c("full","sample_id","condition","cluster")
  grp$cluster <- as.character(grp$cluster)
  
  # Aggregate across samples → mean within (cluster, condition)
  keys <- paste0("Cluster_", grp$cluster, "_", grp$condition)
  idxs <- split(seq_len(ncol(avg)), keys)
  mat_means <- sapply(idxs, function(cols) {
    if (length(cols) == 1L) avg[, cols, drop = FALSE][, 1]
    else rowMeans(avg[, cols, drop = FALSE], na.rm = TRUE)
  })
  mat_means <- as.matrix(mat_means)  # genes x (Cluster_k_{treated|untreated})
  
  # Optional: keep top-N by variance for a cleaner heat map
  if (!is.null(top_n_by_var)) {
    v <- apply(mat_means, 1, var, na.rm = TRUE)
    keep <- names(sort(v, decreasing = TRUE))[seq_len(min(top_n_by_var, length(v)))]
    mat_means <- mat_means[keep, , drop = FALSE]
  }
  
  # Row-wise z-score (per gene across all columns)
  mat_z <- t(scale(t(mat_means)))
  mat_z[!is.finite(mat_z)] <- 0
  if (!is.null(cap_z) && is.finite(cap_z)) {
    mat_z[mat_z >  cap_z] <-  cap_z
    mat_z[mat_z < -cap_z] <- -cap_z
  }
  
  # Hierarchical ordering (use z for ordering, then apply to both)
  d_row <- dist(mat_z, method = "euclidean")
  d_col <- dist(t(mat_z), method = "euclidean")
  o_row <- order.dendrogram(as.dendrogram(hclust(d_row)))
  o_col <- order.dendrogram(as.dendrogram(hclust(d_col)))
  mat_means_ord <- mat_means[o_row, o_col, drop = FALSE]
  mat_z_ord     <- mat_z[o_row, o_col, drop = FALSE]
  
  # Write CSVs
  readr::write_csv(
    data.frame(gene = rownames(mat_means_ord), mat_means_ord, check.names = FALSE),
    out_csv_means
  )
  readr::write_csv(
    data.frame(gene = rownames(mat_z_ord), mat_z_ord, check.names = FALSE),
    out_csv_z
  )
  
  message("Wrote:\n  - ", out_csv_means, " (raw means)\n  - ", out_csv_z, " (Z-scores)")
  invisible(list(means = mat_means_ord, z = mat_z_ord))
}