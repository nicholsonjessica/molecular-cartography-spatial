#==========================
# Spatial plotting function
# Script to be used after analysis pipeline to generate images with color coded ROIs by cluster, transcripts, and DAPI image

project_sample <- function(
    seurat_combined,
    sample_id,
    dapi_path,
    roi_path,                        # .zip or folder of .roi files
    tx_path        = NULL,           # optional transcript table (X Y Z gene) or (X Y gene)
    outdir_figures = NULL,
    outdir_tables = NULL,
    prefix = NULL,

    # >>> ADDED: single gene toggle
    gene_of_interest = NULL,         # e.g., "GAD1"; NULL = don't plot transcripts
    # >>> ADDED: transcript styling
    tx_point_size   = 0.18,
    tx_alpha        = 0.8,
    
    # image / scale
    downsample     = 4,
    tx_delim       = "\t",
    # labels / toggles
    label_col      = "seurat_clusters",
    include_clusters = NULL,
    exclude_clusters = NULL,
    outlines_only  = FALSE,
    # filters for ROIs
    drop_name_regex = NULL,          # e.g., "(?i)outline|tissue|whole|mask"
    area_filter    = TRUE,
    area_quantile  = 0.999,
    area_frac      = 0.001,
    drop_unmatched = TRUE,
    # styling
    fill_alpha     = 0.25,
    outline_lwd    = 0.4,
    # memory helpers
    simplify_tol_px = 0.75           # simplify polygons by ~0.75 px at *downsampled* scale
) {
  suppressPackageStartupMessages({
    library(ijtiff); library(sf); library(dplyr); library(tibble)
    library(stringr); library(ggplot2); library(data.table)
    library(ggnewscale) # >>> ADDED
  })
  sf::sf_use_s2(FALSE)
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  # --- DAPI (grayscale raster, downsampled) -----------------------------------
  read_dapi_gray <- function(path, downsample = 4) {
    img <- ijtiff::read_tif(path)
    dm  <- dim(img)
    mat <- if (length(dm) == 2) img else if (length(dm) == 3) img[,,1, drop=TRUE] else img[,,1,1, drop=TRUE]
    if (!is.null(downsample) && downsample > 1) {
      mat <- mat[seq(1, nrow(mat), by = downsample),
                 seq(1, ncol(mat), by = downsample), drop = FALSE]
    }
    rng  <- range(mat, na.rm = TRUE); if (diff(rng) == 0) rng <- c(0, 1)
    cols <- grDevices::gray((mat - rng[1]) / diff(rng))
    ras  <- as.raster(matrix(cols, nrow = nrow(mat), ncol = ncol(mat)))
    list(raster = ras,
         height = nrow(mat), width = ncol(mat),
         height_full = dm[1], width_full = dm[2],
         ds = ifelse(is.null(downsample), 1, downsample))
  }
  
  # --- ROI I/O (preserve names) ----------------------------------------------
  read_roi_zip_preserve_names <- function(path) {
    if (dir.exists(path)) {
      files <- list.files(path, pattern="\\.roi$", full.names=TRUE, recursive=TRUE)
    } else if (file.exists(path) && grepl("\\.zip$", path, ignore.case=TRUE)) {
      tmp <- tempfile(); dir.create(tmp)
      utils::unzip(path, exdir = tmp)
      files <- list.files(tmp, pattern="\\.roi$", full.names=TRUE, recursive=TRUE)
    } else stop("ROI path not found: ", path)
    stopifnot(length(files) > 0)
    rois <- lapply(files, RImageJROI::read.ijroi)
    names(rois) <- basename(files)   # keep e.g. "123-Cell.roi"
    rois
  }
  
  # Parse polygons from RImageJROI list -> sf; parse idx from "n-Cell.roi"
  roi_to_sf_from_named <- function(rois, img_height_full, sample_id) {
    pls <- list(); nm_out <- character(0)
    add_poly <- function(mat, nm) {
      # ImageJ Y origin at top; flip to bottom-left for ggplot’s raster
      mat[,2] <- img_height_full - mat[,2]
      if (!all(mat[1,] == mat[nrow(mat),])) mat <- rbind(mat, mat[1, , drop=FALSE])
      pls[[length(pls)+1]] <<- sf::st_polygon(list(mat))
      nm_out <<- c(nm_out, nm)
    }
    for (nm in names(rois)) {
      r <- rois[[nm]]
      if (!is.null(r$coords) && nrow(r$coords) >= 3) {
        add_poly(as.matrix(r$coords[, c("x","y")]), nm); next
      }
      if (all(c("left","top","right","bottom") %in% names(r))) {
        box <- rbind(c(r$left,r$top), c(r$right,r$top),
                     c(r$right,r$bottom), c(r$left,r$bottom))
        add_poly(box, nm); next
      }
      if (is.list(r$coords)) {
        for (part in r$coords) {
          if (is.data.frame(part) && nrow(part) >= 3 && all(c("x","y") %in% names(part))) {
            add_poly(as.matrix(part[, c("x","y")]), nm)
          }
        }
      }
    }
    stopifnot(length(pls) > 0)
    sfc <- sf::st_sfc(pls, crs = sf::st_crs(NA))
    sf  <- sf::st_sf(roi_name = nm_out, geometry = sfc)
    
    # Extract leading integer before "-Cell" (robust to spaces/underscores)
    idx_str <- sub("^\\s*([0-9]+)\\s*[-_ ]\\s*Cell.*$", "\\1", sf$roi_name, perl = TRUE)
    # Fallback: last integer anywhere
    bad <- is.na(idx_str) | (idx_str == sf$roi_name)
    if (any(bad)) idx_str[bad] <- sub(".*?(\\d+)(?:\\.roi)?$", "\\1", sf$roi_name[bad], perl = TRUE)
    
    sf$roi_idx    <- suppressWarnings(as.integer(idx_str))
    sf$sample_id  <- sample_id
    sf
  }
  
  # Optional transcripts
  read_tx_xyz_gene <- function(path, img_height_full, delim = "\t") {
    if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
    if (file.info(path)$size == 0) return(NULL)
    dt <- data.table::fread(path, sep = delim, header = FALSE,
                            data.table = FALSE, fill = TRUE, strip.white = TRUE)
    keep <- which(colSums(!(is.na(dt) | dt == "")) > 0)
    dt <- dt[, keep, drop = FALSE]
    if (ncol(dt) < 4) return(NULL)
    out <- dt[, c(1, 2, 4)]
    colnames(out) <- c("x","y","gene")
    out$x <- as.numeric(out$x)
    out$y <- img_height_full - as.numeric(out$y)
    out$gene <- as.character(out$gene)
    out
  }
  
  # -------------------- pipeline --------------------
  dapi     <- read_dapi_gray(dapi_path, downsample = downsample)
  
  rois     <- read_roi_zip_preserve_names(roi_path)
  cells_sf <- roi_to_sf_from_named(rois, img_height_full = dapi$height_full, sample_id = sample_id)
  
  # scale ROI geometry to downsampled pixel space
  sf::st_crs(cells_sf) <- sf::st_crs(NA)
  sf::st_geometry(cells_sf) <- sf::st_geometry(cells_sf) / dapi$ds
  
  # Simplify polygons to save RAM (tolerance in *downsampled* pixels)
  if (!is.null(simplify_tol_px) && simplify_tol_px > 0) {
    cells_sf <- suppressWarnings(sf::st_simplify(cells_sf, dTolerance = simplify_tol_px, preserveTopology = TRUE))
  }
  
  # Seurat metadata (join by sample_id + roi_idx)
  meta <- seurat_combined@meta.data |>
    tibble::as_tibble(rownames = "barcode") |>
    mutate(
      sample_id       = as.character(sample_id),
      sample_id       = .data$sample_id, # ensure column exists
      roi_idx         = suppressWarnings(as.integer(.data$roi_idx)),
      seurat_clusters = as.character(.data$seurat_clusters %||% Idents(seurat_combined))
    )
  
  # Some builds of Seurat may have lost sample_id/roi_idx; if so, recover from barcode
  if (!("sample_id" %in% names(meta)) || !("roi_idx" %in% names(meta))) {
    m <- stringr::str_match(meta$barcode, "^Cell[\\._-]?([0-9]+)_(.+)$")
    if (!("roi_idx" %in% names(meta))) meta$roi_idx <- suppressWarnings(as.integer(m[,2]))
    if (!("sample_id" %in% names(meta))) meta$sample_id <- as.character(meta$sample_id %||% NA)
  }
  
  # Filter meta to this sample only, drop NAs, and deduplicate (one label per ROI)
  meta_keyed <- meta |>
    filter(.data$sample_id == !!sample_id) |>
    filter(!is.na(.data$roi_idx)) |>
    distinct(sample_id, roi_idx, .keep_all = TRUE) |>
    select(sample_id, roi_idx, seurat_clusters)
  
  # Join strictly on (sample_id, roi_idx)
  cells_sf <- cells_sf |>
    left_join(meta_keyed, by = c("sample_id", "roi_idx"))
  
  # optional name filters
  if (!is.null(drop_name_regex) && nzchar(drop_name_regex)) {
    cells_sf <- filter(cells_sf, !grepl(drop_name_regex, .data$roi_name, perl = TRUE))
  }
  
  # area-based filters (optional)
  if (isTRUE(area_filter)) {
    areas <- suppressWarnings(as.numeric(sf::st_area(cells_sf)))
    q_hi  <- stats::quantile(areas, probs = area_quantile, na.rm = TRUE)
    q_lo  <- max(min(areas, na.rm = TRUE) * area_frac, 0)
    keep  <- areas <= q_hi & areas >= q_lo
    cells_sf <- cells_sf[keep, , drop = FALSE]
  }
  
  # Drop unmatched (no cluster)
  if (isTRUE(drop_unmatched)) {
    cells_sf <- filter(cells_sf, !is.na(.data$seurat_clusters))
  }
  
  # Choose label column
  if (!label_col %in% names(cells_sf)) label_col <- "seurat_clusters"
  cells_sf$._label <- as.character(cells_sf[[label_col]])
  
  # include / exclude filters
  if (!is.null(include_clusters)) cells_sf <- filter(cells_sf, ._label %in% include_clusters)
  if (!is.null(exclude_clusters)) cells_sf <- filter(cells_sf, !._label %in% exclude_clusters)
  
  # Factor levels from all clusters in seurat_combined (stable palette across samples)
  all_levels <- sort(unique(as.character(seurat_combined$seurat_clusters)))
  cells_sf$._label <- factor(cells_sf$._label, levels = all_levels)
  # Palette
  global_pal <- setNames(scales::hue_pal()(length(all_levels)), all_levels)
  fill_override <- unname(global_pal[levels(cells_sf$._label)])
  
  # -------------------- plot --------------------
  g <- ggplot() +
    annotation_raster(dapi$raster, xmin = 0, xmax = dapi$width, ymin = 0, ymax = dapi$height)
  
  if (!isTRUE(outlines_only)) { # fill
    g <- g + geom_sf(
      data = cells_sf, aes(fill = ._label),
      color = NA, alpha = fill_alpha, inherit.aes = FALSE, show.legend = FALSE
    )
  }
  
  g <- g + geom_sf( # outlines
    data = cells_sf, aes(color = ._label),
    fill = NA, linewidth = outline_lwd, inherit.aes = FALSE, show.legend = TRUE
  ) +
    scale_color_manual(
      values = global_pal, name = label_col, drop = TRUE,
      guide = guide_legend(override.aes = list(fill = fill_override, linewidth = 0, alpha = 1))
    ) +
    scale_fill_manual(values = global_pal, guide = "none", drop = TRUE) +
    coord_sf(xlim = c(0, dapi$width), ylim = c(0, dapi$height),
             expand = FALSE, default_crs = sf::st_crs(NA)) +
    theme_void(base_size = 12) +
    theme(legend.position = "right") +
    labs(title = paste("Sample", sample_id))
  
  # Transcripts overlay (single gene)
  if (!is.null(tx_path) && !is.null(gene_of_interest) && nzchar(gene_of_interest)) {
    # read: try hinted delimiter first, then whitespace; accept 3 or 4 cols
    dt <- tryCatch(
      data.table::fread(tx_path, sep = tx_delim, header = FALSE,
                        data.table = TRUE, fill = TRUE, strip.white = TRUE),
      error = function(e) NULL
    )
    if (is.null(dt) || ncol(dt) < 3L) {
      dt <- data.table::fread(tx_path, sep = "", header = FALSE,
                              data.table = TRUE, fill = TRUE, strip.white = TRUE)
    }
    # drop empty cols; name columns
    keep <- which(colSums(!(is.na(dt) | dt == "")) > 0)
    if (length(keep)) dt <- dt[, ..keep]
    nc <- ncol(dt)
    if (nc >= 4L) {
      data.table::setnames(dt, c(1,2,nc), c("x","y","gene"))
      dt <- dt[, .(x, y, gene)]
    } else if (nc == 3L) {
      data.table::setnames(dt, c("x","y","gene"))
    } else {
      warning("[transcripts] expected ≥3 columns; found: ", nc)
      dt <- NULL
    }
    
    if (!is.null(dt) && nrow(dt) > 0) {
      dt[, gene := trimws(as.character(gene))]
      dt <- dt[tolower(gene) == tolower(gene_of_interest)]
      if (nrow(dt) > 0) {
        tx_df <- data.frame(
          x = as.numeric(dt$x) / dapi$ds,
          y = (dapi$height_full - as.numeric(dt$y)) / dapi$ds
        )
        tx_df <- subset(tx_df, is.finite(x) & is.finite(y) &
                          x >= 0 & y >= 0 & x <= dapi$width & y <= dapi$height)
        
        if (nrow(tx_df) > 0) {
          # separate color scale so cluster colors stay intact
          g <- g +
            ggnewscale::new_scale_color() +
            geom_point(
              data = tx_df,
              aes(x = x, y = y, color = !!as.name("gene_of_interest")),  # constant legend label
              size = tx_point_size, alpha = tx_alpha, inherit.aes = FALSE, show.legend = TRUE
            ) +
            scale_color_manual(
              name   = "Transcripts",
              values = setNames("#000000", gene_of_interest),  # default black; change if you want
              labels = setNames(gene_of_interest, gene_of_interest),
              breaks = gene_of_interest
            )
        } else {
          message("[transcripts] no points inside bounds after transform")
        }
      } else {
        message("[transcripts] gene not found in file: ", gene_of_interest)
      }
    } else {
      message("[transcripts] could not parse transcript table")
    }
  }
  
  match_rate <- mean(!is.na(cells_sf$seurat_clusters))
  message(sprintf("[project_sample] %s: polygons=%d, matched=%.1f%%",
                  sample_id, nrow(cells_sf), 100 * match_rate))
  
  list(plot = g, cells_sf = cells_sf, dapi = dapi)
}

res_sample <- project_sample(
  sample_id        = "sample",
  dapi_path        = "path/dapi.tiff",
  roi_path         = "sample_rois.zip",
  tx_path          = "path/transcripts.txt",
  seurat_combined  = seurat_combined,
  gene_of_interest = "GENE",   # <- just set this to plot transcripts
  downsample       = 4,
  drop_name_regex  = "(?i)outline|tissue|whole|mask",
  fill_alpha       = 0, # changes ROI fill opacity
  outline_lwd      = 0, # changes ROI outline opacity
  simplify_tol_px  = 0.75,
  tx_point_size    = 0.025,     # transcript dot size
  tx_alpha         = 0.75
)

print(res_sample$plot)

ggsave(
  file.path(DIR_FIGURES, "sample_transcr.png"),
  plot     = res_sample$plot,
  width    = 8, height = 8, dpi = 300,  # adjust size/dpi for full res
  limitsize = FALSE
)
