#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(ggplot2)
  library(glue)
  library(dbscan)
  library(RcppML)
})

source("utils/R/niche.R")

`%||%` <- function(x, y) if (is.null(x) || x == "") y else x

parse_args <- function(argv) {
  args <- list()
  if (length(argv) %% 2 != 0) {
    stop("Arguments must be provided as --key value pairs.")
  }
  if (length(argv) == 0) return(args)
  keys <- argv[seq(1, length(argv), by = 2)]
  vals <- argv[seq(2, length(argv), by = 2)]
  keys <- sub("^--", "", keys)
  for (i in seq_along(keys)) {
    args[[keys[i]]] <- vals[i]
  }
  args
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

as_num <- function(x, default) {
  if (is.null(x)) return(default)
  as.numeric(x)
}

as_int <- function(x, default) {
  if (is.null(x)) return(default)
  as.integer(x)
}

scan_neighbor_counts <- function(coords,
                                 query_cells,
                                 r_values,
                                 max_query_cells = 5000,
                                 use_pbapply = TRUE) {
  if (!all(query_cells %in% rownames(coords))) {
    query_cells <- intersect(query_cells, rownames(coords))
  }
  if (length(query_cells) == 0) {
    stop("No query cells available for r scan.")
  }
  if (length(query_cells) > max_query_cells) {
    set.seed(1024)
    query_cells <- sample(query_cells, max_query_cells)
  }

  Q <- coords[query_cells, , drop = FALSE]
  apply_fun <- if (use_pbapply && requireNamespace("pbapply", quietly = TRUE)) {
    pbapply::pblapply
  } else {
    lapply
  }

  res_list <- apply_fun(r_values, function(r) {
    nbrs <- dbscan::frNN(coords, eps = r, query = Q, sort = FALSE)
    counts <- vapply(nbrs$id, length, integer(1))
    data.frame(
      r = r,
      mean_neighbors = mean(counts),
      median_neighbors = median(counts),
      p90_neighbors = as.numeric(stats::quantile(counts, 0.9)),
      query_n = length(query_cells)
    )
  })
  do.call(rbind, res_list)
}

choose_stable_r <- function(scan_df,
                            stable_delta = 0.1,
                            stable_window = 3,
                            min_neighbors = 30,
                            target_neighbors = NA_real_) {
  if (nrow(scan_df) <= 1) {
    return(list(r = scan_df$r[1], reason = "single_r"))
  }
  median_vals <- scan_df$median_neighbors
  delta <- diff(median_vals) / pmax(median_vals[-length(median_vals)], 1)
  ok <- delta <= stable_delta

  if (length(ok) >= stable_window) {
    for (i in seq_len(length(ok) - stable_window + 1)) {
      if (all(ok[i:(i + stable_window - 1)]) &&
          median_vals[i] >= min_neighbors) {
        return(list(r = scan_df$r[i + stable_window - 1], reason = "stable_window"))
      }
    }
  }

  if (!is.na(target_neighbors)) {
    idx <- which.min(abs(median_vals - target_neighbors))
    return(list(r = scan_df$r[idx], reason = "closest_target"))
  }

  idx <- which(median_vals >= min_neighbors)
  if (length(idx) > 0) {
    return(list(r = scan_df$r[min(idx)], reason = "min_neighbors"))
  }

  list(r = scan_df$r[which.max(median_vals)], reason = "max_neighbors")
}

cluster_average <- function(mat, clusters) {
  cluster_names <- names(clusters)
  clusters <- as.character(clusters)
  names(clusters) <- cluster_names
  split_cells <- split(names(clusters), clusters)
  avg <- sapply(split_cells, function(cells) {
    if (!all(cells %in% rownames(mat))) {
      cells <- intersect(cells, rownames(mat))
    }
    if (length(cells) == 0) {
      rep(NA_real_, ncol(mat))
    } else {
      colMeans(mat[cells, , drop = FALSE])
    }
  })
  out <- t(avg)
  rownames(out) <- names(split_cells)
  colnames(out) <- colnames(mat)
  out
}

plot_r_scan <- function(scan_df, out_file) {
  p <- ggplot(scan_df, aes(r, median_neighbors)) +
    geom_line(linewidth = 0.8, color = "#2C3E50") +
    geom_point(size = 1.6, color = "#2C3E50") +
    theme_minimal(base_size = 11) +
    labs(x = "r", y = "Median neighbor count")
  ggsave(out_file, p, width = 4, height = 3, dpi = 200)
}

plot_spatial_clusters <- function(seu,
                                  caf_clusters,
                                  out_file,
                                  max_background = 100000,
                                  max_caf = 200000) {
  meta <- seu@meta.data
  if (!all(c("Center_X", "Center_Y") %in% colnames(meta))) {
    warning("Center_X/Center_Y not found; skip spatial plot.")
    return(invisible(NULL))
  }

  df <- data.frame(
    cell_id = rownames(meta),
    x = meta$Center_X,
    y = meta$Center_Y,
    cluster = NA_character_
  )
  df$cluster[match(names(caf_clusters), df$cell_id)] <- as.character(caf_clusters)

  bg_idx <- which(is.na(df$cluster))
  caf_idx <- which(!is.na(df$cluster))
  if (length(bg_idx) > max_background) {
    set.seed(1024)
    bg_idx <- sample(bg_idx, max_background)
  }
  if (length(caf_idx) > max_caf) {
    set.seed(1024)
    caf_idx <- sample(caf_idx, max_caf)
  }

  bg_df <- df[bg_idx, , drop = FALSE]
  caf_df <- df[caf_idx, , drop = FALSE]

  p <- ggplot() +
    geom_point(
      data = bg_df,
      aes(x = x, y = y),
      color = "grey80",
      size = 0.12,
      alpha = 0.5
    ) +
    geom_point(
      data = caf_df,
      aes(x = x, y = y, color = cluster),
      size = 0.18,
      alpha = 0.9
    ) +
    coord_equal() +
    theme_void(base_size = 11) +
    theme(legend.position = "right")

  ggsave(out_file, p, width = 6, height = 5, dpi = 200)
}

plot_heatmap <- function(mat, out_file, title = NULL) {
  df <- as.data.frame(mat)
  df$cluster <- rownames(df)
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required for plotting heatmaps.")
  }
  df_long <- tidyr::pivot_longer(df, cols = -cluster, names_to = "cell_type", values_to = "value")

  p <- ggplot(df_long, aes(cell_type, cluster, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#1F78B4") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(x = NULL, y = NULL, title = title)

  ggsave(out_file, p, width = 6.5, height = 4, dpi = 200)
}

plot_stacked_bar <- function(mat, out_file, title = NULL) {
  df <- as.data.frame(mat)
  df$cluster <- rownames(df)
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required for plotting stacked bars.")
  }
  df_long <- tidyr::pivot_longer(df, cols = -cluster, names_to = "cell_type", values_to = "fraction")

  p <- ggplot(df_long, aes(cluster, fraction, fill = cell_type)) +
    geom_col(position = "stack", width = 0.8) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(x = NULL, y = "Fraction", title = title)

  ggsave(out_file, p, width = 6.5, height = 4, dpi = 200)
}

argv <- commandArgs(trailingOnly = TRUE)
opts <- parse_args(argv)

input_dir <- opts$input %||% "data"
output_dir <- opts$output %||% file.path("data", "caf_niche")
pattern <- opts$pattern %||% "\\.qs$"
samples <- opts$samples %||% ""

anno_coarse <- opts$anno_coarse %||% "annotation_coarse"
anno_fine <- opts$anno_fine %||% "annotation_fine"
caf_label <- opts$caf_label %||% "CAF"
reduction <- opts$reduction %||% "spatial"

r_min <- as_num(opts$r_min, 20)
r_max <- as_num(opts$r_max, 200)
r_step <- as_num(opts$r_step, 10)
stable_delta <- as_num(opts$stable_delta, 0.1)
stable_window <- as_int(opts$stable_window, 3)
min_neighbors <- as_int(opts$min_neighbors, 30)
target_neighbors <- as_num(opts$target_neighbors, NA_real_)
max_query_cells <- as_int(opts$max_query_cells, 5000)
use_pbapply <- as_bool(opts$use_pbapply, TRUE)

k_list <- opts$k_list %||% "4,5"
k_list <- as.integer(strsplit(k_list, ",")[[1]])
cluster_resolution <- as_num(opts$cluster_resolution, 0.1)
nmf_cores <- as_int(opts$nmf_cores, 5)
seed <- as_int(opts$seed, 1024)

max_background <- as_int(opts$max_background, 100000)
max_plot_caf <- as_int(opts$max_plot_caf, 200000)
detail_fine <- as_bool(opts$detail_fine, TRUE)
coord_scale <- as_num(opts$coord_scale, 1.0)

input_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
if (samples != "") {
  sample_list <- strsplit(samples, ",")[[1]]
  input_files <- input_files[basename(input_files) %in% paste0(sample_list, ".qs")]
}
if (length(input_files) == 0) {
  stop("No input files found.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

summary_rows <- list()

for (file in input_files) {
  sample_id <- sub("\\.qs$", "", basename(file))
  message(glue("Processing {sample_id} ..."))

  seu <- qs::qread(file)
  if (!anno_coarse %in% colnames(seu@meta.data)) {
    warning(glue("{sample_id}: missing {anno_coarse}; skip."))
    next
  }
  if (!reduction %in% names(seu@reductions)) {
    warning(glue("{sample_id}: reduction {reduction} missing; skip."))
    next
  }

  if (coord_scale != 1.0) {
    coords_scaled <- Seurat::Embeddings(seu, reduction = reduction) * coord_scale
    seu[[reduction]]@cell.embeddings <- coords_scaled
  }

  caf_cells <- rownames(seu@meta.data)[seu@meta.data[[anno_coarse]] == caf_label]
  if (length(caf_cells) == 0) {
    warning(glue("{sample_id}: no CAF cells; skip."))
    next
  }

  coords <- Seurat::Embeddings(seu, reduction = reduction)
  r_values <- seq(r_min, r_max, by = r_step)
  scan_df <- scan_neighbor_counts(
    coords = coords,
    query_cells = caf_cells,
    r_values = r_values,
    max_query_cells = max_query_cells,
    use_pbapply = use_pbapply
  )

  choice <- choose_stable_r(
    scan_df,
    stable_delta = stable_delta,
    stable_window = stable_window,
    min_neighbors = min_neighbors,
    target_neighbors = target_neighbors
  )
  r_choice <- choice$r

  sample_dir <- file.path(output_dir, sample_id)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(scan_df, file.path(sample_dir, "r_scan.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  plot_r_scan(scan_df, file.path(sample_dir, "r_scan.png"))

  niche_res <- CalNichMatrix(
    seu,
    query.cells = caf_cells,
    r = r_choice,
    reduction = reduction,
    anno_col = anno_coarse,
    use_pbapply = use_pbapply,
    return_neighbors = FALSE
  )

  niche_mat_norm <- niche_res$niche.mat.norm
  niche_mat_correct <- niche_res$niche.mat.correct
  niche_mat_correct[!is.finite(niche_mat_correct)] <- 0

  if (detail_fine && anno_fine %in% colnames(seu@meta.data)) {
    fine_res <- CalNichMatrix(
      seu,
      query.cells = caf_cells,
      r = r_choice,
      reduction = reduction,
      anno_col = anno_fine,
      use_pbapply = use_pbapply,
      return_neighbors = FALSE
    )
    fine_norm <- fine_res$niche.mat.norm
  } else {
    fine_norm <- NULL
  }

  seu_caf <- subset(seu, cells = caf_cells)

  for (k in k_list) {
    nmf_res <- RunNMF(niche_mat_correct, k = k, cores = nmf_cores, seed = seed)
    W <- nmf_res$W
    H <- nmf_res$H

    reduction_name <- paste0("nmf_k", k)
    key_name <- paste0("NMF", k, "_")
    seu_caf[[reduction_name]] <- CreateDimReducObject(
      embeddings = W[rownames(seu_caf@meta.data), , drop = FALSE],
      assay = "RNA",
      key = key_name
    )
    seu_caf <- FindNeighbors(seu_caf, reduction = reduction_name, dims = 1:k, verbose = FALSE)

    cluster_name <- paste0("CAF.sub_k", k)
    seu_caf <- FindClusters(
      seu_caf,
      resolution = cluster_resolution,
      algorithm = 4,
      cluster.name = cluster_name,
      random.seed = max(seed, 1),
      verbose = FALSE
    )
    clusters <- seu_caf[[cluster_name]][, 1]
    seu_caf[[cluster_name]] <- factor(paste0("CAF-", clusters))

    cluster_vec <- seu_caf[[cluster_name]][, 1]
    names(cluster_vec) <- rownames(seu_caf@meta.data)

    cluster_counts <- table(cluster_vec)
    cluster_summary <- data.frame(
      sample_id = sample_id,
      k = k,
      cluster = names(cluster_counts),
      n_cells = as.integer(cluster_counts),
      fraction = as.numeric(cluster_counts) / length(cluster_vec)
    )
    utils::write.table(
      cluster_summary,
      file.path(sample_dir, glue("caf_cluster_summary_k{k}.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    cluster_df <- data.frame(
      cell_id = names(cluster_vec),
      caf_cluster = as.character(cluster_vec)
    )
    utils::write.table(
      cluster_df,
      file.path(sample_dir, glue("caf_clusters_k{k}.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    utils::write.table(
      H,
      file.path(sample_dir, glue("nmf_H_k{k}.tsv")),
      sep = "\t",
      quote = FALSE
    )

    niche_avg_norm <- cluster_average(niche_mat_norm, cluster_vec)
    niche_avg_correct <- cluster_average(niche_mat_correct, cluster_vec)

    utils::write.table(
      niche_avg_norm,
      file.path(sample_dir, glue("niche_average_norm_k{k}.tsv")),
      sep = "\t",
      quote = FALSE
    )
    utils::write.table(
      niche_avg_correct,
      file.path(sample_dir, glue("niche_average_correct_k{k}.tsv")),
      sep = "\t",
      quote = FALSE
    )

    plot_spatial_clusters(
      seu,
      cluster_vec,
      out_file = file.path(sample_dir, glue("caf_spatial_k{k}.png")),
      max_background = max_background,
      max_caf = max_plot_caf
    )

    plot_heatmap(
      niche_avg_correct,
      out_file = file.path(sample_dir, glue("niche_heatmap_k{k}.png")),
      title = glue("{sample_id} CAF niche (corrected, k={k})")
    )
    plot_stacked_bar(
      niche_avg_norm,
      out_file = file.path(sample_dir, glue("niche_bar_k{k}.png")),
      title = glue("{sample_id} CAF niche (norm, k={k})")
    )

    if (!is.null(fine_norm)) {
      fine_avg <- cluster_average(fine_norm, cluster_vec)
      utils::write.table(
        fine_avg,
        file.path(sample_dir, glue("niche_fine_average_norm_k{k}.tsv")),
        sep = "\t",
        quote = FALSE
      )
    }

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      sample_id = sample_id,
      k = k,
      r_choice = r_choice,
      r_reason = choice$reason,
      caf_cells = length(caf_cells)
    )
  }

  rm(seu)
  gc()
}

summary_df <- do.call(rbind, summary_rows)
if (!is.null(summary_df) && nrow(summary_df) > 0) {
  utils::write.table(
    summary_df,
    file.path(output_dir, "summary_r_choice.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}
