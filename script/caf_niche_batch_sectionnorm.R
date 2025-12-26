#!/usr/bin/env Rscript

# CAF niche batch (coarse, section-normalized):
# 1) pick r (adaptive vs baseline or fixed)
# 2) compute niche matrices + section normalization
# 3) NMF + clustering
# 4) summaries + plots
suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(ggplot2)
  library(glue)
  library(dbscan)
  library(RcppML)
  library(ggsci)
})

if (requireNamespace("pbapply", quietly = TRUE)) {
  pbapply::pboptions(type = "timer")
}

`%||%` <- function(x, y) if (is.null(x) || x == "") y else x

# Fixed coordinate scaling: pixel -> micron.
# Assumption: 1 pixel = 50 microns.
COORD_SCALE <- 50
# Baseline r (um) to compare against adaptive r.
R_BASELINE <- 80
# Global RNG seed for reproducibility (overwritten by cfg$seed).
RNG_SEED <- 1024

CANONICAL_COARSE <- character(0)

enforce_canonical_levels <- function(seu, anno_col, sample_id) {
  x <- as.character(seu@meta.data[[anno_col]])
  unknown <- sort(unique(x[!is.na(x) & !x %in% CANONICAL_COARSE]))
  if (length(unknown) > 0) {
    stop(glue("{sample_id}: unknown {anno_col} values: {paste(unknown, collapse = ', ')}"))
  }
  seu@meta.data[[anno_col]] <- factor(x, levels = CANONICAL_COARSE)
  seu
}

align_canonical_columns <- function(mat, sample_id, matrix_name) {
  if (is.null(colnames(mat))) {
    stop(glue("{sample_id}: {matrix_name} missing column names."))
  }
  extras <- setdiff(colnames(mat), CANONICAL_COARSE)
  if (length(extras) > 0) {
    stop(glue("{sample_id}: {matrix_name} has unexpected columns: {paste(extras, collapse = ', ')}"))
  }
  missing <- setdiff(CANONICAL_COARSE, colnames(mat))
  if (length(missing) > 0) {
    zeros <- matrix(0, nrow = nrow(mat), ncol = length(missing))
    colnames(zeros) <- missing
    mat <- cbind(mat, zeros)
  }
  mat[, CANONICAL_COARSE, drop = FALSE]
}

compute_section_prop <- function(meta, anno_col, sample_id) {
  x <- meta[[anno_col]]
  x <- factor(as.character(x), levels = CANONICAL_COARSE)
  counts <- table(x)
  total <- sum(counts)
  if (total == 0) {
    prop <- rep(NA_real_, length(CANONICAL_COARSE))
  } else {
    prop <- as.numeric(counts) / total
  }
  names(prop) <- CANONICAL_COARSE
  prop
}

compute_log2_enrichment <- function(niche_counts, section_prop, alpha = 0.5) {
  if (is.null(names(section_prop))) {
    names(section_prop) <- colnames(niche_counts)
  }
  section_prop <- section_prop[colnames(niche_counts)]
  N <- rowSums(niche_counts)
  expected <- outer(N, section_prop, "*")
  rownames(expected) <- rownames(niche_counts)
  colnames(expected) <- colnames(niche_counts)

  log2_enrich <- log2((niche_counts + alpha) / (expected + alpha))
  log2_enrich[!is.finite(log2_enrich) | is.na(log2_enrich)] <- 0
  list(log2_enrich = log2_enrich, expected = expected, total = N)
}

summarize_support <- function(niche_counts, expected, total_neighbors) {
  quantiles <- function(x) {
    stats::quantile(x, probs = c(0.1, 0.5, 0.9), na.rm = TRUE, names = FALSE)
  }
  total_stats <- quantiles(total_neighbors)
  support_rows <- list(
    data.frame(
      cell_type = "neighbor_total",
      x_mean = mean(total_neighbors, na.rm = TRUE),
      x_p10 = total_stats[1],
      x_median = total_stats[2],
      x_p90 = total_stats[3],
      E_mean = NA_real_,
      E_p10 = NA_real_,
      E_median = NA_real_,
      E_p90 = NA_real_
    )
  )
  for (cell_type in colnames(niche_counts)) {
    x <- niche_counts[, cell_type]
    e <- expected[, cell_type]
    x_stats <- quantiles(x)
    e_stats <- quantiles(e)
    support_rows[[length(support_rows) + 1]] <- data.frame(
      cell_type = cell_type,
      x_mean = mean(x, na.rm = TRUE),
      x_p10 = x_stats[1],
      x_median = x_stats[2],
      x_p90 = x_stats[3],
      E_mean = mean(e, na.rm = TRUE),
      E_p10 = e_stats[1],
      E_median = e_stats[2],
      E_p90 = e_stats[3]
    )
  }
  do.call(rbind, support_rows)
}

get_root_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) > 0) {
    return(normalizePath(file.path(dirname(file_arg[1]), "..")))
  }
  normalizePath(getwd())
}

root_dir <- get_root_dir()
source(file.path(root_dir, "utils", "R", "niche.R"))
source(file.path(root_dir, "utils", "R", "coarse_schema.R"))

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

to_bool <- function(x) {
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

get_opt <- function(opts, key, default, cast = NULL) {
  val <- opts[[key]]
  if (is.null(val) || val == "") return(default)
  if (is.null(cast)) return(val)
  cast(val)
}

# CLI parameters and defaults (used by --help).
PARAM_HELP <- list(
  input = "Input directory with Seurat .qs files",
  output = "Output directory for CAF-niche results",
  samples = "Comma-separated sample IDs (e.g., 13P,20P)",
  pattern = "Regex pattern for .qs files",
  anno_coarse = "Coarse annotation column (used for niche composition)",
  anno_fine = "Fine annotation column (optional; off by default)",
  coarse_schema = "Coarse schema: base or myeloid_refined",
  refine_myeloid = "Refine Myeloid into Macrophage/Neutrophil/Mast (true/false)",
  myeloid_fine_cols = "Comma-separated meta columns used to refine Myeloid",
  caf_label = "Label in anno_coarse defining CAF cells",
  reduction = "Seurat reduction with spatial coords (default: spatial)",
  r_min = "Min r (um) for scan; clipped to >=40",
  r_max = "Max r (um) for scan; clipped to <=100",
  r_step = "Step size for r scan (um)",
  r_fixed = "Fixed r (um); skip r scan if set",
  stable_delta = "Max relative change for stable window",
  stable_window = "Number of consecutive steps for stable window",
  min_neighbors = "Minimum median neighbors for stability",
  target_neighbors = "Optional target neighbor count",
  max_query_cells = "Subsample CAFs for r scan only",
  use_pbapply = "Use pbapply progress bars if available",
  k_list = "Comma-separated NMF k values",
  use_section_normalization = "Use section-normalized niche matrix for NMF (true/false)",
  log_enrich_alpha = "Pseudocount alpha for log2 enrichment",
  cluster_resolution = "Leiden resolution for CAF clustering",
  max_caf_clusters = "Upper bound on CAF cluster count",
  min_cluster_resolution = "Lower bound for cluster resolution",
  cluster_resolution_step = "Resolution step when reducing clusters",
  nmf_cores = "RcppML threads",
  seed = "Random seed",
  max_plot_cells = "Max niche cells to plot (proportional sampling)",
  max_plot_all_cells = "Max background cells to plot",
  detail_fine = "Also compute fine niche matrix (default: false)",
  plot_only = "Only redraw plots using existing results (skip r scan/NMF)"
)

print_help <- function(defaults) {
  cat("CAF niche batch (coarse) - parameters\n")
  for (key in names(PARAM_HELP)) {
    def <- defaults[[key]]
    def_str <- if (is.na(def) || is.null(def)) "NA" else as.character(def)
    cat(glue("  --{key} {PARAM_HELP[[key]]} (default: {def_str})\n"))
  }
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
    set.seed(RNG_SEED)
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

select_r_with_baseline <- function(scan_df,
                                   adaptive_choice,
                                   r_baseline = R_BASELINE,
                                   min_neighbors = 30,
                                   target_neighbors = NA_real_) {
  if (nrow(scan_df) == 0) {
    return(list(r = adaptive_choice$r, reason = "adaptive_only"))
  }
  medians <- scan_df$median_neighbors
  r_values <- scan_df$r

  idx_adapt <- which.min(abs(r_values - adaptive_choice$r))
  idx_base <- which.min(abs(r_values - r_baseline))

  median_adapt <- medians[idx_adapt]
  median_base <- medians[idx_base]

  instability <- function(idx) {
    if (length(medians) == 1) return(0)
    if (idx <= 1) {
      delta <- abs(medians[2] - medians[1]) / pmax(medians[1], 1)
      return(delta)
    }
    if (idx >= length(medians)) {
      delta <- abs(medians[length(medians)] - medians[length(medians) - 1]) /
        pmax(medians[length(medians) - 1], 1)
      return(delta)
    }
    delta_prev <- abs(medians[idx] - medians[idx - 1]) / pmax(medians[idx - 1], 1)
    delta_next <- abs(medians[idx + 1] - medians[idx]) / pmax(medians[idx], 1)
    max(delta_prev, delta_next)
  }

  score <- function(idx, median_val) {
    penalty <- if (is.na(median_val) || median_val < min_neighbors) {
      100 + (min_neighbors - median_val)
    } else {
      0
    }
    target_score <- if (!is.na(target_neighbors)) {
      abs(median_val - target_neighbors)
    } else {
      0
    }
    penalty + target_score + instability(idx)
  }

  score_adapt <- score(idx_adapt, median_adapt)
  score_base <- score(idx_base, median_base)

  if (score_base < score_adapt) {
    list(
      r = r_values[idx_base],
      reason = "baseline_better",
      r_adaptive = r_values[idx_adapt],
      r_baseline = r_values[idx_base],
      median_adaptive = median_adapt,
      median_baseline = median_base,
      score_adaptive = score_adapt,
      score_baseline = score_base
    )
  } else {
    list(
      r = r_values[idx_adapt],
      reason = "adaptive_better",
      r_adaptive = r_values[idx_adapt],
      r_baseline = r_values[idx_base],
      median_adaptive = median_adapt,
      median_baseline = median_base,
      score_adaptive = score_adapt,
      score_baseline = score_base
    )
  }
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

find_clusters_with_cap <- function(seu_obj,
                                   cluster_name,
                                   resolution,
                                   max_clusters = 5,
                                   min_resolution = 0.02,
                                   step = 0.02,
                                   seed = 1,
                                   log_prefix = "") {
  res <- resolution
  n_clusters <- Inf
  clusters <- NULL
  repeat {
    seu_obj <- FindClusters(
      seu_obj,
      resolution = res,
      algorithm = 4,
      cluster.name = cluster_name,
      random.seed = max(seed, 1),
      verbose = FALSE
    )
    clusters <- seu_obj[[cluster_name]][, 1]
    n_clusters <- length(unique(clusters))
    if (n_clusters <= max_clusters || res <= min_resolution) break
    next_res <- max(res - step, min_resolution)
    message(glue("{log_prefix} clusters={n_clusters} > {max_clusters}; resolution {res} -> {next_res}"))
    res <- next_res
  }
  list(seu = seu_obj, clusters = clusters, resolution = res, n_clusters = n_clusters)
}

compute_neighbors <- function(coords, query_cells, r) {
  if (!all(query_cells %in% rownames(coords))) {
    query_cells <- intersect(query_cells, rownames(coords))
  }
  if (length(query_cells) == 0) return(list())
  Q <- coords[query_cells, , drop = FALSE]
  nbrs <- dbscan::frNN(coords, eps = r, query = Q, sort = FALSE)
  neighbors <- lapply(nbrs$id, function(idx) rownames(coords)[idx])
  names(neighbors) <- rownames(Q)
  neighbors
}

collect_levels_from_outputs <- function(output_dir, sample_ids, k_value, use_pbapply = TRUE) {
  apply_fun <- if (use_pbapply && requireNamespace("pbapply", quietly = TRUE)) {
    pbapply::pblapply
  } else {
    lapply
  }
  levels_list <- apply_fun(sample_ids, function(sample_id) {
    path <- file.path(output_dir, sample_id, glue("niche_average_raw_k{k_value}.tsv"))
    if (!file.exists(path)) return(character(0))
    df <- read.delim(path, row.names = 1, check.names = FALSE)
    colnames(df)
  })
  unique(unlist(levels_list))
}

collect_levels_from_qs <- function(files, anno_col, use_pbapply = TRUE) {
  apply_fun <- if (use_pbapply && requireNamespace("pbapply", quietly = TRUE)) {
    pbapply::pblapply
  } else {
    lapply
  }
  levels_list <- apply_fun(files, function(file) {
    obj <- qs::qread(file)
    if (!anno_col %in% colnames(obj@meta.data)) {
      rm(obj)
      return(character(0))
    }
    vals <- unique(as.character(obj@meta.data[[anno_col]]))
    vals <- vals[!is.na(vals)]
    rm(obj)
    vals
  })
  unique(unlist(levels_list))
}

order_annotation_levels <- function(levels, canonical_levels = CANONICAL_COARSE) {
  levels <- unique(as.character(levels))
  levels <- levels[!is.na(levels)]
  extras <- levels[!levels %in% canonical_levels]
  extras <- extras[order(tolower(extras))]
  unique(c(canonical_levels, extras))
}

plot_r_scan <- function(scan_df, out_file, r_choice = NA_real_, reason = NULL, r_baseline = NA_real_, r_adaptive = NA_real_) {
  p <- ggplot(scan_df, aes(r, median_neighbors)) +
    geom_line(linewidth = 0.8, color = "#2C3E50") +
    geom_point(size = 1.6, color = "#2C3E50") +
    theme_minimal(base_size = 11) +
    labs(x = "r (um)", y = "Median neighbor count")
  if (!is.na(r_baseline)) {
    p <- p + geom_vline(xintercept = r_baseline, linetype = "dotted", color = "#7F8C8D")
  }
  if (!is.na(r_adaptive)) {
    p <- p + geom_vline(xintercept = r_adaptive, linetype = "dotdash", color = "#2C3E50")
  }
  if (!is.na(r_choice)) {
    p <- p +
      geom_vline(xintercept = r_choice, linetype = "dashed", color = "#E74C3C") +
      annotate("text", x = r_choice, y = max(scan_df$median_neighbors, na.rm = TRUE),
               label = paste0("r=", r_choice, if (!is.null(reason)) paste0(" (", reason, ")") else ""),
               vjust = -0.5, hjust = 0.5, size = 3.2, color = "#E74C3C")
  }
  ggsave(out_file, p, width = 5, height = 4, device = grDevices::cairo_pdf)
}

downsample_niche_proportional <- function(df, max_total) {
  # Downsample while preserving niche proportions.
  n_total <- nrow(df)
  if (n_total <= max_total) return(df)

  counts <- table(df$niche)
  types <- names(counts)
  n_types <- length(types)

  if (max_total < n_types) {
    set.seed(RNG_SEED)
    return(df[sample(seq_len(n_total), max_total), , drop = FALSE])
  }

  props <- counts / sum(counts)
  raw <- props * max_total
  target <- floor(raw)
  target[target == 0] <- 1

  diff <- max_total - sum(target)
  if (diff > 0) {
    frac <- raw - floor(raw)
    order_idx <- order(frac, decreasing = TRUE)
    for (idx in order_idx) {
      if (diff == 0) break
      if (target[idx] < counts[idx]) {
        target[idx] <- target[idx] + 1
        diff <- diff - 1
      }
    }
  } else if (diff < 0) {
    order_idx <- order(target, decreasing = TRUE)
    for (idx in order_idx) {
      if (diff == 0) break
      if (target[idx] > 1) {
        target[idx] <- target[idx] - 1
        diff <- diff + 1
      }
    }
  }
  names(target) <- types

  set.seed(RNG_SEED)
  do.call(
    rbind,
    lapply(types, function(tp) {
      df_tp <- df[df$niche == tp, , drop = FALSE]
      n_keep <- min(target[tp], nrow(df_tp))
      if (n_keep <= 0) return(NULL)
      df_tp[sample(seq_len(nrow(df_tp)), n_keep), , drop = FALSE]
    })
  )
}

plot_niche_pdf <- function(seu,
                           caf_cells,
                           niche_cells,
                           out_file,
                           annotation_col,
                           max_all_cells = 200000,
                           max_niche_cells = 200000,
                           title = NULL,
                           level_order = NULL) {
  meta <- seu@meta.data
  if (!all(c("Center_X", "Center_Y") %in% colnames(meta))) {
    warning("Center_X/Center_Y not found; skip spatial plot.")
    return(invisible(NULL))
  }

  df_all <- data.frame(
    cell_id = rownames(meta),
    x = meta$Center_X * COORD_SCALE,
    y = meta$Center_Y * COORD_SCALE,
    niche = if (annotation_col %in% colnames(meta)) {
      as.character(meta[[annotation_col]])
    } else {
      NA_character_
    }
  )
  if (nrow(df_all) == 0) {
    warning("No cells available for plotting; skip.")
    return(invisible(NULL))
  }

  df_niche <- df_all[df_all$cell_id %in% niche_cells, , drop = FALSE]
  if (nrow(df_niche) == 0) {
    warning("No niche cells to plot; skip.")
    return(invisible(NULL))
  }
  if (nrow(df_all) > max_all_cells) {
    set.seed(RNG_SEED)
    df_all <- df_all[sample(seq_len(nrow(df_all)), max_all_cells), , drop = FALSE]
  }
  if (nrow(df_niche) > max_niche_cells) {
    df_niche <- downsample_niche_proportional(df_niche, max_niche_cells)
  }

  level_vals <- unique(as.character(df_niche$niche))
  level_vals <- level_vals[!is.na(level_vals)]
  if (length(level_vals) > 0) {
    level_order_use <- level_order
    if (is.null(level_order_use) || length(level_order_use) == 0) {
      level_order_use <- CANONICAL_COARSE
    }
    df_niche$niche <- factor(df_niche$niche, levels = level_order_use)
    df_niche <- df_niche[order(df_niche$niche), , drop = FALSE]
  } else {
    level_order_use <- character(0)
  }

  caf_df <- df_all[df_all$cell_id %in% caf_cells, , drop = FALSE]

  # Panel 1: CAF anchors on full slice.
  p_all <- ggplot() +
    geom_point(
      data = df_all,
      aes(x = x, y = y),
      color = "grey85",
      size = 0.01,
      alpha = 0.25
    ) +
    geom_point(
      data = caf_df,
      aes(x = x, y = y),
      color = "#E74C3C",
      size = 0.01,
      alpha = 0.9
    ) +
    coord_equal() +
    theme_void(base_size = 11) +
    labs(title = paste0(title, " - CAF niche anchors"))

  # Panel 2: niche cells on full slice, colored by annotation.
  p_niche <- ggplot() +
    geom_point(
      data = df_all,
      aes(x = x, y = y),
      color = "grey85",
      size = 0.01,
      alpha = 0.35
    ) +
    geom_point(
      data = df_niche,
      aes(x = x, y = y, color = niche),
      shape = 16,
      size = 0.08,
      alpha = 0.6
    ) +
    coord_equal() +
    theme_void(base_size = 11) +
    theme(
      legend.position = "right",
      legend.key.size = grid::unit(0.6, "cm"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11)
    ) +
    labs(title = paste0(title, " - niche cells on slice"), color = annotation_col) +
    guides(color = guide_legend(override.aes = list(size = 3.2, alpha = 1, shape = 16)))

  if (!requireNamespace("ggsci", quietly = TRUE)) {
    stop("Package 'ggsci' is required for D3 palettes.")
  }
  p_niche <- p_niche + ggsci::scale_color_d3(
    "category20",
    limits = level_order_use,
    breaks = level_order_use,
    drop = FALSE
  )

  if (requireNamespace("cowplot", quietly = TRUE)) {
    p <- cowplot::plot_grid(p_all, p_niche, ncol = 2)
    ggsave(out_file, p, width = 14, height = 7, device = grDevices::cairo_pdf)
  } else {
    grDevices::cairo_pdf(out_file, width = 14, height = 7)
    print(p_all)
    print(p_niche)
    grDevices::dev.off()
  }
}

plot_heatmap <- function(mat, out_file, title = NULL, level_order = NULL) {
  df <- as.data.frame(mat)
  df$cluster <- rownames(df)
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required for plotting heatmaps.")
  }
  df_long <- tidyr::pivot_longer(df, cols = -cluster, names_to = "cell_type", values_to = "value")
  level_order_use <- level_order
  if (is.null(level_order_use) || length(level_order_use) == 0) {
    level_order_use <- CANONICAL_COARSE
  }
  df_long$cell_type <- factor(df_long$cell_type, levels = level_order_use)

  p <- ggplot(df_long, aes(cell_type, cluster, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#1F78B4") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(x = NULL, y = NULL, title = title)

  ggsave(out_file, p, width = 8, height = 5, device = grDevices::cairo_pdf)
}

plot_stacked_bar <- function(mat, out_file, title = NULL, level_order = NULL) {
  df <- as.data.frame(mat)
  df$cluster <- rownames(df)
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required for plotting stacked bars.")
  }
  df_long <- tidyr::pivot_longer(df, cols = -cluster, names_to = "cell_type", values_to = "fraction")
  level_order_use <- level_order
  if (is.null(level_order_use) || length(level_order_use) == 0) {
    level_order_use <- CANONICAL_COARSE
  }
  df_long$cell_type <- factor(df_long$cell_type, levels = level_order_use)

  p <- ggplot(df_long, aes(cluster, fraction, fill = cell_type)) +
    geom_col(position = "stack", width = 0.8) +
    theme_minimal(base_size = 15) +
    theme(axis.text = element_text(color = "black")) +
    labs(x = NULL, y = "Fraction of cell types", title = title) +
    Seurat::RotatedAxis()

  if (!requireNamespace("ggsci", quietly = TRUE)) {
    stop("Package 'ggsci' is required for D3 palettes.")
  }
  p <- p + ggsci::scale_fill_d3(
    "category20",
    limits = level_order_use,
    breaks = level_order_use,
    drop = FALSE
  )

  ggsave(out_file, p, width = 8, height = 5, device = grDevices::cairo_pdf)
}

argv <- commandArgs(trailingOnly = TRUE)
opts <- parse_args(argv)

defaults <- list(
  input = file.path(root_dir, "data"),
  output = file.path(root_dir, "data", "caf_niche_sectionnorm"),
  pattern = "\\.qs$",
  samples = "",
  anno_coarse = "annotation_coarse",
  anno_fine = "annotation_fine",
  coarse_schema = "myeloid_refined",
  refine_myeloid = FALSE,
  myeloid_fine_cols = "annotation_fine,cell_type_lvl2,cell_type_lvl1,annotation",
  caf_label = "CAF",
  reduction = "spatial",
  r_min = 40,
  r_max = 80,
  r_step = 5,
  r_fixed = 80,
  stable_delta = 0.1,
  stable_window = 3L,
  min_neighbors = 30L,
  target_neighbors = NA_real_,
  max_query_cells = 5000L,
  use_pbapply = TRUE,
  k_list = "5,6",
  use_section_normalization = TRUE,
  log_enrich_alpha = 0.5,
  cluster_resolution = 0.1,
  max_caf_clusters = 5L,
  min_cluster_resolution = 0.02,
  cluster_resolution_step = 0.02,
  nmf_cores = 5L,
  seed = 1024L,
  max_plot_cells = 200000L,
  max_plot_all_cells = NA_integer_,
  detail_fine = FALSE,
  plot_only = FALSE
)
defaults_display <- defaults
defaults_display$max_plot_all_cells <- "same as max_plot_cells"

if (!is.null(opts$help) || !is.null(opts$h)) {
  print_help(defaults_display)
  quit(status = 0)
}

# Resolve CLI options and validate.
cfg <- list(
  input = get_opt(opts, "input", defaults$input),
  output = get_opt(opts, "output", defaults$output),
  pattern = get_opt(opts, "pattern", defaults$pattern),
  samples = get_opt(opts, "samples", defaults$samples),
  anno_coarse = get_opt(opts, "anno_coarse", defaults$anno_coarse),
  anno_fine = get_opt(opts, "anno_fine", defaults$anno_fine),
  coarse_schema = get_opt(opts, "coarse_schema", defaults$coarse_schema),
  refine_myeloid = get_opt(opts, "refine_myeloid", defaults$refine_myeloid, to_bool),
  myeloid_fine_cols = get_opt(opts, "myeloid_fine_cols", defaults$myeloid_fine_cols),
  caf_label = get_opt(opts, "caf_label", defaults$caf_label),
  reduction = get_opt(opts, "reduction", defaults$reduction),
  r_min = get_opt(opts, "r_min", defaults$r_min, as.numeric),
  r_max = get_opt(opts, "r_max", defaults$r_max, as.numeric),
  r_step = get_opt(opts, "r_step", defaults$r_step, as.numeric),
  r_fixed = get_opt(opts, "r_fixed", defaults$r_fixed, as.numeric),
  stable_delta = get_opt(opts, "stable_delta", defaults$stable_delta, as.numeric),
  stable_window = get_opt(opts, "stable_window", defaults$stable_window, as.integer),
  min_neighbors = get_opt(opts, "min_neighbors", defaults$min_neighbors, as.integer),
  target_neighbors = get_opt(opts, "target_neighbors", defaults$target_neighbors, as.numeric),
  max_query_cells = get_opt(opts, "max_query_cells", defaults$max_query_cells, as.integer),
  use_pbapply = get_opt(opts, "use_pbapply", defaults$use_pbapply, to_bool),
  k_list = get_opt(opts, "k_list", defaults$k_list),
  use_section_normalization = get_opt(opts, "use_section_normalization", defaults$use_section_normalization, to_bool),
  log_enrich_alpha = get_opt(opts, "log_enrich_alpha", defaults$log_enrich_alpha, as.numeric),
  cluster_resolution = get_opt(opts, "cluster_resolution", defaults$cluster_resolution, as.numeric),
  max_caf_clusters = get_opt(opts, "max_caf_clusters", defaults$max_caf_clusters, as.integer),
  min_cluster_resolution = get_opt(opts, "min_cluster_resolution", defaults$min_cluster_resolution, as.numeric),
  cluster_resolution_step = get_opt(opts, "cluster_resolution_step", defaults$cluster_resolution_step, as.numeric),
  nmf_cores = get_opt(opts, "nmf_cores", defaults$nmf_cores, as.integer),
  seed = get_opt(opts, "seed", defaults$seed, as.integer),
  max_plot_cells = get_opt(opts, "max_plot_cells", defaults$max_plot_cells, as.integer),
  max_plot_all_cells = get_opt(opts, "max_plot_all_cells", defaults$max_plot_all_cells, as.integer),
  detail_fine = get_opt(opts, "detail_fine", defaults$detail_fine, to_bool),
  plot_only = get_opt(opts, "plot_only", defaults$plot_only, to_bool)
)

if (is.na(cfg$max_plot_all_cells)) {
  cfg$max_plot_all_cells <- cfg$max_plot_cells
}

cfg$coarse_schema <- tolower(cfg$coarse_schema)
if (cfg$refine_myeloid && cfg$coarse_schema == "base") {
  cfg$coarse_schema <- "myeloid_refined"
}
if (!cfg$coarse_schema %in% c("base", "myeloid_refined")) {
  stop(glue("Unknown coarse_schema: {cfg$coarse_schema}. Use base or myeloid_refined."))
}
CANONICAL_COARSE <- get_canonical_coarse(cfg$coarse_schema, order = "analysis")
cfg$myeloid_fine_cols <- parse_csv_list(cfg$myeloid_fine_cols)
if (length(cfg$myeloid_fine_cols) == 0) {
  cfg$myeloid_fine_cols <- default_myeloid_fine_cols()
}

if (cfg$r_min < 40) {
  warning(glue("r_min {cfg$r_min} < 40; reset to 40"))
  cfg$r_min <- 40
}
if (cfg$r_max > 100) {
  warning(glue("r_max {cfg$r_max} > 100; reset to 100"))
  cfg$r_max <- 100
}
if (cfg$r_max < cfg$r_min) {
  stop(glue("r_max {cfg$r_max} < r_min {cfg$r_min}"))
}

cfg$k_list <- as.integer(trimws(strsplit(cfg$k_list, ",")[[1]]))
if (anyNA(cfg$k_list)) {
  stop("k_list must be a comma-separated list of integers.")
}

RNG_SEED <- cfg$seed

input_files <- list.files(cfg$input, pattern = cfg$pattern, full.names = TRUE)
if (cfg$samples != "") {
  sample_list <- strsplit(cfg$samples, ",")[[1]]
  input_files <- input_files[basename(input_files) %in% paste0(sample_list, ".qs")]
}
if (length(input_files) == 0) {
  stop("No input files found.")
}

dir.create(cfg$output, recursive = TRUE, showWarnings = FALSE)

sample_ids <- sub("\\.qs$", "", basename(input_files))
global_level_order <- CANONICAL_COARSE

summary_rows <- list()
sample_times <- numeric(0)

# Per-sample processing: r scan -> niche -> NMF -> plots.
for (idx in seq_along(input_files)) {
  file <- input_files[idx]
  sample_id <- sub("\\.qs$", "", basename(file))
  sample_start <- Sys.time()
  message(glue("[{idx}/{length(input_files)}] {sample_id}: start"))

  seu <- qs::qread(file)
  if (!cfg$anno_coarse %in% colnames(seu@meta.data)) {
    warning(glue("{sample_id}: missing {cfg$anno_coarse}; skip."))
    next
  }
  if (!cfg$reduction %in% names(seu@reductions)) {
    warning(glue("{sample_id}: reduction {cfg$reduction} missing; skip."))
    next
  }

  if (cfg$refine_myeloid) {
    seu@meta.data[[cfg$anno_coarse]] <- refine_myeloid_in_coarse(
      seu@meta.data,
      anno_coarse = cfg$anno_coarse,
      fine_cols = cfg$myeloid_fine_cols
    )
  }
  seu <- enforce_canonical_levels(seu, cfg$anno_coarse, sample_id)

  message(glue("{sample_id}: scale coords by {COORD_SCALE}"))
  coords_scaled <- Seurat::Embeddings(seu, reduction = cfg$reduction) * COORD_SCALE
  seu[[cfg$reduction]]@cell.embeddings <- coords_scaled

  caf_cells <- rownames(seu@meta.data)[seu@meta.data[[cfg$anno_coarse]] == cfg$caf_label]
  if (length(caf_cells) == 0) {
    warning(glue("{sample_id}: no CAF cells; skip."))
    next
  }

  coords <- Seurat::Embeddings(seu, reduction = cfg$reduction)
  r_values <- seq(cfg$r_min, cfg$r_max, by = cfg$r_step)
  sample_dir <- file.path(cfg$output, sample_id)

  if (cfg$plot_only) {
    message(glue("{sample_id}: plot_only mode"))
    if (!dir.exists(sample_dir)) {
      warning(glue("{sample_id}: output dir missing; skip."))
      next
    }
    choice_path <- file.path(sample_dir, "r_choice.tsv")
    if (file.exists(choice_path)) {
      choice_df <- read.delim(choice_path, stringsAsFactors = FALSE)
      r_choice <- choice_df$r_choice[1]
    } else if (!is.na(cfg$r_fixed)) {
      r_choice <- cfg$r_fixed
    } else {
      warning(glue("{sample_id}: r_choice.tsv missing and r_fixed is NA; skip."))
      next
    }
    message(glue("{sample_id}: compute neighbors (r={r_choice})"))
    t_neighbors <- system.time({
      neighbors <- compute_neighbors(coords, caf_cells, r_choice)
    })
    message(glue("{sample_id}: neighbors ready in {round(t_neighbors['elapsed'], 1)}s"))
    neighbor_counts <- vapply(neighbors, length, integer(1))
    if (length(neighbor_counts) > 0) {
      message(glue("{sample_id}: median CAF neighbors @ r={r_choice} is {median(neighbor_counts)}"))
    }

    for (k in cfg$k_list) {
      cluster_path <- file.path(sample_dir, glue("caf_clusters_k{k}.tsv"))
      if (file.exists(cluster_path)) {
        cluster_df <- read.delim(cluster_path, stringsAsFactors = FALSE)
        cluster_vec <- factor(cluster_df$caf_cluster)
        names(cluster_vec) <- cluster_df$cell_id

        cluster_ids <- sort(unique(as.character(cluster_vec)))
        for (cluster_id in cluster_ids) {
          caf_cluster_cells <- names(cluster_vec)[cluster_vec == cluster_id]
          caf_cluster_cells <- intersect(caf_cluster_cells, names(neighbors))
          niche_cells <- unique(unlist(neighbors[caf_cluster_cells], use.names = FALSE))
          if (length(niche_cells) == 0) {
            next
          }
          message(glue("{sample_id}: niche plot {cluster_id} (cells={length(niche_cells)})"))
          plot_niche_pdf(
            seu,
            caf_cells = caf_cluster_cells,
            niche_cells = niche_cells,
            out_file = file.path(sample_dir, glue("caf_niche_spatial_k{k}_{cluster_id}.pdf")),
            annotation_col = cfg$anno_coarse,
            max_all_cells = cfg$max_plot_all_cells,
            max_niche_cells = cfg$max_plot_cells,
            title = glue("{sample_id} {cluster_id} niche (coarse)"),
            level_order = global_level_order
          )
        }
      } else {
        warning(glue("{sample_id}: {basename(cluster_path)} missing; skip spatial plots for k={k}."))
      }

      raw_path <- file.path(sample_dir, glue("niche_average_raw_k{k}.tsv"))
      if (file.exists(raw_path)) {
        niche_avg_raw <- read.delim(raw_path, row.names = 1, check.names = FALSE)
        niche_avg_raw <- align_canonical_columns(niche_avg_raw, sample_id, "niche_average_raw")
        plot_stacked_bar(
          niche_avg_raw,
          out_file = file.path(sample_dir, glue("niche_bar_raw_k{k}.pdf")),
          title = glue("{sample_id} CAF niche (raw, k={k})"),
          level_order = global_level_order
        )
        plot_heatmap(
          niche_avg_raw,
          out_file = file.path(sample_dir, glue("niche_heatmap_raw_k{k}.pdf")),
          title = glue("{sample_id} CAF niche (raw, k={k})"),
          level_order = global_level_order
        )
      } else {
        warning(glue("{sample_id}: niche_average_raw_k{k}.tsv missing; skip raw plots."))
      }

      corrected_path <- file.path(sample_dir, glue("niche_average_corrected_k{k}.tsv"))
      if (file.exists(corrected_path)) {
        niche_avg_corrected <- read.delim(corrected_path, row.names = 1, check.names = FALSE)
        niche_avg_corrected <- align_canonical_columns(niche_avg_corrected, sample_id, "niche_average_corrected")
        plot_stacked_bar(
          niche_avg_corrected,
          out_file = file.path(sample_dir, glue("niche_bar_corrected_k{k}.pdf")),
          title = glue("{sample_id} CAF niche (corrected, k={k})"),
          level_order = global_level_order
        )
        plot_heatmap(
          niche_avg_corrected,
          out_file = file.path(sample_dir, glue("niche_heatmap_corrected_k{k}.pdf")),
          title = glue("{sample_id} CAF niche (corrected, k={k})"),
          level_order = global_level_order
        )
      } else {
        warning(glue("{sample_id}: niche_average_corrected_k{k}.tsv missing; skip corrected plots."))
      }
    }

    rm(seu)
    gc()

    sample_elapsed <- as.numeric(difftime(Sys.time(), sample_start, units = "mins"))
    sample_times <- c(sample_times, sample_elapsed)
    avg_time <- mean(sample_times)
    remaining <- avg_time * (length(input_files) - idx)
    message(glue("[{idx}/{length(input_files)}] {sample_id}: done in {round(sample_elapsed, 2)} min, ETA {round(remaining, 2)} min"))
    next
  }

  if (!is.na(cfg$r_fixed)) {
    choice <- list(r = cfg$r_fixed, reason = "fixed")
    r_choice <- cfg$r_fixed
    scan_df <- NULL
    message(glue("{sample_id}: r fixed to {r_choice}"))
  } else {
    message(glue("{sample_id}: r scan ({length(r_values)} values, max_query_cells={cfg$max_query_cells})"))
    t_scan <- system.time({
      scan_df <- scan_neighbor_counts(
        coords = coords,
        query_cells = caf_cells,
        r_values = r_values,
        max_query_cells = cfg$max_query_cells,
        use_pbapply = cfg$use_pbapply
      )
    })
    adaptive_choice <- choose_stable_r(
      scan_df,
      stable_delta = cfg$stable_delta,
      stable_window = cfg$stable_window,
      min_neighbors = cfg$min_neighbors,
      target_neighbors = cfg$target_neighbors
    )
    choice <- select_r_with_baseline(
      scan_df,
      adaptive_choice,
      r_baseline = R_BASELINE,
      min_neighbors = cfg$min_neighbors,
      target_neighbors = cfg$target_neighbors
    )
    r_choice <- choice$r
    message(glue("{sample_id}: r adaptive={choice$r_adaptive}, r baseline={choice$r_baseline}, selected={r_choice} ({choice$reason}), scan {round(t_scan['elapsed'], 1)}s"))
  }

  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    data.frame(
      r_choice = r_choice,
      r_reason = choice$reason,
      r_adaptive = choice$r_adaptive %||% NA_real_,
      r_baseline = choice$r_baseline %||% NA_real_,
      median_adaptive = choice$median_adaptive %||% NA_real_,
      median_baseline = choice$median_baseline %||% NA_real_,
      score_adaptive = choice$score_adaptive %||% NA_real_,
      score_baseline = choice$score_baseline %||% NA_real_
    ),
    file.path(sample_dir, "r_choice.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  if (is.na(cfg$r_fixed)) {
    utils::write.table(scan_df, file.path(sample_dir, "r_scan.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    plot_r_scan(
      scan_df,
      file.path(sample_dir, "r_scan.pdf"),
      r_choice = r_choice,
      reason = choice$reason,
      r_baseline = choice$r_baseline,
      r_adaptive = choice$r_adaptive
    )
  }

  message(glue("{sample_id}: niche matrix (coarse, r={r_choice})"))
  t_niche <- system.time({
    niche_res <- CalNichMatrix(
      seu,
      query.cells = caf_cells,
      r = r_choice,
      reduction = cfg$reduction,
      anno_col = cfg$anno_coarse,
      use_pbapply = cfg$use_pbapply,
      return_neighbors = TRUE
    )
  })
  message(glue("{sample_id}: niche matrix done in {round(t_niche['elapsed'], 1)}s"))

  niche_mat <- align_canonical_columns(niche_res$niche.mat, sample_id, "niche_celltype_count")
  niche_mat_norm <- align_canonical_columns(niche_res$niche.mat.norm, sample_id, "niche_celltype_fraction")
  section_prop <- compute_section_prop(seu@meta.data, cfg$anno_coarse, sample_id)
  log2_res <- compute_log2_enrichment(niche_mat, section_prop, alpha = cfg$log_enrich_alpha)
  niche_mat_corrected <- log2_res$log2_enrich
  niche_mat_expected <- log2_res$expected
  neighbor_counts <- log2_res$total
  if (length(neighbor_counts) > 0) {
    message(glue("{sample_id}: median CAF neighbors @ r={r_choice} is {median(neighbor_counts)}"))
  }

  utils::write.table(
    niche_mat,
    file.path(sample_dir, "niche_celltype_count.tsv"),
    sep = "\t",
    quote = FALSE
  )
  utils::write.table(
    niche_mat_norm,
    file.path(sample_dir, "niche_celltype_fraction.tsv"),
    sep = "\t",
    quote = FALSE
  )
  utils::write.table(
    niche_mat_corrected,
    file.path(sample_dir, "niche_celltype_corrected.tsv"),
    sep = "\t",
    quote = FALSE
  )
  utils::write.table(
    niche_mat_expected,
    file.path(sample_dir, "niche_celltype_expected.tsv"),
    sep = "\t",
    quote = FALSE
  )
  utils::write.table(
    data.frame(cell_id = names(neighbor_counts), neighbor_total = as.numeric(neighbor_counts)),
    file.path(sample_dir, "niche_neighbor_total.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  support_summary <- summarize_support(niche_mat, niche_mat_expected, neighbor_counts)
  utils::write.table(
    support_summary,
    file.path(sample_dir, "niche_support_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    data.frame(
      cell_type = CANONICAL_COARSE,
      fraction = as.numeric(section_prop)
    ),
    file.path(sample_dir, "section_celltype_fraction.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  fine_norm <- NULL
  if (cfg$detail_fine && cfg$anno_fine %in% colnames(seu@meta.data)) {
    message(glue("{sample_id}: niche matrix (fine, r={r_choice})"))
    t_fine <- system.time({
      fine_res <- CalNichMatrix(
        seu,
        query.cells = caf_cells,
        r = r_choice,
        reduction = cfg$reduction,
        anno_col = cfg$anno_fine,
        use_pbapply = cfg$use_pbapply,
        return_neighbors = FALSE
      )
    })
    fine_norm <- fine_res$niche.mat.norm
    message(glue("{sample_id}: fine niche matrix done in {round(t_fine['elapsed'], 1)}s"))
  }

  seu_caf <- subset(seu, cells = caf_cells)

  nmf_input <- if (cfg$use_section_normalization) niche_mat_corrected else niche_mat_norm
  if (cfg$use_section_normalization) {
    min_val <- min(nmf_input, na.rm = TRUE)
    if (is.finite(min_val) && min_val < 0) {
      shift_val <- -min_val
      nmf_input <- nmf_input + shift_val
      message(glue("{sample_id}: shift log2 enrichment by {round(shift_val, 3)} for NMF"))
    }
  }
  n_caf_total <- ncol(seu_caf)
  for (k in cfg$k_list) {
    message(glue("{sample_id}: NMF k={k}"))
    t_nmf <- system.time({
      nmf_res <- RunNMF(nmf_input, k = k, cores = cfg$nmf_cores, seed = cfg$seed)
    })
    message(glue("{sample_id}: NMF k={k} done in {round(t_nmf['elapsed'], 1)}s"))
    W <- nmf_res$W
    H <- nmf_res$H

    reduction_name <- paste0("nmf_k", k)
    key_name <- paste0("NMF", k, "_")
    seu_caf[[reduction_name]] <- CreateDimReducObject(
      embeddings = W[rownames(seu_caf@meta.data), , drop = FALSE],
      assay = "RNA",
      key = key_name
    )

    cluster_name <- paste0("CAF.sub_k", k)
    if (n_caf_total < 2) {
      warning(glue("{sample_id}: only {n_caf_total} CAF cell(s); skip Leiden and assign single cluster for k={k}."))
      cluster_vec <- factor(rep("CAF-1", n_caf_total))
      names(cluster_vec) <- rownames(seu_caf@meta.data)
      seu_caf[[cluster_name]] <- cluster_vec
      cluster_res <- list(resolution = cfg$cluster_resolution, n_clusters = 1L)
    } else {
      seu_caf <- FindNeighbors(seu_caf, reduction = reduction_name, dims = 1:k, verbose = FALSE)
      cluster_res <- find_clusters_with_cap(
        seu_caf,
        cluster_name = cluster_name,
        resolution = cfg$cluster_resolution,
        max_clusters = cfg$max_caf_clusters,
        min_resolution = cfg$min_cluster_resolution,
        step = cfg$cluster_resolution_step,
        seed = cfg$seed,
        log_prefix = glue("{sample_id}: k={k}")
      )
      seu_caf <- cluster_res$seu
      if (cluster_res$resolution != cfg$cluster_resolution) {
        message(glue("{sample_id}: k={k} resolution {cfg$cluster_resolution} -> {cluster_res$resolution} (clusters={cluster_res$n_clusters})"))
      }
      if (cluster_res$n_clusters > cfg$max_caf_clusters) {
        warning(glue("{sample_id}: k={k} clusters={cluster_res$n_clusters} > max {cfg$max_caf_clusters} at resolution {cluster_res$resolution}"))
      }
      cluster_vec <- factor(paste0("CAF-", cluster_res$clusters))
      seu_caf[[cluster_name]] <- cluster_vec
      names(cluster_vec) <- rownames(seu_caf@meta.data)
    }

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

    niche_avg_raw <- cluster_average(niche_mat_norm, cluster_vec)
    niche_avg_corrected <- cluster_average(niche_mat_corrected, cluster_vec)

    utils::write.table(
      niche_avg_raw,
      file.path(sample_dir, glue("niche_average_raw_k{k}.tsv")),
      sep = "\t",
      quote = FALSE
    )
    utils::write.table(
      niche_avg_corrected,
      file.path(sample_dir, glue("niche_average_corrected_k{k}.tsv")),
      sep = "\t",
      quote = FALSE
    )

    message(glue("{sample_id}: plots k={k}"))
    neighbors <- niche_res$neighbors
    if (!is.null(neighbors)) {
      if (!is.list(neighbors)) {
        neighbors <- as.list(as.data.frame(neighbors))
      }
      if (is.null(names(neighbors))) {
        names(neighbors) <- rownames(niche_mat)
      }
      cluster_ids <- sort(unique(as.character(cluster_vec)))
      for (cluster_id in cluster_ids) {
        caf_cluster_cells <- names(cluster_vec)[cluster_vec == cluster_id]
        caf_cluster_cells <- intersect(caf_cluster_cells, names(neighbors))
        niche_cells <- unique(unlist(neighbors[caf_cluster_cells], use.names = FALSE))
        if (length(niche_cells) == 0) {
          next
        }
        message(glue("{sample_id}: niche plot {cluster_id} (cells={length(niche_cells)})"))
        plot_niche_pdf(
          seu,
          caf_cells = caf_cluster_cells,
          niche_cells = niche_cells,
          out_file = file.path(sample_dir, glue("caf_niche_spatial_k{k}_{cluster_id}.pdf")),
          annotation_col = cfg$anno_coarse,
          max_all_cells = cfg$max_plot_all_cells,
          max_niche_cells = cfg$max_plot_cells,
          title = glue("{sample_id} {cluster_id} niche (coarse)"),
          level_order = global_level_order
        )
      }
    } else {
      warning(glue("{sample_id}: neighbors not available; skip niche spatial plots."))
    }

    plot_stacked_bar(
      niche_avg_raw,
      out_file = file.path(sample_dir, glue("niche_bar_raw_k{k}.pdf")),
      title = glue("{sample_id} CAF niche (raw, k={k})"),
      level_order = global_level_order
    )
    plot_heatmap(
      niche_avg_raw,
      out_file = file.path(sample_dir, glue("niche_heatmap_raw_k{k}.pdf")),
      title = glue("{sample_id} CAF niche (raw, k={k})"),
      level_order = global_level_order
    )
    plot_stacked_bar(
      niche_avg_corrected,
      out_file = file.path(sample_dir, glue("niche_bar_corrected_k{k}.pdf")),
      title = glue("{sample_id} CAF niche (corrected, k={k})"),
      level_order = global_level_order
    )
    plot_heatmap(
      niche_avg_corrected,
      out_file = file.path(sample_dir, glue("niche_heatmap_corrected_k{k}.pdf")),
      title = glue("{sample_id} CAF niche (corrected, k={k})"),
      level_order = global_level_order
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
      caf_cells = length(caf_cells),
      r_adaptive = choice$r_adaptive %||% NA_real_,
      r_baseline = choice$r_baseline %||% NA_real_,
      median_adaptive = choice$median_adaptive %||% NA_real_,
      median_baseline = choice$median_baseline %||% NA_real_
    )
  }

  rm(seu)
  gc()

  sample_elapsed <- as.numeric(difftime(Sys.time(), sample_start, units = "mins"))
  sample_times <- c(sample_times, sample_elapsed)
  avg_time <- mean(sample_times)
  remaining <- avg_time * (length(input_files) - idx)
  message(glue("[{idx}/{length(input_files)}] {sample_id}: done in {round(sample_elapsed, 2)} min, ETA {round(remaining, 2)} min"))
}

summary_df <- do.call(rbind, summary_rows)
if (!is.null(summary_df) && nrow(summary_df) > 0) {
  utils::write.table(
    summary_df,
    file.path(cfg$output, "summary_r_choice.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}
