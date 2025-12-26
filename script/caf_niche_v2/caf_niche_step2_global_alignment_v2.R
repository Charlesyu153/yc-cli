#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glue)
  library(ggplot2)
  library(pheatmap)
  library(ggsci)
})

`%||%` <- function(x, y) if (is.null(x) || x == "") y else x

parse_args <- function(argv) {
  args <- list()
  if (length(argv) %% 2 != 0) stop("Arguments must be provided as --key value pairs.")
  if (length(argv) == 0) return(args)
  keys <- argv[seq(1, length(argv), by = 2)]
  vals <- argv[seq(2, length(argv), by = 2)]
  keys <- sub("^--", "", keys)
  for (i in seq_along(keys)) args[[keys[i]]] <- vals[i]
  args
}

to_bool <- function(x) tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")

get_opt <- function(opts, key, default, cast = NULL) {
  val <- opts[[key]]
  if (is.null(val) || val == "") return(default)
  if (is.null(cast)) return(val)
  cast(val)
}

get_root_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  start <- if (length(file_arg) > 0) dirname(normalizePath(file_arg[1])) else getwd()
  d <- normalizePath(start)
  for (i in seq_len(10)) {
    if (file.exists(file.path(d, ".gitignore")) &&
      (dir.exists(file.path(d, ".git")) || file.exists(file.path(d, "AGENTS.md")))) {
      return(d)
    }
    parent <- normalizePath(file.path(d, ".."))
    if (parent == d) break
    d <- parent
  }
  normalizePath(getwd())
}

infer_group <- function(sample_id) {
  if (grepl("P$", sample_id)) return("P")
  if (grepl("[LR]$", sample_id)) return("NonP")
  NA_character_
}

align_profile_columns <- function(mat, canonical) {
  mat <- as.matrix(mat)
  if (is.null(colnames(mat))) stop("Profile matrix missing column names.")
  extras <- setdiff(colnames(mat), canonical)
  if (length(extras) > 0) {
    # Allow extra columns (e.g. removed levels like Myeloid_Mast) and drop them.
    mat <- mat[, setdiff(colnames(mat), extras), drop = FALSE]
  }
  missing <- setdiff(canonical, colnames(mat))
  if (length(missing) > 0) {
    zeros <- matrix(0, nrow = nrow(mat), ncol = length(missing))
    colnames(zeros) <- missing
    mat <- cbind(mat, zeros)
  }
  mat[, canonical, drop = FALSE]
}

transform_matrix <- function(x, transform, eps) {
  if (transform == "none") return(x)
  if (transform == "log2p") return(log2(x + eps))
  if (transform == "log1p") return(log1p(x))
  stop(glue("Unknown transform: {transform}"))
}

weighted_centroid <- function(x, w) {
  w <- as.numeric(w)
  w[!is.finite(w) | is.na(w) | w < 0] <- 0
  if (sum(w) <= 0) return(rep(NA_real_, ncol(x)))
  as.numeric(colSums(x * w) / sum(w))
}

within_dispersion <- function(x, labels, weights) {
  labs <- sort(unique(labels))
  disp <- numeric(0)
  for (g in labs) {
    idx <- which(labels == g)
    if (length(idx) == 0) next
    w <- weights[idx]
    cent <- weighted_centroid(x[idx, , drop = FALSE], w)
    if (anyNA(cent)) next
    d <- x[idx, , drop = FALSE] - matrix(cent, nrow = length(idx), ncol = ncol(x), byrow = TRUE)
    w_sum <- sum(w)
    if (w_sum <= 0) next
    disp_g <- sum(w * rowSums(d * d)) / w_sum
    disp <- c(disp, disp_g)
  }
  mean(disp, na.rm = TRUE)
}

subtype_group_stats <- function(labels, meta, weights, cov_min, pur_min, lambda_p, lambda_l) {
  meta$label <- labels
  meta$w <- weights
  meta <- meta[is.finite(meta$w) & !is.na(meta$w) & meta$w > 0, , drop = FALSE]
  labs <- sort(unique(meta$label))
  p_samples <- unique(meta$sample_id[meta$group == "P"])
  np_samples <- unique(meta$sample_id[meta$group == "NonP"])
  n_p <- max(1L, length(p_samples))
  n_np <- max(1L, length(np_samples))

  rows <- list()
  for (g in labs) {
    m <- meta[meta$label == g, , drop = FALSE]
    w_total <- sum(m$w)
    w_p <- sum(m$w[m$group == "P"])
    w_np <- sum(m$w[m$group == "NonP"])
    pur_p <- if (w_total > 0) w_p / w_total else NA_real_
    pur_np <- if (w_total > 0) w_np / w_total else NA_real_

    cov_p <- length(unique(m$sample_id[m$group == "P"])) / n_p
    cov_np <- length(unique(m$sample_id[m$group == "NonP"])) / n_np

    is_universal <- cov_p >= cov_min && cov_np >= cov_min
    is_p_only <- cov_p >= cov_min && cov_np < cov_min
    is_np_only <- cov_np >= cov_min && cov_p < cov_min
    is_valid <- is_universal || is_p_only || is_np_only

    pen_purity <- 0
    pen_leak <- 0
    if (is_p_only) {
      pen_purity <- max(0, pur_min - pur_p)
      pen_leak <- max(0, pur_np - (1 - pur_min))
    } else if (is_np_only) {
      pen_purity <- max(0, pur_min - pur_np)
      pen_leak <- max(0, pur_p - (1 - pur_min))
    }

    rows[[length(rows) + 1]] <- data.frame(
      subtype = g,
      w_total = w_total,
      cov_p = cov_p,
      cov_nonp = cov_np,
      pur_p = pur_p,
      pur_nonp = pur_np,
      is_valid = as.integer(is_valid),
      is_universal = as.integer(is_universal),
      is_p_only = as.integer(is_p_only),
      is_nonp_only = as.integer(is_np_only),
      pen_purity = pen_purity,
      pen_leak = pen_leak,
      penalty = (1 - as.integer(is_valid)) + lambda_p * pen_purity + lambda_l * pen_leak
    )
  }
  do.call(rbind, rows)
}

coverage_score <- function(labels, meta) {
  labs <- sort(unique(labels))
  n_samples <- length(unique(meta$sample_id))
  if (n_samples == 0 || length(labs) == 0) return(NA_real_)
  covs <- sapply(labs, function(g) length(unique(meta$sample_id[labels == g])) / n_samples)
  mean(covs)
}

max_sample_share <- function(labels, meta, weights) {
  out <- stats::setNames(rep(NA_real_, length(unique(labels))), sort(unique(labels)))
  for (g in names(out)) {
    idx <- which(labels == g)
    if (length(idx) == 0) next
    df <- data.frame(sample_id = meta$sample_id[idx], w = weights[idx])
    w_by <- stats::aggregate(w ~ sample_id, df, sum)
    w_tot <- sum(w_by$w)
    out[g] <- if (w_tot > 0) max(w_by$w) / w_tot else NA_real_
  }
  out
}

rotate_leaf_order_by_label <- function(hc, leaf_labels) {
  n <- length(hc$labels)
  if (length(leaf_labels) != n) stop("leaf_labels length must match hc$labels length.")
  merge <- hc$merge

  solve_node <- function(node) {
    if (node < 0) {
      leaf_idx <- -node
      lab <- as.character(leaf_labels[[leaf_idx]])
      return(list(order = leaf_idx, score = 0L, first = lab, last = lab))
    }
    left <- merge[node, 1]
    right <- merge[node, 2]
    a <- solve_node(left)
    b <- solve_node(right)

    score_ab <- a$score + b$score + as.integer(!identical(a$last, b$first))
    score_ba <- b$score + a$score + as.integer(!identical(b$last, a$first))

    if (score_ba < score_ab) {
      return(list(order = c(b$order, a$order), score = score_ba, first = b$first, last = a$last))
    }
    if (score_ab < score_ba) {
      return(list(order = c(a$order, b$order), score = score_ab, first = a$first, last = b$last))
    }
    list(order = c(a$order, b$order), score = score_ab, first = a$first, last = b$last)
  }

  res <- solve_node(n - 1L)
  hc$labels[res$order]
}

count_label_runs <- function(labels) {
  labels <- as.character(labels)
  labels <- labels[!is.na(labels) & labels != ""]
  if (length(labels) == 0) return(data.frame(label = character(0), runs = integer(0)))
  runs <- c(TRUE, labels[-1] != labels[-length(labels)])
  r <- table(labels[runs])
  data.frame(label = names(r), runs = as.integer(r), stringsAsFactors = FALSE)
}

merge_bad_subtypes <- function(x, labels, meta, weights, cov_min, pur_min, lambda_p, lambda_l,
                               p_max = 0.8, low_weight_q = 0.05, pen_soft_max = 0.25) {
  labels <- as.character(labels)
  iter <- 0L
  repeat {
    iter <- iter + 1L
    labs <- sort(unique(labels))
    if (length(labs) <= 1) break

    stats_df <- subtype_group_stats(labels, meta, weights, cov_min, pur_min, lambda_p, lambda_l)
    p_share <- max_sample_share(labels, meta, weights)
    stats_df$p_max_sample <- as.numeric(p_share[stats_df$subtype])
    # NOTE: With K<=5, a "bottom 5% subtype weight" rule degenerates into "always merge the smallest",
    # which can collapse everything into 1 subtype. Only enable this heuristic when subtype count is larger.
    w_thresh <- -Inf
    if (length(labs) >= 6 && low_weight_q > 0) {
      w_thresh <- stats::quantile(stats_df$w_total, probs = low_weight_q, na.rm = TRUE)
    }

    bad <- stats_df[
      stats_df$is_valid == 0 |
        stats_df$p_max_sample > p_max |
        stats_df$w_total <= w_thresh |
        stats_df$pen_purity > pen_soft_max |
        stats_df$pen_leak > pen_soft_max,
      ,
      drop = FALSE
    ]
    if (nrow(bad) == 0) break
    if (iter > 50) break

    bad <- bad[order(bad$is_valid, bad$penalty, bad$p_max_sample, decreasing = TRUE), , drop = FALSE]
    g_bad <- bad$subtype[1]
    others <- setdiff(labs, g_bad)
    if (length(others) == 0) break

    cent_bad <- weighted_centroid(x[labels == g_bad, , drop = FALSE], weights[labels == g_bad])
    best_g <- NA_character_
    best_sim <- -Inf
    for (g in others) {
      cent <- weighted_centroid(x[labels == g, , drop = FALSE], weights[labels == g])
      sim <- suppressWarnings(stats::cor(cent_bad, cent, use = "pairwise.complete.obs"))
      if (!is.finite(sim) || is.na(sim)) next
      if (sim > best_sim) {
        best_sim <- sim
        best_g <- g
      }
    }
    if (is.na(best_g)) break
    labels[labels == g_bad] <- best_g
  }
  factor(labels)
}

PARAM_HELP <- list(
  input_root = "Step1 output root directory (contains per-sample folders)",
  output_root = "Output root directory for Step2 results",
  samples = "Comma-separated sample IDs (blank = auto-detect)",
  profile = "Profile to use: corrected or raw",
  coarse_schema = "Coarse schema: myeloid_refined",
  transform = "Transform: log2p, log1p, none",
  eps = "Epsilon for log2p",
  corr_method = "Correlation method: pearson or spearman",
  hclust_method = "hclust method",
  k_candidates = "Comma-separated K candidates (default: 2,3,4,5)",
  cov_min = "Coverage threshold within group (default: 0.5)",
  pur_min = "Purity reference threshold for penalties (default: 0.7)",
  alpha = "Score weight for coverage",
  beta = "Score weight for group penalty",
  gamma = "Score weight for within-dispersion",
  lambda_p = "Penalty weight for purity gap",
  lambda_l = "Penalty weight for leakage",
  pen_soft_max = "Threshold to trigger merge based on pen_purity/pen_leak",
  p_max_sample = "Threshold to trigger merge based on max sample share",
  low_weight_q = "Quantile threshold to trigger merge based on low total weight",
  write_heatmaps = "Write heatmap PDFs"
)

print_help <- function(defaults) {
  cat("CAF niche Step2 v2: global alignment on Step1 selected local cluster profiles\n")
  for (key in names(PARAM_HELP)) {
    def <- defaults[[key]]
    def_str <- if (is.na(def) || is.null(def)) "NA" else as.character(def)
    cat(glue("  --{key} {PARAM_HELP[[key]]} (default: {def_str})\n"))
  }
}

argv <- commandArgs(trailingOnly = TRUE)
opts <- parse_args(argv)

defaults <- list(
  input_root = file.path(getwd(), "data", "caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2", "step1"),
  output_root = file.path(getwd(), "data", "caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2", "step2", "global_corrected"),
  samples = "",
  profile = "corrected",
  coarse_schema = "myeloid_refined",
  transform = "log2p",
  eps = 1e-6,
  corr_method = "pearson",
  hclust_method = "average",
  k_candidates = "2,3,4,5",
  cov_min = 0.5,
  pur_min = 0.7,
  alpha = 1.0,
  beta = 1.0,
  gamma = 0.5,
  lambda_p = 1.0,
  lambda_l = 1.0,
  pen_soft_max = 0.25,
  p_max_sample = 0.8,
  low_weight_q = 0.05,
  write_heatmaps = TRUE
)

if (!is.null(opts[["help"]]) || !is.null(opts[["h"]])) {
  print_help(defaults)
  quit(status = 0)
}

cfg <- list(
  input_root = get_opt(opts, "input_root", defaults$input_root),
  output_root = get_opt(opts, "output_root", defaults$output_root),
  samples = get_opt(opts, "samples", defaults$samples),
  profile = tolower(get_opt(opts, "profile", defaults$profile)),
  coarse_schema = tolower(get_opt(opts, "coarse_schema", defaults$coarse_schema)),
  transform = tolower(get_opt(opts, "transform", defaults$transform)),
  eps = get_opt(opts, "eps", defaults$eps, as.numeric),
  corr_method = tolower(get_opt(opts, "corr_method", defaults$corr_method)),
  hclust_method = get_opt(opts, "hclust_method", defaults$hclust_method),
  k_candidates = get_opt(opts, "k_candidates", defaults$k_candidates),
  cov_min = get_opt(opts, "cov_min", defaults$cov_min, as.numeric),
  pur_min = get_opt(opts, "pur_min", defaults$pur_min, as.numeric),
  alpha = get_opt(opts, "alpha", defaults$alpha, as.numeric),
  beta = get_opt(opts, "beta", defaults$beta, as.numeric),
  gamma = get_opt(opts, "gamma", defaults$gamma, as.numeric),
  lambda_p = get_opt(opts, "lambda_p", defaults$lambda_p, as.numeric),
  lambda_l = get_opt(opts, "lambda_l", defaults$lambda_l, as.numeric),
  pen_soft_max = get_opt(opts, "pen_soft_max", defaults$pen_soft_max, as.numeric),
  p_max_sample = get_opt(opts, "p_max_sample", defaults$p_max_sample, as.numeric),
  low_weight_q = get_opt(opts, "low_weight_q", defaults$low_weight_q, as.numeric),
  write_heatmaps = get_opt(opts, "write_heatmaps", defaults$write_heatmaps, to_bool)
)

if (!dir.exists(cfg$input_root)) stop(glue("input_root not found: {cfg$input_root}"))
dir.create(cfg$output_root, recursive = TRUE, showWarnings = FALSE)

if (!cfg$profile %in% c("corrected", "raw")) stop("profile must be corrected or raw.")
if (!cfg$transform %in% c("log2p", "log1p", "none")) stop("transform must be log2p/log1p/none.")
if (!cfg$corr_method %in% c("pearson", "spearman")) stop("corr_method must be pearson or spearman.")

cfg$k_candidates <- as.integer(trimws(strsplit(cfg$k_candidates, ",")[[1]]))
cfg$k_candidates <- cfg$k_candidates[is.finite(cfg$k_candidates) & !is.na(cfg$k_candidates)]
cfg$k_candidates <- sort(unique(cfg$k_candidates))
if (length(cfg$k_candidates) == 0) stop("k_candidates is empty.")
if (max(cfg$k_candidates) > 5) stop("k_candidates must be <= 5 for this v2 task.")

root_dir <- get_root_dir()
source(file.path(root_dir, "utils", "R", "coarse_schema.R"))
canonical <- get_canonical_coarse(cfg$coarse_schema, order = "analysis")

sample_dirs <- list.dirs(cfg$input_root, full.names = TRUE, recursive = FALSE)
sample_ids <- basename(sample_dirs)
if (cfg$samples != "") {
  keep <- trimws(strsplit(cfg$samples, ",")[[1]])
  sample_dirs <- sample_dirs[sample_ids %in% keep]
  sample_ids <- basename(sample_dirs)
}
sample_dirs <- sample_dirs[order(sample_ids)]
if (length(sample_dirs) == 0) stop("No sample directories found under input_root.")

profiles <- list()
meta_rows <- list()

for (sample_dir in sample_dirs) {
  sample_id <- basename(sample_dir)
  group <- infer_group(sample_id)
  k_path <- file.path(sample_dir, "k_selected.tsv")
  prof_path <- file.path(sample_dir, glue("niche_average_{cfg$profile}_selected.tsv"))
  summary_path <- file.path(sample_dir, "caf_cluster_summary_selected.tsv")
  if (!file.exists(k_path) || !file.exists(prof_path) || !file.exists(summary_path)) {
    stop(glue("{sample_id}: missing required Step1 outputs in {sample_dir}"))
  }
  k_sel <- read.delim(k_path, stringsAsFactors = FALSE)
  if (nrow(k_sel) != 1) stop(glue("{sample_id}: invalid k_selected.tsv rows: {nrow(k_sel)}"))
  stab <- as.numeric(k_sel$stab_combined[1])
  fit <- as.numeric(k_sel$fit_score[1])
  is_unstable <- as.integer(k_sel$is_unstable[1])

  prof <- read.delim(prof_path, check.names = FALSE)
  prof <- align_profile_columns(prof, canonical)
  if (is.null(rownames(prof))) stop(glue("{sample_id}: profile missing rownames (local clusters)."))

  summ <- read.delim(summary_path, stringsAsFactors = FALSE)
  if (!all(c("local_cluster", "n_cells") %in% colnames(summ))) {
    stop(glue("{sample_id}: caf_cluster_summary_selected.tsv missing local_cluster/n_cells."))
  }
  n_cells <- stats::setNames(as.integer(summ$n_cells), summ$local_cluster)

  for (lc in rownames(prof)) {
    gid <- paste0(sample_id, ":", lc)
    row <- prof[lc, , drop = FALSE]
    rownames(row) <- gid
    profiles[[gid]] <- row
    meta_rows[[gid]] <- data.frame(
      global_local_cluster_id = gid,
      sample_id = sample_id,
      group = group,
      local_cluster = lc,
      n_cells = n_cells[[lc]] %||% NA_integer_,
      stab_combined = stab,
      fit_score = fit,
      is_unstable = is_unstable,
      stringsAsFactors = FALSE
    )
  }
}

X <- do.call(rbind, profiles)
meta <- do.call(rbind, meta_rows)
rownames(meta) <- meta$global_local_cluster_id

meta$n_cells[is.na(meta$n_cells) | meta$n_cells < 0] <- 0L
meta$stab_combined[!is.finite(meta$stab_combined) | is.na(meta$stab_combined)] <- 0
meta$is_unstable[is.na(meta$is_unstable)] <- 1L

weights <- log1p(meta$n_cells) * meta$stab_combined * (1 - 0.5 * meta$is_unstable)
weights[!is.finite(weights) | is.na(weights) | weights < 0] <- 0

utils::write.table(X, file.path(cfg$output_root, glue("global_X_{cfg$profile}.tsv")), sep = "\t", quote = FALSE)
utils::write.table(meta, file.path(cfg$output_root, "global_local_cluster_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(data.frame(global_local_cluster_id = rownames(X), weight = weights), file.path(cfg$output_root, "global_local_cluster_weights.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

X_t <- transform_matrix(X, cfg$transform, cfg$eps)
cor_mat <- suppressWarnings(stats::cor(t(X_t), method = cfg$corr_method, use = "pairwise.complete.obs"))
cor_mat[!is.finite(cor_mat) | is.na(cor_mat)] <- 0
dist_mat <- stats::as.dist(1 - cor_mat)
hc <- stats::hclust(dist_mat, method = cfg$hclust_method)

k_rows <- list()
for (K in cfg$k_candidates) {
  labels <- stats::cutree(hc, k = K)
  cov <- coverage_score(labels, meta)
  stats_df <- subtype_group_stats(labels, meta, weights, cfg$cov_min, cfg$pur_min, cfg$lambda_p, cfg$lambda_l)
  group_pen <- mean(stats_df$penalty, na.rm = TRUE)
  within <- within_dispersion(X_t, labels, weights)
  k_rows[[as.character(K)]] <- data.frame(
    K = K,
    coverage = cov,
    group_penalty = group_pen,
    within_dispersion = within,
    stringsAsFactors = FALSE
  )
}
k_scan <- do.call(rbind, k_rows)
within_vals <- k_scan$within_dispersion
if (length(unique(within_vals[is.finite(within_vals)])) <= 1) {
  k_scan$within_scaled <- 0
} else {
  rng <- range(within_vals, na.rm = TRUE)
  k_scan$within_scaled <- (within_vals - rng[1]) / (rng[2] - rng[1] + 1e-9)
  k_scan$within_scaled[!is.finite(k_scan$within_scaled) | is.na(k_scan$within_scaled)] <- 1
}
k_scan$score <- cfg$alpha * k_scan$coverage - cfg$beta * k_scan$group_penalty - cfg$gamma * k_scan$within_scaled
k_best <- k_scan$K[which.max(k_scan$score)]

utils::write.table(k_scan, file.path(cfg$output_root, "global_k_scan.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message(glue("Selected K_global={k_best}"))

labels0 <- stats::cutree(hc, k = k_best)
labels1 <- merge_bad_subtypes(
  X_t,
  labels0,
  meta,
  weights,
  cov_min = cfg$cov_min,
  pur_min = cfg$pur_min,
  lambda_p = cfg$lambda_p,
  lambda_l = cfg$lambda_l,
  p_max = cfg$p_max_sample,
  low_weight_q = cfg$low_weight_q,
  pen_soft_max = cfg$pen_soft_max
)

labs_final <- levels(labels1)
sub_weight <- sapply(labs_final, function(g) sum(weights[labels1 == g]))
labs_final <- labs_final[order(sub_weight, decreasing = TRUE)]
sub_map <- stats::setNames(paste0("S", seq_along(labs_final)), labs_final)
global_subtype <- as.character(sub_map[as.character(labels1)])

cent_rows <- list()
for (g in sort(unique(global_subtype))) {
  idx <- which(global_subtype == g)
  cent_rows[[g]] <- weighted_centroid(X[idx, , drop = FALSE], weights[idx])
}
cent <- do.call(rbind, cent_rows)
rownames(cent) <- sort(unique(global_subtype))
colnames(cent) <- colnames(X)

mapping <- data.frame(
  global_local_cluster_id = rownames(X),
  sample_id = meta$sample_id,
  group = meta$group,
  local_cluster = meta$local_cluster,
  n_cells = meta$n_cells,
  weight = weights,
  global_subtype = global_subtype,
  stringsAsFactors = FALSE
)

utils::write.table(mapping, file.path(cfg$output_root, "global_subtype_mapping.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(cent, file.path(cfg$output_root, glue("global_subtype_centroids_{cfg$profile}.tsv")), sep = "\t", quote = FALSE)

stats_final <- subtype_group_stats(global_subtype, mapping, mapping$weight, cfg$cov_min, cfg$pur_min, cfg$lambda_p, cfg$lambda_l)
p_share <- max_sample_share(global_subtype, mapping, mapping$weight)
stats_final$p_max_sample <- as.numeric(p_share[stats_final$subtype])
utils::write.table(stats_final, file.path(cfg$output_root, "global_subtype_stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

if (cfg$write_heatmaps) {
  # Original heatmap style (earliest version):
  # - Row/col order comes purely from hierarchical clustering (subtypes may interleave)
  # - Show dendrograms
  # - Annotate by group + global_subtype (no sample_id track)
  ann <- data.frame(group = mapping$group, global_subtype = mapping$global_subtype)
  rownames(ann) <- mapping$global_local_cluster_id
  ann$group <- factor(ann$group, levels = c("P", "NonP"))
  subtype_levels <- sort(unique(mapping$global_subtype))
  subtype_key <- suppressWarnings(as.integer(sub("^S([0-9]+)$", "\\1", subtype_levels)))
  if (!all(is.na(subtype_key))) subtype_levels <- subtype_levels[order(subtype_key, subtype_levels)]
  ann$global_subtype <- factor(ann$global_subtype, levels = subtype_levels)

  group_cols <- c(P = "#4F81BD", NonP = "#C0504D")
  subtype_cols <- stats::setNames(ggsci::pal_d3("category20")(length(subtype_levels)), subtype_levels)
  ann_cols <- list(group = group_cols, global_subtype = subtype_cols)

  grDevices::cairo_pdf(file.path(cfg$output_root, "global_similarity_heatmap.pdf"), width = 10, height = 9)
  pheatmap::pheatmap(
    cor_mat,
    cluster_rows = hc,
    cluster_cols = hc,
    treeheight_row = 25,
    treeheight_col = 25,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_row = ann,
    annotation_col = ann,
    annotation_colors = ann_cols,
    color = grDevices::colorRampPalette(c("#1A1A1A", "#FFFFFF", "#B2182B"))(100),
    border_color = NA
  )
  grDevices::dev.off()

  # Optional: rotated leaf order to reduce subtype fragmentation while preserving the same dendrogram.
  ids_rot <- rotate_leaf_order_by_label(hc, leaf_labels = stats::setNames(mapping$global_subtype, mapping$global_local_cluster_id)[hc$labels])
  cor_rot <- cor_mat[ids_rot, ids_rot, drop = FALSE]
  ann_rot <- ann[ids_rot, , drop = FALSE]
  hc_rot <- hc
  hc_rot$order <- match(ids_rot, hc_rot$labels)

  runs_before <- count_label_runs(stats::setNames(mapping$global_subtype, mapping$global_local_cluster_id)[hc$labels[hc$order]])
  runs_after <- count_label_runs(stats::setNames(mapping$global_subtype, mapping$global_local_cluster_id)[ids_rot])
  runs_out <- merge(runs_before, runs_after, by = "label", all = TRUE, suffixes = c("_hc", "_rot"))
  runs_out$runs_hc[is.na(runs_out$runs_hc)] <- 0L
  runs_out$runs_rot[is.na(runs_out$runs_rot)] <- 0L
  utils::write.table(runs_out, file.path(cfg$output_root, "global_similarity_heatmap_run_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  grDevices::cairo_pdf(file.path(cfg$output_root, "global_similarity_heatmap_rotated_by_subtype.pdf"), width = 10, height = 9)
  pheatmap::pheatmap(
    cor_rot,
    cluster_rows = hc_rot,
    cluster_cols = hc_rot,
    treeheight_row = 25,
    treeheight_col = 25,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_row = ann_rot,
    annotation_col = ann_rot,
    annotation_colors = ann_cols,
    color = grDevices::colorRampPalette(c("#1A1A1A", "#FFFFFF", "#B2182B"))(100),
    border_color = NA
  )
  grDevices::dev.off()

  # subtype x sample weight heatmap
  subtypes <- sort(unique(mapping$global_subtype))
  samples <- sort(unique(mapping$sample_id))
  mat <- matrix(0, nrow = length(subtypes), ncol = length(samples))
  rownames(mat) <- subtypes
  colnames(mat) <- samples
  for (g in subtypes) {
    for (s in samples) {
      mat[g, s] <- sum(mapping$weight[mapping$global_subtype == g & mapping$sample_id == s])
    }
  }
  grDevices::cairo_pdf(file.path(cfg$output_root, "subtype_by_sample_weight_heatmap.pdf"), width = 10, height = 3 + 0.25 * length(subtypes))
  pheatmap::pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = NA,
    color = grDevices::colorRampPalette(c("#FFFFFF", "#2166AC"))(100)
  )
  grDevices::dev.off()
}

message(glue("Done. Output: {cfg$output_root}"))
