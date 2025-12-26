#!/usr/bin/env Rscript

# CAF niche global alignment across samples.
suppressPackageStartupMessages({
  library(ggplot2)
  library(glue)
  library(pheatmap)
  library(ggsci)
})

`%||%` <- function(x, y) if (is.null(x) || x == "") y else x

CANONICAL_COARSE <- character(0)

get_root_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) > 0) {
    return(normalizePath(file.path(dirname(file_arg[1]), "..")))
  }
  normalizePath(getwd())
}

root_dir <- get_root_dir()
source(file.path(root_dir, "utils", "R", "coarse_schema.R"))

infer_group <- function(sample_id) {
  if (grepl("P$", sample_id)) return("Primary (P)")
  if (grepl("[LR]$", sample_id)) return("Metastasis (L/R)")
  NA_character_
}

palette_for_levels <- function(levels, palette = "category20") {
  levels <- levels[!is.na(levels)]
  if (length(levels) == 0) return(c())
  pal <- switch(
    palette,
    category20 = ggsci::pal_d3("category20")(20),
    npg = ggsci::pal_npg("nrc")(10),
    lancet = ggsci::pal_lancet("lanonc")(9),
    nejm = ggsci::pal_nejm("default")(8),
    ggsci::pal_d3("category20")(20)
  )
  pal <- rep(pal, length.out = length(levels))
  stats::setNames(pal, levels)
}

palette_for_group <- function(levels, avoid_colors = character(0)) {
  palettes <- c("npg", "lancet", "nejm", "category20")
  for (pal_name in palettes) {
    pal <- palette_for_levels(levels, palette = pal_name)
    if (length(intersect(tolower(pal), tolower(avoid_colors))) == 0) {
      return(pal)
    }
  }
  palette_for_levels(levels, palette = "category20")
}

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

PARAM_HELP <- list(
  input_root = "Root directory of section-normalized CAF niche outputs",
  output_root = "Root directory to write global alignment outputs",
  samples = "Comma-separated sample IDs (blank = auto-detect)",
  k_local = "Local CAF clustering k",
  profile = "Profile matrix to use: corrected or raw",
  coarse_schema = "Coarse schema: base or myeloid_refined",
  corr_method = "Correlation method: pearson or spearman",
  transform = "Transform: log2p, log1p (log1p(M)), or none",
  eps = "Epsilon for log2p",
  hclust_method = "hclust method (e.g., average)",
  k_global = "Global subtype count (auto or integer)",
  k_global_min = "Min K for auto scan",
  k_global_max = "Max K for auto scan",
  min_cluster_cells = "Exclude local clusters with fewer cells",
  min_subtype_samples = "Minimum distinct samples per global subtype (merge rare subtypes)",
  sample_meta = "Optional TSV with sample_id and annotations",
  write_renamed_clusters = "Write caf_clusters_k{k}_global.tsv per sample"
)

print_help <- function(defaults) {
  cat("CAF niche global alignment - parameters\n")
  for (key in names(PARAM_HELP)) {
    def <- defaults[[key]]
    def_str <- if (is.na(def) || is.null(def)) "NA" else as.character(def)
    cat(glue("  --{key} {PARAM_HELP[[key]]} (default: {def_str})\n"))
  }
}

align_profile_columns <- function(mat, sample_id, path_label) {
  if (is.null(colnames(mat))) {
    stop(glue("{sample_id}: {path_label} missing column names."))
  }
  extras <- setdiff(colnames(mat), CANONICAL_COARSE)
  if (length(extras) > 0) {
    stop(glue("{sample_id}: {path_label} has unexpected columns: {paste(extras, collapse = ', ')}"))
  }
  missing <- setdiff(CANONICAL_COARSE, colnames(mat))
  if (length(missing) > 0) {
    zeros <- matrix(0, nrow = nrow(mat), ncol = length(missing))
    colnames(zeros) <- missing
    mat <- cbind(mat, zeros)
  }
  mat[, CANONICAL_COARSE, drop = FALSE]
}

scan_global_k <- function(cor_mat, hc, k_min, k_max) {
  scan_rows <- list()
  pair_idx <- upper.tri(cor_mat)
  pair_vals <- cor_mat[pair_idx]
  row_idx <- row(cor_mat)[pair_idx]
  col_idx <- col(cor_mat)[pair_idx]

  for (k in k_min:k_max) {
    labels <- cutree(hc, k = k)
    same <- labels[row_idx] == labels[col_idx]
    within_vals <- pair_vals[same]
    between_vals <- pair_vals[!same]
    within_med <- if (length(within_vals) > 0) median(within_vals) else NA_real_
    between_med <- if (length(between_vals) > 0) median(between_vals) else NA_real_
    score <- if (is.na(within_med) || is.na(between_med)) NA_real_ else within_med - between_med
    scan_rows[[length(scan_rows) + 1]] <- data.frame(
      k = k,
      within_median = within_med,
      between_median = between_med,
      score = score
    )
  }
  do.call(rbind, scan_rows)
}

choose_k_from_scan <- function(scan_df, k_min) {
  valid <- scan_df[is.finite(scan_df$score), , drop = FALSE]
  if (nrow(valid) == 0) {
    return(list(k = k_min, reason = "fallback_min"))
  }
  best <- valid[valid$score == max(valid$score), , drop = FALSE]
  best_k <- min(best$k)
  list(k = best_k, reason = "max_score")
}

count_samples_per_label <- function(labels, sample_ids) {
  labs <- sort(unique(labels))
  counts <- sapply(labs, function(k) length(unique(sample_ids[labels == k])))
  stats::setNames(counts, labs)
}

reduce_k_for_min_samples <- function(hc, sample_ids, k_start, min_samples, k_min = 2L) {
  if (min_samples <= 1) {
    return(list(k = k_start, labels = stats::cutree(hc, k = k_start), coverage = NULL))
  }
  if (length(sample_ids) == 0) {
    stop("sample_ids is empty; cannot enforce min_samples.")
  }
  k_start <- as.integer(k_start)
  k_min <- as.integer(k_min)
  if (is.na(k_start) || k_start < k_min) {
    stop("k_start must be >= k_min.")
  }

  best <- list(k = k_start, labels = stats::cutree(hc, k = k_start), coverage = NULL)
  for (k in seq(k_start, k_min, by = -1L)) {
    labels <- stats::cutree(hc, k = k)
    coverage <- count_samples_per_label(labels, sample_ids)
    if (all(coverage >= min_samples)) {
      best <- list(k = k, labels = labels, coverage = coverage)
      break
    }
    best <- list(k = k, labels = labels, coverage = coverage)
  }
  best
}

argv <- commandArgs(trailingOnly = TRUE)
opts <- parse_args(argv)

defaults <- list(
  input_root = file.path(getwd(), "data", "caf_niche_sectionnorm"),
  output_root = file.path(getwd(), "data", "caf_niche_sectionnorm", "_global"),
  samples = "",
  k_local = 5L,
  profile = "corrected",
  coarse_schema = "myeloid_refined",
  corr_method = "pearson",
  transform = "log2p",
  eps = 1e-6,
  hclust_method = "average",
  k_global = "auto",
  k_global_min = 2L,
  k_global_max = 8L,
  min_cluster_cells = 200L,
  min_subtype_samples = 1L,
  sample_meta = "",
  write_renamed_clusters = TRUE
)

if (!is.null(opts$help) || !is.null(opts$h)) {
  print_help(defaults)
  quit(status = 0)
}

cfg <- list(
  input_root = get_opt(opts, "input_root", defaults$input_root),
  output_root = get_opt(opts, "output_root", defaults$output_root),
  samples = get_opt(opts, "samples", defaults$samples),
  k_local = get_opt(opts, "k_local", defaults$k_local, as.integer),
  profile = get_opt(opts, "profile", defaults$profile),
  coarse_schema = get_opt(opts, "coarse_schema", defaults$coarse_schema),
  corr_method = get_opt(opts, "corr_method", defaults$corr_method),
  transform = get_opt(opts, "transform", defaults$transform),
  eps = get_opt(opts, "eps", defaults$eps, as.numeric),
  hclust_method = get_opt(opts, "hclust_method", defaults$hclust_method),
  k_global = get_opt(opts, "k_global", defaults$k_global),
  k_global_min = get_opt(opts, "k_global_min", defaults$k_global_min, as.integer),
  k_global_max = get_opt(opts, "k_global_max", defaults$k_global_max, as.integer),
  min_cluster_cells = get_opt(opts, "min_cluster_cells", defaults$min_cluster_cells, as.integer),
  min_subtype_samples = get_opt(opts, "min_subtype_samples", defaults$min_subtype_samples, as.integer),
  sample_meta = get_opt(opts, "sample_meta", defaults$sample_meta),
  write_renamed_clusters = get_opt(opts, "write_renamed_clusters", defaults$write_renamed_clusters, to_bool)
)

cfg$profile <- tolower(cfg$profile)
if (!cfg$profile %in% c("raw", "corrected")) {
  stop("profile must be 'raw' or 'corrected'.")
}

cfg$coarse_schema <- tolower(cfg$coarse_schema)
if (!cfg$coarse_schema %in% c("base", "myeloid_refined")) {
  stop(glue("Unknown coarse_schema: {cfg$coarse_schema}. Use base or myeloid_refined."))
}
CANONICAL_COARSE <- get_canonical_coarse(cfg$coarse_schema, order = "analysis")

cfg$corr_method <- tolower(cfg$corr_method)
if (!cfg$corr_method %in% c("pearson", "spearman")) {
  stop("corr_method must be 'pearson' or 'spearman'.")
}

cfg$transform <- tolower(cfg$transform)
if (!cfg$transform %in% c("log2p", "log1p", "none")) {
  stop("transform must be 'log2p', 'log1p', or 'none'.")
}

if (cfg$profile == "corrected" && cfg$transform != "none") {
  warning("profile=corrected now uses log2 enrichment; reset transform to 'none'.")
  cfg$transform <- "none"
}

if (is.na(cfg$min_subtype_samples) || cfg$min_subtype_samples < 1) {
  stop("min_subtype_samples must be >= 1.")
}

sample_ids <- character(0)
if (cfg$samples != "") {
  sample_ids <- trimws(strsplit(cfg$samples, ",")[[1]])
} else if (dir.exists(cfg$input_root)) {
  sample_ids <- list.dirs(cfg$input_root, recursive = FALSE, full.names = FALSE)
  sample_ids <- sample_ids[sample_ids != ""]
  sample_ids <- sample_ids[!startsWith(sample_ids, "_")]
}
if (length(sample_ids) == 0) {
  stop("No samples found.")
}

sample_meta_df <- NULL
if (cfg$sample_meta != "") {
  if (!file.exists(cfg$sample_meta)) {
    stop(glue("sample_meta not found: {cfg$sample_meta}"))
  }
  sample_meta_df <- read.delim(cfg$sample_meta, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"sample_id" %in% colnames(sample_meta_df)) {
    stop("sample_meta must contain a 'sample_id' column.")
  }
}

profile_rows <- list()
mapping_rows <- list()

for (sample_id in sample_ids) {
  sample_dir <- file.path(cfg$input_root, sample_id)
  if (!dir.exists(sample_dir)) {
    warning(glue("{sample_id}: sample dir missing; skip."))
    next
  }

  profile_path <- file.path(sample_dir, glue("niche_average_{cfg$profile}_k{cfg$k_local}.tsv"))
  summary_path <- file.path(sample_dir, glue("caf_cluster_summary_k{cfg$k_local}.tsv"))

  if (!file.exists(profile_path)) {
    warning(glue("{sample_id}: {basename(profile_path)} missing; skip."))
    next
  }
  if (!file.exists(summary_path)) {
    warning(glue("{sample_id}: {basename(summary_path)} missing; skip."))
    next
  }

  profile_mat <- read.delim(profile_path, row.names = 1, check.names = FALSE)
  profile_mat <- align_profile_columns(profile_mat, sample_id, basename(profile_path))
  summary_df <- read.delim(summary_path, stringsAsFactors = FALSE)

  local_clusters <- rownames(profile_mat)
  n_cells <- summary_df$n_cells[match(local_clusters, summary_df$cluster)]
  if (any(is.na(n_cells))) {
    missing_clusters <- local_clusters[is.na(n_cells)]
    warning(glue("{sample_id}: missing n_cells for {paste(missing_clusters, collapse = ', ')}"))
  }

  cluster_id <- paste0(sample_id, "::", local_clusters)
  excluded <- is.na(n_cells) | n_cells < cfg$min_cluster_cells

  mapping_rows[[length(mapping_rows) + 1]] <- data.frame(
    cluster_id = cluster_id,
    sample_id = sample_id,
    local_cluster = local_clusters,
    n_cells = n_cells,
    excluded_from_global = excluded,
    global_subtype = NA_integer_,
    global_label = NA_character_,
    stringsAsFactors = FALSE
  )

  if (any(!excluded)) {
    profile_keep <- profile_mat[!excluded, , drop = FALSE]
    rownames(profile_keep) <- cluster_id[!excluded]
    profile_rows[[length(profile_rows) + 1]] <- profile_keep
  }
}

mapping_df <- do.call(rbind, mapping_rows)
if (is.null(mapping_df) || nrow(mapping_df) == 0) {
  stop("No valid clusters found for global alignment.")
}

if (length(profile_rows) == 0) {
  stop("No clusters meet min_cluster_cells; nothing to align.")
}

profiles_mat <- do.call(rbind, profile_rows)
if (nrow(profiles_mat) < 2) {
  stop("Need at least 2 clusters to compute global alignment.")
}

if (cfg$transform == "log2p") {
  profiles_mat <- log2(profiles_mat + cfg$eps)
} else if (cfg$transform == "log1p") {
  profiles_mat <- log1p(profiles_mat)
}

cor_mat <- stats::cor(t(profiles_mat), method = cfg$corr_method, use = "pairwise.complete.obs")
cor_mat[!is.finite(cor_mat)] <- 0
diag(cor_mat) <- 1
cor_mat <- pmax(pmin(cor_mat, 1), -1)

dist_mat <- stats::as.dist(1 - cor_mat)
hc <- stats::hclust(dist_mat, method = cfg$hclust_method)

k_global <- NA_integer_
scan_df <- NULL
if (tolower(cfg$k_global) == "auto") {
  k_min <- max(cfg$k_global_min, 2L)
  k_max <- min(cfg$k_global_max, nrow(profiles_mat))
  if (k_max < k_min) k_min <- k_max
  scan_df <- scan_global_k(cor_mat, hc, k_min, k_max)
  scan_choice <- choose_k_from_scan(scan_df, k_min)
  k_global <- scan_choice$k
} else {
  k_global <- as.integer(cfg$k_global)
  if (is.na(k_global) || k_global < 2) {
    stop("k_global must be 'auto' or an integer >= 2.")
  }
  if (k_global > nrow(profiles_mat)) {
    warning(glue("k_global {k_global} > clusters {nrow(profiles_mat)}; reset to {nrow(profiles_mat)}"))
    k_global <- nrow(profiles_mat)
  }
}

if (tolower(cfg$k_global) == "auto" && !is.null(scan_df)) {
  out_dir <- file.path(cfg$output_root, glue("k{cfg$k_local}"), cfg$profile)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    scan_df,
    file.path(out_dir, "global_k_scan.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  p <- ggplot(scan_df, aes(k, score)) +
    geom_line(linewidth = 0.8, color = "#2C3E50") +
    geom_point(size = 1.8, color = "#2C3E50") +
    theme_minimal(base_size = 11) +
    labs(x = "K", y = "within - between median corr", title = "Global K scan")
  ggsave(file.path(out_dir, "global_k_scan.pdf"), p, width = 5.5, height = 4, device = grDevices::cairo_pdf)
}

cluster_samples <- mapping_df$sample_id[match(rownames(profiles_mat), mapping_df$cluster_id)]
if (any(is.na(cluster_samples))) {
  stop("Missing sample_id for some clusters; cannot enforce min_subtype_samples.")
}

reduced <- reduce_k_for_min_samples(
  hc = hc,
  sample_ids = cluster_samples,
  k_start = k_global,
  min_samples = cfg$min_subtype_samples,
  k_min = 2L
)
labels <- reduced$labels
if (reduced$k != k_global) {
  message(glue("Reduced K to satisfy min_subtype_samples={cfg$min_subtype_samples}: K {k_global} -> {reduced$k}"))
}
k_global <- reduced$k
ordered_ids <- rownames(profiles_mat)[hc$order]
ordered_labels <- labels[ordered_ids]
label_order <- ordered_labels[!duplicated(ordered_labels)]
label_map <- setNames(seq_along(label_order), label_order)
global_subtype <- label_map[as.character(labels)]
global_label <- paste0("s", global_subtype, "-CAFs")
k_global <- length(label_order)

mapping_idx <- match(rownames(profiles_mat), mapping_df$cluster_id)
mapping_df$global_subtype[mapping_idx] <- global_subtype
mapping_df$global_label[mapping_idx] <- global_label

centroids <- sapply(seq_len(k_global), function(k) {
  colMeans(profiles_mat[global_subtype == k, , drop = FALSE])
})
centroids <- t(centroids)
rownames(centroids) <- paste0("s", seq_len(k_global), "-CAFs")
colnames(centroids) <- colnames(profiles_mat)

out_dir <- file.path(cfg$output_root, glue("k{cfg$k_local}"), cfg$profile)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

utils::write.table(
  profiles_mat,
  file.path(out_dir, "profiles.tsv"),
  sep = "\t",
  quote = FALSE
)
utils::write.table(
  cor_mat,
  file.path(out_dir, "correlation.tsv"),
  sep = "\t",
  quote = FALSE
)
utils::write.table(
  data.frame(order = seq_along(ordered_ids), cluster_id = ordered_ids),
  file.path(out_dir, "hclust_order.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  mapping_df,
  file.path(out_dir, "global_subtype_mapping.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  centroids,
  file.path(out_dir, "global_subtype_centroids.tsv"),
  sep = "\t",
  quote = FALSE
)

subtype_list <- sort(unique(mapping_df$global_label[!is.na(mapping_df$global_label)]))
sum_raw <- matrix(0, nrow = length(subtype_list), ncol = length(CANONICAL_COARSE))
sum_corr <- matrix(0, nrow = length(subtype_list), ncol = length(CANONICAL_COARSE))
counts <- setNames(rep(0, length(subtype_list)), subtype_list)
rownames(sum_raw) <- subtype_list
colnames(sum_raw) <- CANONICAL_COARSE
rownames(sum_corr) <- subtype_list
colnames(sum_corr) <- CANONICAL_COARSE

for (sample_id in unique(mapping_df$sample_id)) {
  sample_dir <- file.path(cfg$input_root, sample_id)
  raw_path <- file.path(sample_dir, "niche_celltype_fraction.tsv")
  corr_path <- file.path(sample_dir, "niche_celltype_corrected.tsv")
  cluster_path <- file.path(sample_dir, glue("caf_clusters_k{cfg$k_local}.tsv"))
  if (!file.exists(raw_path) || !file.exists(corr_path) || !file.exists(cluster_path)) {
    warning(glue("{sample_id}: missing niche or cluster files; skip composition."))
    next
  }

  raw_mat <- read.delim(raw_path, row.names = 1, check.names = FALSE)
  corr_mat <- read.delim(corr_path, row.names = 1, check.names = FALSE)
  raw_mat <- align_profile_columns(raw_mat, sample_id, basename(raw_path))
  corr_mat <- align_profile_columns(corr_mat, sample_id, basename(corr_path))

  cluster_df <- read.delim(cluster_path, stringsAsFactors = FALSE)
  map_sub <- mapping_df[mapping_df$sample_id == sample_id, c("local_cluster", "global_label")]
  map_vec <- setNames(map_sub$global_label, map_sub$local_cluster)
  cluster_df$global_label <- map_vec[cluster_df$caf_cluster]
  cluster_df <- cluster_df[!is.na(cluster_df$global_label), , drop = FALSE]
  cluster_df <- cluster_df[cluster_df$cell_id %in% rownames(raw_mat), , drop = FALSE]
  if (nrow(cluster_df) == 0) next

  split_cells <- split(cluster_df$cell_id, cluster_df$global_label)
  for (subtype in names(split_cells)) {
    cells <- split_cells[[subtype]]
    if (length(cells) == 0) next
    sum_raw[subtype, ] <- sum_raw[subtype, ] + colSums(raw_mat[cells, , drop = FALSE])
    sum_corr[subtype, ] <- sum_corr[subtype, ] + colSums(corr_mat[cells, , drop = FALSE])
    counts[subtype] <- counts[subtype] + length(cells)
  }
}

raw_mean <- sum_raw
corr_mean <- sum_corr
for (subtype in subtype_list) {
  if (counts[subtype] > 0) {
    raw_mean[subtype, ] <- sum_raw[subtype, ] / counts[subtype]
    corr_mean[subtype, ] <- sum_corr[subtype, ] / counts[subtype]
  } else {
    raw_mean[subtype, ] <- NA_real_
    corr_mean[subtype, ] <- NA_real_
  }
}

raw_frac <- raw_mean
row_sums <- rowSums(raw_frac, na.rm = TRUE)
raw_frac[row_sums > 0, ] <- raw_frac[row_sums > 0, , drop = FALSE] / row_sums[row_sums > 0]

pie_df <- data.frame(
  subtype = rep(rownames(raw_frac), times = ncol(raw_frac)),
  cell_type = rep(colnames(raw_frac), each = nrow(raw_frac)),
  fraction = as.vector(raw_frac),
  enrichment = as.vector(corr_mean),
  stringsAsFactors = FALSE
)
pie_df <- pie_df[!is.na(pie_df$fraction), , drop = FALSE]
pie_df$subtype <- factor(pie_df$subtype, levels = rownames(raw_frac))
pie_df$cell_type <- factor(pie_df$cell_type, levels = colnames(raw_frac))

label_min_fraction <- 0.08
label_top_n <- 2
label_df <- pie_df[order(pie_df$subtype, pie_df$cell_type), ]
label_df$ypos <- ave(label_df$fraction, label_df$subtype, FUN = function(x) cumsum(x) - x / 2)
label_df$label <- sprintf(
  "%.1f%% %s",
  100 * label_df$fraction,
  gsub("_", " ", label_df$cell_type)
)
label_list <- lapply(split(label_df, label_df$subtype), function(df) {
  df_enriched <- df[!is.na(df$enrichment) & df$enrichment >= 0, , drop = FALSE]
  df_enriched <- df_enriched[order(df_enriched$enrichment, df_enriched$fraction, decreasing = TRUE), , drop = FALSE]
  df_enriched <- df_enriched[df_enriched$fraction >= label_min_fraction, , drop = FALSE]
  if (nrow(df_enriched) == 0) {
    df_enriched <- df[order(df$fraction, decreasing = TRUE), , drop = FALSE]
  }
  head(df_enriched, label_top_n)
})
label_df <- do.call(rbind, label_list)

pie_df$x <- 1
pie_plot <- ggplot(pie_df, aes(x = x, y = fraction, fill = cell_type)) +
  geom_col(width = 1, color = "white", linewidth = 0.2) +
  coord_polar(theta = "y") +
  facet_wrap(~subtype) +
  theme_void(base_size = 12) +
  theme(legend.position = "right") +
  labs(fill = "Cell type", title = "Neighboring cell composition") +
  scale_x_continuous(limits = c(0.5, 1.6))

pie_plot <- pie_plot + ggsci::scale_fill_d3(
  "category20",
  limits = levels(pie_df$cell_type),
  breaks = levels(pie_df$cell_type),
  drop = FALSE
)

if (!is.null(label_df) && nrow(label_df) > 0) {
  pie_plot <- pie_plot + geom_text(
    data = label_df,
    aes(x = 1.2, y = ypos, label = label),
    inherit.aes = FALSE,
    size = 3,
    hjust = 0,
    color = "black",
    check_overlap = TRUE
  )
}

ggsave(
  file.path(out_dir, "global_subtype_feature_pie.pdf"),
  pie_plot,
  width = 7,
  height = 4.5,
  device = grDevices::cairo_pdf
)

group_df <- mapping_df[mapping_df$excluded_from_global == FALSE, , drop = FALSE]
group_df$Group <- vapply(group_df$sample_id, infer_group, character(1))
group_df <- group_df[!is.na(group_df$Group), , drop = FALSE]
if (nrow(group_df) > 0) {
  group_df$subtype <- group_df$global_label
  group_df$Group <- factor(group_df$Group, levels = c("Primary (P)", "Metastasis (L/R)"))
  group_df$subtype <- factor(group_df$subtype, levels = paste0("s", seq_len(k_global), "-CAFs"))

  group_agg <- aggregate(n_cells ~ Group + subtype, data = group_df, sum)
  group_totals <- aggregate(n_cells ~ Group, data = group_agg, sum)
  group_agg <- merge(group_agg, group_totals, by = "Group", suffixes = c("", "_total"))
  group_agg$fraction <- group_agg$n_cells / group_agg$n_cells_total

  p_group_frac <- ggplot(group_agg, aes(Group, fraction, fill = subtype)) +
    geom_col(width = 0.7) +
    theme_minimal(base_size = 12) +
    labs(x = NULL, y = "Fraction", fill = "Subtype", title = "CAF subtype composition by group")
  p_group_frac <- p_group_frac + ggsci::scale_fill_d3(
    "category20",
    limits = levels(group_agg$subtype),
    breaks = levels(group_agg$subtype),
    drop = FALSE
  )
  ggsave(
    file.path(out_dir, "group_subtype_fraction.pdf"),
    p_group_frac,
    width = 6,
    height = 4,
    device = grDevices::cairo_pdf
  )

  p_group_counts <- ggplot(group_agg, aes(Group, n_cells, fill = subtype)) +
    geom_col(width = 0.7) +
    theme_minimal(base_size = 12) +
    labs(x = NULL, y = "CAF cells", fill = "Subtype", title = "CAF subtype counts by group")
  p_group_counts <- p_group_counts + ggsci::scale_fill_d3(
    "category20",
    limits = levels(group_agg$subtype),
    breaks = levels(group_agg$subtype),
    drop = FALSE
  )
  ggsave(
    file.path(out_dir, "group_subtype_counts.pdf"),
    p_group_counts,
    width = 6,
    height = 4,
    device = grDevices::cairo_pdf
  )

  sample_agg <- aggregate(n_cells ~ Group + sample_id + subtype, data = group_df, sum)
  sample_totals <- aggregate(n_cells ~ Group + sample_id, data = sample_agg, sum)
  sample_agg <- merge(sample_agg, sample_totals, by = c("Group", "sample_id"), suffixes = c("", "_total"))
  sample_agg$fraction <- sample_agg$n_cells / sample_agg$n_cells_total

  sample_order <- unique(sample_agg$sample_id[order(sample_agg$Group, sample_agg$sample_id)])
  sample_agg$sample_id <- factor(sample_agg$sample_id, levels = sample_order)
  sample_agg$subtype <- factor(sample_agg$subtype, levels = paste0("s", seq_len(k_global), "-CAFs"))

  p_sample_frac <- ggplot(sample_agg, aes(sample_id, fraction, fill = subtype)) +
    geom_col(width = 0.85) +
    facet_wrap(~Group, scales = "free_x") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = "Fraction", fill = "Subtype", title = "CAF subtype fractions by sample")
  p_sample_frac <- p_sample_frac + ggsci::scale_fill_d3(
    "category20",
    limits = levels(sample_agg$subtype),
    breaks = levels(sample_agg$subtype),
    drop = FALSE
  )
  ggsave(
    file.path(out_dir, "group_subtype_fraction_by_sample.pdf"),
    p_sample_frac,
    width = 10,
    height = 4.5,
    device = grDevices::cairo_pdf
  )
}

annotation_df <- mapping_df[mapping_df$excluded_from_global == FALSE, , drop = FALSE]
annotation_df <- annotation_df[match(rownames(profiles_mat), annotation_df$cluster_id), , drop = FALSE]
annotation_row <- data.frame(
  subtype = annotation_df$global_label,
  Group = vapply(annotation_df$sample_id, infer_group, character(1)),
  stringsAsFactors = FALSE
)
if (!is.null(sample_meta_df)) {
  meta_match <- sample_meta_df[match(annotation_df$sample_id, sample_meta_df$sample_id), , drop = FALSE]
  extra_cols <- setdiff(colnames(meta_match), "sample_id")
  for (col in extra_cols) {
    annotation_row[[col]] <- meta_match[[col]]
  }
}
rownames(annotation_row) <- rownames(profiles_mat)
annotation_row$subtype <- factor(annotation_row$subtype, levels = paste0("s", seq_len(k_global), "-CAFs"))
annotation_row$Group <- factor(annotation_row$Group, levels = c("Primary (P)", "Metastasis (L/R)"))

annotation_colors <- list()
annotation_colors$subtype <- palette_for_levels(levels(annotation_row$subtype), palette = "category20")
annotation_colors$Group <- palette_for_group(
  levels(annotation_row$Group),
  avoid_colors = annotation_colors$subtype
)
if (!is.null(sample_meta_df)) {
  extra_cols <- setdiff(colnames(annotation_row), c("subtype", "Group"))
  for (col in extra_cols) {
    vals <- unique(annotation_row[[col]])
    annotation_colors[[col]] <- palette_for_levels(vals)
  }
}

heatmap_size <- max(7, min(20, nrow(cor_mat) * 0.35))
heatmap_path <- file.path(out_dir, "global_similarity_heatmap.pdf")
grDevices::cairo_pdf(heatmap_path, width = heatmap_size, height = heatmap_size)
pheatmap::pheatmap(
  cor_mat,
  cluster_rows = hc,
  cluster_cols = hc,
  annotation_row = annotation_row,
  annotation_col = annotation_row,
  annotation_colors = annotation_colors,
  border_color = NA,
  main = glue("CAF global alignment (k={k_global}, {cfg$profile})")
)
grDevices::dev.off()

if (cfg$write_renamed_clusters) {
  for (sample_id in sample_ids) {
    sample_dir <- file.path(cfg$input_root, sample_id)
    cluster_path <- file.path(sample_dir, glue("caf_clusters_k{cfg$k_local}.tsv"))
    if (!file.exists(cluster_path)) {
      warning(glue("{sample_id}: {basename(cluster_path)} missing; skip renamed clusters."))
      next
    }
    cluster_df <- read.delim(cluster_path, stringsAsFactors = FALSE)
    map_df <- mapping_df[mapping_df$sample_id == sample_id, c("local_cluster", "global_label")]
    map_vec <- setNames(map_df$global_label, map_df$local_cluster)
    cluster_df$caf_cluster_global <- map_vec[cluster_df$caf_cluster]
    out_path <- file.path(sample_dir, glue("caf_clusters_k{cfg$k_local}_global.tsv"))
    utils::write.table(
      cluster_df,
      out_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
}
