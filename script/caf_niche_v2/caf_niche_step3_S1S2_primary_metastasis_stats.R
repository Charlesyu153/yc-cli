#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(glue)
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

get_opt <- function(opts, key, default, cast = NULL) {
  val <- opts[[key]]
  if (is.null(val) || val == "") return(default)
  if (is.null(cast)) return(val)
  cast(val)
}

parse_csv_list <- function(x) {
  x <- x %||% ""
  x <- trimws(x)
  if (x == "") return(character(0))
  out <- trimws(strsplit(x, ",")[[1]])
  out[out != ""]
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
  if (grepl("QT", sample_id, ignore.case = TRUE)) return(NA_character_)
  if (grepl("P$", sample_id, ignore.case = TRUE)) return("Primary")
  if (grepl("[LR]$", sample_id, ignore.case = TRUE)) return("Metastasis")
  NA_character_
}

read_cell_level_counts <- function(path, sample_id, subtypes) {
  df <- read.delim(path, stringsAsFactors = FALSE)
  if (!all(c("cell_id", "global_subtype") %in% colnames(df))) {
    stop(glue("{basename(path)} missing required columns cell_id/global_subtype"))
  }
  df <- df[!is.na(df$global_subtype) & df$global_subtype != "", , drop = FALSE]
  total <- nrow(df)
  counts <- stats::setNames(rep(0L, length(subtypes)), subtypes)
  if (total > 0) {
    tab <- table(df$global_subtype)
    hit <- intersect(names(tab), subtypes)
    counts[hit] <- as.integer(tab[hit])
  }
  out <- data.frame(
    sample_id = sample_id,
    group = infer_group(sample_id),
    subtype = subtypes,
    n = as.integer(counts),
    n_total = as.integer(total),
    stringsAsFactors = FALSE
  )
  out$prop <- ifelse(out$n_total > 0, out$n / out$n_total, NA_real_)
  out
}

write_pdf <- function(path, plot, width = 12, height = 4.5) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
root_dir <- get_root_dir()

defaults <- list(
  step2_cell_level_root = file.path(root_dir, "data", "caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2", "step2", "cell_level"),
  output_dir = file.path(root_dir, "data", "caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2", "step3", "S1S2_group_stats"),
  subtypes = "S1,S2",
  samples = "",
  font_family = "Times",
  font_base = 9
)

cfg <- list(
  step2_cell_level_root = get_opt(opts, "step2_cell_level_root", defaults$step2_cell_level_root),
  output_dir = get_opt(opts, "output_dir", defaults$output_dir),
  subtypes = parse_csv_list(get_opt(opts, "subtypes", defaults$subtypes)),
  samples = parse_csv_list(get_opt(opts, "samples", defaults$samples)),
  font_family = get_opt(opts, "font_family", defaults$font_family),
  font_base = get_opt(opts, "font_base", defaults$font_base, as.numeric)
)

if (length(cfg$subtypes) == 0) stop("No subtypes specified.")

files <- list.files(cfg$step2_cell_level_root, pattern = "_caf_cells_global_subtype\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop(glue("No cell-level files under {cfg$step2_cell_level_root}"))
sample_ids <- sub("_caf_cells_global_subtype\\.tsv$", "", basename(files))

if (length(cfg$samples) > 0) {
  keep <- sample_ids %in% cfg$samples
  files <- files[keep]
  sample_ids <- sample_ids[keep]
}

rows <- list()
for (i in seq_along(files)) {
  s <- sample_ids[i]
  g <- infer_group(s)
  if (is.na(g)) next
  rows[[length(rows) + 1]] <- read_cell_level_counts(files[i], s, cfg$subtypes)
}
stats_df <- do.call(rbind, rows)
if (is.null(stats_df) || nrow(stats_df) == 0) stop("No samples after filtering (QT excluded / no P/L/R samples found).")

stats_df$group <- factor(stats_df$group, levels = c("Primary", "Metastasis"))
stats_df$subtype <- factor(stats_df$subtype, levels = cfg$subtypes)

sample_order <- with(stats_df, unique(sample_id[order(group, sample_id)]))
stats_df$sample_id <- factor(stats_df$sample_id, levels = sample_order)

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
out_tsv <- file.path(cfg$output_dir, "S1S2_counts_props_by_sample.tsv")
utils::write.table(stats_df, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

subtype_colors <- stats::setNames(c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728")[seq_along(cfg$subtypes)], cfg$subtypes)

base_theme <- theme_minimal(base_size = cfg$font_base, base_family = cfg$font_family) +
  theme(
    text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "top"
  )

p_prop <- ggplot(stats_df, aes(x = sample_id, y = prop, fill = subtype)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, color = NA) +
  facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = subtype_colors, drop = FALSE) +
  scale_y_continuous(labels = function(x) sprintf("%.0f%%", x * 100), limits = c(0, 1)) +
  labs(x = NULL, y = "Fraction of CAF", title = "S1/S2 composition by sample") +
  base_theme

p_count <- ggplot(stats_df, aes(x = sample_id, y = n, fill = subtype)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, color = NA) +
  facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = subtype_colors, drop = FALSE) +
  labs(x = NULL, y = "CAF cell count", title = "S1/S2 counts by sample") +
  base_theme

out_prop <- file.path(cfg$output_dir, "S1S2_fraction_by_sample_primary_vs_metastasis.pdf")
out_count <- file.path(cfg$output_dir, "S1S2_count_by_sample_primary_vs_metastasis.pdf")

write_pdf(out_prop, p_prop, width = max(10, length(sample_order) * 0.35), height = 4.5)
write_pdf(out_count, p_count, width = max(10, length(sample_order) * 0.35), height = 4.5)

group_df <- aggregate(n ~ group + subtype, stats_df, sum)
group_total <- aggregate(n_total ~ group, unique(stats_df[, c("sample_id", "group", "n_total")]), sum)
group_df <- merge(group_df, group_total, by = "group", all.x = TRUE)
group_df$prop <- ifelse(group_df$n_total > 0, group_df$n / group_df$n_total, NA_real_)
group_df$group <- factor(group_df$group, levels = c("Primary", "Metastasis"))
group_df$subtype <- factor(group_df$subtype, levels = cfg$subtypes)

group_labels <- c(Primary = "Primary (P)", Metastasis = "Metastasis (L/R)")

p_group_prop <- ggplot(group_df, aes(x = group, y = prop, fill = subtype)) +
  geom_col(width = 0.65, color = NA) +
  scale_fill_manual(values = subtype_colors, drop = FALSE) +
  scale_x_discrete(labels = group_labels) +
  scale_y_continuous(labels = function(x) sprintf("%.0f%%", x * 100), limits = c(0, 1)) +
  labs(x = NULL, y = "Fraction", title = "CAF subtype composition by group") +
  base_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

out_group_prop <- file.path(cfg$output_dir, "S1S2_fraction_by_group_primary_vs_metastasis.pdf")
write_pdf(out_group_prop, p_group_prop, width = 4.0, height = 3.0)

if (!requireNamespace("ggalluvial", quietly = TRUE)) {
  warning("Missing package 'ggalluvial'; skip Sankey plots.")
} else {
  # Sankey 1: weights by CAF cell counts
  sankey_cells <- group_df
  sankey_cells$axis1 <- sankey_cells$group
  sankey_cells$axis2 <- sankey_cells$subtype
  sankey_cells$weight <- sankey_cells$n

  p_sankey_cells <- ggplot(
    sankey_cells,
    aes(axis1 = axis1, axis2 = axis2, y = weight)
  ) +
    ggalluvial::geom_alluvium(aes(fill = axis2), width = 1 / 12, alpha = 0.9) +
    ggalluvial::geom_stratum(width = 1 / 12, color = "grey30", fill = "white") +
    ggalluvial::stat_stratum(geom = "text", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Group", "Subtype"), expand = c(0.15, 0.05)) +
    scale_fill_manual(values = subtype_colors, drop = FALSE) +
    labs(y = "CAF cell count", x = NULL, title = "S1/S2 Sankey (weights = CAF cell count)") +
    theme_minimal(base_size = cfg$font_base, base_family = cfg$font_family) +
    theme(
      text = element_text(face = "bold"),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "top"
    )

  out_sankey_cells <- file.path(cfg$output_dir, "S1S2_sankey_by_group_cellcount.pdf")
  write_pdf(out_sankey_cells, p_sankey_cells, width = 5.5, height = 3.2)

  # Sankey 2: weights by number of samples with subtype present
  sample_hit <- aggregate(n ~ group + subtype, transform(stats_df, n = as.integer(n > 0)), sum)
  sample_hit$axis1 <- factor(sample_hit$group, levels = c("Primary", "Metastasis"))
  sample_hit$axis2 <- factor(sample_hit$subtype, levels = cfg$subtypes)
  sample_hit$weight <- sample_hit$n

  p_sankey_samples <- ggplot(
    sample_hit,
    aes(axis1 = axis1, axis2 = axis2, y = weight)
  ) +
    ggalluvial::geom_alluvium(aes(fill = axis2), width = 1 / 12, alpha = 0.9) +
    ggalluvial::geom_stratum(width = 1 / 12, color = "grey30", fill = "white") +
    ggalluvial::stat_stratum(geom = "text", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Group", "Subtype"), expand = c(0.15, 0.05)) +
    scale_fill_manual(values = subtype_colors, drop = FALSE) +
    labs(y = "Number of samples", x = NULL, title = "S1/S2 Sankey (weights = sample count)") +
    theme_minimal(base_size = cfg$font_base, base_family = cfg$font_family) +
    theme(
      text = element_text(face = "bold"),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "top"
    )

  out_sankey_samples <- file.path(cfg$output_dir, "S1S2_sankey_by_group_samplecount.pdf")
  write_pdf(out_sankey_samples, p_sankey_samples, width = 5.5, height = 3.2)
}

message(glue("Wrote: {out_tsv}"))
message(glue("Wrote: {out_prop}"))
message(glue("Wrote: {out_count}"))
message(glue("Wrote: {out_group_prop}"))
