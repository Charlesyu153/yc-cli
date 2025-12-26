#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(qs)
  library(ggplot2)
  library(glue)
  library(ggsci)
})

CANONICAL_COARSE <- c(
  "B_cell",
  "CAF",
  "Endothelial",
  "Macrophage",
  "Neutrophil",
  "Myeloid_Mast",
  "Myeloid_other",
  "Stromal_other",
  "T_cell",
  "Tumor"
)

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

get_opt <- function(opts, key, default, cast = NULL) {
  val <- opts[[key]]
  if (is.null(val) || val == "") return(default)
  if (is.null(cast)) return(val)
  cast(val)
}

PARAM_HELP <- list(
  input = "Input directory with Seurat .qs files",
  output = "Output directory for coarse composition summaries",
  samples = "Comma-separated sample IDs (blank = all)",
  pattern = "Regex pattern for .qs files",
  anno_coarse = "Coarse annotation column name"
)

print_help <- function(defaults) {
  cat("Coarse annotation composition summary - parameters\n")
  for (key in names(PARAM_HELP)) {
    def <- defaults[[key]]
    def_str <- if (is.na(def) || is.null(def)) "NA" else as.character(def)
    cat(glue("  --{key} {PARAM_HELP[[key]]} (default: {def_str})\n"))
  }
}

argv <- commandArgs(trailingOnly = TRUE)
opts <- parse_args(argv)

defaults <- list(
  input = file.path(getwd(), "data"),
  output = file.path(getwd(), "data", "coarse_composition"),
  samples = "",
  pattern = "\\.qs$",
  anno_coarse = "annotation_coarse"
)

if (!is.null(opts$help) || !is.null(opts$h)) {
  print_help(defaults)
  quit(status = 0)
}

cfg <- list(
  input = get_opt(opts, "input", defaults$input),
  output = get_opt(opts, "output", defaults$output),
  samples = get_opt(opts, "samples", defaults$samples),
  pattern = get_opt(opts, "pattern", defaults$pattern),
  anno_coarse = get_opt(opts, "anno_coarse", defaults$anno_coarse)
)

input_files <- list.files(cfg$input, pattern = cfg$pattern, full.names = TRUE)
if (cfg$samples != "") {
  sample_list <- trimws(strsplit(cfg$samples, ",")[[1]])
  input_files <- input_files[basename(input_files) %in% paste0(sample_list, ".qs")]
}
input_files <- input_files[basename(input_files) != "conversion_summary.qs"]

if (length(input_files) == 0) {
  stop("No input files found.")
}

dir.create(cfg$output, recursive = TRUE, showWarnings = FALSE)

counts_list <- list()
levels_seen <- character(0)

for (file in input_files) {
  sample_id <- sub("\\.qs$", "", basename(file))
  cat(glue("{sample_id}: reading\n"))
  seu <- qs::qread(file)
  if (!cfg$anno_coarse %in% colnames(seu@meta.data)) {
    warning(glue("{sample_id}: missing {cfg$anno_coarse}; skip."))
    next
  }
  anno <- as.character(seu@meta.data[[cfg$anno_coarse]])
  anno <- anno[!is.na(anno)]
  tab <- table(anno)
  counts_list[[sample_id]] <- tab
  levels_seen <- unique(c(levels_seen, names(tab)))
  rm(seu)
  gc()
}

if (length(counts_list) == 0) {
  stop("No samples with annotation_coarse found.")
}

extras <- setdiff(levels_seen, CANONICAL_COARSE)
level_order <- c(CANONICAL_COARSE, sort(extras))
level_order <- level_order[level_order %in% levels_seen]

sample_ids <- names(counts_list)
count_mat <- matrix(0, nrow = length(sample_ids), ncol = length(level_order))
rownames(count_mat) <- sample_ids
colnames(count_mat) <- level_order

for (sample_id in sample_ids) {
  tab <- counts_list[[sample_id]]
  count_mat[sample_id, names(tab)] <- as.integer(tab)
}

total_counts <- rowSums(count_mat)
frac_mat <- count_mat / ifelse(total_counts == 0, NA_real_, total_counts)

utils::write.table(
  count_mat,
  file.path(cfg$output, "coarse_counts.tsv"),
  sep = "\t",
  quote = FALSE
)
utils::write.table(
  frac_mat,
  file.path(cfg$output, "coarse_fractions.tsv"),
  sep = "\t",
  quote = FALSE
)

long_df <- data.frame(
  sample_id = rep(rownames(count_mat), times = ncol(count_mat)),
  cell_type = rep(colnames(count_mat), each = nrow(count_mat)),
  n = as.vector(count_mat),
  fraction = as.vector(frac_mat)
)

utils::write.table(
  long_df,
  file.path(cfg$output, "coarse_composition.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_df <- long_df
plot_df$sample_id <- factor(plot_df$sample_id, levels = sample_ids)
plot_df$cell_type <- factor(plot_df$cell_type, levels = level_order)

p_fraction <- ggplot(plot_df, aes(sample_id, fraction, fill = cell_type)) +
  geom_col(width = 0.85) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Fraction", fill = cfg$anno_coarse, title = "Coarse annotation fractions")

p_counts <- ggplot(plot_df, aes(sample_id, n, fill = cell_type)) +
  geom_col(width = 0.85) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Cell count", fill = cfg$anno_coarse, title = "Coarse annotation counts")

p_fraction <- p_fraction + ggsci::scale_fill_d3(
  "category20",
  limits = level_order,
  breaks = level_order,
  drop = FALSE
)
p_counts <- p_counts + ggsci::scale_fill_d3(
  "category20",
  limits = level_order,
  breaks = level_order,
  drop = FALSE
)

ggsave(
  file.path(cfg$output, "coarse_composition_fraction.pdf"),
  p_fraction,
  width = 10,
  height = 5,
  device = grDevices::cairo_pdf
)
ggsave(
  file.path(cfg$output, "coarse_composition_counts.pdf"),
  p_counts,
  width = 10,
  height = 5,
  device = grDevices::cairo_pdf
)

cat(glue("Done. Output: {cfg$output}\n"))
