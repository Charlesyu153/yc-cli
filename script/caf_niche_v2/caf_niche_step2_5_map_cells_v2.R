#!/usr/bin/env Rscript

suppressPackageStartupMessages({
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

to_bool <- function(x) tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")

PARAM_HELP <- list(
  step1_root = "Step1 output root directory (per-sample folders with caf_clusters_selected.tsv)",
  step2_root = "Step2 output root directory (contains global_subtype_mapping.tsv)",
  output_root = "Output root directory for Step2.5 cell-level mappings",
  samples = "Comma-separated sample IDs (blank = auto-detect)",
  overwrite = "Overwrite existing outputs"
)

print_help <- function(defaults) {
  cat("CAF niche Step2.5 v2: map CAF cells to global subtypes\n")
  for (key in names(PARAM_HELP)) {
    def <- defaults[[key]]
    def_str <- if (is.na(def) || is.null(def)) "NA" else as.character(def)
    cat(glue("  --{key} {PARAM_HELP[[key]]} (default: {def_str})\n"))
  }
}

argv <- commandArgs(trailingOnly = TRUE)
opts <- parse_args(argv)

defaults <- list(
  step1_root = file.path(getwd(), "data", "caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2", "step1"),
  step2_root = file.path(getwd(), "data", "caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2", "step2", "global_corrected"),
  output_root = file.path(getwd(), "data", "caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2", "step2", "cell_level"),
  samples = "",
  overwrite = FALSE
)

if (!is.null(opts[["help"]]) || !is.null(opts[["h"]])) {
  print_help(defaults)
  quit(status = 0)
}

cfg <- list(
  step1_root = get_opt(opts, "step1_root", defaults$step1_root),
  step2_root = get_opt(opts, "step2_root", defaults$step2_root),
  output_root = get_opt(opts, "output_root", defaults$output_root),
  samples = get_opt(opts, "samples", defaults$samples),
  overwrite = get_opt(opts, "overwrite", defaults$overwrite, to_bool)
)

map_path <- file.path(cfg$step2_root, "global_subtype_mapping.tsv")
if (!file.exists(map_path)) stop(glue("Missing Step2 mapping: {map_path}"))
mapping <- read.delim(map_path, stringsAsFactors = FALSE)
if (!all(c("sample_id", "local_cluster", "global_subtype") %in% colnames(mapping))) {
  stop("global_subtype_mapping.tsv must include sample_id/local_cluster/global_subtype.")
}

dir.create(cfg$output_root, recursive = TRUE, showWarnings = FALSE)

sample_dirs <- list.dirs(cfg$step1_root, full.names = TRUE, recursive = FALSE)
sample_ids <- basename(sample_dirs)
if (cfg$samples != "") {
  keep <- trimws(strsplit(cfg$samples, ",")[[1]])
  sample_dirs <- sample_dirs[sample_ids %in% keep]
  sample_ids <- basename(sample_dirs)
}
sample_dirs <- sample_dirs[order(sample_ids)]
if (length(sample_dirs) == 0) stop("No sample directories found under step1_root.")

for (sample_dir in sample_dirs) {
  sample_id <- basename(sample_dir)
  cl_path <- file.path(sample_dir, "caf_clusters_selected.tsv")
  if (!file.exists(cl_path)) stop(glue("{sample_id}: missing {cl_path}"))

  map_sub <- mapping[mapping$sample_id == sample_id, , drop = FALSE]
  if (nrow(map_sub) == 0) stop(glue("{sample_id}: no rows in global_subtype_mapping.tsv"))
  lc_to_sub <- stats::setNames(map_sub$global_subtype, map_sub$local_cluster)

  cl_df <- read.delim(cl_path, stringsAsFactors = FALSE)
  if (!"cell_id" %in% colnames(cl_df)) {
    stop(glue("{sample_id}: caf_clusters_selected.tsv must include cell_id."))
  }
  lc_col <- if ("local_cluster" %in% colnames(cl_df)) {
    "local_cluster"
  } else if ("caf_cluster" %in% colnames(cl_df)) {
    "caf_cluster"
  } else {
    stop(glue("{sample_id}: caf_clusters_selected.tsv must include local_cluster or caf_cluster."))
  }
  cl_df$local_cluster <- cl_df[[lc_col]]
  cl_df$global_subtype <- lc_to_sub[cl_df$local_cluster]
  cl_df$sample_id <- sample_id

  out_path <- file.path(cfg$output_root, glue("{sample_id}_caf_cells_global_subtype.tsv"))
  if (file.exists(out_path) && !cfg$overwrite) {
    message(glue("{sample_id}: skip (exists) {out_path}"))
    next
  }
  utils::write.table(
    cl_df[, c("cell_id", "sample_id", "local_cluster", "global_subtype")],
    out_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  message(glue("{sample_id}: wrote {out_path}"))
}

message(glue("Done. Output: {cfg$output_root}"))
