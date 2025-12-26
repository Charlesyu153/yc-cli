#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glue)
  library(qs)
  library(Seurat)
})

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

to_bool <- function(x) {
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

get_opt <- function(opts, key, default, cast = NULL) {
  val <- opts[[key]]
  if (is.null(val) || val == "") return(default)
  if (is.null(cast)) return(val)
  cast(val)
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
source(file.path(root_dir, "utils", "R", "coarse_schema.R"))

PARAM_HELP <- list(
  input = "Input directory with Seurat .qs files",
  pattern = "Regex pattern for .qs files",
  samples = "Comma-separated sample IDs (blank = all)",
  include_conversion_summary = "Also update conversion_summary.qs (true/false)",
  anno_coarse = "Column name to overwrite (default: annotation_coarse)",
  anno_coarse_old = "Backup column name to create/preserve (default: annotation_coarse_old)",
  myeloid_fine_cols = "Comma-separated meta columns used to refine Myeloid",
  in_place = "Overwrite input files (true/false)",
  output = "Output directory when in_place=false",
  dry_run = "Print changes but do not write (true/false)"
)

print_help <- function(defaults) {
  cat("Update annotation_coarse: backup old + refine Myeloid\n")
  for (key in names(PARAM_HELP)) {
    def <- defaults[[key]]
    def_str <- if (is.na(def) || is.null(def)) "NA" else as.character(def)
    cat(glue("  --{key} {PARAM_HELP[[key]]} (default: {def_str})\n"))
  }
}

argv <- commandArgs(trailingOnly = TRUE)
opts <- parse_args(argv)

defaults <- list(
  input = file.path(root_dir, "data"),
  pattern = "\\.qs$",
  samples = "",
  include_conversion_summary = FALSE,
  anno_coarse = "annotation_coarse",
  anno_coarse_old = "annotation_coarse_old",
  myeloid_fine_cols = "annotation_fine,cell_type_lvl2,cell_type_lvl1,annotation",
  in_place = TRUE,
  output = file.path(root_dir, "data", "_updated_qs"),
  dry_run = FALSE
)

if (!is.null(opts$help) || !is.null(opts$h)) {
  print_help(defaults)
  quit(status = 0)
}

cfg <- list(
  input = get_opt(opts, "input", defaults$input),
  pattern = get_opt(opts, "pattern", defaults$pattern),
  samples = get_opt(opts, "samples", defaults$samples),
  include_conversion_summary = get_opt(opts, "include_conversion_summary", defaults$include_conversion_summary, to_bool),
  anno_coarse = get_opt(opts, "anno_coarse", defaults$anno_coarse),
  anno_coarse_old = get_opt(opts, "anno_coarse_old", defaults$anno_coarse_old),
  myeloid_fine_cols = get_opt(opts, "myeloid_fine_cols", defaults$myeloid_fine_cols),
  in_place = get_opt(opts, "in_place", defaults$in_place, to_bool),
  output = get_opt(opts, "output", defaults$output),
  dry_run = get_opt(opts, "dry_run", defaults$dry_run, to_bool)
)

cfg$myeloid_fine_cols <- parse_csv_list(cfg$myeloid_fine_cols)
if (length(cfg$myeloid_fine_cols) == 0) {
  cfg$myeloid_fine_cols <- default_myeloid_fine_cols()
}

input_files <- list.files(cfg$input, pattern = cfg$pattern, full.names = TRUE)
if (!cfg$include_conversion_summary) {
  input_files <- input_files[basename(input_files) != "conversion_summary.qs"]
}

if (cfg$samples != "") {
  sample_list <- trimws(strsplit(cfg$samples, ",")[[1]])
  keep <- paste0(sample_list, ".qs")
  input_files <- input_files[basename(input_files) %in% keep]
}

if (length(input_files) == 0) {
  stop("No input files found.")
}

if (!cfg$in_place) {
  dir.create(cfg$output, recursive = TRUE, showWarnings = FALSE)
}

canonical_levels <- get_canonical_coarse("myeloid_refined", order = "analysis")

process_one <- function(file) {
  sample_id <- sub("\\.qs$", "", basename(file))
  cat(glue("{sample_id}: reading\n"))

  seu <- qs::qread(file)
  if (!inherits(seu, "Seurat")) {
    warning(glue("{sample_id}: not a Seurat object; skip."))
    return(invisible(NULL))
  }

  meta <- seu@meta.data
  if (!cfg$anno_coarse %in% colnames(meta)) {
    warning(glue("{sample_id}: missing {cfg$anno_coarse}; skip."))
    return(invisible(NULL))
  }

  if (cfg$anno_coarse_old %in% colnames(meta)) {
    coarse_source <- as.character(meta[[cfg$anno_coarse_old]])
  } else {
    coarse_source <- as.character(meta[[cfg$anno_coarse]])
    meta[[cfg$anno_coarse_old]] <- coarse_source
  }

  meta_tmp <- meta
  meta_tmp[[cfg$anno_coarse]] <- coarse_source
  refined <- refine_myeloid_in_coarse(
    meta_tmp,
    anno_coarse = cfg$anno_coarse,
    fine_cols = cfg$myeloid_fine_cols
  )

  meta[[cfg$anno_coarse]] <- factor(refined, levels = canonical_levels)
  seu@meta.data <- meta

  before <- coarse_source
  after <- as.character(meta[[cfg$anno_coarse]])
  changed <- sum(!is.na(before) & !is.na(after) & before != after)

  cat(glue("{sample_id}: updated (changed={changed})\n"))
  if (cfg$dry_run) return(invisible(NULL))

  out_path <- if (cfg$in_place) file else file.path(cfg$output, basename(file))
  tmp_path <- file.path(dirname(out_path), paste0(basename(out_path), ".tmp"))
  qs::qsave(seu, tmp_path)
  ok <- file.rename(tmp_path, out_path)
  if (!ok) {
    unlink(tmp_path)
    stop(glue("{sample_id}: failed to write {out_path}"))
  }
  invisible(NULL)
}

for (i in seq_along(input_files)) {
  file <- input_files[i]
  cat(glue("[{i}/{length(input_files)}] "))
  process_one(file)
  gc()
}

cat("Done.\n")

