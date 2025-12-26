#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glue)
  library(qs)
  library(Seurat)
})

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
  include_conversion_summary = "Also rename conversion_summary.qs (true/false)",
  anno_coarse = "Coarse annotation column to rename (default: annotation_coarse)",
  in_place = "Overwrite input files (true/false)",
  output = "Output directory when in_place=false",
  dry_run = "Print summary but do not write (true/false)"
)

print_help <- function(defaults) {
  cat("Rename myeloid refined labels in annotation_coarse\n")
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
  in_place = get_opt(opts, "in_place", defaults$in_place, to_bool),
  output = get_opt(opts, "output", defaults$output),
  dry_run = get_opt(opts, "dry_run", defaults$dry_run, to_bool)
)

input_files <- list.files(cfg$input, pattern = cfg$pattern, full.names = TRUE)
if (!cfg$include_conversion_summary) {
  input_files <- input_files[basename(input_files) != "conversion_summary.qs"]
}
if (cfg$samples != "") {
  sample_list <- trimws(strsplit(cfg$samples, ",")[[1]])
  keep <- paste0(sample_list, ".qs")
  input_files <- input_files[basename(input_files) %in% keep]
}
if (length(input_files) == 0) stop("No input files found.")

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

  before <- as.character(meta[[cfg$anno_coarse]])
  after <- rename_myeloid_refined_labels(before)
  changed <- sum(!is.na(before) & !is.na(after) & before != after)

  meta[[cfg$anno_coarse]] <- factor(after, levels = canonical_levels)
  seu@meta.data <- meta

  cat(glue("{sample_id}: renamed (changed={changed})\n"))
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
  cat(glue("[{i}/{length(input_files)}] "))
  process_one(input_files[[i]])
  gc()
}

cat("Done.\n")

