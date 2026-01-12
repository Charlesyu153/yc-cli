#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(glue)
  library(RANN)
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

opts <- parse_args(commandArgs(trailingOnly = TRUE))
root_dir <- get_root_dir()

defaults <- list(
  seurat_root = file.path(root_dir, "data"),
  step2_cell_level_root = file.path(root_dir, "data", "caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2", "step2", "cell_level"),
  output_dir = file.path(root_dir, "data", "caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2", "step3", "S1S2_group_stats"),
  subtype = "S1",
  anno_coarse = "annotation_coarse",
  reduction = "spatial",
  r_um = 80,
  coord_scale = 1,
  include_caf_in_denom = TRUE,
  bcell_label = "B_cell",
  caf_label = "CAF"
)

cfg <- list(
  seurat_root = get_opt(opts, "seurat_root", defaults$seurat_root),
  step2_cell_level_root = get_opt(opts, "step2_cell_level_root", defaults$step2_cell_level_root),
  output_dir = get_opt(opts, "output_dir", defaults$output_dir),
  subtype = get_opt(opts, "subtype", defaults$subtype),
  anno_coarse = get_opt(opts, "anno_coarse", defaults$anno_coarse),
  reduction = get_opt(opts, "reduction", defaults$reduction),
  r_um = get_opt(opts, "r_um", defaults$r_um, as.numeric),
  coord_scale = get_opt(opts, "coord_scale", defaults$coord_scale, as.numeric),
  include_caf_in_denom = get_opt(opts, "include_caf_in_denom", defaults$include_caf_in_denom, to_bool),
  bcell_label = get_opt(opts, "bcell_label", defaults$bcell_label),
  caf_label = get_opt(opts, "caf_label", defaults$caf_label)
)

cell_files <- list.files(cfg$step2_cell_level_root, pattern = "_caf_cells_global_subtype\\.tsv$", full.names = TRUE)
if (length(cell_files) == 0) stop(glue("No cell-level files under {cfg$step2_cell_level_root}"))
sample_ids <- sub("_caf_cells_global_subtype\\.tsv$", "", basename(cell_files))

rows <- list()
for (i in seq_along(cell_files)) {
  sample_id <- sample_ids[i]
  group <- infer_group(sample_id)
  if (is.na(group)) next

  seurat_path <- file.path(cfg$seurat_root, glue("{sample_id}.qs"))
  if (!file.exists(seurat_path)) {
    warning(glue("{sample_id}: missing Seurat file {seurat_path}; skip"))
    next
  }

  sub_df <- read.delim(cell_files[i], stringsAsFactors = FALSE)
  if (!all(c("cell_id", "global_subtype") %in% colnames(sub_df))) {
    warning(glue("{sample_id}: cell-level file missing columns; skip"))
    next
  }

  caf_cells <- intersect(sub_df$cell_id[sub_df$global_subtype == cfg$subtype], sub_df$cell_id)
  caf_cells <- unique(caf_cells)
  if (length(caf_cells) == 0) next

  seu <- qs::qread(seurat_path)
  meta <- seu@meta.data
  if (!cfg$anno_coarse %in% colnames(meta)) stop(glue("{sample_id}: missing {cfg$anno_coarse}"))
  if (!cfg$reduction %in% names(seu@reductions)) stop(glue("{sample_id}: missing reduction {cfg$reduction}"))

  coords <- Seurat::Embeddings(seu, cfg$reduction)[, 1:2, drop = FALSE] * cfg$coord_scale
  caf_cells <- intersect(caf_cells, rownames(coords))
  if (length(caf_cells) == 0) next

  caf_xy <- as.matrix(coords[caf_cells, , drop = FALSE])
  all_xy <- as.matrix(coords[rownames(meta), , drop = FALSE])
  nn <- RANN::nn2(data = caf_xy, query = all_xy, k = 1)
  dist1 <- nn$nn.dists[, 1]
  keep <- dist1 <= cfg$r_um

  anno <- as.character(meta[[cfg$anno_coarse]])
  names(anno) <- rownames(meta)
  anno_keep <- anno[keep]

  denom_mask <- rep(TRUE, length(anno_keep))
  if (!cfg$include_caf_in_denom) denom_mask <- anno_keep != cfg$caf_label
  denom <- sum(denom_mask & !is.na(anno_keep))
  b_n <- sum(denom_mask & anno_keep == cfg$bcell_label, na.rm = TRUE)
  frac <- if (denom > 0) b_n / denom else NA_real_

  caf_n <- length(caf_cells)
  niche_n <- sum(keep)
  rows[[length(rows) + 1]] <- data.frame(
    sample_id = sample_id,
    group = group,
    subtype = cfg$subtype,
    r_um = cfg$r_um,
    coord_scale = cfg$coord_scale,
    n_caf_subtype = caf_n,
    n_niche_cells = niche_n,
    denom_niche = denom,
    b_cell_n = b_n,
    b_cell_frac = frac,
    include_caf_in_denom = cfg$include_caf_in_denom,
    stringsAsFactors = FALSE
  )
  rm(seu)
  gc()
}

out_df <- do.call(rbind, rows)
if (is.null(out_df) || nrow(out_df) == 0) stop(glue("No samples with subtype {cfg$subtype}."))

out_df$group <- factor(out_df$group, levels = c("Primary", "Metastasis"))
out_df <- out_df[order(out_df$group, out_df$sample_id), , drop = FALSE]

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
suffix <- if (cfg$include_caf_in_denom) "inclCAF" else "exclCAF"
out_path <- file.path(cfg$output_dir, glue("{cfg$subtype}_Bcell_in_niche_by_sample_r{cfg$r_um}_{suffix}.tsv"))
utils::write.table(out_df, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

message(glue("Done. Output: {out_path}"))
