#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(glue)
  library(RANN)
})

`%||%` <- function(x, y) if (is.null(x) || x == "") y else x

if (!requireNamespace("cowplot", quietly = TRUE)) {
  stop("Missing package 'cowplot'. Please install it (e.g., install.packages('cowplot')).")
}
if (!requireNamespace("ggrastr", quietly = TRUE)) {
  stop("Missing package 'ggrastr'. Please install it (e.g., install.packages('ggrastr')).")
}

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

root_dir <- get_root_dir()
source(file.path(root_dir, "utils", "R", "coarse_schema.R"))

FONT_FAMILY <- "Times"
FONT_BASE <- 9
FONT_TITLE <- 10
FONT_LEGEND_TITLE <- 9
FONT_LEGEND_TEXT <- 8

PLOT_BG_SIZE <- 0.03
PLOT_BG_ALPHA <- 0.05
PLOT_CAF_SIZE <- 0.07
PLOT_CAF_ALPHA <- 0.85
PLOT_NICHE_SIZE <- 0.025
PLOT_NICHE_ALPHA <- 0.22
RASTER_DPI <- 300

get_celltype_colors <- function(schema = c("base", "myeloid_refined")) {
  schema <- match.arg(schema)
  base <- c(
    Stromal_other = "#4D4D4D",
    CAF = "#E69F00",
    Endothelial = "#0072B2",
    Tumor = "#CC79A7",
    T_cell = "#D55E00",
    B_cell = "#6A3D9A"
  )
  myeloid <- if (schema == "myeloid_refined") {
    c(
      Macrophage = "#009E73",
      Neutrophil = "#56B4E9",
      Myeloid_other = "#984EA3"
    )
  } else {
    c(Myeloid = "#009E73")
  }
  c(base, myeloid)
}

geom_point_rast <- function(..., raster_dpi = RASTER_DPI) {
  ggrastr::geom_point_rast(..., raster.dpi = raster_dpi)
}

pad_limits <- function(x, pad = 0.02) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2])) return(c(-1, 1))
  if (rng[1] == rng[2]) return(rng + c(-1, 1))
  delta <- diff(rng) * pad
  c(rng[1] - delta, rng[2] + delta)
}

plot_pair <- function(df_bg, df_caf, df_niche, celltype_colors, title_left, title_right, xlim, ylim) {
  p_left <- ggplot() +
    geom_point_rast(data = df_bg, aes(x, y), size = PLOT_BG_SIZE, alpha = PLOT_BG_ALPHA, color = "grey50") +
    geom_point_rast(data = df_caf, aes(x, y), size = PLOT_CAF_SIZE, alpha = PLOT_CAF_ALPHA, color = "#D7301F") +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_void(base_size = FONT_BASE) +
    labs(title = title_left) +
    theme(
      text = element_text(family = FONT_FAMILY, face = "bold"),
      plot.title = element_text(size = FONT_TITLE, hjust = 0.5)
    )

  present <- intersect(names(celltype_colors), unique(df_niche$annotation_coarse))
  celltype_colors <- celltype_colors[present]

  p_right <- ggplot() +
    geom_point_rast(data = df_bg, aes(x, y), size = PLOT_BG_SIZE, alpha = PLOT_BG_ALPHA, color = "grey50") +
    geom_point_rast(
      data = df_niche,
      aes(x, y, color = annotation_coarse),
      size = PLOT_NICHE_SIZE,
      alpha = PLOT_NICHE_ALPHA
    ) +
    scale_color_manual(values = celltype_colors, drop = FALSE) +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_void(base_size = FONT_BASE) +
    labs(title = title_right, color = "annotation_coarse") +
    theme(
      text = element_text(family = FONT_FAMILY, face = "bold"),
      plot.title = element_text(size = FONT_TITLE, hjust = 0.5),
      legend.title = element_text(size = FONT_LEGEND_TITLE),
      legend.text = element_text(size = FONT_LEGEND_TEXT)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4.0, alpha = 1)))

  cowplot::plot_grid(p_left, p_right, nrow = 1, rel_widths = c(1, 1.15))
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
defaults <- list(
  seurat_root = file.path(root_dir, "data"),
  step2_cell_level_root = file.path(root_dir, "data", "caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2", "step2", "cell_level"),
  output_dir = file.path(root_dir, "data", "caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2", "step3", "subtype_niche_maps"),
  samples = "",
  subtypes = "",
  r_um = 80,
  coord_scale = 1,
  reduction = "spatial",
  anno_coarse = "annotation_coarse",
  coarse_schema = "myeloid_refined",
  max_bg = 200000L,
  max_niche = 200000L,
  seed = 1L,
  one_pdf_per_sample = TRUE,
  font_family = FONT_FAMILY
)

cfg <- list(
  seurat_root = get_opt(opts, "seurat_root", defaults$seurat_root),
  step2_cell_level_root = get_opt(opts, "step2_cell_level_root", defaults$step2_cell_level_root),
  output_dir = get_opt(opts, "output_dir", defaults$output_dir),
  samples = get_opt(opts, "samples", defaults$samples),
  subtypes = get_opt(opts, "subtypes", defaults$subtypes),
  r_um = get_opt(opts, "r_um", defaults$r_um, as.numeric),
  coord_scale = get_opt(opts, "coord_scale", defaults$coord_scale, as.numeric),
  reduction = get_opt(opts, "reduction", defaults$reduction),
  anno_coarse = get_opt(opts, "anno_coarse", defaults$anno_coarse),
  coarse_schema = get_opt(opts, "coarse_schema", defaults$coarse_schema),
  max_bg = get_opt(opts, "max_bg", defaults$max_bg, as.integer),
  max_niche = get_opt(opts, "max_niche", defaults$max_niche, as.integer),
  seed = get_opt(opts, "seed", defaults$seed, as.integer),
  one_pdf_per_sample = get_opt(opts, "one_pdf_per_sample", defaults$one_pdf_per_sample, to_bool),
  font_family = get_opt(opts, "font_family", defaults$font_family)
)

canonical <- get_canonical_coarse(cfg$coarse_schema, order = "grid")
celltype_colors <- get_celltype_colors(cfg$coarse_schema)

cell_files <- list.files(cfg$step2_cell_level_root, pattern = "_caf_cells_global_subtype\\.tsv$", full.names = TRUE)
if (length(cell_files) == 0) stop(glue("No cell-level files under {cfg$step2_cell_level_root}"))
available_samples <- sub("_caf_cells_global_subtype\\.tsv$", "", basename(cell_files))

sample_ids <- if (cfg$samples == "") available_samples else parse_csv_list(cfg$samples)
sample_ids <- intersect(sample_ids, available_samples)
if (length(sample_ids) == 0) stop("No samples selected/found.")

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

write_pair_pdf <- function(pdf_path, pages) {
  if (!dir.exists(dirname(pdf_path))) dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(pdf_path, width = 10, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (p in pages) print(p)
}

for (sample_id in sample_ids) {
  message(glue("{sample_id}: start"))
  seurat_path <- file.path(cfg$seurat_root, glue("{sample_id}.qs"))
  if (!file.exists(seurat_path)) stop(glue("{sample_id}: missing {seurat_path}"))
  cell_level_path <- file.path(cfg$step2_cell_level_root, glue("{sample_id}_caf_cells_global_subtype.tsv"))
  if (!file.exists(cell_level_path)) stop(glue("{sample_id}: missing {cell_level_path}"))

  FONT_FAMILY <<- cfg$font_family
  seu <- qs::qread(seurat_path)
  meta <- seu@meta.data
  if (!cfg$anno_coarse %in% colnames(meta)) stop(glue("{sample_id}: missing {cfg$anno_coarse}"))
  if (!cfg$reduction %in% names(seu@reductions)) stop(glue("{sample_id}: missing reduction {cfg$reduction}"))

  coords <- Seurat::Embeddings(seu, cfg$reduction)[, 1:2, drop = FALSE] * cfg$coord_scale
  colnames(coords) <- c("x", "y")

  anno <- as.character(meta[[cfg$anno_coarse]])
  unknown <- sort(unique(anno[!is.na(anno) & !anno %in% canonical]))
  if (length(unknown) > 0) stop(glue("{sample_id}: unknown {cfg$anno_coarse} values: {paste(unknown, collapse = ', ')}"))

  df_all <- data.frame(
    cell_id = rownames(meta),
    x = coords[rownames(meta), "x"],
    y = coords[rownames(meta), "y"],
    annotation_coarse = anno,
    stringsAsFactors = FALSE
  )

  df_bg <- df_all
  set.seed(cfg$seed)
  if (nrow(df_bg) > cfg$max_bg) {
    df_bg <- df_bg[sample(seq_len(nrow(df_bg)), cfg$max_bg), , drop = FALSE]
  }

  xlim <- pad_limits(df_all$x)
  ylim <- pad_limits(df_all$y)

  sub_df <- read.delim(cell_level_path, stringsAsFactors = FALSE)
  if (!"global_subtype" %in% colnames(sub_df) || !"cell_id" %in% colnames(sub_df)) {
    stop(glue("{sample_id}: cell-level file must include cell_id/global_subtype"))
  }

  subtype_levels <- unique(sub_df$global_subtype)
  subtype_levels <- subtype_levels[!is.na(subtype_levels) & subtype_levels != ""]
  subtype_levels <- sort(subtype_levels)
  if (cfg$subtypes != "") subtype_levels <- parse_csv_list(cfg$subtypes)
  if (length(subtype_levels) == 0) {
    warning(glue("{sample_id}: no subtypes found"))
    next
  }

  pages <- list()
  for (subtype in subtype_levels) {
    sub_cells <- intersect(sub_df$cell_id[sub_df$global_subtype == subtype], df_all$cell_id)
    if (length(sub_cells) == 0) next

    caf_df <- df_all[df_all$cell_id %in% sub_cells, , drop = FALSE]
    if (nrow(caf_df) == 0) next

    caf_xy <- as.matrix(caf_df[, c("x", "y"), drop = FALSE])
    all_xy <- as.matrix(df_all[, c("x", "y"), drop = FALSE])
    nn <- RANN::nn2(data = caf_xy, query = all_xy, k = 1)
    dist1 <- nn$nn.dists[, 1]
    keep <- dist1 <= cfg$r_um
    niche_df <- df_all[keep, , drop = FALSE]
    if (nrow(niche_df) > cfg$max_niche) {
      niche_df <- niche_df[sample(seq_len(nrow(niche_df)), cfg$max_niche), , drop = FALSE]
    }
    niche_df$annotation_coarse <- factor(niche_df$annotation_coarse, levels = canonical)

    pages[[length(pages) + 1]] <- plot_pair(
      df_bg = df_bg,
      df_caf = caf_df,
      df_niche = niche_df,
      celltype_colors = celltype_colors,
      title_left = glue("{subtype} CAF"),
      title_right = glue("{subtype} CAF niche (r={cfg$r_um} µm)"),
      xlim = xlim,
      ylim = ylim
    )
  }

  if (length(pages) == 0) {
    warning(glue("{sample_id}: no pages generated"))
    next
  }

  if (cfg$one_pdf_per_sample) {
    out_path <- file.path(cfg$output_dir, glue("{sample_id}_subtype_caf_niche_pairs_r{cfg$r_um}.pdf"))
    write_pair_pdf(out_path, pages)
    message(glue("{sample_id}: wrote {out_path} ({length(pages)} pages)"))
  } else {
    for (i in seq_along(pages)) {
      subtype <- subtype_levels[i]
      out_path <- file.path(cfg$output_dir, glue("{sample_id}_{subtype}_caf_vs_niche_r{cfg$r_um}.pdf"))
      write_pair_pdf(out_path, list(pages[[i]]))
      message(glue("{sample_id}: wrote {out_path}"))
    }
  }
}
