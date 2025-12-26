#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(ggplot2)
  library(glue)
  library(dbscan)
})

`%||%` <- function(x, y) if (is.null(x) || x == "") y else x

if (!requireNamespace("cowplot", quietly = TRUE)) {
  stop("Missing package 'cowplot'. Please install it (e.g., install.packages('cowplot')).")
}

COORD_SCALE <- 1
RNG_SEED <- 1024
RASTER_DPI <- 300
PLOT_BG_SIZE <- 0.06
PLOT_BG_ALPHA <- 0.08
PLOT_FG_SIZE <- 0.14
PLOT_FG_ALPHA <- 0.9

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
      Myeloid_Mast = "#F0E442",
      Myeloid_other = "#00A087"
    )
  } else {
    c(Myeloid = "#009E73")
  }
  c(base, myeloid)
}

geom_point_layer <- function(..., raster_dpi = RASTER_DPI) {
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    return(ggrastr::geom_point_rast(..., raster.dpi = raster_dpi))
  }
  ggplot2::geom_point(...)
}

sort_subtypes <- function(subtypes) {
  key <- suppressWarnings(as.integer(sub("^S([0-9]+)$", "\\1", subtypes)))
  if (all(is.na(key))) return(sort(subtypes))
  subtypes[order(key, subtypes)]
}

make_blank_panel <- function(title = "") {
  ggplot() +
    theme_void(base_size = 12) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0, face = "bold", size = 12),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      aspect.ratio = 1
    ) +
    labs(title = title)
}

pad_limits <- function(x, pad = 0.03) {
  if (length(x) == 0 || all(is.na(x))) return(c(0, 1))
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2])) return(c(0, 1))
  if (rng[1] == rng[2]) return(rng + c(-1, 1))
  delta <- diff(rng) * pad
  c(rng[1] - delta, rng[2] + delta)
}

apply_limits <- function(p, xlim, ylim) {
  p + coord_equal(xlim = xlim, ylim = ylim)
}

plot_subtype_panel <- function(df_bg, caf_df, subtype_levels, subtype_colors, title, xlim, ylim) {
  p <- ggplot() +
    geom_point_layer(
      data = df_bg,
      aes(x, y),
      size = PLOT_BG_SIZE,
      alpha = PLOT_BG_ALPHA,
      color = "#BDBDBD"
    ) +
    geom_point_layer(
      data = caf_df,
      aes(x, y, color = subtype),
      size = PLOT_FG_SIZE,
      alpha = PLOT_FG_ALPHA
    ) +
    scale_color_manual(values = subtype_colors, breaks = subtype_levels, drop = FALSE) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 12),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      legend.position = "none",
      aspect.ratio = 1
    ) +
    labs(title = title)
  apply_limits(p, xlim, ylim)
}

plot_neighbor_panel <- function(df_bg, niche_df, celltype_colors, title, xlim, ylim) {
  p <- ggplot() +
    geom_point_layer(
      data = df_bg,
      aes(x, y),
      size = PLOT_BG_SIZE,
      alpha = PLOT_BG_ALPHA,
      color = "#BDBDBD"
    ) +
    geom_point_layer(
      data = niche_df,
      aes(x, y, color = cell_type),
      size = PLOT_FG_SIZE,
      alpha = PLOT_FG_ALPHA
    ) +
    scale_color_manual(values = celltype_colors, drop = FALSE) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 11),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      legend.position = "none",
      aspect.ratio = 1
    ) +
    labs(title = title)
  apply_limits(p, xlim, ylim)
}

plot_pie_panel <- function(frac, celltype_colors, title) {
  df <- data.frame(cell_type = names(frac), frac = as.numeric(frac), stringsAsFactors = FALSE)
  df <- df[df$frac > 0, , drop = FALSE]
  if (nrow(df) == 0) return(make_blank_panel(title))
  ggplot(df, aes(x = "", y = frac, fill = cell_type)) +
    geom_col(width = 1, color = NA) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = celltype_colors, drop = FALSE) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 11),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      legend.position = "none"
    ) +
    labs(title = title)
}

compute_neighbors <- function(coords, query_cells, r) {
  query_cells <- intersect(query_cells, rownames(coords))
  if (length(query_cells) == 0) return(list())
  Q <- coords[query_cells, , drop = FALSE]
  nbrs <- dbscan::frNN(coords, eps = r, query = Q, sort = FALSE)
  nn <- lapply(nbrs$id, function(id) rownames(coords)[id])
  names(nn) <- rownames(Q)
  nn
}

pick_example_samples <- function(cell_level_root, subtype_levels) {
  files <- list.files(cell_level_root, pattern = "_caf_cells_global_subtype\\.tsv$", full.names = TRUE)
  if (length(files) == 0) stop(glue("No cell-level files found under {cell_level_root}"))
  sample_ids <- sub("_caf_cells_global_subtype\\.tsv$", "", basename(files))
  sample_ids <- sort(unique(sample_ids))
  cov <- sapply(sample_ids, function(s) {
    path <- file.path(cell_level_root, glue("{s}_caf_cells_global_subtype.tsv"))
    df <- read.delim(path, stringsAsFactors = FALSE)
    length(intersect(unique(df$global_subtype), subtype_levels))
  })
  full <- names(cov)[which.max(cov)][1]
  missing <- names(cov)[which.min(cov)][1]
  c(full, missing)
}

argv <- commandArgs(trailingOnly = TRUE)
opts <- parse_args(argv)

defaults <- list(
  seurat_root = file.path(getwd(), "data"),
  step1_root = file.path(getwd(), "data", "caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2", "step1"),
  step2_root = file.path(getwd(), "data", "caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2", "step2", "global_corrected"),
  cell_level_root = file.path(getwd(), "data", "caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2", "step2", "cell_level"),
  output = "",
  samples = "",
  subtypes = "",
  example_auto = TRUE,
  anno_coarse = "annotation_coarse",
  coarse_schema = "myeloid_refined",
  reduction = "spatial",
  r_fixed = 80,
  coord_scale = COORD_SCALE,
  max_plot_all_cells = 200000L,
  max_plot_niche_cells = 200000L,
  use_r_choice = TRUE
)

cfg <- list(
  seurat_root = get_opt(opts, "seurat_root", defaults$seurat_root),
  step1_root = get_opt(opts, "step1_root", defaults$step1_root),
  step2_root = get_opt(opts, "step2_root", defaults$step2_root),
  cell_level_root = get_opt(opts, "cell_level_root", defaults$cell_level_root),
  output = get_opt(opts, "output", defaults$output),
  samples = get_opt(opts, "samples", defaults$samples),
  subtypes = get_opt(opts, "subtypes", defaults$subtypes),
  example_auto = get_opt(opts, "example_auto", defaults$example_auto, to_bool),
  anno_coarse = get_opt(opts, "anno_coarse", defaults$anno_coarse),
  coarse_schema = tolower(get_opt(opts, "coarse_schema", defaults$coarse_schema)),
  reduction = get_opt(opts, "reduction", defaults$reduction),
  r_fixed = get_opt(opts, "r_fixed", defaults$r_fixed, as.numeric),
  coord_scale = get_opt(opts, "coord_scale", defaults$coord_scale, as.numeric),
  max_plot_all_cells = get_opt(opts, "max_plot_all_cells", defaults$max_plot_all_cells, as.integer),
  max_plot_niche_cells = get_opt(opts, "max_plot_niche_cells", defaults$max_plot_niche_cells, as.integer),
  use_r_choice = get_opt(opts, "use_r_choice", defaults$use_r_choice, to_bool)
)

if (!cfg$coarse_schema %in% c("base", "myeloid_refined")) {
  stop(glue("Unknown coarse_schema: {cfg$coarse_schema}. Use base or myeloid_refined."))
}
canonical <- get_canonical_coarse(cfg$coarse_schema, order = "grid")
celltype_colors <- get_celltype_colors(cfg$coarse_schema)

map_path <- file.path(cfg$step2_root, "global_subtype_mapping.tsv")
if (!file.exists(map_path)) stop(glue("Missing Step2 mapping: {map_path}"))
map_df <- read.delim(map_path, stringsAsFactors = FALSE)
if (!"global_subtype" %in% colnames(map_df)) stop("global_subtype_mapping.tsv missing global_subtype column.")
subtype_levels <- sort_subtypes(unique(map_df$global_subtype))
subtype_levels <- subtype_levels[!is.na(subtype_levels) & subtype_levels != ""]
if (cfg$subtypes != "") {
  subtype_levels <- parse_csv_list(cfg$subtypes)
}
if (length(subtype_levels) == 0) stop("No global subtypes to plot.")

if (cfg$example_auto) {
  sample_ids <- pick_example_samples(cfg$cell_level_root, subtype_levels)
} else {
  sample_ids <- parse_csv_list(cfg$samples)
}
if (length(sample_ids) == 0) stop("No samples selected.")

subtype_colors <- stats::setNames(ggsci::pal_d3("category20")(max(6, length(subtype_levels) + 1)), c(subtype_levels, "Unmapped"))

read_r_choice <- function(sample_dir, default_r, use_r_choice) {
  if (!use_r_choice) return(default_r)
  path <- file.path(sample_dir, "r_choice.tsv")
  if (!file.exists(path)) return(default_r)
  df <- read.delim(path, stringsAsFactors = FALSE)
  if (!"r_fixed" %in% colnames(df)) return(default_r)
  as.numeric(df$r_fixed[1]) %||% default_r
}

build_sample_row <- function(sample_id) {
  sample_step1 <- file.path(cfg$step1_root, sample_id)
  if (!dir.exists(sample_step1)) stop(glue("{sample_id}: missing Step1 dir {sample_step1}"))
  seurat_path <- file.path(cfg$seurat_root, glue("{sample_id}.qs"))
  if (!file.exists(seurat_path)) stop(glue("{sample_id}: missing Seurat file {basename(seurat_path)}"))
  subtype_path <- file.path(cfg$cell_level_root, glue("{sample_id}_caf_cells_global_subtype.tsv"))
  if (!file.exists(subtype_path)) stop(glue("{sample_id}: missing subtype file {basename(subtype_path)}"))

  seu <- qs::qread(seurat_path)
  meta <- seu@meta.data
  if (!cfg$anno_coarse %in% colnames(meta)) stop(glue("{sample_id}: missing {cfg$anno_coarse}"))
  if (!cfg$reduction %in% names(seu@reductions)) stop(glue("{sample_id}: missing reduction {cfg$reduction}"))
  if (!all(c("Center_X", "Center_Y") %in% colnames(meta))) stop(glue("{sample_id}: missing Center_X/Center_Y"))

  anno <- as.character(meta[[cfg$anno_coarse]])
  unknown <- sort(unique(anno[!is.na(anno) & !anno %in% canonical]))
  if (length(unknown) > 0) stop(glue("{sample_id}: unknown {cfg$anno_coarse} values: {paste(unknown, collapse = ', ')}"))

  df_all_full <- data.frame(
    cell_id = rownames(meta),
    x = meta$Center_X * cfg$coord_scale,
    y = meta$Center_Y * cfg$coord_scale,
    cell_type = anno,
    stringsAsFactors = FALSE
  )

  df_bg <- df_all_full
  if (nrow(df_bg) > cfg$max_plot_all_cells) {
    set.seed(RNG_SEED)
    df_bg <- df_bg[sample(seq_len(nrow(df_bg)), cfg$max_plot_all_cells), , drop = FALSE]
  }
  xlim <- pad_limits(df_all_full$x)
  ylim <- pad_limits(df_all_full$y)

  sub_df <- read.delim(subtype_path, stringsAsFactors = FALSE)
  if (!all(c("cell_id", "global_subtype") %in% colnames(sub_df))) {
    stop(glue("{sample_id}: subtype file must include cell_id/global_subtype"))
  }
  sub_df$subtype <- sub_df$global_subtype
  sub_df$subtype[is.na(sub_df$subtype) | sub_df$subtype == ""] <- "Unmapped"
  caf_df <- merge(sub_df[, c("cell_id", "subtype"), drop = FALSE], df_all_full, by = "cell_id")

  r_use <- read_r_choice(sample_step1, cfg$r_fixed, cfg$use_r_choice)
  coords <- Seurat::Embeddings(seu, reduction = cfg$reduction) * cfg$coord_scale
  caf_cells_all <- unique(sub_df$cell_id)
  neighbors <- compute_neighbors(coords, caf_cells_all, r_use)

  p_left <- if (nrow(caf_df) == 0) {
    make_blank_panel(sample_id)
  } else {
    plot_subtype_panel(df_bg, caf_df, c(subtype_levels, "Unmapped"), subtype_colors, sample_id, xlim, ylim)
  }

  neighbor_panels <- list()
  pie_panels <- list()

  subtype_for_cell <- stats::setNames(sub_df$subtype, sub_df$cell_id)
  subtype_for_cell <- subtype_for_cell[unique(names(subtype_for_cell))]

  for (subtype in subtype_levels) {
    caf_cells <- names(subtype_for_cell)[subtype_for_cell == subtype]
    caf_cells <- intersect(caf_cells, names(neighbors))
    if (length(caf_cells) == 0) {
      neighbor_panels[[subtype]] <- make_blank_panel(subtype)
      pie_panels[[subtype]] <- make_blank_panel(subtype)
      next
    }

    niche_cells <- unique(unlist(neighbors[caf_cells], use.names = FALSE))
    niche_cells <- intersect(niche_cells, df_all_full$cell_id)
    if (length(niche_cells) == 0) {
      neighbor_panels[[subtype]] <- make_blank_panel(subtype)
      pie_panels[[subtype]] <- make_blank_panel(subtype)
      next
    }

    niche_df <- df_all_full[df_all_full$cell_id %in% niche_cells, c("x", "y", "cell_type"), drop = FALSE]
    if (nrow(niche_df) > cfg$max_plot_niche_cells) {
      set.seed(RNG_SEED)
      niche_df <- niche_df[sample(seq_len(nrow(niche_df)), cfg$max_plot_niche_cells), , drop = FALSE]
    }
    niche_df$cell_type <- factor(niche_df$cell_type, levels = canonical)

    neighbor_panels[[subtype]] <- plot_neighbor_panel(df_bg, niche_df, celltype_colors, subtype, xlim, ylim)

    counts <- table(niche_df$cell_type)
    frac <- as.numeric(counts) / sum(counts)
    names(frac) <- names(counts)
    frac <- frac[!is.na(names(frac))]
    pie_panels[[subtype]] <- plot_pie_panel(frac, celltype_colors, subtype)
  }

  neighbor_grid <- cowplot::plot_grid(plotlist = neighbor_panels, nrow = 1, align = "hv")
  pie_grid <- cowplot::plot_grid(plotlist = pie_panels, ncol = 1, align = "hv")

  row <- cowplot::plot_grid(
    p_left,
    neighbor_grid,
    pie_grid,
    nrow = 1,
    rel_widths = c(1.1, length(subtype_levels) * 1.0, 0.9)
  )
  rm(seu)
  gc()
  row
}

rows <- lapply(sample_ids, build_sample_row)
panel <- cowplot::plot_grid(plotlist = rows, ncol = 1, align = "v", rel_heights = rep(1, length(rows)))

header <- cowplot::ggdraw() +
  cowplot::draw_label("CAF subtype", x = 0.02, y = 0.98, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  cowplot::draw_label("Neighboring cells", x = 0.35, y = 0.98, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  cowplot::draw_label("Neighboring composition", x = 0.84, y = 0.98, hjust = 0, vjust = 1, fontface = "bold", size = 14)

final <- cowplot::plot_grid(header, panel, ncol = 1, rel_heights = c(0.12, 1))

if (cfg$output == "") {
  cfg$output <- file.path(cfg$step2_root, "niche_by_sample", "example_niche_by_sample_grid_custom_v2.pdf")
}
dir.create(dirname(cfg$output), recursive = TRUE, showWarnings = FALSE)

png_width <- as.integer(360 * (length(subtype_levels) + 3))
png_height <- as.integer(360 * length(sample_ids) + 180)

ext <- tolower(tools::file_ext(cfg$output))
if (ext == "pdf") {
  w_in <- png_width / 300
  h_in <- png_height / 300
  grDevices::cairo_pdf(cfg$output, width = w_in, height = h_in)
  print(final)
  grDevices::dev.off()
  message(glue("Done. Output: {cfg$output}"))
} else if (ext == "png") {
  if (!requireNamespace("ragg", quietly = TRUE)) {
    stop("Need R package 'ragg' to write PNG output. Please install it.")
  }
  ragg::agg_png(cfg$output, width = png_width, height = png_height, res = 300)
  print(final)
  grDevices::dev.off()
  message(glue("Done. Output: {cfg$output}"))
} else {
  stop("output must end with .pdf or .png")
}
