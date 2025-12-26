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

COORD_SCALE <- 50
RNG_SEED <- 1024
RASTER_DPI <- 300
PLOT_BG_SIZE <- 0.06
PLOT_BG_ALPHA <- 0.08
PLOT_FG_SIZE <- 0.14
PLOT_FG_ALPHA <- 0.9

CANONICAL_COARSE <- character(0)
CELLTYPE_COLORS <- c()

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

SUBTYPE_COLORS <- c(
  "s1-CAFs" = "#7B3294",
  "s2-CAFs" = "#008837",
  "s3-CAFs" = "#D7191C",
  "s4-CAFs" = "#2C7BB6",
  "s5-CAFs" = "#FDAE61",
  "s6-CAFs" = "#1B9E77",
  "s7-CAFs" = "#E7298A",
  "s8-CAFs" = "#66A61E",
  "Unmapped" = "#BDBDBD"
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

to_bool <- function(x) {
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

get_opt <- function(opts, key, default, cast = NULL) {
  val <- opts[[key]]
  if (is.null(val) || val == "") return(default)
  if (is.null(cast)) return(val)
  cast(val)
}

sort_subtypes <- function(subtypes) {
  key <- suppressWarnings(as.integer(sub("^s([0-9]+)-.*$", "\\1", subtypes)))
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

make_blank_panel_limited <- function(title, xlim, ylim) {
  make_blank_panel(title) +
    geom_blank(data = data.frame(x = xlim, y = ylim), aes(x = x, y = y)) +
    coord_equal(xlim = xlim, ylim = ylim)
}

plot_caf_subtype_panel <- function(seu_df, caf_df, subtype_levels, title, xlim, ylim) {
  caf_df$subtype <- factor(caf_df$subtype, levels = subtype_levels)
  caf_df <- caf_df[order(caf_df$subtype), , drop = FALSE]
  p <- ggplot() +
    geom_point_layer(
      data = seu_df,
      aes(x = x, y = y),
      color = "grey85",
      size = PLOT_BG_SIZE,
      alpha = PLOT_BG_ALPHA
    ) +
    geom_point_layer(
      data = caf_df,
      aes(x = x, y = y, color = subtype),
      size = PLOT_FG_SIZE,
      alpha = PLOT_FG_ALPHA
    ) +
    theme_void(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0, face = "bold", size = 12),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      aspect.ratio = 1
    ) +
    scale_color_manual(values = SUBTYPE_COLORS, drop = FALSE) +
    labs(title = title, color = "CAF subtype")
  apply_limits(p, xlim, ylim)
}

plot_subtype_niche_panel <- function(seu_df, niche_df, subtype_levels, title, xlim, ylim) {
  niche_df$subtype <- factor(niche_df$subtype, levels = subtype_levels)
  niche_df <- niche_df[order(niche_df$subtype), , drop = FALSE]
  p <- ggplot() +
    geom_point_layer(
      data = seu_df,
      aes(x = x, y = y),
      color = "grey85",
      size = PLOT_BG_SIZE,
      alpha = PLOT_BG_ALPHA
    ) +
    geom_point_layer(
      data = niche_df,
      aes(x = x, y = y, color = subtype),
      size = PLOT_FG_SIZE,
      alpha = PLOT_FG_ALPHA
    ) +
    theme_void(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0, face = "bold", size = 12),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      aspect.ratio = 1
    ) +
    scale_color_manual(values = SUBTYPE_COLORS, drop = FALSE) +
    labs(title = title, color = "CAF subtype")
  apply_limits(p, xlim, ylim)
}

plot_neighbor_panel <- function(seu_df, niche_df, title, xlim, ylim) {
  niche_df$cell_type <- factor(niche_df$cell_type, levels = CANONICAL_COARSE)
  niche_df <- niche_df[order(niche_df$cell_type), , drop = FALSE]
  p <- ggplot() +
    geom_point_layer(
      data = seu_df,
      aes(x = x, y = y),
      color = "grey85",
      size = PLOT_BG_SIZE,
      alpha = PLOT_BG_ALPHA
    ) +
    geom_point_layer(
      data = niche_df,
      aes(x = x, y = y, color = cell_type),
      size = PLOT_FG_SIZE,
      alpha = PLOT_FG_ALPHA
    ) +
    theme_void(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0, face = "bold", size = 12),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      aspect.ratio = 1
    ) +
    scale_color_manual(values = CELLTYPE_COLORS, drop = FALSE) +
    labs(title = title, color = "Cell type")
  apply_limits(p, xlim, ylim)
}

plot_pie_panel <- function(fraction_named, title) {
  see <- data.frame(
    cell_type = names(fraction_named),
    fraction = as.numeric(fraction_named),
    stringsAsFactors = FALSE
  )
  see <- see[see$fraction > 0 & !is.na(see$fraction), , drop = FALSE]
  if (nrow(see) == 0) return(make_blank_panel(title))
  see$cell_type <- factor(see$cell_type, levels = CANONICAL_COARSE)
  see <- see[order(see$cell_type), , drop = FALSE]
  ggplot(see, aes(x = 1, y = fraction, fill = cell_type)) +
    geom_col(width = 1, color = "white", linewidth = 0.2) +
    coord_polar(theta = "y") +
    theme_void(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0, face = "bold", size = 11),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      aspect.ratio = 1
    ) +
    scale_fill_manual(values = CELLTYPE_COLORS, drop = FALSE) +
    labs(title = title)
}

compute_neighbors <- function(coords, query_cells, r) {
  query_cells <- intersect(query_cells, rownames(coords))
  if (length(query_cells) == 0) return(list())
  Q <- coords[query_cells, , drop = FALSE]
  nbrs <- dbscan::frNN(coords, eps = r, query = Q, sort = FALSE)
  neighbors <- lapply(nbrs$id, function(idx) rownames(coords)[idx])
  names(neighbors) <- rownames(Q)
  neighbors
}

read_r_choice <- function(sample_dir, r_fixed, use_r_choice) {
  if (!use_r_choice) return(r_fixed)
  path <- file.path(sample_dir, "r_choice.tsv")
  if (!file.exists(path)) return(r_fixed)
  df <- read.delim(path, stringsAsFactors = FALSE)
  as.numeric(df$r_choice[1]) %||% r_fixed
}

choose_example_samples <- function(mapping_df, max_subtypes) {
  df <- mapping_df[mapping_df$excluded_from_global == FALSE & !is.na(mapping_df$global_label), , drop = FALSE]
  if (nrow(df) == 0) stop("No mapped clusters available for example selection.")

  sample_to_subtypes <- split(df$global_label, df$sample_id)
  sample_to_subtypes <- lapply(sample_to_subtypes, function(x) unique(sort_subtypes(unique(x))))
  sizes <- vapply(sample_to_subtypes, length, integer(1))
  full_sample <- names(which.max(sizes))[1]

  subtype_list <- sample_to_subtypes[[full_sample]]
  subtype_list <- head(subtype_list, max_subtypes)

  missing_candidates <- setdiff(names(sample_to_subtypes), full_sample)
  missing_counts <- vapply(missing_candidates, function(s) {
    present <- sample_to_subtypes[[s]]
    sum(!subtype_list %in% present)
  }, integer(1))
  missing_candidates <- missing_candidates[missing_counts > 0]
  if (length(missing_candidates) == 0) {
    stop("Cannot find a sample missing any of the chosen subtypes.")
  }
  missing_sample <- missing_candidates[which.max(missing_counts[missing_candidates])][1]
  list(full_sample = full_sample, missing_sample = missing_sample, subtypes = subtype_list)
}

make_legends <- function(subtypes_to_show) {
  subtype_levels <- c(subtypes_to_show, "Unmapped")
  subtype_colors <- SUBTYPE_COLORS[subtype_levels]
  subtype_colors <- subtype_colors[!is.na(subtype_colors)]

  p1 <- ggplot(
    data.frame(x = 1, y = seq_along(names(subtype_colors)), subtype = names(subtype_colors)),
    aes(x = x, y = y, color = subtype)
  ) +
    geom_point(size = 3.2) +
    theme_void(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.justification = "left",
      legend.box.just = "left",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      legend.text.align = 0
    ) +
    scale_color_manual(values = subtype_colors, drop = FALSE) +
    guides(color = guide_legend(title = NULL, nrow = 1, override.aes = list(alpha = 1)))
  leg_subtype <- cowplot::get_legend(p1)

  celltype_labels <- stats::setNames(gsub("_", " ", CANONICAL_COARSE), CANONICAL_COARSE)
  p2 <- ggplot(
    data.frame(x = 1, y = seq_along(CANONICAL_COARSE), cell_type = CANONICAL_COARSE),
    aes(x = x, y = y, color = cell_type)
  ) +
    geom_point(size = 3.2) +
    theme_void(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.justification = "left",
      legend.box.just = "left",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      legend.text.align = 0
    ) +
    scale_color_manual(values = CELLTYPE_COLORS, breaks = CANONICAL_COARSE, labels = celltype_labels, drop = FALSE) +
    guides(color = guide_legend(title = NULL, nrow = 1, override.aes = list(alpha = 1)))
  leg_celltype <- cowplot::get_legend(p2)

  subtype_block <- cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label("CAFsubtype", x = 0, hjust = 0, fontface = "bold", size = 12),
    cowplot::ggdraw() + cowplot::draw_grob(leg_subtype, x = 0, y = 1, hjust = 0, vjust = 1),
    ncol = 1,
    rel_heights = c(0.25, 1)
  )

  celltype_block <- cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label("celltype", x = 0, hjust = 0, fontface = "bold", size = 12),
    cowplot::ggdraw() + cowplot::draw_grob(leg_celltype, x = 0, y = 1, hjust = 0, vjust = 1),
    ncol = 1,
    rel_heights = c(0.25, 1)
  )

  cowplot::plot_grid(subtype_block, celltype_block, ncol = 2, rel_widths = c(1, 2))
}

argv <- commandArgs(trailingOnly = TRUE)
opts <- parse_args(argv)

defaults <- list(
  input_root = file.path(getwd(), "data", "caf_niche_sectionnorm"),
  global_root = file.path(getwd(), "data", "caf_niche_sectionnorm", "_global_noQT"),
  seurat_root = file.path(getwd(), "data"),
  output = "",
  samples = "",
  subtypes = "",
  example_auto = TRUE,
  k_local = 5L,
  profile = "corrected",
  anno_coarse = "annotation_coarse",
  coarse_schema = "myeloid_refined",
  refine_myeloid = FALSE,
  myeloid_fine_cols = "annotation_fine,cell_type_lvl2,cell_type_lvl1,annotation",
  reduction = "spatial",
  r_fixed = 80,
  max_plot_all_cells = 200000L,
  max_plot_niche_cells = 200000L,
  use_r_choice = TRUE,
  max_subtypes = 6L
)

cfg <- list(
  input_root = get_opt(opts, "input_root", defaults$input_root),
  global_root = get_opt(opts, "global_root", defaults$global_root),
  seurat_root = get_opt(opts, "seurat_root", defaults$seurat_root),
  output = get_opt(opts, "output", defaults$output),
  samples = get_opt(opts, "samples", defaults$samples),
  subtypes = get_opt(opts, "subtypes", defaults$subtypes),
  example_auto = get_opt(opts, "example_auto", defaults$example_auto, to_bool),
  k_local = get_opt(opts, "k_local", defaults$k_local, as.integer),
  profile = get_opt(opts, "profile", defaults$profile),
  anno_coarse = get_opt(opts, "anno_coarse", defaults$anno_coarse),
  coarse_schema = get_opt(opts, "coarse_schema", defaults$coarse_schema),
  refine_myeloid = get_opt(opts, "refine_myeloid", defaults$refine_myeloid, to_bool),
  myeloid_fine_cols = get_opt(opts, "myeloid_fine_cols", defaults$myeloid_fine_cols),
  reduction = get_opt(opts, "reduction", defaults$reduction),
  r_fixed = get_opt(opts, "r_fixed", defaults$r_fixed, as.numeric),
  max_plot_all_cells = get_opt(opts, "max_plot_all_cells", defaults$max_plot_all_cells, as.integer),
  max_plot_niche_cells = get_opt(opts, "max_plot_niche_cells", defaults$max_plot_niche_cells, as.integer),
  use_r_choice = get_opt(opts, "use_r_choice", defaults$use_r_choice, to_bool),
  max_subtypes = get_opt(opts, "max_subtypes", defaults$max_subtypes, as.integer)
)

cfg$coarse_schema <- tolower(cfg$coarse_schema)
if (cfg$refine_myeloid && cfg$coarse_schema == "base") {
  cfg$coarse_schema <- "myeloid_refined"
}
if (!cfg$coarse_schema %in% c("base", "myeloid_refined")) {
  stop(glue("Unknown coarse_schema: {cfg$coarse_schema}. Use base or myeloid_refined."))
}
CANONICAL_COARSE <- get_canonical_coarse(cfg$coarse_schema, order = "grid")
CELLTYPE_COLORS <- get_celltype_colors(cfg$coarse_schema)
cfg$myeloid_fine_cols <- parse_csv_list(cfg$myeloid_fine_cols)
if (length(cfg$myeloid_fine_cols) == 0) {
  cfg$myeloid_fine_cols <- default_myeloid_fine_cols()
}

global_dir <- file.path(cfg$global_root, glue("k{cfg$k_local}"), cfg$profile)
mapping_path <- file.path(global_dir, "global_subtype_mapping.tsv")
if (!file.exists(mapping_path)) stop(glue("Missing mapping: {mapping_path}"))
mapping_df <- read.delim(mapping_path, stringsAsFactors = FALSE)

if (cfg$example_auto) {
  picked <- choose_example_samples(mapping_df, cfg$max_subtypes)
  sample_ids <- c(picked$full_sample, picked$missing_sample)
  subtype_levels <- picked$subtypes
} else if (cfg$samples != "") {
  sample_ids <- trimws(strsplit(cfg$samples, ",")[[1]])
  if (length(sample_ids) == 0) stop("No samples provided.")
  if (cfg$subtypes != "") {
    subtype_levels <- trimws(strsplit(cfg$subtypes, ",")[[1]])
  } else {
    sub_df <- mapping_df[mapping_df$sample_id %in% sample_ids & mapping_df$excluded_from_global == FALSE & !is.na(mapping_df$global_label), , drop = FALSE]
    subtype_levels <- sort_subtypes(unique(sub_df$global_label))
  }
} else {
  stop("Provide --samples or set --example_auto true.")
}

if (length(sample_ids) != 2) {
  stop("This script currently generates a 2-sample example. Provide exactly 2 samples or use --example_auto true.")
}

subtype_levels <- head(sort_subtypes(unique(subtype_levels)), cfg$max_subtypes)
subtype_levels <- subtype_levels[subtype_levels != "" & !is.na(subtype_levels)]
if (length(subtype_levels) == 0) stop("No subtypes available to plot.")

subtype_slots <- subtype_levels

if (cfg$output == "") {
  cfg$output <- file.path(global_dir, "niche_by_sample", "example_niche_by_sample_grid.pdf")
}
dir.create(dirname(cfg$output), recursive = TRUE, showWarnings = FALSE)

build_sample_row <- function(sample_id) {
  sample_dir <- file.path(cfg$input_root, sample_id)
  cluster_path <- file.path(sample_dir, glue("caf_clusters_k{cfg$k_local}.tsv"))
  if (!file.exists(cluster_path)) stop(glue("{sample_id}: missing {basename(cluster_path)}"))

  seurat_path <- file.path(cfg$seurat_root, glue("{sample_id}.qs"))
  if (!file.exists(seurat_path)) stop(glue("{sample_id}: missing Seurat file {basename(seurat_path)}"))

  cluster_df <- read.delim(cluster_path, stringsAsFactors = FALSE)
  map_sub <- mapping_df[mapping_df$sample_id == sample_id, c("local_cluster", "global_label")]
  map_vec <- setNames(map_sub$global_label, map_sub$local_cluster)
  cluster_df$subtype <- map_vec[cluster_df$caf_cluster]
  cluster_df$subtype[is.na(cluster_df$subtype)] <- "Unmapped"

  seu <- qs::qread(seurat_path)
  meta <- seu@meta.data
  if (!cfg$anno_coarse %in% colnames(meta)) stop(glue("{sample_id}: missing {cfg$anno_coarse}"))
  if (!cfg$reduction %in% names(seu@reductions)) stop(glue("{sample_id}: missing reduction {cfg$reduction}"))
  if (!all(c("Center_X", "Center_Y") %in% colnames(meta))) stop(glue("{sample_id}: missing Center_X/Center_Y"))

  if (cfg$refine_myeloid) {
    meta[[cfg$anno_coarse]] <- refine_myeloid_in_coarse(
      meta,
      anno_coarse = cfg$anno_coarse,
      fine_cols = cfg$myeloid_fine_cols
    )
  }
  x <- as.character(meta[[cfg$anno_coarse]])
  unknown <- sort(unique(x[!is.na(x) & !x %in% CANONICAL_COARSE]))
  if (length(unknown) > 0) {
    stop(glue("{sample_id}: unknown {cfg$anno_coarse} values: {paste(unknown, collapse = ', ')}"))
  }

  df_all_full <- data.frame(
    cell_id = rownames(meta),
    x = meta$Center_X * COORD_SCALE,
    y = meta$Center_Y * COORD_SCALE,
    cell_type = x,
    stringsAsFactors = FALSE
  )
  df_bg <- df_all_full
  if (nrow(df_bg) > cfg$max_plot_all_cells) {
    set.seed(RNG_SEED)
    df_bg <- df_bg[sample(seq_len(nrow(df_bg)), cfg$max_plot_all_cells), , drop = FALSE]
  }
  xlim <- pad_limits(df_all_full$x)
  ylim <- pad_limits(df_all_full$y)

  caf_df <- cluster_df[, c("cell_id", "subtype"), drop = FALSE]
  caf_df <- merge(caf_df, df_all_full, by = "cell_id")

  r_use <- read_r_choice(sample_dir, cfg$r_fixed, cfg$use_r_choice)
  coords <- Seurat::Embeddings(seu, reduction = cfg$reduction) * COORD_SCALE

  caf_cells_all <- unique(cluster_df$cell_id)
  neighbors <- compute_neighbors(coords, caf_cells_all, r_use)

  frac_path <- file.path(sample_dir, "niche_celltype_fraction.tsv")
  count_path <- file.path(sample_dir, "niche_celltype_count.tsv")
  frac_mat <- NULL
  count_mat <- NULL
  if (file.exists(frac_path)) {
    frac_mat <- read.delim(frac_path, row.names = 1, check.names = FALSE)
    missing_cols <- setdiff(CANONICAL_COARSE, colnames(frac_mat))
    if (length(missing_cols) > 0) {
      for (mc in missing_cols) frac_mat[[mc]] <- 0
    }
    frac_mat <- frac_mat[, CANONICAL_COARSE, drop = FALSE]
  }
  if (file.exists(count_path)) {
    count_mat <- read.delim(count_path, row.names = 1, check.names = FALSE)
    missing_cols <- setdiff(CANONICAL_COARSE, colnames(count_mat))
    if (length(missing_cols) > 0) {
      for (mc in missing_cols) count_mat[[mc]] <- 0
    }
    count_mat <- count_mat[, CANONICAL_COARSE, drop = FALSE]
  }

  subtype_for_cell <- setNames(cluster_df$subtype, cluster_df$cell_id)
  subtype_for_cell <- subtype_for_cell[unique(names(subtype_for_cell))]

  subtype_niche_cells <- list()
  for (subtype in subtype_levels) {
    caf_cells <- names(subtype_for_cell)[subtype_for_cell == subtype]
    caf_cells <- intersect(caf_cells, names(neighbors))
    if (length(caf_cells) == 0) {
      subtype_niche_cells[[subtype]] <- character(0)
      next
    }
    subtype_niche_cells[[subtype]] <- unique(unlist(neighbors[caf_cells], use.names = FALSE))
  }
  all_niche_cells <- unique(unlist(subtype_niche_cells, use.names = FALSE))
  niche_subtype_df <- df_all_full[df_all_full$cell_id %in% all_niche_cells, c("cell_id", "x", "y"), drop = FALSE]
  niche_subtype_df$subtype <- NA_character_
  for (subtype in subtype_levels) {
    cells <- subtype_niche_cells[[subtype]]
    if (is.null(cells) || length(cells) == 0) next
    idx <- niche_subtype_df$cell_id %in% cells & is.na(niche_subtype_df$subtype)
    niche_subtype_df$subtype[idx] <- subtype
  }
  niche_subtype_df <- niche_subtype_df[!is.na(niche_subtype_df$subtype), , drop = FALSE]

  p_left <- if (nrow(niche_subtype_df) == 0) {
    make_blank_panel_limited(sample_id, xlim, ylim)
  } else {
    plot_subtype_niche_panel(
      seu_df = df_bg,
      niche_df = niche_subtype_df,
      subtype_levels = subtype_levels,
      title = sample_id,
      xlim = xlim,
      ylim = ylim
    )
  }

  neighbor_panels <- list()
  pie_panels <- list()

  for (subtype in subtype_slots) {
    if (is.na(subtype) || subtype == "") {
      neighbor_panels[[paste0("blank_", length(neighbor_panels) + 1)]] <- make_blank_panel("")
      pie_panels[[paste0("blank_", length(pie_panels) + 1)]] <- make_blank_panel("")
      next
    }
    caf_cells <- names(subtype_for_cell)[subtype_for_cell == subtype]
    caf_cells <- intersect(caf_cells, names(neighbors))
    if (length(caf_cells) == 0) {
      neighbor_panels[[subtype]] <- make_blank_panel_limited(subtype, xlim, ylim)
      pie_panels[[subtype]] <- make_blank_panel(subtype)
      next
    }

    niche_cells <- unique(unlist(neighbors[caf_cells], use.names = FALSE))
    niche_cells <- intersect(niche_cells, df_all_full$cell_id)
    if (length(niche_cells) == 0) {
      neighbor_panels[[subtype]] <- make_blank_panel_limited(subtype, xlim, ylim)
    } else {
      niche_df <- df_all_full[df_all_full$cell_id %in% niche_cells, c("x", "y", "cell_type"), drop = FALSE]
      if (nrow(niche_df) > cfg$max_plot_niche_cells) {
        set.seed(RNG_SEED)
        niche_df <- niche_df[sample(seq_len(nrow(niche_df)), cfg$max_plot_niche_cells), , drop = FALSE]
      }
      neighbor_panels[[subtype]] <- plot_neighbor_panel(df_bg, niche_df, subtype, xlim, ylim)
    }

    if (!is.null(count_mat)) {
      caf_cells_in_mat <- intersect(caf_cells, rownames(count_mat))
      if (length(caf_cells_in_mat) == 0) {
        pie_panels[[subtype]] <- make_blank_panel(subtype)
      } else {
        counts <- colSums(count_mat[caf_cells_in_mat, , drop = FALSE], na.rm = TRUE)
        frac <- if (sum(counts) > 0) counts / sum(counts) else counts
        pie_panels[[subtype]] <- plot_pie_panel(frac, subtype)
      }
    } else if (!is.null(frac_mat)) {
      caf_cells_in_mat <- intersect(caf_cells, rownames(frac_mat))
      if (length(caf_cells_in_mat) == 0) {
        pie_panels[[subtype]] <- make_blank_panel(subtype)
      } else {
        frac_mean <- colMeans(frac_mat[caf_cells_in_mat, , drop = FALSE], na.rm = TRUE)
        pie_panels[[subtype]] <- plot_pie_panel(frac_mean, subtype)
      }
    } else {
      pie_panels[[subtype]] <- make_blank_panel(subtype)
    }
  }

  neighbor_grid <- cowplot::plot_grid(
    plotlist = neighbor_panels,
    nrow = 1,
    align = "hv",
    rel_widths = rep(1, length(neighbor_panels))
  )
  pie_grid <- cowplot::plot_grid(
    plotlist = pie_panels,
    ncol = 2,
    align = "hv"
  )

  row <- cowplot::plot_grid(
    p_left,
    neighbor_grid,
    pie_grid,
    nrow = 1,
    rel_widths = c(1.05, length(subtype_slots) * 1.0, 1.05)
  )
  rm(seu)
  gc()
  row
}

rows <- lapply(sample_ids, build_sample_row)
panel <- cowplot::plot_grid(plotlist = rows, ncol = 1, align = "v", rel_heights = rep(1, length(rows)))

left_w <- 1.05
mid_w <- length(subtype_slots) * 1.0
right_w <- 1.05
total_w <- left_w + mid_w + right_w
x_neighbor <- (left_w / total_w) + 0.01
x_pie <- ((left_w + mid_w) / total_w) + 0.01

header <- cowplot::ggdraw() +
  cowplot::draw_label("CAF subtype", x = 0.02, y = 0.98, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  cowplot::draw_label("Neighboring cells", x = x_neighbor, y = 0.98, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  cowplot::draw_label("Neighboring cell composition", x = x_pie, y = 0.98, hjust = 0, vjust = 1, fontface = "bold", size = 14)

legends <- make_legends(subtype_levels)

final <- cowplot::plot_grid(
  header,
  panel,
  legends,
  ncol = 1,
  rel_heights = c(0.12, 1, 0.22)
)

png_width <- as.integer(360 * (length(subtype_slots) + 4))
png_height <- as.integer(360 * length(sample_ids) + 320)

ext <- tolower(tools::file_ext(cfg$output))
if (ext == "png") {
  if (!requireNamespace("ragg", quietly = TRUE)) {
    stop("Need R package 'ragg' to write PNG output. Please install it.")
  }
  ragg::agg_png(cfg$output, width = png_width, height = png_height, res = 300)
  print(final)
  grDevices::dev.off()
  message(glue("Done. Output: {cfg$output}"))
} else if (ext == "pdf") {
  w_in <- png_width / 300
  h_in <- png_height / 300
  grDevices::cairo_pdf(cfg$output, width = w_in, height = h_in)
  print(final)
  grDevices::dev.off()
  message(glue("Done. Output: {cfg$output}"))
} else {
  stop("output must end with .png or .pdf")
}
