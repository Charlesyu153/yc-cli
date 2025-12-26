#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(ggplot2)
  library(glue)
  library(dbscan)
  library(ggsci)
})

`%||%` <- function(x, y) if (is.null(x) || x == "") y else x

COORD_SCALE <- 50
RNG_SEED <- 1024

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
  input_root = "CAF niche output root (per-sample dirs)",
  global_root = "Global alignment root (contains global_subtype_mapping.tsv)",
  seurat_root = "Directory containing Seurat .qs files",
  output_dir = "Directory for per-sample PDF outputs",
  samples = "Comma-separated sample IDs (blank = from mapping)",
  subtypes = "Comma-separated subtypes (blank = s1-CAFs,s2-CAFs if present)",
  k_local = "Local CAF clustering k",
  profile = "Global alignment profile: corrected or raw",
  anno_coarse = "Coarse annotation column name",
  coarse_schema = "Coarse schema: base or myeloid_refined",
  refine_myeloid = "Refine Myeloid into Macrophage/Neutrophil/Mast (true/false)",
  myeloid_fine_cols = "Comma-separated meta columns used to refine Myeloid",
  reduction = "Seurat reduction name for spatial coords",
  r_fixed = "Fallback r when r_choice.tsv missing",
  max_plot_cells = "Max niche cells to plot per subtype",
  max_plot_all_cells = "Max background cells to plot",
  use_r_choice = "Use r_choice.tsv when available (true/false)"
)

print_help <- function(defaults) {
  cat("CAF global subtype spatial niche plots - parameters\n")
  for (key in names(PARAM_HELP)) {
    def <- defaults[[key]]
    def_str <- if (is.na(def) || is.null(def)) "NA" else as.character(def)
    cat(glue("  --{key} {PARAM_HELP[[key]]} (default: {def_str})\n"))
  }
}

enforce_canonical_levels <- function(seu, anno_col, sample_id) {
  x <- as.character(seu@meta.data[[anno_col]])
  unknown <- sort(unique(x[!is.na(x) & !x %in% CANONICAL_COARSE]))
  if (length(unknown) > 0) {
    stop(glue("{sample_id}: unknown {anno_col} values: {paste(unknown, collapse = ', ')}"))
  }
  seu@meta.data[[anno_col]] <- factor(x, levels = CANONICAL_COARSE)
  seu
}

downsample_niche_proportional <- function(df, max_total) {
  n_total <- nrow(df)
  if (n_total <= max_total) return(df)

  counts <- table(df$niche)
  types <- names(counts)
  n_types <- length(types)

  if (max_total < n_types) {
    set.seed(RNG_SEED)
    return(df[sample(seq_len(n_total), max_total), , drop = FALSE])
  }

  props <- counts / sum(counts)
  raw <- props * max_total
  target <- floor(raw)
  target[target == 0] <- 1

  diff <- max_total - sum(target)
  if (diff > 0) {
    frac <- raw - floor(raw)
    order_idx <- order(frac, decreasing = TRUE)
    for (idx in order_idx) {
      if (diff == 0) break
      if (target[idx] < counts[idx]) {
        target[idx] <- target[idx] + 1
        diff <- diff - 1
      }
    }
  } else if (diff < 0) {
    order_idx <- order(target, decreasing = TRUE)
    for (idx in order_idx) {
      if (diff == 0) break
      if (target[idx] > 1) {
        target[idx] <- target[idx] - 1
        diff <- diff + 1
      }
    }
  }
  names(target) <- types

  set.seed(RNG_SEED)
  do.call(
    rbind,
    lapply(types, function(tp) {
      df_tp <- df[df$niche == tp, , drop = FALSE]
      n_keep <- min(target[tp], nrow(df_tp))
      if (n_keep <= 0) return(NULL)
      df_tp[sample(seq_len(nrow(df_tp)), n_keep), , drop = FALSE]
    })
  )
}

plot_niche_panels <- function(seu,
                              caf_cells,
                              niche_cells,
                              annotation_col,
                              max_all_cells = 200000,
                              max_niche_cells = 200000,
                              title = NULL,
                              level_order = CANONICAL_COARSE) {
  meta <- seu@meta.data
  if (!all(c("Center_X", "Center_Y") %in% colnames(meta))) {
    warning("Center_X/Center_Y not found; skip spatial plot.")
    return(NULL)
  }

  df_all <- data.frame(
    cell_id = rownames(meta),
    x = meta$Center_X * COORD_SCALE,
    y = meta$Center_Y * COORD_SCALE,
    niche = if (annotation_col %in% colnames(meta)) {
      as.character(meta[[annotation_col]])
    } else {
      NA_character_
    }
  )
  if (nrow(df_all) == 0) {
    warning("No cells available for plotting; skip.")
    return(NULL)
  }

  df_niche <- df_all[df_all$cell_id %in% niche_cells, , drop = FALSE]
  if (nrow(df_niche) == 0) {
    warning("No niche cells to plot; skip.")
    return(NULL)
  }
  if (nrow(df_all) > max_all_cells) {
    set.seed(RNG_SEED)
    df_all <- df_all[sample(seq_len(nrow(df_all)), max_all_cells), , drop = FALSE]
  }
  if (nrow(df_niche) > max_niche_cells) {
    df_niche <- downsample_niche_proportional(df_niche, max_niche_cells)
  }

  df_niche$niche <- factor(df_niche$niche, levels = level_order)
  df_niche <- df_niche[order(df_niche$niche), , drop = FALSE]

  caf_df <- df_all[df_all$cell_id %in% caf_cells, , drop = FALSE]

  p_all <- ggplot() +
    geom_point(
      data = df_all,
      aes(x = x, y = y),
      color = "grey85",
      size = 0.01,
      alpha = 0.25
    ) +
    geom_point(
      data = caf_df,
      aes(x = x, y = y),
      color = "#E74C3C",
      size = 0.01,
      alpha = 0.9
    ) +
    coord_equal() +
    theme_void(base_size = 11) +
    labs(title = paste0(title, " - CAF niche anchors"))

  p_niche <- ggplot() +
    geom_point(
      data = df_all,
      aes(x = x, y = y),
      color = "grey85",
      size = 0.01,
      alpha = 0.35
    ) +
    geom_point(
      data = df_niche,
      aes(x = x, y = y, color = niche),
      shape = 16,
      size = 0.08,
      alpha = 0.6
    ) +
    coord_equal() +
    theme_void(base_size = 11) +
    theme(
      legend.position = "right",
      legend.key.size = grid::unit(0.6, "cm"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11)
    ) +
    labs(title = paste0(title, " - niche cells on slice"), color = annotation_col) +
    guides(color = guide_legend(override.aes = list(size = 3.2, alpha = 1, shape = 16)))

  p_niche <- p_niche + ggsci::scale_color_d3(
    "category20",
    limits = level_order,
    breaks = level_order,
    drop = FALSE
  )

  if (requireNamespace("cowplot", quietly = TRUE)) {
    return(cowplot::plot_grid(p_all, p_niche, ncol = 2))
  }

  list(p_all = p_all, p_niche = p_niche)
}

plot_coarse_slice <- function(seu,
                              annotation_col,
                              max_all_cells = 200000,
                              title = NULL,
                              level_order = CANONICAL_COARSE) {
  meta <- seu@meta.data
  if (!all(c("Center_X", "Center_Y") %in% colnames(meta))) {
    warning("Center_X/Center_Y not found; skip coarse slice plot.")
    return(NULL)
  }

  df_all <- data.frame(
    cell_id = rownames(meta),
    x = meta$Center_X * COORD_SCALE,
    y = meta$Center_Y * COORD_SCALE,
    niche = if (annotation_col %in% colnames(meta)) {
      as.character(meta[[annotation_col]])
    } else {
      NA_character_
    }
  )
  df_all <- df_all[!is.na(df_all$niche), , drop = FALSE]
  if (nrow(df_all) == 0) {
    warning("No annotated cells available for plotting; skip.")
    return(NULL)
  }
  if (nrow(df_all) > max_all_cells) {
    df_all <- downsample_niche_proportional(df_all, max_all_cells)
  }
  df_all$niche <- factor(df_all$niche, levels = level_order)
  df_all <- df_all[order(df_all$niche), , drop = FALSE]

  p <- ggplot(df_all, aes(x = x, y = y, color = niche)) +
    geom_point(size = 0.08, alpha = 0.6) +
    coord_equal() +
    theme_void(base_size = 11) +
    theme(
      legend.position = "right",
      legend.key.size = grid::unit(0.6, "cm"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11)
    ) +
    labs(title = title, color = annotation_col) +
    guides(color = guide_legend(override.aes = list(size = 3.2, alpha = 1, shape = 16)))

  p <- p + ggsci::scale_color_d3(
    "category20",
    limits = level_order,
    breaks = level_order,
    drop = FALSE
  )
  p
}

compute_neighbors <- function(coords, query_cells, r) {
  if (!all(query_cells %in% rownames(coords))) {
    query_cells <- intersect(query_cells, rownames(coords))
  }
  if (length(query_cells) == 0) return(list())
  Q <- coords[query_cells, , drop = FALSE]
  nbrs <- dbscan::frNN(coords, eps = r, query = Q, sort = FALSE)
  neighbors <- lapply(nbrs$id, function(idx) rownames(coords)[idx])
  names(neighbors) <- rownames(Q)
  neighbors
}

argv <- commandArgs(trailingOnly = TRUE)
opts <- parse_args(argv)

defaults <- list(
  input_root = file.path(getwd(), "data", "caf_niche_sectionnorm"),
  global_root = file.path(getwd(), "data", "caf_niche_sectionnorm", "_global_noQT"),
  seurat_root = file.path(getwd(), "data"),
  output_dir = "",
  samples = "",
  subtypes = "",
  k_local = 5L,
  profile = "corrected",
  anno_coarse = "annotation_coarse",
  coarse_schema = "myeloid_refined",
  refine_myeloid = FALSE,
  myeloid_fine_cols = "annotation_fine,cell_type_lvl2,cell_type_lvl1,annotation",
  reduction = "spatial",
  r_fixed = 80,
  max_plot_cells = 200000L,
  max_plot_all_cells = 200000L,
  use_r_choice = TRUE
)

if (!is.null(opts$help) || !is.null(opts$h)) {
  print_help(defaults)
  quit(status = 0)
}

cfg <- list(
  input_root = get_opt(opts, "input_root", defaults$input_root),
  global_root = get_opt(opts, "global_root", defaults$global_root),
  seurat_root = get_opt(opts, "seurat_root", defaults$seurat_root),
  output_dir = get_opt(opts, "output_dir", defaults$output_dir),
  samples = get_opt(opts, "samples", defaults$samples),
  subtypes = get_opt(opts, "subtypes", defaults$subtypes),
  k_local = get_opt(opts, "k_local", defaults$k_local, as.integer),
  profile = get_opt(opts, "profile", defaults$profile),
  anno_coarse = get_opt(opts, "anno_coarse", defaults$anno_coarse),
  coarse_schema = get_opt(opts, "coarse_schema", defaults$coarse_schema),
  refine_myeloid = get_opt(opts, "refine_myeloid", defaults$refine_myeloid, to_bool),
  myeloid_fine_cols = get_opt(opts, "myeloid_fine_cols", defaults$myeloid_fine_cols),
  reduction = get_opt(opts, "reduction", defaults$reduction),
  r_fixed = get_opt(opts, "r_fixed", defaults$r_fixed, as.numeric),
  max_plot_cells = get_opt(opts, "max_plot_cells", defaults$max_plot_cells, as.integer),
  max_plot_all_cells = get_opt(opts, "max_plot_all_cells", defaults$max_plot_all_cells, as.integer),
  use_r_choice = get_opt(opts, "use_r_choice", defaults$use_r_choice, to_bool)
)

cfg$coarse_schema <- tolower(cfg$coarse_schema)
if (cfg$refine_myeloid && cfg$coarse_schema == "base") {
  cfg$coarse_schema <- "myeloid_refined"
}
if (!cfg$coarse_schema %in% c("base", "myeloid_refined")) {
  stop(glue("Unknown coarse_schema: {cfg$coarse_schema}. Use base or myeloid_refined."))
}
CANONICAL_COARSE <- get_canonical_coarse(cfg$coarse_schema, order = "grid")
cfg$myeloid_fine_cols <- parse_csv_list(cfg$myeloid_fine_cols)
if (length(cfg$myeloid_fine_cols) == 0) {
  cfg$myeloid_fine_cols <- default_myeloid_fine_cols()
}

global_dir <- file.path(cfg$global_root, glue("k{cfg$k_local}"), cfg$profile)
mapping_path <- file.path(global_dir, "global_subtype_mapping.tsv")
if (!file.exists(mapping_path)) {
  stop(glue("global_subtype_mapping.tsv not found: {mapping_path}"))
}
mapping_df <- read.delim(mapping_path, stringsAsFactors = FALSE)
mapping_df <- mapping_df[mapping_df$excluded_from_global == FALSE, , drop = FALSE]

if (cfg$samples != "") {
  sample_ids <- trimws(strsplit(cfg$samples, ",")[[1]])
} else {
  sample_ids <- unique(mapping_df$sample_id)
}
sample_ids <- sample_ids[!is.na(sample_ids)]

if (cfg$subtypes != "") {
  subtype_list <- trimws(strsplit(cfg$subtypes, ",")[[1]])
} else {
  subtype_list <- c("s1-CAFs", "s2-CAFs")
  present <- intersect(subtype_list, unique(mapping_df$global_label))
  if (length(present) > 0) {
    subtype_list <- present
  } else {
    subtype_list <- unique(mapping_df$global_label)
  }
}

if (cfg$output_dir == "") {
  cfg$output_dir <- file.path(global_dir, "niche_by_sample")
}
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

for (sample_id in sample_ids) {
  sample_dir <- file.path(cfg$input_root, sample_id)
  if (!dir.exists(sample_dir)) {
    warning(glue("{sample_id}: output dir missing; skip."))
    next
  }

  cluster_path <- file.path(sample_dir, glue("caf_clusters_k{cfg$k_local}.tsv"))
  if (!file.exists(cluster_path)) {
    warning(glue("{sample_id}: {basename(cluster_path)} missing; skip."))
    next
  }

  mapping_sub <- mapping_df[mapping_df$sample_id == sample_id, c("local_cluster", "global_label")]
  if (nrow(mapping_sub) == 0) {
    warning(glue("{sample_id}: no global mapping rows; skip."))
    next
  }
  map_vec <- setNames(mapping_sub$global_label, mapping_sub$local_cluster)

  cluster_df <- read.delim(cluster_path, stringsAsFactors = FALSE)
  cluster_df$global_label <- map_vec[cluster_df$caf_cluster]
  cluster_df <- cluster_df[!is.na(cluster_df$global_label), , drop = FALSE]
  if (nrow(cluster_df) == 0) {
    warning(glue("{sample_id}: no CAF cells with global subtype; skip."))
    next
  }

  seurat_path <- file.path(cfg$seurat_root, glue("{sample_id}.qs"))
  if (!file.exists(seurat_path)) {
    warning(glue("{sample_id}: missing Seurat file; skip."))
    next
  }

  seu <- qs::qread(seurat_path)
  if (!cfg$anno_coarse %in% colnames(seu@meta.data)) {
    warning(glue("{sample_id}: missing {cfg$anno_coarse}; skip."))
    next
  }
  if (!cfg$reduction %in% names(seu@reductions)) {
    warning(glue("{sample_id}: reduction {cfg$reduction} missing; skip."))
    next
  }
  if (cfg$refine_myeloid) {
    seu@meta.data[[cfg$anno_coarse]] <- refine_myeloid_in_coarse(
      seu@meta.data,
      anno_coarse = cfg$anno_coarse,
      fine_cols = cfg$myeloid_fine_cols
    )
  }
  seu <- enforce_canonical_levels(seu, cfg$anno_coarse, sample_id)

  coords_scaled <- Seurat::Embeddings(seu, reduction = cfg$reduction) * COORD_SCALE
  seu[[cfg$reduction]]@cell.embeddings <- coords_scaled
  coords <- Seurat::Embeddings(seu, reduction = cfg$reduction)

  coarse_out <- file.path(cfg$output_dir, glue("{sample_id}_coarse_annotation.pdf"))
  coarse_plot <- plot_coarse_slice(
    seu,
    annotation_col = cfg$anno_coarse,
    max_all_cells = cfg$max_plot_all_cells,
    title = glue("{sample_id} coarse annotation"),
    level_order = CANONICAL_COARSE
  )
  if (!is.null(coarse_plot)) {
    ggsave(
      coarse_out,
      coarse_plot,
      width = 7,
      height = 7,
      device = grDevices::cairo_pdf
    )
  }

  r_choice <- cfg$r_fixed
  if (cfg$use_r_choice) {
    choice_path <- file.path(sample_dir, "r_choice.tsv")
    if (file.exists(choice_path)) {
      choice_df <- read.delim(choice_path, stringsAsFactors = FALSE)
      r_choice <- choice_df$r_choice[1]
    }
  }

  caf_cells_all <- cluster_df$cell_id
  message(glue("{sample_id}: compute neighbors (r={r_choice}, CAF={length(caf_cells_all)})"))
  neighbors <- compute_neighbors(coords, caf_cells_all, r_choice)

  out_file <- file.path(cfg$output_dir, glue("{sample_id}_global_subtype_niche.pdf"))
  grDevices::cairo_pdf(out_file, width = 14, height = 7)

  for (subtype in subtype_list) {
    subtype_cells <- cluster_df$cell_id[cluster_df$global_label == subtype]
    subtype_cells <- intersect(subtype_cells, names(neighbors))
    if (length(subtype_cells) == 0) {
      warning(glue("{sample_id}: no CAF cells for {subtype}; skip."))
      next
    }
    niche_cells <- unique(unlist(neighbors[subtype_cells], use.names = FALSE))
    if (length(niche_cells) == 0) {
      warning(glue("{sample_id}: no niche cells for {subtype}; skip."))
      next
    }
    title <- glue("{sample_id} {subtype} niche (coarse)")
    p <- plot_niche_panels(
      seu,
      caf_cells = subtype_cells,
      niche_cells = niche_cells,
      annotation_col = cfg$anno_coarse,
      max_all_cells = cfg$max_plot_all_cells,
      max_niche_cells = cfg$max_plot_cells,
      title = title,
      level_order = CANONICAL_COARSE
    )
    if (is.null(p)) next
    if (is.list(p) && !inherits(p, "gg")) {
      print(p$p_all)
      print(p$p_niche)
    } else {
      print(p)
    }
  }

  grDevices::dev.off()
  rm(seu)
  gc()
}

message(glue("Done. Output: {cfg$output_dir}"))
