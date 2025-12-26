#!/usr/bin/env Rscript

# CAF niche Step1 (v2): per-sample auto-K with stability-first selection.
suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(glue)
  library(dbscan)
  library(RcppML)
})

`%||%` <- function(x, y) if (is.null(x) || x == "") y else x

if (requireNamespace("pbapply", quietly = TRUE)) {
  pbapply::pboptions(type = "timer")
}

# NOTE: QC on current `.qs` indicates spatial coords are already in micron-like units:
# with coord_scale=1 and r=80, CAF neighborhoods contain ~O(10^2) cells;
# with coord_scale=50 and r=80, neighborhoods collapse to self-only (count=1).
COORD_SCALE <- 1
R_FIXED <- 80
RNG_SEED <- 1024

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
  if (grepl("P$", sample_id)) return("P")
  if (grepl("[LR]$", sample_id)) return("M")
  NA_character_
}

align_canonical_columns <- function(mat, canonical) {
  if (is.null(colnames(mat))) stop("Matrix missing column names.")
  missing <- setdiff(canonical, colnames(mat))
  if (length(missing) > 0) {
    zeros <- matrix(0, nrow = nrow(mat), ncol = length(missing))
    colnames(zeros) <- missing
    mat <- cbind(mat, zeros)
  }
  extras <- setdiff(colnames(mat), canonical)
  if (length(extras) > 0) {
    mat <- mat[, setdiff(colnames(mat), extras), drop = FALSE]
  }
  mat[, canonical, drop = FALSE]
}

compute_section_prop <- function(meta, anno_col, canonical) {
  x <- factor(as.character(meta[[anno_col]]), levels = canonical)
  tab <- table(x)
  total <- sum(tab)
  prop <- if (total == 0) rep(NA_real_, length(canonical)) else as.numeric(tab) / total
  names(prop) <- canonical
  prop
}

cluster_average <- function(mat, cluster_vec) {
  clusters <- unique(as.character(cluster_vec))
  clusters <- clusters[order(clusters)]
  out <- matrix(0, nrow = length(clusters), ncol = ncol(mat))
  rownames(out) <- clusters
  colnames(out) <- colnames(mat)
  for (cl in clusters) {
    cells <- names(cluster_vec)[as.character(cluster_vec) == cl]
    if (length(cells) == 0) next
    out[cl, ] <- colMeans(mat[cells, , drop = FALSE], na.rm = TRUE)
  }
  out
}

find_clusters_with_cap <- function(seu_obj,
                                  cluster_name,
                                  resolution,
                                  max_clusters,
                                  min_resolution,
                                  step,
                                  seed,
                                  log_prefix = "") {
  res <- resolution
  repeat {
    seu_obj <- FindClusters(
      seu_obj,
      resolution = res,
      algorithm = 4,
      cluster.name = cluster_name,
      random.seed = max(seed, 1),
      verbose = FALSE
    )
    clusters <- seu_obj[[cluster_name]][, 1]
    n_clusters <- length(unique(clusters))
    if (n_clusters <= max_clusters || res <= min_resolution + 1e-12) {
      return(list(seu = seu_obj, clusters = clusters, n_clusters = n_clusters, resolution = res))
    }
    new_res <- res - step
    if (new_res < min_resolution) new_res <- min_resolution
    if (log_prefix != "") {
      message(glue("{log_prefix}: clusters={n_clusters} > {max_clusters}; resolution {res} -> {new_res}"))
    }
    if (abs(new_res - res) < 1e-12) {
      return(list(seu = seu_obj, clusters = clusters, n_clusters = n_clusters, resolution = res))
    }
    res <- new_res
  }
}

calc_fit_within_dispersion <- function(mat, clusters) {
  clusters_chr <- as.character(clusters)
  if (!is.null(names(clusters))) {
    names(clusters_chr) <- names(clusters)
  }
  ok <- !is.na(clusters_chr) & clusters_chr != ""
  clusters_chr <- clusters_chr[ok]
  if (length(clusters_chr) < 3) return(NA_real_)

  common_cells <- intersect(names(clusters_chr), rownames(mat))
  if (length(common_cells) < 3) return(NA_real_)
  clusters_chr <- clusters_chr[common_cells]
  mat <- mat[common_cells, , drop = FALSE]

  cents <- cluster_average(mat, clusters_chr)
  disp <- 0
  n <- 0
  for (cl in rownames(cents)) {
    cells <- names(clusters_chr)[clusters_chr == cl]
    if (length(cells) == 0) next
    d <- mat[cells, , drop = FALSE] - matrix(cents[cl, ], nrow = length(cells), ncol = ncol(mat), byrow = TRUE)
    disp <- disp + sum(rowSums(d * d))
    n <- n + length(cells)
  }
  if (n == 0) return(NA_real_)
  disp / n
}

profile_similarity_one_pair <- function(p1, p2) {
  p1 <- as.matrix(p1)
  p2 <- as.matrix(p2)
  if (nrow(p1) == 0 || nrow(p2) == 0) return(NA_real_)
  cor_mat <- suppressWarnings(stats::cor(t(p1), t(p2), method = "pearson", use = "pairwise.complete.obs"))
  cor_mat[!is.finite(cor_mat)] <- -1
  n1 <- nrow(cor_mat)
  n2 <- ncol(cor_mat)
  n <- max(n1, n2)
  pad <- matrix(-1, nrow = n, ncol = n)
  pad[seq_len(n1), seq_len(n2)] <- cor_mat
  assign <- clue::solve_LSAP(as.matrix(1 - pad))
  idx1 <- seq_len(n1)
  idx2 <- as.integer(assign)[seq_len(n1)]
  keep <- idx2 >= 1 & idx2 <= n2
  if (!any(keep)) return(NA_real_)
  mean(pad[cbind(idx1[keep], idx2[keep])])
}

compute_stability_metrics <- function(run_profiles, run_assignments) {
  runs <- sort(unique(names(run_profiles)))
  if (length(runs) < 2) {
    return(list(stab_profile = NA_real_, stab_assign = NA_real_, stab_combined = NA_real_))
  }
  pair_idx <- utils::combn(runs, 2, simplify = FALSE)
  prof_scores <- numeric(0)
  ari_scores <- numeric(0)
  for (pair in pair_idx) {
    r1 <- pair[1]
    r2 <- pair[2]
    sim <- profile_similarity_one_pair(run_profiles[[r1]], run_profiles[[r2]])
    prof_scores <- c(prof_scores, sim)

    a1 <- run_assignments[[r1]]
    a2 <- run_assignments[[r2]]
    common <- intersect(names(a1), names(a2))
    if (length(common) < 3) {
      ari_scores <- c(ari_scores, NA_real_)
    } else {
      ari <- mclust::adjustedRandIndex(as.character(a1[common]), as.character(a2[common]))
      ari_scores <- c(ari_scores, ari)
    }
  }
  prof_mean <- mean(prof_scores, na.rm = TRUE)
  ari_mean <- mean(ari_scores, na.rm = TRUE)
  stab_profile <- max(0, min(1, (prof_mean + 1) / 2))
  stab_assign <- max(0, min(1, (ari_mean + 1) / 2))
  stab_combined <- 0.7 * stab_profile + 0.3 * stab_assign
  list(stab_profile = stab_profile, stab_assign = stab_assign, stab_combined = stab_combined)
}

select_k_stable <- function(k_df, delta = 0.08, fit_min = 0.65) {
  stab_max <- max(k_df$stab_combined, na.rm = TRUE)
  k_stable <- k_df$k[k_df$stab_combined >= (stab_max - delta)]
  k_good <- k_df$k[k_df$k %in% k_stable & k_df$fit_score >= fit_min]
  if (length(k_good) > 0) return(min(k_good))
  k_top <- k_df[order(k_df$stab_combined, decreasing = TRUE), , drop = FALSE]
  k_top <- head(k_top, 3)
  k_top$fit_score[is.na(k_top$fit_score)] <- -Inf
  k_top$k[which.max(k_top$fit_score)]
}

build_consensus_profile <- function(run_profiles, ref_run = NULL) {
  runs <- sort(unique(names(run_profiles)))
  if (length(runs) == 0) stop("No run profiles provided.")
  if (is.null(ref_run) || !ref_run %in% runs) ref_run <- runs[1]
  ref <- run_profiles[[ref_run]]
  ref <- as.matrix(ref)
  if (nrow(ref) == 0) stop("Reference run has no clusters.")
  ref_ids <- rownames(ref)
  ref_ids <- ref_ids[order(ref_ids)]
  ref <- ref[ref_ids, , drop = FALSE]

  aligned <- list()
  aligned[[ref_run]] <- ref

  maps <- list()
  maps[[ref_run]] <- stats::setNames(ref_ids, ref_ids)

  for (r in setdiff(runs, ref_run)) {
    p <- as.matrix(run_profiles[[r]])
    if (nrow(p) == 0) next
    cor_mat <- suppressWarnings(stats::cor(t(ref), t(p), method = "pearson", use = "pairwise.complete.obs"))
    cor_mat[!is.finite(cor_mat)] <- -1
    n1 <- nrow(cor_mat)
    n2 <- ncol(cor_mat)
    n <- max(n1, n2)
    pad <- matrix(-1, nrow = n, ncol = n)
    pad[seq_len(n1), seq_len(n2)] <- cor_mat
    assign <- clue::solve_LSAP(as.matrix(1 - pad))
    idx2 <- as.integer(assign)[seq_len(n1)]
    map <- rep(NA_character_, nrow(p))
    names(map) <- rownames(p)
    for (i in seq_len(n1)) {
      j <- idx2[i]
      if (j >= 1 && j <= n2) {
        map[colnames(cor_mat)[j]] <- rownames(cor_mat)[i]
      }
    }
    maps[[r]] <- map
  }

  consensus <- matrix(NA_real_, nrow = nrow(ref), ncol = ncol(ref))
  rownames(consensus) <- ref_ids
  colnames(consensus) <- colnames(ref)
  for (ref_cl in ref_ids) {
    mats <- list()
    for (r in runs) {
      p <- as.matrix(run_profiles[[r]])
      if (nrow(p) == 0) next
      map <- maps[[r]]
      if (is.null(map)) next
      src <- names(map)[map == ref_cl]
      if (length(src) != 1) next
      mats[[r]] <- p[src, , drop = FALSE]
    }
    if (length(mats) == 0) next
    stack <- do.call(rbind, mats)
    consensus[ref_cl, ] <- colMeans(stack, na.rm = TRUE)
  }
  list(profile = consensus, maps = maps)
}

build_consensus_assignment <- function(run_assignments, run_to_ref_maps) {
  runs <- sort(unique(names(run_assignments)))
  all_cells <- unique(unlist(lapply(run_assignments, names), use.names = FALSE))
  out <- rep(NA_character_, length(all_cells))
  names(out) <- all_cells
  for (cell in all_cells) {
    votes <- character(0)
    for (r in runs) {
      a <- run_assignments[[r]]
      if (is.null(a) || !cell %in% names(a)) next
      cl <- as.character(a[[cell]])
      map <- run_to_ref_maps[[r]]
      if (is.null(map) || is.na(map[[cl]]) || map[[cl]] == "") next
      votes <- c(votes, map[[cl]])
    }
    if (length(votes) == 0) next
    tab <- sort(table(votes), decreasing = TRUE)
    out[cell] <- names(tab)[1]
  }
  out
}

PARAM_HELP <- list(
  input = "Input directory with Seurat .qs files",
  output_root = "Output root directory for Step1 auto-K results",
  samples = "Comma-separated sample IDs (blank = all, excluding QT and conversion_summary)",
  pattern = "Regex pattern for .qs files",
  anno_coarse = "Coarse annotation column name",
  caf_label = "Label in anno_coarse defining CAF cells",
  reduction = "Seurat reduction with spatial coords (default: spatial)",
  r_fixed = "Neighborhood radius in microns (default: 80)",
  coord_scale = "Pixel-to-micron scale (default: 50)",
  k_list = "Comma-separated candidate NMF k values",
  n_runs = "Number of repeated runs per k",
  subsample_frac = "CAF subsample fraction per run (default: 0.8)",
  max_caf_run_cells = "Optional cap on CAF cells per run (blank = no cap)",
  profile = "Matrix for clustering/profiles: corrected or raw",
  max_caf_clusters = "Upper bound on CAF cluster count per run",
  cluster_resolution = "Leiden resolution (will be decreased if clusters exceed cap)",
  cluster_resolution_step = "Resolution decrement step",
  min_cluster_resolution = "Minimum resolution",
  stab_delta = "Delta for stable k set",
  fit_min = "Minimum fit_score to accept k",
  seed = "Base random seed"
)

print_help <- function(defaults) {
  cat("CAF niche Step1 v2: auto-K (stability-first)\n")
  for (key in names(PARAM_HELP)) {
    def <- defaults[[key]]
    def_str <- if (is.na(def) || is.null(def)) "NA" else as.character(def)
    cat(glue("  --{key} {PARAM_HELP[[key]]} (default: {def_str})\n"))
  }
}

argv <- commandArgs(trailingOnly = TRUE)
opts <- parse_args(argv)

root_dir <- get_root_dir()
source(file.path(root_dir, "utils", "R", "niche.R"))
source(file.path(root_dir, "utils", "R", "coarse_schema.R"))

defaults <- list(
  input = file.path(root_dir, "data"),
  output_root = file.path(root_dir, "data", "caf_niche_v2_step1"),
  samples = "",
  pattern = "\\.qs$",
  anno_coarse = "annotation_coarse",
  caf_label = "CAF",
  reduction = "spatial",
  r_fixed = R_FIXED,
  coord_scale = COORD_SCALE,
  k_list = "3,4,5,6,7,8",
  n_runs = 5L,
  subsample_frac = 0.8,
  max_caf_run_cells = NA_integer_,
  profile = "corrected",
  max_caf_clusters = 5L,
  cluster_resolution = 0.1,
  cluster_resolution_step = 0.02,
  min_cluster_resolution = 0.02,
  stab_delta = 0.08,
  fit_min = 0.65,
  seed = RNG_SEED
)

if (!is.null(opts$help) || !is.null(opts$h)) {
  print_help(defaults)
  quit(status = 0)
}

cfg <- list(
  input = get_opt(opts, "input", defaults$input),
  output_root = get_opt(opts, "output_root", defaults$output_root),
  samples = get_opt(opts, "samples", defaults$samples),
  pattern = get_opt(opts, "pattern", defaults$pattern),
  anno_coarse = get_opt(opts, "anno_coarse", defaults$anno_coarse),
  caf_label = get_opt(opts, "caf_label", defaults$caf_label),
  reduction = get_opt(opts, "reduction", defaults$reduction),
  r_fixed = get_opt(opts, "r_fixed", defaults$r_fixed, as.numeric),
  coord_scale = get_opt(opts, "coord_scale", defaults$coord_scale, as.numeric),
  k_list = get_opt(opts, "k_list", defaults$k_list),
  n_runs = get_opt(opts, "n_runs", defaults$n_runs, as.integer),
  subsample_frac = get_opt(opts, "subsample_frac", defaults$subsample_frac, as.numeric),
  max_caf_run_cells = get_opt(opts, "max_caf_run_cells", defaults$max_caf_run_cells, as.integer),
  profile = tolower(get_opt(opts, "profile", defaults$profile)),
  max_caf_clusters = get_opt(opts, "max_caf_clusters", defaults$max_caf_clusters, as.integer),
  cluster_resolution = get_opt(opts, "cluster_resolution", defaults$cluster_resolution, as.numeric),
  cluster_resolution_step = get_opt(opts, "cluster_resolution_step", defaults$cluster_resolution_step, as.numeric),
  min_cluster_resolution = get_opt(opts, "min_cluster_resolution", defaults$min_cluster_resolution, as.numeric),
  stab_delta = get_opt(opts, "stab_delta", defaults$stab_delta, as.numeric),
  fit_min = get_opt(opts, "fit_min", defaults$fit_min, as.numeric),
  seed = get_opt(opts, "seed", defaults$seed, as.integer)
)

cfg$k_list <- as.integer(trimws(strsplit(cfg$k_list, ",")[[1]]))
if (anyNA(cfg$k_list) || length(cfg$k_list) == 0) stop("Invalid k_list")
if (cfg$n_runs < 1) stop("n_runs must be >= 1")
if (cfg$subsample_frac <= 0 || cfg$subsample_frac > 1) stop("subsample_frac must be in (0,1]")
if (!cfg$profile %in% c("corrected", "raw")) stop("profile must be corrected or raw")

dir.create(cfg$output_root, recursive = TRUE, showWarnings = FALSE)

input_files <- list.files(cfg$input, pattern = cfg$pattern, full.names = TRUE)
input_files <- input_files[basename(input_files) != "conversion_summary.qs"]
if (cfg$samples != "") {
  sample_list <- trimws(strsplit(cfg$samples, ",")[[1]])
  input_files <- input_files[basename(input_files) %in% paste0(sample_list, ".qs")]
} else {
  # Exclude QT by default.
  input_files <- input_files[!grepl("QT\\.qs$", basename(input_files))]
}
if (length(input_files) == 0) stop("No input files found.")

canonical <- get_canonical_coarse("myeloid_refined", order = "analysis")

run_one_sample <- function(file, idx, total) {
  sample_id <- sub("\\.qs$", "", basename(file))
  group <- infer_group(sample_id)
  message(glue("[{idx}/{total}] {sample_id}: start (group={group %||% 'NA'})"))
  sample_dir <- file.path(cfg$output_root, sample_id)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

  seu <- qs::qread(file)
  if (!cfg$anno_coarse %in% colnames(seu@meta.data)) stop(glue("{sample_id}: missing {cfg$anno_coarse}"))
  if (!cfg$reduction %in% names(seu@reductions)) stop(glue("{sample_id}: missing reduction {cfg$reduction}"))

  x <- as.character(seu@meta.data[[cfg$anno_coarse]])
  unknown <- sort(unique(x[!is.na(x) & !x %in% canonical]))
  if (length(unknown) > 0) {
    stop(glue("{sample_id}: unknown {cfg$anno_coarse} values: {paste(unknown, collapse = ', ')}"))
  }
  seu@meta.data[[cfg$anno_coarse]] <- factor(x, levels = canonical)

  coords_scaled <- Seurat::Embeddings(seu, reduction = cfg$reduction) * cfg$coord_scale
  seu[[cfg$reduction]]@cell.embeddings <- coords_scaled

  caf_cells_all <- rownames(seu@meta.data)[seu@meta.data[[cfg$anno_coarse]] == cfg$caf_label]
  if (length(caf_cells_all) == 0) stop(glue("{sample_id}: no CAF cells"))

  section_prop <- compute_section_prop(seu@meta.data, cfg$anno_coarse, canonical)
  section_df <- data.frame(cell_type = canonical, fraction = as.numeric(section_prop))
  utils::write.table(section_df, file.path(sample_dir, "section_celltype_fraction.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(data.frame(r_fixed = cfg$r_fixed, coord_scale = cfg$coord_scale), file.path(sample_dir, "r_choice.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  message(glue("{sample_id}: CalNichMatrix (r={cfg$r_fixed}um, CAF={length(caf_cells_all)})"))
  t_niche <- system.time({
    niche_res <- CalNichMatrix(
      seu,
      query.cells = caf_cells_all,
      r = cfg$r_fixed,
      reduction = cfg$reduction,
      anno_col = cfg$anno_coarse,
      use_pbapply = TRUE,
      return_neighbors = FALSE
    )
  })
  message(glue("{sample_id}: niche matrix done in {round(t_niche['elapsed'], 1)}s"))
  niche_count <- align_canonical_columns(niche_res$niche.mat, canonical)
  niche_frac <- align_canonical_columns(niche_res$niche.mat.norm, canonical)
  neigh_n <- rowSums(niche_count)
  neigh_summary <- data.frame(
    sample_id = sample_id,
    group = group,
    coord_scale = cfg$coord_scale,
    r_fixed = cfg$r_fixed,
    n_query = nrow(niche_count),
    n_min = min(neigh_n),
    n_p25 = as.numeric(stats::quantile(neigh_n, 0.25)),
    n_median = as.numeric(stats::quantile(neigh_n, 0.5)),
    n_mean = mean(neigh_n),
    n_p75 = as.numeric(stats::quantile(neigh_n, 0.75)),
    n_max = max(neigh_n)
  )
  utils::write.table(
    neigh_summary,
    file.path(sample_dir, "niche_neighbor_count_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  if (neigh_summary$n_median <= 2) {
    warning(glue("{sample_id}: median neighbor count is {neigh_summary$n_median} (very low). Check coord_scale/r_fixed units."))
  }

  # section-normalized enrichment-like matrix
  denom <- section_prop[colnames(niche_frac)]
  denom[denom == 0] <- NA_real_
  niche_corr <- niche_frac / matrix(denom, nrow = nrow(niche_frac), ncol = length(denom), byrow = TRUE)
  niche_corr[!is.finite(niche_corr) | is.na(niche_corr)] <- 0

  nmf_input <- if (cfg$profile == "corrected") niche_corr else niche_frac

  fit_disp <- numeric(0)
  k_rows <- list()

  for (k in cfg$k_list) {
    message(glue("{sample_id}: k={k} runs={cfg$n_runs} subsample={cfg$subsample_frac}"))
    run_profiles <- list()
    run_assign <- list()
    run_disp <- numeric(0)

    for (run_i in seq_len(cfg$n_runs)) {
      seed_i <- cfg$seed + k * 1000L + run_i * 37L
      set.seed(seed_i)
      caf_cells <- caf_cells_all
      if (cfg$subsample_frac < 1) {
        n_take <- max(2, floor(length(caf_cells_all) * cfg$subsample_frac))
        caf_cells <- sample(caf_cells_all, n_take)
      }
      if (!is.na(cfg$max_caf_run_cells) && length(caf_cells) > cfg$max_caf_run_cells) {
        caf_cells <- sample(caf_cells, cfg$max_caf_run_cells)
      }
      caf_cells <- sort(unique(caf_cells))
      if (length(caf_cells) < 2) next

      seu_caf <- subset(seu, cells = caf_cells)
      nmf_mat <- nmf_input[caf_cells, , drop = FALSE]

      nmf_res <- RunNMF(nmf_mat, k = k, cores = 5L, seed = seed_i)
      W <- nmf_res$W
      H <- nmf_res$H

      reduction_name <- paste0("nmf_k", k)
      key_name <- paste0("NMF", k, "_")
      seu_caf[[reduction_name]] <- CreateDimReducObject(
        embeddings = W[rownames(seu_caf@meta.data), , drop = FALSE],
        assay = "RNA",
        key = key_name
      )

      cluster_name <- paste0("CAF.sub_k", k)
      seu_caf <- FindNeighbors(seu_caf, reduction = reduction_name, dims = 1:k, verbose = FALSE)
      cluster_res <- find_clusters_with_cap(
        seu_caf,
        cluster_name = cluster_name,
        resolution = cfg$cluster_resolution,
        max_clusters = cfg$max_caf_clusters,
        min_resolution = cfg$min_cluster_resolution,
        step = cfg$cluster_resolution_step,
        seed = seed_i,
        log_prefix = glue("{sample_id}: k={k} run={run_i}")
      )
      seu_caf <- cluster_res$seu
      clusters <- cluster_res$clusters
      cluster_vec <- factor(paste0("CAF-", clusters))
      names(cluster_vec) <- rownames(seu_caf@meta.data)

      run_tag <- paste0("run", run_i)
      cluster_df <- data.frame(cell_id = names(cluster_vec), caf_cluster = as.character(cluster_vec))
      utils::write.table(
        cluster_df,
        file.path(sample_dir, glue("caf_clusters_k{k}_{run_tag}.tsv")),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )

      counts <- table(cluster_vec)
      cluster_summary <- data.frame(
        sample_id = sample_id,
        group = group,
        k = k,
        run = run_i,
        cluster = names(counts),
        n_cells = as.integer(counts),
        fraction = as.numeric(counts) / length(cluster_vec)
      )
      utils::write.table(
        cluster_summary,
        file.path(sample_dir, glue("caf_cluster_summary_k{k}_{run_tag}.tsv")),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )

      avg <- cluster_average(nmf_mat, cluster_vec)
      utils::write.table(
        avg,
        file.path(sample_dir, glue("niche_average_{cfg$profile}_k{k}_{run_tag}.tsv")),
        sep = "\t",
        quote = FALSE
      )

      run_profiles[[run_tag]] <- avg
      run_assign[[run_tag]] <- cluster_vec
      run_disp <- c(run_disp, calc_fit_within_dispersion(nmf_mat, cluster_vec))
    }

    stab <- compute_stability_metrics(run_profiles, run_assign)
    fit_raw <- mean(run_disp, na.rm = TRUE)
    fit_disp <- c(fit_disp, fit_raw)
    k_rows[[as.character(k)]] <- data.frame(
      sample_id = sample_id,
      group = group,
      k = k,
      stab_profile = stab$stab_profile,
      stab_assign = stab$stab_assign,
      stab_combined = stab$stab_combined,
      fit_dispersion = fit_raw,
      stringsAsFactors = FALSE
    )
  }

  k_df <- do.call(rbind, k_rows)
  if (is.null(k_df) || nrow(k_df) == 0) stop(glue("{sample_id}: no k results"))

  # rescale fit_score to [0,1] (higher is better)
  fit_vals <- k_df$fit_dispersion
  if (all(is.na(fit_vals)) || length(unique(fit_vals[is.finite(fit_vals)])) <= 1) {
    k_df$fit_score <- 1
  } else {
    rng <- range(fit_vals, na.rm = TRUE)
    k_df$fit_score <- (rng[2] - fit_vals) / (rng[2] - rng[1] + 1e-9)
    k_df$fit_score[!is.finite(k_df$fit_score) | is.na(k_df$fit_score)] <- 0
  }
  utils::write.table(k_df, file.path(sample_dir, "k_scan.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  selected_k <- select_k_stable(k_df, delta = cfg$stab_delta, fit_min = cfg$fit_min)
  row_sel <- k_df[k_df$k == selected_k, , drop = FALSE]
  is_unstable <- as.integer(is.na(row_sel$stab_combined) || row_sel$stab_combined < 0.2 || row_sel$fit_score < cfg$fit_min)

  message(glue("{sample_id}: selected_k={selected_k} (stab={round(row_sel$stab_combined,3)}, fit={round(row_sel$fit_score,3)}, unstable={is_unstable})"))

  # Build consensus using run files for selected_k.
  run_profiles <- list()
  run_assign <- list()
  for (run_i in seq_len(cfg$n_runs)) {
    run_tag <- paste0("run", run_i)
    prof_path <- file.path(sample_dir, glue("niche_average_{cfg$profile}_k{selected_k}_{run_tag}.tsv"))
    clus_path <- file.path(sample_dir, glue("caf_clusters_k{selected_k}_{run_tag}.tsv"))
    if (!file.exists(prof_path) || !file.exists(clus_path)) next
    run_profiles[[run_tag]] <- read.delim(prof_path, row.names = 1, check.names = FALSE)
    cl_df <- read.delim(clus_path, stringsAsFactors = FALSE)
    vec <- setNames(cl_df$caf_cluster, cl_df$cell_id)
    run_assign[[run_tag]] <- vec
  }
  cons <- build_consensus_profile(run_profiles)
  cons_profile <- cons$profile
  cons_maps <- cons$maps

  # Rename local clusters to stable LC-1..LC-n based on reference order.
  lc_levels <- rownames(cons_profile)
  lc_new <- paste0("LC-", seq_along(lc_levels))
  rownames(cons_profile) <- lc_new
  for (r in names(cons_maps)) {
    m <- cons_maps[[r]]
    if (is.null(m)) next
    m2 <- m
    for (i in seq_along(lc_levels)) {
      m2[m == lc_levels[i]] <- lc_new[i]
    }
    cons_maps[[r]] <- m2
  }

  utils::write.table(
    cons_profile,
    file.path(sample_dir, glue("niche_average_{cfg$profile}_selected.tsv")),
    sep = "\t",
    quote = FALSE
  )

  cons_assign <- build_consensus_assignment(run_assign, cons_maps)
  cons_assign <- cons_assign[!is.na(cons_assign)]
  if (length(cons_assign) == 0) stop(glue("{sample_id}: consensus assignment empty"))

  caf_sel_df <- data.frame(cell_id = names(cons_assign), local_cluster = as.character(cons_assign), stringsAsFactors = FALSE)
  utils::write.table(caf_sel_df, file.path(sample_dir, "caf_clusters_selected.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  sel_counts <- table(cons_assign)
  sel_summary <- data.frame(
    sample_id = sample_id,
    group = group,
    selected_k = selected_k,
    local_cluster = names(sel_counts),
    n_cells = as.integer(sel_counts),
    fraction = as.numeric(sel_counts) / length(cons_assign),
    stringsAsFactors = FALSE
  )
  utils::write.table(sel_summary, file.path(sample_dir, "caf_cluster_summary_selected.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  k_selected <- data.frame(
    sample_id = sample_id,
    group = group,
    selected_k = selected_k,
    stab_profile = row_sel$stab_profile[1],
    stab_assign = row_sel$stab_assign[1],
    stab_combined = row_sel$stab_combined[1],
    fit_score = row_sel$fit_score[1],
    is_unstable = is_unstable,
    profile = cfg$profile,
    profile_file = glue("niche_average_{cfg$profile}_selected.tsv"),
    cluster_file = "caf_clusters_selected.tsv",
    stringsAsFactors = FALSE
  )
  utils::write.table(k_selected, file.path(sample_dir, "k_selected.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  rm(seu)
  gc()
  message(glue("{sample_id}: done"))
}

for (i in seq_along(input_files)) {
  run_one_sample(input_files[[i]], i, length(input_files))
}

message(glue("Done. Output root: {cfg$output_root}"))
