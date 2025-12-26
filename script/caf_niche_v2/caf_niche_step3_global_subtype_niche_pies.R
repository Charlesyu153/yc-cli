#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glue)
  library(ggplot2)
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

argv <- commandArgs(trailingOnly = TRUE)
opts <- parse_args(argv)

defaults <- list(
  input = file.path(getwd(), "data", "caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2", "step2", "global_corrected", "global_subtype_centroids_corrected.tsv"),
  output = file.path(getwd(), "data", "caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2", "step3", "global_subtype_niche_pies.pdf"),
  coarse_schema = "myeloid_refined",
  zero_negative = TRUE,
  drop_zeros = TRUE,
  title = "Global CAF subtype niche composition (normalized centroid)"
)

if (!is.null(opts[["help"]]) || !is.null(opts[["h"]])) {
  cat("Plot global CAF subtype niche composition pies from Step2 centroids.\n")
  cat(glue("  --input {defaults$input}\n"))
  cat(glue("  --output {defaults$output}\n"))
  cat(glue("  --coarse_schema {defaults$coarse_schema}\n"))
  cat(glue("  --zero_negative {defaults$zero_negative}\n"))
  cat(glue("  --drop_zeros {defaults$drop_zeros}\n"))
  cat(glue("  --title {defaults$title}\n"))
  quit(status = 0)
}

cfg <- list(
  input = get_opt(opts, "input", defaults$input),
  output = get_opt(opts, "output", defaults$output),
  coarse_schema = tolower(get_opt(opts, "coarse_schema", defaults$coarse_schema)),
  zero_negative = get_opt(opts, "zero_negative", defaults$zero_negative, to_bool),
  drop_zeros = get_opt(opts, "drop_zeros", defaults$drop_zeros, to_bool),
  title = get_opt(opts, "title", defaults$title)
)

if (!file.exists(cfg$input)) stop(glue("Missing input: {cfg$input}"))
dir.create(dirname(cfg$output), recursive = TRUE, showWarnings = FALSE)

root_dir <- get_root_dir()
source(file.path(root_dir, "utils", "R", "coarse_schema.R"))
canonical <- get_canonical_coarse(cfg$coarse_schema, order = "grid")
celltype_colors <- get_celltype_colors(cfg$coarse_schema)

cent <- read.delim(cfg$input, check.names = FALSE)
if (ncol(cent) == 0) stop("Empty centroid file.")
if (is.null(rownames(cent))) {
  if ("global_subtype" %in% colnames(cent)) {
    rownames(cent) <- cent$global_subtype
    cent$global_subtype <- NULL
  } else {
    stop("Centroid file must have rownames (subtypes).")
  }
}

present <- intersect(canonical, colnames(cent))
if (length(present) == 0) stop("No canonical coarse columns found in centroid file.")
cent <- as.matrix(cent[, present, drop = FALSE])

long <- list()
for (subtype in rownames(cent)) {
  vals <- as.numeric(cent[subtype, ])
  names(vals) <- colnames(cent)
  if (cfg$zero_negative) vals[vals < 0] <- 0
  if (cfg$drop_zeros) {
    keep <- vals > 0
    vals <- vals[keep]
  }
  total <- sum(vals)
  if (!is.finite(total) || total <= 0) next
  frac <- vals / total
  df <- data.frame(
    global_subtype = subtype,
    cell_type = names(frac),
    fraction = as.numeric(frac),
    stringsAsFactors = FALSE
  )
  long[[subtype]] <- df
}
df_all <- do.call(rbind, long)
if (is.null(df_all) || nrow(df_all) == 0) stop("No non-zero centroid entries to plot.")

df_all$global_subtype <- factor(df_all$global_subtype, levels = sort(unique(df_all$global_subtype)))
df_all$cell_type <- factor(df_all$cell_type, levels = canonical)
colors <- celltype_colors[levels(df_all$cell_type)]

p <- ggplot(df_all, aes(x = "", y = fraction, fill = cell_type)) +
  geom_col(width = 1, color = NA) +
  coord_polar(theta = "y") +
  facet_wrap(~global_subtype, nrow = 1) +
  scale_fill_manual(values = colors, drop = FALSE) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  ) +
  labs(title = cfg$title)

ggsave(cfg$output, p, width = 3.2 * nlevels(df_all$global_subtype) + 3.5, height = 3.6, device = grDevices::cairo_pdf)

out_tsv <- sub("\\.pdf$", ".tsv", cfg$output)
utils::write.table(df_all, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message(glue("Done. Output: {cfg$output}"))
message(glue("Data: {out_tsv}"))
