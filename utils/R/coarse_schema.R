default_myeloid_fine_cols <- function() {
  c(
    "annotation_fine",
    "cell_type_lvl2",
    "cell_type_lvl1",
    "annotation"
  )
}

parse_csv_list <- function(x) {
  if (is.null(x)) return(character(0))
  x <- as.character(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(character(0))
  if (length(x) == 1) {
    out <- trimws(strsplit(x, ",")[[1]])
  } else {
    out <- trimws(x)
  }
  out <- out[out != ""]
  unique(out)
}

get_canonical_coarse <- function(schema = c("base", "myeloid_refined"),
                                order = c("analysis", "grid")) {
  schema <- match.arg(schema)
  order <- match.arg(order)

  myeloid_levels <- if (schema == "myeloid_refined") {
    c("Macrophage", "Neutrophil", "Myeloid_Mast", "Myeloid_other")
  } else {
    "Myeloid"
  }

  if (order == "grid") {
    return(c(
      "Stromal_other",
      "CAF",
      "Endothelial",
      "Tumor",
      myeloid_levels,
      "T_cell",
      "B_cell"
    ))
  }

  c(
    "B_cell",
    "CAF",
    "Endothelial",
    myeloid_levels,
    "Stromal_other",
    "T_cell",
    "Tumor"
  )
}

refine_myeloid_in_coarse <- function(meta,
                                     anno_coarse = "annotation_coarse",
                                     fine_cols = default_myeloid_fine_cols()) {
  if (!anno_coarse %in% colnames(meta)) {
    stop(sprintf("Column '%s' not found in meta.", anno_coarse))
  }

  coarse <- as.character(meta[[anno_coarse]])
  out <- coarse

  idx <- which(!is.na(coarse) & coarse == "Myeloid")
  if (length(idx) == 0) return(out)

  fine_cols <- intersect(parse_csv_list(fine_cols), colnames(meta))
  if (length(fine_cols) == 0) return(out)

  chosen <- rep(NA_character_, length(idx))
  for (col in fine_cols) {
    v <- as.character(meta[[col]][idx])
    take <- is.na(chosen) & !is.na(v) & v != ""
    chosen[take] <- v[take]
  }

  low <- tolower(chosen)
  mapped <- rep(NA_character_, length(chosen))
  is_valid <- !is.na(low) & low != ""

  mapped[is_valid & grepl("mast", low)] <- "Myeloid_Mast"
  mapped[is.na(mapped) & is_valid & (grepl("neut", low) | grepl("\\bpmn\\b", low))] <- "Neutrophil"
  mapped[is.na(mapped) & is_valid & (grepl("macroph", low) | grepl("macro", low) | grepl("monocyte", low) | grepl("\\bmono\\b", low))] <- "Macrophage"

  keep <- !is.na(mapped)
  out[idx[keep]] <- mapped[keep]
  out[idx[!keep]] <- "Myeloid_other"
  out
}

rename_myeloid_refined_labels <- function(x) {
  x <- as.character(x)
  x[x == "Myeloid_Macrophage"] <- "Macrophage"
  x[x == "Myeloid_Neutrophil"] <- "Neutrophil"
  x[x == "Myeloid"] <- "Myeloid_other"
  x
}
