CalQCMetrics <- function(seu) {
    # Load required packages
    suppressPackageStartupMessages({
      library(SingleCellExperiment)
      library(celda)
    })

  log_msg <- function(msg) {
    message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  }

  log_msg("Step 1: Calculating cell cycle score ...")
  seu <- NormalizeData(seu)
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes)

  log_msg("Step 2: Calculating doublet score (scrublet) ...")
  sc <- reticulate::import("scanpy")
  counts <- LayerData(seu, assay = "RNA", layer = "counts")
  obs <- seu@meta.data
  var <- data.frame(row.names = rownames(counts), gene_symbols = rownames(counts))
  adata <- sc$AnnData(X = t(counts), obs = obs, var = var)
  sc$external$pp$scrublet(adata)
  meta.data <- adata$obs[, c("doublet_score", "predicted_doublet")]
  seu <- AddMetaData(seu, meta.data)

  log_msg("Step 3: Calculating mitochondrial gene expression ...")
  seu[["mito_percent"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

  log_msg("Step 4: Calculating ambient RNA contamination (decontX) ...")
  counts <- LayerData(seu, assay = "RNA", layer = "counts")
  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = obs,
    rowData = var
  )
  sce <- decontX(sce)
  meta.data2 <- colData(sce)
  seu[["decontX_contamination"]] <- meta.data2[, "decontX_contamination"]

  log_msg("All steps finished.")
  return(seu)
}
