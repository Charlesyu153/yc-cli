#' CalNichMatrix: Compute neighborhood cell-type composition around query cells (Visium data)
#'
#' @param seu Seurat object. 需要包含 GetTissueCoordinates 和细胞类型权重矩阵。
#' @param query.cells character vector，待评估邻域的细胞（必须在 seu 中）。
#' @param r numeric，邻域半径（与坐标同单位，例：400 μm）。
#' @param anno_col character vector，细胞类型权重矩阵的列名（如 c("Tumor", "CAF", ...)）。
#' @param use_pbapply logical，是否用 pbapply 显示进度条；若未安装将自动回退为 lapply。
#' @param return_neighbors logical，是否在结果里返回 frNN 的邻居索引。
#'
#' @return list，包含
#'   - niche.mat: 邻域内各类型平均权重矩阵（行=query cells，列=cell types）
#'   - niche.mat.norm: niche.mat 的行归一化（各类型比例，行和=1）
#'   - type.levels: 细胞类型名称（列顺序）
#'   - neighbors (可选): frNN 返回的邻居索引列表
#'
#' @examples
#' # Visium 数据
#' ref.celltypes <- colnames(seu@meta.data)[4:14]
#' res <- CalNichMatrix(seu,
#'   query.cells = query.cells, r = 400,
#'   anno_col = ref.celltypes
#' )
CalNichMatrix <- function(seu,
                          query.cells,
                          r = 400,
                          anno_col,
                          use_pbapply = TRUE,
                          return_neighbors = FALSE) {
  # -------------- 基本检查 --------------
  stopifnot(is.character(query.cells), length(query.cells) > 0)
  stopifnot(is.character(anno_col), length(anno_col) > 0)

  # 检查所有列是否存在
  missing_cols <- setdiff(anno_col, colnames(seu@meta.data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Columns not found in seu@meta.data: %s", paste(missing_cols, collapse = ", ")))
  }

  # -------------- 提取坐标 --------------
  coords <- Seurat::GetTissueCoordinates(seu)[, c("x", "y"), drop = FALSE]
  ref.cells <- rownames(coords)
  if (is.null(ref.cells)) {
    # 如果 GetTissueCoordinates 返回的没有 rownames，使用 meta.data 的 rownames
    ref.cells <- rownames(seu@meta.data)
    rownames(coords) <- ref.cells
  }

  # 仅保留在对象中的 query cells
  query.cells <- intersect(query.cells, ref.cells)
  if (length(query.cells) == 0) stop("No query.cells found in the object.")

  # -------------- 提取细胞类型权重矩阵 --------------
  weights <- seu@meta.data[, anno_col, drop = FALSE]
  type.levels <- colnames(weights)

  # -------------- 邻域查询（半径 r） --------------
  Q <- coords[query.cells, , drop = FALSE]
  nbrs <- dbscan::frNN(coords, eps = r, query = Q, sort = FALSE)

  # -------------- 计算各邻域的平均权重 --------------
  # 进度条函数准备
  apply_fun <- if (use_pbapply && requireNamespace("pbapply", quietly = TRUE)) {
    pbapply::pblapply
  } else {
    lapply
  }

  niche.list <- apply_fun(nbrs$id, function(nn) {
    if (length(nn) == 0) {
      # 无邻居时给全零向量
      rep(0, length(type.levels))
    } else {
      colMeans(weights[ref.cells[nn], , drop = FALSE])
    }
  })

  niche.mat <- do.call(rbind, niche.list)
  rownames(niche.mat) <- rownames(Q)
  colnames(niche.mat) <- type.levels

  # -------------- 归一化 --------------
  # 行归一化：邻域内部组成（比例）
  row_sums <- rowSums(niche.mat)
  # 避免除零：对行和为0的行保持 0
  niche.mat.norm <- niche.mat
  nz <- row_sums > 0
  niche.mat.norm[nz, ] <- niche.mat[nz, , drop = FALSE] / row_sums[nz]

  # -------------- 返回结果 --------------
  out <- list(
    niche.mat = niche.mat,
    niche.mat.norm = niche.mat.norm,
    type.levels = type.levels
  )
  if (return_neighbors) {
    nn <- sapply(nbrs$id, function(id) {
      ref.cells[id]
    })
    names(nn) <- names(nbrs$id)
    out$neighbors <- nn
  }
  return(out)
}


RunNMF <- function(niche.mat, k = 5, cores = 5, seed = 1024) {
  # RcppML::getRcppMLthreads()
  RcppML::setRcppMLthreads(cores) ## 使用5个线程
  model <- RcppML::nmf(niche.mat, k = k, verbose = F, seed = seed) ## 技术细节：k的数量？

  H <- model$h
  rownames(H) <- paste0("factor_", 1:nrow(H))
  colnames(H) <- colnames(niche.mat)
  # round(H, 4)

  W <- model$w
  rownames(W) <- rownames(niche.mat)
  colnames(W) <- rownames(H)
  # head(W)

  return(list(W = W, H = H))
}
