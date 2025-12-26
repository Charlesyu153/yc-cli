#' CalNichMatrix: Compute neighborhood cell-type composition around query cells
#'
#' @param seu Seurat object. 需要包含 reduction (默认 "spatial") 的二维坐标
#'            和细胞类型注释（默认列名 "annotation"）。
#' @param query.cells character vector，待评估邻域的细胞（必须在 seu 中）。
#' @param r numeric，邻域半径（与 coords 同单位，例：80）。
#' @param reduction character，用于邻域查询的坐标降维名，默认 "spatial"。
#' @param anno_col character，细胞类型注释所在的 meta.data 列名，默认 "annotation"。
#' @param use_pbapply logical，是否用 pbapply 显示进度条；若未安装将自动回退为 lapply。
#' @param return_neighbors logical，是否在结果里返回 frNN 的邻居索引。
#'
#' @return list，包含
#'   - niche.mat: 邻域内各类型计数矩阵（行=query cells，列=cell types）
#'   - niche.mat.norm: niche.mat 的行归一化（各类型比例，行和=1）
#'   - niche.mat.correct: 在 niche.mat.norm 基础上再除以“整体组成”得到的校正矩阵
#'   - type.levels: 因子水平（列顺序）
#'   - neighbors (可选): frNN 返回的邻居索引列表
#'
#' @examples
#' # 假设 Idents(seu) 已设置且你有 query.cells 向量
#' res <- CalNichMatrix(seu, query.cells = CellsByIdentities(seu)$CAF, r = 80)
CalNichMatrix <- function(seu,
                          query.cells,
                          r = 80,
                          reduction = "spatial",
                          anno_col = "annotation",
                          use_pbapply = TRUE,
                          return_neighbors = FALSE) {
  # -------------- 基本检查 --------------
  stopifnot(is.character(query.cells), length(query.cells) > 0)
  if (!anno_col %in% colnames(seu@meta.data)) {
    stop(sprintf("Column '%s' not found in seu@meta.data.", anno_col))
  }
  if (!reduction %in% names(seu@reductions)) {
    stop(sprintf("Reduction '%s' not found in seu@reductions.", reduction))
  }
  
  # -------------- 提取坐标与注释 --------------
  coords <- Seurat::Embeddings(seu, reduction = reduction)
  ref.cells <- rownames(coords)
  if (is.null(ref.cells)) stop("Embeddings must have rownames (cell names).")
  
  # 仅保留在对象中的 query cells
  query.cells <- intersect(query.cells, ref.cells)
  if (length(query.cells) == 0) stop("No query.cells found in the object.")
  
  annotation <- seu@meta.data[[anno_col]]
  # 转为因子并固定水平顺序
  if (!is.factor(annotation)) annotation <- factor(annotation)
  type.levels <- levels(annotation)
  
  # 映射到整数（为 tabulate 高效统计做准备）
  # names(type.int) = cell names（与 ref.cells 索引相容）
  type.int <- setNames(match(as.character(annotation), type.levels),
                       rownames(seu@meta.data))
  
  # -------------- 邻域查询（半径 r） --------------
  Q <- coords[query.cells, , drop = FALSE]
  nbrs <- dbscan::frNN(coords, eps = r, query = Q, sort = FALSE)
  
  # -------------- 统计各邻域的 cell-type 计数 --------------
  # 进度条函数准备
  apply_fun <- if (use_pbapply && requireNamespace("pbapply", quietly = TRUE)) {
    pbapply::pblapply
  } else {
    lapply
  }
  
  # 对每个 query cell 的邻居索引，统计各类型计数（tabulate 基于整数标签）
  # 注意：ref.cells[nn] 将邻居索引映射回细胞名，再取 type.int
  niche.list <- apply_fun(nbrs$id, function(nn) {
    if (length(nn) == 0) {
      # 无邻居时给全零向量
      integer(length(type.levels))
    } else {
      tabulate(type.int[ref.cells[nn]], nbins = length(type.levels))
    }
  })
  
  niche.mat <- do.call(rbind, niche.list)
  rownames(niche.mat) <- rownames(Q)
  colnames(niche.mat) <- type.levels
  
  # -------------- 归一化与样本层级校正 --------------
  # 行归一化：邻域内部组成（比例）
  row_sums <- rowSums(niche.mat)
  # 避免除零：对行和为0的行保持 0
  niche.mat.norm <- niche.mat
  nz <- row_sums > 0
  niche.mat.norm[nz, ] <- niche.mat[nz, , drop = FALSE] / row_sums[nz]
  
  # “整张切片”的总体组成（与原脚本一致：先把所有邻域计数求和后再归一化）
  niche.mat.section <- colSums(niche.mat)
  total <- sum(niche.mat.section)
  if (total == 0) {
    # 极端情况：所有邻域都没有邻居
    section.prop <- rep(NA_real_, length(type.levels))
  } else {
    section.prop <- niche.mat.section / total
  }
  names(section.prop) <- type.levels
  
  # 校正矩阵：每行比例 / 全局组成比例
  # 对全局比例为0的列，结果设为 NA 以避免除零
  denom <- section.prop
  denom[denom == 0] <- NA_real_
  # -------------- 返回结果 --------------
  out <- list(
    niche.mat = niche.mat,
    niche.mat.norm = niche.mat.norm,
    type.levels = type.levels
  )
  if (return_neighbors) {
    nn = sapply(nbrs$id, function(id) {
      rownames(coords)[id]
    })
    names(nn) = names(nbrs$id)
    out$neighbors <- nn
  }
  return(out)
}


RunNMF = function(niche.mat, k=5, cores=5, seed=1024) {
  # RcppML::getRcppMLthreads()
  RcppML::setRcppMLthreads(cores) ## 使用5个线程
  model <- RcppML::nmf(niche.mat, k = k, verbose = F, seed = seed) ## 技术细节：k的数量？
  
  H = model$h
  rownames(H) = paste0("factor_",1:nrow(H))
  colnames(H) = colnames(niche.mat)
  # round(H, 4)
  
  W = model$w
  rownames(W) = rownames(niche.mat)
  colnames(W) = rownames(H)
  # head(W)
  
  return(list(W=W, H=H))
}
