#!/usr/bin/env Rscript
# 测试h5ad到Seurat的转换功能

# 加载必要的包
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(reticulate)
  library(glue)
})

# 配置Python环境
cat("检查Python环境...\n")
reticulate::py_config()

# 导入scanpy
cat("导入scanpy...\n")
sc <- reticulate::import("scanpy")

# 定义转换函数
h5ad_to_seurat <- function(h5ad_path, output_path = NULL) {
  
  cat(glue("正在读取h5ad文件: {h5ad_path}\n"))
  adata <- sc$read_h5ad(h5ad_path)
  
  # 提取counts矩阵
  counts <- adata$X
  
  # 设置行名（细胞名）和列名（基因名）
  rownames(counts) <- adata$obs_names$to_list()
  colnames(counts) <- adata$var_names$to_list()
  
  # 提取元数据
  meta.data <- adata$obs
  
  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(
    counts = t(counts),
    meta.data = meta.data
  )
  
  # 添加降维结果（如果存在）
  obsm_keys <- tryCatch({
    builtins <- reticulate::import_builtins()
    reticulate::py_to_r(builtins$list(adata$obsm$keys()))
  }, error = function(e) {
    character(0)
  })
  obsm_keys <- as.character(obsm_keys)
  
  for (obsm_key in obsm_keys) {
    embedding <- tryCatch({
      reticulate::py_to_r(adata$obsm[[obsm_key]])
    }, error = function(e) {
      NULL
    })
    
    if (is.null(embedding)) {
      next
    }
    if (!is.matrix(embedding)) {
      embedding <- tryCatch(as.matrix(embedding), error = function(e) NULL)
    }
    if (is.null(embedding) || length(dim(embedding)) != 2 || nrow(embedding) == 0) {
      next
    }
    
    rownames(embedding) <- adata$obs_names$to_list()
    reduction_name <- gsub("^X_", "", tolower(obsm_key))
    n_dims <- ncol(embedding)
    colnames(embedding) <- paste0(toupper(reduction_name), "_", 1:n_dims)
    
    seurat_obj[[reduction_name]] <- CreateDimReducObject(
      embeddings = embedding,
      assay = "RNA"
    )
    
    cat(glue("已添加降维结果: {reduction_name}\n"))
  }
  
  if (!is.null(output_path)) {
    saveRDS(seurat_obj, output_path)
    cat(glue("已保存rds文件: {output_path}\n"))
  }
  
  cat(glue("成功创建Seurat对象，包含 {ncol(seurat_obj)} 个细胞和 {nrow(seurat_obj)} 个基因\n"))
  
  return(seurat_obj)
}

# 测试单个样本（选择最小的文件）
test_h5ad <- "/home/jacekyu/PCF/spatial/20P.annotated.h5ad"

if (file.exists(test_h5ad)) {
  cat("\n=== 开始测试单个样本转换 ===\n")
  
  # 创建输出目录
  output_dir <- "/home/jacekyu/PCF/data"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 转换
  test_seurat <- h5ad_to_seurat(
    test_h5ad,
    output_path = file.path(output_dir, "test_20P.rds")
  )
  
  # 检查结果
  cat("\n=== 转换结果检查 ===\n")
  cat(glue("细胞数量: {ncol(test_seurat)}\n"))
  cat(glue("基因数量: {nrow(test_seurat)}\n"))
  cat(glue("元数据列数: {ncol(test_seurat@meta.data)}\n"))
  cat(glue("降维结果: {paste(Reductions(test_seurat), collapse = ', ')}\n"))
  
  # 显示元数据预览
  cat("\n元数据预览:\n")
  print(head(test_seurat@meta.data, 3))
  
  # 显示基因和细胞名预览
  cat("\n基因名预览（前5个）:\n")
  print(head(rownames(test_seurat), 5))
  
  cat("\n细胞名预览（前5个）:\n")
  print(head(colnames(test_seurat), 5))
  
  cat("\n✓ 测试成功完成！\n")
  
} else {
  cat(glue("错误: 测试文件不存在: {test_h5ad}\n"))
  quit(status = 1)
}
