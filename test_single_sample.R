#!/usr/bin/env Rscript
# 测试单个样本的h5ad到Seurat转换

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(reticulate)
  library(glue)
  library(qs)
})

# 配置Python环境（使用spatial conda环境）
cat("检查Python环境...\n")
reticulate::use_python("~/miniconda3/envs/spatial/bin/python3.11", required = TRUE)
reticulate::py_config()

# 导入scanpy
cat("导入scanpy...\n")
sc <- reticulate::import("scanpy")

# 定义转换函数（从主脚本复制）
h5ad_to_seurat <- function(h5ad_path, output_path = NULL) {
  
  if (!file.exists(h5ad_path)) {
    stop(glue("文件不存在: {h5ad_path}"))
  }
  
  cat(glue("正在读取h5ad文件: {h5ad_path}\n"))
  adata <- sc$read_h5ad(h5ad_path)
  
  # 提取counts矩阵
  counts <- adata$X
  
  # 设置行名（细胞名）和列名（基因名）
  rownames(counts) <- adata$obs_names$to_list()
  colnames(counts) <- adata$var_names$to_list()
  
  cat(glue("  矩阵维度: {nrow(counts)} 细胞 x {ncol(counts)} 基因\n"))
  
  # 提取元数据
  meta.data <- adata$obs
  
  # 创建Seurat对象（只使用counts数据）
  cat("正在创建Seurat对象（仅counts数据）...\n")
  seurat_obj <- CreateSeuratObject(
    counts = t(counts),
    meta.data = meta.data
  )
  
  # 确保只保留counts layer，移除其他layer（如果有）
  # Seurat v5中，layers是list，需要检查并移除
  if (length(seurat_obj@assays$RNA@layers) > 0) {
    layer_names <- names(seurat_obj@assays$RNA@layers)
    for (layer_name in layer_names) {
      if (layer_name != "counts") {
        seurat_obj@assays$RNA@layers[[layer_name]] <- NULL
      }
    }
  }
  # 移除scale.data（Seurat v5中可能不存在，需要安全检查）
  tryCatch({
    if (length(seurat_obj@assays$RNA@scale.data) > 0) {
      seurat_obj@assays$RNA@scale.data <- matrix(nrow = 0, ncol = 0)
    }
  }, error = function(e) {
    # scale.data slot不存在，忽略
  })
  
  # 处理obsm中的内容（降维结果和spatial坐标）
  obsm_keys <- names(adata$obsm)
  
  if (length(obsm_keys) > 0) {
    cat(glue("  发现 {length(obsm_keys)} 个obsm键\n"))
    
    for (obsm_key in obsm_keys) {
      tryCatch({
        embedding <- adata$obsm[[obsm_key]]
        
        # 特别处理spatial坐标（CODEX/CosMx等空间数据）
        if (obsm_key == "spatial") {
          if (!is.null(embedding)) {
            embedding_dim <- tryCatch({
              dim(embedding)
            }, error = function(e) NULL)
            
            if (!is.null(embedding_dim) && length(embedding_dim) == 2 && 
                !is.na(embedding_dim[1]) && !is.na(embedding_dim[2]) &&
                embedding_dim[1] > 0 && embedding_dim[2] >= 2) {
              
              # 转换为R矩阵
              spatial_coords <- as.matrix(embedding)
              rownames(spatial_coords) <- adata$obs_names$to_list()
              
              # 对于Seurat V5，将spatial坐标添加到metadata
              # 检查是否已有Center_X和Center_Y列（可能来自原始数据）
              if (!"Center_X" %in% colnames(seurat_obj@meta.data)) {
                seurat_obj@meta.data[["Center_X"]] <- spatial_coords[, 1]
              }
              if (!"Center_Y" %in% colnames(seurat_obj@meta.data)) {
                seurat_obj@meta.data[["Center_Y"]] <- spatial_coords[, 2]
              }
              
              cat(glue("  ✓ 已添加spatial坐标到metadata: Center_X, Center_Y\n"))
            }
          }
        } else {
          # 其他obsm键作为降维结果处理（如X_umap, X_pca等）
          # 安全地检查embedding是否有效
          if (!is.null(embedding)) {
            embedding_dim <- tryCatch({
              dim(embedding)
            }, error = function(e) NULL)
            
            if (!is.null(embedding_dim) && length(embedding_dim) == 2 && 
                !is.na(embedding_dim[1]) && !is.na(embedding_dim[2]) &&
                embedding_dim[1] > 0 && embedding_dim[2] > 0) {
              
              # 转换为R矩阵
              embedding_matrix <- as.matrix(embedding)
              rownames(embedding_matrix) <- adata$obs_names$to_list()
              reduction_name <- gsub("^X_", "", tolower(obsm_key))
              n_dims <- ncol(embedding_matrix)
              colnames(embedding_matrix) <- paste0(toupper(reduction_name), "_", 1:n_dims)
              
              seurat_obj[[reduction_name]] <- CreateDimReducObject(
                embeddings = embedding_matrix,
                assay = "RNA"
              )
              
              cat(glue("  ✓ 已添加降维结果: {reduction_name} ({n_dims}维)\n"))
            }
          }
        }
      }, error = function(e) {
        cat(glue("  ✗ 跳过 {obsm_key}: {e$message}\n"))
      })
    }
  }
  
  if (!is.null(output_path)) {
    qs::qsave(seurat_obj, output_path)
    cat(glue("  ✓ 已保存qs文件: {output_path}\n"))
  }
  
  cat(glue("  ✓ 成功创建Seurat对象: {ncol(seurat_obj)} 细胞, {nrow(seurat_obj)} 基因\n"))
  
  return(seurat_obj)
}

# 测试单个样本（选择最小的文件20P）
test_h5ad <- "/home/jacekyu/PCF/spatial/20P.annotated.h5ad"
output_dir <- "/home/jacekyu/PCF/data"

if (!file.exists(test_h5ad)) {
  cat(glue("错误: 测试文件不存在: {test_h5ad}\n"))
  quit(status = 1)
}

cat("\n=== 开始测试单个样本转换 ===\n")
cat(glue("测试文件: {test_h5ad}\n\n"))

# 创建输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 转换
test_seurat <- h5ad_to_seurat(
  test_h5ad,
  output_path = file.path(output_dir, "test_20P.qs")
)

# 详细检查结果
cat("\n=== 转换结果详细检查 ===\n")
cat(glue("细胞数量: {ncol(test_seurat)}\n"))
cat(glue("基因数量: {nrow(test_seurat)}\n"))
cat(glue("元数据列数: {ncol(test_seurat@meta.data)}\n"))
cat(glue("降维结果: {paste(Reductions(test_seurat), collapse = ', ')}\n"))

# 检查spatial坐标
cat("\n=== Spatial坐标检查 ===\n")
if ("Center_X" %in% colnames(test_seurat@meta.data) && 
    "Center_Y" %in% colnames(test_seurat@meta.data)) {
  cat("✓ Spatial坐标已保留在metadata中\n")
  x_min <- min(test_seurat@meta.data$Center_X, na.rm=TRUE)
  x_max <- max(test_seurat@meta.data$Center_X, na.rm=TRUE)
  y_min <- min(test_seurat@meta.data$Center_Y, na.rm=TRUE)
  y_max <- max(test_seurat@meta.data$Center_Y, na.rm=TRUE)
  cat(glue("  Center_X范围: {round(x_min, 2)} - {round(x_max, 2)}\n"))
  cat(glue("  Center_Y范围: {round(y_min, 2)} - {round(y_max, 2)}\n"))
  cat("  前5个细胞的坐标:\n")
  print(head(test_seurat@meta.data[, c("Center_X", "Center_Y")], 5))
} else {
  cat("✗ 警告: 未找到spatial坐标（Center_X, Center_Y）\n")
  cat("  可能原因: h5ad文件中obsm['spatial']不存在或为空\n")
}

# 显示元数据预览
cat("\n元数据预览（前3行）:\n")
print(head(test_seurat@meta.data, 3))

# 显示基因和细胞名预览
cat("\n基因名预览（前5个）:\n")
print(head(rownames(test_seurat), 5))

cat("\n细胞名预览（前5个）:\n")
print(head(colnames(test_seurat), 5))

# 检查counts矩阵
cat("\nCounts矩阵统计:\n")
counts_data <- GetAssayData(test_seurat, assay = "RNA", layer = "counts")
cat(glue("  矩阵类型: {class(counts_data)}\n"))
cat(glue("  非零值数量: {sum(counts_data != 0)}\n"))
cat(glue("  总元素数: {length(counts_data)}\n"))
cat(glue("  稀疏度: {round((1 - sum(counts_data != 0) / length(counts_data)) * 100, 2)}%\n"))

cat("\n✓ 测试成功完成！\n")
cat(glue("转换后的Seurat对象已保存到: {file.path(output_dir, 'test_20P.qs')}\n"))

