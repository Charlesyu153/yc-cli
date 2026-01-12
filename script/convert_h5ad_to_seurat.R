#!/usr/bin/env Rscript
# h5ad到Seurat格式转换脚本
# 用法: Rscript convert_h5ad_to_seurat.R [input_dir] [output_dir]

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(reticulate)
  library(glue)
  library(qs)
})

# scanpy 模块句柄（显式声明，避免 R CMD check / lintr 关于全局变量的告警）
sc <- NULL

# 配置Python环境（默认使用 spatial conda 环境；也可通过环境变量覆盖）
# - PCF_PYTHON=/path/to/python3  (推荐)
cat("检查Python环境...\n")
tryCatch({
  py_default <- "~/miniconda3/envs/spatial/bin/python3.11"
  py_path <- Sys.getenv("PCF_PYTHON", unset = py_default)
  py_path <- path.expand(py_path)

  if (file.exists(py_path)) {
    reticulate::use_python(py_path, required = TRUE)
  } else {
    # 回退到系统 python3
    py_sys <- Sys.which("python3")
    if (nzchar(py_sys)) {
      reticulate::use_python(py_sys, required = TRUE)
    } else {
      stop(glue("未找到可用的 Python：默认路径不存在且系统无 python3。\n  默认: {py_default}\n  可设环境变量 PCF_PYTHON=/path/to/python3"))
    }
  }

  reticulate::py_config()
}, error = function(e) {
  stop(glue("Python环境配置失败：{e$message}\n建议：确认 reticulate 可用，并设置 PCF_PYTHON 指向包含 scanpy 的 python。"))
})

# 导入scanpy
cat("导入scanpy...\n")
tryCatch({
  sc <<- reticulate::import("scanpy")
}, error = function(e) {
  stop("无法导入scanpy，请确保已安装: pip install scanpy")
})

#' 将h5ad格式转换为Seurat对象
#' 
#' @param h5ad_path h5ad文件路径
#' @param output_path 可选，输出R对象路径（如.qs）
#' @return 返回Seurat对象（仅包含counts数据）
h5ad_to_seurat <- function(h5ad_path, output_path = NULL) {
  
  if (!file.exists(h5ad_path)) {
    stop(glue("文件不存在: {h5ad_path}"))
  }
  
  cat(glue("正在读取h5ad文件: {h5ad_path}\n"))
  
  tryCatch({
    adata <- sc$read_h5ad(h5ad_path)
  }, error = function(e) {
    stop(glue("读取h5ad文件失败: {e$message}"))
  })
  
  # 提取 counts 矩阵（AnnData: cells x genes）
  # 注意：adata$X 常为 scipy sparse；不要在 python 对象上直接设置 rownames/colnames。
  counts_py <- adata$X

  # 若为 scipy sparse，尽量转为 CSC（更适合下游转置/列操作）
  tryCatch({
    sp <- reticulate::import("scipy.sparse", delay_load = TRUE)
    if (sp$issparse(counts_py)) {
      counts_py <- counts_py$tocsc()
    }
  }, error = function(e) {
    # scipy 不可用时忽略；后续尝试直接 py_to_r
  })

  # 转为 R 对象（优先保留稀疏矩阵，避免爆内存）
  counts <- tryCatch({
    reticulate::py_to_r(counts_py)
  }, error = function(e) {
    stop(glue("无法将 adata$X 转换为 R 矩阵（可能过大或类型不支持）：{e$message}"))
  })

  # 设置行名（细胞名）和列名（基因名）
  cell_names <- tryCatch(adata$obs_names$to_list(), error = function(e) NULL)
  gene_names <- tryCatch(adata$var_names$to_list(), error = function(e) NULL)
  if (is.null(cell_names) || is.null(gene_names)) {
    stop("无法获取 obs_names/var_names（h5ad 可能损坏或 scanpy 版本不兼容）")
  }

  if (inherits(counts, "Matrix")) {
    Matrix::Dimnames(counts) <- list(cell_names, gene_names)
  } else {
    rownames(counts) <- cell_names
    colnames(counts) <- gene_names
  }
  
  cat(glue("  矩阵维度: {nrow(counts)} 细胞 x {ncol(counts)} 基因\n"))
  
  # 提取元数据（确保为 R data.frame，并以细胞名为行名）
  meta.data <- tryCatch({
    reticulate::py_to_r(adata$obs)
  }, error = function(e) {
    stop(glue("无法将 adata$obs 转换为 R data.frame：{e$message}"))
  })
  meta.data <- as.data.frame(meta.data, stringsAsFactors = FALSE)
  rownames(meta.data) <- cell_names
  
  # 创建Seurat对象（注意：Seurat中counts是genes x cells，所以需要转置）
  # 只使用counts数据，不包含其他layer
  cat("正在创建Seurat对象（仅counts数据）...\n")
  seurat_obj <- CreateSeuratObject(
    counts = Matrix::t(counts),
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
  # 首先特别处理spatial坐标（CODEX/CosMx等空间数据）
  tryCatch({
    spatial_embedding <- adata$obsm[["spatial"]]
    if (!is.null(spatial_embedding)) {
      # 转换为R矩阵
      spatial_coords <- tryCatch({
        reticulate::py_to_r(spatial_embedding)
      }, error = function(e) {
        as.matrix(spatial_embedding)
      })
      
      # 确保是矩阵格式
      if (!is.matrix(spatial_coords)) {
        spatial_coords <- as.matrix(spatial_coords)
      }
      
      # 设置行名
      rownames(spatial_coords) <- cell_names
      
      # 检查维度
      if (ncol(spatial_coords) >= 2 && nrow(spatial_coords) > 0) {
        # 对于Seurat V5，将spatial坐标添加到metadata
        # 检查是否已有Center_X和Center_Y列（可能来自原始数据）
        if (!"Center_X" %in% colnames(seurat_obj@meta.data)) {
          seurat_obj@meta.data[["Center_X"]] <- spatial_coords[, 1]
        }
        if (!"Center_Y" %in% colnames(seurat_obj@meta.data)) {
          seurat_obj@meta.data[["Center_Y"]] <- spatial_coords[, 2]
        }
        
        # 创建spatial reduction以便用DimPlot可视化
        # 确保列名符合Seurat规范
        colnames(spatial_coords) <- c("spatial_1", "spatial_2")
        
        # 创建DimReducObject
        seurat_obj[["spatial"]] <- CreateDimReducObject(
          embeddings = spatial_coords,
          assay = "RNA",
          key = "spatial_"
        )
        
        cat(glue("  ✓ 已添加spatial坐标到metadata: Center_X, Center_Y\n"))
        cat(glue("  ✓ 已创建spatial reduction用于DimPlot可视化\n"))
      }
    }
  }, error = function(e) {
    # spatial键不存在或无法访问，忽略
    cat(glue("  未找到spatial坐标（跳过）\n"))
  })
  
  # 处理其他obsm键（降维结果，如X_umap, X_pca等）
  # 注意：由于无法直接获取所有键，我们只处理已知的常见降维结果
  # 如果需要处理其他降维结果，可以在h5ad文件中明确指定
  
  # 保存R对象（使用qs格式）
  if (!is.null(output_path)) {
    tryCatch({
      qs::qsave(seurat_obj, output_path)
      cat(glue("  ✓ 已保存qs文件: {output_path}\n"))
    }, error = function(e) {
      warning(glue("保存文件失败: {e$message}"))
    })
  }
  
  cat(glue("  ✓ 成功创建Seurat对象: {ncol(seurat_obj)} 细胞, {nrow(seurat_obj)} 基因\n"))
  
  return(seurat_obj)
}

#' 批量将h5ad文件转换为Seurat对象
#' 
#' @param input_dir 输入目录（包含h5ad文件）
#' @param output_dir 输出目录
#' @param pattern 文件名模式，默认为"*.annotated.h5ad"
#' @param save_qs 是否保存R对象（qs格式），默认为TRUE
#' @return 返回转换结果列表
batch_h5ad_to_seurat <- function(input_dir, output_dir, pattern = "*.annotated.h5ad", save_qs = TRUE) {
  
  if (!dir.exists(input_dir)) {
    stop(glue("输入目录不存在: {input_dir}"))
  }
  
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(glue("已创建输出目录: {output_dir}\n"))
  }
  
  # 查找所有h5ad文件
  h5ad_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  
  if (length(h5ad_files) == 0) {
    stop(glue("在 {input_dir} 中未找到匹配 {pattern} 的文件"))
  }
  
  cat(glue("\n找到 {length(h5ad_files)} 个h5ad文件\n"))
  cat(paste0(rep("=", 60), collapse = ""), "\n")
  
  # 存储结果
  results <- list()
  
  # 逐个转换
  for (i in seq_along(h5ad_files)) {
    h5ad_file <- h5ad_files[i]
    file_name <- basename(h5ad_file)
    sample_name <- gsub("\\.annotated\\.h5ad$|\\.h5ad$", "", file_name)
    
    cat(glue("\n[{i}/{length(h5ad_files)}] 正在处理: {file_name}\n"))
    cat(paste0(rep("-", 60), collapse = ""), "\n")
    
    tryCatch({
      # 确定输出路径
      if (save_qs) {
        output_path <- file.path(output_dir, paste0(sample_name, ".qs"))
      } else {
        output_path <- NULL
      }
      
      # 转换
      seurat_obj <- h5ad_to_seurat(h5ad_file, output_path = output_path)
      
      results[[sample_name]] <- list(
        status = "success",
        seurat_obj = seurat_obj,
        output_path = output_path,
        n_cells = ncol(seurat_obj),
        n_genes = nrow(seurat_obj)
      )
      
      cat(glue("✓ 成功转换: {sample_name}\n"))
      
    }, error = function(e) {
      cat(glue("✗ 转换失败: {sample_name}\n"))
      cat(glue("  错误信息: {e$message}\n"))
      results[[sample_name]] <<- list(
        status = "error",
        error = e$message
      )
    })
  }
  
  # 汇总结果
  cat("\n", paste0(rep("=", 60), collapse = ""), "\n")
  cat("转换完成汇总\n")
  cat(paste0(rep("=", 60), collapse = ""), "\n")
  
  success_count <- sum(sapply(results, function(x) x$status == "success"))
  error_count <- length(results) - success_count
  
  cat(glue("成功: {success_count}/{length(h5ad_files)}\n"))
  if (error_count > 0) {
    cat(glue("失败: {error_count}/{length(h5ad_files)}\n"))
  }
  
  # 显示详细信息
  if (success_count > 0) {
    cat("\n成功转换的样本:\n")
    for (sample_name in names(results)) {
      if (results[[sample_name]]$status == "success") {
        r <- results[[sample_name]]
        cat(glue("  ✓ {sample_name}: {r$n_cells} 细胞, {r$n_genes} 基因\n"))
      }
    }
  }
  
  if (error_count > 0) {
    cat("\n转换失败的样本:\n")
    for (sample_name in names(results)) {
      if (results[[sample_name]]$status == "error") {
        r <- results[[sample_name]]
        cat(glue("  ✗ {sample_name}: {r$error}\n"))
      }
    }
  }
  
  return(results)
}

# 主程序
main <- function() {
  # 从命令行参数获取输入输出目录，或使用默认值
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) >= 2) {
    input_dir <- args[1]
    output_dir <- args[2]
  } else {
    # 默认值
    input_dir <- "/home/jacekyu/PCF/spatial"
    output_dir <- "/home/jacekyu/PCF/data"
  }
  
  cat(glue("\n=== h5ad到Seurat批量转换工具 ===\n"))
  cat(glue("输入目录: {input_dir}\n"))
  cat(glue("输出目录: {output_dir}\n\n"))
  
  # 执行批量转换
  results <- batch_h5ad_to_seurat(
    input_dir = input_dir,
    output_dir = output_dir,
    pattern = "*.annotated.h5ad",
    save_qs = TRUE
  )
  
  # 保存转换结果摘要
  summary_path <- file.path(output_dir, "conversion_summary.qs")
  qs::qsave(results, summary_path)
  cat(glue("\n转换结果摘要已保存到: {summary_path}\n"))
}

# 如果直接运行此脚本，执行主程序
if (!interactive()) {
  main()
}

