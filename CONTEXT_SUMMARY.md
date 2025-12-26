# h5ad到Seurat V5转换工具 - 上下文总结

## 项目概述

将CODEX空间蛋白组多重染色数据（h5ad格式）转换为Seurat V5对象，保留spatial坐标信息。

## 已完成功能

### 1. 核心转换函数
- **`h5ad_to_seurat()`**: 将h5ad文件转换为Seurat V5对象
  - 保留counts数据（仅原始数据）
  - 正确处理spatial坐标（`obsm["spatial"]` → `Center_X`, `Center_Y` in metadata）
  - 保留降维结果（UMAP, PCA等）
  - 兼容Seurat V5和空转数据格式

- **`seurat_to_h5ad()`**: 将Seurat对象转换为h5ad格式（可选）

- **`batch_h5ad_to_seurat()`**: 批量转换函数

### 2. 关键特性
- ✅ 使用`.qs`格式保存（`qs::qsave`）
- ✅ 仅保留counts数据（移除data、scale.data等layer）
- ✅ 正确保留spatial坐标到metadata
- ✅ 兼容Seurat V5数据结构
- ✅ Python环境：spatial conda环境 (`~/miniconda3/envs/spatial/bin/python3.11`)
- ✅ R环境：RStudio R (`/opt/R/R-4.4.3/lib/R/bin/R`)

## 文件位置

### 主要脚本
- `/home/jacekyu/PCF/convert_h5ad_to_seurat.R` - 批量转换脚本
- `/home/jacekyu/PCF/test_single_sample.R` - 单样本测试脚本
- `/home/jacekyu/GCPM/GC&PM/0seurat<-->h5ad_optimized.Rmd` - R Markdown文档

### 数据路径
- 输入：`/home/jacekyu/PCF/spatial/*.annotated.h5ad`
- 输出：`/home/jacekyu/PCF/data/*.qs`

## 使用方法

### 测试单个样本
```bash
cd /home/jacekyu/PCF
export PATH="/opt/R/R-4.4.3/lib/R/bin:$PATH"
Rscript test_single_sample.R
```

### 批量转换
```bash
cd /home/jacekyu/PCF
export PATH="/opt/R/R-4.4.3/lib/R/bin:$PATH"
Rscript convert_h5ad_to_seurat.R
```

## 技术细节

### Spatial坐标处理
- h5ad中的`obsm["spatial"]`被提取并添加到Seurat对象的metadata中
- 列名：`Center_X`, `Center_Y`
- 如果metadata中已存在这些列，不会覆盖（保留原始数据）

### 数据格式
- 输入：AnnData (h5ad) - CODEX数据，包含21个marker基因
- 输出：Seurat V5对象，仅包含counts layer
- 保存格式：`.qs` (使用`qs::qsave`)

### 环境配置
- R路径已添加到`~/.bashrc`：`export PATH="/opt/R/R-4.4.3/lib/R/bin:$PATH"`
- Python环境在代码中自动配置为spatial conda环境

## 测试结果

✅ 测试样本：`20P.annotated.h5ad`
- 细胞数：74,444
- 基因数：21
- Spatial坐标：已正确保留（Center_X, Center_Y）
- 坐标范围：X: -5301.4 到 -114.49, Y: 7736.78 到 13016.05

## 注意事项

1. 确保已安装必要的R包：`Seurat`, `SeuratObject`, `reticulate`, `glue`, `qs`
2. 确保Python环境中有`scanpy`包
3. 转换后的对象仅包含counts数据，不包含其他layer
4. Spatial坐标存储在metadata中，不是作为降维结果

## 数据来源

h5ad文件由`build_spatialdata_batch.py`生成，包含：
- 21个marker基因的表达矩阵
- 细胞注释信息（major_lineage, cell_type_lvl1, cell_type_lvl2等）
- Spatial坐标（`obsm["spatial"]`）
- 各种状态标签（st_*）

