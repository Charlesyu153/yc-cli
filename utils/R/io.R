slim_save <- function(seu, file, layers = c("counts"),
                      slots = c("pca", "umap", "spatial")) {
  # 构建一个新的精简版 Seurat 对象
  seu.export <- CreateSeuratObject(
    counts = LayerData(seu, layer = layers[1]),
    meta.data = seu@meta.data
  )
  
  # 复制指定的 slots
  for (slot in slots) {
    if (slot %in% names(seu)) {
      seu.export[[slot]] <- seu[[slot]]
    } else {
      warning(paste("Slot", slot, "not found in input Seurat object."))
    }
  }
  
  # 保存为 .qs 文件
  qs::qsave(seu.export, file)
  message("Saved slim Seurat object to: ", file)
  
  invisible(seu.export)
}
