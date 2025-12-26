# script/ 脚本分类（快速索引）

本目录主要放“可直接运行的分析脚本”和少量 Rmd。下面按用途做一个最小可用的索引：每个条目给出**做什么**、**主要输入**、**主要输出**（输出路径以脚本 defaults 为准，可用参数覆盖）。

## CAF niche（legacy：固定 k 的旧流程）

- `script/caf_niche_batch.R`
  - 作用：单/多样本 CAF niche（未做 section-normalization 的旧版本；含 r 扫描、NMF、聚类与基础可视化）。
  - 输入：`data/<sample>.qs`
  - 输出：脚本内 `output_dir` 指定的目录（按样本分子目录）。

- `script/caf_niche_batch_sectionnorm.R`
  - 作用：逐样本 CAF niche（section-normalized），固定 `k_local`（默认 5）；输出 per-sample niche 矩阵、聚类与图。
  - 输入：`data/<sample>.qs`
  - 输出：`data/caf_niche_sectionnorm/<sample>/`（例如 `caf_clusters_k5.tsv`、`niche_celltype_count.tsv`、`niche_celltype_fraction.tsv`、`niche_celltype_corrected.tsv` 等）。

- `script/caf_niche_global_alignment.R`
  - 作用：把各样本的 local CAF cluster profiles 做全局对齐/聚类（默认用 `corrected`），生成全局 subtype 映射与相似度热图。
  - 输入：`data/caf_niche_sectionnorm/<sample>/` 的 `niche_average_*` 与 `caf_cluster_summary_*` 等（固定 k）。
  - 输出：`data/caf_niche_sectionnorm/_global*/k5/corrected/`（例如 `global_subtype_mapping.tsv`、`global_similarity_heatmap.pdf`）。

- `script/caf_niche_global_subtype_spatial.R`
  - 作用：把全局 subtype 映射回空间坐标并输出空间可视化（配合 global alignment 的输出）。
  - 输入：`data/caf_niche_sectionnorm/` + `_global*/k*/.../global_subtype_mapping.tsv`
  - 输出：脚本 `output_dir` 下的空间图/汇总文件。

- `script/caf_niche_niche_by_sample_grid.R`
  - 作用：niche-by-sample 网格图（旧流程；脚本当前默认做 “2-sample example”），左侧 CAFsubtype，中央邻域细胞，右侧 composition。
  - 输入：`data/caf_niche_sectionnorm/<sample>/` + `_global*/.../global_subtype_mapping.tsv` + `data/<sample>.qs`
  - 输出：`.../niche_by_sample/example_niche_by_sample_grid*.pdf`

## 注释/数据更新

- `script/update_annotation_coarse_refine_myeloid.R`
  - 作用：备份原 `annotation_coarse` 到 `annotation_coarse_old`，并在 `annotation_coarse` 内细分 Myeloid（从 fine 注释列推断）。
  - 输入：`data/*.qs`
  - 输出：默认 in-place 覆盖写回 `data/*.qs`（也支持写到新目录）。

- `script/rename_annotation_coarse_myeloid_labels.R`
  - 作用：重命名 coarse 中的 myeloid 细分类标签（例如 `Myeloid_Macrophage -> Macrophage` 等）。
  - 输入：`data/*.qs`
  - 输出：默认 in-place 覆盖写回 `data/*.qs`（也支持写到新目录）。

## QC/汇总

- `script/summary_coarse_composition.R`
  - 作用：汇总每个样本的 `annotation_coarse` 组成（count + fraction）并输出堆叠柱状图。
  - 输入：`data/*.qs`
  - 输出：默认 `data/coarse_composition/`（`coarse_counts.tsv`、`coarse_fractions.tsv`、`*.pdf`）。

## v2 pipeline（本次“auto-K + 全局 subtype”重新发现）

本次任务相关的**可复现脚本副本**、以及所有输出都放在：

- `data/caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2/scripts/`
- `data/caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2/step1/`、`step2/`、`step3/`

