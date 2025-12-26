# CAF Niche Log2 Enrichment (Counts-Based) / CAF 邻域 Log2 富集（基于计数）

This document describes the log2 enrichment algorithm used by
`script/caf_niche_batch_sectionnorm.R` and the support metrics written
alongside the CAF niche outputs.  
本文档说明 `script/caf_niche_batch_sectionnorm.R` 使用的 log2 富集算法，
以及输出的支撑性统计。

## Summary / 总览

We replace simple fold-change correction with **count-based log2 enrichment**
and a pseudocount to stabilize rare types (e.g., B_cell) and make effect
sizes interpretable.  
我们用**计数版 log2 富集**替代简单倍数校正，并加入伪计数，稳定稀有类型的
估计并提升可解释性。

- `log2_enrich = 0` means observed ~= expected / 观察≈期望
- `log2_enrich = 1` means ~2x expected / 约 2× 期望
- `log2_enrich = -1` means ~0.5x expected / 约 0.5× 期望

## Algorithm / 算法

For each sample and each CAF cell (query cell), define a circular neighborhood
with radius `r` (microns). For cell type `t` in neighborhood `i`:

```
x_i,t = count of type t in neighborhood i
N_i   = total neighbor count in neighborhood i
```

Section (whole-slice) composition:

```
section_prop_t = (number of cells of type t in slice) / (total cells in slice)
```

Expected counts per neighborhood:

```
E_i,t = N_i * section_prop_t
```

Log2 enrichment with pseudocount `alpha`:

```
log2_enrich_i,t = log2((x_i,t + alpha) / (E_i,t + alpha))
```

Default `alpha = 0.5` (Jeffreys-style smoothing), overridable by
`--log_enrich_alpha`.

### 公式细节 / Formula details

Let `C` be all cells in a slice, `I(condition)` an indicator function:

```
x_i,t = sum_{c in C} I(c in neighborhood i AND type(c)=t)
N_i   = sum_t x_i,t
section_prop_t = (sum_{c in C} I(type(c)=t)) / |C|
E_i,t = N_i * section_prop_t
```

The log2 enrichment uses base-2 logarithm. `alpha` can be 0.5 (Jeffreys)
or 1.0 (Laplace) depending on desired smoothing.

### Why this approach / 为什么这样算

- **Stability for rare types / 稳定稀有类型**:
  When `E` is small, pseudocount prevents division by ~0 and huge ratios.
- **Interpretability / 可解释性**:
  Log2 scale gives symmetric enrichment/depletion.
- **Counts-based expectation / 基于计数的期望**:
  Uses neighborhood size `N_i`, so large neighborhoods are not over-penalized.

## Outputs (per sample) / 输出（每样本）

Raw counts and composition / 原始计数与组成:

- `niche_celltype_count.tsv` (CAF x celltype): `x_i,t`
- `niche_celltype_fraction.tsv` (CAF x celltype): `x_i,t / N_i`
- `section_celltype_fraction.tsv` (celltype): `section_prop_t`

Expected counts and log2 enrichment / 期望计数与 log2 富集:

- `niche_celltype_expected.tsv` (CAF x celltype): `E_i,t`
- `niche_celltype_corrected.tsv` (CAF x celltype): `log2_enrich_i,t`

Support statistics (effect size + evidence) / 支撑性统计（效应量 + 证据）:

- `niche_neighbor_total.tsv`: `N_i` per CAF
- `niche_support_summary.tsv`:
  - `neighbor_total`: mean / p10 / median / p90 for `N_i`
  - for each cell type `t`: mean / p10 / median / p90 for `x_i,t` and `E_i,t`

## NMF and clustering / NMF 与聚类

When `--use_section_normalization true`, NMF uses `log2_enrich`.
Because log2 enrichment can be negative, the script shifts the matrix
by `-min_value` to keep it non-negative before NMF, and logs the shift.
  
当启用 `--use_section_normalization true` 时，NMF 使用 log2_enrich。
由于 log2 富集可能为负值，脚本会按最小值平移为非负并记录平移量。

## Global alignment / 全局对齐

Global alignment uses the corrected profile matrix directly. Since the
corrected matrix is already log2 enrichment, the alignment script
auto-sets `transform = none` for `profile=corrected`.  
全局对齐直接使用 log2 enrichment，因此 `profile=corrected` 时会自动将
`transform` 设为 `none`，避免重复取 log。

## Notes / 注意事项

- Interpret enrichment together with `N_i`, `x_i,t`, and `E_i,t`. A high
  log2 enrichment with tiny `N_i` is less reliable than the same effect
  supported by large `N_i`.
- You can change `alpha` (e.g., 1.0) via `--log_enrich_alpha`.
