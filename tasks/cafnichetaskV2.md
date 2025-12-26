# CAF Niche Analysis Pipeline v2.0（auto-K + 全局 CAF subtype）

本文件描述一套可由 AI/脚本实现的 CAF niche 分析流程（Plan v2.0），包含输入输出、目录结构、关键步骤与指标定义。  
目标：在多样本空间数据上，得到 ≤5 个、跨样本可复现的 CAF 全局 subtype，并在每个样本上做空间可视化。

---

## 0. 数据与目录约定

### 0.1 输入样本

- 输入文件目录：`data/`
- 样本文件格式：`data/<sample>.qs`
- 从 `data/` 中选取样本：
  - 包含：所有 `*.qs`
  - 排除：
    - `13QT.qs`
    - `19QT.qs`
    - `conversion_summary.qs`
- 分组（用于“稀有 subtype”判断与合并）：
  - `P`（原发灶）：样本名以 `P` 结尾，例如 `13P.qs`
  - `M`（卵巢转移灶）：样本名以 `L` 或 `R` 结尾，例如 `19L.qs`、`7R.qs`
  - 为seurat对象加入group信息
  - 坐标单位（重要）：本项目现有 `.qs` 的 `spatial`/`Center_X/Y` 坐标在 **`coord_scale=1`** 时，
    以 `r=80` 得到的 CAF 邻域规模约为 O(10^2) 细胞（符合“niche 周围 <200 细胞”的预期）；若使用 `coord_scale=50` 且仍 `r=80`，
    邻域会退化到几乎只有自身（计数≈1），会导致 niche profile 失真。因此本 v2 默认按 `coord_scale=1` 运行。
- 每个 `<sample>.qs` 至少包含：
  - 细胞级坐标（如 `x`, `y`）
  - CAF 细胞的标记
  - `annotation_coarse`（粗 cell-type 注释，Myeloid 已细分并重命名为 `Macrophage/Neutrophil/Myeloid_other`）
  - （可选）其它元信息：`sample_id` 等

### 0.2 输出根目录

- 新输出根目录（本 pipeline 的所有结果）：  
  `data/caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2/`

所有后续步骤的输出都以此为根目录，按样本和步骤分子目录。

---

## 1. Step1：逐样本本地 CAF niche 聚类（auto-K，稳定性优先）

### 1.1 目标

对每个样本独立地：
- 在 CAF 细胞上做本地 niche clustering（local CAF clusters）
- 针对一组候选 k 值进行重复运行，评估：
  - 聚类稳定性（stability）
  - 拟合/解释力（fit）
- 自动选择一个最合适的 `selected_k`（稳定性优先，兼顾 fit）
- 输出该样本的：
  - 共识本地 cluster profile（供全局使用）
  - 本地 CAF cell → local_cluster 的赋值
  - auto-K 决策信息（k_selected.tsv）

### 1.2 候选 k 设置

- 对每个样本定义候选 k 集合，例如：
  - `k_list = {3, 4, 5, 6, 7, 8}`
- `k_list` 仅用于比较，**并不强制所有样本用同一个 k**。

### 1.3 重复运行设置（用于稳定性评估）

- 对每个样本、每个候选 k：
  - 运行 R 次聚类（例如 R = 5）
  - 每次运行应具有随机性：
    - 不同的随机 seed
    - 或/且 对 CAF 细胞进行子采样（例如随机保留 80% 细胞）
  - 每次运行输出：
    - `caf_clusters_k{k}_run{r}.tsv`
      - 内容：cell_id → local_cluster_id（基于该 run 的 clustering）
    - `niche_average_<profile>_k{k}_run{r}.tsv`
      - 内容：每个 local_cluster 的 profile
      - 列：coarse cell types（含 `Macrophage/Neutrophil/Myeloid_other`）的比例/频数

### 1.4 稳定性指标定义

对每个样本、每个 k，定义两个稳定性指标和一个综合稳定性分数。

#### 1.4.1 profile-level stability（`stab_profile(k)`）

- 输入：
  - 同一样本、同一 k 的多次运行的 `niche_average_<profile>_k{k}_run{r}.tsv`。
- 做法概述：
  1. 对任意两次运行 r1, r2：
     - 计算 cluster profile 相似度矩阵：每个 run 中的每个 cluster 与另一个 run 中每个 cluster 的相关性（如 Pearson）。
     - 采用匈牙利算法或其它最优化匹配算法，得到 cluster 间的最佳匹配对。
     - 对所有匹配对的相关性取平均，得到一次 (r1, r2) 的 profile 相似度。
  2. 对所有 (r1, r2) 组合的 profile 相似度取平均，得到 `stab_profile(k)`。
- 范围：
  - 将其 rescale 到 [0, 1] 区间，越大越稳定。

#### 1.4.2 assignment-level stability（`stab_assign(k)`）

- 输入：
  - 同一 k 的多次 run 输出的 `caf_clusters_k{k}_run{r}.tsv`（cell → cluster）。
- 做法概述：
  1. 对每对运行 r1, r2：
     - 只考虑在两个 run 中都被抽样到的 cell（若有子采样）。
     - 比较这些 cell 的 cluster label 一致性，计算一个聚类一致性指数，如：
       - ARI（Adjusted Rand Index）
       - 或 Fowlkes–Mallows Index
  2. 对所有 (r1, r2) 的指数取平均，得到 `stab_assign(k)`。
- 同样 rescale 到 [0, 1] 区间，越大越稳定。

#### 1.4.3 综合稳定性分数（`stab_combined(k)`）

- 定义：
  \[
  \text{stab\_combined}(k) = 0.7 \cdot \text{stab\_profile}(k) + 0.3 \cdot \text{stab\_assign}(k)
  \]
- 用于后续 auto-K 决策。

### 1.5 拟合/解释力指标 `fit_score(k)`

- 目的是区分：
  - k 太小 → 聚类过粗，拟合差
  - k 适中 → 聚类合理
  - 不追求 k 极大带来的微小拟合提升（避免过拟合）
- 可选定义（择一实现即可）：
  1. silhouette coefficient，基于 cluster profile 的空间；
  2. Calinski–Harabasz 指标；
  3. 基于 within-cluster variance 的简单评分。
- 将 `fit_score(k)` rescale 到 [0, 1] 区间，数值越大表示拟合越好。

### 1.6 auto-K 决策规则（逐样本）

对每个样本，在候选 k 中计算完 `stab_combined(k)` 和 `fit_score(k)` 后，按以下规则选 `selected_k`：

1. 找出该样本的最大稳定性：
   - `stab_max = max_k stab_combined(k)`
2. 定义“稳定候选集合”：
   - 给定阈值 δ（例如 0.05~0.1）
   - `K_stable = { k | stab_combined(k) >= stab_max - δ }`
3. 在 `K_stable` 中筛选拟合合理的 k：
   - 给定最小拟合阈值 `fit_min`（例如 0.6~0.7）
   - `K_good = { k ∈ K_stable | fit_score(k) >= fit_min }`
4. 如果 `K_good` 非空：
   - 选择 **最小的 k** 作为 `selected_k`（“最早稳定点”）。
5. 如果 `K_good` 为空（所有稳定 k 的拟合都低于阈值）：
   - fallback 逻辑：
     - 在所有 k 中选出 `stab_combined(k)` 排名前 2–3 的候选；
     - 在这几者中选 `fit_score(k)` 最大的一个作为 `selected_k`；
     - 同时标记该样本 `is_unstable = 1`，供后续全局加权/过滤。

### 1.7 共识 profile 与 cluster assignment（为后续步骤提供单一版本）

对每个样本、已选定的 `selected_k`：

1. **共识 profile：**
   - 收集该样本 `selected_k` 下的所有 run 的 cluster profiles；
   - 通过运行间 cluster 匹配（同 1.4.1 中使用的策略）对齐 cluster；
   - 对对应 cluster 的 profile 取平均，得到共识 `local_cluster_profile`。
   - 输出文件：
     - `niche_average_<profile>_selected.tsv`
       - 行：local_cluster_id
       - 列：coarse cell types（含 `Macrophage/Neutrophil/Myeloid_other`）的比例/频数

2. **共识 cluster assignment（可选但推荐）：**
   - 对 cell → cluster 的多次 run 结果做投票/最大后验：
     - 对每个 cell，统计其在 R 个 run 中最常见的 cluster label；
     - 作为最终的 `local_cluster_id`。
   - 输出文件（每样本）：
     - `caf_clusters_selected.tsv`
       - 内容：cell_id → consensus local_cluster_id

### 1.8 Step1 per-sample 结果：`k_selected.tsv`

每个样本输出一份 `k_selected.tsv`，至少包含以下字段：

- `sample_id`
- `selected_k`
- `stab_profile`（即 `stab_profile(selected_k)`）
- `stab_assign`（即 `stab_assign(selected_k)`）
- `stab_combined`（即 `stab_combined(selected_k)`）
- `fit_score`（即 `fit_score(selected_k)`）
- `is_unstable`（0 或 1）
- `profile_file`（例如 `niche_average_<profile>_selected.tsv` 的相对路径）
- （可选）`cluster_file`（例如 `caf_clusters_selected.tsv` 的相对路径）

---

## 2. Step1.5：构建全局输入矩阵与权重

### 2.1 汇总全样本本地 cluster profiles

遍历所有样本：

- 读取每个样本的：
  - `niche_average_<profile>_selected.tsv`
  - `k_selected.tsv`
  - CAF cluster summary（例如已有的 `caf_cluster_summary`，用于获取 `n_cells`）

构建如下两个全局对象：

1. **矩阵 X（local_cluster × coarse cell types）：**
   - 每一行对应一个“本地 cluster”，建议用组合 ID 标识：
     - `global_local_cluster_id = sample_id + ":" + local_cluster_id`
   - 每一个元素为该本地 cluster 在某 coarse cell-type（含 `Macrophage/Neutrophil/Myeloid_other`）上的比例/频数。

2. **metadata 表 M：**
   - 行与 X 一一对应，列包括：
     - `global_local_cluster_id`
     - `sample_id`
     - `local_cluster_id`
     - `selected_k`（该样本的 k）
     - `n_cells`（该 local cluster 中的 CAF cell 数）
     - `stab_combined`（从对应样本的 `k_selected.tsv` 继承）
     - `fit_score`
     - `is_unstable`（0 / 1）

### 2.2 定义全局权重 w_i

对每一行（local cluster）定义一个权重 `w_i` 用于全局聚类时加权：

\[
w_i = \log(1 + n\_cells) \times \text{stab\_combined}_i \times (1 - 0.5 \cdot is\_unstable_i)
\]

- `n_cells`：本地 cluster 的 CAF cell 数
- `stab_combined_i`：该 cluster 所在样本的稳定性分数
- `is_unstable_i` = 1 时，将权重打折（乘以 0.5）

---

## 3. Step2：全局 CAF subtype 聚类（K_global ≤ 5）

### 3.1 目标

- 在所有样本的 local cluster profiles（矩阵 X）上，得到 **≤5 个全局 CAF subtype**。
- 优先标准：
  - 每个 subtype 尽可能跨样本出现（覆盖样本多）
  - 减少“单样本独占”的 subtype
  - subtype 内部的 profile 相似度合理（类内离散度不要太大）
- 对“稀有/不稳定 subtype”的判定与合并：优先满足下述三种“组可复现”形态（见 3.3.2、3.5）：
  - Universal：同时覆盖 `P` 与 `NonP`
  - P-only：主要存在于 `P`
  - NonP-only：主要存在于 `NonP`

### 3.2 聚类方法

推荐使用带权重的层级聚类（或其它支持样本权重的聚类方法）：

- 输入：矩阵 X，权重 w_i。
- 距离度量：例如 `1 - 相关性` 或 Cosine 距离。
- 对需要评估的 K_global 取值范围：
  - `K_global_candidate = {2, 3, 4, 5}`

对于每个候选 K_global：

- 对聚类树进行 `cut` 得到 K_global 个 cluster（全局 CAF subtype）
- 计算以下指标以形成一个总目标函数。

### 3.3 全局目标函数 Score(K)

对每个 K_global，定义：

\[
\mathrm{Score}(K) = 
\alpha \cdot \mathrm{coverage}(K) 
- \beta \cdot \mathrm{group\_validity\_penalty}(K)
- \gamma \cdot \mathrm{within\_dispersion}(K)
\]

其中：

#### 3.3.1 coverage(K)

- 衡量 subtype 覆盖样本的广度。
- 对每个全局 subtype g：
  - 统计该 subtype 中涉及到的样本数 `S_g`；
  - 并可按权重 w_i 加权。
- 定义：
  \[
  \mathrm{coverage}(K) = \frac{1}{K} \sum_{g=1}^{K} \frac{S_g}{N_{\text{sample}}}
  \]
- 取值范围约在 [0,1]，越大表示越多 subtype 跨样本。

#### 3.3.2 group_validity_penalty(K)

- 惩罚“不符合组可复现形态”的 subtype，并用“软惩罚”引导 subtype 形成清晰的 P-only / NonP-only。
- 先对每个 subtype g 统计（均建议用权重 `w_i`，见 2.2）：
  - `cov_P(g)`：在 `P` 组中出现的样本占比（至少有非零权重的样本数 / `P` 样本总数）
  - `cov_NonP(g)`：在 `NonP` 组中出现的样本占比
  - `pur_P(g)`：该 subtype 的总权重中来自 `P` 组的占比（`pur_NonP(g)=1-pur_P(g)`）
- purity 口径（写死为加权口径）：
  \[
  \mathrm{pur\_P}(g) = \frac{\sum_{i \in g, \mathrm{sample}(i)\in P} w_i}{\sum_{i \in g} w_i}
  \]
- 规则参数（可调）：
  - 覆盖阈值：`COV_MIN = 0.5`（至少覆盖组内一半以上样本）
  - 纯度参考阈值（用于软惩罚，不是硬门槛）：`PUR_MIN = 0.7`
- 定义“形态标签”（至少满足其一，否则视为“不合法/需合并”）：
  - Universal：`cov_P(g) >= COV_MIN` 且 `cov_NonP(g) >= COV_MIN`
  - P-only-candidate：`cov_P(g) >= COV_MIN` 且 `cov_NonP(g) < COV_MIN`
  - NonP-only-candidate：`cov_NonP(g) >= COV_MIN` 且 `cov_P(g) < COV_MIN`
- 软惩罚项（连续）：
  - purity 缺口惩罚 `pen_purity(g)`：
    - 若 `g` 是 P-only-candidate：`max(0, PUR_MIN - pur_P(g))`
    - 若 `g` 是 NonP-only-candidate：`max(0, PUR_MIN - pur_NonP(g))`
    - 其它情况：`0`
  - leakage 惩罚 `pen_leak(g)`（防止“对侧权重强驱动但组内覆盖不足”的坏例子）：
    - 若 `g` 是 P-only-candidate：`max(0, pur_NonP(g) - (1 - PUR_MIN))`
    - 若 `g` 是 NonP-only-candidate：`max(0, pur_P(g) - (1 - PUR_MIN))`
    - 其它情况：`0`
- 定义整体惩罚（示例实现，包含一个硬有效性项 + 两个软惩罚项）：
  \[
  \mathrm{group\_validity\_penalty}(K) = \frac{1}{K} \sum_{g=1}^{K} \Big(\mathbf{1}\{\text{subtype } g \text{ is not Universal/P-only-candidate/NonP-only-candidate}\} + \lambda_p \cdot \mathrm{pen\_purity}(g) + \lambda_l \cdot \mathrm{pen\_leak}(g)\Big)
  \]
  - 建议取值：`λ_p = 1.0`，`λ_l = 1.0`（覆盖优先，但推动形态更清晰、减少边界抖动）

#### 3.3.3 within_dispersion(K)

- 衡量 subtype 内部的 profile 离散程度。
- 对每个 subtype g：
  - 计算 cluster 内所有 pairwise 距离的加权平均（权重 w_i）；
- 再对 g 取平均，得到 `within_dispersion(K)`。
- 越大表示类内差异越大。

#### 3.3.4 系数权重建议

- α = 1.0（覆盖样本优先）
- β = 1.0（强烈惩罚样本特异性）
- γ = 0.5（兼顾类内一致性）

### 3.4 选择 K_global

- 对 K = 2, 3, 4, 5 分别计算 `Score(K)`；
- 选取 `Score(K)` 最大的 K 作为 `K_global_selected`；
- 若多个 K 的得分接近，可在此基础上再结合生物学解释进行人工微调。

### 3.5 稀有/样本特异 subtype 合并策略

在已经选定 `K_global_selected` 后，对 subtype 进一步后处理：

1. 对每个 subtype g，计算 `cov_P/cov_NonP/pur_P/pen_purity/pen_leak`（同 3.3.2，purity 口径用权重 `w_i`）。
2. 判定 “需合并/不稳定” 的 subtype（满足任一条即可）：
   - 不满足 Universal / P-only-candidate / NonP-only-candidate 的任一形态（按 `COV_MIN=0.5` 判定）；
   - 或极端单样本独占（例如最大样本权重占比 `p_g > 0.8`）；
   - 或总权重（∑w_i）极低（例如处于全局权重分布的最下 5%）；
   - 或 `pen_purity/pen_leak` 很大（提示其“看起来属于某组但被对侧强驱动且覆盖不足”，优先合并以获得更清晰形态）。
3. 合并方式（重复迭代直到满足约束）：
   - 计算该 subtype centroid 与其它 subtype centroid 的相似度（相关性）；
   - 将其合并到“最相近”的 subtype（相关性最高/距离最小）；
   - 合并后重新计算 `cov_* / pur_* / pen_* / p_g`，直到所有 subtype 都能归入三类形态，且总数 ≤5。
4. 若合并导致类内离散度显著变差（或出现明显强分离被硬并），可回退并改为保留该 subtype；但必须保证最终 subtype 数目仍 ≤5。

---

## 4. Step2.5：全局 CAF subtype 反哺本地（半监督 refine）

### 4.1 目标

- 利用 Step2 得出的全局 CAF subtype，对每个样本的 CAF 细胞进行统一标注；
- 在保留每个样本本地信息的同时，提高跨样本可比性；
- 修正单个样本中由于局部噪声造成的不理想聚类。

### 4.2 输入

- 全局 subtype 的 centroid：
  - `C_global[g, ·]`，g=1..K_global（在 coarse cell-type 空间中）
- 每个样本的本地 cluster profile：
  - `L_i`：第 i 个 local cluster 的 profile
- 每个样本的 CAF 细胞 → local_cluster 的共识赋值：
  - `cell_id → local_cluster_id`

### 4.3 local cluster → global subtype 映射

对每个样本，逐个 local cluster：

1. 计算该 local cluster 与各 global subtype 的相似度：
   - `sim(i, g) = corr(L_i, C_global[g])`（或其它相似度）
2. 选择最大相似度的 global subtype 作为该 local cluster 的主要 subtype：
   - `g* = argmax_g sim(i, g)`
3. 得到映射：
   - `local_cluster_id i → global_subtype g*`
   - （必要时也可保留 soft assignment `sim(i, g)`）

### 4.4 CAF cell → global subtype 标注

对每个 CAF 细胞：

- 通过其 `local_cluster_id` 找到对应的 `global_subtype g*`；
- 把这个 `g*` 作为 CAF cell 的全局 subtype 标签；
- 输出每个样本的 cell-level 结果，例如：
  - `caf_cells_global_subtype.tsv`
    - 字段：`cell_id`, `sample_id`, `local_cluster_id`, `global_subtype`

### 4.5 可选局部 refine（如需要）

- 在每张切片内，可以基于全局 centroid 做一个小范围的局部优化：
  - 例如对每个 CAF cell，比较其局部表达/邻域信息与各 global subtype 的距离；
  - 对明显偏离当前 subtype 的 cell 进行重新归属。
- 此步为可选增强，不是必须。

---

## 5. Step3：展示与验收

### 5.1 质量控制与验收标准

#### 5.1.1 验收标准

- 全局 CAF subtype 数目：`K_global_final ≤ 5`
- 每个 subtype 尽可能跨样本出现：
  - 大部分 subtype 在多数样本中有非零权重
- Myeloid 细分（当前命名：`Macrophage/Neutrophil/Myeloid_other`）在：
  - local niche profiles
  - global subtype profiles
  中均有合理信号、具有可解释性。

#### 5.1.2 推荐 QC 图

1. **Subtype × Sample 权重热图**
   - 行：global subtype
   - 列：sample
   - 单元格：该 subtype 在该样本中对应的总权重（∑w_i）
   - 用于评估 subtype 的跨样本覆盖和样本特异性。
2. **Subtype 内 profile 热图**
   - 对每个 subtype 内的 local cluster profile 做聚合，可视化：
     - 各 coarse cell-type（含 `Macrophage/Neutrophil/Myeloid_other`）的平均比例；
   - 用于理解 CAF–免疫–其它 cell type 的组合模式。
3. **稳定性 vs 拟合度散点图**
   - 对每个样本的 `selected_k`，画点：
     - x 轴：`stab_combined`
     - y 轴：`fit_score`
   - 可以快速识别出“极不稳定”或“拟合很差”的样本。

### 5.2 空间展示（先输出示例确认后进行全量可视化）

基于 Step2.5 得到的 CAF cell → global_subtype 的标注：

- 对每个样本生成空间图（grid / image overlay 等）：
  - 每个 CAF cell 按 global_subtype 着色；
  - 结合原始组织切片（如 H&E）进行展示。
  - 参考脚本caf_niche_niche_by_sample_grid进行分样本展示。



---
