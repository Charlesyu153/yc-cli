明白。既然你的 **Seurat 对象就是由 `*_objects.tsv` 构建的**，那意味着：

* Seurat 内的细胞 ID、坐标（`Center X/Y`）、23 蛋白强度列、Positivity 列、以及 8 类细胞注释，理论上都已经在对象里了；
* 你现在新增的“公司交付的切片边界分割文件”（mask 或 polygon）主要用于两件事：

  1. 更高质量、可审稿的可视化（用真实细胞轮廓而非点）
  2. 更精确的“接触/邻接图”（cell-cell contact），从而增强 CN 耦合/边界混合等分析的可信度。

下面是更新后的规划（以“已有 Seurat”为起点）。

---

# 更新版 Plan（以 Seurat 为主数据源 + 外部边界文件增强）

## 0. 统一输入与命名规范（先把工程问题一次性定死）

对每张切片（每个 Seurat 对象）准备并记录：

* `sample_id`：从文件名或 Seurat 的 metadata 字段解析（包含 P / L / R 信息）
* `group`：`Primary`（P） vs `Met`（L/R）
* `cell_id`：Seurat 的细胞名（`Cells(seu)`），它应当与 objects.tsv 的 `Name` 一致
* `x,y`：Seurat 中的细胞坐标（来源于 `Center X/Center Y`）

**验收标准：**

* Seurat 中每个细胞都有唯一 `cell_id`，且坐标非空。
* `group` 能对每张切片唯一确定。

---

## 1. 导入并对齐“边界分割文件”（mask 或 polygon）

你的 Seurat 已经有质心坐标，因此这一步的核心不是“重建细胞”，而是把“边界几何”正确挂到同一个 `cell_id` 上。

### 1.1 如果交付的是 polygon（优先）

* 读取 polygon 文件，确保每个 polygon 有一个明确的 `cell_id` 字段（或可映射字段）。
* 将 polygon 表整理成：
  `cell_id`, `geometry`（多边形），可选 `area`, `perimeter` 等派生几何属性。

### 1.2 如果交付的是 label mask

* 读取 mask（整型图像，像素值=cell_id）。
* 将每个 cell_id 的像素区域矢量化成 polygon（可同时保存 mask 便于追溯）。
* 输出同样的 polygon 表：`cell_id`, `geometry`（以及几何属性）。

### 1.3 坐标系对齐（必须做一次抽样验证）

对每张切片随机抽 100 个细胞，验证：

* polygon 的 centroid 与 Seurat 的 `(x,y)` 是否重合（允许小误差）
* 若不重合，修正：缩放（像素↔微米）、平移（offset）、或 Y 轴翻转

**验收标准：**

* 修正后，多数抽样细胞的 centroid 与 Seurat 坐标误差应远小于典型细胞直径（例如 < 10% 细胞直径）。

---

## 2. 形成“主细胞表”（Master Cell Table）并做几何 QC

以 Seurat 为主，外部 polygon 作为附加字段。

### 2.1 主表字段（每个细胞一行）

从 Seurat 取：

* `cell_id`
* `sample_id`, `group`
* `x,y`
* `cell_type`（8 类）
* 23 个 marker 的连续强度列
* Positivity 列（0/1）

从 polygon 取（若有）：

* `geometry`（细胞边界）
* `area`, `perimeter`, `aspect_ratio`（建议计算）

### 2.2 几何 QC（只影响“接触/可视化”，不强制影响 CN kNN）

* 标记并可选剔除：面积极小碎片、面积极大合并体、异常细长形状、自交 polygon
* 记录每张切片的 QC 剔除比例作为质量指标（用于解释下游差异是否受分割质量影响）

---

## 3. 蛋白表达标准化（仅用于功能比较，不用于 CN 定义）

CN 定义使用“细胞类型组成”即可；蛋白用于后续解释与机制。

### 3.1 连续强度做 arcsinh 变换

对 23 个 marker 强度列统一做：

* `marker_asinh = asinh(marker / cofactor)`
  cofactor 先固定一个（例如 5 或 10），全队列一致。

### 3.2 Positivity 列不做变换

Positivity 用于“阳性比例/状态标签”，不做标准化。

### 3.3样本内鲁棒缩放用于可视化

如果你要做跨切片热图或聚类展示，可以对每张切片每个 marker 的 asinh 值做 robust z-score。
但用于组间推断时，建议以“样本/患者汇总 + 模型控制”替代强行拉齐分布。

---

## 4. CN（Cellular Neighborhood）定义：以 Schürch 的 kNN-window 为主线

你现在的优势是：Seurat 已经提供了所有细胞的坐标和细胞类型，CN 可以直接算。

### 4.1 计算 kNN window（每张切片内）

* 对每个细胞，用 `(x,y)` 找 **k=10** 个最近邻（建议包含自身，保持一致即可）。
* 得到每个细胞一个“局部窗口”。

### 4.2 计算窗口细胞组成向量

* 在窗口内统计 8 类 cell_type 的比例，得到长度为 8 的组成向量（和=1）。

### 4.3 聚类得到 CN

* 对所有细胞的组成向量做 KMeans/MiniBatchKMeans 聚类。
* CN 数量 K：建议从 8–12 扫描，选择标准：

  * 每个 CN 占比不至于过小（例如 >1% 或 >N 细胞）
  * 跨抽样重复聚类稳定（ARI/NMI 高）
  * CN 组成具有可解释结构（如 Tumor-enriched、CAF-stroma、Immune、Myeloid、Endothelial/vascular、Neutrophil-hot 等）

### 4.4 将 CN 标签写回 Seurat

* 把每个细胞的 CN 作为 `seu@meta.data$CN`，后续所有比较都在 Seurat 框架里完成。

---

## 5. 利用“真实边界”增强空间相互作用：接触图与边界混合

这部分是你拿到边界文件后，最值得做的“增值项”。

### 5.1 构建 cell-cell 接触图（优先用 polygon）

* 定义“接触”：两细胞 polygon 共享边界，或 polygon 间最短距离小于阈值（阈值可设为 0 或非常小）。
* 得到 adjacency graph（边表示接触）。

### 5.2 计算三类指标（用于 P vs L/R 对比）

1. **CN-CN 接触富集**：CN 之间接触频率相对随机置换（shuffle CN 或 shuffle cell_type）的富集倍数
2. **边界混合指数**：一个细胞的邻居中“不同 CN”的比例；按 CN 汇总
3. **“耦合/重路由”**：每个 CN 的主要接触对象在 P 与 L/R 是否变化（例如 Immune CN 原本主要接触 Boundary CN，转移灶变成主要接触 Myeloid CN）

---

## 6. 组间比较（Primary P vs Met L/R）

分析单位建议优先按“切片/患者”汇总，避免单细胞数量差导致虚假显著。

### 6.1 CN 频率差异

* 每张切片计算每个 CN 的占比
* 比较 P vs L/R（若能解析 patient_id，优先用配对或混合效应模型）

### 6.2 CN 内部组成重排

* 在同一 CN 内比较 8 类细胞比例的变化（P vs L/R）
* 用于解释“CN 频率变化背后是什么细胞在驱动”。

### 6.3 CN 分层下的功能蛋白差异

两种输出都做，便于形成论文式论证：

* “整体差异”：不分 CN，比较 P vs L/R
* “CN-特异差异”：在每个 CN 内比较 P vs L/R
  重点关注你面板里典型功能轴：PD-1 / PD-L1 / TIGIT / CTLA-4 / Ki67 / GZMB 以及 P5CS/PYCR1 等。

同时输出两种度量：

* 连续表达（asinh 强度）
* 阳性比例（Positivity）

---

## 7. 最终交付物（建议你按这个顺序做图和写结果）

1. 每张切片：用 polygon 按 `cell_type` 着色的组织图（代表性样本即可）
2. 每张切片：用 polygon 按 `CN` 着色的组织图（展示 CN 的空间“域”）
3. CN × cell_type 的组成热图（用来命名 CN）
4. CN 频率在 P vs L/R 的差异图（box/violin + 效应量）
5. CN-CN 接触网络或差异热图（体现耦合与重路由）
6. CN 分层的功能 marker 差异（最好给出“整体不显著但 CN-特异显著”的案例）

---

## 你现在只需要补充一项信息，就能把计划完全落地

公司交付的“边界分割文件”具体是哪一种格式：

* label mask（TIFF/OME-TIFF）
  还是
* polygon（GeoJSON/Parquet/CSV 顶点）

你告诉我格式后，我可以把第 1 步“对齐与验收”部分进一步具体化成非常明确的执行清单（包括：如何确认 cell_id 对应、如何判断需要缩放/翻转、以及如何定义接触阈值）。
