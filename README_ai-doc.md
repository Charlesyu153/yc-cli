# PCF 项目结构文档

## 项目概述

空间转录组学分析项目，专注于CAF Niche分析。
- 18个样本（P/QT/L/R类型）
- Python数据预处理 + R下游分析
- 支持h5ad/zarr/Seurat格式

## 快速导航

| 功能 | 路径 | 说明 |
|------|------|------|
| AI协作文档 | `ai-docs/` | 任务管理和协作规范 |
| Python预处理 | `PCF/*.py` | SpatialData构建和清理 |
| R分析脚本 | `script/*.R` | CAF Niche批量分析 |
| CAF Niche v2 | `script/caf_niche_v2/` | 最新分析流程 |
| 参考脚本 | `script_ref/*.Rmd` | CosMx分析参考 |
| 原始数据 | `rawdata/` | 18个样本原始数据 |
| 分析结果 | `data/caf_niche*/` | CAF Niche分析输出 |
| Seurat对象 | `data/*.qs` | R分析对象 |
| 空间数据 | `spatial/*.h5ad` | Python空间数据 |

## 目录结构

### 1. 文档管理

```
ai-docs/                                # AI协作文档管理系统
├── README.md, QUICKREF.md, USAGE.md    # 使用指南
├── current/                            # 进行中任务
├── archive/                            # 已完成任务
└── templates/                          # 模板和工具脚本
    ├── ai-guidelines.md                # AI协作规范
    ├── task-template.md                # 任务模板
    └── *.sh                            # 任务管理脚本

tasks/                                  # 历史任务（待迁移到ai-docs）
├── 0.TASKS.md
├── 1.cafnichetask.md
├── cafnichetaskV2.md
└── task_CNanalysis.md
```

### 2. 代码脚本

```
PCF/                                    # Python预处理脚本
├── 001prepare_csd(latest).py           # 数据准备
├── build_spatialdata_*.py              # SpatialData构建（批量/单样本）
├── clean_*_objects.py                  # 数据清理
├── check_missing_markers.py            # 标记检查
├── qc_summary.py                       # 质控汇总
└── *.md                                # 说明文档

script/                                 # R分析脚本
├── README_scripts.md                   # 脚本说明
├── caf_niche_batch*.R                  # CAF Niche批量分析
├── caf_niche_global_*.R                # 全局分析（比对/亚型空间）
├── caf_niche_niche_by_sample_grid.R    # 样本网格分析
├── caf_niche_v2/                       # CAF Niche v2流程（推荐）
├── *_annotation_coarse_*.R             # 细胞类型注释
├── summary_coarse_composition.R        # 组成汇总
└── Codex_*.Rmd                         # 单样本分析流程

script_ref/                             # CosMx参考流程（12个脚本）
├── CosMx_01_annotation_*.Rmd           # 注释流程
├── CosMx_02_cell_niche_*.Rmd           # Niche分析
├── CosMx_03_integration_*.{Rmd,ipynb}  # 整合分析
├── CosMx_04_visualization.Rmd          # 可视化
└── CosMx_05_cell_niche_across_*.Rmd    # 跨切片分析

utils/R/                                # R工具函数库
```

### 3. 数据存储

```
rawdata/                                # 原始数据（18个样本）
├── {1,2,3,6,7,8}P/                     # P类样本（原发灶）
├── {13,17,19,20,21}P/                  # P类样本
├── {6,7,21}R/                          # R类样本（复发）
├── {2,19}L/                            # L类样本（淋巴转移）
└── {13,19}QT/                          # QT类样本

rawdata1/                               # 原始数据备份
cleandata/                              # 清理后数据（同样本结构）

spatial/                                # Python空间数据
├── *.annotated.h5ad                    # AnnData格式（18个）
└── *.sdata.zarr/                       # SpatialData格式（18个）

data/                                   # R分析数据和结果
├── *.qs                                # Seurat对象（18个）
├── *.{geojson,tif}                     # 空间几何和影像
├── caf_niche/                          # CAF Niche v1结果
├── caf_niche_sectionnorm/              # 切片标准化结果
├── caf_niche_sectionnorm_myeloidrefined_redisco_autoK_stable_v2/
├── caf_niche_v2run_YYYYMMDD_HHMMSS/    # v2运行结果（时间戳）
├── coarse_composition/                 # 粗分类组成
├── qc_summary_{regions,samples}.tsv    # 质控汇总表
└── conversion_summary.qs               # 格式转换汇总

backup/step3_caf_niche/                 # 备份文件
```

### 4. 配置和工具

```
convert_h5ad_to_seurat.R                # h5ad -> Seurat转换
test_conversion.R                       # 转换测试
test_single_sample.R                    # 单样本测试
batch_process.log                       # 批处理日志
datapre.txt                             # 数据预处理说明
README_*.md                             # 各类说明文档
AGENTS.md                               # Agent配置
CONTEXT_SUMMARY.md                      # 上下文总结
```

## 数据处理流程

```
rawdata (18 samples)
    |
    v  [PCF/Python]
cleandata
    |
    v  [build_spatialdata_*.py]
spatial/*.{h5ad, sdata.zarr}
    |
    v  [convert_h5ad_to_seurat.R]
data/*.qs (Seurat objects)
    |
    v  [script/*.R]
data/caf_niche*/ (Analysis results)
```

## 样本信息

18个样本分类：
- P类（原发）：1P, 2P, 3P, 6P, 7P, 8P, 13P, 17P, 19P, 20P, 21P
- R类（复发）：6R, 7R, 21R
- L类（淋巴转移）：2L, 19L
- QT类：13QT, 19QT

配对样本：
- 6P-6R, 7P-7R, 21P-21R（原发-复发配对）
- 2P-2L, 19P-19L（原发-淋巴转移配对）
- 13P-13QT, 19P-19QT（原发-QT配对）

## AI协作规范

详见：`ai-docs/templates/ai-guidelines.md`

核心原则：
1. 无表情符号
2. 步骤完成即commit
3. 任务完成即归档
4. 代码简洁不超过500行

启动新任务：
```bash
cd ai-docs/templates
./new-task.sh "任务名称" "负责人"
```

归档已完成任务：
```bash
cd ai-docs/templates
./archive-task.sh "任务名称"
```

---

文档版本：v2.0
最后更新：2026-01-12