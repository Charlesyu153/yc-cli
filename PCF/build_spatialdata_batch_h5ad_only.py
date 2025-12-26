#!/usr/bin/env python3
"""
批量处理脚本：遍历 cleandata 下所有样本，
构建 AnnData 对象，并基于 marker 规则进行细胞谱系与状态注释。

用法：
    python build_spatialdata_batch_h5ad_only.py

输出：
    - spatial/<sample_id>.annotated.h5ad
    - 不生成 .sdata.zarr
"""

from __future__ import annotations

import logging
import os
import warnings
import re
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import anndata as ad

# ============================================================================
# 配置参数
# ============================================================================

ROOT = Path("/home/jacekyu/PCF")
CLEANDATA_ROOT = ROOT / "cleandata"
OUT_DIR = ROOT / "spatial"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 标准 marker 列表
CANONICAL_MARKERS = [
    "CD163", "CD20", "CD31", "CD3e", "CD4", "CD66", "CD68", "CD8", "CD86",
    "CTLA-4", "CLDN18.2", "FOXP3", "TIGIT", "GZMB", "Ki67", "P5CS",
    "PD-1", "PD-L1", "PYCR1", "PanCK", "SMA", "Vimentin", "DAPI"
]

# Positivity 列前缀
POSITIVITY_PREFIX = "Positivity-"

# 坐标列名
COORD_X = "Center X"
COORD_Y = "Center Y"

# Negative/Artifact 判断：所有 marker 强度总和阈值（可调参数）
# 若细胞所有 marker 强度加和 < 此阈值，则视为 Negative_or_Artifact
NEGATIVE_THRESHOLD = 1.0

# Positivity 缺失时的兜底策略（默认关闭，避免引入未经校准的阈值判断）
# - PCF_INFER_POSITIVITY_FROM_INTENSITY=1: 启用强度阈值推断
# - PCF_INTENSITY_POSITIVITY_THRESHOLD: 阈值（默认 0.0）
INFER_POSITIVITY_FROM_INTENSITY = os.environ.get("PCF_INFER_POSITIVITY_FROM_INTENSITY") == "1"
INTENSITY_POSITIVITY_THRESHOLD = float(os.environ.get("PCF_INTENSITY_POSITIVITY_THRESHOLD", "0.0"))

# 设置日志
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# ============================================================================
# 数据加载
# ============================================================================

def load_sample_objects(path: Path) -> pd.DataFrame:
    """
    读取样本的 objects.tsv 文件。
    
    Parameters
    ----------
    path : Path
        tsv 文件路径
    
    Returns
    -------
    DataFrame
        原始数据 DataFrame
    """
    if not path.exists():
        raise FileNotFoundError(f"文件不存在: {path}")
    
    df = pd.read_csv(path, sep="\t")
    logger.info(f"已加载 {len(df)} 个细胞，共 {len(df.columns)} 列")
    return df


def find_marker_columns(df: pd.DataFrame) -> dict:
    """
    在 DataFrame 中查找 marker 强度列和 positivity 列。
    
    Returns
    -------
    dict
        {
            "intensity": {marker_name: column_name, ...},
            "positivity": {marker_name: column_name, ...},
            "missing_intensity": [marker_name, ...],
            "missing_positivity": [marker_name, ...]
        }
    """
    result = {
        "intensity": {},
        "positivity": {},
        "missing_intensity": [],
        "missing_positivity": []
    }
    
    cols = set(df.columns)
    
    for marker in CANONICAL_MARKERS:
        # 查找强度列（直接 marker 名）
        if marker in cols:
            result["intensity"][marker] = marker
        else:
            result["missing_intensity"].append(marker)
        
        # 查找 positivity 列
        pos_col = f"{POSITIVITY_PREFIX}{marker}"
        if pos_col in cols:
            result["positivity"][marker] = pos_col
        else:
            result["missing_positivity"].append(marker)
    
    # 打印警告
    if result["missing_intensity"]:
        logger.warning(f"缺少强度列: {result['missing_intensity']}")
    if result["missing_positivity"]:
        logger.warning(f"缺少 positivity 列: {result['missing_positivity']}")
    
    return result


# ============================================================================
# Marker Positivity 推断（若无 positivity 列，可用阈值法推断）
# ============================================================================

def infer_marker_positivity(df: pd.DataFrame, marker_config: dict) -> pd.DataFrame:
    """
    为 DataFrame 添加 marker 的布尔阳性列（bool_xxx 格式）。
    优先使用已有的 Positivity 列，否则可用强度列 + 阈值推断。
    
    Parameters
    ----------
    df : DataFrame
        输入数据
    marker_config : dict
        find_marker_columns 的返回结果
    
    Returns
    -------
    DataFrame
        添加了 bool_<marker> 列的 DataFrame
    """
    df = df.copy()
    
    for marker in CANONICAL_MARKERS:
        bool_col = f"bool_{marker}"
        
        # 优先使用 positivity 列
        pos_col = marker_config["positivity"].get(marker)
        if pos_col and pos_col in df.columns:
            # Positivity 列是 0/1 整数，转换为布尔值
            df[bool_col] = (df[pos_col] == 1)
        else:
            intensity_col = marker_config["intensity"].get(marker)
            if INFER_POSITIVITY_FROM_INTENSITY and intensity_col and intensity_col in df.columns:
                # 基于强度阈值推断（可选）
                df[bool_col] = (pd.to_numeric(df[intensity_col], errors="coerce").fillna(0.0) > INTENSITY_POSITIVITY_THRESHOLD)
                logger.warning(
                    f"Marker {marker} 缺少 {POSITIVITY_PREFIX}{marker}，已用强度阈值推断："
                    f"{intensity_col} > {INTENSITY_POSITIVITY_THRESHOLD}"
                )
            else:
                # 保守：无法确定则默认 False
                df[bool_col] = False
                if intensity_col and intensity_col in df.columns:
                    logger.warning(
                        f"Marker {marker} 缺少 {POSITIVITY_PREFIX}{marker}，已默认 False。"
                        f"（可设 PCF_INFER_POSITIVITY_FROM_INTENSITY=1 启用强度阈值推断）"
                    )
    
    return df


# ============================================================================
# 细胞注释核心逻辑
# ============================================================================

def is_positive(row: pd.Series, marker: str) -> bool:
    """检查细胞某个 marker 是否阳性。"""
    return bool(row.get(f"bool_{marker}", False))


def _set_annotation_fields(result: dict, coarse: str, fine: str) -> None:
    """设置 annotation/annotation_coarse/annotation_fine 字段（annotation 默认使用 coarse）。"""
    result["annotation_coarse"] = coarse
    result["annotation_fine"] = fine
    result["annotation"] = coarse


def annotate_cell(row: pd.Series, marker_intensity_cols: list, vimentin_observable: bool) -> dict:
    """
    根据 marker 规则对单个细胞进行谱系与状态注释。
    
    Parameters
    ----------
    row : Series
        单个细胞的行数据（包含 bool_<marker> 列）
    marker_intensity_cols : list
        marker 强度列名列表，用于计算总强度
    vimentin_observable : bool
        Vimentin 是否在该数据中“可观测”（存在 positivity 列，或在启用强度阈值推断时存在强度列）。
        若不可观测，则不应使用 pos("Vimentin") 作为硬性门控条件。
    
    Returns
    -------
    dict
        包含 major_lineage, cell_type_lvl1, cell_type_lvl2, annotation_coarse/fine 及各 st_* 状态标签
    """
    result = {
        "major_lineage": "Unknown",
        "cell_type_lvl1": "Unknown",
        "cell_type_lvl2": "Unknown",
        "annotation": "Unknown",
        "annotation_coarse": "Unknown",
        "annotation_fine": "Unknown",
    }
    
    # 初始化所有状态标签为 False
    state_labels = [
        "st_proliferating", "st_checkpoint_PD1", "st_checkpoint_PDL1",
        "st_checkpoint_CTLA4", "st_checkpoint_TIGIT",
        "st_metabolic_P5CS_high", "st_metabolic_PYCR1_high",
        "st_Treg", "st_CD8_effector", "st_CD8_exhausted", "st_CD4_exhausted",
        "st_B_proliferating", "st_B_exhausted",
        "st_M1_like", "st_M2_like", "st_M_mixed", "st_M_exhausted", "st_M_proliferating",
        "st_Neutrophil_activated", "st_Neutrophil_checkpoint",
        "st_CAF_PD-L1_pos", "st_CAF_metabolic_P5CS", "st_CAF_metabolic_PYCR1"
    ]
    for st in state_labels:
        result[st] = False
    
    # ========================================================================
    # 1. Negative / Artifact 检查
    # ========================================================================
    total_intensity = 0.0
    for col in marker_intensity_cols:
        val = row.get(col, 0)
        if pd.notna(val):
            total_intensity += float(val)
    
    if total_intensity < NEGATIVE_THRESHOLD:
        result["major_lineage"] = "Negative_or_Artifact"
        result["cell_type_lvl1"] = "Negative_or_Artifact"
        result["cell_type_lvl2"] = "Negative_or_Artifact"
        _set_annotation_fields(result, "Negative_or_Artifact", "Negative_or_Artifact")
        return result
    
    # 快捷函数
    def pos(marker: str) -> bool:
        return is_positive(row, marker)
    
    # ========================================================================
    # 2. Endothelial（血管内皮）: CD31+ → 归入 Stromal
    # ========================================================================
    if pos("CD31"):
        result["major_lineage"] = "Stromal"
        result["cell_type_lvl1"] = "Endothelial"
        if pos("SMA"):
            result["cell_type_lvl2"] = "Vessel_mature"
        else:
            result["cell_type_lvl2"] = "Vessel_immature"
        _set_annotation_fields(result, "Endothelial", "Endothelial")
        # 添加通用状态标签后返回
        _add_common_states(result, row, pos)
        return result
    
    # ========================================================================
    # 3. Tumor / Epithelial: PanCK+ AND NOT CD31+
    # ========================================================================
    if pos("PanCK") and not pos("CD31"):
        result["major_lineage"] = "Tumor/Epithelial"
        result["cell_type_lvl1"] = "Tumor_epithelial"
        
        # CLDN18.2 分层
        if pos("CLDN18.2"):
            prefix = "Tumor_CLDN18p"
        else:
            prefix = "Tumor_CLDN18n"
        
        # P5CS 分层（用于 cell_type_lvl2）
        if pos("P5CS"):
            result["cell_type_lvl2"] = f"{prefix}_P5CSp"
        else:
            result["cell_type_lvl2"] = f"{prefix}_P5CSn"

        _set_annotation_fields(result, "Tumor", result["cell_type_lvl2"])
        
        _add_common_states(result, row, pos)
        return result
    
    # ========================================================================
    # 4. Immune（免疫谱系）: PanCK- AND CD31- 且满足免疫 marker 条件
    # ========================================================================
    is_panck_neg = not pos("PanCK")
    is_cd31_neg = not pos("CD31")
    
    has_immune_marker = (
        pos("CD3e") or pos("CD20") or 
        pos("CD68") or pos("CD163") or pos("CD86") or 
        pos("CD66")
    )
    
    if is_panck_neg and is_cd31_neg and has_immune_marker:
        result["major_lineage"] = "Immune"
        
        # ====================================================================
        # 4.1 T 细胞: CD3e+
        # ====================================================================
        if pos("CD3e"):
            result["cell_type_lvl1"] = "T_cell"
            
            cd4_pos = pos("CD4")
            cd8_pos = pos("CD8")
            
            if cd4_pos and not cd8_pos:
                result["cell_type_lvl1"] = "CD4_T"
            elif cd8_pos and not cd4_pos:
                result["cell_type_lvl1"] = "CD8_T"
            elif cd4_pos and cd8_pos:
                result["cell_type_lvl1"] = "CD4_CD8_DP_T"
            else:
                result["cell_type_lvl1"] = "T_DN_or_gamma_delta"
            
            result["cell_type_lvl2"] = result["cell_type_lvl1"]
            _set_annotation_fields(result, "T_cell", result["cell_type_lvl1"])
            
            # T 细胞特有状态
            if result["cell_type_lvl1"] == "CD4_T":
                if pos("FOXP3") or pos("CTLA-4"):
                    result["st_Treg"] = True
                if pos("PD-1") or pos("PD-L1") or pos("TIGIT"):
                    result["st_CD4_exhausted"] = True
            
            if result["cell_type_lvl1"] == "CD8_T":
                if pos("GZMB"):
                    result["st_CD8_effector"] = True
                if pos("PD-1") or pos("PD-L1") or pos("TIGIT"):
                    result["st_CD8_exhausted"] = True
            
            _add_common_states(result, row, pos)
            return result
        
        # ====================================================================
        # 4.2 B 细胞: CD20+ 且 CD3e-
        # ====================================================================
        if pos("CD20") and not pos("CD3e"):
            result["cell_type_lvl1"] = "B_cell"
            result["cell_type_lvl2"] = "B_cell"
            _set_annotation_fields(result, "B_cell", "B_cell")
            
            # B 细胞状态
            if pos("Ki67"):
                result["st_B_proliferating"] = True
            if pos("PD-1") or pos("PD-L1"):
                result["st_B_exhausted"] = True
            
            _add_common_states(result, row, pos)
            return result
        
        # ====================================================================
        # 4.3 中性粒细胞: CD66+ AND CD3e- AND CD20- AND CD68-
        # ====================================================================
        if pos("CD66") and not pos("CD3e") and not pos("CD20") and not pos("CD68"):
            result["cell_type_lvl1"] = "Neutrophil"
            result["cell_type_lvl2"] = "Neutrophil"
            _set_annotation_fields(result, "Myeloid", "Neutrophil")
            
            # Neutrophil 状态
            if pos("GZMB") or pos("Ki67"):
                result["st_Neutrophil_activated"] = True
            if pos("PD-1") or pos("PD-L1"):
                result["st_Neutrophil_checkpoint"] = True
            
            _add_common_states(result, row, pos)
            return result
        
        # ====================================================================
        # 4.4 单核-巨噬 / 髓系: CD3e- AND CD20- AND CD66- AND (CD68+ OR CD163+ OR CD86+)
        # ====================================================================
        if (not pos("CD3e") and not pos("CD20") and not pos("CD66") and
            (pos("CD68") or pos("CD163") or pos("CD86"))):
            result["cell_type_lvl1"] = "Monocyte_macrophage"
            
            cd163_pos = pos("CD163")
            cd86_pos = pos("CD86")
            cd68_pos = pos("CD68")
            
            if cd163_pos and not cd86_pos:
                result["cell_type_lvl2"] = "Macrophage_M2_like"
                result["st_M2_like"] = True
            elif cd86_pos and not cd163_pos:
                result["cell_type_lvl2"] = "Macrophage_M1_like"
                result["st_M1_like"] = True
            elif cd163_pos and cd86_pos:
                result["cell_type_lvl2"] = "Macrophage_M_mixed"
                result["st_M_mixed"] = True
            elif cd68_pos and not cd163_pos and not cd86_pos:
                result["cell_type_lvl2"] = "Macrophage_unspecified"
            else:
                result["cell_type_lvl2"] = "Monocyte_macrophage_other"
            
            # 髓系状态
            if pos("PD-1") or pos("PD-L1"):
                result["st_M_exhausted"] = True
            if pos("Ki67"):
                result["st_M_proliferating"] = True
            
            _set_annotation_fields(result, "Myeloid", result["cell_type_lvl2"])
            _add_common_states(result, row, pos)
            return result
        
        # ====================================================================
        # 4.5 NK-like: CD3e- AND CD8+ AND GZMB+ AND PanCK- AND CD31-
        # ====================================================================
        if (not pos("CD3e") and pos("CD8") and pos("GZMB") and 
            not pos("PanCK") and not pos("CD31")):
            result["cell_type_lvl1"] = "NK_like"
            result["cell_type_lvl2"] = "NK_like"
            _set_annotation_fields(result, "T_cell", "NK_like")
            _add_common_states(result, row, pos)
            return result
        
        # 未分类的免疫细胞
        result["cell_type_lvl1"] = "Immune_other"
        result["cell_type_lvl2"] = "Immune_other"
        _set_annotation_fields(result, "Myeloid", "Immune_other")
        _add_common_states(result, row, pos)
        return result
    
    # ========================================================================
    # 5. Stromal / Mesenchymal 分支
    # 注意：Endothelial (CD31+) 已在前面处理并返回
    # ========================================================================

    # 5.1 Pericyte 判别已禁用：按需求优先判别 CAF
    # 5.2 CAF (Fibroblast) / 其他 Stromal:
    # - Vimentin 可观测时：CAF 需要 Vimentin+ & SMA+
    # - Vimentin 不可观测时：允许用 SMA+ 直接判别 CAF（避免缺失 marker 导致误分型/丢弃）
    if not pos("PanCK") and not pos("CD31") and not has_immune_marker:
        # CAF: (Vimentin+ 或 Vimentin 不可观测) AND SMA+
        if pos("SMA") and (pos("Vimentin") or (not vimentin_observable)):
            result["major_lineage"] = "Stromal"
            result["cell_type_lvl1"] = "CAF"
            # CAF 按 PYCR1 细分 cell_type_lvl2
            if pos("PYCR1"):
                result["cell_type_lvl2"] = "CAF_PYCR1p"
            else:
                result["cell_type_lvl2"] = "CAF_PYCR1n"
            _set_annotation_fields(result, "CAF", result["cell_type_lvl2"])

            # CAF 状态
            if pos("PD-L1"):
                result["st_CAF_PD-L1_pos"] = True
            if pos("P5CS"):
                result["st_CAF_metabolic_P5CS"] = True
            if pos("PYCR1"):
                result["st_CAF_metabolic_PYCR1"] = True

            _add_common_states(result, row, pos)
            return result

        # 兜底：Vimentin+ 但 SMA- 的其他 Stromal 类型（仅在 Vimentin 可观测/为阳性时使用）
        if pos("Vimentin"):
            result["major_lineage"] = "Stromal"
            result["cell_type_lvl1"] = "Stromal_other"
            result["cell_type_lvl2"] = "Stromal_other"
            _set_annotation_fields(result, "Stromal_other", "Stromal_other")
            _add_common_states(result, row, pos)
            return result
    
    # ========================================================================
    # 6. Unknown（这些细胞将在后续被删除）
    # ========================================================================
    result["major_lineage"] = "Unknown"
    result["cell_type_lvl1"] = "Unknown"
    result["cell_type_lvl2"] = "Unknown"
    _set_annotation_fields(result, "Unknown", "Unknown")
    _add_common_states(result, row, pos)
    return result


def _add_common_states(result: dict, row: pd.Series, pos) -> None:
    """添加通用状态标签。"""
    # 通用状态
    if pos("Ki67"):
        result["st_proliferating"] = True
    if pos("PD-1"):
        result["st_checkpoint_PD1"] = True
    if pos("PD-L1"):
        result["st_checkpoint_PDL1"] = True
    if pos("CTLA-4"):
        result["st_checkpoint_CTLA4"] = True
    if pos("TIGIT"):
        result["st_checkpoint_TIGIT"] = True
    if pos("P5CS"):
        result["st_metabolic_P5CS_high"] = True
    if pos("PYCR1"):
        result["st_metabolic_PYCR1_high"] = True


def annotate_dataframe(df: pd.DataFrame, marker_config: dict) -> pd.DataFrame:
    """
    对整个 DataFrame 执行细胞注释。
    
    Parameters
    ----------
    df : DataFrame
        包含 bool_<marker> 列的 DataFrame
    marker_config : dict
        find_marker_columns 的返回结果
    
    Returns
    -------
    DataFrame
        添加了注释列的 DataFrame
    """
    # 获取 marker 强度列名列表
    marker_intensity_cols = list(marker_config["intensity"].values())
    
    # Vimentin 是否可观测：存在 positivity 列，或（启用阈值推断时）存在强度列
    has_vim_positivity = "Vimentin" in marker_config.get("positivity", {})
    has_vim_intensity = "Vimentin" in marker_config.get("intensity", {})
    vimentin_observable = bool(has_vim_positivity or (INFER_POSITIVITY_FROM_INTENSITY and has_vim_intensity))

    logger.info(f"开始细胞注释... (Vimentin_observable={vimentin_observable})")
    
    # 使用 apply 进行注释
    annotations = df.apply(
        lambda row: annotate_cell(row, marker_intensity_cols, vimentin_observable),
        axis=1, 
        result_type="expand"
    )
    
    # 合并注释结果
    for col in annotations.columns:
        df[col] = annotations[col]
    
    logger.info(f"注释完成，添加了 {len(annotations.columns)} 个注释列")
    
    # 删除 major_lineage=Unknown 的细胞
    unknown_count = (df["major_lineage"] == "Unknown").sum()
    if unknown_count > 0:
        logger.info(f"删除 {unknown_count} 个 Unknown 细胞（major_lineage='Unknown'）")
        df = df[df["major_lineage"] != "Unknown"].copy()
    else:
        logger.info("未发现 Unknown 细胞")
    
    # 打印统计
    logger.info("细胞类型统计:")
    logger.info(df["major_lineage"].value_counts().to_string())
    
    return df


# ============================================================================
# 构建 AnnData
# ============================================================================

def sanitize_column_name(name: str) -> str:
    """
    清洗列名，将空格和特殊字符替换为下划线。
    只保留字母、数字、下划线、点和连字符。
    """
    # 先将空格替换为下划线
    name = name.replace(" ", "_")
    # 移除不允许的字符（只保留字母数字、下划线、点、连字符）
    name = re.sub(r"[^a-zA-Z0-9_.\-]", "_", name)
    # 合并连续的下划线
    name = re.sub(r"_+", "_", name)
    # 去除首尾下划线
    name = name.strip("_")
    return name


def _make_unique(names: list[str]) -> list[str]:
    """将列名列表改为唯一（同名自动追加 _1/_2/...），避免 AnnData/SpatialData 因重复列名报错。"""
    seen: dict[str, int] = {}
    out: list[str] = []
    for n in names:
        if n not in seen:
            seen[n] = 0
            out.append(n)
        else:
            seen[n] += 1
            out.append(f"{n}_{seen[n]}")
    return out


def _select_obs_columns(df: pd.DataFrame, marker_cols: set[str]) -> list[str]:
    """
    选择写入 obs 的列（尽量只保留元数据/positivity/注释结果），避免把超大的 MX 色彩通道塞进 obs。
    """
    keep: list[str] = []
    for c in df.columns:
        if c in marker_cols:
            continue

        if c in {COORD_X, COORD_Y, "Name", "Image", "LayerData", "MX - Type"}:
            keep.append(c)
            continue
        if c.startswith("Study level"):
            keep.append(c)
            continue
        if c.startswith("Object info -"):
            keep.append(c)
            continue

        if c.startswith(POSITIVITY_PREFIX):
            keep.append(c)
            continue

        if c.startswith("bool_") or c in {"major_lineage", "cell_type_lvl1", "cell_type_lvl2", "annotation", "annotation_coarse", "annotation_fine"} or c.startswith("st_"):
            keep.append(c)
            continue

        if re.match(r"^MX\s*-\s*.*(resolution\s*#\d+\s*-\s*)?(Red|Green|Blue|Hematoxylin)\s*$", c, flags=re.IGNORECASE):
            continue

        continue

    return keep


def build_anndata(df_annot: pd.DataFrame, marker_config: dict, sample_id: str) -> ad.AnnData:
    """
    根据注释后的 DataFrame 构建 AnnData 对象。
    
    Parameters
    ----------
    df_annot : DataFrame
        注释后的 DataFrame
    marker_config : dict
        find_marker_columns 的返回结果
    sample_id : str
        样本 ID
    
    Returns
    -------
    AnnData
        构建的 AnnData 对象
    """
    # 构建表达矩阵 X（使用 marker 强度列）
    marker_cols = list(marker_config["intensity"].values())
    marker_names = list(marker_config["intensity"].keys())
    
    X = df_annot[marker_cols].values.astype(np.float32)
    
    # 构建 obs（尽量只保留元数据/positivity/注释列，避免输出过大）
    exclude_cols = set(marker_cols)
    obs_cols = _select_obs_columns(df_annot, exclude_cols)
    obs = df_annot[obs_cols].copy()
    
    # 清洗 obs 列名（SpatialData 要求）
    obs.columns = _make_unique([sanitize_column_name(c) for c in obs.columns])
    
    # 添加 cell_id
    obs["cell_id"] = [f"{sample_id}_cell_{i}" for i in range(len(obs))]
    obs["region"] = sample_id
    obs.index = obs["cell_id"]
    
    # 构建 var
    var = pd.DataFrame(index=marker_names)
    var.index.name = "marker"
    
    # 创建 AnnData
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    # 添加空间坐标到 obsm
    if COORD_X in df_annot.columns and COORD_Y in df_annot.columns:
        spatial_coords = df_annot[[COORD_X, COORD_Y]].values.astype(np.float64)
        adata.obsm["spatial"] = spatial_coords
    else:
        logger.warning(f"缺少坐标列 {COORD_X} 或 {COORD_Y}")
    
    logger.info(f"构建 AnnData: {adata.n_obs} cells × {adata.n_vars} markers")
    
    return adata


# ============================================================================
# 构建 SpatialData
# ============================================================================

def build_spatialdata_from_adata(adata: ad.AnnData, sample_id: str) -> sd.SpatialData:
    """
    根据 AnnData 构建 SpatialData 对象。
    
    Parameters
    ----------
    adata : AnnData
        AnnData 对象
    sample_id : str
        样本 ID
    
    Returns
    -------
    SpatialData
        构建的 SpatialData 对象
    """
    import geopandas as gpd
    from shapely.geometry import Point
    import spatialdata as sd
    from spatialdata.models import ShapesModel, TableModel

    # 从 obsm["spatial"] 获取坐标
    if "spatial" not in adata.obsm:
        raise ValueError("AnnData 缺少 obsm['spatial'] 坐标信息")
    
    coords = adata.obsm["spatial"]
    
    # 创建点几何（使用 Point）
    geometries = [Point(x, y) for x, y in coords]
    
    # 创建 GeoDataFrame，添加 radius 列（ShapesModel 需要）
    gdf = gpd.GeoDataFrame(
        {
            "geometry": geometries,
            "radius": 1.0  # 默认半径为 1
        },
        index=adata.obs.index
    )
    
    # 使用 ShapesModel 解析
    shapes = ShapesModel.parse(gdf, transformations=None)
    
    # 准备 table 的 AnnData（确保有 region 和 instance_key）
    adata_table = adata.copy()
    adata_table.obs["region"] = sample_id
    adata_table.obs["cell_id"] = adata_table.obs.index
    
    # 使用 TableModel 解析
    table = TableModel.parse(
        adata_table,
        region=sample_id,
        region_key="region",
        instance_key="cell_id"
    )
    
    # 构建 SpatialData
    sdata = sd.SpatialData(
        shapes={sample_id: shapes},
        tables={"table": table}
    )
    
    logger.info(f"构建 SpatialData: shapes={sample_id}, table={adata.n_obs} cells")
    
    return sdata


# ============================================================================
# Niche Analysis 占位函数
# ============================================================================

def run_niche_analysis(adata: ad.AnnData):
    """
    Placeholder for future niche analysis.
    
    Parameters
    ----------
    adata : AnnData
        注释后的 AnnData 对象
    """
    pass


# ============================================================================
# 主程序
# ============================================================================

def process_single_sample(sample_dir: Path, outdir: Path) -> Optional[ad.AnnData]:
    """
    处理单个样本。
    
    Parameters
    ----------
    sample_dir : Path
        样本目录路径
    outdir : Path
        输出目录
    
    Returns
    -------
    AnnData 或 None 如果处理失败
    """
    sample_id = sample_dir.name
    logger.info(f"=" * 60)
    logger.info(f"处理样本: {sample_id}")
    logger.info(f"=" * 60)
    
    # 查找 objects.tsv 文件
    tsv_files = list(sample_dir.glob("*_objects.tsv"))
    if not tsv_files:
        logger.error(f"样本 {sample_id} 下未找到 *_objects.tsv 文件")
        return None
    
    tsv_path = tsv_files[0]
    logger.info(f"读取文件: {tsv_path}")
    
    try:
        # 1. 加载数据
        df = load_sample_objects(tsv_path)
        
        # 2. 查找 marker 列
        marker_config = find_marker_columns(df)
        
        # 3. 推断 marker positivity
        df = infer_marker_positivity(df, marker_config)
        
        # 4. 细胞注释
        df = annotate_dataframe(df, marker_config)
        
        # 5. 构建 AnnData
        adata = build_anndata(df, marker_config, sample_id)
        
        # 6. 保存输出
        outdir.mkdir(parents=True, exist_ok=True)
        
        h5ad_path = outdir / f"{sample_id}.annotated.h5ad"
        
        # 保存 h5ad
        adata.write_h5ad(h5ad_path)
        logger.info(f"已保存 AnnData: {h5ad_path}")
        
        return adata
        
    except Exception as e:
        logger.error(f"处理样本 {sample_id} 时出错: {e}")
        import traceback
        traceback.print_exc()
        return None


def main(cleandata_root: str = None, outdir: str = None):
    """
    主入口函数（批量处理）。
    
    Parameters
    ----------
    cleandata_root : str
        cleandata 根目录路径
    outdir : str
        输出目录路径
    """
    if cleandata_root is None:
        cleandata_root = CLEANDATA_ROOT
    else:
        cleandata_root = Path(cleandata_root)
    
    if outdir is None:
        outdir = OUT_DIR
    else:
        outdir = Path(outdir)
    
    # 发现所有样本目录
    sample_dirs = sorted([
        d for d in cleandata_root.iterdir() 
        if d.is_dir() and list(d.glob("*_objects.tsv"))
    ])
    
    if not sample_dirs:
        raise FileNotFoundError(f"在 {cleandata_root} 下未找到任何样本目录")
    
    logger.info(f"发现 {len(sample_dirs)} 个样本目录")
    logger.info(f"样本列表: {[d.name for d in sample_dirs]}")
    
    # 处理统计
    success_count = 0
    failed_samples = []
    
    for sample_dir in sample_dirs:
        adata = process_single_sample(sample_dir, outdir)
        if adata is not None:
            success_count += 1
            # 预留 niche 分析接口
            run_niche_analysis(adata)
        else:
            failed_samples.append(sample_dir.name)
    
    # 打印汇总
    logger.info("=" * 60)
    logger.info("批处理完成！")
    logger.info(f"成功处理: {success_count}/{len(sample_dirs)} 个样本")
    if failed_samples:
        logger.warning(f"失败样本: {failed_samples}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
