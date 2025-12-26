#!/usr/bin/env python3
"""
单样本示例脚本：清洗 rawdata 下某个 objects.tsv 文件的 marker 列名。
逻辑：
- 去除 marker 列名中的 UUID 前缀与 (MX)/(DAB) 后缀，只保留标准 marker 名。
- 裸列名为 DAB 时视为 CLDN18.2。
- 同一 marker 若同时存在 DAB/MX，优先保留 DAB 数据。
- positive_* 前缀改写为 positivity-<marker>（细胞注释结果）。
- 非 marker 列保持原名。
输出：写入 cleandata/ 下与原路径对应的位置。
"""

from __future__ import annotations

import re
from pathlib import Path
import pandas as pd

# 根目录（可根据需要调整）
ROOT = Path("/home/jacekyu/PCF")
RAW = ROOT / "rawdata"
OUT = ROOT / "cleandata"
OUT.mkdir(parents=True, exist_ok=True)

# 允许的标准 marker 名（增加 DAPI 作为特例保留）
CANONICAL = {
    "CD163",
    "CD20",
    "CD31",
    "CD3e",
    "CD4",
    "CD66",
    "CD68",
    "CD8",
    "CD86",
    "CTLA-4",
    "CLDN18.2",
    "FOXP3",
    "TIGIT",
    "GZMB",
    "Ki67",
    "P5CS",
    "PD-1",
    "PD-L1",
    "PYCR1",
    "PanCK",
    "SMA",
    "Vimentin",
    "DAPI",
}

ALIAS = {
    "pan-cytokeratin": "PanCK",
    "pan cytokeratin": "PanCK",
    "granzyme b": "GZMB",
    "granzyme-b": "GZMB",
    "dab cldn18.2": "CLDN18.2",
}


def _extract_resolution_suffix(text: str) -> str | None:
    """
    从长列名中提取 marker token：
    - 形如 "... - resolution #1-CD163 (MX)" -> "CD163 (MX)"
    - 形如 "... - resolution #1-DAB" -> "DAB"
    若不存在该模式则返回 None（兼容旧数据：UUID-<marker> 等）。
    """
    m = re.search(r"resolution\s*#\d+\s*-\s*(.+)$", text, flags=re.IGNORECASE)
    if not m:
        return None
    return m.group(1).strip()


def normalize_marker(token: str) -> str | None:
    """标准化 marker 名；若不在白名单，返回 None。"""
    extracted = _extract_resolution_suffix(token)
    if extracted:
        token = extracted

    # 去掉前缀 "MX -", "DAB -" 等
    token = re.sub(r"^\s*(mx|dab)\s*-\s*", "", token, flags=re.IGNORECASE)

    lower = token.lower().strip()

    # 裸 DAB 视为 CLDN18.2
    if lower == "dab":
        return "CLDN18.2"

    # 去掉 (MX)/(DAB) 及其他括号链后缀：如 CD86 (AF 750) (MX)、CD3e (Cy5) (MX)
    cleaned = token.strip()
    while True:
        stripped = re.sub(r"\s*\([^)]*\)\s*$", "", cleaned).strip()
        if stripped == cleaned:
            break
        cleaned = stripped

    # 去掉 UUID 前缀（8-4-4-4-12 hex）
    cleaned = re.sub(
        r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}-",
        "",
        cleaned,
        flags=re.IGNORECASE,
    )

    # 去掉前置的 "MX - " 再次尝试
    cleaned = re.sub(r"^\s*(mx|dab)\s*-\s*", "", cleaned, flags=re.IGNORECASE)

    cleaned_lower = cleaned.lower().strip()

    # 去除 UUID 后若变成 DAB，视为 CLDN18.2
    if cleaned_lower == "dab":
        return "CLDN18.2"

    # alias 归一
    alias = ALIAS.get(cleaned.lower())
    if alias:
        cleaned = alias
        cleaned_lower = cleaned.lower()

    for canon in CANONICAL:
        if canon.lower() == cleaned_lower:
            # CLDN18.2 仅接受含 DAB 信息的列（或裸 DAB 已提前返回）
            if canon == "CLDN18.2":
                if "dab" not in token.lower() and "dab" not in cleaned.lower():
                    return None
            return canon
    return None


def parse_column(col: str):
    """
    返回 (kind, name):
    - kind == "positivity": positivity 列，name 是标准 marker 或 None
    - kind == "marker": marker 数据列，name 是标准 marker 或 None
    - kind == "other": 非 marker 列
    """
    lower = col.lower().strip()

    # 若包含 CLDN18.2 但不含 DAB，跳过（不用作数据列，也不保留）
    # 注意：某些样本的 DAB 列名可能在文件名前缀中包含 "CLDN18.2-IHC"，但列本身也包含 "DAB"，不会被 skip。
    if "cldn18.2" in lower and "dab" not in lower:
        return "skip_cldn", None, None

    # 兼容 "Positivity - xxx" / "positive_xxx"
    if lower.startswith("positive_") or lower.startswith("positivity"):
        if lower.startswith("positive_"):
            remainder = col[len("positive_") :]
        else:
            remainder = re.sub(r"^\s*positivity\s*[-_\s]*", "", col, flags=re.IGNORECASE)
        remainder = remainder.strip()
        name = normalize_marker(remainder)
        return "positivity", name, remainder

    name = normalize_marker(col)
    if name:
        return "marker", name, None
    return "other", None, None


def clean_one_file(path: Path):
    df = pd.read_csv(path, sep="\t")

    # 先移除现有包含 CLDN18.2 的列（包括 positivity / marker），避免混入旧列
    drop_cols = [c for c in df.columns if "cldn18.2" in c.lower()]
    if drop_cols:
        df = df.drop(columns=drop_cols)

    # 选择每个 marker 的数据列，优先 DAB；CLDN18.2 仅允许 DAB 来源
    chosen_for_marker = {}  # name -> source column
    for col in df.columns:
        kind, name, _ = parse_column(col)
        if kind != "marker" or not name:
            continue
        if name == "CLDN18.2":
            # 仅接受含 DAB 的列或裸 DAB
            if "dab" not in col.lower():
                continue
        if name not in chosen_for_marker:
            chosen_for_marker[name] = col
        else:
            # 若已有，优先含 DAB 的列
            current = chosen_for_marker[name]
            if "dab" in col.lower() and "dab" not in current.lower():
                chosen_for_marker[name] = col

    out_df = pd.DataFrame(index=df.index)

    for col in df.columns:
        kind, name, raw_remainder = parse_column(col)
        if kind == "skip_cldn":
            continue
        if kind == "other":
            out_df[col] = df[col]
        elif kind == "positivity":
            if name:
                new_col = f"Positivity-{name}"
            else:
                # 无法标准化时，尽量用清洗后的 remainder
                new_col = f"Positivity-{raw_remainder}"
            out_df[new_col] = df[col]
        elif kind == "marker":
            if not name:
                # 未标准化的 marker 列，保持原名
                out_df[col] = df[col]
                continue
            source = chosen_for_marker.get(name)
            if source == col:
                out_df[name] = df[col]

    rel = path.relative_to(RAW)
    out_path = OUT / rel
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"已输出: {out_path}")


def main():
    # 这里示例只处理一个文件，请按需修改 input_path
    input_path = RAW / "1P" / "1P_objects.tsv"  # 示例路径
    if not input_path.exists():
        raise FileNotFoundError(f"未找到示例文件: {input_path}")
    clean_one_file(input_path)


if __name__ == "__main__":
    main()

