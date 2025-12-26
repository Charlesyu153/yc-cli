#!/usr/bin/env python3
"""
批处理：遍历 rawdata 下所有 *_objects.tsv，使用与 clean_one_objects_example.py 相同逻辑进行清洗，
输出到 cleandata/ 对应路径。
"""

from __future__ import annotations

import os
import re
from pathlib import Path
import pandas as pd

# 直接复用原脚本的设置
ROOT = Path("/home/jacekyu/PCF")
RAW = ROOT / "rawdata"
OUT = ROOT / "cleandata"
OUT.mkdir(parents=True, exist_ok=True)

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


def _regex_for_token(token: str) -> str:
    """
    构建宽松的正则：允许空格/连字符/下划线替换，并做单词边界。
    """
    escaped = re.escape(token.lower())
    # 替换为可选的分隔符，需要双反斜杠避免替换串中的转义报错
    pattern = re.sub(r"[-\s_]+", r"[-\\s_]*", escaped)
    return rf"\b{pattern}\b"


def detect_marker_from_text(text: str) -> str | None:
    """
    在任意字符串中搜索已知 marker/别名，命中返回标准名。
    """
    lower = text.lower()

    # 先匹配别名，保证 granzyme b / pan cytokeratin 等被转换
    for alias, canon in ALIAS.items():
        if re.search(_regex_for_token(alias), lower):
            return canon

    # 再匹配标准名
    for canon in CANONICAL:
        if re.search(_regex_for_token(canon), lower):
            return canon

    return None


def normalize_marker(token: str) -> str | None:
    # 对新数据：优先从 "resolution #n-<Marker>" 抽取末尾 token，避免在超长字符串里靠全文搜索碰运气
    extracted = _extract_resolution_suffix(token)
    if extracted:
        token = extracted

    token = re.sub(r"^\s*(mx|dab)\s*-\s*", "", token, flags=re.IGNORECASE)
    lower = token.lower().strip()

    if lower == "dab":
        return "CLDN18.2"

    cleaned = re.sub(r"\s*\((mx|dab)\)\s*$", "", token, flags=re.IGNORECASE)
    # 兼容荧光染料等后缀括号：如 CD86 (AF 750) (MX)、CD3e (Cy5) (MX)
    # 循环剥离末尾括号块，直到末尾不再是括号。
    cleaned = cleaned.strip()
    while True:
        stripped = re.sub(r"\s*\([^)]*\)\s*$", "", cleaned).strip()
        if stripped == cleaned:
            break
        cleaned = stripped
    cleaned = re.sub(
        r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}-",
        "",
        cleaned,
        flags=re.IGNORECASE,
    )
    cleaned = re.sub(r"^\s*(mx|dab)\s*-\s*", "", cleaned, flags=re.IGNORECASE)

    cleaned_lower = cleaned.lower().strip()
    if cleaned_lower == "dab":
        return "CLDN18.2"

    alias = ALIAS.get(cleaned.lower())
    if alias:
        cleaned = alias
        cleaned_lower = cleaned.lower()

    for canon in CANONICAL:
        if canon.lower() == cleaned_lower:
            return canon

    # 宽松匹配：在整段字符串里找 marker 关键词（处理长 MX 列名）
    detected = detect_marker_from_text(cleaned_lower)
    if detected:
        return detected

    return None


def parse_column(col: str):
    lower = col.lower().strip()

    # 规则：CLDN18.2 仅接受 DAB 来源。若列名包含 CLDN18.2 但不含 DAB，则直接跳过（不输出、不参与 marker 选择）。
    # 这同时避免了 20P 这类列名包含 "CLDN18.2-IHC" 导致被误删的历史问题（DAB 列自身包含 "dab"，不会被 skip）。
    if "cldn18.2" in lower and "dab" not in lower:
        return "skip", None, None

    # 兼容多种 positivity 前缀：Positivity-xxx / Positivity - xxx / positive_xxx
    if lower.startswith("positive_") or lower.startswith("positivity"):
        if lower.startswith("positive_"):
            remainder = col[len("positive_") :]
        else:
            # 去掉前缀“positivity”后的分隔符（-/_/空格）
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

    chosen_for_marker = {}
    for col in df.columns:
        kind, name, _ = parse_column(col)
        if kind != "marker" or not name:
            continue
        if name not in chosen_for_marker:
            chosen_for_marker[name] = col
        else:
            current = chosen_for_marker[name]
            if "dab" in col.lower() and "dab" not in current.lower():
                chosen_for_marker[name] = col

    out_df = pd.DataFrame(index=df.index)

    for col in df.columns:
        kind, name, raw_remainder = parse_column(col)
        if kind == "skip":
            continue
        if kind == "other":
            out_df[col] = df[col]
        elif kind == "positivity":
            if name:
                new_col = f"Positivity-{name}"
            else:
                new_col = f"Positivity-{raw_remainder}"
            out_df[new_col] = df[col]
        elif kind == "marker":
            if not name:
                out_df[col] = df[col]
                continue
            source = chosen_for_marker.get(name)
            if source == col:
                out_df[name] = df[col]

    # 自检：若原始数据存在 CD86 相关列，但输出缺失 CD86/Positivity-CD86，则告警（可选严格报错）
    raw_cd86_cols = [c for c in df.columns if "cd86" in c.lower()]
    if raw_cd86_cols:
        missing = []
        if "CD86" not in out_df.columns:
            missing.append("CD86")
        if "Positivity-CD86" not in out_df.columns:
            missing.append("Positivity-CD86")
        if missing:
            sample_id = path.parent.name
            msg = (
                f"[WARN] sample={sample_id} cleaned output missing {missing}; "
                f"raw CD86-like columns={raw_cd86_cols}"
            )
            if os.environ.get("PCF_STRICT_MARKERS") == "1":
                raise RuntimeError(msg)
            print(msg)

    # 自检：若原始数据存在 CD163 相关列，但输出缺失 CD163/Positivity-CD163，则告警（可选严格报错）
    raw_cd163_cols = [c for c in df.columns if "cd163" in c.lower()]
    if raw_cd163_cols:
        missing = []
        if "CD163" not in out_df.columns:
            missing.append("CD163")
        if "Positivity-CD163" not in out_df.columns:
            missing.append("Positivity-CD163")
        if missing:
            sample_id = path.parent.name
            msg = (
                f"[WARN] sample={sample_id} cleaned output missing {missing}; "
                f"raw CD163-like columns={raw_cd163_cols}"
            )
            if os.environ.get("PCF_STRICT_MARKERS") == "1":
                raise RuntimeError(msg)
            print(msg)

    rel = path.relative_to(RAW)
    out_path = OUT / rel
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"已输出: {out_path}")


def main():
    # 遍历 rawdata 下的所有 *_objects.tsv
    files = sorted(RAW.rglob("*_objects.tsv"))
    if not files:
        raise FileNotFoundError(f"未在 {RAW} 下找到 *_objects.tsv")
    for f in files:
        print(f"处理: {f}")
        clean_one_file(f)


if __name__ == "__main__":
    main()