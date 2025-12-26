#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Normalize phenotype tables under rawdata by standardizing marker names."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd

# ---------------------------------------------------------------------------
# 配置
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parents[1]
RAW_DATA_DIR = BASE_DIR / "rawdata"
OUTPUT_SUFFIX = "_markers.tsv"

CANONICAL_MARKERS: List[str] = [
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
]
CANONICAL_SET = set(CANONICAL_MARKERS)



UUID_PREFIX_RE = re.compile(r"^[0-9a-f]{4,}(?:-[0-9a-f]{2,}){1,}-", re.IGNORECASE)


def normalize_marker_label(label: str) -> str:
    """Strip UUID/Positivity/dye suffixes to extract the marker token."""
    cleaned = label.strip()
    cleaned = re.sub(r"^Positivity\s*-\s*", "", cleaned, flags=re.IGNORECASE)
    cleaned = UUID_PREFIX_RE.sub("", cleaned)
    cleaned = re.sub(r"\s*\([^)]*\)", "", cleaned)
    cleaned = re.sub(r"\s+", " ", cleaned)
    return cleaned.strip()


def split_marker_tokens(text: str) -> List[str]:
    """Split 'Name' strings into marker tokens terminated by +/-."""
    tokens: List[str] = []
    buff: List[str] = []
    for char in text:
        buff.append(char)
        if char in "+-":
            token = "".join(buff).strip()
            if token:
                tokens.append(token)
            buff.clear()
    if buff:
        tail = "".join(buff).strip()
        if tail:
            tokens.append(tail)
    return tokens


def normalize_row_name(value: object) -> Optional[str]:
    """Convert a raw Name entry into canonical marker +/- combos."""
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text:
        return None
    if text.lower() == "negative":
        return "Negative"

    normalized_tokens: List[str] = []
    for token in split_marker_tokens(text):
        sign = token[-1] if token and token[-1] in "+-" else ""
        marker_body = token[:-1] if sign else token
        normalized = normalize_marker_label(marker_body)
        final_name = MARKER_ALIAS.get(normalized)
        if final_name and final_name in CANONICAL_SET:
            normalized_tokens.append(f"{final_name}{sign}")

    if not normalized_tokens:
        return None
    return " ".join(normalized_tokens)


def standardize_marker_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename marker columns to canonical names and drop unsupported ones."""
    rename_map: Dict[str, str] = {}
    for col in df.columns:
        normalized = normalize_marker_label(col)
        final_name = MARKER_ALIAS.get(normalized)
        if final_name:
            rename_map[col] = final_name
    df = df.rename(columns=rename_map)
    keep_cols = ["Name", "Count"] + [m for m in CANONICAL_MARKERS if m in df.columns]
    missing_core = [col for col in ["Name", "Count"] if col not in df.columns]
    if missing_core:
        raise ValueError(f"Missing required columns {missing_core}")
    return df.loc[:, keep_cols]


def normalize_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Apply column and row normalization to a phenotype dataframe."""
    df = standardize_marker_columns(df)
    df["Name"] = df["Name"].map(normalize_row_name)
    df = df[df["Name"].notna()].copy()
    df.reset_index(drop=True, inplace=True)
    return df


def find_phenotype_files(root: Path) -> Iterable[Path]:
    """Yield all phenotype TSV files under a root directory."""
    return sorted(root.rglob("*Phenotypes*.tsv"))


def process_file(tsv_path: Path) -> Path:
    """Normalize a single phenotype TSV and write alongside the source."""
    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    normalized = normalize_dataframe(df)
    output_path = tsv_path.with_name(f"{tsv_path.stem}{OUTPUT_SUFFIX}")
    normalized.to_csv(output_path, sep="\t", index=False)
    return output_path


def main() -> None:
    if not RAW_DATA_DIR.exists():
        raise SystemExit(f"rawdata directory not found: {RAW_DATA_DIR}")

    files = list(find_phenotype_files(RAW_DATA_DIR))
    if not files:
        raise SystemExit(f"No phenotype TSV files found under: {RAW_DATA_DIR}")

    for src in files:
        dst = process_file(src)
        rel_src = src.relative_to(BASE_DIR)
        rel_dst = dst.relative_to(BASE_DIR)
        print(f"Normalized {rel_src} -> {rel_dst}")


if __name__ == "__main__":
    main()
MARKER_ALIAS: Dict[str, str] = {
    "CD163": "CD163",
    "CD20": "CD20",
    "CD31": "CD31",
    "CD3e": "CD3e",
    "CD4": "CD4",
    "CD66": "CD66",
    "CD68": "CD68",
    "CD8": "CD8",
    "CD86": "CD86",
    "CTLA-4": "CTLA-4",
    "CLDN18.2": "CLDN18.2",
    "DAB CLDN18.2": "CLDN18.2",
    "DAB": "CLDN18.2",
    "FOXP3": "FOXP3",
    "TIGIT": "TIGIT",
    "Granzyme B": "GZMB",
    "GZMB": "GZMB",
    "Ki67": "Ki67",
    "P5CS": "P5CS",
    "PD-1": "PD-1",
    "PD-L1": "PD-L1",
    "PYCR1": "PYCR1",
    "Pan-Cytokeratin": "PanCK",
    "Pan Cytokeratin": "PanCK",
    "PanCK": "PanCK",
    "SMA": "SMA",
    "Vimentin": "Vimentin",
}