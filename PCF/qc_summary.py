#!/usr/bin/env python3
"""
Generate QC summaries for annotated h5ad files.

Outputs:
  - data/qc_summary_samples.tsv
  - data/qc_summary_regions.tsv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# ============================================================================
# Configuration
# ============================================================================

ROOT = Path("/home/jacekyu/PCF")
DEFAULT_INPUT_DIR = ROOT / "spatial"
DEFAULT_OUTPUT_SAMPLES = ROOT / "data" / "qc_summary_samples.tsv"
DEFAULT_OUTPUT_REGIONS = ROOT / "data" / "qc_summary_regions.tsv"

COORD_X = "Center_X"
COORD_Y = "Center_Y"

ANNOTATION = "annotation"
ANNOTATION_COARSE = "annotation_coarse"
ANNOTATION_FINE = "annotation_fine"
MAJOR_LINEAGE = "major_lineage"
REGION = "region"

NEGATIVE_LABELS = {"Negative_or_Artifact", "Unknown"}
KEY_POSITIVITY_MARKERS = ["P5CS", "Vimentin", "CD20", "CD3e"]


def _sanitize_label(label: str) -> str:
    return re.sub(r"[^0-9A-Za-z]+", "_", str(label)).strip("_")


def _make_filter_mask(obs: pd.DataFrame) -> pd.Series:
    mask = pd.Series(False, index=obs.index)
    if ANNOTATION in obs.columns:
        mask |= obs[ANNOTATION].isin(NEGATIVE_LABELS)
    if MAJOR_LINEAGE in obs.columns:
        mask |= obs[MAJOR_LINEAGE].isin(NEGATIVE_LABELS)
    return mask


def _coord_summary(obs: pd.DataFrame) -> dict:
    summary = {
        "coord_missing": None,
        "coord_missing_pct": None,
        "coord_min_x": None,
        "coord_max_x": None,
        "coord_range_x": None,
        "coord_min_y": None,
        "coord_max_y": None,
        "coord_range_y": None,
        "flag_coord_missing": True,
        "flag_coord_zero_range": True,
    }

    if COORD_X not in obs.columns or COORD_Y not in obs.columns:
        return summary

    x = pd.to_numeric(obs[COORD_X], errors="coerce")
    y = pd.to_numeric(obs[COORD_Y], errors="coerce")
    missing = ~np.isfinite(x) | ~np.isfinite(y)
    summary["coord_missing"] = int(missing.sum())
    summary["coord_missing_pct"] = float(missing.mean()) if len(obs) else 0.0
    summary["flag_coord_missing"] = summary["coord_missing"] > 0

    finite = ~missing
    if finite.any():
        x_f = x[finite]
        y_f = y[finite]
        summary["coord_min_x"] = float(x_f.min())
        summary["coord_max_x"] = float(x_f.max())
        summary["coord_range_x"] = float(x_f.max() - x_f.min())
        summary["coord_min_y"] = float(y_f.min())
        summary["coord_max_y"] = float(y_f.max())
        summary["coord_range_y"] = float(y_f.max() - y_f.min())
        summary["flag_coord_zero_range"] = (
            summary["coord_range_x"] <= 0.0 or summary["coord_range_y"] <= 0.0
        )

    return summary


def _annotation_missing_counts(obs: pd.DataFrame) -> dict:
    summary = {
        "annotation_missing": None,
        "annotation_coarse_missing": None,
        "annotation_fine_missing": None,
    }
    if ANNOTATION in obs.columns:
        summary["annotation_missing"] = int(obs[ANNOTATION].isna().sum())
    if ANNOTATION_COARSE in obs.columns:
        summary["annotation_coarse_missing"] = int(obs[ANNOTATION_COARSE].isna().sum())
    if ANNOTATION_FINE in obs.columns:
        summary["annotation_fine_missing"] = int(obs[ANNOTATION_FINE].isna().sum())
    return summary


def _coarse_counts(obs: pd.DataFrame) -> tuple[dict, dict]:
    if ANNOTATION_COARSE not in obs.columns:
        return {}, {}

    counts = obs[ANNOTATION_COARSE].value_counts(dropna=False)
    raw_counts = {str(k): int(v) for k, v in counts.items()}
    sanitized = {
        f"coarse_{_sanitize_label(k)}": int(v) for k, v in raw_counts.items()
    }
    return raw_counts, sanitized


def _positivity_missing(obs: pd.DataFrame) -> dict:
    summary = {}
    for marker in KEY_POSITIVITY_MARKERS:
        col = f"Positivity-{marker}"
        summary[f"missing_pos_{marker}"] = col not in obs.columns
    return summary


def summarize_obs(obs: pd.DataFrame, sample_id: str, region: str | None = None) -> dict:
    mask = _make_filter_mask(obs)
    filtered = obs.loc[~mask]

    summary = {
        "sample_id": sample_id,
        "region": region if region is not None else sample_id,
        "total_cells_raw": int(len(obs)),
        "total_cells_filtered": int(len(filtered)),
        "filtered_out": int(mask.sum()),
    }
    summary.update(_annotation_missing_counts(obs))
    summary.update(_coord_summary(obs))
    summary.update(_positivity_missing(obs))

    raw_counts, sanitized_counts = _coarse_counts(filtered)
    summary.update(sanitized_counts)

    b_count = raw_counts.get("B_cell", 0)
    t_count = raw_counts.get("T_cell", 0)
    summary["B_cell_count"] = b_count
    summary["T_cell_count"] = t_count
    summary["B_cell_pct"] = (
        float(b_count) / summary["total_cells_filtered"] * 100.0
        if summary["total_cells_filtered"] > 0
        else 0.0
    )
    summary["T_cell_pct"] = (
        float(t_count) / summary["total_cells_filtered"] * 100.0
        if summary["total_cells_filtered"] > 0
        else 0.0
    )
    summary["flag_low_B_cell"] = b_count < 30
    summary["flag_low_T_cell"] = t_count < 20
    summary["flag_tls_not_evaluable"] = summary["flag_low_B_cell"] or summary["flag_low_T_cell"]

    return summary


def summarize_file(path: Path) -> tuple[dict, list[dict]]:
    sample_id = path.name.replace(".annotated.h5ad", "")
    adata = ad.read_h5ad(path, backed="r")
    obs = adata.obs

    sample_summary = summarize_obs(obs, sample_id)

    region_summaries: list[dict] = []
    if REGION in obs.columns:
        for region, group in obs.groupby(REGION, sort=True, observed=False):
            region_summaries.append(summarize_obs(group, sample_id, region=str(region)))
    else:
        region_summaries.append(summarize_obs(obs, sample_id, region=sample_id))

    try:
        adata.file.close()
    except Exception:
        pass

    return sample_summary, region_summaries


def _normalize_table(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    count_cols = [c for c in df.columns if c.startswith("coarse_")]
    df[count_cols] = df[count_cols].fillna(0).astype(int)
    return df


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate QC summaries for h5ad files.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT_DIR))
    parser.add_argument("--output-samples", default=str(DEFAULT_OUTPUT_SAMPLES))
    parser.add_argument("--output-regions", default=str(DEFAULT_OUTPUT_REGIONS))
    parser.add_argument("--pattern", default="*.annotated.h5ad")
    args = parser.parse_args()

    input_dir = Path(args.input)
    files = sorted(input_dir.glob(args.pattern))
    if not files:
        raise FileNotFoundError(f"No files matched in {input_dir} with pattern {args.pattern}")

    sample_rows: list[dict] = []
    region_rows: list[dict] = []

    for path in files:
        print(f"Processing {path.name} ...")
        sample_summary, region_summaries = summarize_file(path)
        sample_rows.append(sample_summary)
        region_rows.extend(region_summaries)

    samples_df = _normalize_table(pd.DataFrame(sample_rows))
    regions_df = _normalize_table(pd.DataFrame(region_rows))

    output_samples = Path(args.output_samples)
    output_regions = Path(args.output_regions)
    output_samples.parent.mkdir(parents=True, exist_ok=True)
    output_regions.parent.mkdir(parents=True, exist_ok=True)

    samples_df.to_csv(output_samples, sep="\t", index=False)
    regions_df.to_csv(output_regions, sep="\t", index=False)
    print(f"Wrote sample QC: {output_samples}")
    print(f"Wrote region QC: {output_regions}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
