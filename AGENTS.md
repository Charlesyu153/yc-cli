# Repository Guidelines

## Project Structure & Module Organization
- `PCF/`: Python pipelines that build annotated `.h5ad` and `.sdata.zarr` from `cleandata/` (e.g., `build_spatialdata_batch.py`, `build_spatialdata_single.py`).
- Root R scripts (`convert_h5ad_to_seurat.R`, `test_single_sample.R`, `test_conversion.R`) handle h5ad to Seurat conversion.
- `utils/R/`: reusable R functions (QC, niche analysis).
- `script/` and `script_ref/`: analysis notebooks and references; not required for core conversion.
- Data folders: `cleandata/` + `rawdata*/` inputs, `spatial/` Python outputs, `data/` Seurat outputs.

## Build, Test, and Development Commands
- `python PCF/build_spatialdata_batch.py` - process all samples under `cleandata/`.
- `python PCF/build_spatialdata_single.py` - process a single sample into `spatial/`.
- `Rscript convert_h5ad_to_seurat.R [input_dir] [output_dir]` - batch convert `.h5ad` to `.qs`.
- `Rscript test_single_sample.R` - full check on `spatial/20P.annotated.h5ad`.
- `Rscript test_conversion.R` - quick conversion sanity check (writes `.rds`).

## Coding Style & Naming Conventions
- Python: PEP 8, 4-space indents, `snake_case` for functions/vars; keep config in module-level constants (see `ROOT`, `CANONICAL_MARKERS`).
- R: 2-space indents, `snake_case` function names; keep console output via `cat()`/`glue()`.
- Prefer `Path`/`file.path` for paths; if you change the hard-coded `ROOT` in Python, update other scripts consistently.

## Testing Guidelines
- No formal test framework; rely on `test_*.R` scripts and verify outputs in `spatial/` and `data/`.
- Avoid committing large binary outputs (`.h5ad`, `.zarr`, `.qs`, `.rds`) unless explicitly requested.

## Commit & Pull Request Guidelines
- No established commit history; use short imperative messages (e.g., "Add batch h5ad conversion").
- PRs should describe data inputs/outputs, list commands run, and note any new dependencies or environment variables.

## Configuration & Environment
- R packages: `Seurat`, `SeuratObject`, `reticulate`, `qs`, `glue`.
- Python env should include `scanpy`, `spatialdata`, `anndata`; R scripts honor `PCF_PYTHON=/path/to/python3`.
- Optional flags for positivity inference: `PCF_INFER_POSITIVITY_FROM_INTENSITY=1` and `PCF_INTENSITY_POSITIVITY_THRESHOLD`.
