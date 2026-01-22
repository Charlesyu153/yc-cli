---
title: "QC threshold decisions are stored as pointers"
id: "example_decision"
type: decision
created: "2026-01-22"
tags: ["qc", "threshold", "process"]
status: done
priority: medium
last_verified: "2026-01-22"
supersedes: ""
superseded_by: ""
related_files:
  - path: "pipelines/qc.nf"
    git_rev: ""
    locator:
      rg_pattern: "min_mapq\\s*="
    why: "MAPQ filtering threshold lives here; pattern is stable across edits"
---

## One-Line Conclusion
Record decisions in `insights/*.md`, but link to code with `related_files` pointers instead of copying the script.

## Background / Context
Bioinformatics pipelines change frequently. Line numbers drift, and copying code into notes becomes stale quickly.

## Decision / Insight
- The vector index should include only `MANIFEST.md` and `insights/*.md`.
- Each insight that depends on code must include `related_files` pointers with a robust locator (prefer `rg_pattern`).

## How to Verify (Repro Checklist)
1. Search for this insight by keyword.
2. Read the insight file.
3. Use the locator to jump into the relevant script:
   - `rg -n "min_mapq\\s*=" pipelines/qc.nf -C 3`

