#!/usr/bin/env python3
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass
class Finding:
    level: str  # "ERROR" | "WARN"
    file: Path
    message: str


def extract_frontmatter(text: str) -> tuple[dict[str, Any], str] | None:
    if not text.startswith("---\n"):
        return None

    end = text.find("\n---\n", 4)
    if end == -1:
        return None

    fm_text = text[4:end]
    body = text[end + 5 :]
    data = yaml.safe_load(fm_text) or {}
    if not isinstance(data, dict):
        raise ValueError("Frontmatter must be a YAML mapping (dict).")
    return data, body


def looks_like_git_rev(value: str) -> bool:
    if not value:
        return False
    if not (7 <= len(value) <= 40):
        return False
    return all(c in "0123456789abcdef" for c in value.lower())


def validate_related_files(
    repo_root: Path, file_path: Path, frontmatter: dict[str, Any]
) -> list[Finding]:
    findings: list[Finding] = []

    related_files = frontmatter.get("related_files") or []
    if related_files is None:
        related_files = []
    if not isinstance(related_files, list):
        findings.append(
            Finding("ERROR", file_path, "`related_files` must be a list (or omitted).")
        )
        return findings

    for i, entry in enumerate(related_files):
        if not isinstance(entry, dict):
            findings.append(
                Finding("ERROR", file_path, f"`related_files[{i}]` must be a mapping.")
            )
            continue

        path = entry.get("path")
        if not isinstance(path, str) or not path.strip():
            findings.append(
                Finding("ERROR", file_path, f"`related_files[{i}].path` is required.")
            )
            continue

        locator = entry.get("locator")
        if not isinstance(locator, dict):
            findings.append(
                Finding(
                    "ERROR",
                    file_path,
                    f"`related_files[{i}].locator` must be a mapping with `rg_pattern` or `symbol`.",
                )
            )
            continue

        rg_pattern = locator.get("rg_pattern")
        symbol = locator.get("symbol")
        line = locator.get("line")

        has_rg = isinstance(rg_pattern, str) and rg_pattern.strip()
        has_symbol = isinstance(symbol, str) and symbol.strip()
        has_line = isinstance(line, int) and line > 0

        if not (has_rg or has_symbol):
            findings.append(
                Finding(
                    "ERROR",
                    file_path,
                    f"`related_files[{i}].locator` must include non-empty `rg_pattern` or `symbol`.",
                )
            )
        elif has_line and not (has_rg or has_symbol):
            findings.append(
                Finding(
                    "WARN",
                    file_path,
                    f"`related_files[{i}].locator.line` is fragile; prefer `rg_pattern`.",
                )
            )

        git_rev = entry.get("git_rev")
        if isinstance(git_rev, str) and git_rev.strip() and not looks_like_git_rev(git_rev.strip()):
            findings.append(
                Finding(
                    "WARN",
                    file_path,
                    f"`related_files[{i}].git_rev` does not look like a git commit hash: {git_rev!r}",
                )
            )

        # path existence check (warn-only)
        resolved = (repo_root / path).resolve()
        if not resolved.exists():
            findings.append(
                Finding(
                    "WARN",
                    file_path,
                    f"`related_files[{i}].path` not found at {path!r} (repo-relative).",
                )
            )

    return findings


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    insights_glob = repo_root / "ai-docs" / "current"
    insights = sorted(insights_glob.glob("*/insights/*.md"))

    findings: list[Finding] = []

    for file_path in insights:
        try:
            text = file_path.read_text(encoding="utf-8")
        except Exception as exc:
            findings.append(Finding("ERROR", file_path, f"Failed to read: {exc}"))
            continue

        try:
            extracted = extract_frontmatter(text)
        except Exception as exc:
            findings.append(Finding("ERROR", file_path, f"Invalid frontmatter: {exc}"))
            continue

        if extracted is None:
            findings.append(Finding("ERROR", file_path, "Missing or unterminated YAML frontmatter."))
            continue

        frontmatter, _ = extracted

        insight_id = frontmatter.get("id")
        if not isinstance(insight_id, str) or not insight_id.strip():
            findings.append(Finding("ERROR", file_path, "Frontmatter `id` is required."))

        findings.extend(validate_related_files(repo_root, file_path, frontmatter))

    for f in findings:
        print(f"{f.level}: {f.file.relative_to(repo_root)}: {f.message}")

    has_error = any(f.level == "ERROR" for f in findings)
    if has_error:
        return 1

    print(f"OK: validated {len(insights)} insight file(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

