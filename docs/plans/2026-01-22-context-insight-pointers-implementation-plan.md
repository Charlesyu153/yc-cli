# Context Insight Pointers Implementation Plan

> **For Codex:** REQUIRED SUB-SKILL: Use `superpowers:executing-plans` to implement this plan task-by-task.

**Goal:** Make `insights/*.md` records robust for frequently-changing bioinformatics projects by storing code pointers (file + git rev + locator), without indexing the whole codebase.

**Architecture:** Keep Markdown files as the source of truth. Only `MANIFEST.md` + `insights/*.md` are indexed. Insights link to scripts via frontmatter `related_files`, and retrieval uses those pointers to read targeted code snippets on demand.

**Tech Stack:** Markdown/YAML frontmatter, optional vector index (Chroma) later, simple Python validation script.

---

### Task 1: Add templates and example docs

**Files:**
- Create: `ai-docs/templates/insight-template.md`
- Create: `ai-docs/templates/manifest-template.md`
- Create: `ai-docs/templates/insights-index-template.txt`
- Create: `ai-docs/current/default/MANIFEST.md`
- Create: `ai-docs/current/default/insights/example_decision.md`
- Create: `ai-docs/.context/insights/default.txt`

**Step 1: Add templates**

Create templates with required frontmatter fields:
- `related_files[].path` (repo-relative)
- `related_files[].git_rev`
- `related_files[].locator.rg_pattern` (preferred) or `locator.symbol`
- `last_verified`, `status`, `supersedes`, `superseded_by`

**Step 2: Add example insight**

Add a short example that demonstrates pointers without copying code.

**Step 3: Verify structure**

Run: `find ai-docs -maxdepth 4 -type f -print`
Expected: templates + default task tree exist.

---

### Task 2: Update `context-insight` skill to require pointers

**Files:**
- Modify: `skills/context-insight/SKILL.md`
- Modify: `skills/SKILL.md`
- Modify: `skills.md`

**Step 1: Update insight format**

Document the required `related_files` structure and locator priority:
`git_rev + rg_pattern` > `symbol` > `line` (line numbers are optional hints only).

**Step 2: Update workflow**

Add a retrieval micro-workflow:
1) search insights
2) read the chosen insight
3) follow pointers to read only relevant code snippets using `rg -n -C 3`

**Step 3: Verify docs consistency**

Run: `rg -n "~/.claude" -S skills.md skills/SKILL.md skills/context-insight/SKILL.md`
Expected: no misleading references to an unrelated directory layout.

---

### Task 3: Update RAG plan docs to align with pointer-based approach

**Files:**
- Modify: `docs/plans/2026-01-22-rag-system-design.md`
- Modify: `docs/plans/2026-01-22-rag-session-startup.md`
- Modify: `docs/plans/2026-01-22-rag-implementation-steps.md`

**Step 1: Clarify scope**

Explicitly state the vector index is built from `MANIFEST.md` + `insights/*.md` only.

**Step 2: Fix ID / update semantics**

In code snippets, change document IDs to be stable (e.g., `sha1(file_path)`), and recommend upsert. Avoid timestamp-based IDs.

**Step 3: Include pointers in retrieval**

Document that search results should expose `related_files` so the AI can fetch precise code context on demand.

---

### Task 4: Add a validator for insights

**Files:**
- Create: `scripts/validate_insights.py`

**Step 1: Implement validator**

The validator should:
- Walk `ai-docs/current/**/insights/*.md`
- Parse YAML frontmatter
- Ensure `related_files` entries have `path`, `locator` (`rg_pattern` or `symbol`), and optionally `git_rev`
- Verify paths exist relative to repo root (warn, don’t hard-fail, if absent)
- Exit non-zero on malformed frontmatter

**Step 2: Run validator**

Run: `python scripts/validate_insights.py`
Expected: exits 0; prints warnings only if example paths don’t exist.

