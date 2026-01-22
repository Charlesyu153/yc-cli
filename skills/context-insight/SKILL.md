---
name: context-insight
description: Use this skill when user says "context skill", needs to record insights/lessons for future sessions, or manage project state. This skill provides a lightweight alternative to JSON-based context systems, using Markdown + .txt index format for efficient knowledge retrieval.
---

# Context Insight Management System

## Quick Start

When user says "context skill" or asks to record something:

1. **Read this skill** for the format
2. **Create insight** following the template
3. **Update index** in `.txt` file

## When to Use

| User says | Do this |
|----------|---------|
| "context skill" | Read this skill, follow format |
| "记录这个教训" | Create insight in `insights/` |
| "更新上下文" | Ask: ~/.claude/context-insight or ai-docs/.context/? |
| "用 skill" | Always check ~/.claude/skills/ first |

## Core Principle

**Separate key records (brief index) from detailed records (full content).**

AI reads index first, then loads details on demand. This saves tokens and keeps sessions efficient.

## File Structure

```
ai-docs/current/{task}/
├── MANIFEST.md          # L0: Hot Data (current state, decisions, next actions)
├── {task}.md            # L2: Full document
└── insights/            # L3: Deep insights
    ├── insight_01.md
    └── ...

ai-docs/.context/insights/
└── {task}.txt           # Key Records Index (brief, searchable)
```

## Index Format (.txt)

```text
# {Task Name} - Key Insight Index
# Last updated: YYYY-MM-DD

[insight_id]
Title: Brief title
Summary: One-liner takeaway
Keywords: tag1, tag2
Detail: current/{task}/insights/{insight_id}.md
Created: YYYY-MM-DD
Related: {optional file path}
Status: draft|done|todo
Priority: high|medium|low
```

## Insight Format (.md)

```markdown
---
title: Insight Title
id: unique_id
created: YYYY-MM-DD
tags: tag1, tag2
status: draft|done|todo
priority: high|medium|low
---

## Core Question/Conclusion
{One sentence}

## Background
{Why this matters}

## Detailed Content
{Full explanation}

## Key Findings
1. Finding 1
2. Finding 2

## Related Files
- {path}: {purpose}
```

## MANIFEST.md Format

```markdown
# Project Manifest: {Project Name}

**Last Updated**: YYYY-MM-DD HH:MM
**Session Context**: {Brief description}

---

## Current State
| Aspect | Status |
|--------|--------|
| Phase | Planning/Implementation/etc. |
| Focus | What we're doing NOW |
| Blocker | Any blockers |

---

## Active Decisions (Last 7 Days)
| Date | Decision | Rationale |
|------|----------|-----------|

---

## Current Questions
| ID | Question | Status |

---

## Next Actions
- [ ] Action 1
- [ ] Action 2
```

## Workflow

### Creating a New Insight

1. Copy template: `cp ai-docs/templates/insight-template.md ai-docs/current/{task}/insights/{new_insight}.md`
2. Fill content
3. Update index: `vim ai-docs/.context/insights/{task}.txt`
4. Commit: `git add ai-docs/ && git commit -m "docs: add insight {title}"`

### Session Recovery

1. Read MANIFEST.md: `Read ai-docs/current/{task}/MANIFEST.md`
2. Read index: `Read ai-docs/.context/insights/{task}.txt`
3. Load specific insights on demand

## Rules

1. **NO INSIGHT WITHOUT PROPER STRUCTURING FIRST** - If it's worth remembering, structure it immediately
2. **Hot vs Cold Data** - MANIFEST.md is hot (read every session), insights are cold (load on demand)
3. **Incremental Loading** - Don't load all insights at start, load when topic is mentioned

## What This Skill Is NOT

- ❌ NOT the `ai-docs/.context/scripts/` system (JSON-based, shell scripts)
- ❌ NOT for file-level code indexing (that's separate)
- ✅ This is for project insights, decisions, and state tracking

## References

See `references/` directory for:
- Template examples
- Advanced patterns
- Troubleshooting
