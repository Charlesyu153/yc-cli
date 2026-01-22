---
name: skills-index
description: Central index for all ~/.claude/skills. Read this FIRST when user mentions "skill", "context skill", or asks what skills are available. This is the entry point for discovering and invoking the correct skill.
---

# ~/.claude Skills Index

## Available Skills

| Skill | Location | When to Use |
|-------|----------|-------------|
| **context-insight** | `skills/context-insight/` | User says "context skill", needs to record insights/lessons, manage project state |

## Quick Decision Tree

```
User mentions "skill"?
├─ "context skill" → Use skills/context-insight/
├─ "record this" → Use skills/context-insight/
├─ "lesson learned" → Use skills/context-insight/
└─ Other → Ask user to specify
```

## Skill Structure

All skills follow this format:
```
skills/{skill-name}/
├── SKILL.md           # Required: Main skill file
└── references/        # Optional: Detailed docs
    └── *.md
```

## Usage

1. **User says "context skill"**:
   ```bash
   Read ~/.claude/skills/context-insight/SKILL.md
   ```

2. **User wants to record something**:
   - Follow context-insight format
   - Create insight in `ai-docs/current/{task}/insights/`
   - Update index in `ai-docs/.context/insights/{task}.txt`

## Important Notes

- This system (`~/.claude/skills/`) is **separate** from `ai-docs/.context/`
- When user says "context skill", default to `~/.claude/skills/context-insight/`
- If unclear, ask: "Use ~/.claude context-insight or ai-docs/.context/?"
