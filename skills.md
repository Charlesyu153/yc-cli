---
name: skills-index
description: Central index for all ~/.claude skills. Read this FIRST when user mentions "skill" or "context skill".
---

# ~/.claude Skills Index

**重要**: 用户说 "context skill" 时，只用 `~/.claude/context-insight.md`，不要与 `ai-docs/.context/` 的 JSON 系统混淆。

---

## 可用 Skills

| Skill | 文件 | 何时调用 |
|-------|------|----------|
| **context-insight** | `~/.claude/context-insight.md` | 记录洞察、教训、项目状态管理 |

---

## context-insight 调用方式

### 1. 读取 Skill 定义

用户说 "用 context skill" 或 "记录这个教训" 时：
```bash
Read ~/.claude/context-insight.md
```

### 2. 会话开始时读取 Hot Data

```bash
# 1. 读取 MANIFEST.md（当前状态）
Read ai-docs/current/{任务}/MANIFEST.md

# 2. 读取洞察索引（按需）
Read ai-docs/.context/insights/{任务}.txt
```

### 3. 创建新 Insight

```bash
# 模板位置
ai-docs/templates/insight-template.md

# 索引位置
ai-docs/.context/insights/{任务}.txt

# 详细记录位置
ai-docs/current/{任务}/insights/{insight_id}.md
```

---

## 禁止混淆

| 禁止 | 替代 |
|------|------|
| `ai-docs/.context/scripts/discover-context.sh` | 用户说 "context skill" 时不要用这个 |
| `ai-docs/.context/scripts/update-context.sh` | 用户说 "context skill" 时不要用这个 |
| JSON 格式的文件上下文 | 用户说 "context skill" 时不要读这些 |

---

## 简单规则

```
用户说 "context skill" → Read ~/.claude/context-insight.md
用户说 "记录这个"  → 按照 context-insight.md 格式创建 insight
用户说 "更新上下文"  → 默认使用 ~/.claude/context-insight.md，如果用户有特殊需求，再询问是用 ~/.claude/context-insight.md 还是 ai-docs/.context/
```

---

**更新日期**: 2026-01-21
