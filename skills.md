---
name: skills-index
description: Central index for this repo's skills. Read this FIRST when user mentions "skill" or "context skill".
---

# Skills Index

**重要**: 用户说 "context skill" 时，只用 `skills/context-insight/SKILL.md`，不要与其它 JSON/脚本类上下文系统混淆。

---

## 可用 Skills

| Skill | 文件 | 何时调用 |
|-------|------|----------|
| **context-insight** | `skills/context-insight/SKILL.md` | 记录洞察、教训、项目状态管理 |
| **rag-management** | `skills/rag-management/SKILL.md` | RAG 系统管理、语义搜索、自动记录生成 |

---

## context-insight 调用方式

### 1. 读取 Skill 定义

用户说 "用 context skill" 或 "记录这个教训" 时：
```bash
Read skills/context-insight/SKILL.md
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

## rag-management 调用方式

### 1. 读取 Skill 定义

用户说 "用 rag skill" 或 "搜索文档" 时：
```bash
Read skills/rag-management/SKILL.md
```

### 2. CLI 工具

```bash
# 状态
ragctl status

# 重建索引
ragctl index

# 搜索
ragctl search "decision"

# 创建记录
ragctl record

# 导出/导入
ragctl export backup.tar.gz
ragctl import backup.tar.gz
```

### 3. Python API

```python
from rag.indexer import index_all
from rag.vector_store import VectorStore

# 索引���有文件
results = index_all()

# 搜索文档
results = VectorStore.search("decision", top_k=5)
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
用户说 "context skill" → Read skills/context-insight/SKILL.md
用户说 "记录这个"     → 按照 templates 创建 insight + 更新索引
用户说 "更新上下文"   → 默认走 context-insight（索引/按需加载），必要时再讨论是否需要 RAG 服务
用户说 "rag skill"    → Read skills/rag-management/SKILL.md
用户说 "搜索文档"     → 使用 RAG 语义搜索
```

---

**更新日期**: 2026-01-23
