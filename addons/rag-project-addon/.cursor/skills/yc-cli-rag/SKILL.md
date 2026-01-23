---
name: yc-cli-rag
description: Use this skill when the user asks to search past decisions/lessons/insights or manage project RAG on localhost:8733. Provides a consistent workflow using yc-cli RAG HTTP API first, with a wrapper script fallback.
---

# yc-cli RAG (localhost:8733)

## When to use

Use this when the user asks to:
- 搜索历史决策/教训/洞察（ai-docs）
- 查“之前怎么做的/为什么这么做”
- 更新或重建 RAG 索引
- 生成新的 decision/lesson/failure 记录

## Default workflow

1) Prefer HTTP API (fast, tool-agnostic):
- Health: `GET /health`
- Status: `GET /api/status`
- Search: `GET /api/search?q=<query>&top_k=<n>`
- Index: `POST /api/index`
- Record: `POST /api/record` (JSON body)

2) If HTTP is not reachable:
- Ask user to run `./scripts/yc_rag.sh start` in the project root (after `mamba activate rag_cli`)
- Then retry HTTP

## Output rules

- Return: top hits with `file_path` and a short 1-2 line summary per hit.
- Do NOT paste large document bodies into the chat.
- If user needs details, open the referenced file(s) and cite only the necessary lines.
