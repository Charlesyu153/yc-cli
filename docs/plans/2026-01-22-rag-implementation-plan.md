# RAG 系统实现计划

**关联设计**: [2026-01-22-rag-system-design.md](./2026-01-22-rag-system-design.md)

---

## 阶段划分

```
┌─────────────────────────────────────────────────────────────┐
│  阶段 1        │  阶段 2        │  阶段 3        │  阶段 4   │
│  核心功能      │  AI 集成       │  自动记录      │  可选功能 │
│  (必须)        │  (必须)        │  (推荐)        │  (未来)   │
├─────────────────────────────────────────────────────────────┤
│  ┌─────────┐  │  ┌─────────┐  │  ┌─────────┐  │  ┌───────┐ │
│  │ 索引构建 │  │  │ Skill更新│  │  │ 记录器  │  │  │ WebUI │ │
│  │ 语义搜索 │  │  │ 会话钩子 │  │  │ 自动触发 │  │  │ 定制化│ │
│  │ 文件监听 │  │  └─────────┘  │  └─────────┘  │  └───────┘ │
│  └─────────┘  │                 │               │            │
└─────────────────────────────────────────────────────────────┘
```

---

## 阶段 1：核心功能

### 1.1 环境准备
```bash
# 创建虚拟环境
python -m venv venv
source venv/bin/activate

# 安装依赖
pip install fastapi uvicorn chromadb sentence-transformers watchdog
```

### 1.2 核心模块

| 文件 | 功能 | 代码行数 (估计) |
|------|------|----------------|
| `rag/config.py` | 配置管理 | ~30 |
| `rag/embeddings.py` | 嵌入生成 | ~40 |
| `rag/vector_store.py` | ChromaDB 封装 | ~80 |
| `rag/indexer.py` | 文件索引构建 | ~60 |
| `rag/server.py` | FastAPI 服务 | ~100 |
| `rag/file_watcher.py` | 文件监听 | ~50 |

### 1.3 API 端点

```
POST /api/index          # 重建整个索引
GET  /api/search?q=xx&n=5 # 语义搜索
GET  /api/status          # 索引状态
POST /api/record          # 保存 AI 生成的记录
```

### 1.4 验收标准
- [ ] 服务能正常启动 (`python -m rag.server`)
- [ ] 能索引现有 ai-docs 目录
- [ ] 能搜索并返回相关文件路径
- [ ] 搜索结果按相关性排序
- [ ] 文件修改后索引自动更新

---

## 阶段 2：AI 集成

### 2.1 Skill 更新

更新 `skills/context-insight/SKILL.md`，添加 RAG 调用流程：

```markdown
## 会话开始时

1. GET /api/status - 检查索引状态
2. 如需要，POST /api/index - 更新索引
3. GET /api/search?q={当前任务} - 获取相关上下文
4. Read {返回的文件路径} - 读取内容
```

### 2.2 新增 RAG Skill

创建 `skills/rag-retrieval/SKILL.md`：

```markdown
---
name: rag-retrieval
description: 自动从向量数据库检索相关项目上下文
---

## 用法

当需要项目上下文时，调用：
GET http://localhost:8733/api/search?q={关键词}&top_k=5

然后读取返回的文件路径。
```

### 2.3 验收标准
- [ ] Skill 文档更新完成
- [ ] AI 能正确调用 API
- [ ] 搜索���果能正确被 AI 使用
- [ ] 会话开始时自动检索相关上下文

### 2.4 新增自动记录 Skill

创建 `skills/auto-record/SKILL.md`：

```markdown
---
name: auto-record
description: 任务完成时自动生成简洁记录并保存
---

## 触发条件

完成以下工作时，生成记录：
- 修复了 bug
- 做了技术选择/决策
- 发现了性能问题/陷阱
- 理解了复杂逻辑

## 记录格式

```markdown
---
title: {简短标题}
type: lesson|failure|decision
tags: {标签}
created: {时间}
---

## 结论
{一句话}

## 背景
{为什么}

## 关键点
1. ...
2. ...
```

## 保存流程

1. 生成 1-3 条记录
2. 显示预览，询问用户确认
3. 确认后调用 POST /api/record 保存
```

---

## 阶段 3：自动记录

### 3.1 记录器模块

| 文件 | 功能 | 代码行数 (估计) |
|------|------|----------------|
| `rag/recorder.py` | 记录生成与保存 | ~80 |
| `rag/templates/` | 记录模板 | ~30 |

### 3.2 API 端点

```
POST /api/record           # 保存记录
  body: {title, type, tags, content, task}
  → 创建 insight 文件 + 更新索引
```

### 3.3 验收标准
- [ ] 能按模板创建 insight 文件
- [ ] 创建后自动更新索引
- [ ] AI 能正确调用记录 API
- [ ] 记录可被搜索检索到

---

## 阶段 4：可选功能

| 功能 | 优先级 | 说明 |
|------|--------|------|
| Web UI | 低 | 浏览器查看索引、手动搜索 |
| 多模型支持 | 低 | 支持切换不同嵌入模型 |
| 跨项目索引 | 低 | 同时索引多个项目的 ai-docs |
| 导入导出 | 低 | 索引的备份与恢复 |

---

## 实现顺序建议

```
1. config.py          # 配置先定好
2. embeddings.py      # 依赖配置
3. vector_store.py    # 依赖 embeddings
4. indexer.py         # 依赖 vector_store
5. file_watcher.py    # 文件监听
6. server.py          # 依赖上面所有
7. recorder.py        # 记录器
8. Skill 更新         # 与 server 并行
```

---

## 风险与应对

| 风险 | 应对 |
|------|------|
| ChromaDB 兼容性问题 | 使用 docker 封装 |
| 嵌入模型下载慢 | 首次运行提供镜像/缓存 |
| 端口冲突 | 配置文件可修改端口 |
| AI 生成记录质量 | 预览确认机制，用户可拒绝 |
| 跨平台迁移 | .rag/ 和 ai-docs/ 独立，复制即用 |

