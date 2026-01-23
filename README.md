# YC-CLI RAG

> 轻量级的知识管理系统，通过语义搜索快速找到项目中的决策、教训和洞察。

## 概述

YC-CLI RAG 是一个轻量级的知识管理系统，帮助你通过语义搜索快速找到项目中的决策、教训和洞察。

**核心价值**
- 只索引重要文档，不索引代码库
- 语义搜索：理解意图，而非关键词匹配
- 自动记录：按照模板生成结构化文档

---

## 快速开始

### 1. 安装

```bash
# 克隆项目
git clone <repository-url>
cd yc-cli

# 安装依赖
pip install -r requirements.txt
```

### 2. 配置

配置文件位于 `.rag/config.toml`，默认配置已可使用：

```toml
[server]
host = "127.0.0.1"
port = 8733

[docs]
base_dir = "ai-docs/current"
auto_index = true

[task]
name = "your-project"
```

### 3. 第一次使用

```bash
# 1. 创建项目文档目录
mkdir -p ai-docs/current/your-project/insights

# 2. 创建项目清单文件
cat > ai-docs/current/your-project/MANIFEST.md << 'EOF'
# Project Manifest: your-project

**Last Updated**: 2025-01-23
**Focus**: 项目描述

## Current State
| Aspect | Status |
|---|---|
| Phase | 初始化 |
| Blocker | 无 |
EOF

# 3. 构建索引
python -m rag.ragctl index
```

---

## 基础操作

### 搜索文档

```bash
# 搜索关键词
python -m rag.ragctl search "决策"

# 搜索并返回更多结果
python -m rag.ragctl search "性能优化" -n 10
```

### 查看状态

```bash
# 查看系统状态
python -m rag.ragctl status
```

### 重建索引

```bash
# 当文档更新后，重建索引
python -m rag.ragctl index
```

---

## 创建记录

### 使用交互式命令

```bash
python -m rag.ragctl record
```

按照提示输入：
1. Title：记录标题
2. Type：记录类型（lesson/failure/decision）
3. Conclusion：一句话总结
4. Background：背景信息
5. Key points：关键点（每行一个）

### 手动创建记录

在 `ai-docs/current/your-project/insights/` 目录下创建 `.md` 文件：

```markdown
---
title: 记录标题
type: decision
created: 2025-01-23
tags: [架构, 性能]
status: done
priority: high
related_files:
  - path: src/main.py
    locator:
      rg_pattern: "class.*Service"
    why: 此文件包含核心逻辑
---

## One-Line Conclusion
一句话总结这个决策。

## Background / Context
为什么需要做这个决策，当时的背景是什么。

## Decision / Insight
详细的决策内容或洞察。

## Key Points
- 关键点一
- 关键点二
```

---

## 文档结构

```
项目根目录/
├── ai-docs/
│   └── current/
│       └── your-project/          # 你的项目
│           ├── MANIFEST.md         # 项目清单（当前状态）
│           └── insights/           # 洞察记录
│               ├── decision_01.md
│               └── lesson_01.md
├── .rag/
│   ├── config.toml                # 配置文件
│   ├── chroma/                    # 向量索引
│   └── cache/                     # 嵌入缓存
└── rag/                           # 系统代码
```

---

## 记录类型

| 类型 | 用途 | 模板 |
|------|------|------|
| `decision` | 记录重要决策及其理由 | 决策记录模板 |
| `lesson` | 记录学到的经验教训 | 经验总结模板 |
| `failure` | 记录失败分析 | 复盘分析模板 |

---

## 高级功能

### 备份与迁移

```bash
# 导出备份
python -m rag.ragctl export backup-$(date +%Y%m%d).tar.gz

# 导入到其他项目
python -m rag.ragctl import backup-20250123.tar.gz
python -m rag.ragctl index
```

### 查看日志

```bash
# 查看最近日志
python -m rag.ragctl logs

# 实时跟踪日志
python -m rag.ragctl logs -f
```

### 启动 HTTP 服务

```bash
# 启动服务
python -m rag.ragctl start

# 在另一个终端测试
curl "http://localhost:8733/api/search?q=决策&top_k=5"
```

---

## 常见问题

### 搜索结果为空？

1. 确认已构建索引：`python -m rag.ragctl status`
2. 如果 indexed documents 为 0，运行：`python -m rag.ragctl index`
3. 确认文档在 `ai-docs/current/your-project/` 目录下

### 想搜索代码怎么办？

RAG 系统只索引文档，不索引代码。如果记录中需要引用代码：

1. 在记录中使用 `related_files` 字段指定文件路径
2. 使用 `rg_pattern` 提供稳定的搜索模式
3. 搜索到记录后，根据指针跳转到代码

### 如何在多个项目间共享？

每个项目维护自己的文档目录。如果需要共享：

1. 将通用洞察放在独立项目
2. 通过导出/导入在不同项目间迁移
3. 或者将 `docs.base_dir` 配置指向共享目录

---

## 最佳实践

1. **及时记录**：决策和教训要及时记录，避免遗忘
2. **结构清晰**：使用模板确保记录格式一致
3. **指针准确**：引用代码时使用稳定的定位方式
4. **定期备份**：定期导出备份，防止数据丢失
5. **渐进式索引**：系统会自动检测文件变化，无需频繁重建索引

---

## 获取帮助

- 查看 `skills/rag-management/SKILL.md` 了解更多技术细节
- 查看 `skills/context-insight/SKILL.md` 了解记录格式规范
- 查看日志文件 `.rag/rag.log` 排查问题

---

**文档版本**: 1.0
**最后更新**: 2025-01-23
