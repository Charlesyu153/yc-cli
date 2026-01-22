# AI 会话启动模板

**关联**: [2026-01-22-rag-system-design.md](./2026-01-22-rag-system-design.md)

---

## 会话启动检查清单

每次新会话开始时，AI 按以下顺序执行：

```
1. [ ] 检查服务状态
2. [ ] 检查索引状态
3. [ ] 获取当前任务
4. [ ] 检索相关上下文
5. [ ] 读取关键文件
6. [ ] 开始工作
```

---

## 自动执行流程

### Step 1: 检查服务

```bash
curl http://127.0.0.1:8733/health
```

| 结果 | 操作 |
|------|------|
| 200 OK | 继续 |
| 连接失败 | 提示用户: `服务未启动，运行 ragctl start` |

---

### Step 2: 检查索引

```bash
curl http://127.0.0.1:8733/api/status
```

| 结果 | 操作 |
|------|------|
| count > 0 | 继续 |
| count = 0 | 提示用户: `索引为空，运行 ragctl index` |

---

### Step 3: 获取当前任务

从配置读取：
```bash
# .rag/config.toml
[task]
name = "当前任务名"
```

如果用户指定了任务，使用用户的任务名。

---

### Step 4: 检索相关上下文

```bash
curl "http://127.0.0.1:8733/api/search?q={当前任务}&top_k=3"
```

获取与当前任务最相关的 3 个文档。

---

### Step 5: 读取关键文件

对 Step 4 返回的每个文件：

```python
for result in search_results:
    Read(result["file_path"])
```

注意：
- 启动阶段只读取 `MANIFEST.md` / `insights/*.md`（不要直接把脚本/全代码库加载进来）。
- 当需要代码上下文时：先读 insight，再根据 `related_files` 的 locator 用 `rg -n -C 3` 按需定位到脚本片段。

---

### Step 6: 向用户报告

```markdown
## 会话就绪

- 服务: ✓ 运行中
- 索引: ✓ {count} 个文档
- 当前任务: {task_name}
- 已加载上下文:
  1. {doc1_title}
  2. {doc2_title}
  3. {doc3_title}

准备就绪，请告诉我需要做什么。
```

---

## AI Skill 模板

创建 `skills/session-startup/SKILL.md`：

```markdown
---
name: session-startup
description: AI 会话启动时自动执行的初始化流程
---

## 执行步骤

### 1. 检查服务

```bash
curl http://127.0.0.1:8733/health
```

- 失败 → 提示: `服务未启动，请先运行 ragctl start`
- 成功 → 继续

### 2. 检查索引

```bash
curl http://127.0.0.1:8733/api/status
```

- count = 0 → 提示: `索引为空，请先运行 ragctl index`
- count > 0 → 继续

### 3. 获取当前任务

从 `.rag/config.toml` 读取 `task.name`，或询问用户。

### 4. 检索上下文

```bash
curl "http://127.0.0.1:8733/api/search?q={任务名}&top_k=3"
```

### 5. 读取文件

```python
for result in search_results:
    Read(result["file_path"])
```

说明：
- 这里默认返回的是 `MANIFEST.md` / `insights/*.md`。
- 如果后续需要代码定位：读取 insight 后，跟随 `related_files` 做按需 snippet 加载。

### 6. 报告状态

向用户展示加载的上下文，确认准备就绪。

---

## 快捷方式

用户也可以直接说：

- "加载项目上下文" → 执行上述流程
- "搜索 {关键词}" → 直接执行 Step 4
- "重建索引" → 执行 `ragctl index`
```

---

## 用户快捷指令

| 用户说 | AI 做 |
|--------|-------|
| "加载项目上下文" | 执行完整启动流程 |
| "搜索 {关键词}" | 调用 /api/search |
| "相关上下文" | 基于当前对话搜索 |
| "状态检查" | 显示服务 + 索引状态 |
