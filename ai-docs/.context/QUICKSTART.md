# Filesystem Context 快速使用指南

快速参考指南，帮助 AI 和开发者高效使用上下文系统。

## 快速开始

### AI 会话开始时

```
请读取项目上下文：
1. ai-docs/.context/index.json - 项目概览
2. ai-docs/current/{任务名}/{任务名}.md - 当前任务

我需要：[具体需求]
```

### 发现相关上下文

```bash
cd ai-docs/.context/scripts

# 搜索与 agents 相关的所有上下文
./discover-context.sh "agents"

# 仅搜索模块
./discover-context.sh "部署" -m

# 详细输出
./discover-context.sh "配置" -v
```

### 更新上下文

```bash
cd ai-docs/.context/scripts

# 更新所有上下文
./update-context.sh --all

# 更新特定模块
./update-context.sh --module agents

# 更新当前任务
./update-context.sh --task "任务名"
```

### 创建快照

```bash
cd ai-docs/.context/scripts
./snapshot-context.sh
```

## 常见场景

### 场景 1: 开始新任务

```bash
# 1. 创建任务（自动创建任务上下文）
cd ai-docs/templates
./new-task.sh "新功能开发"

# 2. AI 读取上下文
# "请读取：
#  - ai-docs/.context/index.json
#  - ai-docs/current/新功能开发/新功能开发.md"

# 3. 开发过程中，每次 commit 自动更新上下文
git commit -m "feat: 完成某功能"
```

### 场景 2: 探索项目

```bash
# 1. 搜索相关上下文
./discover-context.sh "agents" -v

# 2. 读取推荐的模块上下文
# AI: "请读取 ai-docs/.context/modules/agents.json"
```

### 场景 3: 任务完成

```bash
# 1. 填写任务完成总结
vim ai-docs/current/任务名/任务名.md

# 2. 归档任务（自动更新任务上下文为 completed）
cd ai-docs/templates
./archive-task.sh "任务名"
```

## 上下文文件结构

### 项目索引 (index.json)

```json
{
  "project": "项目名",
  "modules": [...],
  "active_tasks": [...],
  "recent_changes": [...]
}
```

用途: 快速了解项目整体结构和最新动态

### 模块上下文 (modules/*.json)

```json
{
  "module": "模块名",
  "description": "模块描述",
  "workflow": [...],
  "key_files": [...]
}
```

用途: 了解模块功能、工作流程和关键文件

### 任务上下文 (tasks/*.json)

```json
{
  "task": "任务名",
  "status": "in_progress|completed",
  "document": "任务文档路径",
  "commits": [...]
}
```

用途: 跟踪任务状态和相关提交

## AI 使用模式

### 模式 1: 渐进式加载

```
1. 读取 index.json 获取项目概览
2. 根据关键词匹配相关模块
3. 读取模块上下文了解详情
4. 按需读取文件上下文
```

优势: 节省 token，只加载需要的信息

### 模式 2: 任务驱动

```
1. 读取当前任务文档
2. 读取任务上下文获取相关 commits
3. 读取涉及模块的上下文
4. 开始工作
```

优势: 快速进入任务状态

## 维护建议

### 定期维护

```bash
# 每周创建一次快照
cd ai-docs/.context/scripts
./snapshot-context.sh

# 检查上下文是否最新
./update-context.sh --all
```

## 故障排查

### 问题: 上下文文件不存在

```bash
# 初始化所有模块上下文
cd ai-docs/.context/scripts
./update-context.sh --all
```

### 问题: 搜索无结果

```bash
# 检查索引是否最新
./update-context.sh --index

# 手动搜索文件
grep -r "关键词" ai-docs/.context/
```

## 与 ai-docs 规范的关系

Filesystem Context 是 ai-docs 系统的补充：

- **ai-docs/current/** - 任务文档（人类可读）
- **ai-docs/.context/** - 结构化上下文（机器可读）

两者配合使用：
- 任务文档记录开发过程和决策
- 上下文系统提供快速索引和发现

## 参考资料

- 完整架构: ai-docs/.context/README.md
- AI 协作规范: ai-docs/templates/ai-guidelines.md
