# Filesystem Context 持久化系统

基于 filesystem-context 模式的项目上下文动态管理系统，与 ai-docs 规范深度集成。

## 核心理念

将 AI 上下文从内存转移到文件系统，实现：
1. 动态发现 - AI 按需读取相关上下文
2. 持久化存储 - 上下文永久保存，跨会话可用
3. 自动更新 - 代码变更自动触发上下文更新
4. 分层索引 - 快速定位所需信息

## 目录结构

```
ai-docs/.context/
├── README.md                    # 本文档
├── QUICKSTART.md                # 快速使用指南
├── index.json                   # 主索引文件
├── modules/                     # 模块级上下文
│   ├── agents.json             # AI Agents 模块
│   └── ai-docs.json            # AI 文档系统模块
├── tasks/                       # 任务级上下文
│   └── {任务名}.json           # 任务元数据和进度快照
├── files/                       # 文件级上下文（按需创建）
├── snapshots/                   # 定期快照
│   └── {日期}/
│       └── full-context.json
└── scripts/                     # 自动化脚本
    ├── discover-context.sh      # 发现相关上下文
    ├── update-context.sh        # 更新上下文
    ├── snapshot-context.sh      # 创建快照
    └── hooks/                   # Git hooks
        └── post-commit          # 提交后更新
```

## 上下文层级

### 1. 项目级 (index.json)
```json
{
  "project": "YC-CLI",
  "description": "AI 辅助开发工具集",
  "updated": "2026-01-14T10:00:00Z",
  "modules": ["agents", "ai-docs"],
  "active_tasks": [],
  "key_files": ["README.md", "config/project-info.json"]
}
```

### 2. 模块级 (modules/*.json)
```json
{
  "module": "agents",
  "path": "agents/",
  "description": "AI Agents 配置系统",
  "key_files": ["agents.json", "code-architect.json"],
  "workflow": ["定义角色", "配置工具", "部署使用"]
}
```

### 3. 任务级 (tasks/*.json)
```json
{
  "task": "任务名",
  "status": "in_progress|completed",
  "document": "ai-docs/current/任务名/任务名.md",
  "commits": ["abc123"]
}
```

## 动态发现机制

### 发现脚本使用

```bash
# 发现与 "agents" 相关的上下文
./scripts/discover-context.sh "agents"

# 发现与 "部署" 相关的上下文
./scripts/discover-context.sh "部署" -v
```

## 自动更新机制

### 手动更新

```bash
# 更新整个项目上下文
./scripts/update-context.sh --all

# 更新特定模块
./scripts/update-context.sh --module agents

# 更新当前任务
./scripts/update-context.sh --task "任务名"
```

### Git Hooks 集成

post-commit hook 会在每次提交后自动更新相关上下文。

## 与 ai-docs 工作流集成

### 任务创建时

```bash
# new-task.sh 自动创建任务上下文
./new-task.sh "新任务名"
# 自动生成 ai-docs/.context/tasks/新任务名.json
```

### 任务完成时

```bash
# archive-task.sh 归档任务上下文
./archive-task.sh "任务名"
# 自动标记任务状态为 completed
```

## AI 使用指南

### 会话开始时

```
请读取项目上下文：
1. ai-docs/.context/index.json - 项目概览
2. ai-docs/current/{任务名}/{任务名}.md - 当前任务详情

我需要：[具体需求]
```

### 探索代码时

```
请使用 discover-context.sh 查找与 "agents" 相关的上下文，
然后读取相关模块的上下文 JSON。
```

## 优势

| 传统方式 | Filesystem Context |
|---------|-------------------|
| 上下文在内存中，会话结束丢失 | 持久化到文件，永久保存 |
| 需要重复读取大量文件 | 按需加载精简的上下文 |
| 手动维护文档 | 自动更新上下文 |
| 难以快速定位相关信息 | 分层索引，快速发现 |

## 最佳实践

1. **会话开始时先读索引** - 快速了解项目状态
2. **使用 discover-context.sh** - 避免盲目搜索
3. **及时更新上下文** - 代码变更后立即更新
4. **定期创建快照** - 保留历史状态

## 参考资料

- AI 协作规范: ai-docs/templates/ai-guidelines.md
- 快速使用指南: ai-docs/.context/QUICKSTART.md
