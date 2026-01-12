# YC-CLI

AI 辅助开发工具集，包含 AI Agents 配置和文档管理系统。

## 项目简介

YC-CLI 是一个用于提升 AI 辅助开发效率的工具集，包含：

1. **AI Agents 配置** - 三个专用 AI agents（架构师、探索者、审查者）
2. **AI 文档管理系统** - 任务管理和协作规范工具

## 快速开始

### 部署 AI Agents

```bash
# 方式 1: 使用部署脚本
./agents/deploy-agents.sh /path/to/your-project/.claude/agents

# 方式 2: 直接复制
cp -r agents /path/to/your-project/.claude/
```

### 部署 AI 文档管理系统

```bash
# 复制到项目根目录
cp -r ai-docs /path/to/your-project/
```

## 目录结构

```
yc-cli/
├── agents/              # AI Agents 配置
│   ├── code-architect.json
│   ├── code-explorer.json
│   ├── code-reviewer.json
│   ├── agents.json
│   ├── deploy-agents.sh
│   ├── README.md
│   ├── USAGE_EXAMPLES.md
│   └── DEPLOY.md
├── ai-docs/             # AI 文档管理系统
│   ├── templates/
│   ├── current/
│   ├── archive/
│   ├── README.md
│   └── *.md
└── README.md            # 本文档
```

## AI Agents 使用

### 三个 Agent 角色

1. **code-architect** - 代码架构师
   - 架构设计、技术选型、实施规划

2. **code-explorer** - 代码探索者
   - 代码搜索、结构分析、文档生成

3. **code-reviewer** - 代码审查者
   - 质量检查、问题发现、改进建议

### 基本用法

在 Claude Code 中使用以下格式调用 agent：

```
请以 code-explorer 的角色，查找所有处理用户认证的代码
```

### Feature-Dev 工作流

完整的功能开发流程：

1. **探索** (code-explorer) - 理解现有代码
2. **设计** (code-architect) - 设计实现方案
3. **实现** - 编写代码
4. **审查** (code-reviewer) - 质量检查
5. **提交** - 更新文档并提交

详细使用示例请查看：`agents/USAGE_EXAMPLES.md`

## AI 文档管理系统

### 核心功能

- 任务文档管理（创建、归档）
- 上下文恢复（长期任务连续性）
- AI 协作规范（统一的工作流程）
- 工具脚本（自动化任务管理）

### 快速使用

```bash
# 创建新任务
cd ai-docs/templates
./new-task.sh "任务名称" "负责人"

# 归档已完成任务
./archive-task.sh "任务名称"

# 导出上下文
./export-context.sh "任务名称"
```

### 核心规范

1. 禁止使用表情符号
2. 每完成一步就 commit
3. 任务完成后归档
4. 代码文件不超过 500 行

详细说明请查看：`ai-docs/README.md`

## 文档导航

| 文档 | 说明 |
|------|------|
| `agents/README.md` | AI Agents 快速开始 |
| `agents/USAGE_EXAMPLES.md` | 详细使用示例 |
| `agents/DEPLOY.md` | 部署指南 |
| `ai-docs/README.md` | AI 文档管理系统说明 |
| `ai-docs/QUICKREF.md` | 快速参考 |
| `ai-docs/USAGE.md` | 完整使用指南 |

## 许可证

MIT License

## 贡献

欢迎提交 Issue 和 Pull Request。

---

**版本**: v1.0.0
**创建日期**: 2026-01-12
