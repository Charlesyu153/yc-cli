# YC-CLI

AI 辅助开发工具集，包含 AI Agents 配置和文档管理系统。

## 项目简介

YC-CLI 提供两个核心子系统：

1. **AI Agents 配置** - 三个专用 AI agents（架构师、探索者、审查者）
2. **AI 文档管理系统** - 任务管理和协作规范工具

## 快速开始

### 部署

```bash
# 完整部署（安全模式，默认）
bash scripts/deploy.sh /path/to/your-project

# 备份模式：备份已存在的文件后更新
bash scripts/deploy.sh --backup /path/to/your-project

# 强制模式：直接覆盖所有文件
bash scripts/deploy.sh --force /path/to/your-project

# 仅部署 AI Agents
bash scripts/deploy.sh --agents-only /path/to/your-project

# 仅部署 AI 文档系统
bash scripts/deploy.sh --docs-only /path/to/your-project
```

**部署模式说明**：
- **安全模式**（默认）：跳过已存在的文件，保护用户修改
- **备份模式**：备份已存在的文件后覆盖
- **强制模式**：直接覆盖所有文件

**特别保护**：`ai-docs/current/` 和 `ai-docs/archive/` 目录的用户数据永远不会被覆盖

### 使用 AI Agents

在 Claude Code 中使用以下格式调用 agent：

```
请以 code-explorer 的角色，查找所有处理用户认证的代码
```

**三个 Agent 角色**：
- **code-architect**: 架构设计、技术选型、实施规划
- **code-explorer**: 代码搜索、结构分析、文档生成
- **code-reviewer**: 质量检查、问题发现、改进建议

### 使用 AI 文档系统

```bash
# 创建新任务
cd ai-docs/templates
./new-task.sh "任务名称" "负责人"

# 归档已完成任务
./archive-task.sh "任务名称"
```

### Skills 集成

优先使用可用的 Skills 来完成任务：
- muratcankoylan/Agent-Skills-for-Context-Engineering（通用开发技能）
- claude-scientific-skills（科学计算和数据分析技能）

**安装 Skills**：

```bash
# 安装通用开发技能
claude mcp install github --owner muratcankoylan --repo Agent-Skills-for-Context-Engineering

# 安装科学计算技能
claude mcp install github --owner K-Dense-AI --repo claude-scientific-skills
```

详细使用指南请查看：`docs/SKILLS-GUIDE.md`

## 核心规范

1. 禁止使用表情符号
2. 每完成一步就 commit
3. 任务完成后归档
4. 代码文件不超过 500 行

详细规范请查看：`docs/STANDARDS.md`

## 文档导航

| 文档 | 说明 |
|------|------|
| `docs/GETTING-STARTED.md` | 快速开始指南 |
| `docs/DEPLOYMENT.md` | 统一部署指南 |
| `docs/STANDARDS.md` | 项目规范 |
| `docs/AGENTS-GUIDE.md` | AI Agents 使用指南 |
| `docs/AI-DOCS-GUIDE.md` | AI 文档系统使用指南 |
| `docs/SKILLS-GUIDE.md` | Skills 使用指南 |

## 目录结构

```
yc-cli/
├── README.md                  # 本文档
├── config/                    # 统一配置
│   ├── standards.json        # 项目规范配置
│   └── project-info.json     # 项目元信息
├── docs/                      # 统一文档中心
│   ├── GETTING-STARTED.md
│   ├── DEPLOYMENT.md
│   ├── STANDARDS.md
│   ├── AGENTS-GUIDE.md
│   ├── AI-DOCS-GUIDE.md
│   └── SKILLS-GUIDE.md
├── scripts/                   # 统一部署脚本
│   ├── deploy.sh
│   ├── deploy-agents.sh
│   ├── deploy-ai-docs.sh
│   ├── pack-agents.sh
│   └── pack-ai-docs.sh
├── agents/                    # AI Agents 配置
│   ├── agents.json
│   ├── code-architect.json
│   ├── code-explorer.json
│   ├── code-reviewer.json
│   └── README.md
└── ai-docs/                   # AI 文档管理系统
    ├── current/               # 进行中的任务
    ├── archive/               # 已完成任务
    ├── templates/             # 模板和工具
    ├── README.md
    └── QUICKREF.md
```

## 许可证

MIT License

## 贡献

欢迎提交 Issue 和 Pull Request。

---

**版本**: v1.0.0
**创建日期**: 2026-01-12
**更新日期**: 2026-01-14
