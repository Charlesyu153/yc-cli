# YC-CLI 快速开始

5 分钟快速部署和使用 YC-CLI。

## 一、快速部署

### 方式 1: 完整部署（推荐）

```bash
# 安全模式（默认）：跳过已存在的文件
bash scripts/deploy.sh /path/to/your-project

# 备份模式：备份已存在的文件后更新
bash scripts/deploy.sh --backup /path/to/your-project

# 强制模式：直接覆盖所有文件
bash scripts/deploy.sh --force /path/to/your-project

# 验证部署
ls /path/to/your-project/.claude/agents/
ls /path/to/your-project/ai-docs/
ls /path/to/your-project/config/
ls /path/to/your-project/docs/
```

**部署模式说明**：
- **安全模式**（默认）：跳过已存在的文件，保护用户修改
- **备份模式**：备份已存在的文件后覆盖，备份文件名格式：`文件名.backup.YYYYMMDD-HHMMSS`
- **强制模式**：直接覆盖所有文件，用于完全重置配置

**特别保护**：`ai-docs/current/` 和 `ai-docs/archive/` 目录的用户数据永远不会被覆盖

### 方式 2: 仅部署 AI Agents

```bash
bash scripts/deploy.sh --agents-only /path/to/your-project
```

### 方式 3: 仅部署 AI 文档系统

```bash
bash scripts/deploy.sh --docs-only /path/to/your-project
```

## 二、基本使用

### 使用 AI Agents

在 Claude Code 中使用以下格式调用 agent：

```
请以 code-explorer 的角色，查找所有处理用户认证的代码
```

**三个 Agent 角色：**
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

# 导出上下文
./export-context.sh "任务名称"
```

### 使用 Skills

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

## 三、核心规范

1. 禁止使用表情符号
2. 每完成一步就 commit
3. 任务完成后归档
4. 代码文件不超过 500 行

详细规范请查看：`docs/STANDARDS.md`

## 四、下一步

- 详细 Agents 指南: `docs/AGENTS-GUIDE.md`
- 详细 AI-Docs 指南: `docs/AI-DOCS-GUIDE.md`
- Skills 使用指南: `docs/SKILLS-GUIDE.md`
- 部署指南: `docs/DEPLOYMENT.md`
- 项目规范: `docs/STANDARDS.md`

## 五、常见问题

### 如何选择部署方式？

- 新项目：完整部署
- 只需要 AI 角色：仅部署 Agents
- 只需要任务管理：仅部署 AI-Docs

### 如何更新配置？

修改 `config/standards.json` 即可，所有子系统会自动引用。

### 如何自定义 Agent？

编辑 `agents/code-*.json` 配置文件，修改 system_prompt 和 capabilities。

---

**版本**: v1.0.0
**更新日期**: 2026-01-14
