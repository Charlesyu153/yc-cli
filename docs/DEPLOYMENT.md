# YC-CLI 部署指南

本文档说明如何在新环境中部署 YC-CLI 的 AI Agents 和 AI 文档管理系统。

---

## 部署方案总览

| 方案 | 命令 | 适用场景 |
|------|------|----------|
| 安全模式（默认） | `./scripts/deploy.sh /path/to/project` | 更新现有安装，保护用户修改 |
| 备份模式 | `./scripts/deploy.sh --backup /path/to/project` | 更新并保留备份 |
| 强制模式 | `./scripts/deploy.sh --force /path/to/project` | 完全重置配置 |
| 仅 Agents | `./scripts/deploy.sh --agents-only /path/to/project` | 只需要 AI 角色配置 |
| 仅 AI-Docs | `./scripts/deploy.sh --docs-only /path/to/project` | 只需要任务管理系统 |
| 打包部署 | 先打包再传输解压 | 跨服务器部署 |

**部署模式说明**：
- **安全模式**（默认）：跳过已存在的文件，保护用户修改
- **备份模式**：备份已存在的文件后覆盖，备份文件名格式：`文件名.backup.YYYYMMDD-HHMMSS`
- **强制模式**：直接覆盖所有文件，用于完全重置配置

**特别保护**：`ai-docs/current/` 和 `ai-docs/archive/` 目录的用户数据永远不会被覆盖

---

## 一、统一部署（推荐）

### 1.1 安全模式部署（默认）

适用于更新现有安装，保护用户修改的配置和文档。

```bash
# 在 yc-cli 项目目录执行
bash scripts/deploy.sh /path/to/your-project

# 或显式指定安全模式
bash scripts/deploy.sh --safe /path/to/your-project

# 验证部署
ls -la /path/to/your-project/.claude/agents/
ls -la /path/to/your-project/ai-docs/
ls -la /path/to/your-project/config/
ls -la /path/to/your-project/docs/
```

**行为**：
- 如果文件已存在，跳过并提示
- 只复制不存在的文件
- 完全保护用户的修改

**输出示例**：
```
✓ 新建: agents.json
⊘ 跳过（已存在）: code-architect.json
✓ 新建: README.md
⊘ 保留: current/ 目录（包含用户任务数据）
⊘ 保留: archive/ 目录（包含用户历史数据）
```

### 1.2 备份模式部署

适用于需要更新配置，但希望保留备份以便回滚。

```bash
bash scripts/deploy.sh --backup /path/to/your-project
```

**行为**：
- 如果文件已存在，先备份（添加 `.backup.YYYYMMDD-HHMMSS` 后缀）
- 然后复制新文件
- 可以通过备份文件回滚

**输出示例**：
```
✓ 备份: agents.json -> agents.json.backup.20260114-105030
✓ 更新: agents.json
✓ 新建: README.md
```

**回滚方法**：
```bash
# 恢复备份文件
cp agents.json.backup.20260114-105030 agents.json
```

### 1.3 强制模式部署

适用于完全重置配置，覆盖所有文件。

```bash
bash scripts/deploy.sh --force /path/to/your-project
```

**行为**：
- 直接覆盖所有文件
- 不创建备份
- 用户修改会丢失（除了 current/ 和 archive/）

**警告**：此模式会覆盖用户的配置修改，请谨慎使用！

**输出示例**：
```
✓ 覆盖: agents.json
✓ 覆盖: code-architect.json
✓ 新建: README.md
```

### 1.4 完整部署

部署 AI Agents + AI 文档系统 + 配置 + 文档：

```bash
# 安全模式（默认）
bash scripts/deploy.sh /path/to/your-project

# 备份模式
bash scripts/deploy.sh --backup /path/to/your-project

# 强制模式
bash scripts/deploy.sh --force /path/to/your-project
```

**部署内容**：
- `.claude/agents/`: AI Agents 配置文件
- `ai-docs/`: AI 文档管理系统（保护 current/ 和 archive/）
- `config/`: 共享配置文件
- `docs/`: 完整文档

### 1.2 仅部署 AI Agents

```bash
bash scripts/deploy.sh --agents-only /path/to/your-project

# 验证部署
ls -la /path/to/your-project/.claude/agents/
```

**部署内容**：
- `agents.json`: 主配置文件
- `code-architect.json`: 架构师配置
- `code-explorer.json`: 探索者配置
- `code-reviewer.json`: 审查者配置
- `README.md`: 简要说明

### 1.3 仅部署 AI 文档系统

```bash
bash scripts/deploy.sh --docs-only /path/to/your-project

# 验证部署
ls -la /path/to/your-project/ai-docs/
```

**部署内容**：
- `current/`: 进行中的任务目录
- `archive/`: 已完成任务归档目录
- `templates/`: 模板和工具脚本
- `README.md`: 简要说明
- `QUICKREF.md`: 快速参考

---

## 二、打包部署

适用于跨服务器部署或离线部署。

### 2.1 打包 AI Agents

```bash
# 在 yc-cli 项目目录执行
bash scripts/pack-agents.sh

# 生成文件: yc-cli-agents-YYYYMMDD-HHMMSS.tar.gz
```

**打包内容**：
- agents/*.json
- agents/README.md
- config/standards.json
- config/project-info.json
- docs/AGENTS-GUIDE.md
- docs/SKILLS-GUIDE.md
- docs/STANDARDS.md

### 2.2 打包 AI 文档系统

```bash
# 打包模板（推荐）
bash scripts/pack-ai-docs.sh template

# 打包完整系统（包含现有任务）
bash scripts/pack-ai-docs.sh full

# 生成文件: yc-cli-ai-docs-template-YYYYMMDD-HHMMSS.tar.gz
```

### 2.3 解压和部署

```bash
# 传输到目标服务器
scp yc-cli-*.tar.gz user@server:/path/to/project/

# 在目标服务器解压
cd /path/to/project
tar -xzf yc-cli-agents-*.tar.gz
tar -xzf yc-cli-ai-docs-*.tar.gz
```

---

## 三、手动部署

适用于需要深度定制的场景。

### 3.1 部署 AI Agents

```bash
# 1. 创建目录
mkdir -p /path/to/project/.claude/agents

# 2. 复制配置文件
cp agents/*.json /path/to/project/.claude/agents/
cp agents/README.md /path/to/project/.claude/agents/

# 3. 复制共享配置
mkdir -p /path/to/project/config
cp config/*.json /path/to/project/config/
```

### 3.2 部署 AI 文档系统

```bash
# 1. 创建目录结构
mkdir -p /path/to/project/ai-docs/{current,archive,templates}

# 2. 复制模板和脚本
cp ai-docs/templates/*.md /path/to/project/ai-docs/templates/
cp ai-docs/templates/*.sh /path/to/project/ai-docs/templates/
chmod +x /path/to/project/ai-docs/templates/*.sh

# 3. 复制文档
cp ai-docs/README.md /path/to/project/ai-docs/
cp ai-docs/QUICKREF.md /path/to/project/ai-docs/

# 4. 创建 .gitkeep
touch /path/to/project/ai-docs/current/.gitkeep
touch /path/to/project/ai-docs/archive/.gitkeep
```

---

## 四、自定义配置

### 4.1 修改项目规范

编辑 `config/standards.json`：

```json
{
  "core_rules": {
    "max_lines": {
      "limit": 500  // 修改行数限制
    }
  },
  "coding_standards": {
    "python": {
      "indent": "4 spaces"  // 修改缩进规范
    }
  }
}
```

### 4.2 自定义 Agent 配置

编辑 `agents/code-*.json`：

```json
{
  "system_prompt": "你是一个专业的代码架构师...",  // 修改提示词
  "capabilities": ["architecture-design", "..."],  // 修改能力
  "tools": ["Glob", "Grep", "Read", "..."]  // 修改可用工具
}
```

### 4.3 自定义文档模板

编辑 `ai-docs/templates/task-template.md`：

```markdown
# 任务名称

## 任务目标
[自定义模板内容]
```

---

## 五、验证部署

### 5.1 验证 AI Agents

```bash
# 检查配置文件
ls -la .claude/agents/
cat .claude/agents/agents.json

# 验证 JSON 格式
jq . .claude/agents/agents.json
```

### 5.2 验证 AI 文档系统

```bash
# 检查目录结构
tree ai-docs/

# 测试工具脚本
cd ai-docs/templates
./new-task.sh "测试任务" "测试人员"
ls -la ../current/测试任务/
```

### 5.3 验证共享配置

```bash
# 检查配置文件
ls -la config/
jq . config/standards.json
jq . config/project-info.json
```

---

## 六、常见问题

### Q1: 部署脚本执行失败

**原因**：可能缺少执行权限

**解决**：
```bash
chmod +x scripts/*.sh
```

### Q2: 配置文件引用错误

**原因**：路径不正确

**解决**：
- 确保 `config/` 目录在项目根目录
- 检查 `agents.json` 中的 `config_reference` 路径

### Q3: 工具脚本无法执行

**原因**：缺少执行权限或路径错误

**解决**：
```bash
chmod +x ai-docs/templates/*.sh
cd ai-docs/templates
./new-task.sh "任务名" "负责人"
```

### Q4: 如何更新已部署的配置

**方案 1**：重新部署（覆盖）
```bash
bash scripts/deploy.sh --full /path/to/project
```

**方案 2**：手动更新
```bash
# 只更新配置文件
cp config/*.json /path/to/project/config/
cp agents/*.json /path/to/project/.claude/agents/
```

---

## 七、团队部署方案

### 7.1 个人新项目

```bash
# 直接完整部署
bash scripts/deploy.sh /path/to/new-project
```

### 7.2 团队推广

```bash
# 1. 打包模板
bash scripts/pack-agents.sh
bash scripts/pack-ai-docs.sh template

# 2. 分发给团队成员
# 团队成员解压后使用
```

### 7.3 现有项目迁移

```bash
# 直接复制整个 yc-cli 目录
cp -r yc-cli /path/to/existing-project/

# 或使用部署脚本
cd yc-cli
bash scripts/deploy.sh /path/to/existing-project
```

---

## 八、部署检查清单

### 部署前
- [ ] 确认目标路径正确
- [ ] 确认有写入权限
- [ ] 备份现有配置（如有）

### 部署中
- [ ] 执行部署脚本
- [ ] 观察输出信息
- [ ] 确认无错误提示

### 部署后
- [ ] 验证目录结构
- [ ] 验证配置文件格式
- [ ] 测试工具脚本
- [ ] 阅读快速开始文档

---

**相关文档**：
- 快速开始：`docs/GETTING-STARTED.md`
- 项目规范：`docs/STANDARDS.md`
- Agents 指南：`docs/AGENTS-GUIDE.md`
- AI-Docs 指南：`docs/AI-DOCS-GUIDE.md`

**版本**：v1.0.0
**更新日期**：2026-01-14
