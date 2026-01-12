# AI Agents 快速部署指南

本文档说明如何在新环境中快速部署 PCF 项目的 AI agents 配置。

## 部署方式

### 方式 1: 直接复制（推荐）

适用于：已有项目，需要添加 agents 配置

```bash
# 1. 复制整个 agents 目录到新项目
cp -r .claude/agents /path/to/new-project/.claude/

# 2. 验证文件
ls -la /path/to/new-project/.claude/agents/
```

### 方式 2: Git 克隆

适用于：从 PCF 项目克隆

```bash
# 1. 克隆项目
git clone <repository-url> new-project
cd new-project

# 2. agents 配置已包含在 .claude/agents/ 目录中
ls -la .claude/agents/
```

### 方式 3: 手动创建

适用于：需要自定义配置

```bash
# 1. 创建目录
mkdir -p .claude/agents

# 2. 复制配置文件
# 从 PCF 项目复制以下文件：
# - agents.json
# - code-architect.json
# - code-explorer.json
# - code-reviewer.json
# - README.md
# - USAGE_EXAMPLES.md
```

---

## 一键部署脚本

### 创建部署脚本

在 PCF 项目中创建部署脚本：

```bash
# 创建脚本文件
cat > .claude/agents/deploy-agents.sh << 'EOF'
#!/bin/bash

# AI Agents 一键部署脚本
# 用法: ./deploy-agents.sh [目标路径]

set -e

SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARGET_DIR="${1:-.claude/agents}"

echo "=========================================="
echo "AI Agents 配置部署工具"
echo "=========================================="
echo ""
echo "源目录: $SOURCE_DIR"
echo "目标目录: $TARGET_DIR"
echo ""

# 创建目标目录
mkdir -p "$TARGET_DIR"

# 复制配置文件
echo "正在复制配置文件..."
cp "$SOURCE_DIR/agents.json" "$TARGET_DIR/"
cp "$SOURCE_DIR/code-architect.json" "$TARGET_DIR/"
cp "$SOURCE_DIR/code-explorer.json" "$TARGET_DIR/"
cp "$SOURCE_DIR/code-reviewer.json" "$TARGET_DIR/"
cp "$SOURCE_DIR/README.md" "$TARGET_DIR/"
cp "$SOURCE_DIR/USAGE_EXAMPLES.md" "$TARGET_DIR/"
cp "$SOURCE_DIR/DEPLOY.md" "$TARGET_DIR/"

echo "✓ 配置文件复制完成"
echo ""

# 验证部署
echo "验证部署..."
if [ -f "$TARGET_DIR/agents.json" ]; then
    echo "✓ agents.json"
fi
if [ -f "$TARGET_DIR/code-architect.json" ]; then
    echo "✓ code-architect.json"
fi
if [ -f "$TARGET_DIR/code-explorer.json" ]; then
    echo "✓ code-explorer.json"
fi
if [ -f "$TARGET_DIR/code-reviewer.json" ]; then
    echo "✓ code-reviewer.json"
fi

echo ""
echo "=========================================="
echo "部署完成！"
echo "=========================================="
echo ""
echo "快速开始："
echo "  cat $TARGET_DIR/README.md"
echo ""
echo "查看使用示例："
echo "  cat $TARGET_DIR/USAGE_EXAMPLES.md"
echo ""
EOF

# 添加执行权限
chmod +x .claude/agents/deploy-agents.sh
```

### 使用部署脚本

```bash
# 部署到当前项目
cd /path/to/pcf-project
./.claude/agents/deploy-agents.sh

# 部署到其他项目
./.claude/agents/deploy-agents.sh /path/to/new-project/.claude/agents
```

---

## 配置验证

部署完成后，验证配置是否正确：

```bash
# 1. 检查文件是否存在
ls -la .claude/agents/

# 应该看到以下文件：
# - agents.json (主配置)
# - code-architect.json
# - code-explorer.json
# - code-reviewer.json
# - README.md
# - USAGE_EXAMPLES.md
# - DEPLOY.md

# 2. 验证 JSON 配置文件格式
cat .claude/agents/agents.json | python3 -m json.tool > /dev/null && echo "✓ agents.json 格式正确"

# 3. 查看快速开始指南
cat .claude/agents/README.md
```

---

## 自定义配置

### 修改项目特定信息

如果部署到非 PCF 项目，需要修改以下配置：

#### 1. 更新 agents.json

```bash
# 编辑主配置文件
vim .claude/agents/agents.json

# 修改以下字段：
# - project.name: 项目名称
# - project.description: 项目描述
# - project.tech_stack: 技术栈
# - project.context_files: 项目文档路径
```

#### 2. 更新 context_files

每个 agent 配置文件中的 `context_files` 需要指向正确的项目文档：

```json
"context_files": [
  "README.md",              // 项目主文档
  "docs/GUIDELINES.md",     // 项目规范
  "docs/ARCHITECTURE.md"    // 架构文档
]
```

#### 3. 调整项目规范

根据项目实际情况，修改 `project_standards` 部分：

```json
"project_standards": {
  "python": {
    "style": "PEP 8",
    "indent": "4 spaces",
    "max_lines": 500
  },
  "javascript": {
    "style": "ESLint",
    "indent": "2 spaces",
    "max_lines": 300
  }
}
```

---

## 常见问题

### Q1: 部署后如何测试 agents 是否工作？

```bash
# 在 Claude Code 中测试
# 输入以下提示语：
请以 code-explorer 的角色，列出当前项目的主要目录结构
```

### Q2: 如何更新已部署的配置？

```bash
# 重新运行部署脚本即可覆盖
./.claude/agents/deploy-agents.sh /path/to/target
```

### Q3: 配置文件可以放在其他位置吗？

可以，但建议统一放在 `.claude/agents/` 目录下，便于管理和版本控制。

### Q4: 如何备份现有配置？

```bash
# 备份当前配置
cp -r .claude/agents .claude/agents.backup.$(date +%Y%m%d)
```

---

## 快速部署总结

### 最简单的方式（3 步）

```bash
# 1. 复制配置目录
cp -r /path/to/pcf/.claude/agents /path/to/new-project/.claude/

# 2. 验证部署
ls -la /path/to/new-project/.claude/agents/

# 3. 查看使用说明
cat /path/to/new-project/.claude/agents/README.md
```

### 使用部署脚本（推荐）

```bash
# 1. 创建并运行部署脚本
cd /path/to/pcf
./.claude/agents/deploy-agents.sh /path/to/new-project/.claude/agents

# 2. 开始使用
# 在 Claude Code 中输入：
# 请以 code-explorer 的角色，帮我了解项目结构
```

---

**文档版本**: v1.0.0
**最后更新**: 2026-01-12
