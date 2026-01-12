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
