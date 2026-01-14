#!/bin/bash

# AI Agents 部署脚本
# 用法: ./deploy-agents.sh [目标路径] [模式]
# 模式: safe (默认) | backup | force

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
SOURCE_DIR="$PROJECT_ROOT/agents"
TARGET_DIR="${1:-.claude/agents}"
COPY_MODE="${2:-safe}"

# 智能复制文件函数
# 参数: $1=源文件, $2=目标文件, $3=模式(safe/backup/force)
smart_copy() {
    local src="$1"
    local dst="$2"
    local mode="${3:-safe}"

    if [ ! -f "$src" ]; then
        echo "警告: 源文件不存在: $src"
        return 1
    fi

    if [ -f "$dst" ]; then
        case "$mode" in
            safe)
                echo "⊘ 跳过（已存在）: $(basename "$dst")"
                return 0
                ;;
            backup)
                local backup="${dst}.backup.$(date +%Y%m%d-%H%M%S)"
                cp "$dst" "$backup"
                echo "✓ 备份: $(basename "$dst") -> $(basename "$backup")"
                cp "$src" "$dst"
                echo "✓ 更新: $(basename "$dst")"
                ;;
            force)
                cp "$src" "$dst"
                echo "✓ 覆盖: $(basename "$dst")"
                ;;
        esac
    else
        cp "$src" "$dst"
        echo "✓ 新建: $(basename "$dst")"
    fi
}

echo "部署 AI Agents..."
echo "源目录: $SOURCE_DIR"
echo "目标目录: $TARGET_DIR"
echo "复制模式: $COPY_MODE"
echo ""

# 创建目标目录
mkdir -p "$TARGET_DIR"

# 复制配置文件（使用智能复制）
echo "复制配置文件..."
smart_copy "$SOURCE_DIR/agents.json" "$TARGET_DIR/agents.json" "$COPY_MODE"
smart_copy "$SOURCE_DIR/code-architect.json" "$TARGET_DIR/code-architect.json" "$COPY_MODE"
smart_copy "$SOURCE_DIR/code-explorer.json" "$TARGET_DIR/code-explorer.json" "$COPY_MODE"
smart_copy "$SOURCE_DIR/code-reviewer.json" "$TARGET_DIR/code-reviewer.json" "$COPY_MODE"
smart_copy "$SOURCE_DIR/README.md" "$TARGET_DIR/README.md" "$COPY_MODE"

echo ""
echo "✓ AI Agents 部署完成"
