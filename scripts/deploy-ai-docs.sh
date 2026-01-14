#!/bin/bash

# AI 文档系统部署脚本
# 用法: ./deploy-ai-docs.sh [目标路径] [模式]
# 模式: safe (默认) | backup | force

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
SOURCE_DIR="$PROJECT_ROOT/ai-docs"
TARGET_DIR="${1:-ai-docs}"
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

echo "部署 AI 文档系统..."
echo "源目录: $SOURCE_DIR"
echo "目标目录: $TARGET_DIR"
echo "复制模式: $COPY_MODE"
echo ""

# 创建目录结构（但不复制 current 和 archive 的内容）
mkdir -p "$TARGET_DIR"/{current,archive,templates}

# 保护用户数据：只创建 .gitkeep（如果不存在）
echo "保护用户数据目录..."
if [ ! -f "$TARGET_DIR/current/.gitkeep" ]; then
    touch "$TARGET_DIR/current/.gitkeep"
    echo "✓ 创建: current/.gitkeep"
else
    echo "⊘ 保留: current/ 目录（包含用户任务数据）"
fi

if [ ! -f "$TARGET_DIR/archive/.gitkeep" ]; then
    touch "$TARGET_DIR/archive/.gitkeep"
    echo "✓ 创建: archive/.gitkeep"
else
    echo "⊘ 保留: archive/ 目录（包含用户历史数据）"
fi
echo ""

# 复制模板文件（使用智能复制）
echo "复制模板文件..."
for file in "$SOURCE_DIR/templates/"*.md; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        smart_copy "$file" "$TARGET_DIR/templates/$filename" "$COPY_MODE"
    fi
done
echo ""

# 复制脚本文件（使用智能复制）
echo "复制工具脚本..."
for file in "$SOURCE_DIR/templates/"*.sh; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        smart_copy "$file" "$TARGET_DIR/templates/$filename" "$COPY_MODE"
        chmod +x "$TARGET_DIR/templates/$filename"
    fi
done
echo ""

# 复制文档（使用智能复制）
echo "复制文档..."
smart_copy "$SOURCE_DIR/README.md" "$TARGET_DIR/README.md" "$COPY_MODE"
smart_copy "$SOURCE_DIR/QUICKREF.md" "$TARGET_DIR/QUICKREF.md" "$COPY_MODE"

echo ""
echo "✓ AI 文档系统部署完成"
echo ""
echo "注意: current/ 和 archive/ 目录的用户数据已被完全保护"
