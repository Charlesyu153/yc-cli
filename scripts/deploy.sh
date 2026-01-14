#!/bin/bash

# YC-CLI 统一部署脚本
# 用法: ./deploy.sh [选项] [目标路径]

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# 显示帮助信息
show_help() {
    cat << EOF
YC-CLI 部署工具

用法:
    ./deploy.sh [选项] [目标路径]

选项:
    -a, --agents-only       仅部署 AI Agents
    -d, --docs-only         仅部署 AI 文档系统
    -f, --full              完整部署（默认）
    -s, --safe              安全模式：跳过已存在的文件（默认）
    -b, --backup            备份模式：备份已存在的文件后覆盖
    --force                 强制模式：直接覆盖所有文件
    -h, --help              显示帮助信息

示例:
    ./deploy.sh /path/to/project                    # 完整部署（安全模式）
    ./deploy.sh --agents-only /path/to/project      # 仅部署 Agents
    ./deploy.sh --backup /path/to/project           # 备份模式部署
    ./deploy.sh --force /path/to/project            # 强制覆盖部署

EOF
}

# 智能复制文件函数（用于 config 和 docs）
smart_copy() {
    local src="$1"
    local dst="$2"
    local mode="$3"

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

# 解析参数
MODE="full"
COPY_MODE="safe"
TARGET_DIR=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -a|--agents-only)
            MODE="agents"
            shift
            ;;
        -d|--docs-only)
            MODE="docs"
            shift
            ;;
        -f|--full)
            MODE="full"
            shift
            ;;
        -s|--safe)
            COPY_MODE="safe"
            shift
            ;;
        -b|--backup)
            COPY_MODE="backup"
            shift
            ;;
        --force)
            COPY_MODE="force"
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            TARGET_DIR="$1"
            shift
            ;;
    esac
done

# 验证目标路径
if [ -z "$TARGET_DIR" ]; then
    echo "错误: 请指定目标路径"
    show_help
    exit 1
fi

echo "=========================================="
echo "YC-CLI 部署工具"
echo "=========================================="
echo ""
echo "部署模式: $MODE"
echo "复制模式: $COPY_MODE"
echo "目标路径: $TARGET_DIR"
echo ""

# 执行部署
case $MODE in
    agents)
        bash "$SCRIPT_DIR/deploy-agents.sh" "$TARGET_DIR/.claude/agents" "$COPY_MODE"
        ;;
    docs)
        bash "$SCRIPT_DIR/deploy-ai-docs.sh" "$TARGET_DIR/ai-docs" "$COPY_MODE"
        ;;
    full)
        echo "开始完整部署..."
        echo ""

        # 部署 Agents
        bash "$SCRIPT_DIR/deploy-agents.sh" "$TARGET_DIR/.claude/agents" "$COPY_MODE"
        echo ""

        # 部署 AI-Docs
        bash "$SCRIPT_DIR/deploy-ai-docs.sh" "$TARGET_DIR/ai-docs" "$COPY_MODE"
        echo ""

        # 复制共享配置
        echo "复制共享配置..."
        mkdir -p "$TARGET_DIR/config"
        for file in "$PROJECT_ROOT/config/"*.json; do
            if [ -f "$file" ]; then
                filename=$(basename "$file")
                smart_copy "$file" "$TARGET_DIR/config/$filename" "$COPY_MODE"
            fi
        done
        echo ""

        # 复制文档
        echo "复制文档..."
        mkdir -p "$TARGET_DIR/docs"
        for file in "$PROJECT_ROOT/docs/"*.md; do
            if [ -f "$file" ]; then
                filename=$(basename "$file")
                smart_copy "$file" "$TARGET_DIR/docs/$filename" "$COPY_MODE"
            fi
        done
        echo ""

        echo "✓ 完整部署完成"
        ;;
esac

echo ""
echo "=========================================="
echo "部署完成！"
echo "=========================================="
echo ""
echo "快速开始: cat $TARGET_DIR/docs/GETTING-STARTED.md"
echo ""
