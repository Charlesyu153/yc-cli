#!/bin/bash

# AI 文档系统打包脚本
# 用法: ./pack-ai-docs.sh [模式]
# 模式: template | full

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
MODE="${1:-template}"
TIMESTAMP=$(date +%Y%m%d-%H%M%S)

case "$MODE" in
    template)
        OUTPUT_FILE="yc-cli-ai-docs-template-${TIMESTAMP}.tar.gz"
        echo "打包 AI 文档系统模板..."

        cd "$PROJECT_ROOT"
        tar -czf "$OUTPUT_FILE" \
            ai-docs/templates/ \
            ai-docs/README.md \
            ai-docs/QUICKREF.md \
            config/standards.json \
            docs/AI-DOCS-GUIDE.md \
            docs/STANDARDS.md
        ;;
    full)
        OUTPUT_FILE="yc-cli-ai-docs-full-${TIMESTAMP}.tar.gz"
        echo "打包完整 AI 文档系统..."

        cd "$PROJECT_ROOT"
        tar -czf "$OUTPUT_FILE" ai-docs/
        ;;
    *)
        echo "错误: 未知模式 '$MODE'"
        echo "可用模式: template | full"
        exit 1
        ;;
esac

echo "✓ 打包完成: $OUTPUT_FILE"
