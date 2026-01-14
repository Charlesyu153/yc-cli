#!/bin/bash

# AI Agents 打包脚本
# 用法: ./pack-agents.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
TIMESTAMP=$(date +%Y%m%d-%H%M%S)
OUTPUT_FILE="yc-cli-agents-${TIMESTAMP}.tar.gz"

echo "打包 AI Agents..."

cd "$PROJECT_ROOT"
tar -czf "$OUTPUT_FILE" \
    agents/*.json \
    agents/README.md \
    config/standards.json \
    config/project-info.json \
    docs/AGENTS-GUIDE.md \
    docs/SKILLS-GUIDE.md \
    docs/STANDARDS.md

echo "✓ 打包完成: $OUTPUT_FILE"
