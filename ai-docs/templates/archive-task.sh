#!/bin/bash

# AI文档管理 - 任务归档脚本
# 使用方法: ./archive-task.sh "任务名称"

set -e

# 颜色定义
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 获取脚本所在目录的父目录（ai-docs目录）
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AI_DOCS_DIR="$(dirname "$SCRIPT_DIR")"

# 检查参数
if [ $# -lt 1 ]; then
    echo -e "${RED}错误: 缺少任务名称${NC}"
    echo "使用方法: $0 \"任务名称\""
    echo "示例: $0 \"添加用户认证功能\""
    exit 1
fi

TASK_NAME="$1"
CURRENT_DIR="$AI_DOCS_DIR/current/$TASK_NAME"
ARCHIVE_DIR="$AI_DOCS_DIR/archive"
TARGET_DIR="$ARCHIVE_DIR/$TASK_NAME"

# 检查任务目录是否存在
if [ ! -d "$CURRENT_DIR" ]; then
    echo -e "${RED}错误: 任务目录不存在: $CURRENT_DIR${NC}"
    echo "当前可用的任务："
    ls "$AI_DOCS_DIR/current/" 2>/dev/null || echo "  (无)"
    exit 1
fi

# 检查任务文档是否存在
TASK_FILE="$CURRENT_DIR/$TASK_NAME.md"
if [ ! -f "$TASK_FILE" ]; then
    echo -e "${YELLOW}警告: 找不到任务文档: $TASK_FILE${NC}"
    read -p "是否继续归档? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "已取消"
        exit 1
    fi
fi

# 检查完成总结是否已填写
if [ -f "$TASK_FILE" ]; then
    if grep -q "\[简述最终实现的内容\]" "$TASK_FILE" || \
       grep -q "\[记录值得注意的经验，供未来参考\]" "$TASK_FILE"; then
        echo -e "${YELLOW}警告: 任务文档的完成总结似乎未填写${NC}"
        echo "建议先填写完成总结再归档"
        read -p "是否继续归档? (y/N): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "已取消，请先填写完成总结"
            exit 1
        fi
    fi
fi

# 检查状态是否为"已完成"
if [ -f "$TASK_FILE" ]; then
    if ! grep -q "已完成" "$TASK_FILE" | head -5; then
        echo -e "${YELLOW}警告: 任务状态可能不是'已完成'${NC}"
        read -p "是否继续归档? (y/N): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "已取消"
            exit 1
        fi
    fi
fi

# 检查目标目录是否已存在
if [ -d "$TARGET_DIR" ]; then
    echo -e "${RED}错误: 归档目录已存在: $TARGET_DIR${NC}"
    read -p "是否覆盖? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm -rf "$TARGET_DIR"
    else
        echo "已取消"
        exit 1
    fi
fi

# 创建归档目录（如果不存在）
mkdir -p "$ARCHIVE_DIR"

# 执行归档
echo -e "${BLUE}正在归档任务: $TASK_NAME${NC}"
mv "$CURRENT_DIR" "$TARGET_DIR"

# 统计文件数量
FILE_COUNT=$(find "$TARGET_DIR" -type f | wc -l)
DIR_SIZE=$(du -sh "$TARGET_DIR" | cut -f1)

echo -e "${GREEN}归档成功！${NC}"
echo ""
echo "任务名称: $TASK_NAME"
echo "原路径: $CURRENT_DIR"
echo "归档路径: $TARGET_DIR"
echo "文件数量: $FILE_COUNT"
echo "目录大小: $DIR_SIZE"
echo ""
echo -e "${BLUE}归档内容:${NC}"
tree "$TARGET_DIR" -L 2 2>/dev/null || ls -lh "$TARGET_DIR"
echo ""
echo -e "${GREEN}下一步:${NC}"
echo "1. 查看归档: ls $TARGET_DIR"
echo "2. 如需要，可以将归档提交到git:"
echo "   git add $TARGET_DIR"
echo "   git commit -m \"docs: 归档任务 $TASK_NAME\""
