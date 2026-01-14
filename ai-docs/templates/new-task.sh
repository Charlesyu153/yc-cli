#!/bin/bash

# AI文档管理 - 新任务初始化脚本
# 使用方法: ./new-task.sh "任务名称" ["负责人"]

set -e

# 颜色定义
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# 获取脚本所在目录的父目录（ai-docs目录）
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AI_DOCS_DIR="$(dirname "$SCRIPT_DIR")"

# 检查参数
if [ $# -lt 1 ]; then
    echo -e "${RED}错误: 缺少任务名称${NC}"
    echo "使用方法: $0 \"任务名称\" [\"负责人\"]"
    echo "示例: $0 \"添加用户认证功能\" \"张三\""
    exit 1
fi

TASK_NAME="$1"
OWNER="${2:-待指定}"
CURRENT_DATE=$(date +%Y-%m-%d)

# 创建任务目录
TASK_DIR="$AI_DOCS_DIR/current/$TASK_NAME"

if [ -d "$TASK_DIR" ]; then
    echo -e "${RED}错误: 任务目录已存在: $TASK_DIR${NC}"
    exit 1
fi

echo -e "${BLUE}正在创建任务: $TASK_NAME${NC}"
mkdir -p "$TASK_DIR"

# 复制并填充模板
TEMPLATE_FILE="$AI_DOCS_DIR/templates/task-template.md"
TARGET_FILE="$TASK_DIR/$TASK_NAME.md"

if [ ! -f "$TEMPLATE_FILE" ]; then
    echo -e "${RED}错误: 模板文件不存在: $TEMPLATE_FILE${NC}"
    exit 1
fi

# 使用sed替换模板中的占位符
sed -e "s/\[任务名称\]/$TASK_NAME/g" \
    -e "s/\[姓名\]/$OWNER/g" \
    -e "s/\[YYYY-MM-DD\]/$CURRENT_DATE/g" \
    -e "s/进行中 \/ 已完成/进行中/g" \
    "$TEMPLATE_FILE" > "$TARGET_FILE"

echo -e "${GREEN}✓ 任务文档已创建: $TARGET_FILE${NC}"
echo -e "${BLUE}任务目录: $TASK_DIR${NC}"

# 创建任务上下文
CONTEXT_SCRIPT="$AI_DOCS_DIR/.context/scripts/update-context.sh"
if [ -f "$CONTEXT_SCRIPT" ]; then
    echo -e "${BLUE}正在创建任务上下文...${NC}"
    "$CONTEXT_SCRIPT" --task "$TASK_NAME" 2>/dev/null && \
        echo -e "${GREEN}✓ 任务上下文已创建${NC}" || \
        echo -e "${RED}! 任务上下文创建失败（可忽略）${NC}"
fi

echo ""
echo -e "${GREEN}下一步:${NC}"
echo "1. 编辑任务文档: $TARGET_FILE"
echo "2. 填写任务目标和技术方案"
echo "3. 开始开发并更新进度追踪"
