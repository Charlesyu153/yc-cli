#!/bin/bash

# AI文档管理系统 - 打包脚本
# 将ai-docs打包以便在新环境中部署
# 使用方法: ./pack-ai-docs.sh [模式]
# 模式: template (仅模板) | full (包含任务) | deploy (包含部署脚本)

set -e

# 颜色定义
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  AI文档管理系统 - 打包工具${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# 获取脚本所在目录（ai-docs目录）
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AI_DOCS_DIR="$SCRIPT_DIR"

# 获取打包模式
MODE="${1:-template}"
TIMESTAMP=$(date +%Y%m%d-%H%M%S)

case "$MODE" in
    template)
        OUTPUT_FILE="ai-docs-template-${TIMESTAMP}.tar.gz"
        echo -e "${BLUE}模式: 仅打包模板和工具（适用于新项目）${NC}"
        ;;
    full)
        OUTPUT_FILE="ai-docs-full-${TIMESTAMP}.tar.gz"
        echo -e "${BLUE}模式: 打包完整系统（包含所有任务）${NC}"
        ;;
    deploy)
        OUTPUT_FILE="ai-docs-deploy-${TIMESTAMP}.tar.gz"
        echo -e "${BLUE}模式: 打包部署版（仅模板+部署脚本）${NC}"
        ;;
    *)
        echo -e "${RED}错误: 未知模式 '$MODE'${NC}"
        echo "可用模式: template | full | deploy"
        exit 1
        ;;
esac

echo ""
echo -e "${BLUE}准备打包...${NC}"

# 创建临时目录
TEMP_DIR=$(mktemp -d)
PACK_DIR="$TEMP_DIR/ai-docs"
mkdir -p "$PACK_DIR"

# 根据模式复制文件
case "$MODE" in
    template)
        # 只打包模板和工具
        echo "复制模板和工具..."
        mkdir -p "$PACK_DIR/templates"
        mkdir -p "$PACK_DIR/current"
        mkdir -p "$PACK_DIR/archive"

        cp "$AI_DOCS_DIR/templates/"*.md "$PACK_DIR/templates/" 2>/dev/null || true
        cp "$AI_DOCS_DIR/templates/"*.sh "$PACK_DIR/templates/" 2>/dev/null || true

        # 复制文档
        cp "$AI_DOCS_DIR/QUICKREF.md" "$PACK_DIR/" 2>/dev/null || true
        cp "$AI_DOCS_DIR/README.md" "$PACK_DIR/" 2>/dev/null || true
        cp "$AI_DOCS_DIR/USAGE.md" "$PACK_DIR/" 2>/dev/null || true

        # 创建.gitkeep保持空目录
        touch "$PACK_DIR/current/.gitkeep"
        touch "$PACK_DIR/archive/.gitkeep"
        ;;

    full)
        # 打包完整系统
        echo "复制完整系统..."
        cp -r "$AI_DOCS_DIR/"* "$PACK_DIR/"

        # 排除一些不需要的文件
        rm -f "$PACK_DIR/"*.tar.gz 2>/dev/null || true
        rm -f "$PACK_DIR/pack-ai-docs.sh" 2>/dev/null || true
        ;;

    deploy)
        # 打包部署版本
        echo "复制部署文件..."
        mkdir -p "$PACK_DIR/templates"

        cp "$AI_DOCS_DIR/templates/"*.md "$PACK_DIR/templates/" 2>/dev/null || true
        cp "$AI_DOCS_DIR/templates/"*.sh "$PACK_DIR/templates/" 2>/dev/null || true

        # 复制文档
        cp "$AI_DOCS_DIR/QUICKREF.md" "$PACK_DIR/" 2>/dev/null || true
        cp "$AI_DOCS_DIR/README.md" "$PACK_DIR/" 2>/dev/null || true
        cp "$AI_DOCS_DIR/DEPLOY.md" "$PACK_DIR/" 2>/dev/null || true

        # 复制部署脚本
        cp "$AI_DOCS_DIR/setup-ai-docs.sh" "$PACK_DIR/" 2>/dev/null || true
        ;;
esac

# 打包
echo -e "${BLUE}正在压缩...${NC}"
cd "$TEMP_DIR"
tar -czf "$OUTPUT_FILE" ai-docs/

# 移动到当前目录
mv "$OUTPUT_FILE" "$AI_DOCS_DIR/"
cd "$AI_DOCS_DIR"

# 清理临时目录
rm -rf "$TEMP_DIR"

# 显示结果
FILE_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  打包完成！${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo -e "${BLUE}输出文件:${NC} $OUTPUT_FILE"
echo -e "${BLUE}文件大小:${NC} $FILE_SIZE"
echo -e "${BLUE}文件路径:${NC} $AI_DOCS_DIR/$OUTPUT_FILE"
echo ""
echo -e "${BLUE}在新环境中部署:${NC}"

case "$MODE" in
    template)
        echo "1. 传输文件到新环境"
        echo "2. tar -xzf $OUTPUT_FILE"
        echo "3. cd ai-docs/templates"
        echo "4. ./new-task.sh \"任务名\" \"负责人\""
        ;;
    full)
        echo "1. 传输文件到新环境"
        echo "2. tar -xzf $OUTPUT_FILE"
        echo "3. 根据需要清理 current/ 和 archive/ 目录"
        echo "4. cd ai-docs/templates"
        echo "5. 验证脚本权限: chmod +x *.sh"
        ;;
    deploy)
        echo "1. 传输文件到新环境"
        echo "2. tar -xzf $OUTPUT_FILE"
        echo "3. bash ai-docs/setup-ai-docs.sh [目标路径]"
        ;;
esac

echo ""
echo -e "${YELLOW}提示:${NC}"
echo "- 使用 'scp $OUTPUT_FILE user@remote:~/' 传输到远程服务器"
echo "- 或上传到云存储供团队成员下载"
