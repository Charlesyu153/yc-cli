#!/bin/bash
# discover-context.sh - 根据关键词发现相关上下文

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTEXT_DIR="$(dirname "$SCRIPT_DIR")"
PROJECT_ROOT="$(dirname "$(dirname "$CONTEXT_DIR")")"

# 颜色定义
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

usage() {
    echo "用法: $0 <关键词> [选项]"
    echo ""
    echo "选项:"
    echo "  -m, --modules-only    仅搜索模块"
    echo "  -t, --tasks-only      仅搜索任务"
    echo "  -v, --verbose         详细输出"
    echo ""
    echo "示例:"
    echo "  $0 agents              # 搜索所有与 agents 相关的上下文"
    echo "  $0 部署 -m             # 仅搜索模块"
    echo "  $0 配置 -v             # 详细输出"
    exit 1
}

# 参数解析
KEYWORD=""
MODULES_ONLY=false
TASKS_ONLY=false
VERBOSE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -m|--modules-only)
            MODULES_ONLY=true
            shift
            ;;
        -t|--tasks-only)
            TASKS_ONLY=true
            shift
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            if [ -z "$KEYWORD" ]; then
                KEYWORD="$1"
            else
                echo "错误: 未知参数 $1"
                usage
            fi
            shift
            ;;
    esac
done

if [ -z "$KEYWORD" ]; then
    echo "错误: 必须提供关键词"
    usage
fi

echo -e "${BLUE}=== 搜索关键词: $KEYWORD ===${NC}\n"

# 搜索项目索引
echo -e "${GREEN}[项目概览]${NC}"
if [ -f "$CONTEXT_DIR/index.json" ]; then
    grep -i "$KEYWORD" "$CONTEXT_DIR/index.json" 2>/dev/null || echo "  未在项目索引中找到匹配"
else
    echo "  警告: 项目索引文件不存在"
fi
echo ""

# 搜索模块
if [ "$TASKS_ONLY" = false ]; then
    echo -e "${GREEN}[相关模块]${NC}"
    found_modules=false
    if [ -d "$CONTEXT_DIR/modules" ]; then
        for module_file in "$CONTEXT_DIR/modules"/*.json; do
            if [ -f "$module_file" ]; then
                if grep -qi "$KEYWORD" "$module_file"; then
                    module_name=$(basename "$module_file" .json)
                    echo -e "  ${YELLOW}模块:${NC} $module_name"

                    if [ "$VERBOSE" = true ]; then
                        echo "    文件: $module_file"
                        description=$(jq -r '.description // "无描述"' "$module_file" 2>/dev/null || echo "无法解析")
                        echo "    描述: $description"
                        echo "    匹配内容:"
                        grep -i "$KEYWORD" "$module_file" | head -3 | sed 's/^/      /'
                    fi

                    found_modules=true
                fi
            fi
        done
    fi

    if [ "$found_modules" = false ]; then
        echo "  未找到相关模块"
    fi
    echo ""
fi

# 搜索任务
if [ "$MODULES_ONLY" = false ]; then
    echo -e "${GREEN}[相关任务]${NC}"
    found_tasks=false

    # 搜索任务上下文
    if [ -d "$CONTEXT_DIR/tasks" ]; then
        for task_file in "$CONTEXT_DIR/tasks"/*.json; do
            if [ -f "$task_file" ]; then
                if grep -qi "$KEYWORD" "$task_file"; then
                    task_name=$(basename "$task_file" .json)
                    echo -e "  ${YELLOW}任务:${NC} $task_name"

                    if [ "$VERBOSE" = true ]; then
                        status=$(jq -r '.status // "unknown"' "$task_file" 2>/dev/null || echo "unknown")
                        echo "    状态: $status"
                        echo "    文件: $task_file"
                    fi

                    found_tasks=true
                fi
            fi
        done
    fi

    # 搜索当前任务文档
    if [ -d "$PROJECT_ROOT/current" ]; then
        for task_dir in "$PROJECT_ROOT/current"/*; do
            if [ -d "$task_dir" ]; then
                task_name=$(basename "$task_dir")
                task_doc="$task_dir/$task_name.md"

                if [ -f "$task_doc" ] && grep -qi "$KEYWORD" "$task_doc"; then
                    echo -e "  ${YELLOW}当前任务:${NC} $task_name"

                    if [ "$VERBOSE" = true ]; then
                        echo "    文档: $task_doc"
                        echo "    匹配内容:"
                        grep -i "$KEYWORD" "$task_doc" | head -3 | sed 's/^/      /'
                    fi

                    found_tasks=true
                fi
            fi
        done
    fi

    if [ "$found_tasks" = false ]; then
        echo "  未找到相关任务"
    fi
    echo ""
fi

# 推荐阅读
echo -e "${GREEN}[推荐阅读]${NC}"
echo "  1. 项目索引: ai-docs/.context/index.json"

# 根据关键词推荐模块
case "$KEYWORD" in
    *agent*|*Agent*|*角色*|*architect*|*explorer*|*reviewer*)
        echo "  2. Agents 模块: ai-docs/.context/modules/agents.json"
        echo "  3. Agents 指南: docs/AGENTS-GUIDE.md"
        ;;
    *doc*|*文档*|*任务*|*task*|*归档*)
        echo "  2. AI-docs 模块: ai-docs/.context/modules/ai-docs.json"
        echo "  3. AI-docs 指南: docs/AI-DOCS-GUIDE.md"
        ;;
    *部署*|*deploy*|*安装*)
        echo "  2. 部署指南: docs/DEPLOYMENT.md"
        echo "  3. 部署脚本: scripts/deploy.sh"
        ;;
    *)
        echo "  2. 查看所有模块: ls ai-docs/.context/modules/"
        ;;
esac

echo ""
echo -e "${BLUE}=== 搜索完成 ===${NC}"
