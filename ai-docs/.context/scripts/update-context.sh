#!/bin/bash
# update-context.sh - 更新项目上下文

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTEXT_DIR="$(dirname "$SCRIPT_DIR")"
PROJECT_ROOT="$(dirname "$(dirname "$CONTEXT_DIR")")"

# 颜色定义
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

usage() {
    echo "用法: $0 [选项]"
    echo ""
    echo "选项:"
    echo "  --all                 更新所有上下文"
    echo "  --module <名称>       更新指定模块"
    echo "  --task <任务名>       更新指定任务"
    echo "  --index               仅更新项目索引"
    echo ""
    echo "示例:"
    echo "  $0 --all                              # 更新所有"
    echo "  $0 --module agents                    # 更新 agents 模块"
    echo "  $0 --task 新功能开发                   # 更新任务"
    exit 1
}

# 更新项目索引
update_index() {
    echo -e "${BLUE}[更新项目索引]${NC}"

    local index_file="$CONTEXT_DIR/index.json"
    local temp_file=$(mktemp)

    # 获取最新 commit 信息
    cd "$PROJECT_ROOT"
    local latest_commit=$(git log -1 --format="%h" 2>/dev/null || echo "unknown")
    local commit_date=$(git log -1 --format="%ci" 2>/dev/null || date -Iseconds)
    local commit_msg=$(git log -1 --format="%s" 2>/dev/null || echo "")

    # 读取现有索引
    if [ -f "$index_file" ]; then
        # 更新时间戳和最新变更
        jq --arg date "$(date -Iseconds)" \
           --arg commit "$latest_commit" \
           --arg commit_date "$commit_date" \
           --arg msg "$commit_msg" \
           '.updated = $date |
            .recent_changes = [{
                date: $commit_date,
                commit: $commit,
                summary: $msg
            }] + .recent_changes[:4]' \
           "$index_file" > "$temp_file"

        mv "$temp_file" "$index_file"
        echo -e "  ${GREEN}✓${NC} 项目索引已更新"
    else
        echo -e "  ${RED}✗${NC} 项目索引不存在"
        return 1
    fi
}

# 更新模块上下文
update_module() {
    local module_name="$1"
    echo -e "${BLUE}[更新模块: $module_name]${NC}"

    local module_file="$CONTEXT_DIR/modules/${module_name}.json"

    if [ ! -f "$module_file" ]; then
        echo -e "  ${RED}✗${NC} 模块上下文不存在: $module_file"
        return 1
    fi

    # 更新最后修改时间
    local temp_file=$(mktemp)
    jq --arg date "$(date +%Y-%m-%d)" \
       '.last_modified = $date' \
       "$module_file" > "$temp_file"

    mv "$temp_file" "$module_file"
    echo -e "  ${GREEN}✓${NC} 模块 $module_name 已更新"
}

# 更新任务上下文
update_task() {
    local task_name="$1"
    echo -e "${BLUE}[更新任务: $task_name]${NC}"

    local task_doc="$PROJECT_ROOT/current/$task_name/$task_name.md"
    local task_context="$CONTEXT_DIR/tasks/${task_name}.json"

    if [ ! -f "$task_doc" ]; then
        # 检查是否已归档
        local archived_doc="$PROJECT_ROOT/archive/$task_name/$task_name.md"
        if [ -f "$archived_doc" ]; then
            echo -e "  ${YELLOW}!${NC} 任务已归档"
            task_doc="$archived_doc"
        else
            echo -e "  ${RED}✗${NC} 任务文档不存在"
            return 1
        fi
    fi

    # 创建或更新任务上下文
    mkdir -p "$CONTEXT_DIR/tasks"

    # 提取任务状态
    local status="in_progress"
    if [ -f "$archived_doc" ]; then
        status="completed"
    fi

    # 获取最近的 commits
    cd "$PROJECT_ROOT"
    local commits=$(git log --all --grep="$task_name" --format="%h" 2>/dev/null | head -5 | jq -R . | jq -s . || echo "[]")

    # 创建/更新任务上下文
    cat > "$task_context" <<EOF
{
  "task": "$task_name",
  "status": "$status",
  "document": "ai-docs/current/$task_name/$task_name.md",
  "updated": "$(date -Iseconds)",
  "commits": $commits
}
EOF

    echo -e "  ${GREEN}✓${NC} 任务上下文已更新"
}

# 更新所有上下文
update_all() {
    echo -e "${BLUE}=== 更新所有上下文 ===${NC}\n"

    # 更新项目索引
    update_index
    echo ""

    # 更新所有模块
    if [ -d "$CONTEXT_DIR/modules" ]; then
        for module_file in "$CONTEXT_DIR/modules"/*.json; do
            if [ -f "$module_file" ]; then
                module_name=$(basename "$module_file" .json)
                update_module "$module_name"
            fi
        done
        echo ""
    fi

    # 更新所有当前任务
    if [ -d "$PROJECT_ROOT/current" ]; then
        for task_dir in "$PROJECT_ROOT/current"/*; do
            if [ -d "$task_dir" ]; then
                task_name=$(basename "$task_dir")
                update_task "$task_name"
            fi
        done
        echo ""
    fi

    echo -e "${GREEN}=== 更新完成 ===${NC}"
}

# 参数解析
if [ $# -eq 0 ]; then
    usage
fi

case "$1" in
    --all)
        update_all
        ;;
    --module)
        if [ -z "$2" ]; then
            echo "错误: 必须指定模块名称"
            usage
        fi
        update_module "$2"
        ;;
    --task)
        if [ -z "$2" ]; then
            echo "错误: 必须指定任务名称"
            usage
        fi
        update_task "$2"
        ;;
    --index)
        update_index
        ;;
    -h|--help)
        usage
        ;;
    *)
        echo "错误: 未知选项 $1"
        usage
        ;;
esac
