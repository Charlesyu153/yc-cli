#!/bin/bash
# snapshot-context.sh - 创建上下文快照

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTEXT_DIR="$(dirname "$SCRIPT_DIR")"
PROJECT_ROOT="$(dirname "$(dirname "$CONTEXT_DIR")")"

# 颜色定义
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

SNAPSHOT_DATE=$(date +%Y-%m-%d)
SNAPSHOT_DIR="$CONTEXT_DIR/snapshots/$SNAPSHOT_DATE"

echo -e "${BLUE}=== 创建上下文快照 ===${NC}\n"

# 创建快照目录
mkdir -p "$SNAPSHOT_DIR"

# 复制索引
if [ -f "$CONTEXT_DIR/index.json" ]; then
    cp "$CONTEXT_DIR/index.json" "$SNAPSHOT_DIR/"
    echo -e "${GREEN}✓${NC} 已备份项目索引"
fi

# 复制所有模块
if [ -d "$CONTEXT_DIR/modules" ]; then
    cp -r "$CONTEXT_DIR/modules" "$SNAPSHOT_DIR/"
    echo -e "${GREEN}✓${NC} 已备份模块上下文"
fi

# 复制所有任务
if [ -d "$CONTEXT_DIR/tasks" ] && [ "$(ls -A "$CONTEXT_DIR/tasks" 2>/dev/null)" ]; then
    cp -r "$CONTEXT_DIR/tasks" "$SNAPSHOT_DIR/"
    echo -e "${GREEN}✓${NC} 已备份任务上下文"
fi

# 创建快照元数据
cat > "$SNAPSHOT_DIR/snapshot.json" <<EOF
{
  "date": "$SNAPSHOT_DATE",
  "timestamp": "$(date -Iseconds)",
  "git_commit": "$(cd "$PROJECT_ROOT" && git rev-parse HEAD 2>/dev/null || echo 'unknown')",
  "git_branch": "$(cd "$PROJECT_ROOT" && git branch --show-current 2>/dev/null || echo 'unknown')"
}
EOF

echo -e "${GREEN}✓${NC} 已创建快照元数据"
echo ""
echo -e "${BLUE}快照位置:${NC} $SNAPSHOT_DIR"
echo -e "${GREEN}=== 快照完成 ===${NC}"
