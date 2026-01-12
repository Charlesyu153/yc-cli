#!/bin/bash

# AI文档管理 - 上下文恢复文档生成脚本
# 当会话上下文超限时，使用此脚本生成新会话用的进度摘要
# 使用方法: ./export-context.sh "任务名称"

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
CURRENT_DATE=$(date +%Y-%m-%d)
TIMESTAMP=$(date +%H:%M:%S)

# 查找任务文档
TASK_DIR="$AI_DOCS_DIR/current/$TASK_NAME"
TASK_FILE="$TASK_DIR/$TASK_NAME.md"

if [ ! -f "$TASK_FILE" ]; then
    echo -e "${RED}错误: 找不到任务文档: $TASK_FILE${NC}"
    exit 1
fi

# 生成上下文恢复文档
CONTEXT_FILE="$TASK_DIR/context-resume-${CURRENT_DATE}.md"

echo -e "${BLUE}正在生成上下文恢复文档...${NC}"

cat > "$CONTEXT_FILE" << EOF
# 上下文恢复 - $TASK_NAME

> **生成时间**: $CURRENT_DATE $TIMESTAMP
> **用途**: 为新会话提供快速上下文恢复

---

## 📋 任务概览

EOF

# 提取基本信息
echo "### 基本信息" >> "$CONTEXT_FILE"
echo "" >> "$CONTEXT_FILE"
sed -n '/## 基本信息/,/---/p' "$TASK_FILE" | sed '1d;$d' >> "$CONTEXT_FILE"
echo "" >> "$CONTEXT_FILE"

# 提取任务目标
echo "### 任务目标" >> "$CONTEXT_FILE"
echo "" >> "$CONTEXT_FILE"
sed -n '/## 任务目标/,/---/p' "$TASK_FILE" | sed '1d;$d' >> "$CONTEXT_FILE"
echo "" >> "$CONTEXT_FILE"
echo "---" >> "$CONTEXT_FILE"
echo "" >> "$CONTEXT_FILE"

# 提取进度追踪
echo "## 🎯 当前进度状态" >> "$CONTEXT_FILE"
echo "" >> "$CONTEXT_FILE"
sed -n '/## 进度追踪/,/---/p' "$TASK_FILE" | sed '1d;$d' >> "$CONTEXT_FILE"
echo "" >> "$CONTEXT_FILE"
echo "---" >> "$CONTEXT_FILE"
echo "" >> "$CONTEXT_FILE"

# 提取技术方案
echo "## 🔧 技术方案摘要" >> "$CONTEXT_FILE"
echo "" >> "$CONTEXT_FILE"
sed -n '/## 技术方案/,/---/p' "$TASK_FILE" | sed '1d;$d' >> "$CONTEXT_FILE"
echo "" >> "$CONTEXT_FILE"
echo "---" >> "$CONTEXT_FILE"
echo "" >> "$CONTEXT_FILE"

# 提取关键决策
echo "## 💡 关键决策记录" >> "$CONTEXT_FILE"
echo "" >> "$CONTEXT_FILE"
sed -n '/## 关键决策/,/---/p' "$TASK_FILE" | sed '1d;$d' >> "$CONTEXT_FILE"
echo "" >> "$CONTEXT_FILE"
echo "---" >> "$CONTEXT_FILE"
echo "" >> "$CONTEXT_FILE"

# 提取变更文件清单
echo "## 📝 变更文件清单" >> "$CONTEXT_FILE"
echo "" >> "$CONTEXT_FILE"
sed -n '/## 变更文件清单/,/---/p' "$TASK_FILE" | sed '1d;$d' >> "$CONTEXT_FILE"
echo "" >> "$CONTEXT_FILE"
echo "---" >> "$CONTEXT_FILE"
echo "" >> "$CONTEXT_FILE"

# 添加新会话指引
cat >> "$CONTEXT_FILE" << 'EOF'

## 🚀 新会话快速启动指引

### 立即开始的步骤

1. **阅读"下一步"部分** - 了解当前需要做什么
2. **检查"进行中"任务** - 继续未完成的工作
3. **查看"关键决策"** - 理解已做的重要选择
4. **参考"变更文件清单"** - 定位需要修改的文件

### 常用命令参考

```bash
# 查看任务目录
cd ai-docs/current/[任务名称]

# 更新主文档
vim [任务名称].md

# 查看代码变更
git status
git diff
```

### 更新进度提示

每次完成一个步骤后，记得：
- 在主文档中更新"进度追踪"部分
- 在"变更文件清单"中记录新的文件变更
- 更新"下一步"说明

---

**原始文档**: `ai-docs/current/[任务名称]/[任务名称].md`
EOF

echo -e "${GREEN}✓ 上下文恢复文档已生成: $CONTEXT_FILE${NC}"
echo ""
echo -e "${YELLOW}使用建议:${NC}"
echo "1. 在新会话中打开此文档: $CONTEXT_FILE"
echo "2. 将此文档内容作为上下文提供给AI"
echo "3. 继续在原文档更新进度: $TASK_FILE"
echo ""
echo -e "${BLUE}新会话开始提示语示例:${NC}"
echo -e "${GREEN}\"请阅读以下上下文恢复文档，继续之前的任务...\"${NC}"
