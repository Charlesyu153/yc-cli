# AI文档管理系统 - 使用指南

> 一个专为AI辅助开发设计的任务文档管理系统
> 帮助在长期项目中保持上下文连续性，特别是处理会话超限问题

---

## 📁 目录结构

```
ai-docs/
├── current/                    # 当前进行中的任务
│   └── {需求名}/               # 文件夹名 = 需求名
│       ├── {需求名}.md         # 主文档 (禁止使用 README.md)
│       ├── context-resume-*.md # 上下文恢复文档（自动生成）
│       └── *.md                # 补充文档（可选）
├── archive/                    # 已完成的任务
│   └── {需求名}/               # 完成后从 current 移动到此
└── templates/                  # 模板和工具
    ├── task-template.md        # 任务文档模板
    ├── context-resume-template.md  # 上下文恢复模板
    ├── ai-guidelines.md        # AI协作预设指令
    ├── new-task.sh             # 新建任务脚本
    ├── export-context.sh       # 导出上下文脚本
    └── archive-task.sh         # 归档任务脚本
```

---

## 🚀 快速开始

### 1. 创建新任务

使用自动化脚本（推荐）：

```bash
cd ai-docs/templates
./new-task.sh "任务名称" "负责人"
```

示例：
```bash
./new-task.sh "CAF亚型分析优化" "张三"
```

手动创建：
```bash
# 创建任务目录
mkdir -p ai-docs/current/任务名称

# 复制并重命名模板
cp ai-docs/templates/task-template.md ai-docs/current/任务名称/任务名称.md

# 编辑文档
vim ai-docs/current/任务名称/任务名称.md
```

### 2. 开发过程中更新文档

随着任务推进，持续更新以下部分：

- **进度追踪**: 标记已完成、进行中、待完成的步骤
- **下一步**: 明确写出下一步要做什么
- **变更文件清单**: 记录新增/修改的文件
- **关键决策**: 记录重要的技术选择和原因

### 3. 任务完成后归档

使用自动化脚本（推荐）：

```bash
cd ai-docs/templates
./archive-task.sh "任务名称"
```

手动归档：

```bash
# 填写"完成总结"部分
vim ai-docs/current/任务名称/任务名称.md

# 移动到归档目录
mv ai-docs/current/任务名称 ai-docs/archive/
```

---

## AI协作预设指令

为确保AI严格遵循项目规范，在开始新会话时请提供预设指令。

### 会话开始模板

```
请遵循以下规范进行协作：

1. 阅读AI协作指令: ai-docs/templates/ai-guidelines.md
2. 阅读任务文档: ai-docs/current/{任务名}/{任务名}.md

关键规范：
- 禁止使用表情符号
- 每完成一个步骤后：更新进度文档 + git commit
- 所有任务完成后：填写完成总结 + 归档到archive目录
- 所有代码文件不超过500行，超过即模块化重构

当前需要：[具体描述需求]
```

### 核心规范说明

**1. 禁止使用表情符号**
- 所有文档和代码中不使用emoji
- 保持专业、纯文本的沟通风格

**2. 完成步骤后的标准流程**
```bash
# a. 更新进度文档
vim ai-docs/current/{任务名}/{任务名}.md

# b. Git提交
git add .
git commit -m "feat: 完成XXX功能

- 具体完成的工作
- 修改的关键文件

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>"
```

**3. 任务完成后的归档流程**
```bash
# a. 填写完成总结
vim ai-docs/current/{任务名}/{任务名}.md

# b. 最终提交
git add ai-docs/current/{任务名}/{任务名}.md
git commit -m "docs: 完成{任务名}任务总结"

# c. 归档任务
cd ai-docs/templates
./archive-task.sh "{任务名}"

# d. 提交归档
git add ai-docs/archive/{任务名}
git commit -m "docs: 归档任务 {任务名}"
```

**4. 代码质量规范**
- 单个代码文件或脚本不超过500行
- 超过500行时必须进行模块化重构
- 重构原则：
  - 按功能拆分成多个模块
  - 使用清晰的函数/类结构
  - 保持代码可读性和可维护性
  - 在任务文档中记录重构决策

详细规范请查看：`ai-docs/templates/ai-guidelines.md`

---

## 🔄 处理上下文超限

当AI会话上下文即将超限时，使用以下方案：

### 方案A: 自动导出（推荐）

```bash
cd ai-docs/templates
./export-context.sh "任务名称"
```

这将自动生成 `context-resume-{日期}.md` 文件，包含：
- 任务概览
- 当前进度
- 技术方案摘要
- 关键决策
- 变更文件清单
- 新会话启动指引

### 方案B: 手动创建

```bash
# 复制模板
cp ai-docs/templates/context-resume-template.md \
   ai-docs/current/任务名称/context-resume-$(date +%Y-%m-%d).md

# 填写内容
vim ai-docs/current/任务名称/context-resume-$(date +%Y-%m-%d).md
```

### 在新会话中使用

1. **启动新会话**
2. **提供上下文恢复文档**：
   ```
   请阅读以下上下文恢复文档，继续之前的开发任务：

   [粘贴 context-resume-*.md 的内容]

   我当前需要：[具体需求]
   ```

3. **继续开发**，AI会根据文档快速理解当前状态

---

## 💡 使用建议

### 文档更新频率

- **高频更新**：进度追踪、下一步
- **中频更新**：变更文件清单
- **低频更新**：技术方案、关键决策

### 与AI协作的最佳实践

1. **任务开始时**
   - 让AI阅读任务文档
   - 明确当前要做的步骤

2. **开发过程中**
   - 完成步骤后立即更新进度
   - 记录重要的技术决策
   - 使用 `path/to/file.kt:line` 格式引用代码

3. **会话切换时**
   - 导出上下文恢复文档
   - 在新会话中提供恢复文档
   - 继续更新主任务文档

4. **任务完成时**
   - 填写完成总结
   - 记录经验教训
   - 归档到 archive 目录

### 代码引用规范

在文档中引用代码时使用以下格式：
- 文件引用：`path/to/file.R`
- 类/函数引用：`ClassName` or `functionName()`
- 精确定位：`path/to/file.R:42` (第42行)

这种格式可以让AI和其他开发者快速定位代码位置。

---

## 📊 文档模板说明

### task-template.md

完整的任务文档模板，包含：
- 基本信息
- 任务目标
- 进度追踪
- 技术方案
- 关键决策
- 变更文件清单
- 完成总结

### context-resume-template.md

精简的上下文恢复模板，专注于：
- 快速概览当前状态
- 明确下一步行动
- 关键技术信息
- 新会话启动指引

---

## 🛠️ 脚本工具说明

### new-task.sh

**功能**: 自动创建新任务结构

**使用**:
```bash
./new-task.sh "任务名称" ["负责人"]
```

**输出**:
- 创建任务目录
- 生成填充好基本信息的任务文档
- 显示下一步提示

### export-context.sh

**功能**: 从主任务文档导出上下文恢复文档

**使用**:
```bash
./export-context.sh "任务名称"
```

**输出**:
- 自动提取关键信息
- 生成带时间戳的恢复文档
- 添加新会话启动指引

### archive-task.sh

**功能**: 归档已完成的任务到archive目录

**使用**:
```bash
./archive-task.sh "任务名称"
```

**输出**:
- 检查完成总结是否已填写
- 将任务从current移动到archive
- 显示归档统计信息
- 提供git提交命令建议

---

## 📝 常见场景示例

### 场景1: 启动新的数据分析任务

```bash
# 1. 创建任务
./new-task.sh "单细胞数据CAF亚型分析" "研究员A"

# 2. 编辑任务文档，填写目标和方案
vim ai-docs/current/单细胞数据CAF亚型分析/单细胞数据CAF亚型分析.md

# 3. 开始与AI对话
# "请阅读 ai-docs/current/单细胞数据CAF亚型分析/单细胞数据CAF亚型分析.md，
#  帮我完成第一步数据预处理"
```

### 场景2: 会话上下文超限

```bash
# 1. 导出当前进度
./export-context.sh "单细胞数据CAF亚型分析"

# 2. 打开新会话，提供恢复文档
# "请阅读以下上下文恢复文档继续任务：
#  [粘贴 context-resume-2026-01-12.md 的内容]"

# 3. 继续开发，同时更新主文档
```

### 场景3: 任务完成归档

```bash
# 1. 填写完成总结部分
vim ai-docs/current/单细胞数据CAF亚型分析/单细胞数据CAF亚型分析.md

# 2. 提交最终文档
git add ai-docs/current/单细胞数据CAF亚型分析/单细胞数据CAF亚型分析.md
git commit -m "docs: 完成单细胞数据CAF亚型分析任务总结"

# 3. 归档任务
cd ai-docs/templates
./archive-task.sh "单细胞数据CAF亚型分析"

# 4. 提交归档
git add ai-docs/archive/单细胞数据CAF亚型分析
git commit -m "docs: 归档任务 单细胞数据CAF亚型分析"
```

### 场景4: 团队协作交接

```bash
# 1. 更新完成总结部分
vim ai-docs/current/任务名称/任务名称.md

# 2. 导出最新上下文
./export-context.sh "任务名称"

# 3. 将上下文文档发送给接手人员
```

---

## ⚠️ 注意事项

### 禁止使用 README.md

不要使用 `README.md` 作为任务文档名，因为：
- GitHub会自动渲染README，可能造成混淆
- 任务名称更具描述性
- 保持文件名与文件夹名一致

### 保持文档同步

- 主任务文档是单一数据源
- 上下文恢复文档是快照，用完即可丢弃
- 所有更新都应反映到主文档

### 版本控制

建议将 `ai-docs/` 目录纳入 Git 版本控制：

```bash
git add ai-docs/
git commit -m "docs: 更新任务进度"
```

---

## 🎯 总结

这个文档管理系统的核心价值：

1. **持续上下文** - 长期任务不丢失进度
2. **快速恢复** - 新会话快速进入状态
3. **知识沉淀** - 记录决策和经验
4. **团队协作** - 便于交接和协作
5. **规范流程** - AI协作遵循统一规范

### 关键文档速查

- **使用指南**: `ai-docs/USAGE.md` (本文档)
- **AI协作规范**: `ai-docs/templates/ai-guidelines.md`
- **任务模板**: `ai-docs/templates/task-template.md`
- **上下文恢复模板**: `ai-docs/templates/context-resume-template.md`

### 工具脚本速查

- **创建任务**: `./new-task.sh "任务名" "负责人"`
- **导出上下文**: `./export-context.sh "任务名"`
- **归档任务**: `./archive-task.sh "任务名"`

### 核心规范速查

1. **无表情** - 禁止使用表情符号
2. **步骤完成即commit** - 每完成一个步骤后更新文档并git commit
3. **任务完成即归档** - 所有任务完成后归档到archive目录
4. **代码简洁不超过500行** - 单文件不超过500行，超过即模块化重构

通过规范的文档管理和AI协作指令，让AI辅助开发更高效！

