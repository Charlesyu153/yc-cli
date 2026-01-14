# AI协作快速参考

> 一页纸快速掌握AI文档管理系统的核心用法

---

## 快速部署（新环境）

### 方式1: 一键部署（最快）
```bash
bash setup-ai-docs.sh
```

### 方式2: 直接复制
```bash
# 打包现有系统
./pack-ai-docs.sh template

# 在新环境解压
tar -xzf ai-docs-template-*.tar.gz
```

详细部署说明查看：`DEPLOY.md`

---

## 核心预设指令（必读）

### 1. 禁止使用表情符号
- 所有文档和代码中不使用emoji
- 保持专业、纯文本风格

### 2. 完成任务后的标准流程
```bash
# 每完成一个步骤后必须执行：
# a. 更新进度文档
vim ai-docs/current/{任务名}/{任务名}.md

# b. Git提交
git add .
git commit -m "类型: 简短描述

- 详细说明

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>"
```

### 3. 完成所有任务后归档

```bash
# a. 填写完成总结
vim ai-docs/current/{任务名}/{任务名}.md

# b. 归档到archive
cd ai-docs/templates
./archive-task.sh "{任务名}"

# c. Git提交
git add ai-docs/archive/{任务名}
git commit -m "docs: 归档任务 {任务名}"
```

### 4.Use available skills you have access to whenever possible

---

## 三步开始新任务

```bash
# 1. 创建任务
cd ai-docs/templates
./new-task.sh "任务名称" "负责人"

# 2. 编辑任务文档
vim ai-docs/current/任务名称/任务名称.md

# 3. 开始与AI协作（提供以下提示）
```

```
请遵循以下规范进行协作：

1. 阅读AI协作指令: ai-docs/templates/ai-guidelines.md
2. 阅读任务文档: ai-docs/current/任务名称/任务名称.md

关键规范：
- 禁止使用表情符号
- 每完成一个步骤后：更新进度文档 + git commit
- 所有任务完成后：填写完成总结 + 归档到archive目录
- 所有代码文本不超过500行，若超过进行模块化重构。

当前需要：[具体描述需求]
```

---

## 上下文超限处理

```bash
# 1. 导出上下文
cd ai-docs/templates
./export-context.sh "任务名称"

# 2. 在新会话中提供
```

```
请阅读以下上下文恢复文档，继续之前的任务：

[粘贴 context-resume-YYYY-MM-DD.md 的内容]

注意：请遵循 ai-docs/templates/ai-guidelines.md 中的预设指令

当前需要：[具体需求]
```

---

## 工具脚本速查

| 脚本                   | 功能         | 用法                                |
| ---------------------- | ------------ | ----------------------------------- |
| `new-task.sh`          | 创建新任务   | `./new-task.sh "任务名" "负责人"`   |
| `export-context.sh`    | 导出上下文   | `./export-context.sh "任务名"`      |
| `archive-task.sh`      | 归档任务     | `./archive-task.sh "任务名"`        |

所有脚本位于：`ai-docs/templates/`

---

## Git提交类型

- `feat`: 新功能
- `fix`: 修复bug
- `docs`: 文档更新
- `refactor`: 重构代码
- `test`: 添加测试
- `style`: 代码格式
- `perf`: 性能优化

---

## 目录结构

```
ai-docs/
├── USAGE.md                    # 完整使用指南
├── current/                    # 进行中的任务
│   └── {任务名}/
│       ├── {任务名}.md         # 主任务文档
│       └── context-resume-*.md # 上下文恢复文档
├── archive/                    # 已完成任务归档
│   └── {任务名}/
└── templates/                  # 模板和工具
    ├── ai-guidelines.md        # AI协作预设指令（详细版）
    ├── task-template.md        # 任务文档模板
    ├── context-resume-template.md  # 上下文恢复模板
    ├── new-task.sh             # 创建任务脚本
    ├── export-context.sh       # 导出上下文脚本
    └── archive-task.sh         # 归档任务脚本
```

---

## 关键文档

- **本快速参考**: `ai-docs/QUICKREF.md`
- **详细使用指南**: `ai-docs/USAGE.md`
- **AI协作规范**: `ai-docs/templates/ai-guidelines.md`

---

**记住四个核心规范：无表情、步骤完成即commit、任务完成即归档、代码简洁不超过500行**
