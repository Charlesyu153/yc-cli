# AI文档管理系统

一个专为AI辅助开发设计的任务文档管理系统，帮助在长期项目中保持上下文连续性。

---

## 快速开始

### 在当前项目使用

```bash
# 1. 创建新任务
cd ai-docs/templates
./new-task.sh "任务名称" "负责人"

# 2. 查看快速参考
cat ../QUICKREF.md
```

### 部署到新环境

```bash
# 方式1: 一键部署
bash ai-docs/setup-ai-docs.sh

# 方式2: 打包后部署
cd ai-docs
./pack-ai-docs.sh template
# 将生成的tar.gz传输到新环境后解压
```

详细部署说明：`DEPLOY.md`

---

## 文档导航

| 文档 | 说明 |
|------|------|
| `QUICKREF.md` | 快速参考指南（一页纸） |
| `USAGE.md` | 完整使用指南（详细） |
| `DEPLOY.md` | 部署指南（新环境部署） |
| `templates/ai-guidelines.md` | AI协作预设指令 |

---

## 目录结构

```
ai-docs/
├── README.md              # 本文档
├── QUICKREF.md            # 快速参考
├── USAGE.md               # 完整使用指南
├── DEPLOY.md              # 部署指南
├── setup-ai-docs.sh       # 一键部署脚本
├── pack-ai-docs.sh        # 打包工具
├── current/               # 进行中的任务
├── archive/               # 已完成任务归档
└── templates/
    ├── task-template.md            # 任务文档模板
    ├── context-resume-template.md  # 上下文恢复模板
    ├── ai-guidelines.md            # AI协作预设指令
    ├── new-task.sh                 # 创建任务脚本
    ├── export-context.sh           # 导出上下文脚本
    └── archive-task.sh             # 归档任务脚本
```

---

## 四个核心规范

1. **无表情** - 禁止使用表情符号
2. **步骤完成即commit** - 每完成一个步骤后更新文档并git commit
3. **任务完成即归档** - 所有任务完成后归档到archive目录
4. **代码简洁不超过500行** - 单文件不超过500行，超过即模块化重构

---

## 工具脚本

| 脚本 | 功能 | 使用 |
|------|------|------|
| `new-task.sh` | 创建新任务 | `./new-task.sh "任务名" "负责人"` |
| `export-context.sh` | 导出上下文 | `./export-context.sh "任务名"` |
| `archive-task.sh` | 归档任务 | `./archive-task.sh "任务名"` |
| `setup-ai-docs.sh` | 部署系统 | `bash setup-ai-docs.sh [目标路径]` |
| `pack-ai-docs.sh` | 打包系统 | `./pack-ai-docs.sh [template|full|deploy]` |

---

## 与AI协作的提示语模板

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

---

## 核心价值

1. **持续上下文** - 长期任务不丢失进度
2. **快速恢复** - 新会话快速进入状态
3. **知识沉淀** - 记录决策和经验
4. **团队协作** - 便于交接和协作
5. **规范流程** - AI协作遵循统一规范

---

## 部署方案

### 个人新项目
```bash
bash ai-docs/setup-ai-docs.sh ./my-project/ai-docs
```

### 团队推广
```bash
# 打包模板
cd ai-docs
./pack-ai-docs.sh template

# 分发给团队成员
# 团队成员解压后使用
```

### 现有项目迁移
直接复制整个 `ai-docs/` 目录到新项目即可。

详细部署方案查看：`DEPLOY.md`

---

**通过规范的文档管理和AI协作指令，让AI辅助开发更高效！**
