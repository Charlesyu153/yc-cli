# AI Agents 配置

PCF 项目的 AI agents 配置，参考 feature-dev 插件功能模式。

## 快速开始

### 三个 Agent 角色

1. **code-architect** - 代码架构师
   - 架构设计、技术选型、实施规划

2. **code-explorer** - 代码探索者
   - 代码搜索、结构分析、文档生成

3. **code-reviewer** - 代码审查者
   - 质量检查、问题发现、改进建议

### 基本使用

在与 Claude Code 对话时，使用以下格式：

```
请以 [agent-name] 的角色，[具体任务描述]
```

示例：
```
请以 code-explorer 的角色，查找所有处理 CAF 细胞类型的代码
```

## 配置文件

```
.claude/agents/
├── README.md                  # 本文档
├── agents.json                # 主配置文件
├── code-architect.json        # 架构师配置
├── code-explorer.json         # 探索者配置
├── code-reviewer.json         # 审查者配置
└── USAGE_EXAMPLES.md          # 详细使用示例
```

## Feature-Dev 工作流

完整的功能开发流程：

1. **探索** (code-explorer) - 理解现有代码
2. **设计** (code-architect) - 设计实现方案
3. **实现** - 编写代码
4. **审查** (code-reviewer) - 质量检查
5. **提交** - 更新文档并提交

## 文档导航

- **详细使用示例**: `USAGE_EXAMPLES.md`
- **项目规范**: `../../docs/AGENTS.md`
- **AI 协作规范**: `../../ai-docs/templates/ai-guidelines.md`

## 项目规范

所有 agents 遵循的规范：
- 禁止使用表情符号
- 代码文件不超过 500 行
- Python: PEP 8, 4-space 缩进, snake_case
- R: 2-space 缩进, snake_case, file.path()

---

**版本**: v1.0.0
**创建日期**: 2026-01-12
