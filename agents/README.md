# AI Agents 配置

YC-CLI 的 AI Agents 配置系统，提供三个专用 AI 角色。

## 三个 Agent 角色

1. **code-architect** - 代码架构师
   - 架构设计、技术选型、实施规划

2. **code-explorer** - 代码探索者
   - 代码搜索、结构分析、文档生成

3. **code-reviewer** - 代码审查者
   - 质量检查、问题发现、改进建议

## 基本使用

在与 Claude Code 对话时，使用以下格式：

```
请以 [agent-name] 的角色，[具体任务描述]
```

示例：
```
请以 code-explorer 的角色，查找所有处理用户认证的代码
```

## Skills 集成

优先使用可用的 Skills 来完成任务：
- muratcankoylan/Agent-Skills-for-Context-Engineering（通用开发技能）
- claude-scientific-skills（科学计算和数据分析技能）

## 详细文档

完整的使用指南和示例请查看：

- **详细使用指南**: `../docs/AGENTS-GUIDE.md`
- **Skills 使用指南**: `../docs/SKILLS-GUIDE.md`
- **项目规范**: `../docs/STANDARDS.md`

---

**版本**: v1.0.0
**更新日期**: 2026-01-14
