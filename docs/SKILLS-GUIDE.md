# Skills 使用指南

本文档详细说明如何在 YC-CLI 中使用 Skills 来提升开发效率。

---

## 概述

YC-CLI 集成了两个强大的 Skills 库，提供丰富的专业能力。

**核心原则**：优先使用可用的 Skills 来完成任务

---

## 可用的 Skills 库

### 1. muratcankoylan/Agent-Skills-for-Context-Engineering

**类型**：通用开发技能集

**主要能力**：
- 代码分析和质量检查
- 代码重构和优化
- 测试生成和执行
- 文档生成和维护
- 项目结构分析

**适用场景**：
- 各类软件开发任务
- 代码质量提升
- 自动化测试
- 文档维护

**安装方式**：

```bash
# 方式 1: 使用 Claude Code CLI 安装
claude mcp install github --owner muratcankoylan --repo Agent-Skills-for-Context-Engineering

# 方式 2: 手动配置
# 编辑 ~/.config/claude/claude_desktop_config.json
# 添加以下配置：
{
  "mcpServers": {
    "agent-skills": {
      "command": "npx",
      "args": ["-y", "@muratcankoylan/agent-skills-mcp"]
    }
  }
}
```

**仓库地址**：https://github.com/muratcankoylan/Agent-Skills-for-Context-Engineering

### 2. claude-scientific-skills

**类型**：科学计算和数据分析技能集

**主要能力**：
- 生物信息学数据处理
- 机器学习模型训练
- 数据可视化
- 统计分析
- 科学计算

**适用场景**：
- 数据分析项目
- 科学研究
- 机器学习应用
- 生物信息学处理

**安装方式**：

```bash
# 方式 1: 使用 Claude Code CLI 安装
claude mcp install github --owner K-Dense-AI --repo claude-scientific-skills

# 方式 2: 手动配置
# 编辑 ~/.config/claude/claude_desktop_config.json
# 添加以下配置：
{
  "mcpServers": {
    "scientific-skills": {
      "command": "npx",
      "args": ["-y", "@k-dense-ai/scientific-skills-mcp"]
    }
  }
}
```

**仓库地址**：https://github.com/K-Dense-AI/claude-scientific-skills

### 验证安装

安装完成后，重启 Claude Code 并验证：

```
请列出当前可用的 Skills，以及它们的主要功能。
```

如果安装成功，Claude 会列出已安装的 Skills 及其功能。

---

## Skills 使用原则

### 核心原则

1. **优先使用**：在执行任务前，先检查是否有相关的 Skill 可用
2. **充分利用**：充分利用 Skills 的专业能力，提升开发效率
3. **协同工作**：Skills 可以与 Agents 配合使用，形成更强大的工作流

### 使用时机

- 代码分析和重构
- 测试和质量检查
- 文档生成
- 数据分析和可视化
- 生物信息学处理
- 机器学习任务

---

## Skills 使用示例

### 示例 1: 使用 Skills 进行代码分析

**场景**：分析 Python 代码的复杂度和质量

**提示语**：
```
请使用可用的代码分析 Skill，分析 src/data_processor.py 的代码质量。

重点关注：
1. 代码复杂度
2. 潜在的性能问题
3. 可维护性评估
4. 改进建议
```

**预期输出**：
- 代码复杂度报告
- 性能瓶颈分析
- 可维护性评分
- 具体改进建议

### 示例 2: 使用 Skills 进行代码重构

**场景**：重构超过 500 行的代码文件

**提示语**：
```
请使用可用的代码重构 Skill，重构 src/legacy_module.py

要求：
1. 分析当前代码结构
2. 识别可以模块化的部分
3. 提供重构方案
4. 生成重构后的代码
```

**预期输出**：
- 代码结构分析
- 模块化方案
- 重构后的代码
- 重构说明文档

### 示例 3: 使用 Skills 进行数据分析

**场景**：使用科学计算 Skills 处理生物信息学数据

**提示语**：
```
请使用 claude-scientific-skills 中的生物信息学工具，分析单细胞数据。

需求：
1. 读取 h5ad 格式数据
2. 进行质控和标准化
3. 生成可视化报告
4. 输出处理后的数据
```

**预期输出**：
- 数据质控报告
- 标准化后的数据
- 可视化图表
- 分析总结

### 示例 4: 使用 Skills 生成测试

**场景**：为现有代码生成单元测试

**提示语**：
```
请使用可用的测试生成 Skill，为 src/utils.py 生成单元测试。

要求：
1. 分析函数功能
2. 生成测试用例
3. 包含边界情况测试
4. 使用 pytest 框架
```

**预期输出**：
- 完整的测试文件
- 测试用例说明
- 测试覆盖率报告

### 示例 5: 使用 Skills 生成文档

**场景**：为代码生成 API 文档

**提示语**：
```
请使用可用的文档生成 Skill，为 src/api/ 目录生成 API 文档。

要求：
1. 提取所有公共函数和类
2. 生成函数签名和参数说明
3. 包含使用示例
4. 输出 Markdown 格式
```

**预期输出**：
- 完整的 API 文档
- 函数签名和参数说明
- 使用示例
- 格式化的 Markdown 文件

---

## Skills 与 Agents 协作

### 协作模式 1: Explore + Skills

**场景**：使用 code-explorer 和 Skills 分析代码库

**提示语**：
```
请以 code-explorer 的角色，使用可用的 Skills，分析项目的代码结构。

要求：
1. 使用 Skills 分析代码复杂度
2. 识别代码组织模式
3. 生成项目结构文档
4. 提供改进建议
```

### 协作模式 2: Architect + Skills

**场景**：使用 code-architect 和 Skills 设计方案

**提示语**：
```
请以 code-architect 的角色，使用可用的 Skills，设计一个数据处理流程。

要求：
1. 使用 Skills 分析现有代码模式
2. 设计新的架构方案
3. 提供实施步骤
4. 生成架构文档
```

### 协作模式 3: Reviewer + Skills

**场景**：使用 code-reviewer 和 Skills 审查代码

**提示语**：
```
请以 code-reviewer 的角色，使用可用的 Skills，审查 src/core/ 目录。

要求：
1. 使用 Skills 进行代码质量分析
2. 识别潜在问题
3. 提供改进建议
4. 生成审查报告
```

---

## Skills 最佳实践

### 1. 任务前检查

在开始任务前，先问自己：
- 是否有相关的 Skill 可以使用？
- 使用 Skill 能否提升效率？
- 如何将 Skill 与 Agent 结合使用？

### 2. 明确需求

使用 Skills 时，明确说明：
- 要完成什么任务
- 期望得到什么输出
- 有什么特殊要求

### 3. 结合 Agents

- 使用 code-explorer 探索代码，然后用 Skills 分析
- 使用 code-architect 设计方案，然后用 Skills 实现
- 使用 code-reviewer 审查代码，然后用 Skills 优化

### 4. 验证结果

使用 Skills 后，验证：
- 输出是否符合预期
- 是否需要进一步调整
- 是否需要人工审查

---

## 常见使用场景

### 场景 1: 代码质量提升

```
1. 使用代码分析 Skill 识别问题
2. 使用代码重构 Skill 优化代码
3. 使用测试生成 Skill 增加测试覆盖率
4. 使用 code-reviewer 验证改进效果
```

### 场景 2: 数据分析项目

```
1. 使用 claude-scientific-skills 处理数据
2. 使用可视化 Skill 生成图表
3. 使用文档生成 Skill 创建报告
4. 使用 code-architect 设计分析流程
```

### 场景 3: 新功能开发

```
1. 使用 code-explorer 和 Skills 分析现有代码
2. 使用 code-architect 和 Skills 设计方案
3. 使用 Skills 实现功能
4. 使用测试生成 Skill 创建测试
5. 使用 code-reviewer 和 Skills 审查代码
```

### 场景 4: 文档维护

```
1. 使用文档生成 Skill 创建 API 文档
2. 使用 Skills 生成使用示例
3. 使用 Skills 更新 README
4. 使用 code-explorer 验证文档完整性
```

---

## 技巧和窍门

### 技巧 1: 组合使用

不要只使用单个 Skill，尝试组合使用：
```
请使用代码分析 Skill 分析代码，然后使用重构 Skill 优化，最后使用测试生成 Skill 创建测试。
```

### 技巧 2: 迭代改进

使用 Skills 进行迭代改进：
```
1. 第一轮：使用 Skill 生成初始版本
2. 第二轮：使用 Skill 分析并改进
3. 第三轮：使用 Skill 优化和完善
```

### 技巧 3: 明确约束

明确告诉 Skills 你的约束条件：
```
请使用代码重构 Skill 优化代码，要求：
- 保持向后兼容
- 不改变公共 API
- 遵循项目编码规范
```

### 技巧 4: 请求解释

要求 Skills 解释其操作：
```
请使用代码分析 Skill 分析代码，并解释：
- 为什么识别出这些问题
- 建议的改进方案的原理
- 预期的改进效果
```

---

## 常见问题

### Q1: 如何知道有哪些 Skills 可用？

**答**：在任务开始前，询问 Claude：
```
请列出当前可用的 Skills，以及它们的主要功能。
```

### Q2: Skills 和 Agents 有什么区别？

**答**：
- **Agents**：定义 AI 的角色和工作方式（架构师、探索者、审查者）
- **Skills**：提供具体的专业能力（代码分析、数据处理等）
- **关系**：Agents 可以使用 Skills 来完成任务

### Q3: 什么时候应该使用 Skills？

**答**：
- 需要专业能力时（代码分析、数据处理等）
- 需要自动化重复任务时
- 需要提升效率时
- 需要标准化输出时

### Q4: Skills 的输出可靠吗？

**答**：
- Skills 的输出通常是可靠的
- 但建议结合 code-reviewer 进行验证
- 对于关键任务，建议人工审查

---

## 快速参考

### Skills 使用模板

```
请使用可用的 [Skill 类型] Skill，[具体任务]。

要求：
1. [要求 1]
2. [要求 2]
3. [要求 3]

期望输出：
- [输出 1]
- [输出 2]
```

### Skills 与 Agents 协作模板

```
请以 [agent-name] 的角色，使用可用的 Skills，[具体任务]。

要求：
1. 使用 Skills 进行 [具体操作]
2. [其他要求]

期望输出：
- [输出内容]
```

---

**相关文档**：
- Agents 使用指南：`docs/AGENTS-GUIDE.md`
- 项目规范：`docs/STANDARDS.md`
- 快速开始：`docs/GETTING-STARTED.md`

**配置文件**：`config/standards.json`

**版本**：v1.0.0
**更新日期**：2026-01-14
