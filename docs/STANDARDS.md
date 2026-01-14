# YC-CLI 项目规范

本文档定义 YC-CLI 项目的统一规范，所有子系统（agents、ai-docs）都遵循这些规范。

配置文件：`config/standards.json`

---

## 核心规范（四条）

### 1. 禁止使用表情符号

**规则**：
- 所有文档和代码中禁止使用 emoji 表情
- 使用纯文本表达，保持专业性
- 例外：用户明确要求使用表情时

**原因**：保持文档的专业性和跨平台兼容性

### 2. 步骤完成即 commit

**规则**：
- 每完成一个步骤后，必须更新文档并 git commit
- 不要批量提交多个步骤的修改

**工作流程**：
```
1. 执行任务
2. 更新进度文档
3. git add .
4. git commit -m "描述性的提交信息"
5. 继续下一个任务
```

**原因**：保持提交历史清晰，便于回滚和追踪

### 3. 任务完成即归档

**规则**：
- 所有任务完成后，必须归档到 `ai-docs/archive/` 目录
- 使用 `archive-task.sh` 脚本进行归档

**归档步骤**：
```bash
cd ai-docs/templates
./archive-task.sh "任务名称"
```

**原因**：保持工作区整洁，便于知识沉淀和回顾

### 4. 代码文件不超过 500 行

**规则**：
- 单个代码文件或脚本不超过 500 行
- 超过 500 行时必须进行模块化重构

**重构原则**：
- 按功能拆分成多个模块
- 使用清晰的函数/类结构
- 保持代码可读性和可维护性
- 在任务文档中记录重构决策

**原因**：提高代码可维护性，降低复杂度

---

## 编码规范

### Python 规范

- **风格**：PEP 8
- **缩进**：4 spaces
- **命名**：snake_case
- **最大行数**：500 行
- **导入**：按标准库、第三方库、本地模块分组

**示例**：
```python
import os
import sys

import numpy as np
import pandas as pd

from my_module import my_function


def process_data(input_file):
    """处理数据文件"""
    pass
```

### R 规范

- **缩进**：2 spaces
- **命名**：snake_case
- **路径处理**：使用 `file.path()`
- **最大行数**：500 行
- **日志**：适当的日志输出

**示例**：
```r
process_data <- function(input_file) {
  # 处理数据文件
  output_path <- file.path("output", "result.csv")
  return(output_path)
}
```

---

## Git 提交规范

### 提交信息格式

```
<type>: <subject>

<body>

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>
```

### Type 类型

| Type | 说明 |
|------|------|
| feat | 新功能 |
| fix | 修复 bug |
| docs | 文档更新 |
| refactor | 重构代码 |
| test | 添加测试 |
| style | 代码格式调整 |
| perf | 性能优化 |

### 示例

```bash
git commit -m "feat: 添加 CAF 亚型统计分析功能

- 实现按样本分组的统计检验
- 添加可视化函数 plot_subtype_stats()
- 更新主分析脚本调用新功能

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>"
```

---

## Skills 使用原则

### 核心原则

**优先使用可用的 Skills 来完成任务**

- 在执行任务前，先检查是否有相关的 Skill 可用
- 充分利用 Skills 的专业能力，提升开发效率
- Skills 可以与 Agents 配合使用

### 可用的 Skills 库

1. **muratcankoylan/Agent-Skills-for-Context-Engineering**
   - 类型：通用开发技能
   - 能力：代码分析、重构、测试、文档生成

2. **claude-scientific-skills**
   - 类型：科学计算和数据分析技能
   - 能力：生物信息学、机器学习、数据可视化

详细使用指南请查看：`docs/SKILLS-GUIDE.md`

---

## 文档规范

### 格式

- **格式**：Markdown
- **风格**：简洁、专业、无表情符号
- **代码引用**：使用 `path/to/file.ext:line_number` 格式

### 代码引用示例

- 文件引用：`path/to/file.R`
- 精确定位：`path/to/file.R:42`（第 42 行）
- 函数引用：`functionName()`
- 类引用：`ClassName`

---

## 检查清单

### 开始任务时
- [ ] 已阅读任务文档
- [ ] 理解当前进度和下一步
- [ ] 理解关键决策

### 编写代码时
- [ ] 检查文件行数是否超过 500 行
- [ ] 如超过 500 行，进行模块化重构
- [ ] 在任务文档中记录重构决策
- [ ] 优先使用可用的 Skills

### 完成每个步骤后
- [ ] 更新进度文档
- [ ] 更新变更文件清单
- [ ] 执行 git commit
- [ ] 未使用表情符号

### 完成所有任务后
- [ ] 填写完成总结
- [ ] 执行最终 git commit
- [ ] 归档到 archive 目录
- [ ] 验证归档成功

---

**配置文件**：`config/standards.json`
**版本**：v1.0.0
**更新日期**：2026-01-14
