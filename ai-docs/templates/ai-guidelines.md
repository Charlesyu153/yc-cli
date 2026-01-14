# AI协作预设指令

> 在与AI开始协作时，将此文档内容提供给AI，确保遵循项目规范
>
> 本文档遵循项目统一规范，详见: ../../config/standards.json

---

## 📋 基本规范

### 1. 沟通风格

**禁止使用表情符号**
- 所有文档和代码中禁止使用emoji表情
- 使用纯文本表达，保持专业性
- 例外：用户明确要求使用表情时

### 2. 文档管理

**任务文档位置**
- 当前任务文档：`ai-docs/current/{任务名}/{任务名}.md`
- 始终优先阅读和更新该文档

**禁止使用 README.md**
- 任务文档必须使用任务名称命名
- 不要创建或使用 README.md

### 3. 代码质量规范

**代码行数限制**
- 单个代码文件或脚本不超过500行
- 超过500行时必须进行模块化重构
- 重构原则：
  - 按功能拆分成多个模块
  - 使用清晰的函数/类结构
  - 保持代码可读性和可维护性
  - 在任务文档中记录重构决策

---

## 🔄 工作流程规范

### 阶段1: 任务启动

1. **阅读任务文档**
   ```
   请阅读任务文档: ai-docs/current/{任务名}/{任务名}.md
   ```

2. **理解当前状态**
   - 查看"进度追踪"了解已完成和待完成的工作
   - 查看"下一步"明确当前要做的事
   - 查看"关键决策"理解重要的技术选择

### 阶段2: 执行任务

**每完成一个步骤后，必须执行以下操作：**

1. **更新进度文档**
   - 标记已完成的步骤为 `[x]`
   - 更新"进行中"部分
   - 更新"下一步"说明
   - 在"变更文件清单"中记录修改的文件

2. **提交代码到本地仓库**
   ```bash
   git add .
   git commit -m "描述性的提交信息
   
   - 完成的具体工作
   - 修改的关键文件
   
   Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>"
   ```

**完整的步骤模板：**
```
1. 执行任务A（如：实现某个功能）
2. 更新进度文档 ai-docs/current/{任务名}/{任务名}.md
3. git add .
4. git commit -m "feat: 完成任务A"
5. 继续下一个任务
```

### 阶段3: 任务完成

**当所有任务都完成后，执行以下操作：**

1. **填写完成总结**
   - 在任务文档的"完成总结"部分填写：
     - 实现概述
     - 经验教训

2. **最终提交**
   ```bash
   git add ai-docs/current/{任务名}/{任务名}.md
   git commit -m "docs: 完成{任务名}任务总结"
   ```

3. **归档任务**
   ```bash
   cd ai-docs/templates
   ./archive-task.sh "{任务名}"
   ```

   或手动归档：
   ```bash
   mv ai-docs/current/{任务名} ai-docs/archive/
   ```

4. **验证归档**
   ```bash
   # 确认任务已移动到archive
   ls ai-docs/archive/{任务名}/
   ```

---

## 📝 文档更新规范

### 进度追踪格式

```markdown
### 已完成
- [x] 步骤1: 描述
- [x] 步骤2: 描述

### 进行中
- [ ] 步骤3: 描述 (当前进度: 正在实现XX功能)

### 待完成
- [ ] 步骤4: 描述
- [ ] 步骤5: 描述

### 下一步
> 接下来需要：实现步骤3的XX功能，涉及文件 path/to/file.R:line
```

### 变更文件清单格式

```markdown
### 新增文件

| **文件**                          | **说明**     |
| --------------------------------- | ------------ |
| `script/analysis/new_function.R`  | 新增分析函数 |

### 修改文件

| **文件**                             | **修改内容**       |
| ------------------------------------ | ------------------ |
| `script/main/process.R:42`           | 优化数据处理逻辑   |
| `script/visualization/plot.R:156`    | 修改绘图参数       |
```

### 代码引用格式

- 文件引用：`path/to/file.R`
- 精确定位：`path/to/file.R:42` (第42行)
- 函数引用：`functionName()`
- 类引用：`ClassName`

---

## 💡 Git提交规范

### 提交信息格式

```
<type>: <subject>

<body>

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>
```

### Type类型

- `feat`: 新功能
- `fix`: 修复bug
- `docs`: 文档更新
- `refactor`: 重构代码
- `test`: 添加测试
- `style`: 代码格式调整
- `perf`: 性能优化

### 示例

```bash
git commit -m "feat: 添加CAF亚型统计分析功能

- 实现按样本分组的统计检验
- 添加可视化函数plot_subtype_stats()
- 更新主分析脚本调用新功能

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>"
```

---

## ⚠️ 重要提醒

### 必须执行的操作

1. **完成每个步骤后**

   - 更新进度文档
   - 执行 git commit

2. **完成所有任务后**
   - 填写完成总结
   - 归档到 archive 目录

3. **编写代码时**
   - 检查文件行数，单文件不超过500行
   - 超过500行立即进行模块化重构
   - 在任务文档中记录重构决策

4. ```txt
   Use available skills you have access to whenever possible.
   ```

### 4. Skills 使用规范

**优先使用 Skills**
- 在执行任务前，检查是否有相关的 Skill 可用
- 充分利用 Skills 的专业能力，提升开发效率
- Skills 可以与 Agents 配合使用

**可用的 Skills 库**
- muratcankoylan/Agent-Skills-for-Context-Engineering（通用开发技能）
- claude-scientific-skills（科学计算和数据分析技能）

**Skills 使用场景**
- 代码分析和重构
- 测试和质量检查
- 文档生成
- 数据分析和可视化
- 生物信息学处理

### 禁止的操作

1. 不要使用表情符号
2. 不要创建 README.md 作为任务文档
3. 不要忘记更新进度文档
4. 不要忘记 git commit
5. 不要在任务完成后忘记归档
6. 不要编写超过500行的单文件代码

---

## 🔄 上下文超限处理

当会话上下文即将超限时：

1. **导出上下文**
   ```bash
   cd ai-docs/templates
   ./export-context.sh "{任务名}"
   ```

2. **提供给新会话**
   ```
   请阅读以下上下文恢复文档，继续之前的任务...
   [粘贴 context-resume-{日期}.md 的内容]
   
   注意：请遵循 ai-docs/templates/ai-guidelines.md 中的预设指令
   ```

---

## 📋 检查清单

使用以下清单确保遵循规范：

**开始任务时**
- [ ] 已阅读任务文档
- [ ] 理解当前进度和下一步
- [ ] 理解关键决策

**编写代码时**
- [ ] 检查文件行数是否超过500行
- [ ] 如超过500行，进行模块化重构
- [ ] 在任务文档中记录重构决策
- [ ] 优先使用可用的 Skills

**完成每个步骤后**
- [ ] 更新进度文档
- [ ] 更新变更文件清单
- [ ] 执行 git commit
- [ ] 未使用表情符号

**完成所有任务后**
- [ ] 填写完成总结
- [ ] 执行最终 git commit
- [ ] 归档到 archive 目录
- [ ] 验证归档成功

---

## 🎯 快速参考

### 常用命令

```bash
# 查看当前任务
ls ai-docs/current/

# 编辑任务文档
vim ai-docs/current/{任务名}/{任务名}.md

# 导出上下文
cd ai-docs/templates && ./export-context.sh "{任务名}"

# 归档任务
cd ai-docs/templates && ./archive-task.sh "{任务名}"

# Git操作
git status
git add .
git commit -m "feat: 描述"
git log --oneline -5
```

### 会话开始提示语模板

```
请遵循以下规范进行协作：

1. 阅读任务文档: ai-docs/current/{任务名}/{任务名}.md
2. 遵循 ai-docs/templates/ai-guidelines.md 中的预设指令
3. 特别注意：
   - 禁止使用表情符号
   - 每完成一个步骤后更新文档并 git commit
   - 所有任务完成后归档到 archive 目录
   - 所有代码文件不超过500行，超过即模块化重构
   - 优先使用可用的 Skills

当前需要：[具体描述需求]
```

---

**最后提醒：这些规范是为了保持项目的连续性和可追溯性，请严格遵守！**

