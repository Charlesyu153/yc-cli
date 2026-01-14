# Plan模式 + AI-DOCS协作工作流程指南

## 完整工作流程

### 阶段1：接收任务
```
1. 理解用户需求
2. 判断是否需要进入Plan模式
3. 如需Plan模式，使用EnterPlanMode工具
```

### 阶段2：Plan模式（如适用）
```
1. 探索代码库（Glob/Grep/Read）
2. 设计实现方案
3. 创建ai-docs任务文档（在plan模式下）
   - 位置：ai-docs/current/任务名称/任务名称.md
   - 遵循task-template.md格式
   - 不使用emoji
4. 将计划写入plan文件
5. 使用ExitPlanMode请求批准
```

### 阶段3：实现阶段
```
1. 创建/更新TodoWrite列表
2. 按步骤实现功能
3. 每完成一个步骤：
   a. 更新TodoWrite状态
   b. 更新ai-docs任务文档
   c. Git commit（遵循commit规范）
4. 保持代码文件不超过500行
```

### 阶段4：测试阶段
```
1. 运行测试
2. 记录测试结果到ai-docs
3. 修复问题并commit
```

### 阶段5：完成阶段
```
1. 填写ai-docs完成总结
2. 归档任务文档到ai-docs/archive/
3. 清理TodoWrite列表
```

## Plan模式与AI-DOCS配合

### 在Plan模式下创建任务文档
```markdown
# 正确做法
1. 进入Plan模式
2. 探索代码库
3. 创建ai-docs/current/任务名称/任务名称.md
4. 在plan文件中引用任务文档位置
5. ExitPlanMode请求批准

# 错误做法
1. 进入Plan模式
2. 直接写plan文件，不创建ai-docs任务文档
3. 或者在实现阶段才创建ai-docs文档
```

### 任务文档与Plan文件的关系
- **Plan文件**：临时的实现计划，用于获得用户批准
- **AI-DOCS任务文档**：持久的任务追踪，记录整个生命周期
- **关系**：Plan文件内容应该同步到AI-DOCS任务文档的"技术方案"部分

## 实际案例：TME分析任务

### 步骤1：Plan模式
```
1. 用户请求：设计PCF项目TME分析
2. 判断：需要Plan模式（新功能、多文件、架构决策）
3. 进入Plan模式
4. 探索：读取recommend.md、现有脚本、工具函数
5. 设计：CLDN18.2验证 + 全细胞30um邻域分析
6. 创建：ai-docs/current/TME分析/TME分析.md（不使用emoji）
7. 写入：plan.md详细实现计划
8. 退出：ExitPlanMode请求批准
```

### 步骤2：实现阶段
```
1. 创建TodoWrite：
   - 创建ai-docs任务文档结构
   - 实现CLDN18.2分类验证分析脚本
   - 实现全细胞30um邻域分析脚本

2. 实现第一个脚本：
   - 创建script/cldn18_validation/目录
   - 实现cldn18_niche_comparison.R（275行，<500行）
   - 更新TodoWrite状态为completed
   - 更新ai-docs进度追踪
   - Git commit: "feat: 实现CLDN18.2分类验证分析脚本"

3. 实现第二个脚本：
   - 创建script/all_cell_niche/目录
   - 实现all_cell_niche_step1.R（261行，<500行）
   - 更新TodoWrite状态为completed
   - 更新ai-docs进度追踪
   - Git commit: "feat: 实现全细胞30um邻域分析脚本"
```

### 步骤3：测试阶段
```
1. 测试20P样本
2. 发现热图绘制错误
3. 修复零方差列问题
4. Git commit: "fix: 修复热图绘制零方差列问题"
5. 重新测试成功
6. 记录测试结果到ai-docs
```

### 步骤4：完成阶段
```
1. 更新ai-docs任务文档：
   - 标记所有步骤为已完成
   - 填写测试结果部分
   - 记录问题与解决方案
2. 准备归档（如任务完全结束）
```

## 快速检查清单

### Plan模式检查
- [ ] 任务需要Plan模式吗？
- [ ] 已充分探索代码库？
- [ ] 已创建ai-docs任务文档？
- [ ] Plan文件包含足够技术细节？
- [ ] 已使用ExitPlanMode请求批准？

### AI-DOCS检查
- [ ] 文档位置正确（在子目录中）？
- [ ] 没有使用emoji？
- [ ] 遵循task-template.md格式？
- [ ] 进度追踪及时更新？
- [ ] 每完成一步就commit？
- [ ] 代码文件不超过500行？
- [ ] Commit消息格式正确？

### 实现阶段检查
- [ ] TodoWrite列表已创建？
- [ ] 每完成一步更新状态？
- [ ] 每完成一步更新ai-docs？
- [ ] 每完成一步git commit？
- [ ] 代码符合项目规范？

## 常见问题

### Q1: 什么时候创建ai-docs任务文档？
A: 在Plan模式下创建，作为计划的一部分。不要等到实现阶段。

### Q2: Plan文件和ai-docs任务文档有什么区别？
A: Plan文件是临时的实现计划，ai-docs是持久的任务追踪记录。

### Q3: 如果任务很简单，还需要Plan模式和ai-docs吗？
A: 简单任务（单行修复、明确需求）不需要Plan模式，但如果涉及多个步骤，仍建议使用ai-docs追踪。

### Q4: 可以在一个commit中提交多个功能吗？
A: 不可以。每完成一个功能模块就应该commit一次。

### Q5: 代码文件超过500行怎么办？
A: 立即重构，按功能拆分成多个文件。

### Q6: ai-docs目录被.gitignore忽略，为什么还要维护？
A: ai-docs用于AI协作追踪，虽然不提交到git，但对任务管理和进度追踪很重要。
