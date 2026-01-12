# AI文档管理系统 - 部署指南

> 在新环境中快速部署AI文档管理系统的完整指南

---

## 部署方案总览

提供三种部署方案，根据实际情况选择：

| 方案 | 适用场景 | 优点 | 缺点 |
|------|----------|------|------|
| 方案1: 一键脚本部署 | 新项目快速搭建 | 最快速，自动化 | 需要bash环境 |
| 方案2: 直接复制目录 | 已有完整配置的ai-docs | 保留所有自定义内容 | 需要手动调整 |
| 方案3: 手动逐步创建 | 需要深度定制 | 完全可控 | 最耗时 |

---

## 方案1: 一键脚本部署（推荐）

### 适用场景
- 新项目需要从零搭建AI文档管理系统
- 需要标准的目录结构和工具
- 快速在多个项目中复制相同的文档管理系统

### 部署步骤

#### 1. 获取部署脚本

方式A - 从现有项目复制：
```bash
# 从已有项目复制脚本
cp /path/to/existing/ai-docs/setup-ai-docs.sh ~/setup-ai-docs.sh
```

方式B - 下载脚本（如果托管在git仓库）：
```bash
# 下载单个脚本文件
curl -O https://your-repo/ai-docs/setup-ai-docs.sh
# 或
wget https://your-repo/ai-docs/setup-ai-docs.sh
```

#### 2. 执行部署

```bash
# 在项目根目录执行
cd /path/to/new/project
bash ~/setup-ai-docs.sh

# 或指定自定义目录
bash ~/setup-ai-docs.sh ./docs/ai-docs
```

#### 3. 验证部署

```bash
# 检查目录结构
tree ai-docs -L 2

# 测试创建任务
cd ai-docs/templates
./new-task.sh "测试任务" "你的名字"

# 查看创建的任务
ls -la ../current/测试任务/
```

#### 4. 添加到版本控制

```bash
git add ai-docs/
git commit -m "docs: 初始化AI文档管理系统"
```

### 部署后的目录结构

```
ai-docs/
├── README.md                       # 系统说明
├── QUICKREF.md                     # 快速参考
├── current/                        # 进行中的任务（空）
├── archive/                        # 已完成任务（空）
└── templates/
    ├── task-template.md            # 任务文档模板
    ├── context-resume-template.md  # 上下文恢复模板
    ├── ai-guidelines.md            # AI协作预设指令
    ├── new-task.sh                 # 创建任务脚本
    ├── export-context.sh           # 导出上下文脚本
    └── archive-task.sh             # 归档任务脚本
```

---

## 方案2: 直接复制目录

### 适用场景
- 已有配置完善的ai-docs系统
- 需要保留所有自定义内容和历史任务
- 在团队内统一文档管理规范

### 部署步骤

#### 1. 打包现有ai-docs

在源项目中：
```bash
cd /path/to/source/project

# 方式A: 只打包模板和工具（推荐新项目）
tar -czf ai-docs-template.tar.gz \
  ai-docs/templates/ \
  ai-docs/QUICKREF.md \
  ai-docs/USAGE.md \
  ai-docs/README.md \
  ai-docs/setup-ai-docs.sh

# 方式B: 打包完整系统（包含历史任务）
tar -czf ai-docs-full.tar.gz ai-docs/

# 方式C: 使用提供的打包脚本
cd ai-docs
./pack-ai-docs.sh
```

#### 2. 传输到新环境

```bash
# 本地复制
cp ai-docs-template.tar.gz /path/to/new/project/

# 远程传输
scp ai-docs-template.tar.gz user@remote:/path/to/new/project/

# 云存储
# 上传到云盘后在新环境下载
```

#### 3. 解压部署

```bash
cd /path/to/new/project

# 解压
tar -xzf ai-docs-template.tar.gz

# 如果需要，清空current和archive目录
rm -rf ai-docs/current/*
rm -rf ai-docs/archive/*

# 验证
tree ai-docs -L 2
```

#### 4. 初始化新环境

```bash
# 验证脚本权限
cd ai-docs/templates
ls -la *.sh

# 如果权限不对，重新设置
chmod +x *.sh

# 测试创建任务
./new-task.sh "测试任务" "你的名字"
```

#### 5. 添加到版本控制

```bash
git add ai-docs/
git commit -m "docs: 部署AI文档管理系统"
```

### 注意事项

1. **清理历史任务**：如果打包了完整系统，记得清理不需要的历史任务
2. **检查脚本权限**：某些压缩格式可能不保留执行权限
3. **更新文档**：根据新项目需求更新ai-guidelines.md中的项目特定规范

---

## 方案3: 手动逐步创建

### 适用场景
- 需要深度定制文档系统
- 只需要部分功能
- 学习和理解系统结构

### 部署步骤

#### 1. 创建目录结构

```bash
cd /path/to/project
mkdir -p ai-docs/{current,archive,templates}
```

#### 2. 创建核心模板文件

```bash
cd ai-docs/templates

# 创建任务模板
touch task-template.md

# 创建上下文恢复模板
touch context-resume-template.md

# 创建AI协作指令
touch ai-guidelines.md
```

#### 3. 填写模板内容

参考现有项目或USAGE.md中的模板内容，手动填写各个文件。

#### 4. 创建工具脚本

```bash
# 创建三个核心脚本
touch new-task.sh export-context.sh archive-task.sh
chmod +x *.sh

# 编写脚本内容（参考现有脚本）
```

#### 5. 创建文档

```bash
cd ai-docs
touch README.md QUICKREF.md
# 填写文档内容
```

---

## 快速部署检查清单

部署完成后，使用以下检查清单确保系统正常：

### 基础检查

- [ ] 目录结构正确（current/archive/templates）
- [ ] 所有脚本有执行权限（new-task.sh等）
- [ ] 文档文件存在（QUICKREF.md等）
- [ ] 模板文件完整（task-template.md等）

### 功能检查

- [ ] 可以创建新任务：`./new-task.sh "测试" "名字"`
- [ ] 任务文档正常生成
- [ ] 脚本可以正常运行
- [ ] Git可以正常追踪文件

### 内容检查

- [ ] ai-guidelines.md包含四个核心规范
- [ ] 脚本中的路径引用正确
- [ ] 模板中的占位符格式正确

---

## 团队部署方案

### 场景：在团队内推广使用

#### 1. 创建共享资源

```bash
# 在共享位置创建标准模板
mkdir -p /shared/ai-docs-template
cp -r ai-docs/* /shared/ai-docs-template/

# 或创建Git仓库
git init ai-docs-template
cd ai-docs-template
cp -r /path/to/ai-docs/* .
git add .
git commit -m "feat: AI文档管理系统模板"
```

#### 2. 团队成员使用

```bash
# 方式A: 从共享位置复制
cp -r /shared/ai-docs-template /path/to/my/project/ai-docs

# 方式B: 从Git仓库克隆
cd /path/to/my/project
git clone /shared/repos/ai-docs-template ai-docs
```

#### 3. 版本更新管理

```bash
# 更新共享模板
cd /shared/ai-docs-template
# 更新文件...
git add .
git commit -m "docs: 更新模板"

# 团队成员同步更新
cd /path/to/my/project/ai-docs/templates
cp /shared/ai-docs-template/templates/*.md .
cp /shared/ai-docs-template/templates/*.sh .
chmod +x *.sh
```

---

## 常见问题

### Q1: 脚本执行权限问题

```bash
# 如果提示 Permission denied
cd ai-docs/templates
chmod +x *.sh
```

### Q2: 路径引用错误

脚本使用相对路径，确保：
- 脚本位于 `ai-docs/templates/` 目录
- 执行时在 `templates` 目录内

### Q3: 如何更新已部署的系统

```bash
# 只更新模板和脚本，不影响现有任务
cd ai-docs/templates
cp /path/to/new/templates/*.md .
cp /path/to/new/templates/*.sh .
chmod +x *.sh
```

### Q4: 如何迁移现有任务到新系统

```bash
# 1. 部署新系统
bash setup-ai-docs.sh ./ai-docs-new

# 2. 迁移任务
cp -r ai-docs-old/current/* ai-docs-new/current/
cp -r ai-docs-old/archive/* ai-docs-new/archive/

# 3. 替换旧系统
mv ai-docs-old ai-docs-backup
mv ai-docs-new ai-docs
```

---

## 部署后的首次使用

### 1. 创建第一个任务

```bash
cd ai-docs/templates
./new-task.sh "AI文档系统测试" "你的名字"
```

### 2. 编辑任务文档

```bash
vim ../current/AI文档系统测试/AI文档系统测试.md
```

### 3. 与AI开始协作

提供给AI：
```
请遵循以下规范进行协作：

1. 阅读AI协作指令: ai-docs/templates/ai-guidelines.md
2. 阅读任务文档: ai-docs/current/AI文档系统测试/AI文档系统测试.md

关键规范：
- 禁止使用表情符号
- 每完成一个步骤后：更新进度文档 + git commit
- 所有任务完成后：填写完成总结 + 归档到archive目录
- 所有代码文件不超过500行，超过即模块化重构

当前需要：测试AI文档管理系统的功能
```

---

## 推荐部署方案

根据不同场景的推荐：

1. **个人新项目** → 方案1（一键脚本）
2. **团队统一规范** → 方案2（打包复制）+ Git仓库管理
3. **定制化需求** → 方案3（手动创建）
4. **快速测试** → 方案1（一键脚本）

---

## 总结

- **最快方式**：使用setup-ai-docs.sh一键部署
- **最灵活**：直接复制已配置的ai-docs目录
- **最可控**：手动创建并定制
- **团队推荐**：建立Git仓库统一管理模板

选择适合你的方案，快速开始使用AI文档管理系统！
