# Git Push 前安全审查清单

**重要**: 每次 `git push` 前必须执行此检查。

---

## 1. 审查待推送内容

```bash
# 查看未推送的提交
git log origin/main..HEAD --oneline

# 查看每个提交的详细差异
git log origin/main..HEAD -p
```

---

## 2. 安全/隐私检查项

| 检查项 | 命令/方法 | 说明 |
|--------|-----------|------|
| **敏感信息** | `git diff -G "password\|secret\|key\|token" origin/main..HEAD` | 检查是否包含密码、密钥等 |
| **API密钥** | `git diff -G "sk-\|api_key\|API_KEY" origin/main..HEAD` | 检查 API 密钥泄露 |
| **个人路径** | `git diff -G "/home/" origin/main..HEAD` | 检查硬编码的本地路径 |
| **内网地址** | `git diff -G "192\.168\|10\.|172\." origin/main..HEAD` | 检查内网IP地址 |
| **配置文件** | 检查 `.rag/config.toml`, `.env` 等是否被提交 | 配置文件应被 gitignore |
| **大文件** | `git log --stat origin/main..HEAD` | 检查是否意外提交大文件 |

---

## 3. 配置文件审查

确保以下文件在 `.gitignore` 中：

```
.rag/config.toml    # 项目特定配置
.env                # 环境变量
*.log               # 日志文件
```

---

## 4. 推送前最终检查

```bash
# 一键检查脚本
git push --dry-run -v  # 预演推送，不实际执行
```

---

## 典型问题案例

### ❌ 不应推送
```toml
# .rag/config.toml
base_dir = "/home/jacekyu/PCF/ai-docs/current"  # 个人路径
name = "pcf"                                      # 项目特定
```

### ✅ 应该是
```toml
# .rag/config.toml (默认模板)
base_dir = "ai-docs/current"  # 相对路径
name = "default"              # 通用默认值
```

---

**记住**: 一旦推送到 GitHub，删除敏感信息很难彻底（git历史保留）。
