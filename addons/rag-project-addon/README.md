# yc-cli RAG Project Addon

本目录用于把 yc-cli RAG（FastAPI, 默认 8733）接入到一个**新项目**的根目录。

## 你会得到什么

- `.rag/config.toml`：项目级配置（建议相对路径）
- `scripts/yc_rag.sh`：项目级 wrapper（start/status/index/search/record）
- `.cursor/skills/yc-cli-rag/`：Cursor 项目级 Skill（让 AI 自动使用 RAG）
- `ai-docs/current/your-project/`：最小 ai-docs 模板（无任何真实数据）
- `gitignore.append.txt`：建议追加到项目 `.gitignore` 的忽略规则

## 使用步骤（新项目）

1) 把本目录内容复制到你的新项目根目录（保持目录结构不变）。

2) 在新项目根目录修改 `.rag/config.toml`：

- 把 `task.name` 改成你的项目名
- `docs.base_dir` 保持 `ai-docs/current`（推荐）

3) 在新项目根目录（先激活 mamba 环境）：

```bash
mamba activate rag_cli
./scripts/yc_rag.sh start
./scripts/yc_rag.sh index
./scripts/yc_rag.sh search "决策" -n 5
```

4) 把 `gitignore.append.txt` 的内容追加到你的项目 `.gitignore`，避免提交索引与日志。

## 依赖约定

- yc-cli 代码默认放在 `~/yc-cli`。
- 如不在该路径，请在运行前设置：`export YC_CLI_HOME=/path/to/yc-cli`。
