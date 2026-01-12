# Repository Guidelines

## AI Agents Configuration

本项目配置了三个专用 AI agents，参考 feature-dev 插件的功能模式，用于辅助开发工作流程。

### 配置文件位置

- 主配置：`.claude/agents/agents.json`
- Agent 配置：`.claude/agents/`
  - `code-architect.json` - 代码架构师
  - `code-explorer.json` - 代码探索者
  - `code-reviewer.json` - 代码审查者

### Agent 角色说明

#### 1. code-architect (代码架构师)

负责架构设计、技术选型和实施规划。

使用场景：
- 设计新功能的技术架构
- 评估不同实现方案的优劣
- 制定详细的实施步骤
- 重构超过 500 行的代码
- 优化代码组织结构

调用方式：
```
请以 code-architect 的角色，帮我 [具体需求]
```

#### 2. code-explorer (代码探索者)

负责快速查找、理解代码结构和生成文档。

使用场景：
- 查找特定功能的实现位置
- 理解项目的代码结构
- 生成代码文档
- 追踪函数调用链
- 分析模块依赖关系

调用方式：
```
请以 code-explorer 的角色，帮我 [具体需求]
```

#### 3. code-reviewer (代码审查者)

负责代码质量检查、发现问题和提供改进建议。

使用场景：
- 审查代码质量
- 检查是否符合项目规范
- 发现潜在的 bug 和问题
- 提供代码改进建议
- 验证文件行数限制（500 行）

调用方式：
```
请以 code-reviewer 的角色，审查 [文件/目录]
```

### Feature-Dev 工作流程

完整的功能开发流程（类似 feature-dev 插件）：

1. **探索阶段** (code-explorer)
   - 理解现有代码结构和相关功能
   - 查找相关代码位置

2. **设计阶段** (code-architect)
   - 设计新功能的架构和实施方案
   - 制定详细的实施步骤

3. **实现阶段**
   - 按照设计方案实现功能代码
   - 遵循项目编码规范

4. **审查阶段** (code-reviewer)
   - 审查代码质量和规范符合性
   - 发现并修复问题

5. **提交阶段**
   - 根据审查意见修改
   - 更新文档并提交代码

### 使用示例

详细使用示例请参考：`.claude/agents/USAGE_EXAMPLES.md`

---

## Project Structure & Module Organization
- `PCF/`: Python pipelines that build annotated `.h5ad` and `.sdata.zarr` from `cleandata/` (e.g., `build_spatialdata_batch.py`, `build_spatialdata_single.py`).
- Root R scripts (`convert_h5ad_to_seurat.R`, `test_single_sample.R`, `test_conversion.R`) handle h5ad to Seurat conversion.
- `utils/R/`: reusable R functions (QC, niche analysis).
- `script/` and `script_ref/`: analysis notebooks and references; not required for core conversion.
- Data folders: `cleandata/` + `rawdata*/` inputs, `spatial/` Python outputs, `data/` Seurat outputs.

## Build, Test, and Development Commands
- `python PCF/build_spatialdata_batch.py` - process all samples under `cleandata/`.
- `python PCF/build_spatialdata_single.py` - process a single sample into `spatial/`.
- `Rscript convert_h5ad_to_seurat.R [input_dir] [output_dir]` - batch convert `.h5ad` to `.qs`.
- `Rscript test_single_sample.R` - full check on `spatial/20P.annotated.h5ad`.
- `Rscript test_conversion.R` - quick conversion sanity check (writes `.rds`).

## Coding Style & Naming Conventions
- Python: PEP 8, 4-space indents, `snake_case` for functions/vars; keep config in module-level constants (see `ROOT`, `CANONICAL_MARKERS`).
- R: 2-space indents, `snake_case` function names; keep console output via `cat()`/`glue()`.
- Prefer `Path`/`file.path` for paths; if you change the hard-coded `ROOT` in Python, update other scripts consistently.

## Testing Guidelines
- No formal test framework; rely on `test_*.R` scripts and verify outputs in `spatial/` and `data/`.
- Avoid committing large binary outputs (`.h5ad`, `.zarr`, `.qs`, `.rds`) unless explicitly requested.

## Commit & Pull Request Guidelines
- No established commit history; use short imperative messages (e.g., "Add batch h5ad conversion").
- PRs should describe data inputs/outputs, list commands run, and note any new dependencies or environment variables.

## Configuration & Environment
- R packages: `Seurat`, `SeuratObject`, `reticulate`, `qs`, `glue`.
- Python env should include `scanpy`, `spatialdata`, `anndata`; R scripts honor `PCF_PYTHON=/path/to/python3`.
- Optional flags for positivity inference: `PCF_INFER_POSITIVITY_FROM_INTENSITY=1` and `PCF_INTENSITY_POSITIVITY_THRESHOLD`.
