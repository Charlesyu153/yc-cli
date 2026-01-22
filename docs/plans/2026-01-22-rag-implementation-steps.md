# RAG 系统实现步骤

**关联**: [2026-01-22-rag-implementation-plan.md](./2026-01-22-rag-implementation-plan.md)

---

## 步骤 1：项目初始化

```bash
# 1. 创建目录
mkdir -p rag/.rag/cache
mkdir -p ai-docs/current/default/insights

# 2. 创建 requirements.txt
cat > requirements.txt << 'EOF'
fastapi>=0.104.0
uvicorn[standard]>=0.24.0
chromadb>=0.4.18
sentence-transformers>=2.2.2
watchdog>=3.0.0
python-dotenv>=1.0.0
toml>=0.10.2
pyyaml>=6.0
requests>=2.31.0
EOF

# 3. 创建默认配置
mkdir -p .rag
cat > .rag/config.toml << 'EOF'
[server]
host = "127.0.0.1"
port = 8733

[embedding]
model = "sentence-transformers/all-MiniLM-L6-v2"
cache_size_mb = 500

[docs]
base_dir = "ai-docs/current"
auto_index = true

[task]
name = "default"

[log]
level = "INFO"
file = ".rag/rag.log"
EOF

# 4. 安装依赖
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

---

## 步骤 2：配置模块 (config.py)

详见: [2026-01-22-rag-config-design.md](./2026-01-22-rag-config-design.md)

---

## 步骤 3：嵌入模块 (embeddings.py)

```python
# rag/embeddings.py
from sentence_transformers import SentenceTransformer
import pickle
import hashlib
import logging

from .config import cfg, RAG_CACHE_DIR

logger = logging.getLogger(__name__)


class Embeddings:
    _model = None

    @classmethod
    def get_model(cls):
        if cls._model is None:
            logger.info(f"Loading embedding model: {cfg.embedding_model}")
            cls._model = SentenceTransformer(cfg.embedding_model)
        return cls._model

    @classmethod
    def embed(cls, texts: list[str]) -> list[list[float]]:
        return cls.get_model().encode(texts, convert_to_tensor=False).tolist()

    @classmethod
    def embed_with_cache(cls, text: str) -> list[float]:
        """带缓存的嵌入生成"""
        cache_key = hashlib.md5(text.encode()).hexdigest()
        cache_file = RAG_CACHE_DIR / f"{cache_key}.pkl"

        if cache_file.exists():
            with open(cache_file, "rb") as f:
                return pickle.load(f)

        embedding = cls.embed([text])[0]
        with open(cache_file, "wb") as f:
            pickle.dump(embedding, f)
        return embedding
```

---

## 步骤 4：向量存储模块 (vector_store.py)

```python
# rag/vector_store.py
import chromadb
from chromadb.config import Settings
from datetime import datetime
import logging

from .config import cfg, CHROMA_DIR
from .embeddings import Embeddings

logger = logging.getLogger(__name__)


class VectorStore:
    _client = None
    _collection = None

    @classmethod
    def get_collection(cls):
        if cls._client is None:
            cls._client = chromadb.PersistentClient(
                path=str(CHROMA_DIR),
                settings=Settings(anonymized_telemetry=False)
            )
            cls._collection = cls._client.get_or_create_collection(
                name="yc_insights",
                metadata={"hnsw:space": "cosine"}
            )
        return cls._collection

    @classmethod
    def upsert(cls, id: str, content: str, metadata: dict):
        """添加/更新单条记录（稳定 id，避免重复与过期记忆）"""
        embedding = Embeddings.embed_with_cache(content)
        cls.get_collection().upsert(
            ids=[id],
            embeddings=[embedding],
            documents=[content],
            metadatas=[metadata]
        )

    @classmethod
    def delete(cls, id: str):
        """删除单条记录（用于文件删除/重命名）"""
        cls.get_collection().delete(ids=[id])

    @classmethod
    def search(cls, query: str, top_k: int = 5) -> list[dict]:
        """语义搜索"""
        embedding = Embeddings.embed_with_cache(query)
        results = cls.get_collection().query(
            query_embeddings=[embedding],
            n_results=top_k
        )

        return [
            {
                "id": results["ids"][0][i],
                "score": 1 - results["distances"][0][i],
                "content": results["documents"][0][i],
                "metadata": results["metadatas"][0][i]
            }
            for i in range(len(results["ids"][0]))
        ]

    @classmethod
    def clear(cls):
        """清空索引"""
        if cls._client:
            cls._client.delete_collection("yc_insights")
            cls._collection = None

    @classmethod
    def status(cls) -> dict:
        """索引状态"""
        coll = cls.get_collection()
        return {
            "count": coll.count(),
            "last_updated": datetime.now().isoformat()
        }
```

---

## 步骤 5：索引模块 (indexer.py)

```python
# rag/indexer.py
import re
from pathlib import Path
from datetime import datetime
import logging
import yaml

from .config import cfg
from .vector_store import VectorStore

logger = logging.getLogger(__name__)


def parse_frontmatter(content: str) -> tuple[dict, str]:
    """解析 Markdown YAML frontmatter（支持 list/dict，用于 related_files 指针）"""
    match = re.match(r'^---\n(.*?)\n---\n(.*)', content, re.DOTALL)
    if not match:
        return {}, content

    frontmatter_text, body = match.groups()
    metadata = yaml.safe_load(frontmatter_text) or {}
    return metadata, body


def doc_id_for(file_path: Path) -> str:
    """稳定 ID：使用 docs_dir 的相对路径（跨机器可复现，支持 upsert/delete）"""
    return file_path.relative_to(cfg.docs_dir).as_posix()


def index_file(file_path: Path) -> dict:
    """索引单个文件"""
    try:
        content = file_path.read_text()
    except Exception as e:
        logger.warning(f"Failed to read {file_path}: {e}")
        return {"error": str(e)}

    metadata, body = parse_frontmatter(content)

    # 构建索引内容
    related_files = metadata.get("related_files") or []
    related_paths = [rf.get("path") for rf in related_files if isinstance(rf, dict) and rf.get("path")]
    related_paths_str = ",".join(related_paths)
    index_content = f"{metadata.get('title', '')}\n\nPointers: {related_paths_str}\n\n{body[:8000]}"

    # 生成稳定 ID（避免 timestamp 造成重复与过期）
    doc_id = doc_id_for(file_path)

    # 添加元数据
    metadata.update({
        "file_path": str(file_path),
        "task": file_path.parent.parent.name,
        "created": metadata.get("created") or datetime.now().isoformat(),
        "has_pointers": bool(related_paths),
        "related_paths": related_paths_str,
    })

    VectorStore.upsert(doc_id, index_content, metadata)
    logger.debug(f"Indexed: {file_path}")
    return {"id": doc_id, "file": str(file_path)}


def index_all():
    """索引所有文档"""
    VectorStore.clear()

    results = []
    docs_dir = cfg.docs_dir

    if not docs_dir.exists():
        logger.warning(f"Docs directory not found: {docs_dir}")
        return results

    for task_dir in docs_dir.iterdir():
        if not task_dir.is_dir():
            continue

        # 索引 MANIFEST.md
        manifest = task_dir / "MANIFEST.md"
        if manifest.exists():
            results.append(index_file(manifest))

        # 索引 insights/
        insights_dir = task_dir / "insights"
        if insights_dir.exists():
            for insight_file in insights_dir.glob("*.md"):
                results.append(index_file(insight_file))

    logger.info(f"Indexed {len(results)} files")
    return results
```

---

## 步骤 6：文件监听模块 (file_watcher.py)

```python
# rag/file_watcher.py
import time
from pathlib import Path
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import logging

from .config import cfg
from .indexer import index_file, doc_id_for
from .vector_store import VectorStore

logger = logging.getLogger(__name__)


class InsightFileHandler(FileSystemEventHandler):
    """监听 insight 文件变化"""

    def on_modified(self, event):
        if event.src_path.endswith('.md'):
            time.sleep(cfg.index_debounce)
            file_path = Path(event.src_path)
            if file_path.exists():
                index_file(file_path)

    def on_created(self, event):
        if event.src_path.endswith('.md'):
            time.sleep(cfg.index_debounce)
            file_path = Path(event.src_path)
            if file_path.exists():
                index_file(file_path)

    def on_deleted(self, event):
        if event.src_path.endswith('.md'):
            try:
                VectorStore.delete(doc_id_for(Path(event.src_path)))
            except Exception:
                # 删除失败不阻塞；下次全量 index 会自愈
                pass


def start_watcher():
    """启动文件监听"""
    if not cfg.docs_dir.exists():
        logger.warning(f"Docs directory not found: {cfg.docs_dir}")
        return None

    event_handler = InsightFileHandler()
    observer = Observer()
    observer.schedule(event_handler, str(cfg.docs_dir), recursive=True)
    observer.start()
    logger.info("File watcher started")
    return observer
```

---

## 步骤 7：记录器模块 (recorder.py)

```python
# rag/recorder.py
from pathlib import Path
from datetime import datetime
import logging

from .config import cfg
from .indexer import index_file

logger = logging.getLogger(__name__)


RECORD_TEMPLATE = """---
title: {title}
type: {record_type}
tags: {tags}
created: {created}
---

## 结论
{conclusion}

## 背景
{background}

## 关键点
{key_points}
"""


def save_record(
    task: str,
    title: str,
    record_type: str,
    conclusion: str,
    background: str,
    key_points: list[str],
    tags: str = ""
) -> dict:
    """保存记录到文件并更新索引"""

    task = task or cfg.task_name
    task_dir = cfg.docs_dir / task
    insights_dir = task_dir / "insights"
    insights_dir.mkdir(parents=True, exist_ok=True)

    # 生成文件名
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{record_type}_{timestamp}.md"
    file_path = insights_dir / filename

    # 格式化关键点
    key_points_text = "\n".join(f"{i+1}. {p}" for i, p in enumerate(key_points))

    # 填充模板
    content = RECORD_TEMPLATE.format(
        title=title,
        record_type=record_type,
        tags=tags or "auto-generated",
        created=datetime.now().isoformat(),
        conclusion=conclusion,
        background=background,
        key_points=key_points_text
    )

    # 写入文件
    file_path.write_text(content)
    logger.info(f"Record saved: {file_path}")

    # 更新索引
    result = index_file(file_path)

    return {
        "file_path": str(file_path),
        "id": result["id"],
        "title": title
    }
```

---

## 步骤 8：服务模块 (server.py)

详见: [2026-01-22-rag-error-handling.md](./2026-01-22-rag-error-handling.md)

---

## 步骤 9：CLI 工具 (ragctl.py)

详见: [2026-01-22-rag-cli-tool.md](./2026-01-22-rag-cli-tool.md)

---

## 步骤 10：创建入口

```bash
# rag/__init__.py
# 空文件

# 创建启动脚本
cat > run_server.sh << 'EOF'
#!/bin/bash
cd "$(dirname "$0")"
source venv/bin/activate
python -m rag.server
EOF
chmod +x run_server.sh

# 安装 CLI 工具
chmod +x rag/ragctl.py
ln -sf $(pwd)/rag/ragctl.py ~/.local/bin/ragctl
```

---

## 步骤 11：测试

```bash
# 启动服务
./run_server.sh

# 或使用 CLI
ragctl start

# 另一个终端测试
ragctl status
ragctl index
ragctl search "数据库" -n 3

# 创建记录
ragctl record -t my-project

# 查看日志
ragctl logs -f
```

---

## 文件清单

```
yc-cli/
├── rag/
│   ├── __init__.py
│   ├── config.py           # 配置模块
│   ├── embeddings.py       # 嵌入生成
│   ├── vector_store.py     # 向量存储
│   ├── indexer.py          # 文件索引
│   ├── file_watcher.py     # 文件监听
│   ├── recorder.py         # 记录器
│   ├── server.py           # HTTP 服务
│   └── ragctl.py           # CLI 工具
├── .rag/
│   ├── config.toml         # 配置文件
│   ├── cache/              # 嵌入缓存
│   ├── chroma/             # 向量数据库
│   └── rag.log             # 日志文件
├── ai-docs/
│   └── current/{task}/
│       ├── MANIFEST.md
│       └── insights/
├── requirements.txt
├── run_server.sh
└── docs/plans/
    ├── 2026-01-22-rag-system-design.md
    ├── 2026-01-22-rag-implementation-plan.md
    ├── 2026-01-22-rag-implementation-steps.md
    ├── 2026-01-22-rag-config-design.md
    ├── 2026-01-22-rag-error-handling.md
    ├── 2026-01-22-rag-cli-tool.md
    └── 2026-01-22-rag-session-startup.md
```

---

## 快速开始

```bash
# 1. 初始化
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt

# 2. 启动服务
./run_server.sh

# 3. 建立索引
ragctl index

# 4. 搜索测试
ragctl search "关键词"
```

---

## 迁移

```bash
# 打包
tar czf rag-backup.tar.gz .rag/ ai-docs/

# 恢复
tar xzf rag-backup.tar.gz
ragctl index
```
