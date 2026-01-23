---
name: rag-management
description: Use this skill when user says "rag skill", needs to manage RAG index, search documents, or generate records. This skill provides a comprehensive interface for managing the RAG system.
---

# RAG Management Skill

## Quick Start

When user says "rag skill" or asks to manage RAG system:

1. **Read this skill** for the format
2. **Use the RAG CLI** or **Python API** to manage the system
3. **Generate records** using the automatic record generation feature

## When to Use

| User says | Do this |
|----------|---------|
| "rag skill" | Read this skill, follow format |
| "搜索文档" | Use RAG search functionality |
| "生成记录" | Use RAG record generation feature |
| "重建索引" | Use RAG index rebuild functionality |

## Core Principle

**RAG system provides semantic search and automatic record generation for project insights.**

- **Index**: Only index `MANIFEST.md` and `insights/*.md` files
- **Search**: Semantic search using embeddings
- **Record**: Automatic record generation with templates

## File Structure

```
rag/
├── config.py          # Configuration management
├── embeddings.py      # Embedding generation
├── vector_store.py    # Vector storage (ChromaDB)
├── indexer.py         # File indexing
├── file_watcher.py    # File watching
├── recorder.py        # Record generation
├── server.py          # HTTP service
└── ragctl.py          # CLI tool
```

## Usage

### 1. CLI Tool

```bash
# Status
ragctl status

# Rebuild index
ragctl index

# Search
ragctl search "decision"

# Create record
ragctl record

# Export/Import
ragctl export backup.tar.gz
ragctl import backup.tar.gz

# Logs
ragctl logs -f

# Start/Stop service
ragctl start
ragctl stop
```

### 2. Python API

```python
from rag.indexer import index_all, index_file
from rag.vector_store import VectorStore
from rag.recorder import save_record

# Index all files
results = index_all()

# Index single file
result = index_file(file_path)

# Search documents
results = VectorStore.search("decision", top_k=5)

# Generate record
result = save_record(
    task="default",
    title="Decision title",
    record_type="decision",
    conclusion="One-line conclusion",
    background="Background information",
    key_points=["Point 1", "Point 2"],
    tags="tag1,tag2"
)
```

### 3. HTTP API

```bash
# Health check
curl http://localhost:8733/health

# Status
curl http://localhost:8733/api/status

# Rebuild index
curl -X POST http://localhost:8733/api/index

# Search
curl "http://localhost:8733/api/search?q=decision&top_k=5"

# Create record
curl -X POST http://localhost:8733/api/record \
  -H "Content-Type: application/json" \
  -d '{
    "task": "default",
    "title": "Decision title",
    "type": "decision",
    "conclusion": "One-line conclusion",
    "background": "Background information",
    "key_points": ["Point 1", "Point 2"]
  }'
```

## Configuration

Configuration file: `.rag/config.toml`

```toml
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
```

## Record Types

| Type | Description | Template |
|------|-------------|----------|
| `lesson` | Key lessons learned | `ai-docs/templates/lesson-template.md` |
| `failure` | Failure analysis | `ai-docs/templates/failure-template.md` |
| `decision` | Decision records | `ai-docs/templates/decision-template.md` |

## Optimization

- **Model caching**: Embedding model is loaded once and cached
- **Embedding caching**: Embeddings are cached to avoid recomputation
- **Batch indexing**: Use batch operations for better performance
- **Incremental indexing**: Only index changed files

## Workflow

### Creating a New Insight

1. Create insight file: `ai-docs/current/{task}/insights/{insight_id}.md`
2. Run index: `ragctl index`
3. Verify: `ragctl search "keyword"`

### Generating Records

1. Use CLI: `ragctl record`
2. Fill in the required fields
3. Review generated file: `ai-docs/current/{task}/insights/{insight_id}.md`

### Backing Up Data

1. Export: `ragctl export backup.tar.gz`
2. Import: `ragctl import backup.tar.gz`
3. Rebuild index: `ragctl index`

## Rules

1. **Only index MANIFEST.md and insights/*.md files** - Keep the index clean and focused
2. **Use pointers for code references** - Don't index code directly, use `related_files` pointers
3. **Run validation before indexing** - Ensure insight files are valid
4. **Regular backups** - Export data regularly to avoid data loss

## What This Skill Is NOT

- ❌ NOT for indexing entire codebase
- ❌ NOT for file-level code indexing
- ✅ This is for project insights, decisions, and state tracking (with semantic search)

## References

See `skills/context-insight/SKILL.md` for more details on insight format and structure.
