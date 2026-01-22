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
    """稳定 ID：使�� docs_dir 的相对路径（跨机器可复现，支持 upsert/delete）"""
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

    # 准备元数据（ChromaDB 只支持 str/int/float/bool/None）
    # 将 tags 等复杂类型转换为字符串
    clean_metadata = {
        "title": str(metadata.get('title', '')),
        "type": str(metadata.get('type', '')),
        "status": str(metadata.get('status', '')),
        "priority": str(metadata.get('priority', '')),
        "file_path": str(file_path),
        "task": file_path.parent.parent.name,
        "created": str(metadata.get("created") or datetime.now().isoformat()),
        "has_pointers": bool(related_paths),
        "related_paths": related_paths_str,
    }

    # 添加 tags 作为字符串
    tags = metadata.get("tags")
    if tags:
        if isinstance(tags, list):
            clean_metadata["tags"] = ",".join(str(t) for t in tags)
        else:
            clean_metadata["tags"] = str(tags)

    VectorStore.upsert(doc_id, index_content, clean_metadata)
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
