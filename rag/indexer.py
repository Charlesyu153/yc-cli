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
    """索引单个文件（包含验证）"""
    try:
        content = file_path.read_text()
    except Exception as e:
        logger.warning(f"Failed to read {file_path}: {e}")
        return {"error": str(e)}

    metadata, body = parse_frontmatter(content)

    # 验证 insight 文件格式
    if file_path.name.endswith(".md") and file_path.parent.name == "insights":
        try:
            import subprocess
            from .config import BASE_DIR
            script_path = BASE_DIR / "scripts" / "validate_insights.py"
            result = subprocess.run(["python", str(script_path), str(file_path)], capture_output=True, text=True)
            if result.returncode != 0:
                logger.warning(f"Insight file validation failed: {file_path}\n{result.stderr}")
        except Exception as e:
            logger.warning(f"Failed to validate insight file: {file_path}: {e}")

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
        "modified": str(datetime.fromtimestamp(file_path.stat().st_mtime).isoformat()),
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
    """索引所有文档（优化版：增量索引 + 批量处理）"""
    results = []
    docs_dir = cfg.docs_dir

    if not docs_dir.exists():
        logger.warning(f"Docs directory not found: {docs_dir}")
        return results

    # 获取当前所有需要索引的文件
    all_files = []
    for task_dir in docs_dir.iterdir():
        if not task_dir.is_dir():
            continue
        manifest = task_dir / "MANIFEST.md"
        if manifest.exists():
            all_files.append(manifest)
        insights_dir = task_dir / "insights"
        if insights_dir.exists():
            all_files.extend(insights_dir.rglob("*.md"))

    # 获取当前索引中的所有文档 ID 和元数据
    current_docs = {}
    try:
        all_docs = VectorStore.search("", top_k=10000)
        for doc in all_docs:
            current_docs[doc["id"]] = doc["metadata"]
    except Exception as e:
        logger.warning(f"Failed to get current index status: {e}")
        # 如果获取索引状态失败，我们就重新索引所有文件
        docs_to_index = []
        for file_path in all_files:
            metadata, body = parse_frontmatter(file_path.read_text())
            related_files = metadata.get("related_files") or []
            related_paths = [rf.get("path") for rf in related_files if isinstance(rf, dict) and rf.get("path")]
            related_paths_str = ",".join(related_paths)
            index_content = f"{metadata.get('title', '')}\n\nPointers: {related_paths_str}\n\n{body[:8000]}"
            doc_id = doc_id_for(file_path)
            # 准备元数据
            clean_metadata = {
                "title": str(metadata.get('title', '')),
                "type": str(metadata.get('type', '')),
                "status": str(metadata.get('status', '')),
                "priority": str(metadata.get('priority', '')),
                "file_path": str(file_path),
                "task": file_path.parent.parent.name,
                "created": str(metadata.get("created") or datetime.now().isoformat()),
                "modified": str(datetime.fromtimestamp(file_path.stat().st_mtime).isoformat()),
                "has_pointers": bool(related_paths),
                "related_paths": related_paths_str,
            }
            tags = metadata.get("tags")
            if tags:
                if isinstance(tags, list):
                    clean_metadata["tags"] = ",".join(str(t) for t in tags)
                else:
                    clean_metadata["tags"] = str(tags)
            docs_to_index.append({"id": doc_id, "content": index_content, "metadata": clean_metadata})
        VectorStore.clear()
        VectorStore.upsert_batch(docs_to_index)
        logger.info(f"Indexed {len(docs_to_index)} files (full reindex)")
        return [{"id": doc["id"], "file": doc["metadata"]["file_path"]} for doc in docs_to_index]

    # 比较当前文件和索引中的文档
    files_to_index = []
    files_to_delete = []

    # 检查需要索引的文件
    for file_path in all_files:
        doc_id = doc_id_for(file_path)
        if doc_id not in current_docs:
            files_to_index.append(file_path)
        else:
            # 检查文件是否已修改
            file_modified = datetime.fromtimestamp(file_path.stat().st_mtime).isoformat()
            index_modified = current_docs[doc_id].get("modified", "")
            if file_modified > index_modified:
                files_to_index.append(file_path)

    # 检查需要删除的文档
    for doc_id in current_docs:
        file_path = cfg.docs_dir / doc_id
        if not file_path.exists():
            files_to_delete.append(doc_id)

    # 执行索引和删除操作
    if files_to_index:
        docs_to_index = []
        for file_path in files_to_index:
            try:
                content = file_path.read_text()
            except Exception as e:
                logger.warning(f"Failed to read {file_path}: {e}")
                results.append({"error": str(e)})
                continue
            metadata, body = parse_frontmatter(content)
            related_files = metadata.get("related_files") or []
            related_paths = [rf.get("path") for rf in related_files if isinstance(rf, dict) and rf.get("path")]
            related_paths_str = ",".join(related_paths)
            index_content = f"{metadata.get('title', '')}\n\nPointers: {related_paths_str}\n\n{body[:8000]}"
            doc_id = doc_id_for(file_path)
            # 准备元数据
            clean_metadata = {
                "title": str(metadata.get('title', '')),
                "type": str(metadata.get('type', '')),
                "status": str(metadata.get('status', '')),
                "priority": str(metadata.get('priority', '')),
                "file_path": str(file_path),
                "task": file_path.parent.parent.name,
                "created": str(metadata.get("created") or datetime.now().isoformat()),
                "modified": str(datetime.fromtimestamp(file_path.stat().st_mtime).isoformat()),
                "has_pointers": bool(related_paths),
                "related_paths": related_paths_str,
            }
            tags = metadata.get("tags")
            if tags:
                if isinstance(tags, list):
                    clean_metadata["tags"] = ",".join(str(t) for t in tags)
                else:
                    clean_metadata["tags"] = str(tags)
            docs_to_index.append({"id": doc_id, "content": index_content, "metadata": clean_metadata})
        VectorStore.upsert_batch(docs_to_index)
        for doc in docs_to_index:
            results.append({"id": doc["id"], "file": doc["metadata"]["file_path"]})

    if files_to_delete:
        VectorStore.delete_batch(files_to_delete)
        for doc_id in files_to_delete:
            logger.debug(f"Deleted: {doc_id}")

    logger.info(f"Indexed {len(files_to_index)} files, deleted {len(files_to_delete)} files")
    return results
