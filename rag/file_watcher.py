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
