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
