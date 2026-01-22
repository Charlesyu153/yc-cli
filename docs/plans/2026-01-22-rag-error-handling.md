# RAG 系统错误处理

**关联**: [2026-01-22-rag-system-design.md](./2026-01-22-rag-system-design.md)

---

## 错误分类与处理

### 1. 启动时错误

| 错误 | 处理 |
|------|------|
| 端口被占用 | 自动尝试 +1 ~ +5 端口，失败则退出 |
| 模型下载失败 | 提示使用 `--no-model` 模式或检查网络 |
| 配置文件损坏 | 使用默认配置，记录警告日志 |
| 目录无权限 | 提示并退出 |

### 2. 运行时错误

| 错误 | 处理 |
|------|------|
| 索引为空 | 返回空结果 + 提示运行 `POST /api/index` |
| 文件不存在 | 跳过，记录日志 |
| 嵌入超时 | 重试 3 次，失败后跳过 |
| ChromaDB 错误 | 返回 500 + 错误详情 |

### 3. API 错误响应格式

```json
{
    "error": {
        "code": "INDEX_EMPTY",
        "message": "No documents indexed. Run POST /api/index first.",
        "details": {}
    }
}
```

---

## 错误代码

| 代码 | HTTP | 说明 |
|------|------|------|
| `INDEX_EMPTY` | 404 | 索引为空 |
| `MODEL_NOT_LOADED` | 503 | 模型未加载 |
| `FILE_NOT_FOUND` | 404 | 文件不存在 |
| `INVALID_INPUT` | 400 | 输入参数无效 |
| `EMBEDDING_TIMEOUT` | 504 | 嵌入生成超时 |
| `PORT_OCCUPIED` | - | 端口被占用（启动时） |

---

## 更新后的 server.py（错误处理）

```python
# rag/server.py
from fastapi import FastAPI, Query, HTTPException, status
from pydantic import BaseModel
from typing import List, Optional
import logging

from .config import cfg, PORT
from .indexer import index_all
from .vector_store import VectorStore
from .recorder import save_record
from .file_watcher import start_watcher

# 配置日志
logging.basicConfig(
    level=getattr(logging, cfg.log_level),
    filename=cfg.log_file,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

app = FastAPI(title="YC-CLI RAG Service")
_watcher = None


class ErrorResponse(BaseModel):
    error: dict


@app.exception_handler(Exception)
async def global_exception_handler(request, exc):
    logger.error(f"Unhandled exception: {exc}", exc_info=True)
    return {
        "error": {
            "code": "INTERNAL_ERROR",
            "message": str(exc)
        }
    }


@app.on_event("startup")
def startup_event():
    global _watcher

    # 检查模型
    try:
        from .embeddings import Embeddings
        Embeddings.get_model()
        logger.info("Embedding model loaded")
    except Exception as e:
        logger.error(f"Failed to load model: {e}")
        raise RuntimeError(f"Failed to load embedding model: {e}")

    # 启动文件监听
    if cfg.auto_index:
        _watcher = start_watcher()
        logger.info("File watcher started")

    logger.info(f"RAG service started on port {cfg.port}")


@app.on_event("shutdown")
def shutdown_event():
    if _watcher:
        _watcher.stop()
    logger.info("RAG service stopped")


class IndexResponse(BaseModel):
    count: int
    files: list[dict]


class RecordRequest(BaseModel):
    task: str
    title: str
    type: str  # lesson | failure | decision
    conclusion: str
    background: str
    key_points: List[str]
    tags: str = ""


@app.post("/api/index", response_model=IndexResponse)
def rebuild_index() -> IndexResponse:
    """重建索引"""
    try:
        results = index_all()
        logger.info(f"Index rebuilt: {len(results)} files")
        return IndexResponse(count=len(results), files=results)
    except Exception as e:
        logger.error(f"Index rebuild failed: {e}")
        raise HTTPException(
            status_code=500,
            detail={"code": "INDEX_FAILED", "message": str(e)}
        )


@app.get("/api/search")
def search(
    q: str = Query(..., description="搜索查询"),
    top_k: int = Query(5, ge=1, le=20, description="返回数量")
) -> list[dict]:
    """语义搜索"""
    # 检查索引状态
    status_info = VectorStore.status()
    if status_info["count"] == 0:
        raise HTTPException(
            status_code=404,
            detail={
                "code": "INDEX_EMPTY",
                "message": "No documents indexed. Run POST /api/index first."
            }
        )

    try:
        results = VectorStore.search(q, top_k)
        return [
            {
                "id": r["id"],
                "score": round(r["score"], 3),
                "title": r["metadata"].get("title", "Untitled"),
                "file_path": r["metadata"]["file_path"]
            }
            for r in results
        ]
    except Exception as e:
        logger.error(f"Search failed: {e}")
        raise HTTPException(
            status_code=500,
            detail={"code": "SEARCH_FAILED", "message": str(e)}
        )


@app.post("/api/record")
def record(req: RecordRequest) -> dict:
    """保存 AI 生成的记录"""
    try:
        result = save_record(
            task=req.task or cfg.task_name,
            title=req.title,
            record_type=req.type,
            conclusion=req.conclusion,
            background=req.background,
            key_points=req.key_points,
            tags=req.tags
        )
        logger.info(f"Record saved: {result['file_path']}")
        return {"success": True, "result": result}
    except Exception as e:
        logger.error(f"Record save failed: {e}")
        raise HTTPException(
            status_code=400,
            detail={"code": "RECORD_FAILED", "message": str(e)}
        )


@app.get("/api/status")
def status() -> dict:
    """索引状态"""
    return VectorStore.status()


@app.get("/health")
def health() -> dict:
    """健康检查"""
    return {"status": "ok", "port": cfg.port}


if __name__ == "__main__":
    import uvicorn

    # 端口占用检测
    for offset in range(6):  # 尝试原端口 + 5
        try:
            uvicorn.run(app, host=cfg.host, port=cfg.port + offset)
            break
        except OSError as e:
            if "address already in use" in str(e).lower():
                if offset == 5:
                    logger.error("No available ports")
                    raise
                continue
            raise
```
