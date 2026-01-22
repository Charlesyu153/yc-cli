# rag/server.py
from fastapi import FastAPI, Query, HTTPException, status
from pydantic import BaseModel
from typing import List, Optional
import logging

from .config import cfg
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
