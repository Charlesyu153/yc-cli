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
        # 重新创建集合
        cls.get_collection()

    @classmethod
    def status(cls) -> dict:
        """索引状态"""
        coll = cls.get_collection()
        return {
            "count": coll.count(),
            "last_updated": datetime.now().isoformat()
        }
