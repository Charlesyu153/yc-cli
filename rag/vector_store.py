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
    def upsert_batch(cls, docs: list[dict]):
        """批量添加/更新记录（提高效率）"""
        if not docs:
            return
        ids = [doc["id"] for doc in docs]
        contents = [doc["content"] for doc in docs]
        metadatas = [doc["metadata"] for doc in docs]
        embeddings = Embeddings.embed_batch_with_cache(contents)
        cls.get_collection().upsert(
            ids=ids,
            embeddings=embeddings,
            documents=contents,
            metadatas=metadatas
        )

    @classmethod
    def delete(cls, id: str):
        """删除单条记录（用于文件删除/重命名）"""
        cls.get_collection().delete(ids=[id])

    @classmethod
    def delete_batch(cls, ids: list[str]):
        """批量删除记录（提高效率）"""
        if not ids:
            return
        cls.get_collection().delete(ids=ids)

    @classmethod
    def search(cls, query: str, top_k: int = 5, where: dict = None) -> list[dict]:
        """语义搜索（支持元数据过滤）"""
        embedding = Embeddings.embed_with_cache(query)
        results = cls.get_collection().query(
            query_embeddings=[embedding],
            n_results=top_k,
            where=where
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

    @classmethod
    def get_all_ids(cls) -> list[str]:
        """获取所有文档的 ID（用于增量索引）"""
        try:
            # 使用一个非常通用的查询来获取所有文档
            all_docs = cls.search("", top_k=10000)
            return [doc["id"] for doc in all_docs]
        except Exception as e:
            logger.warning(f"Failed to get all document IDs: {e}")
            return []
