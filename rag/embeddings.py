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
        """批量嵌入（提高效率）"""
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

    @classmethod
    def embed_batch_with_cache(cls, texts: list[str]) -> list[list[float]]:
        """批量带缓存的嵌入生成（提高效率）"""
        embeddings = []
        uncached_texts = []
        uncached_indices = []

        for i, text in enumerate(texts):
            cache_key = hashlib.md5(text.encode()).hexdigest()
            cache_file = RAG_CACHE_DIR / f"{cache_key}.pkl"
            if cache_file.exists():
                with open(cache_file, "rb") as f:
                    embeddings.append(pickle.load(f))
            else:
                uncached_texts.append(text)
                uncached_indices.append(i)

        if uncached_texts:
            uncached_embeddings = cls.embed(uncached_texts)
            for i, j in enumerate(uncached_indices):
                embeddings.insert(j, uncached_embeddings[i])
                # 缓存嵌入
                cache_key = hashlib.md5(uncached_texts[i].encode()).hexdigest()
                cache_file = RAG_CACHE_DIR / f"{cache_key}.pkl"
                with open(cache_file, "wb") as f:
                    pickle.dump(uncached_embeddings[i], f)

        return embeddings
