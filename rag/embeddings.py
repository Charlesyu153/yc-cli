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
