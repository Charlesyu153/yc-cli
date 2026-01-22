# rag/config.py
from pathlib import Path
import tomllib  # Python 3.11+
import os
from dataclasses import dataclass


@dataclass
class Config:
    # 服务
    host: str = "127.0.0.1"
    port: int = 8733

    # 嵌入
    embedding_model: str = "sentence-transformers/all-MiniLM-L6-v2"
    cache_size_mb: int = 500
    embedding_timeout: int = 30

    # 文档
    docs_dir: str = "ai-docs/current"
    auto_index: bool = True
    index_debounce: float = 1.0

    # 任务
    task_name: str = "default"
    auto_create_task: bool = True

    # 日志
    log_level: str = "INFO"
    log_file: str = ".rag/rag.log"

    @classmethod
    def load(cls) -> "Config":
        """加载配置：环境变量 > config.toml > 默认值"""
        config_file = BASE_DIR / ".rag" / "config.toml"
        values = {}

        # 1. 读取配置文件
        if config_file.exists():
            with open(config_file, "rb") as f:
                toml_config = tomllib.load(f)

            if "server" in toml_config:
                values["host"] = toml_config["server"].get("host", cls.host)
                values["port"] = toml_config["server"].get("port", cls.port)

            if "embedding" in toml_config:
                values["embedding_model"] = toml_config["embedding"].get("model", cls.embedding_model)
                values["cache_size_mb"] = toml_config["embedding"].get("cache_size_mb", cls.cache_size_mb)

            if "docs" in toml_config:
                values["docs_dir"] = toml_config["docs"].get("base_dir", cls.docs_dir)
                values["auto_index"] = toml_config["docs"].get("auto_index", cls.auto_index)

            if "task" in toml_config:
                values["task_name"] = toml_config["task"].get("name", cls.task_name)

            if "log" in toml_config:
                values["log_level"] = toml_config["log"].get("level", cls.log_level)
                values["log_file"] = toml_config["log"].get("file", cls.log_file)

        # 2. 环境变量覆盖
        values["port"] = int(os.getenv("RAG_PORT", values.get("port", cls.port)))
        values["log_level"] = os.getenv("RAG_LOG_LEVEL", values.get("log_level", cls.log_level))
        values["task_name"] = os.getenv("RAG_TASK_NAME", values.get("task_name", cls.task_name))

        # 3. 解析路径
        values["docs_dir"] = BASE_DIR / values["docs_dir"]
        values["log_file"] = BASE_DIR / values["log_file"]

        return cls(**{k: v for k, v in values.items() if v is not None})


# 全局配置实例
BASE_DIR = Path(__file__).parent.parent
RAG_DIR = BASE_DIR / ".rag"
RAG_CACHE_DIR = RAG_DIR / "cache"
CHROMA_DIR = RAG_DIR / "chroma"

# 确保目录存在
for d in [RAG_DIR, RAG_CACHE_DIR, CHROMA_DIR]:
    d.mkdir(parents=True, exist_ok=True)

cfg = Config.load()

# 兼容性：导出常量
HOST = cfg.host
PORT = cfg.port
EMBEDDING_MODEL = cfg.embedding_model
DOCS_DIR = cfg.docs_dir
