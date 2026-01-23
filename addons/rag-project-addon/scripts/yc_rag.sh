#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

YC_CLI_HOME="${YC_CLI_HOME:-$HOME/yc-cli}"
export PYTHONPATH="${YC_CLI_HOME}${PYTHONPATH:+:${PYTHONPATH}}"

cd "${PROJECT_ROOT}"

usage() {
  cat <<'EOF'
Usage:
  ./scripts/yc_rag.sh status
  ./scripts/yc_rag.sh index
  ./scripts/yc_rag.sh search "<query>" [-n <top_k>]
  ./scripts/yc_rag.sh record [-t <task>]
  ./scripts/yc_rag.sh start

Notes:
  - Please activate mamba env first: mamba activate rag_cli
  - Default yc-cli path is ~/yc-cli; override via YC_CLI_HOME
EOF
}

if [[ $# -lt 1 ]]; then
  usage
  exit 1
fi

cmd="$1"
shift || true

if [[ ! -d "${YC_CLI_HOME}" ]]; then
  echo "[ERR] YC_CLI_HOME not found: ${YC_CLI_HOME}" >&2
  echo "[INFO] Set YC_CLI_HOME=/path/to/yc-cli" >&2
  exit 2
fi

python -m rag.ragctl "${cmd}" "$@"
