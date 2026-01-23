#!/usr/bin/env python3
"""
RAG CLI Tool

Usage:
    ragctl status
    ragctl index
    ragctl search <query>
    ragctl record
    ragctl export [output.tar.gz]
    ragctl import <input.tar.gz>
    ragctl logs
    ragctl start
    ragctl stop
"""
import sys
import argparse
import requests
import subprocess
import os
from pathlib import Path
from colorama import Fore, Style

# 服务地址 - 从配置系统读取
def get_api_base():
    """获取API地址，优先使用环境变量"""
    try:
        from rag.config import cfg
        port = os.getenv("RAG_PORT", str(cfg.port))
    except (ImportError, AttributeError):
        port = os.getenv("RAG_PORT", "8733")
    return f"http://127.0.0.1:{port}"

API_BASE = get_api_base()

# 颜色输出函数
def print_success(msg):
    print(f"{Fore.GREEN}✓ {msg}{Style.RESET_ALL}")

def print_error(msg):
    print(f"{Fore.RED}✗ {msg}{Style.RESET_ALL}")

def print_warning(msg):
    print(f"{Fore.YELLOW}⚠️  {msg}{Style.RESET_ALL}")

def print_info(msg):
    print(f"{Fore.BLUE}ℹ️  {msg}{Style.RESET_ALL}")


def cmd_status(args):
    """查看状态"""
    try:
        r = requests.get(f"{API_BASE}/health", timeout=2)
        if r.status_code == 200:
            print_success("Service is running")

            r = requests.get(f"{API_BASE}/api/status")
            data = r.json()
            print(f"  Indexed documents: {data['count']}")
            print(f"  Last updated: {data['last_updated']}")
        else:
            print_error("Service is not healthy")
    except requests.exceptions.ConnectionError:
        print_error("Service is not running")
        print_info("Start with: ragctl start")
    return 0


def cmd_index(args):
    """重建索引"""
    print_info("Rebuilding index...")
    try:
        r = requests.post(f"{API_BASE}/api/index")
        data = r.json()
        print_success(f"Indexed {data['count']} files")
    except requests.exceptions.ConnectionError:
        print_error("Service is not running")
        return 1
    return 0


def cmd_search(args):
    """搜索"""
    try:
        r = requests.get(f"{API_BASE}/api/search", params={"q": args.query, "top_k": args.n})
        results = r.json()

        if not results:
            print_info("No results found")
        else:
            print_info(f"Found {len(results)} results:\n")
            for i, item in enumerate(results, 1):
                print(f"{i}. [{item['score']}] {item['title']}")
                print(f"   {item['file_path']}\n")
    except requests.exceptions.ConnectionError:
        print_error("Service is not running")
        return 1
    return 0


def cmd_record(args):
    """交互式创建记录"""
    print_info("Creating new record (Ctrl+C to cancel)\n")

    try:
        title = input("Title: ")
        rtype = input("Type (lesson/failure/decision): ") or "lesson"
        conclusion = input("Conclusion (one line): ")
        background = input("Background: ")

        print("\nKey points (one per line, empty line to finish):")
        points = []
        while True:
            p = input("  - ")
            if not p:
                break
            points.append(p)

        payload = {
            "task": args.task or "default",
            "title": title,
            "type": rtype,
            "conclusion": conclusion,
            "background": background,
            "key_points": points
        }

        r = requests.post(f"{API_BASE}/api/record", json=payload)
        if r.status_code == 200:
            print_success(f"Record saved: {r.json()['result']['file_path']}")
        else:
            print_error(f"Failed: {r.text}")
    except KeyboardInterrupt:
        print("\nCancelled")
    except requests.exceptions.ConnectionError:
        print_error("Service is not running")
        return 1
    return 0


def cmd_export(args):
    """导出备份"""
    import tarfile
    from datetime import datetime

    output = args.output or f"rag-backup-{datetime.now().strftime('%Y%m%d')}.tar.gz"

    print_info(f"Exporting to {output}...")
    with tarfile.open(output, "w:gz") as tar:
        tar.add(".rag/", arcname=".rag/")
        tar.add("ai-docs/", arcname="ai-docs/")
    print_success(f"Exported to {output}")
    return 0


def cmd_import(args):
    """导入备份"""
    import tarfile

    print_info(f"Importing from {args.input}...")
    with tarfile.open(args.input, "r:gz") as tar:
        tar.extractall()
    print_success("Imported")
    print_info("Run: ragctl index")
    return 0


def cmd_logs(args):
    """查看日志"""
    import os

    log_file = Path(".rag/rag.log")
    if not log_file.exists():
        print_info("No log file found")
        return 1

    if args.follow:
        # tail -f
        subprocess.run(["tail", "-f", str(log_file)])
    else:
        # cat last N lines
        n = args.n or 50
        subprocess.run(["tail", "-n", str(n), str(log_file)])
    return 0


def cmd_start(args):
    """启动服务"""
    print_info("Starting RAG service...")
    subprocess.run(["python", "-m", "rag.server"])
    return 0


def cmd_stop(args):
    """停止服务"""
    try:
        r = requests.get(f"{API_BASE}/health")
        print_error("Stop not implemented yet")
        print_info("Use: pkill -f 'rag.server'")
        return 1
    except requests.exceptions.ConnectionError:
        print_error("Service is not running")
        return 1


def main():
    parser = argparse.ArgumentParser(description="RAG CLI Tool")
    subparsers = parser.add_subparsers(dest="command")

    # status
    subparsers.add_parser("status", help="Show service status")

    # index
    subparsers.add_parser("index", help="Rebuild index")

    # search
    p_search = subparsers.add_parser("search", help="Search documents")
    p_search.add_argument("query", help="Search query")
    p_search.add_argument("-n", type=int, default=5, help="Number of results")

    # record
    p_record = subparsers.add_parser("record", help="Create new record")
    p_record.add_argument("-t", "--task", help="Task name")

    # export
    p_export = subparsers.add_parser("export", help="Export backup")
    p_export.add_argument("output", nargs="?", help="Output file")

    # import
    p_import = subparsers.add_parser("import", help="Import backup")
    p_import.add_argument("input", help="Input file")

    # logs
    p_logs = subparsers.add_parser("logs", help="View logs")
    p_logs.add_argument("-f", "--follow", action="store_true", help="Follow logs")
    p_logs.add_argument("-n", type=int, help="Number of lines")

    # start/stop
    subparsers.add_parser("start", help="Start service")
    subparsers.add_parser("stop", help="Stop service")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 1

    commands = {
        "status": cmd_status,
        "index": cmd_index,
        "search": cmd_search,
        "record": cmd_record,
        "export": cmd_export,
        "import": cmd_import,
        "logs": cmd_logs,
        "start": cmd_start,
        "stop": cmd_stop,
    }

    return commands[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
