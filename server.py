#!/usr/bin/env python3
"""Local web UI server for editing config.yaml and running the pipeline."""

from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
import threading
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import pipeline


ROOT = Path(__file__).resolve().parent
STATIC_DIR = ROOT / "web"


class Job:
    def __init__(self, config_path: Path, log_path: Path):
        self.config_path = config_path
        self.log_path = log_path
        self.process: subprocess.Popen[bytes] | None = None
        self.returncode: int | None = None
        self.lock = threading.Lock()

    def start(self) -> None:
        with self.lock:
            if self.process is not None and self.process.poll() is None:
                raise RuntimeError("パイプラインは既に実行中です。")
            pipeline.load_config(self.config_path)
            self.log_path.parent.mkdir(parents=True, exist_ok=True)
            log = self.log_path.open("wb")
            try:
                self.process = subprocess.Popen(
                    [sys.executable, "-u", str(ROOT / "pipeline.py"), "--config", str(self.config_path)],
                    cwd=ROOT,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                )
            finally:
                log.close()
            self.returncode = None
            threading.Thread(target=self._wait, daemon=True).start()

    def _wait(self) -> None:
        process = self.process
        if process is None:
            return
        code = process.wait()
        with self.lock:
            self.returncode = code

    def stop(self) -> None:
        with self.lock:
            process = self.process
            if process is None or process.poll() is not None:
                raise RuntimeError("実行中のパイプラインはありません。")
            os.killpg(process.pid, signal.SIGTERM)

    def status(self) -> dict[str, Any]:
        with self.lock:
            process = self.process
            running = process is not None and process.poll() is None
            code = None if running else (self.returncode if process is not None else None)
            pid = process.pid if running else None
        log = ""
        if self.log_path.exists():
            with self.log_path.open("rb") as fp:
                fp.seek(0, os.SEEK_END)
                size = fp.tell()
                fp.seek(max(0, size - 200_000))
                log = fp.read().decode("utf-8", errors="replace")
        return {"running": running, "returncode": code, "pid": pid, "log": log}


class Handler(BaseHTTPRequestHandler):
    server_version = "ScanpyBNPreprocessUI/1.0"

    @property
    def app(self) -> "AppServer":
        return self.server  # type: ignore[return-value]

    def _json(self, data: Any, status: int = HTTPStatus.OK) -> None:
        body = json.dumps(data, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _body(self) -> dict[str, Any]:
        try:
            length = int(self.headers.get("Content-Length", "0"))
        except ValueError as exc:
            raise ValueError("Content-Length が不正です。") from exc
        if length > 1_000_000:
            raise ValueError("リクエストが大きすぎます。")
        data = json.loads(self.rfile.read(length) or b"{}")
        if not isinstance(data, dict):
            raise ValueError("JSON object が必要です。")
        return data

    def do_GET(self) -> None:  # noqa: N802
        path = urlparse(self.path).path
        if path == "/api/config":
            try:
                config = pipeline.load_config(self.app.config_path)
                self._json({"config": config, "path": str(self.app.config_path)})
            except pipeline.ConfigError as exc:
                self._json({"error": str(exc)}, HTTPStatus.BAD_REQUEST)
            return
        if path == "/api/status":
            self._json(self.app.job.status())
            return
        if path == "/api/plan":
            try:
                self._json({"plan": pipeline.command_plan(self.app.config_path)})
            except pipeline.ConfigError as exc:
                self._json({"error": str(exc)}, HTTPStatus.BAD_REQUEST)
            return
        if path in {"/", "/index.html"}:
            self._serve(STATIC_DIR / "index.html", "text/html; charset=utf-8")
            return
        self.send_error(HTTPStatus.NOT_FOUND)

    def do_PUT(self) -> None:  # noqa: N802
        if urlparse(self.path).path != "/api/config":
            self.send_error(HTTPStatus.NOT_FOUND)
            return
        try:
            body = self._body()
            config = body.get("config")
            if not isinstance(config, dict):
                raise ValueError("config object が必要です。")
            normalized = pipeline.normalize_config(pipeline._deep_merge(pipeline.DEFAULTS, config))
            if pipeline.yaml is None:
                raise pipeline.ConfigError("PyYAML がありません。")
            text = pipeline.yaml.safe_dump(normalized, allow_unicode=True, sort_keys=False)
            temp = self.app.config_path.with_suffix(self.app.config_path.suffix + ".tmp")
            temp.write_text(text, encoding="utf-8")
            temp.replace(self.app.config_path)
            self._json({"ok": True, "config": normalized})
        except (ValueError, pipeline.ConfigError) as exc:
            self._json({"error": str(exc)}, HTTPStatus.BAD_REQUEST)

    def do_POST(self) -> None:  # noqa: N802
        path = urlparse(self.path).path
        try:
            if path == "/api/run":
                self.app.job.start()
                self._json({"ok": True, **self.app.job.status()}, HTTPStatus.ACCEPTED)
                return
            if path == "/api/stop":
                self.app.job.stop()
                self._json({"ok": True})
                return
            self.send_error(HTTPStatus.NOT_FOUND)
        except (RuntimeError, pipeline.ConfigError) as exc:
            self._json({"error": str(exc)}, HTTPStatus.CONFLICT)

    def _serve(self, path: Path, content_type: str) -> None:
        try:
            body = path.read_bytes()
        except FileNotFoundError:
            self.send_error(HTTPStatus.NOT_FOUND)
            return
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, fmt: str, *args: Any) -> None:
        sys.stderr.write("[%s] %s\n" % (self.log_date_time_string(), fmt % args))


class AppServer(ThreadingHTTPServer):
    def __init__(self, address: tuple[str, int], config_path: Path, log_path: Path):
        super().__init__(address, Handler)
        self.config_path = config_path.resolve()
        self.job = Job(self.config_path, log_path.resolve())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1", help="既定はローカル接続のみ")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--config", type=Path, default=ROOT / "config.yaml")
    parser.add_argument("--log", type=Path, default=ROOT / ".pipeline-ui.log")
    args = parser.parse_args()
    try:
        pipeline.load_config(args.config)
    except pipeline.ConfigError as exc:
        print(f"config error: {exc}", file=sys.stderr)
        return 2
    server = AppServer((args.host, args.port), args.config, args.log)
    print(f"ScanpyBNPreprocess UI: http://{args.host}:{args.port}")
    print(f"config: {args.config.resolve()}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nserver stopped")
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

