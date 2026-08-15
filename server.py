#!/usr/bin/env python3
"""Local web console for configuring, running, and tracking the pipeline."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import mimetypes
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import uuid
from datetime import datetime
from email.utils import formatdate
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlparse

import pipeline


ROOT = Path(__file__).resolve().parent
STATIC_DIR = ROOT / "web"
RUN_ID_RE = re.compile(r"^[0-9]{8}-[0-9]{6}-[0-9a-f]{8}$")
MAX_LOG_BYTES = 500_000
MAX_RESULT_FILES = 10_000


def _now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def _read_json(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        return {}
    return data if isinstance(data, dict) else {}


def _write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    temp.replace(path)


def _tail(path: Path, limit: int = MAX_LOG_BYTES) -> str:
    if not path.exists():
        return ""
    with path.open("rb") as fp:
        fp.seek(0, os.SEEK_END)
        size = fp.tell()
        fp.seek(max(0, size - limit))
        return fp.read().decode("utf-8", errors="replace")


def _parse_time(value: str | None) -> datetime | None:
    if not value:
        return None
    try:
        return datetime.fromisoformat(value)
    except ValueError:
        return None


def _resolved(base: Path, value: str) -> str:
    return pipeline._resolved(base, value)


def _absolute_config(config: dict[str, Any], config_path: Path) -> dict[str, Any]:
    """Make path fields absolute so a run snapshot is location independent."""
    result = copy.deepcopy(config)
    base = config_path.resolve().parent
    result["dataset"]["input"] = _resolved(base, result["dataset"]["input"])
    result["preprocess"]["output_dir"] = _resolved(base, result["preprocess"]["output_dir"])
    result["aggregation"]["input_glob"] = _resolved(base, result["aggregation"]["input_glob"])
    result["aggregation"]["output_dir"] = _resolved(base, result["aggregation"]["output_dir"])
    result["discrete"]["merge_output"] = _resolved(base, result["discrete"]["merge_output"])
    result["discrete"]["discretize_output_dir"] = _resolved(
        base, result["discrete"]["discretize_output_dir"]
    )
    result["discrete"]["prepare_output_dir"] = _resolved(
        base, result["discrete"]["prepare_output_dir"]
    )
    return result


class RunManager:
    """Own one active process and a persistent directory of run records."""

    def __init__(self, config_path: Path, runs_dir: Path):
        self.config_path = config_path.resolve()
        self.runs_dir = runs_dir.resolve()
        self.runs_dir.mkdir(parents=True, exist_ok=True)
        self.lock = threading.RLock()
        self.process: subprocess.Popen[bytes] | None = None
        self.active_id: str | None = None
        self.stop_requested = False
        self._recover_interrupted_runs()

    def _run_dir(self, run_id: str) -> Path:
        if not RUN_ID_RE.fullmatch(run_id):
            raise ValueError("run id が不正です。")
        path = (self.runs_dir / run_id).resolve()
        if path.parent != self.runs_dir:
            raise ValueError("run id が不正です。")
        return path

    def _metadata_path(self, run_id: str) -> Path:
        return self._run_dir(run_id) / "metadata.json"

    def _recover_interrupted_runs(self) -> None:
        for metadata_path in self.runs_dir.glob("*/metadata.json"):
            metadata = _read_json(metadata_path)
            if metadata.get("status") in {"running", "stopping"}:
                metadata["status"] = "interrupted"
                metadata["finished_at"] = _now()
                metadata["pid"] = None
                metadata["message"] = "server restart detected"
                _write_json(metadata_path, metadata)

    def _new_id(self) -> str:
        stamp = datetime.now().astimezone().strftime("%Y%m%d-%H%M%S")
        return f"{stamp}-{uuid.uuid4().hex[:8]}"

    def start(self) -> dict[str, Any]:
        with self.lock:
            if self.process is not None and self.process.poll() is None:
                raise RuntimeError("パイプラインは既に実行中です。")

            config = pipeline.load_config(self.config_path)
            plan = pipeline.command_plan(self.config_path)
            run_id = self._new_id()
            run_dir = self._run_dir(run_id)
            run_dir.mkdir(parents=True)
            snapshot_path = run_dir / "config.yaml"
            runtime_config = _absolute_config(config, self.config_path)
            assert pipeline.yaml is not None
            snapshot_path.write_text(
                pipeline.yaml.safe_dump(runtime_config, allow_unicode=True, sort_keys=False),
                encoding="utf-8",
            )
            metadata = {
                "id": run_id,
                "status": "starting",
                "created_at": _now(),
                "started_at": None,
                "finished_at": None,
                "pid": None,
                "returncode": None,
                "stop_requested": False,
                "config_path": str(self.config_path),
                "config_snapshot": str(snapshot_path),
                "dataset": config["dataset"],
                "aggregation": {
                    "method": config["aggregation"]["method"],
                    "level": config["aggregation"]["level"],
                },
                "stages": [
                    {"name": item["stage"], "status": "pending", "index": index + 1}
                    for index, item in enumerate(plan)
                ],
                "total_stages": len(plan),
                "completed_stages": 0,
                "progress": 0,
                "current_stage": None,
                "result_count": 0,
                "results_truncated": False,
            }
            _write_json(run_dir / "metadata.json", metadata)

            log = (run_dir / "pipeline.log").open("wb")
            try:
                process = subprocess.Popen(
                    [
                        sys.executable,
                        "-u",
                        str(ROOT / "pipeline.py"),
                        "--config",
                        str(snapshot_path),
                        "--events",
                        str(run_dir / "events.jsonl"),
                    ],
                    cwd=ROOT,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                )
            except Exception:
                metadata["status"] = "failed"
                metadata["finished_at"] = _now()
                metadata["message"] = "failed to start process"
                _write_json(run_dir / "metadata.json", metadata)
                raise
            finally:
                log.close()

            self.process = process
            self.active_id = run_id
            self.stop_requested = False
            metadata["status"] = "running"
            metadata["started_at"] = _now()
            metadata["pid"] = process.pid
            _write_json(run_dir / "metadata.json", metadata)
            threading.Thread(target=self._wait, args=(run_id, process), daemon=True).start()
            return self.detail(run_id)

    def _wait(self, run_id: str, process: subprocess.Popen[bytes]) -> None:
        returncode = process.wait()
        with self.lock:
            metadata_path = self._metadata_path(run_id)
            metadata = self._with_events(run_id, _read_json(metadata_path))
            stopped = bool(metadata.get("stop_requested") or self.stop_requested)
            metadata["returncode"] = returncode
            metadata["pid"] = None
            metadata["finished_at"] = _now()
            metadata["status"] = "stopped" if stopped else ("succeeded" if returncode == 0 else "failed")
            metadata["current_stage"] = None
            try:
                manifest, truncated = self._scan_results(run_id)
                metadata["result_count"] = len(manifest)
                metadata["results_truncated"] = truncated
                _write_json(
                    self._run_dir(run_id) / "results.json",
                    {"files": manifest, "truncated": truncated},
                )
            except Exception as exc:
                metadata["results_error"] = str(exc)
            _write_json(metadata_path, metadata)
            if self.process is process:
                self.process = None
                self.stop_requested = False

    def stop(self) -> None:
        with self.lock:
            process = self.process
            if process is None or process.poll() is not None or self.active_id is None:
                raise RuntimeError("実行中のパイプラインはありません。")
            self.stop_requested = True
            metadata_path = self._metadata_path(self.active_id)
            metadata = _read_json(metadata_path)
            metadata["status"] = "stopping"
            metadata["stop_requested"] = True
            metadata["message"] = "termination requested"
            _write_json(metadata_path, metadata)
            os.killpg(process.pid, signal.SIGTERM)

    def _events(self, run_id: str) -> list[dict[str, Any]]:
        path = self._run_dir(run_id) / "events.jsonl"
        if not path.exists():
            return []
        events: list[dict[str, Any]] = []
        with path.open(encoding="utf-8", errors="replace") as fp:
            for line in fp:
                try:
                    item = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if isinstance(item, dict):
                    events.append(item)
        return events

    def _with_events(self, run_id: str, metadata: dict[str, Any]) -> dict[str, Any]:
        result = copy.deepcopy(metadata)
        stages = {item["name"]: item for item in result.get("stages", [])}
        for event in self._events(run_id):
            name = event.get("stage")
            stage = stages.get(name)
            kind = event.get("event")
            if kind == "stage_started" and stage is not None:
                stage["status"] = "running"
                stage["started_at"] = event.get("time")
                stage["command"] = event.get("command", [])
                result["current_stage"] = name
            elif kind == "stage_completed" and stage is not None:
                stage["status"] = "succeeded"
                stage["finished_at"] = event.get("time")
                stage["elapsed_seconds"] = event.get("elapsed_seconds")
                if result.get("current_stage") == name:
                    result["current_stage"] = None
            elif kind == "stage_failed" and stage is not None:
                stage["status"] = "failed"
                stage["finished_at"] = event.get("time")
                stage["elapsed_seconds"] = event.get("elapsed_seconds")
                stage["returncode"] = event.get("returncode")
            elif kind == "pipeline_completed":
                result["current_stage"] = None
        completed = sum(item.get("status") == "succeeded" for item in stages.values())
        total = int(result.get("total_stages") or len(stages))
        result["completed_stages"] = completed
        result["progress"] = round(completed * 100 / total) if total else 0
        if result.get("status") == "succeeded":
            result["progress"] = 100
        started = _parse_time(result.get("started_at"))
        finished = _parse_time(result.get("finished_at"))
        if started:
            end = finished or datetime.now().astimezone()
            result["elapsed_seconds"] = max(0, round((end - started).total_seconds()))
        return result

    def _summary(self, metadata: dict[str, Any]) -> dict[str, Any]:
        keys = {
            "id", "status", "created_at", "started_at", "finished_at", "pid", "returncode",
            "dataset", "aggregation", "stages", "total_stages", "completed_stages", "progress",
            "current_stage", "elapsed_seconds", "result_count", "results_truncated", "message",
        }
        return {key: value for key, value in metadata.items() if key in keys}

    def list_runs(self) -> list[dict[str, Any]]:
        items: list[dict[str, Any]] = []
        for metadata_path in self.runs_dir.glob("*/metadata.json"):
            metadata = _read_json(metadata_path)
            run_id = str(metadata.get("id", ""))
            if RUN_ID_RE.fullmatch(run_id):
                items.append(self._summary(self._with_events(run_id, metadata)))
        return sorted(items, key=lambda item: item.get("created_at") or "", reverse=True)

    def detail(
        self, run_id: str, include_log: bool = True, include_files: bool = True
    ) -> dict[str, Any]:
        metadata_path = self._metadata_path(run_id)
        metadata = _read_json(metadata_path)
        if not metadata:
            raise FileNotFoundError("run が見つかりません。")
        result = self._with_events(run_id, metadata)
        result["running"] = result.get("status") in {"starting", "running", "stopping"}
        if include_log:
            result["log"] = _tail(self._run_dir(run_id) / "pipeline.log")
        if include_files:
            manifest = _read_json(self._run_dir(run_id) / "results.json")
            result["files"] = [
                {key: value for key, value in item.items() if key != "absolute_path"}
                for item in manifest.get("files", [])
                if isinstance(item, dict)
            ]
        return result

    def status(self) -> dict[str, Any]:
        with self.lock:
            run_id = self.active_id
            if run_id is None:
                runs = self.list_runs()
                run_id = runs[0]["id"] if runs else None
            if run_id is None:
                return {
                    "running": False, "returncode": None, "pid": None, "log": "",
                    "run": None,
                }
            detail = self.detail(run_id, include_files=False)
            return {
                "running": detail["running"],
                "returncode": detail.get("returncode"),
                "pid": detail.get("pid"),
                "log": detail.get("log", ""),
                "run": detail,
            }

    def _result_roots(self, run_id: str) -> list[tuple[str, Path]]:
        config = pipeline.load_config(self._run_dir(run_id) / "config.yaml")
        roots: list[tuple[str, Path]] = []
        if config["pipeline"]["preprocess"]:
            roots.append(("preprocess", Path(config["preprocess"]["output_dir"])))
        if config["pipeline"]["aggregation"]:
            roots.append(("aggregation", Path(config["aggregation"]["output_dir"])))
        if config["pipeline"]["discrete"]:
            roots.extend([
                ("merged", Path(config["discrete"]["merge_output"])),
                ("discretized", Path(config["discrete"]["discretize_output_dir"])),
                ("prepared", Path(config["discrete"]["prepare_output_dir"])),
            ])
        return roots

    def _scan_results(self, run_id: str) -> tuple[list[dict[str, Any]], bool]:
        metadata = _read_json(self._metadata_path(run_id))
        started = _parse_time(metadata.get("started_at"))
        started_timestamp = started.timestamp() if started else 0
        files: list[dict[str, Any]] = []
        seen: set[Path] = set()
        truncated = False
        for group, root in self._result_roots(run_id):
            candidates = [root] if root.is_file() else (root.rglob("*") if root.is_dir() else [])
            for candidate in candidates:
                if not candidate.is_file():
                    continue
                resolved = candidate.resolve()
                if resolved in seen:
                    continue
                seen.add(resolved)
                stat = resolved.stat()
                try:
                    display = str(resolved.relative_to(ROOT))
                except ValueError:
                    display = str(resolved)
                file_id = hashlib.sha256(str(resolved).encode()).hexdigest()[:16]
                files.append({
                    "id": file_id,
                    "path": display,
                    "name": resolved.name,
                    "group": group,
                    "size": stat.st_size,
                    "modified_at": datetime.fromtimestamp(
                        stat.st_mtime, tz=datetime.now().astimezone().tzinfo
                    ).isoformat(timespec="seconds"),
                    "updated_in_run": stat.st_mtime >= started_timestamp - 1,
                    "absolute_path": str(resolved),
                })
                if len(files) >= MAX_RESULT_FILES:
                    truncated = True
                    return files, truncated
        files.sort(key=lambda item: (item["group"], item["path"]))
        return files, truncated

    def refresh_results(self, run_id: str) -> dict[str, Any]:
        if self.detail(run_id, include_log=False)["running"]:
            raise RuntimeError("実行中は結果一覧を更新できません。")
        files, truncated = self._scan_results(run_id)
        _write_json(self._run_dir(run_id) / "results.json", {"files": files, "truncated": truncated})
        metadata_path = self._metadata_path(run_id)
        metadata = _read_json(metadata_path)
        metadata["result_count"] = len(files)
        metadata["results_truncated"] = truncated
        _write_json(metadata_path, metadata)
        return self.detail(run_id)

    def result_file(self, run_id: str, file_id: str) -> Path:
        manifest = _read_json(self._run_dir(run_id) / "results.json")
        for item in manifest.get("files", []):
            if isinstance(item, dict) and item.get("id") == file_id:
                path = Path(str(item.get("absolute_path", ""))).resolve()
                if path.is_file():
                    return path
                raise FileNotFoundError("結果ファイルがありません。")
        raise FileNotFoundError("結果ファイルが見つかりません。")

    def delete(self, run_id: str) -> None:
        with self.lock:
            if self.active_id == run_id and self.process is not None and self.process.poll() is None:
                raise RuntimeError("実行中の履歴は削除できません。")
            run_dir = self._run_dir(run_id)
            if not (run_dir / "metadata.json").exists():
                raise FileNotFoundError("run が見つかりません。")
            shutil.rmtree(run_dir)
            if self.active_id == run_id:
                self.active_id = None


class Handler(BaseHTTPRequestHandler):
    server_version = "ScanpyBNPreprocessUI/2.0"

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

    def _run_route(self, path: str) -> tuple[str | None, str | None, str | None]:
        parts = [unquote(part) for part in path.strip("/").split("/")]
        if len(parts) >= 3 and parts[:2] == ["api", "runs"]:
            run_id = parts[2]
            action = parts[3] if len(parts) >= 4 else None
            value = parts[4] if len(parts) >= 5 else None
            return run_id, action, value
        return None, None, None

    def do_GET(self) -> None:  # noqa: N802
        path = urlparse(self.path).path
        try:
            if path == "/api/config":
                config = pipeline.load_config(self.app.config_path)
                self._json({"config": config, "path": str(self.app.config_path)})
                return
            if path == "/api/status":
                self._json(self.app.runs.status())
                return
            if path == "/api/plan":
                self._json({"plan": pipeline.command_plan(self.app.config_path)})
                return
            if path == "/api/runs":
                self._json({"runs": self.app.runs.list_runs()})
                return
            run_id, action, value = self._run_route(path)
            if run_id and action is None:
                self._json({"run": self.app.runs.detail(run_id)})
                return
            if run_id and action == "files" and value:
                self._serve_result(self.app.runs.result_file(run_id, value))
                return
            if run_id and action == "log" and value is None:
                self._serve_result(self.app.runs._run_dir(run_id) / "pipeline.log")
                return
            if path in {"/", "/index.html"}:
                self._serve(STATIC_DIR / "index.html", "text/html; charset=utf-8")
                return
            self.send_error(HTTPStatus.NOT_FOUND)
        except pipeline.ConfigError as exc:
            self._json({"error": str(exc)}, HTTPStatus.BAD_REQUEST)
        except (ValueError, FileNotFoundError) as exc:
            self._json({"error": str(exc)}, HTTPStatus.NOT_FOUND)

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
                run = self.app.runs.start()
                self._json({"ok": True, "run": run}, HTTPStatus.ACCEPTED)
                return
            if path == "/api/stop":
                self.app.runs.stop()
                self._json({"ok": True})
                return
            run_id, action, _ = self._run_route(path)
            if run_id and action == "refresh-results":
                self._json({"ok": True, "run": self.app.runs.refresh_results(run_id)})
                return
            self.send_error(HTTPStatus.NOT_FOUND)
        except (RuntimeError, pipeline.ConfigError) as exc:
            self._json({"error": str(exc)}, HTTPStatus.CONFLICT)
        except (ValueError, FileNotFoundError) as exc:
            self._json({"error": str(exc)}, HTTPStatus.NOT_FOUND)

    def do_DELETE(self) -> None:  # noqa: N802
        path = urlparse(self.path).path
        run_id, action, _ = self._run_route(path)
        if not run_id or action is not None:
            self.send_error(HTTPStatus.NOT_FOUND)
            return
        try:
            self.app.runs.delete(run_id)
            self._json({"ok": True})
        except RuntimeError as exc:
            self._json({"error": str(exc)}, HTTPStatus.CONFLICT)
        except (ValueError, FileNotFoundError) as exc:
            self._json({"error": str(exc)}, HTTPStatus.NOT_FOUND)

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

    def _serve_result(self, path: Path) -> None:
        stat = path.stat()
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", mimetypes.guess_type(path.name)[0] or "application/octet-stream")
        self.send_header("Content-Length", str(stat.st_size))
        self.send_header("Last-Modified", formatdate(stat.st_mtime, usegmt=True))
        self.send_header("Content-Disposition", f'attachment; filename="{path.name.replace(chr(34), "")}"')
        self.end_headers()
        with path.open("rb") as fp:
            shutil.copyfileobj(fp, self.wfile)

    def log_message(self, fmt: str, *args: Any) -> None:
        sys.stderr.write("[%s] %s\n" % (self.log_date_time_string(), fmt % args))


class AppServer(ThreadingHTTPServer):
    def __init__(self, address: tuple[str, int], config_path: Path, runs_dir: Path):
        super().__init__(address, Handler)
        self.config_path = config_path.resolve()
        self.runs = RunManager(self.config_path, runs_dir)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1", help="既定はローカル接続のみ")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--config", type=Path, default=ROOT / "config.yaml")
    parser.add_argument("--runs-dir", type=Path, default=ROOT / ".pipeline-runs")
    parser.add_argument("--log", type=Path, default=None, help=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        pipeline.load_config(args.config)
    except pipeline.ConfigError as exc:
        print(f"config error: {exc}", file=sys.stderr)
        return 2
    server = AppServer((args.host, args.port), args.config, args.runs_dir)
    print(f"ScanpyBNPreprocess UI: http://{args.host}:{args.port}")
    print(f"config: {args.config.resolve()}")
    print(f"run history: {args.runs_dir.resolve()}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nserver stopped")
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
