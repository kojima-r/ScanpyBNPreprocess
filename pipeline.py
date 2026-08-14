#!/usr/bin/env python3
"""Run the preprocessing pipeline from a YAML configuration file."""

from __future__ import annotations

import argparse
import copy
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:  # pragma: no cover - exercised by installation state
    yaml = None


ROOT = Path(__file__).resolve().parent
MODES = {"bbknn", "facs", "droplet"}
LEVELS = {"tissue", "age", "batch"}
METHODS = {"resample", "pseudo_bulk"}

DEFAULTS: dict[str, Any] = {
    "dataset": {"mode": "bbknn", "target": "bbknn", "input": ""},
    "pipeline": {
        "preprocess": True,
        "aggregation": True,
        "continuous": True,
        "discrete": True,
    },
    "preprocess": {"output_dir": ""},
    "aggregation": {
        "method": "resample",
        "level": "batch",
        "workers": 4,
        "input_glob": "",
        "output_dir": "",
        "resample": {"count": 10, "size": "same"},
    },
    "continuous": {
        "transpose_output_dir": "",
        "merge_output": "",
        "min_variance": 0.0,
    },
    "discrete": {
        "merge_output": "",
        "min_variance": 0.0,
        "discretize_output_dir": "",
        "min_entropy": 0.0,
        "prepare_output_dir": "",
    },
}


class ConfigError(ValueError):
    pass


def _deep_merge(base: dict[str, Any], supplied: dict[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(base)
    for key, value in supplied.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def load_config(path: Path) -> dict[str, Any]:
    if yaml is None:
        raise ConfigError("PyYAML がありません。`python -m pip install -r requirements.txt` を実行してください。")
    try:
        raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ConfigError(f"設定ファイルがありません: {path}") from exc
    except yaml.YAMLError as exc:
        raise ConfigError(f"YAML の構文エラー: {exc}") from exc
    if not isinstance(raw, dict):
        raise ConfigError("config のルートはマッピングである必要があります。")
    return normalize_config(_deep_merge(DEFAULTS, raw))


def normalize_config(config: dict[str, Any]) -> dict[str, Any]:
    errors: list[str] = []
    dataset = config.get("dataset", {})
    flags = config.get("pipeline", {})
    aggregation = config.get("aggregation", {})
    continuous = config.get("continuous", {})
    discrete = config.get("discrete", {})

    mode = dataset.get("mode")
    target = str(dataset.get("target", "")).strip()
    method = aggregation.get("method")
    level = aggregation.get("level")
    if mode not in MODES:
        errors.append(f"dataset.mode は {sorted(MODES)} のいずれかです。")
    if not target or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in target):
        errors.append("dataset.target は英数字、_、- のみ使用できます。")
    if method not in METHODS:
        errors.append(f"aggregation.method は {sorted(METHODS)} のいずれかです。")
    if level not in LEVELS:
        errors.append(f"aggregation.level は {sorted(LEVELS)} のいずれかです。")
    for key in ("preprocess", "aggregation", "continuous", "discrete"):
        if not isinstance(flags.get(key), bool):
            errors.append(f"pipeline.{key} は true または false にしてください。")
    try:
        workers = int(aggregation.get("workers"))
        if workers < 1:
            raise ValueError
        aggregation["workers"] = workers
    except (TypeError, ValueError):
        errors.append("aggregation.workers は 1 以上の整数です。")
    try:
        count = int(aggregation.get("resample", {}).get("count"))
        if count < 1:
            raise ValueError
        aggregation["resample"]["count"] = count
    except (TypeError, ValueError):
        errors.append("aggregation.resample.count は 1 以上の整数です。")
    if aggregation.get("resample", {}).get("size") not in {"same", "root"}:
        errors.append("aggregation.resample.size は same または root です。")
    for section, key in ((continuous, "min_variance"), (discrete, "min_variance"), (discrete, "min_entropy")):
        try:
            section[key] = float(section.get(key))
            if section[key] < 0:
                raise ValueError
        except (TypeError, ValueError):
            errors.append(f"{key} は 0 以上の数値です。")
    if flags.get("continuous") or flags.get("discrete"):
        if not flags.get("aggregation"):
            errors.append("continuous/discrete を実行するには pipeline.aggregation を true にしてください。")
    if errors:
        raise ConfigError("\n".join(errors))

    source = "r" if method == "resample" else "p"
    dataset["input"] = dataset.get("input") or f"data/tabula-muris-senis-{mode}-processed-official-annotations.h5ad"
    config["preprocess"]["output_dir"] = config["preprocess"].get("output_dir") or f"data01_{target}"
    aggregation["input_glob"] = aggregation.get("input_glob") or f"data01_{target}/*.txt"
    aggregation["output_dir"] = aggregation.get("output_dir") or f"data02_{target}_{source}_{level}"
    continuous["transpose_output_dir"] = continuous.get("transpose_output_dir") or f"data02_{target}_{source}_{level}_t"
    continuous["merge_output"] = continuous.get("merge_output") or f"data03_{target}_{source}_{level}_t/all.txt"
    discrete["merge_output"] = discrete.get("merge_output") or f"data03_{target}_{source}_{level}/all.txt"
    discrete["discretize_output_dir"] = discrete.get("discretize_output_dir") or str(Path(discrete["merge_output"]).parent)
    discrete["prepare_output_dir"] = discrete.get("prepare_output_dir") or f"data04_{target}_{source}_{level}_disc"
    return config


def _resolved(config_dir: Path, value: str) -> str:
    path = Path(value).expanduser()
    return str(path if path.is_absolute() else (config_dir / path).resolve())


def build_commands(config: dict[str, Any], config_path: Path) -> list[tuple[str, list[str]]]:
    base = config_path.resolve().parent
    py = sys.executable
    dataset = config["dataset"]
    flags = config["pipeline"]
    agg = config["aggregation"]
    continuous = config["continuous"]
    discrete = config["discrete"]
    target, level = dataset["target"], agg["level"]
    source = "r" if agg["method"] == "resample" else "p"
    commands: list[tuple[str, list[str]]] = []

    if flags["preprocess"]:
        commands.append(("01 preprocess", [py, str(ROOT / "01preprocess.py"), "--mode", dataset["mode"],
            "--input", _resolved(base, dataset["input"]), "--out-dir", _resolved(base, config["preprocess"]["output_dir"])]))
    if flags["aggregation"]:
        script = "02resample.py" if source == "r" else "02pseudo_bulk.py"
        cmd = [py, str(ROOT / script), "--target", target, "--input-glob", _resolved(base, agg["input_glob"]),
               "--level", level, "--out-dir", _resolved(base, agg["output_dir"]), "--workers", str(agg["workers"])]
        if source == "r":
            cmd += ["--n-resamples", str(agg["resample"]["count"]), "--size", agg["resample"]["size"]]
        commands.append(("02 aggregation", cmd))
    if flags["continuous"]:
        commands.append(("03 transpose", [py, str(ROOT / "03transpose.py"), "--target", target, "--source", source,
            "--level", level, "--input-glob", _resolved(base, f"{agg['output_dir']}/*.txt"),
            "--out-dir", _resolved(base, continuous["transpose_output_dir"])]))
        commands.append(("04 merge continuous", [py, str(ROOT / "04merge.py"), "--target", target, "--source", source,
            "--level", level, "--transposed", "--input-glob", _resolved(base, f"{continuous['transpose_output_dir']}/*.txt"),
            "--out", _resolved(base, continuous["merge_output"]), "--min-var", str(continuous["min_variance"])]))
    if flags["discrete"]:
        commands.append(("04 merge discrete", [py, str(ROOT / "04merge.py"), "--target", target, "--source", source,
            "--level", level, "--input-glob", _resolved(base, f"{agg['output_dir']}/*.txt"),
            "--out", _resolved(base, discrete["merge_output"]), "--min-var", str(discrete["min_variance"])]))
        commands.append(("05 discretize", [py, str(ROOT / "05disc.py"), "--target", target, "--source", source,
            "--level", level, "--input", _resolved(base, discrete["merge_output"]),
            "--out-dir", _resolved(base, discrete["discretize_output_dir"]), "--min-entropy", str(discrete["min_entropy"])]))
        disc_dir = Path(discrete["discretize_output_dir"])
        commands.append(("06 prepare discrete", [py, str(ROOT / "06prep_disc.py"), "--target", target, "--source", source,
            "--level", level, "--bin-input", _resolved(base, str(disc_dir / "all_disc.txt")),
            "--tri-input", _resolved(base, str(disc_dir / "all_disc_tri.txt")),
            "--out-dir", _resolved(base, discrete["prepare_output_dir"])]))
    if not commands:
        raise ConfigError("実行対象がありません。pipeline の項目を1つ以上 true にしてください。")
    return commands


def command_plan(config_path: Path) -> list[dict[str, Any]]:
    config = load_config(config_path)
    return [{"stage": stage, "command": command} for stage, command in build_commands(config, config_path)]


def run(config_path: Path, dry_run: bool = False) -> int:
    config = load_config(config_path)
    commands = build_commands(config, config_path)
    print(f"config: {config_path.resolve()}", flush=True)
    for stage, command in commands:
        print(f"\n=== {stage} ===", flush=True)
        print("command: " + json.dumps(command, ensure_ascii=False), flush=True)
        if dry_run:
            continue
        started = time.monotonic()
        completed = subprocess.run(command, cwd=ROOT)
        if completed.returncode:
            print(f"ERROR: {stage} failed (exit={completed.returncode})", file=sys.stderr, flush=True)
            return completed.returncode
        print(f"completed in {time.monotonic() - started:.1f}s", flush=True)
    print("\nPipeline completed." if not dry_run else "\nDry run completed.", flush=True)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config.yaml")
    parser.add_argument("--dry-run", action="store_true", help="設定を検証し、実行コマンドだけ表示")
    parser.add_argument("--print-plan", action="store_true", help="実行計画を JSON で表示")
    args = parser.parse_args()
    try:
        if args.print_plan:
            print(json.dumps(command_plan(args.config), ensure_ascii=False, indent=2))
            return 0
        return run(args.config, args.dry_run)
    except ConfigError as exc:
        print(f"config error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

