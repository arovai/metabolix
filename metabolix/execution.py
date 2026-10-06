"""Safe discovery and execution of FSL-MRS command-line tools."""

from __future__ import annotations

import json
import logging
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def find_executable(name: str, binary_dir: str | None = None) -> str | None:
    if binary_dir:
        candidate = Path(binary_dir) / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate.resolve())
    return shutil.which(name)


def dependency_report(config: dict[str, Any]) -> dict[str, str | None]:
    binary_dir = config["basis"].get("binary_dir")
    return {
        "fsl_mrs_proc": find_executable("fsl_mrs_proc", binary_dir),
        "fsl_mrsi": find_executable("fsl_mrsi", binary_dir),
    }


def executable_version(name: str, binary_dir: str | None = None) -> str:
    executable = find_executable(name, binary_dir)
    if not executable:
        return "unavailable"
    completed = subprocess.run([executable, "--version"], text=True, capture_output=True, check=False)
    value = (completed.stdout or completed.stderr).strip()
    return value.splitlines()[0] if value and completed.returncode == 0 else "unreported"


def run_command(
    argv: list[str],
    stage_dir: Path,
    label: str,
    logger: logging.Logger,
    dry_run: bool = False,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "label": label,
        "argv": argv,
        "started": datetime.now(timezone.utc).isoformat(),
        "dry_run": dry_run,
    }
    if dry_run:
        logger.info("DRY RUN: %s", argv)
        record.update({"exit_status": None, "stdout_file": None, "stderr_file": None})
        return record
    logger.info("Running %s", label)
    completed = subprocess.run(argv, text=True, capture_output=True, check=False)
    log_dir = stage_dir.parent / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = log_dir / f"{label}_stdout.txt"
    stderr_path = log_dir / f"{label}_stderr.txt"
    stdout_path.write_text(completed.stdout)
    stderr_path.write_text(completed.stderr)
    record.update(
        {
            "finished": datetime.now(timezone.utc).isoformat(),
            "exit_status": completed.returncode,
            "stdout_file": str(stdout_path),
            "stderr_file": str(stderr_path),
        }
    )
    if completed.stdout.strip():
        logger.debug("%s stdout:\n%s", label, completed.stdout.rstrip())
    if completed.stderr.strip():
        logger.warning("%s stderr:\n%s", label, completed.stderr.rstrip())
    combined_output = f"{completed.stdout}\n{completed.stderr}"
    if "not have enough samples" in combined_output.lower() or "100,000 samples" in combined_output:
        record["warnings"] = [line.strip() for line in combined_output.splitlines() if "sample" in line.lower() or "covariance" in line.lower()]
    return record


def save_manifest(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")