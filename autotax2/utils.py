"""General utilities."""
from __future__ import annotations

import hashlib
import logging
import shutil
import subprocess
import time
from pathlib import Path
from typing import Iterable, Sequence

logger = logging.getLogger("autotax2.utils")


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def require_file(path: str | Path, label: str = "file") -> Path:
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        raise FileNotFoundError(f"Required {label} not found or empty: {p}")
    return p


def require_executable(name: str) -> str:
    exe = shutil.which(name)
    if not exe:
        raise FileNotFoundError(f"Executable not found in PATH: {name}")
    return exe


def run_command(cmd: Sequence[str], log_path: str | Path | None = None, dry_run: bool = False) -> None:
    """Run an external command with detailed logging."""
    cmd_str = " ".join(map(str, cmd))
    logger.info("Running command: %s", cmd_str)
    if dry_run:
        logger.info("Dry run enabled; command not executed.")
        return

    start = time.time()
    stdout_target = subprocess.PIPE
    stderr_target = subprocess.STDOUT
    with subprocess.Popen(
        list(map(str, cmd)), stdout=stdout_target, stderr=stderr_target, text=True, bufsize=1
    ) as proc:
        log_fh = open(log_path, "a", encoding="utf-8") if log_path else None
        try:
            assert proc.stdout is not None
            for line in proc.stdout:
                line = line.rstrip("\n")
                logger.debug("[external] %s", line)
                if log_fh:
                    log_fh.write(line + "\n")
            ret = proc.wait()
        finally:
            if log_fh:
                log_fh.close()
    elapsed = time.time() - start
    if ret != 0:
        raise RuntimeError(f"Command failed with exit code {ret}: {cmd_str}")
    logger.info("Command completed in %.2f seconds", elapsed)


def md5_text(text: str) -> str:
    return hashlib.md5(text.encode("utf-8")).hexdigest()


def md5_file(path: str | Path, block_size: int = 1024 * 1024) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(block_size), b""):
            h.update(block)
    return h.hexdigest()


def strip_gaps(seq: str) -> str:
    return seq.replace("-", "").replace(".", "")


def chunked(iterable: Iterable, size: int):
    chunk = []
    for item in iterable:
        chunk.append(item)
        if len(chunk) == size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk
