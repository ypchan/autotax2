"""Dated audit log helpers for autotax2 commands."""

from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path


def write_event_log(
    build: str | Path,
    command: str,
    details: dict[str, object] | None = None,
) -> Path:
    """Write a small dated command audit log and return its path."""
    build_dir = Path(build)
    logs_dir = build_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().astimezone()
    safe_command = re.sub(r"[^A-Za-z0-9_.-]+", "_", command.strip()) or "command"
    stamp = timestamp.strftime("%Y%m%d%H%M%S")
    path = logs_dir / f"{safe_command}_date{stamp}.log"
    counter = 2
    while path.exists():
        path = logs_dir / f"{safe_command}_date{stamp}_{counter}.log"
        counter += 1

    rows = {
        "timestamp": timestamp.isoformat(timespec="seconds"),
        "command": command,
    }
    if details:
        rows.update({key: str(value) for key, value in details.items()})
    path.write_text("".join(f"{key}={value}\n" for key, value in rows.items()), encoding="utf-8")
    return path
