"""Configuration loading for autotax2."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class Autotax2Config:
    """Lightweight configuration wrapper."""

    path: Path
    data: dict[str, Any] = field(default_factory=dict)


def load_config(path: str | Path) -> Autotax2Config:
    """Load a YAML configuration file."""
    config_path = Path(path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with config_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}

    if not isinstance(data, dict):
        raise ValueError("Configuration root must be a mapping.")

    return Autotax2Config(path=config_path, data=data)
