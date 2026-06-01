"""Logging helpers for AutoTax2."""
from __future__ import annotations

import logging
from pathlib import Path
from rich.logging import RichHandler


def setup_logging(log_file: str | Path | None = None, debug: bool = False) -> logging.Logger:
    """Configure terminal and optional file logging.

    Parameters
    ----------
    log_file:
        Optional path to a log file. Parent directories are created.
    debug:
        If true, set DEBUG level; otherwise INFO.
    """
    level = logging.DEBUG if debug else logging.INFO
    handlers: list[logging.Handler] = [RichHandler(rich_tracebacks=True, show_time=True)]

    if log_file:
        path = Path(log_file)
        path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(path, mode="a", encoding="utf-8")
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")
        )
        handlers.append(file_handler)

    logging.basicConfig(
        level=level,
        format="%(message)s",
        handlers=handlers,
        force=True,
    )
    logger = logging.getLogger("autotax2")
    logger.debug("Logging initialized. debug=%s log_file=%s", debug, log_file)
    return logger
