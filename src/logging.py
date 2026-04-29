from __future__ import annotations

import shlex
import sys
import time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Optional, Sequence, TextIO


@dataclass(frozen=True)
class LogConfig:
    quiet: bool = False
    verbose: bool = False
    stream: TextIO = sys.stderr


_DEFAULT_CONFIG = LogConfig()


def _write(message: str, *, stream: Optional[TextIO] = None) -> None:
    handle = stream or _DEFAULT_CONFIG.stream
    print(message, file=handle)


def info(message: str, *, stream: Optional[TextIO] = None) -> None:
    _write(f"[autotax2] {message}", stream=stream)


def step(message: str, *, stream: Optional[TextIO] = None) -> None:
    _write(f"[autotax2] {message}", stream=stream)


def success(message: str, *, stream: Optional[TextIO] = None) -> None:
    _write(f"OK {message}", stream=stream)


def warning(message: str, *, stream: Optional[TextIO] = None) -> None:
    _write(f"WARNING {message}", stream=stream)


def warn(message: str, *, stream: Optional[TextIO] = None) -> None:
    warning(message, stream=stream)


def error(message: str, *, stream: Optional[TextIO] = None) -> None:
    _write(f"ERROR {message}", stream=stream)


def debug(
    message: str,
    *,
    enabled: bool = False,
    stream: Optional[TextIO] = None,
) -> None:
    if enabled:
        _write(f"DEBUG {message}", stream=stream)


def section(title: str, *, stream: Optional[TextIO] = None) -> None:
    _write("", stream=stream)
    _write(title, stream=stream)


def key_value(
    key: str,
    value: object,
    *,
    stream: Optional[TextIO] = None,
) -> None:
    _write(f"  {key:<28} {value}", stream=stream)


def command_to_string(command: Sequence[str | Path]) -> str:
    return " ".join(shlex.quote(str(part)) for part in command)


def command(
    command_parts: Sequence[str | Path],
    *,
    stream: Optional[TextIO] = None,
) -> None:
    _write(f"[command] {command_to_string(command_parts)}", stream=stream)


def print_command(
    command_parts: Sequence[str | Path],
    *,
    stream: Optional[TextIO] = None,
) -> None:
    command(command_parts, stream=stream)


def print_outputs(
    outputs: dict[str, object],
    *,
    stream: Optional[TextIO] = None,
) -> None:
    for key, value in outputs.items():
        _write(f"{key}\t{value}", stream=stream)


@contextmanager
def timed_step(
    message: str,
    *,
    stream: Optional[TextIO] = None,
) -> Iterator[None]:
    start = time.time()
    step(message, stream=stream)

    try:
        yield
    except Exception:
        elapsed = time.time() - start
        error(f"Failed after {elapsed:.1f}s: {message}", stream=stream)
        raise
    else:
        elapsed = time.time() - start
        success(f"{message} ({elapsed:.1f}s)", stream=stream)


class ProgressCounter:
    """Simple text progress reporter for streaming FASTA/TSV processing."""

    def __init__(
        self,
        label: str,
        *,
        interval: int = 100000,
        stream: Optional[TextIO] = None,
    ) -> None:
        if interval < 1:
            raise ValueError("Progress interval must be at least 1.")

        self.label = label
        self.interval = interval
        self.stream = stream
        self.count = 0
        self.start_time = time.time()
        self.last_report_time = self.start_time

    def update(self, n: int = 1) -> None:
        self.count += n

        if self.count % self.interval == 0:
            self.report()

    def report(self) -> None:
        now = time.time()
        elapsed = max(now - self.start_time, 1e-9)
        rate = self.count / elapsed

        _write(
            f"[autotax2] processed {self.count:,} {self.label} "
            f"({rate:,.1f}/s)",
            stream=self.stream,
        )

        self.last_report_time = now

    def finish(self) -> None:
        elapsed = time.time() - self.start_time

        if self.count == 0:
            _write(f"[autotax2] processed 0 {self.label}", stream=self.stream)
            return

        rate = self.count / max(elapsed, 1e-9)

        _write(
            f"[autotax2] processed {self.count:,} {self.label} "
            f"in {elapsed:.1f}s ({rate:,.1f}/s)",
            stream=self.stream,
        )


def summarize_paths(title: str, paths: Iterable[tuple[str, str | Path]]) -> None:
    section(title)

    for key, value in paths:
        key_value(key, value)