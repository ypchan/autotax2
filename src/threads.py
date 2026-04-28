from __future__ import annotations

import os
from typing import Optional

from .logging import warn


def available_cpus() -> int:
    try:
        return os.cpu_count() or 1
    except Exception:
        return 1


def resolve_threads(value: Optional[int | str], reserve: int = 1, min_threads: int = 1) -> int:
    """Resolve a user thread request into a safe integer.

    Accepted values:
    - None or "auto": max(1, os.cpu_count() - reserve)
    - positive integer string/int: requested thread count
    - values <= 0: treated as auto

    The result is capped at os.cpu_count() to avoid accidental oversubscription.
    """
    cpus = available_cpus()
    auto_threads = max(min_threads, cpus - max(0, reserve))

    if value is None:
        requested = auto_threads
    elif isinstance(value, str) and value.lower() == "auto":
        requested = auto_threads
    else:
        try:
            requested = int(value)  # type: ignore[arg-type]
        except Exception as exc:
            raise ValueError(f"Invalid --threads value: {value}") from exc
        if requested <= 0:
            requested = auto_threads

    if requested > cpus:
        warn(f"--threads {requested} is larger than detected CPUs {cpus}; capping to {cpus}.")
        requested = cpus
    if requested < min_threads:
        requested = min_threads
    return requested


def per_job_threads(total_threads: int, jobs: int = 1, min_threads: int = 1) -> int:
    """Return threads per job for workflows that run multiple jobs.

    AutoTax2 currently runs VSEARCH jobs sequentially, but this helper is used
    to keep future parallel multi-level execution from oversubscribing CPUs.
    """
    jobs = max(1, int(jobs))
    return max(min_threads, total_threads // jobs)
