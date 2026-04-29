from __future__ import annotations

import os


DEFAULT_THREADS = 4


class ThreadValueError(ValueError):
    """Raised when a thread value is not a positive integer."""


def available_cpu_count() -> int:
    """Return the number of CPUs visible to Python.

    This is informational only. AutoTax2 does not use CPU count as the
    default thread number.
    """

    count = os.cpu_count()
    if count is None or count < 1:
        return 1
    return count


def validate_threads(threads: int) -> int:
    """Validate and return a positive integer thread count.

    AutoTax2 only accepts explicit positive integers.

    Valid:
        1
        4
        12

    Invalid:
        "auto"
        "default"
        "4"
        0
        -1
        None
    """

    if not isinstance(threads, int):
        raise ThreadValueError(
            f"Thread count must be a positive integer, not {type(threads).__name__}."
        )

    if threads < 1:
        raise ThreadValueError("Thread count must be at least 1.")

    return threads


def parse_threads(text: str) -> int:
    """Parse a thread count from config text.

    This function is only for reading config files. It accepts numeric strings,
    but rejects 'auto' and all non-integer values.
    """

    value = text.strip()

    if not value:
        raise ThreadValueError("Thread count is empty.")

    try:
        threads = int(value)
    except ValueError as exc:
        raise ThreadValueError(
            f"Thread count must be a positive integer, not {text!r}."
        ) from exc

    return validate_threads(threads)


def default_threads() -> int:
    """Return the AutoTax2 default thread count."""

    return DEFAULT_THREADS