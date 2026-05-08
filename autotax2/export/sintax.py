"""SINTAX export helpers."""

from __future__ import annotations

from collections.abc import Sequence

from autotax2.export.formats import format_taxonomy_sintax


SINTAX_FORMAT = "sintax-reference"


def format_sintax_header(sequence_id: str, taxonomy: str | Sequence[str]) -> str:
    """Format a SINTAX-compatible FASTA header."""
    return f"{sequence_id};tax={format_taxonomy_sintax(taxonomy)};"
