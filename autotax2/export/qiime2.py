"""QIIME2 export helpers."""

from __future__ import annotations

from collections.abc import Sequence

from autotax2.export.formats import format_taxonomy_qiime2


QIIME2_REFERENCE_SEQUENCES_FORMAT = "qiime2-reference-sequences"
QIIME2_TAXONOMY_FORMAT = "qiime2-taxonomy"
QIIME2_FORMATS = (QIIME2_REFERENCE_SEQUENCES_FORMAT, QIIME2_TAXONOMY_FORMAT)


def format_qiime2_taxonomy_row(sequence_id: str, taxonomy: str | Sequence[str]) -> str:
    """Format a QIIME2 taxonomy TSV row."""
    return f"{sequence_id}\t{format_taxonomy_qiime2(taxonomy)}"
