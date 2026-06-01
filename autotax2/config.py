"""Configuration constants for AutoTax2."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

RANKS_FINE_TO_COARSE: List[str] = [
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
]

RANKS_COARSE_TO_FINE: List[str] = list(reversed(RANKS_FINE_TO_COARSE))

RANK_PREFIX = {
    "domain": "d",
    "phylum": "p",
    "class": "c",
    "order": "o",
    "family": "f",
    "genus": "g",
    "species": "s",
}

DEFAULT_THRESHOLDS_PERCENT: Dict[str, float] = {
    "species": 97.2,
    "genus": 90.1,
    "family": 80.1,
    "order": 72.9,
    "class": 72.2,
    "phylum": 69.6,
}

DEFAULT_THRESHOLDS_FRACTION: Dict[str, float] = {
    rank: value / 100.0 for rank, value in DEFAULT_THRESHOLDS_PERCENT.items()
}

DEFAULT_REQUIRED_TAX_COLUMNS = [
    "seq_id",
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

@dataclass
class RuntimeConfig:
    """Runtime configuration persisted into database config.yaml."""

    db: Path
    ref_fa: Path | None = None
    ref_tax: Path | None = None
    ref_arb: Path | None = None
    threads: int = 1
    thresholds_percent: Dict[str, float] = field(default_factory=lambda: DEFAULT_THRESHOLDS_PERCENT.copy())
    iddef: int = 2
    min_align_quality: float = 0.0
    min_aligned_len: int = 0
    max_head_cutoff: int | None = None
    max_tail_cutoff: int | None = None
    min_qcov: float = 0.80
    min_tcov: float = 0.80
    placeholder_digits: int = 6
    vsearch_bin: str = "vsearch"
    sina_bin: str = "sina"
