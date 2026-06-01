"""Taxonomy parsing and placeholder helpers."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

from .config import RANK_PREFIX, RANKS_COARSE_TO_FINE, RANKS_FINE_TO_COARSE


@dataclass
class AnchorCall:
    seq_id: str
    anchor_rank: str | None
    anchor_taxonomy: Dict[str, str | None]
    align_ident: float | None
    align_quality: float | None
    reason: str = ""


def normalize_taxon(rank: str, value: str | None) -> str | None:
    if value is None:
        return None
    value = value.strip()
    if not value or value.lower() in {"na", "nan", "none", "unclassified", "unclassified;"}:
        return None
    prefix = RANK_PREFIX.get(rank)
    if prefix and not value.startswith(f"{prefix}__"):
        if ":" in value:
            value = value.split(":", 1)[1]
        value = f"{prefix}__{value}"
    return value.rstrip(";")


def parse_tax_string(tax: str | None) -> Dict[str, str | None]:
    """Parse semicolon taxonomy into rank dictionary.

    Accepts GTDB/SILVA-like strings such as:
    d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;...
    """
    result: Dict[str, str | None] = {rank: None for rank in ["domain"] + RANKS_COARSE_TO_FINE}
    if tax is None:
        return result
    tax = tax.strip().strip(";")
    if not tax or tax.lower() == "unclassified":
        return result
    parts = [p.strip() for p in tax.split(";") if p.strip()]
    prefix_to_rank = {v: k for k, v in RANK_PREFIX.items()}
    for part in parts:
        if "__" in part:
            pref = part.split("__", 1)[0]
            rank = prefix_to_rank.get(pref)
            if rank:
                result[rank] = normalize_taxon(rank, part)
        elif ":" in part:
            pref, name = part.split(":", 1)
            rank = prefix_to_rank.get(pref)
            if rank:
                result[rank] = normalize_taxon(rank, name)
    return result


def taxonomy_to_string(tax: Dict[str, str | None], include_domain: bool = True) -> str:
    ranks = (["domain"] if include_domain else []) + RANKS_COARSE_TO_FINE
    return ";".join(tax.get(rank) or "" for rank in ranks)


def infer_anchor_rank(
    align_ident_percent: float | None,
    lca_tax: Dict[str, str | None],
    thresholds_percent: Dict[str, float],
) -> str | None:
    """Return the finest trusted anchor rank based on SINA identity and LCA taxonomy."""
    if align_ident_percent is None:
        return None
    for rank in RANKS_FINE_TO_COARSE:
        if align_ident_percent >= thresholds_percent[rank] and lca_tax.get(rank):
            return rank
    return None


def inherited_ranks(anchor_rank: str | None) -> list[str]:
    """Ranks inherited from backbone for an anchor rank."""
    if anchor_rank is None:
        return []
    idx = RANKS_COARSE_TO_FINE.index(anchor_rank)
    return RANKS_COARSE_TO_FINE[: idx + 1]


def novel_ranks_below_anchor(anchor_rank: str | None) -> list[str]:
    """Ranks that need novel clustering/placeholder construction."""
    if anchor_rank is None:
        return RANKS_COARSE_TO_FINE.copy()
    idx = RANKS_COARSE_TO_FINE.index(anchor_rank)
    return RANKS_COARSE_TO_FINE[idx + 1 :]


def make_placeholder(rank: str, prefix: str, number: int, digits: int = 6) -> str:
    letter = RANK_PREFIX[rank]
    return f"{letter}__{prefix}_{letter}{number:0{digits}d}"
