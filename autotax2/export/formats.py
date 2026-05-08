"""Taxonomy formatting helpers for classifier exports."""

from __future__ import annotations

from collections.abc import Sequence
import re


RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
RANK_PREFIXES = ("d", "p", "c", "o", "f", "g", "s")


def validate_taxonomy_7rank(taxonomy: str | Sequence[str]) -> tuple[str, ...]:
    """Return a normalized seven-rank taxonomy or raise a clear validation error."""
    if isinstance(taxonomy, str):
        ranks = tuple(part.strip() for part in taxonomy.split(";") if part.strip())
    else:
        ranks = tuple(str(part).strip() for part in taxonomy if str(part).strip())

    if len(ranks) != 7:
        raise ValueError(f"Taxonomy must contain exactly 7 ranks; observed {len(ranks)}: {taxonomy}")

    normalized: list[str] = []
    for prefix, value in zip(RANK_PREFIXES, ranks, strict=True):
        if not value:
            raise ValueError(f"Taxonomy contains an empty {prefix} rank: {taxonomy}")
        if re.match(r"^[a-z]__", value):
            normalized.append(value)
        else:
            normalized.append(f"{prefix}__{value}")
    return tuple(normalized)


def strip_rank_prefix(value: str) -> str:
    """Strip a rank prefix such as ``g__`` from a taxonomy value."""
    return re.sub(r"^[a-z]__", "", value.strip())


def format_taxonomy_sintax(taxonomy: str | Sequence[str]) -> str:
    """Format taxonomy for VSEARCH SINTAX FASTA headers."""
    ranks = validate_taxonomy_7rank(taxonomy)
    return ",".join(
        f"{prefix}:{strip_rank_prefix(value)}"
        for prefix, value in zip(RANK_PREFIXES, ranks, strict=True)
    )


def format_taxonomy_qiime2(taxonomy: str | Sequence[str]) -> str:
    """Format taxonomy for QIIME2 reference taxonomy TSV."""
    return "; ".join(validate_taxonomy_7rank(taxonomy))


def format_taxonomy_dada2_genus(taxonomy: str | Sequence[str]) -> str:
    """Format taxonomy to genus level for DADA2 assignTaxonomy training FASTA."""
    ranks = validate_taxonomy_7rank(taxonomy)
    return ";".join(strip_rank_prefix(value) for value in ranks[:6])


def format_taxonomy_dada2_species(taxonomy: str | Sequence[str]) -> str:
    """Format genus/species label for DADA2 assignSpecies FASTA headers."""
    ranks = validate_taxonomy_7rank(taxonomy)
    genus = strip_rank_prefix(ranks[5])
    species = strip_rank_prefix(ranks[6])
    if not genus or not species:
        raise ValueError("DADA2 assignSpecies requires non-empty genus and species ranks.")
    if species.startswith(f"{genus} "):
        return species
    return f"{genus} {species}"
