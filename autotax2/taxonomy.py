"""Taxonomy primitives for rank-aware reference building."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class Rank(str, Enum):
    """Supported taxonomic ranks."""

    DOMAIN = "d"
    PHYLUM = "p"
    CLASS = "c"
    ORDER = "o"
    FAMILY = "f"
    GENUS = "g"
    SPECIES = "s"


RANK_ORDER = (
    Rank.DOMAIN,
    Rank.PHYLUM,
    Rank.CLASS,
    Rank.ORDER,
    Rank.FAMILY,
    Rank.GENUS,
    Rank.SPECIES,
)
RANK_INDEX = {rank: index for index, rank in enumerate(RANK_ORDER)}


@dataclass(frozen=True)
class Taxon:
    """A taxon label at a specific rank."""

    rank: Rank | str
    name: str

    def __post_init__(self) -> None:
        rank = Rank(self.rank)
        name = self.name.strip()
        if not name:
            raise ValueError("Taxon name must not be empty.")
        object.__setattr__(self, "rank", rank)
        object.__setattr__(self, "name", name)

    @property
    def prefixed_name(self) -> str:
        """Return a rank-prefixed taxonomy label."""
        return f"{self.rank.value}__{self.name}"


@dataclass(frozen=True)
class TaxonomyPath:
    """A parsed, rank-aware taxonomy path."""

    taxa: tuple[Taxon, ...]

    @property
    def labels(self) -> tuple[str, ...]:
        """Return rank-prefixed labels in path order."""
        return tuple(taxon.prefixed_name for taxon in self.taxa)

    def by_rank(self) -> dict[Rank, Taxon]:
        """Return taxa keyed by rank."""
        return {taxon.rank: taxon for taxon in self.taxa}

    def get(self, rank: Rank | str) -> Taxon | None:
        """Return a taxon by rank if present."""
        return self.by_rank().get(Rank(rank))

    def require(self, rank: Rank | str) -> Taxon:
        """Return a taxon by rank or raise an error."""
        parsed_rank = Rank(rank)
        taxon = self.get(parsed_rank)
        if taxon is None:
            raise KeyError(f"Taxonomy path lacks rank: {parsed_rank.value}")
        return taxon

    def to_semicolon_string(self, trailing_semicolon: bool = False) -> str:
        """Render the taxonomy path as a semicolon-delimited string."""
        rendered = ";".join(self.labels)
        if trailing_semicolon and rendered:
            return f"{rendered};"
        return rendered


def split_taxonomy_path(path: str) -> list[str]:
    """Split a semicolon-delimited taxonomy string into non-empty labels."""
    return [part.strip() for part in path.split(";") if part.strip()]


def rank_prefix(label: str) -> str:
    """Return the rank prefix from a taxonomy label."""
    return parse_taxon_label(label).rank.value


def parse_taxon_label(label: str) -> Taxon:
    """Parse a single rank-prefixed taxonomy label."""
    normalized = label.strip()
    if "__" not in normalized:
        raise ValueError(f"Taxonomy label lacks rank prefix: {label}")

    prefix, name = normalized.split("__", maxsplit=1)
    try:
        rank = Rank(prefix)
    except ValueError as exc:
        raise ValueError(f"Unsupported taxonomy rank prefix: {prefix}") from exc

    return Taxon(rank=rank, name=name)


def parse_taxonomy_path(path: str, strict_order: bool = True) -> TaxonomyPath:
    """Parse a semicolon-delimited taxonomy path into rank-aware taxa."""
    taxa = tuple(parse_taxon_label(label) for label in split_taxonomy_path(path))
    seen: set[Rank] = set()
    previous_index = -1

    for taxon in taxa:
        if taxon.rank in seen:
            raise ValueError(f"Duplicate taxonomy rank: {taxon.rank.value}")
        seen.add(taxon.rank)

        index = RANK_INDEX[taxon.rank]
        if strict_order and index <= previous_index:
            raise ValueError("Taxonomy ranks are out of canonical order.")
        previous_index = index

    return TaxonomyPath(taxa=taxa)
