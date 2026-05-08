"""Placeholder ID helpers."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum

from autotax2.validate import validate_dataset_prefix


PLACEHOLDER_ID_RE = re.compile(
    r"^(?P<rank>[cofgs])__(?P<prefix>[A-Za-z][A-Za-z0-9]*)(?P=rank)(?P<ordinal>[0-9]{6})$"
)


class PlaceholderRank(str, Enum):
    """Ranks that currently receive generated placeholder IDs."""

    CLASS = "c"
    ORDER = "o"
    FAMILY = "f"
    GENUS = "g"
    SPECIES = "s"


@dataclass(frozen=True)
class PlaceholderId:
    """A parsed placeholder identifier."""

    rank: PlaceholderRank
    dataset_prefix: str
    ordinal: int

    def as_string(self) -> str:
        """Render the placeholder ID."""
        return make_placeholder_id(self.rank, self.dataset_prefix, self.ordinal)


@dataclass
class PlaceholderAllocator:
    """Allocate stable placeholder IDs without reusing deprecated IDs."""

    issued_ids: set[str] = field(default_factory=set)
    deprecated_ids: set[str] = field(default_factory=set)
    _next_ordinals: dict[tuple[PlaceholderRank, str], int] = field(
        default_factory=dict,
        init=False,
        repr=False,
    )

    def __post_init__(self) -> None:
        self.issued_ids = set(self.issued_ids)
        self.deprecated_ids = set(self.deprecated_ids)

        for placeholder_id in self.issued_ids | self.deprecated_ids:
            parsed = parse_placeholder_id(placeholder_id)
            self._advance_next(parsed)

    def allocate(self, rank: PlaceholderRank | str, dataset_prefix: str) -> str:
        """Allocate the next active placeholder for a rank and dataset prefix."""
        parsed_rank = PlaceholderRank(rank)
        normalized_prefix = validate_dataset_prefix(dataset_prefix)
        key = (parsed_rank, normalized_prefix)
        ordinal = self._next_ordinals.get(key, 1)

        while True:
            placeholder_id = make_placeholder_id(parsed_rank, normalized_prefix, ordinal)
            if placeholder_id not in self.issued_ids and placeholder_id not in self.deprecated_ids:
                self.issued_ids.add(placeholder_id)
                self._next_ordinals[key] = ordinal + 1
                return placeholder_id
            ordinal += 1

    def reserve(self, placeholder_id: str) -> None:
        """Reserve an existing placeholder ID as active."""
        if placeholder_id in self.deprecated_ids:
            raise ValueError(f"Cannot reserve deprecated placeholder ID: {placeholder_id}")
        parsed = parse_placeholder_id(placeholder_id)
        self.issued_ids.add(parsed.as_string())
        self._advance_next(parsed)

    def deprecate(self, placeholder_id: str) -> None:
        """Mark a placeholder ID as deprecated so it can never be reused."""
        parsed = parse_placeholder_id(placeholder_id)
        rendered = parsed.as_string()
        self.issued_ids.discard(rendered)
        self.deprecated_ids.add(rendered)
        self._advance_next(parsed)

    def next_ordinal(self, rank: PlaceholderRank | str, dataset_prefix: str) -> int:
        """Return the next ordinal that would be considered for allocation."""
        parsed_rank = PlaceholderRank(rank)
        normalized_prefix = validate_dataset_prefix(dataset_prefix)
        return self._next_ordinals.get((parsed_rank, normalized_prefix), 1)

    def is_deprecated(self, placeholder_id: str) -> bool:
        """Return whether a placeholder ID has been deprecated."""
        return parse_placeholder_id(placeholder_id).as_string() in self.deprecated_ids

    def _advance_next(self, parsed: PlaceholderId) -> None:
        key = (parsed.rank, parsed.dataset_prefix)
        self._next_ordinals[key] = max(self._next_ordinals.get(key, 1), parsed.ordinal + 1)


def make_placeholder_id(rank: PlaceholderRank | str, dataset_prefix: str, ordinal: int) -> str:
    """Create a placeholder ID such as ``g__D20g000001``."""
    parsed_rank = PlaceholderRank(rank)
    if ordinal < 1:
        raise ValueError("Ordinal must be >= 1.")

    normalized_prefix = validate_dataset_prefix(dataset_prefix)

    return f"{parsed_rank.value}__{normalized_prefix}{parsed_rank.value}{ordinal:06d}"


def parse_placeholder_id(placeholder_id: str) -> PlaceholderId:
    """Parse a placeholder ID created by :func:`make_placeholder_id`."""
    match = PLACEHOLDER_ID_RE.fullmatch(placeholder_id.strip())
    if match is None:
        raise ValueError(f"Invalid placeholder ID: {placeholder_id}")

    return PlaceholderId(
        rank=PlaceholderRank(match.group("rank")),
        dataset_prefix=validate_dataset_prefix(match.group("prefix")),
        ordinal=int(match.group("ordinal")),
    )
