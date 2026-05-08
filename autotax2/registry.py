"""Registry models for stable internal IDs and sequence de-duplication."""

from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import md5
from typing import Sequence

from autotax2.io import FastaRecord, normalize_sequence
from autotax2.validate import validate_dataset_prefix


@dataclass(frozen=True)
class DatasetPrefixRecord:
    """A frozen dataset prefix assignment."""

    prefix: str
    dataset_name: str | None = None
    frozen: bool = True


@dataclass
class DatasetPrefixRegistry:
    """Register dataset prefixes and keep assignments frozen."""

    prefixes: dict[str, DatasetPrefixRecord] = field(default_factory=dict)

    def register(self, prefix: str, dataset_name: str | None = None) -> DatasetPrefixRecord:
        """Register a dataset prefix, or return the existing frozen assignment."""
        normalized_prefix = validate_dataset_prefix(prefix)
        normalized_name = dataset_name.strip() if dataset_name is not None else None
        if normalized_name == "":
            normalized_name = None

        existing = self.prefixes.get(normalized_prefix)
        if existing is not None:
            if normalized_name is not None and existing.dataset_name not in (None, normalized_name):
                raise ValueError(
                    f"Dataset prefix {normalized_prefix} is already frozen for "
                    f"{existing.dataset_name!r}."
                )
            return existing

        record = DatasetPrefixRecord(
            prefix=normalized_prefix,
            dataset_name=normalized_name,
            frozen=True,
        )
        self.prefixes[normalized_prefix] = record
        return record

    def require(self, prefix: str) -> DatasetPrefixRecord:
        """Return a registered prefix or raise an error."""
        normalized_prefix = validate_dataset_prefix(prefix)
        try:
            return self.prefixes[normalized_prefix]
        except KeyError as exc:
            raise KeyError(f"Dataset prefix is not registered: {normalized_prefix}") from exc

    def is_registered(self, prefix: str) -> bool:
        """Return whether a dataset prefix has been registered."""
        return validate_dataset_prefix(prefix) in self.prefixes


@dataclass(frozen=True)
class SequenceRecord:
    """A registered sequence and its stable internal metadata."""

    internal_id: str
    original_id: str
    sequence: str
    sequence_md5: str


@dataclass(frozen=True)
class AssignedSequence:
    """A sequence remapped to an internal autotax2 ID."""

    internal_seq_id: str
    original_seq_id: str
    original_header: str
    dataset: str
    prefix: str
    sequence: str
    sequence_md5: str
    sequence_length: int


@dataclass(frozen=True)
class UniqueSequence:
    """One exact normalized sequence represented once in exports."""

    unique_seq_id: str
    sequence_md5: str
    representative_internal_seq_id: str
    sequence: str
    length: int
    first_seen_dataset: str


@dataclass(frozen=True)
class SequenceMembership:
    """Membership of an internal sequence ID in an exact-sequence group."""

    internal_seq_id: str
    original_seq_id: str
    dataset: str
    prefix: str
    sequence_md5: str
    unique_seq_id: str
    is_duplicate_sequence: bool


@dataclass
class SequenceRegistry:
    """Track unique sequences and duplicate memberships by normalized MD5."""

    unique_sequences: dict[str, UniqueSequence] = field(default_factory=dict)
    memberships: dict[str, SequenceMembership] = field(default_factory=dict)
    md5_index: dict[str, str] = field(default_factory=dict)

    def add(self, assigned: AssignedSequence) -> SequenceMembership:
        """Add one assigned sequence and record its unique-sequence membership."""
        if assigned.internal_seq_id in self.memberships:
            raise ValueError(f"Internal sequence ID already exists: {assigned.internal_seq_id}")

        unique_seq_id = self.md5_index.get(assigned.sequence_md5)
        is_duplicate = unique_seq_id is not None

        if unique_seq_id is None:
            unique_seq_id = self._next_unique_seq_id()
            self.unique_sequences[unique_seq_id] = UniqueSequence(
                unique_seq_id=unique_seq_id,
                sequence_md5=assigned.sequence_md5,
                representative_internal_seq_id=assigned.internal_seq_id,
                sequence=assigned.sequence,
                length=assigned.sequence_length,
                first_seen_dataset=assigned.dataset,
            )
            self.md5_index[assigned.sequence_md5] = unique_seq_id

        membership = SequenceMembership(
            internal_seq_id=assigned.internal_seq_id,
            original_seq_id=assigned.original_seq_id,
            dataset=assigned.dataset,
            prefix=assigned.prefix,
            sequence_md5=assigned.sequence_md5,
            unique_seq_id=unique_seq_id,
            is_duplicate_sequence=is_duplicate,
        )
        self.memberships[assigned.internal_seq_id] = membership
        return membership

    def add_many(self, assigned_records: Sequence[AssignedSequence]) -> list[SequenceMembership]:
        """Add assigned sequences in order."""
        return [self.add(assigned) for assigned in assigned_records]

    def unique_records(self) -> list[UniqueSequence]:
        """Return unique sequence records in stable ID order."""
        return [self.unique_sequences[key] for key in sorted(self.unique_sequences)]

    def membership_records(self) -> list[SequenceMembership]:
        """Return membership records in internal sequence ID order."""
        return [self.memberships[key] for key in sorted(self.memberships)]

    def representative_internal_ids(self) -> list[str]:
        """Return one representative internal ID per unique sequence."""
        return [
            unique.representative_internal_seq_id
            for unique in self.unique_records()
        ]

    def _next_unique_seq_id(self) -> str:
        return f"U_{len(self.unique_sequences) + 1:06d}"


@dataclass
class Registry:
    """In-memory placeholder registry used until durable storage is designed."""

    records: dict[str, SequenceRecord] = field(default_factory=dict)
    md5_index: dict[str, str] = field(default_factory=dict)
    dataset_prefixes: DatasetPrefixRegistry = field(default_factory=DatasetPrefixRegistry)
    _next_ordinals: dict[str, int] = field(default_factory=dict, init=False, repr=False)

    def add_sequence(
        self,
        dataset_prefix: str,
        ordinal: int,
        original_id: str,
        sequence: str,
    ) -> SequenceRecord:
        """Add a sequence with a deterministic internal ID."""
        prefix_record = self.dataset_prefixes.register(dataset_prefix)
        internal_id = remap_sequence_id(prefix_record.prefix, ordinal)
        if internal_id in self.records:
            raise ValueError(f"Internal sequence ID already exists: {internal_id}")

        normalized_sequence = normalize_sequence(sequence)
        digest = sequence_md5(normalized_sequence)
        record = SequenceRecord(
            internal_id=internal_id,
            original_id=original_id,
            sequence=normalized_sequence,
            sequence_md5=digest,
        )
        self.records[internal_id] = record
        self.md5_index.setdefault(digest, internal_id)
        self._next_ordinals[prefix_record.prefix] = max(
            self._next_ordinals.get(prefix_record.prefix, 1),
            ordinal + 1,
        )
        return record

    def add_sequence_auto(
        self,
        dataset_prefix: str,
        original_id: str,
        sequence: str,
    ) -> SequenceRecord:
        """Add a sequence using the next ordinal for the dataset prefix."""
        prefix_record = self.dataset_prefixes.register(dataset_prefix)
        ordinal = self._next_ordinals.get(prefix_record.prefix, 1)
        return self.add_sequence(prefix_record.prefix, ordinal, original_id, sequence)

    def register_dataset_prefix(
        self,
        prefix: str,
        dataset_name: str | None = None,
    ) -> DatasetPrefixRecord:
        """Freeze a dataset prefix assignment."""
        return self.dataset_prefixes.register(prefix, dataset_name)

    def exported_internal_ids(self) -> list[str]:
        """Return one internal ID per exact sequence MD5."""
        return sorted(self.md5_index.values())


def remap_sequence_id(dataset_prefix: str, ordinal: int) -> str:
    """Map an input sequence to an internal ID such as ``D20_000001``."""
    if ordinal < 1:
        raise ValueError("Ordinal must be >= 1.")
    normalized_prefix = validate_dataset_prefix(dataset_prefix)
    return f"{normalized_prefix}_{ordinal:06d}"


def assign_internal_ids(
    records: Sequence[FastaRecord],
    dataset_name: str,
    prefix: str,
) -> list[AssignedSequence]:
    """Assign deterministic internal IDs to FASTA records from one dataset."""
    normalized_prefix = validate_dataset_prefix(prefix)
    normalized_dataset = dataset_name.strip()
    if not normalized_dataset:
        raise ValueError("Dataset name must not be empty.")

    assigned: list[AssignedSequence] = []
    for ordinal, record in enumerate(records, start=1):
        sequence = normalize_sequence(record.sequence)
        assigned.append(
            AssignedSequence(
                internal_seq_id=remap_sequence_id(normalized_prefix, ordinal),
                original_seq_id=record.seq_id,
                original_header=record.header,
                dataset=normalized_dataset,
                prefix=normalized_prefix,
                sequence=sequence,
                sequence_md5=sequence_md5(sequence),
                sequence_length=len(sequence),
            )
        )
    return assigned


def sequence_md5(sequence: str) -> str:
    """Return the MD5 hex digest for a normalized sequence string."""
    normalized = normalize_sequence(sequence)
    return md5(normalized.encode("ascii")).hexdigest()
