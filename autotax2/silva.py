"""SILVA backbone initialization helpers."""

from __future__ import annotations

import csv
import gzip
import re
import hashlib
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, TextIO

import yaml

from autotax2.placeholders import (
    PlaceholderAllocator,
    PlaceholderRank,
    make_placeholder_id,
    parse_placeholder_id,
)
from autotax2.io import FastaRecord, read_fasta, write_fasta
from autotax2.registry import sequence_md5
from autotax2.thresholds import (
    DEFAULT_CLASS_ID,
    DEFAULT_FAMILY_ID,
    DEFAULT_GENUS_ID,
    DEFAULT_ORDER_ID,
    DEFAULT_PHYLUM_ID,
    DEFAULT_SPECIES_ID,
)
from autotax2.vsearch import (
    parse_uc_records,
    read_uc_assignments,
    run_vsearch_cluster,
    write_singleton_uc,
)


SILVA_PLACEHOLDER_PREFIX = "SILVA"
SILVA_SOURCE = "SILVA138.2_NR99"
SILVA_RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
BUILD_MANIFEST_FIELDS = [
    "build_dir",
    "silva_fasta",
    "type_strain_metadata",
    "gtdb_taxonomy",
    "threads",
    "metadata_contract",
    "strict_legal_names",
    "organelle_policy",
]
LEGAL_NAME_CATALOG_FIELDS = ["rank", "name", "source", "canonical_species_name"]
RANK_PREFIXES = {
    "domain": "d",
    "phylum": "p",
    "class": "c",
    "order": "o",
    "family": "f",
    "genus": "g",
    "species": "s",
}
PLACEHOLDER_RANKS = ("phylum", "class", "order", "family", "genus", "species")
PLACEHOLDER_RANK_BY_NAME = {
    "phylum": PlaceholderRank.PHYLUM,
    "class": PlaceholderRank.CLASS,
    "order": PlaceholderRank.ORDER,
    "family": PlaceholderRank.FAMILY,
    "genus": PlaceholderRank.GENUS,
    "species": PlaceholderRank.SPECIES,
}
RESOLVE_RANKS = ("phylum", "class", "order", "family", "genus", "species")
RESOLVE_EVIDENCE_FIELDS = [
    "source_stage",
    "dataset",
    "seq_id",
    "rank",
    "parent_rank",
    "parent_taxon",
    "input_rank_value",
    "input_rank_status",
    "threshold",
    "candidate_scope",
    "anchor_count",
    "best_anchor_taxon",
    "best_anchor_identity",
    "passing_anchor_count",
    "competing_anchor_taxa",
    "residual_cluster_key",
    "residual_cluster_size",
    "decision",
    "output_taxon",
    "reason",
    "job_size",
    "job_threads",
]
MIXED_PARENT_WARNING = "mixed_silva_parent_unresolved_cluster"
SILVA_REJECTED_FIELDS = ["seq_id", "original_taxonomy", "reject_reason"]
SILVA_INVALID_SEQUENCE_REASON = "invalid_sequence_characters"
SILVA_ATGC_RE = re.compile(r"^[ATGC]*$")
DEFAULT_FLOOR_ID = DEFAULT_PHYLUM_ID
UNRESOLVED_TOKENS = (
    "unidentified",
    "unclassified",
    "uncultured",
    "uncultivated",
    "unknown",
    "environmental",
    "metagenome",
    "metagenomic",
    "sample",
    "clone",
    "bacterium",
    "archaeon",
    "organism",
)
STANDALONE_UNRESOLVED_VALUES = {
    "candidatus",
    "incertae",
    "incertae sedis",
}
ENVIRONMENTAL_DESCRIPTOR_VALUES = {
    "aquatic",
    "biofilm",
    "deep sea",
    "fecal",
    "faecal",
    "freshwater",
    "gut",
    "hydrothermal",
    "marine",
    "oral",
    "rumen",
    "sediment",
    "sludge",
    "soil",
    "thermal",
    "wastewater",
}
ORGANELLE_MARKERS = {
    "chloroplast": "Chloroplast",
    "mitochondria": "Mitochondria",
    "mitochondrion": "Mitochondria",
    "plastid": "Plastid",
    "apicoplast": "Apicoplast",
    "cyanelle": "Cyanelle",
}
GTDB_PREFIX_TO_RANK = {
    "d": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
MISSING_VALUE_PREFIX = "unresolved_"


@dataclass(frozen=True)
class SilvaTaxonomy:
    """Seven-rank SILVA taxonomy parsed from a FASTA header."""

    domain: str
    phylum: str
    class_name: str
    order: str
    family: str
    genus: str
    species: str
    original_taxonomy: str

    @property
    def values(self) -> tuple[str, str, str, str, str, str, str]:
        """Return taxonomy values in canonical rank order."""
        return (
            self.domain,
            self.phylum,
            self.class_name,
            self.order,
            self.family,
            self.genus,
            self.species,
        )

    @property
    def taxonomy_7rank(self) -> str:
        """Return the seven-rank taxonomy string used in TSV outputs."""
        return ";".join(self.values)

    def as_rank_dict(self) -> dict[str, str]:
        """Return taxonomy values keyed by TSV rank names."""
        return dict(zip(SILVA_RANKS, self.values, strict=True))


@dataclass(frozen=True)
class UnresolvedClassification:
    """Rank-aware unresolved taxonomy classification."""

    lowest_reliable_rank: str
    unresolved_ranks: tuple[str, ...]
    unresolved_reason: str

    @property
    def is_unresolved(self) -> bool:
        """Return whether any rank is unresolved."""
        return bool(self.unresolved_ranks)


@dataclass(frozen=True)
class SilvaRecord:
    """A SILVA FASTA record with parsed taxonomy and classification."""

    fasta_record: FastaRecord
    taxonomy: SilvaTaxonomy
    classification: UnresolvedClassification

    @property
    def seq_id(self) -> str:
        """Return the SILVA sequence ID."""
        return self.fasta_record.seq_id

    @property
    def is_unresolved(self) -> bool:
        """Return whether the record has unresolved taxonomy."""
        return self.classification.is_unresolved


@dataclass(frozen=True)
class TypeStrainMetadata:
    """Optional type-strain metadata attached by SILVA sequence ID."""

    seq_id: str
    is_type_strain: str
    species_name: str
    strain_id: str
    source: str
    evidence: str


@dataclass(frozen=True)
class LegalNameCatalog:
    """Taxon names accepted by SILVA type strains or GTDB taxonomy."""

    names_by_rank: dict[str, frozenset[str]]
    canonical_species_by_name: dict[str, str]
    sources_by_rank_name: dict[tuple[str, str], frozenset[str]]
    display_names_by_rank_name: dict[tuple[str, str], str]

    @property
    def has_names(self) -> bool:
        """Return whether the catalog contains any type-strain-derived names."""
        return any(self.names_by_rank.get(rank, frozenset()) for rank in SILVA_RANKS)

    def accepts(self, rank: str, value: str) -> bool:
        """Return whether a rank value is accepted by the legal-name catalog."""
        if rank == "species":
            if self.canonical_species(value):
                return True
            return _normalize_taxon_value(value) in self.names_by_rank.get(rank, frozenset())
        return _normalize_taxon_value(value) in self.names_by_rank.get(rank, frozenset())

    def canonical_species(self, value: str) -> str:
        """Return the accepted binomial species name, without strain suffixes."""
        canonical = _canonical_species_name(value)
        if not canonical:
            return ""
        return self.canonical_species_by_name.get(_normalize_taxon_value(canonical), "")


@dataclass(frozen=True)
class ClusterAssignment:
    """Resolved cluster assignment for one sequence."""

    cluster_key: str
    cluster_id: str
    representative_seq_id: str


@dataclass(frozen=True)
class ResolveAnchor:
    """Known same-parent child used as an anchor during SILVA resolve."""

    job_id: str
    seq_id: str
    rank: str
    taxon_name: str
    taxon_label: str
    sequence: str
    source: str


@dataclass(frozen=True)
class ResolveNamedAnchorRecord:
    """Named SILVA sequence that can anchor same-parent unresolved records."""

    seq_id: str
    values: tuple[str, str, str, str, str, str, str]
    sequence: str


@dataclass
class ResolveRecordState:
    """Mutable per-record state while resolving one SILVA unresolved sequence."""

    record: FastaRecord
    source_row: dict[str, str]
    unresolved_ranks: tuple[str, ...]
    resolved_labels: dict[str, str]
    placeholder_taxa: list[str] = field(default_factory=list)
    rank_cluster_keys: dict[str, str] = field(default_factory=dict)
    rank_placeholders: dict[str, str] = field(default_factory=dict)
    warnings: set[str] = field(default_factory=set)
    blocked: bool = False


@dataclass(frozen=True)
class ResolveParentJob:
    """One same-parent, same-rank SILVA resolve job."""

    rank: str
    parent_key: tuple[str, ...]
    candidates: list[ResolveRecordState]
    anchors: list[ResolveAnchor]
    threshold: float
    cluster_dir: Path
    job_threads: int
    vsearch_bin: str
    iddef: int


@dataclass(frozen=True)
class ResolveParentJobResult:
    """Parsed result from a parent-local resolve job."""

    rank: str
    parent_key: tuple[str, ...]
    candidates: list[ResolveRecordState]
    anchors: list[ResolveAnchor]
    threshold: float
    cluster_members: dict[str, list[str]]
    job_id_to_anchor: dict[str, ResolveAnchor]
    candidate_anchor_identities: dict[tuple[str, str], str]
    job_size: int
    job_threads: int


@dataclass(frozen=True)
class ResolveSilvaSummary:
    """Summary from resolving SILVA unresolved records."""

    unresolved_records: int
    placeholder_taxa: int
    dry_run: bool


def named_backbone_is_mutable() -> bool:
    """Return whether the named SILVA backbone may be mutated."""
    return False


def parse_silva_header(header: str) -> tuple[str, SilvaTaxonomy]:
    """Parse a SILVA-style FASTA header into sequence ID and seven-rank taxonomy."""
    normalized_header = header.strip()
    if not normalized_header:
        raise ValueError("SILVA FASTA header must not be empty.")

    parts = normalized_header.split(maxsplit=1)
    seq_id = parts[0]
    original_taxonomy = parts[1].strip() if len(parts) > 1 else ""
    taxonomy_parts = [part.strip() for part in original_taxonomy.split(";")]
    padded_parts = list(taxonomy_parts[: len(SILVA_RANKS)])
    for rank in SILVA_RANKS[len(padded_parts) :]:
        padded_parts.append(f"{MISSING_VALUE_PREFIX}{rank}")

    return seq_id, SilvaTaxonomy(
        domain=padded_parts[0],
        phylum=padded_parts[1],
        class_name=padded_parts[2],
        order=padded_parts[3],
        family=padded_parts[4],
        genus=padded_parts[5],
        species=padded_parts[6],
        original_taxonomy=original_taxonomy,
    )


def classify_unresolved(
    taxonomy: SilvaTaxonomy,
    legal_names: LegalNameCatalog | None = None,
) -> UnresolvedClassification:
    """Detect unresolved ranks in a seven-rank SILVA taxonomy."""
    direct_unresolved: list[str] = []
    reasons: list[str] = []

    for rank, value in zip(SILVA_RANKS, taxonomy.values, strict=True):
        reason = _unresolved_reason(rank, value)
        if not reason and legal_names is not None and legal_names.has_names:
            if not legal_names.accepts(rank, value):
                reason = "not in legal name catalog"
        if reason:
            direct_unresolved.append(rank)
            reasons.append(f"{rank}: {reason}")

    unresolved_ranks: tuple[str, ...]
    if direct_unresolved:
        first_unresolved_index = min(SILVA_RANKS.index(rank) for rank in direct_unresolved)
        unresolved_ranks = SILVA_RANKS[first_unresolved_index:]
        lowest_reliable_rank = (
            SILVA_RANKS[first_unresolved_index - 1] if first_unresolved_index > 0 else ""
        )
    else:
        unresolved_ranks = ()
        lowest_reliable_rank = SILVA_RANKS[-1]

    return UnresolvedClassification(
        lowest_reliable_rank=lowest_reliable_rank,
        unresolved_ranks=unresolved_ranks,
        unresolved_reason="; ".join(reasons),
    )


def read_silva_records(path: str | Path) -> list[SilvaRecord]:
    """Read all valid-domain SILVA FASTA records."""
    records, _ = _read_silva_records_with_rejections(path)
    return records


def _read_silva_records_with_rejections(
    path: str | Path,
    metadata: dict[str, TypeStrainMetadata] | None = None,
) -> tuple[list[SilvaRecord], list[dict[str, str]]]:
    records: list[SilvaRecord] = []
    rejected_rows: list[dict[str, str]] = []
    type_metadata = metadata or {}

    for fasta_record in read_fasta(path):
        seq_id, taxonomy = parse_silva_header(fasta_record.header)
        if seq_id != fasta_record.seq_id:
            raise ValueError(
                f"Header parser ID mismatch for {fasta_record.header!r}: "
                f"{seq_id!r} != {fasta_record.seq_id!r}"
            )
        reject_reason = _silva_record_reject_reason(fasta_record, taxonomy, type_metadata)
        if reject_reason:
            rejected_rows.append(
                {
                    "seq_id": seq_id,
                    "original_taxonomy": taxonomy.original_taxonomy,
                    "reject_reason": reject_reason,
                }
            )
            continue
        records.append(
            SilvaRecord(
                fasta_record=fasta_record,
                taxonomy=taxonomy,
                classification=classify_unresolved(taxonomy),
            )
        )

    return records, rejected_rows


def initialize_silva_build(
    silva_fasta: str | Path,
    outdir: str | Path,
    type_strain_metadata: str | Path,
    gtdb_taxonomies: Iterable[str | Path],
    threads: int = 1,
) -> dict[str, int]:
    """Initialize an autotax2 build from SILVA FASTA records."""
    output_dir = Path(outdir)
    registry_dir = output_dir / "registry"
    silva_dir = output_dir / "silva"
    logs_dir = output_dir / "logs"
    for directory in (registry_dir, silva_dir, logs_dir):
        directory.mkdir(parents=True, exist_ok=True)

    gtdb_paths = tuple(Path(path) for path in gtdb_taxonomies)
    if len(gtdb_paths) < 2:
        raise ValueError(
            "GTDB ar53 and bac120 taxonomy inputs are required for SILVA initialization."
        )
    metadata = read_type_strain_metadata(type_strain_metadata) if type_strain_metadata else {}
    records, rejected_rows = _read_silva_records_with_rejections(silva_fasta, metadata)
    records, legal_catalog = _curate_silva_records_with_legal_names(records, metadata, gtdb_paths)
    named_records = [record for record in records if not record.is_unresolved]
    unresolved_records = [record for record in records if record.is_unresolved]
    taxon_rows = _taxon_node_rows(named_records)
    species_taxon_by_taxonomy = _species_taxon_id_by_taxonomy(taxon_rows)

    write_fasta(
        [record.fasta_record for record in named_records],
        silva_dir / "silva_named_backbone.fa",
    )
    write_fasta(
        [record.fasta_record for record in unresolved_records],
        silva_dir / "silva_unresolved.fa",
    )
    _write_tsv(
        _named_taxonomy_rows(named_records),
        silva_dir / "silva_named_backbone.tax.tsv",
        ["seq_id", "taxonomy_7rank", *SILVA_RANKS],
    )
    _write_tsv(
        _unresolved_rows(unresolved_records),
        silva_dir / "silva_unresolved.tsv",
        [
            "seq_id",
            "original_taxonomy",
            *SILVA_RANKS,
            "lowest_reliable_rank",
            "unresolved_ranks",
            "unresolved_reason",
        ],
    )
    _write_tsv(
        rejected_rows,
        silva_dir / "silva_rejected.tsv",
        SILVA_REJECTED_FIELDS,
    )
    _write_tsv(
        taxon_rows,
        registry_dir / "taxon_nodes.tsv",
        [
            "taxon_id",
            "rank",
            "name",
            "parent_taxon_id",
            "taxonomy_path",
            "protected",
            "is_silva_named",
            "source",
        ],
    )
    _write_tsv(
        _sequence_registry_rows(records, metadata, species_taxon_by_taxonomy),
        registry_dir / "sequence_registry.tsv",
        [
            "seq_id",
            "taxon_id",
            "source",
            "sequence_md5",
            "sequence_length",
            "protected",
            "is_silva_named",
            "is_silva_unresolved",
            "original_taxonomy",
            "taxonomy_7rank",
            "is_type_strain",
            "type_species_name",
            "type_strain_id",
            "type_source",
            "type_evidence",
        ],
    )
    _write_tsv(
        _legal_name_catalog_rows(legal_catalog),
        registry_dir / "legal_name_catalog.tsv",
        LEGAL_NAME_CATALOG_FIELDS,
    )
    _write_tsv(
        [
            {
                "dataset_id": SILVA_SOURCE,
                "source": SILVA_SOURCE,
                "description": "SILVA 138.2 SSURef NR99 taxonomy FASTA",
                "protected": "true",
            }
        ],
        registry_dir / "dataset_registry.tsv",
        ["dataset_id", "source", "description", "protected"],
    )
    _write_tsv(
        [
            {
                "build_dir": str(output_dir),
                "silva_fasta": str(silva_fasta),
                "type_strain_metadata": str(type_strain_metadata),
                "gtdb_taxonomy": ";".join(str(path) for path in gtdb_paths),
                "threads": str(threads),
                "metadata_contract": "SILVA full_metadata + GTDB r232 taxonomy",
                "strict_legal_names": "true",
                "organelle_policy": "standardize_organelle_lineages_as_7rank",
            }
        ],
        registry_dir / "build_manifest.tsv",
        BUILD_MANIFEST_FIELDS,
    )
    _write_tsv(
        _name_index_rows(taxon_rows),
        registry_dir / "name_index.tsv",
        ["name", "rank", "taxon_id", "protected", "source"],
    )
    _write_tsv(
        [],
        registry_dir / "cluster_to_taxon.tsv",
        ["cluster_key", "taxon_id", "rank", "name", "status", "source_prefix"],
    )
    _write_tsv(
        _representative_registry_rows(named_records, metadata, species_taxon_by_taxonomy),
        registry_dir / "representative_registry.tsv",
        [
            "representative_seq_id",
            "taxon_id",
            "dataset",
            "source_category",
            "status",
            "protected",
            "source",
            "is_type_strain",
            "representative_reason",
        ],
    )
    _write_placeholder_counters(registry_dir, {SILVA_PLACEHOLDER_PREFIX: {rank: 1 for rank in PLACEHOLDER_RANKS}})
    _write_tsv(
        [
            {
                "taxon_id": row["taxon_id"],
                "rank": row["rank"],
                "name": row["name"],
                "parent_taxon_id": row["parent_taxon_id"],
            }
            for row in taxon_rows
        ],
        registry_dir / "protected_taxa_snapshot.tsv",
        ["taxon_id", "rank", "name", "parent_taxon_id"],
    )
    return {
        "records": len(records),
        "rejected": len(rejected_rows),
        "named": len(named_records),
        "unresolved": len(unresolved_records),
        "threads": threads,
    }


def read_type_strain_metadata(path: str | Path) -> dict[str, TypeStrainMetadata]:
    """Read SILVA official full_metadata keyed by SILVA sequence ID."""
    with _open_metadata_text(path) as handle:
        reader = csv.DictReader(handle, delimiter=chr(9))
        if reader.fieldnames is None:
            raise ValueError("Type strain metadata must have a header row.")
        _validate_silva_full_metadata_header(reader.fieldnames)
        return _read_silva_full_metadata(reader)


def _curate_silva_records_with_legal_names(
    records: Iterable[SilvaRecord],
    metadata: dict[str, TypeStrainMetadata],
    gtdb_taxonomies: Iterable[Path],
) -> tuple[list[SilvaRecord], LegalNameCatalog]:
    record_list = list(records)
    catalog = _legal_name_catalog(record_list, metadata, gtdb_taxonomies)
    if not catalog.has_names:
        return record_list, catalog

    curated_records: list[SilvaRecord] = []
    for record in record_list:
        taxonomy = _curate_taxonomy_with_legal_names(record.taxonomy, catalog)
        organelle_label = _organelle_label(taxonomy)
        if organelle_label:
            taxonomy = _standardize_organelle_taxonomy(taxonomy, organelle_label)
            classification = _resolved_classification()
        else:
            classification = classify_unresolved(taxonomy, catalog)
        curated_records.append(
            SilvaRecord(
                fasta_record=record.fasta_record,
                taxonomy=taxonomy,
                classification=classification,
            )
        )
    return curated_records, catalog


def _legal_name_catalog(
    records: Iterable[SilvaRecord],
    metadata: dict[str, TypeStrainMetadata],
    gtdb_taxonomies: Iterable[Path],
) -> LegalNameCatalog:
    (
        names_by_rank,
        canonical_species_by_name,
        sources_by_rank_name,
        display_names_by_rank_name,
    ) = _type_strain_name_catalog(records, metadata)
    _merge_gtdb_names(
        names_by_rank,
        canonical_species_by_name,
        sources_by_rank_name,
        display_names_by_rank_name,
        gtdb_taxonomies,
    )
    return LegalNameCatalog(
        names_by_rank={rank: frozenset(values) for rank, values in names_by_rank.items()},
        canonical_species_by_name=canonical_species_by_name,
        sources_by_rank_name={
            key: frozenset(values) for key, values in sources_by_rank_name.items()
        },
        display_names_by_rank_name=display_names_by_rank_name,
    )


def _type_strain_name_catalog(
    records: Iterable[SilvaRecord],
    metadata: dict[str, TypeStrainMetadata],
) -> tuple[
    dict[str, set[str]],
    dict[str, str],
    dict[tuple[str, str], set[str]],
    dict[tuple[str, str], str],
]:
    names_by_rank: dict[str, set[str]] = {rank: set() for rank in SILVA_RANKS}
    canonical_species_by_name: dict[str, str] = {}
    sources_by_rank_name: dict[tuple[str, str], set[str]] = {}
    display_names_by_rank_name: dict[tuple[str, str], str] = {}

    for record in records:
        type_metadata = _metadata_for_seq_id(record.seq_id, metadata)
        if not _is_type_strain_metadata(type_metadata):
            continue

        species_name = _canonical_species_name(type_metadata.species_name)
        if not species_name:
            species_name = _canonical_species_name(record.taxonomy.species)
        if not species_name:
            continue

        for rank, value in zip(SILVA_RANKS, record.taxonomy.values, strict=True):
            if rank == "species":
                continue
            clean_value = value.strip()
            if clean_value and not _unresolved_reason(rank, clean_value):
                _add_legal_name(
                    names_by_rank,
                    sources_by_rank_name,
                    display_names_by_rank_name,
                    rank,
                    clean_value,
                    "SILVA full_metadata type material",
                )

        species_key = _normalize_taxon_value(species_name)
        _add_legal_name(
            names_by_rank,
            sources_by_rank_name,
            display_names_by_rank_name,
            "species",
            species_name,
            "SILVA full_metadata type material",
        )
        canonical_species_by_name[species_key] = species_name

    return names_by_rank, canonical_species_by_name, sources_by_rank_name, display_names_by_rank_name


def _merge_gtdb_names(
    names_by_rank: dict[str, set[str]],
    canonical_species_by_name: dict[str, str],
    sources_by_rank_name: dict[tuple[str, str], set[str]],
    display_names_by_rank_name: dict[tuple[str, str], str],
    paths: Iterable[Path],
) -> None:
    for path in paths:
        for rank, name in _iter_gtdb_taxon_names(path):
            normalized = _normalize_taxon_value(name)
            if not normalized or _is_reserved_autotax_placeholder_name(rank, name):
                continue
            _add_legal_name(
                names_by_rank,
                sources_by_rank_name,
                display_names_by_rank_name,
                rank,
                name,
                "GTDB r232 taxonomy",
            )
            if rank == "species":
                canonical = _canonical_species_name(name)
                if canonical:
                    canonical_species_by_name[_normalize_taxon_value(canonical)] = canonical


def _add_legal_name(
    names_by_rank: dict[str, set[str]],
    sources_by_rank_name: dict[tuple[str, str], set[str]],
    display_names_by_rank_name: dict[tuple[str, str], str],
    rank: str,
    name: str,
    source: str,
) -> None:
    normalized = _normalize_taxon_value(name)
    if not normalized:
        return
    names_by_rank.setdefault(rank, set()).add(normalized)
    sources_by_rank_name.setdefault((rank, normalized), set()).add(source)
    display_names_by_rank_name.setdefault((rank, normalized), name.strip())


def _iter_gtdb_taxon_names(path: Path) -> Iterable[tuple[str, str]]:
    with _open_metadata_text(path) as handle:
        reader = csv.reader(handle, delimiter=chr(9))
        for row in reader:
            if len(row) < 2:
                continue
            taxonomy = row[1].strip()
            if "__" not in taxonomy:
                continue
            for part in taxonomy.split(";"):
                rank, name = _parse_gtdb_taxon(part)
                if rank and name:
                    yield rank, name


def _parse_gtdb_taxon(value: str) -> tuple[str, str]:
    clean = value.strip()
    if len(clean) < 4 or clean[1:3] != "__":
        return "", ""
    rank = GTDB_PREFIX_TO_RANK.get(clean[0].lower(), "")
    name = clean[3:].strip()
    return rank, name


def _is_reserved_autotax_placeholder_name(rank: str, name: str) -> bool:
    rank_prefix = RANK_PREFIXES.get(rank, "")
    if not rank_prefix:
        return False
    return re.fullmatch(rf"[A-Za-z][A-Za-z0-9]*{rank_prefix}[0-9]{{6}}", name) is not None


def _curate_taxonomy_with_legal_names(
    taxonomy: SilvaTaxonomy,
    catalog: LegalNameCatalog,
) -> SilvaTaxonomy:
    values = list(taxonomy.values)
    species = catalog.canonical_species(taxonomy.species)
    if species:
        values[SILVA_RANKS.index("species")] = species
    return _taxonomy_with_values(taxonomy, values)


def _taxonomy_with_values(taxonomy: SilvaTaxonomy, values: list[str]) -> SilvaTaxonomy:
    return SilvaTaxonomy(
        domain=values[0],
        phylum=values[1],
        class_name=values[2],
        order=values[3],
        family=values[4],
        genus=values[5],
        species=values[6],
        original_taxonomy=taxonomy.original_taxonomy,
    )


def _organelle_label(taxonomy: SilvaTaxonomy) -> str:
    normalized_values = [_normalize_taxon_value(value) for value in taxonomy.values]
    for normalized_value in normalized_values:
        for marker, label in ORGANELLE_MARKERS.items():
            if marker in normalized_value:
                return label
    return ""


def _standardize_organelle_taxonomy(taxonomy: SilvaTaxonomy, organelle_label: str) -> SilvaTaxonomy:
    values = list(taxonomy.values)
    marker_index = _organelle_marker_index(taxonomy)
    if marker_index is None:
        return taxonomy
    for index in range(marker_index, len(values)):
        values[index] = organelle_label
    return _taxonomy_with_values(taxonomy, values)


def _organelle_marker_index(taxonomy: SilvaTaxonomy) -> int | None:
    for index, value in enumerate(taxonomy.values):
        normalized_value = _normalize_taxon_value(value)
        if any(marker in normalized_value for marker in ORGANELLE_MARKERS):
            return index
    return None


def _resolved_classification() -> UnresolvedClassification:
    return UnresolvedClassification(
        lowest_reliable_rank=SILVA_RANKS[-1],
        unresolved_ranks=(),
        unresolved_reason="",
    )


def _read_silva_full_metadata(reader: csv.DictReader) -> dict[str, TypeStrainMetadata]:
    field_map = _metadata_field_map(reader.fieldnames or [])
    metadata: dict[str, TypeStrainMetadata] = {}
    for row in reader:
        seq_id = _metadata_value(
            row,
            field_map,
            "seqid",
            "sequenceid",
            "primaryaccession",
            "accession",
            "accessionnumber",
            "acc",
            "id",
        )
        if not seq_id:
            continue
        type_field, type_value = _type_material_evidence(row)
        if not _is_type_strain_value(type_value):
            continue
        metadata[seq_id] = TypeStrainMetadata(
            seq_id=seq_id,
            is_type_strain="true",
            species_name=_metadata_species_name(row, field_map),
            strain_id=_metadata_value(
                row,
                field_map,
                "strainid",
                "strain",
                "isolate",
                "culturescollection",
                "culturecollection",
            ),
            source="SILVA full_metadata",
            evidence=f"{type_field}={type_value}" if type_field else type_value,
        )
    return metadata


def _validate_silva_full_metadata_header(fieldnames: Iterable[str]) -> None:
    field_map = _metadata_field_map(fieldnames)
    has_accession = any(
        candidate in field_map
        for candidate in ("primaryaccession", "accession", "accessionnumber", "acc")
    )
    has_type_material = any("type" in normalized for normalized in field_map)
    if not has_accession or not has_type_material:
        raise ValueError(
            "type-strain metadata must be an official SILVA full_metadata TSV/TSV.gz "
            "with accession and type-material fields."
        )


def _metadata_field_map(fieldnames: Iterable[str]) -> dict[str, str]:
    return {_normalize_metadata_field(field): field for field in fieldnames}


def _normalize_metadata_field(field: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", field.lower())


def _metadata_value(
    row: dict[str, str],
    field_map: dict[str, str],
    *candidates: str,
) -> str:
    for candidate in candidates:
        field = field_map.get(candidate)
        if field:
            value = (row.get(field) or "").strip()
            if value:
                return value
    return ""


def _metadata_species_name(row: dict[str, str], field_map: dict[str, str]) -> str:
    direct = _metadata_value(
        row,
        field_map,
        "speciesname",
        "organismname",
        "organism",
        "species",
        "taxspecies",
    )
    if direct:
        return direct
    for field, value in row.items():
        normalized = _normalize_metadata_field(field)
        if "tax" not in normalized:
            continue
        parts = [part.strip() for part in (value or "").split(";") if part.strip()]
        if parts:
            return parts[-1]
    return ""


def _type_material_evidence(row: dict[str, str]) -> tuple[str, str]:
    preferred_fields = (
        "is_type_strain",
        "type_strain",
        "type_material",
        "typematerial",
        "type",
    )
    normalized_lookup = {_normalize_metadata_field(field): field for field in row}
    for preferred in preferred_fields:
        field = normalized_lookup.get(_normalize_metadata_field(preferred))
        if field:
            value = (row.get(field) or "").strip()
            if value:
                return field, value
    for field, value in row.items():
        if "type" not in _normalize_metadata_field(field):
            continue
        clean_value = (value or "").strip()
        if clean_value:
            return field, clean_value
    return "", ""


def _is_type_strain_value(value: str) -> bool:
    normalized = value.strip().lower()
    if not normalized:
        return False
    if normalized in {"false", "0", "no", "n", "none", "na", "n/a"}:
        return False
    if "non-type" in normalized or "not type" in normalized or "not a type" in normalized:
        return False
    if normalized in {"true", "1", "yes", "y", "t"}:
        return True
    compact = re.sub(r"[^a-z0-9]+", "", normalized)
    return (
        "typestrain" in compact
        or "typematerial" in compact
        or "typespecies" in compact
    )


def _open_metadata_text(path: str | Path) -> TextIO:
    metadata_path = Path(path)
    if metadata_path.suffix == ".gz":
        return gzip.open(metadata_path, "rt", encoding="utf-8", newline="")
    return metadata_path.open("r", encoding="utf-8", newline="")


def resolve_silva_unresolved_build(
    build: str | Path,
    threads: int = 4,
    species_id: float = DEFAULT_SPECIES_ID,
    genus_id: float = DEFAULT_GENUS_ID,
    family_id: float = DEFAULT_FAMILY_ID,
    order_id: float = DEFAULT_ORDER_ID,
    class_id: float = DEFAULT_CLASS_ID,
    phylum_id: float = DEFAULT_PHYLUM_ID,
    floor_id: float | None = None,
    vsearch_bin: str = "vsearch",
    iddef: int = 2,
    dry_run: bool = False,
) -> ResolveSilvaSummary:
    """Resolve SILVA unresolved records into mutable placeholder taxa."""
    if floor_id is not None:
        phylum_id = floor_id

    build_dir = Path(build)
    silva_dir = build_dir / "silva"
    registry_dir = build_dir / "registry"
    cluster_dir = silva_dir / "silva_unresolved_clusters"
    cluster_dir.mkdir(parents=True, exist_ok=True)
    registry_dir.mkdir(parents=True, exist_ok=True)

    unresolved_fasta_path = silva_dir / "silva_unresolved.fa"
    unresolved_tsv_path = silva_dir / "silva_unresolved.tsv"
    records = read_fasta(unresolved_fasta_path) if unresolved_fasta_path.exists() else []
    unresolved_rows = _read_tsv(unresolved_tsv_path) if unresolved_tsv_path.exists() else []
    rows_by_seq_id = {row["seq_id"]: row for row in unresolved_rows}
    records = [record for record in records if record.seq_id in rows_by_seq_id]

    output_paths = _resolve_output_paths(silva_dir, dry_run)
    if not records:
        _write_resolve_outputs(
            taxa_rows=[],
            member_rows=[],
            mapping_rows=[],
            evidence_rows=[],
            resolved_records=[],
            output_paths=output_paths,
        )
        return ResolveSilvaSummary(unresolved_records=0, placeholder_taxa=0, dry_run=dry_run)

    state = _load_placeholder_state(registry_dir)
    allocator = PlaceholderAllocator(
        issued_ids=state["issued_placeholders"],
        deprecated_ids=state["deprecated_placeholders"],
    )
    active_cluster_to_taxon = state["active_cluster_to_taxon"]

    thresholds = {
        "phylum": phylum_id,
        "class": class_id,
        "order": order_id,
        "family": family_id,
        "genus": genus_id,
        "species": species_id,
    }
    states = _initial_resolve_states(records, rows_by_seq_id)
    named_anchors = _load_named_silva_anchors(silva_dir)
    taxa_by_cluster_key: dict[str, dict[str, str]] = {}
    evidence_rows: list[dict[str, str]] = []

    for rank in RESOLVE_RANKS:
        _resolve_rank_with_parent_scope(
            rank=rank,
            states=states,
            named_anchors=named_anchors,
            threshold=thresholds[rank],
            cluster_dir=cluster_dir,
            total_threads=threads,
            vsearch_bin=vsearch_bin,
            iddef=iddef,
            allocator=allocator,
            active_cluster_to_taxon=active_cluster_to_taxon,
            taxa_by_cluster_key=taxa_by_cluster_key,
            evidence_rows=evidence_rows,
        )

    taxa_rows = list(taxa_by_cluster_key.values())
    member_rows = _resolve_member_rows(states)
    mapping_rows = _resolve_mapping_rows(states)
    resolved_records = _resolve_fasta_records(states)
    _write_resolve_outputs(
        taxa_rows=taxa_rows,
        member_rows=member_rows,
        mapping_rows=mapping_rows,
        evidence_rows=evidence_rows,
        resolved_records=resolved_records,
        output_paths=output_paths,
    )

    if not dry_run:
        _update_resolve_registry(
            registry_dir=registry_dir,
            taxa_rows=taxa_rows,
            allocator=allocator,
        )

    return ResolveSilvaSummary(
        unresolved_records=len(records),
        placeholder_taxa=len(taxa_rows),
        dry_run=dry_run,
    )


def _ensure_uc(
    fasta_path: Path,
    uc_path: Path,
    identity: float,
    records: list[FastaRecord],
    threads: int,
    vsearch_bin: str,
    iddef: int,
) -> None:
    if uc_path.exists():
        return
    try:
        run_vsearch_cluster(
            fasta_path=fasta_path,
            uc_path=uc_path,
            identity=identity,
            threads=threads,
            vsearch_bin=vsearch_bin,
            iddef=iddef,
        )
    except (FileNotFoundError, OSError, subprocess.CalledProcessError):
        write_singleton_uc([record.seq_id for record in records], uc_path)


def _assignments_by_seq_id(
    uc_path: Path,
    records: list[FastaRecord],
    level: str,
) -> dict[str, ClusterAssignment]:
    raw_assignments = read_uc_assignments(uc_path)
    assignments: dict[str, ClusterAssignment] = {
        assignment.seq_id: ClusterAssignment(
            cluster_key=f"{level}:{assignment.cluster_id}",
            cluster_id=assignment.cluster_id,
            representative_seq_id=assignment.representative_seq_id,
        )
        for assignment in raw_assignments
    }
    used_cluster_ids = {
        int(assignment.cluster_id)
        for assignment in assignments.values()
        if assignment.cluster_id.isdigit()
    }
    next_cluster_id = (max(used_cluster_ids) + 1) if used_cluster_ids else 0
    for record in records:
        if record.seq_id in assignments:
            continue
        cluster_id = str(next_cluster_id)
        next_cluster_id += 1
        assignments[record.seq_id] = ClusterAssignment(
            cluster_key=f"{level}:{cluster_id}",
            cluster_id=cluster_id,
            representative_seq_id=record.seq_id,
        )
    return assignments


def _cluster_counts(assignments: dict[str, ClusterAssignment]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for assignment in assignments.values():
        counts[assignment.cluster_key] = counts.get(assignment.cluster_key, 0) + 1
    return counts


def _initial_resolve_states(
    records: list[FastaRecord],
    rows_by_seq_id: dict[str, dict[str, str]],
) -> list[ResolveRecordState]:
    states: list[ResolveRecordState] = []
    for record in records:
        source_row = rows_by_seq_id[record.seq_id]
        values = _row_taxonomy_values(source_row)
        resolved_labels = {"domain": _label_for_taxon("domain", values[0])}
        states.append(
            ResolveRecordState(
                record=record,
                source_row=source_row,
                unresolved_ranks=_split_unresolved_ranks(source_row.get("unresolved_ranks", "")),
                resolved_labels=resolved_labels,
            )
        )
    return states


def _load_named_silva_anchors(silva_dir: Path) -> list[ResolveNamedAnchorRecord]:
    fasta_by_id = {
        record.seq_id: record
        for record in read_fasta(silva_dir / "silva_named_backbone.fa")
    } if (silva_dir / "silva_named_backbone.fa").exists() else {}
    anchors: list[ResolveNamedAnchorRecord] = []
    for row in _read_tsv(silva_dir / "silva_named_backbone.tax.tsv"):
        seq_id = row.get("seq_id", "")
        record = fasta_by_id.get(seq_id)
        if record is None:
            continue
        anchors.append(
            ResolveNamedAnchorRecord(
                seq_id=seq_id,
                values=_row_taxonomy_values(row),
                sequence=record.sequence,
            )
        )
    return anchors


def _resolve_rank_with_parent_scope(
    rank: str,
    states: list[ResolveRecordState],
    named_anchors: list[ResolveNamedAnchorRecord],
    threshold: float,
    cluster_dir: Path,
    total_threads: int,
    vsearch_bin: str,
    iddef: int,
    allocator: PlaceholderAllocator,
    active_cluster_to_taxon: dict[str, dict[str, str]],
    taxa_by_cluster_key: dict[str, dict[str, str]],
    evidence_rows: list[dict[str, str]],
) -> None:
    parent_groups = _resolve_parent_groups(rank, states)
    parent_items = list(parent_groups.items())
    job_sizes = [
        len([state for state in group_states if _state_needs_rank_resolution(state, rank)])
        for _, group_states in parent_items
    ]
    job_threads_by_parent = dict(
        zip(
            [parent_key for parent_key, _ in parent_items],
            _allocate_parent_job_threads(total_threads, job_sizes),
            strict=True,
        )
    )
    jobs: list[ResolveParentJob] = []

    for parent_key, group_states in parent_items:
        candidates = [
            state
            for state in group_states
            if _state_needs_rank_resolution(state, rank)
        ]
        anchors = _same_parent_anchors(rank, parent_key, group_states, named_anchors)
        job_threads = job_threads_by_parent.get(parent_key, 1)

        for state in group_states:
            if state.blocked:
                _set_unplaced_rank(state, rank)
                evidence_rows.append(
                    _resolve_evidence_row(
                        state=state,
                        rank=rank,
                        threshold=threshold,
                        parent_key=parent_key,
                        input_rank_status="blocked",
                        anchor_count=len(anchors),
                        decision="unplaced",
                        output_taxon=state.resolved_labels.get(rank, ""),
                        reason="upstream_rank_ambiguous",
                        job_size=len(candidates),
                        job_threads=job_threads,
                    )
                )
            elif not _state_needs_rank_resolution(state, rank):
                label = _label_for_taxon(rank, _row_taxonomy_value(state.source_row, rank))
                state.resolved_labels[rank] = label
                evidence_rows.append(
                    _resolve_evidence_row(
                        state=state,
                        rank=rank,
                        threshold=threshold,
                        parent_key=parent_key,
                        input_rank_status="known",
                        anchor_count=len(anchors),
                        decision="keep_known",
                        output_taxon=label,
                        reason="rank_already_reliable_in_silva",
                        job_size=len(candidates),
                        job_threads=job_threads,
                    )
                )

        if not candidates:
            continue

        jobs.append(
            ResolveParentJob(
                rank=rank,
                parent_key=parent_key,
                candidates=candidates,
                anchors=anchors,
                threshold=threshold,
                cluster_dir=cluster_dir,
                job_threads=job_threads,
                vsearch_bin=vsearch_bin,
                iddef=iddef,
            )
        )

    for result in _run_parent_resolve_jobs(jobs, total_threads=total_threads):
        _apply_candidate_group_result(
            rank=rank,
            parent_key=result.parent_key,
            candidates=result.candidates,
            anchors=result.anchors,
            threshold=result.threshold,
            cluster_members=result.cluster_members,
            job_id_to_anchor=result.job_id_to_anchor,
            candidate_anchor_identities=result.candidate_anchor_identities,
            job_threads=result.job_threads,
            allocator=allocator,
            active_cluster_to_taxon=active_cluster_to_taxon,
            taxa_by_cluster_key=taxa_by_cluster_key,
            evidence_rows=evidence_rows,
        )


def _state_needs_rank_resolution(state: ResolveRecordState, rank: str) -> bool:
    return not state.blocked and rank in state.unresolved_ranks and rank in PLACEHOLDER_RANKS


def _resolve_parent_groups(
    rank: str,
    states: list[ResolveRecordState],
) -> dict[tuple[str, ...], list[ResolveRecordState]]:
    groups: dict[tuple[str, ...], list[ResolveRecordState]] = {}
    for state in states:
        parent_key = _state_parent_key(state, rank)
        groups.setdefault(parent_key, []).append(state)
    return groups


def _state_parent_key(state: ResolveRecordState, rank: str) -> tuple[str, ...]:
    rank_index = SILVA_RANKS.index(rank)
    labels: list[str] = []
    for parent_rank in SILVA_RANKS[:rank_index]:
        label = state.resolved_labels.get(parent_rank)
        if not label:
            label = _label_for_taxon(parent_rank, _row_taxonomy_value(state.source_row, parent_rank))
            state.resolved_labels[parent_rank] = label
        labels.append(label)
    return tuple(labels)


def _same_parent_anchors(
    rank: str,
    parent_key: tuple[str, ...],
    group_states: list[ResolveRecordState],
    named_anchors: list[ResolveNamedAnchorRecord],
) -> list[ResolveAnchor]:
    anchors: list[ResolveAnchor] = []
    seen: set[tuple[str, str]] = set()

    for named in named_anchors:
        if _values_parent_key(named.values, rank) != parent_key:
            continue
        anchor = _make_resolve_anchor(
            rank=rank,
            taxon_name=named.values[SILVA_RANKS.index(rank)],
            seq_id=named.seq_id,
            sequence=named.sequence,
            source="named_silva",
            index=len(anchors) + 1,
        )
        key = (anchor.taxon_label, anchor.seq_id)
        if key not in seen:
            anchors.append(anchor)
            seen.add(key)

    for state in group_states:
        if rank in state.unresolved_ranks:
            continue
        value = _row_taxonomy_value(state.source_row, rank)
        if _unresolved_reason(rank, value):
            continue
        anchor = _make_resolve_anchor(
            rank=rank,
            taxon_name=value,
            seq_id=state.record.seq_id,
            sequence=state.record.sequence,
            source="silva_unresolved_known_rank",
            index=len(anchors) + 1,
        )
        key = (anchor.taxon_label, anchor.seq_id)
        if key not in seen:
            anchors.append(anchor)
            seen.add(key)

    return anchors


def _make_resolve_anchor(
    rank: str,
    taxon_name: str,
    seq_id: str,
    sequence: str,
    source: str,
    index: int,
) -> ResolveAnchor:
    return ResolveAnchor(
        job_id=f"ANCHOR_{index:06d}",
        seq_id=seq_id,
        rank=rank,
        taxon_name=taxon_name,
        taxon_label=_label_for_taxon(rank, taxon_name),
        sequence=sequence,
        source=source,
    )


def _values_parent_key(
    values: tuple[str, str, str, str, str, str, str],
    rank: str,
) -> tuple[str, ...]:
    rank_index = SILVA_RANKS.index(rank)
    return tuple(
        _label_for_taxon(parent_rank, values[index])
        for index, parent_rank in enumerate(SILVA_RANKS[:rank_index])
    )


def _run_parent_resolve_jobs(
    jobs: list[ResolveParentJob],
    total_threads: int,
) -> list[ResolveParentJobResult]:
    if not jobs:
        return []
    results: list[ResolveParentJobResult] = []
    for batch in _parent_job_batches(jobs, total_threads):
        with ThreadPoolExecutor(max_workers=len(batch)) as executor:
            results.extend(executor.map(_run_parent_resolve_job, batch))
    return results


def _parent_job_batches(
    jobs: list[ResolveParentJob],
    total_threads: int,
) -> list[list[ResolveParentJob]]:
    capacity = max(1, total_threads)
    batches: list[list[ResolveParentJob]] = []
    batch: list[ResolveParentJob] = []
    used_threads = 0
    for job in jobs:
        job_threads = max(1, min(capacity, job.job_threads))
        if batch and used_threads + job_threads > capacity:
            batches.append(batch)
            batch = []
            used_threads = 0
        batch.append(job)
        used_threads += job_threads
    if batch:
        batches.append(batch)
    return batches


def _run_parent_resolve_job(job: ResolveParentJob) -> ResolveParentJobResult:
    job_fasta, uc_path, job_id_to_anchor = _write_parent_job_fasta(
        rank=job.rank,
        parent_key=job.parent_key,
        candidates=job.candidates,
        anchors=job.anchors,
        threshold=job.threshold,
        cluster_dir=job.cluster_dir,
    )
    job_records = read_fasta(job_fasta)
    _ensure_uc(
        fasta_path=job_fasta,
        uc_path=uc_path,
        identity=job.threshold,
        records=job_records,
        threads=job.job_threads,
        vsearch_bin=job.vsearch_bin,
        iddef=job.iddef,
    )
    return ResolveParentJobResult(
        rank=job.rank,
        parent_key=job.parent_key,
        candidates=job.candidates,
        anchors=job.anchors,
        threshold=job.threshold,
        cluster_members=_cluster_members_by_id(uc_path, job_records),
        job_id_to_anchor=job_id_to_anchor,
        candidate_anchor_identities=_candidate_anchor_identities(
            uc_path=uc_path,
            candidate_ids={state.record.seq_id for state in job.candidates},
            job_id_to_anchor=job_id_to_anchor,
        ),
        job_size=len(job.candidates),
        job_threads=job.job_threads,
    )


def _candidate_anchor_identities(
    uc_path: Path,
    candidate_ids: set[str],
    job_id_to_anchor: dict[str, ResolveAnchor],
) -> dict[tuple[str, str], str]:
    identities: dict[tuple[str, str], str] = {}
    for record in parse_uc_records(uc_path):
        if record.record_type != "H":
            continue
        query_id = record.query_label
        target_id = record.target_label
        identity = _uc_percent_identity_to_fraction(record.percent_identity)
        if not identity:
            continue
        if query_id in candidate_ids and target_id in job_id_to_anchor:
            key = (query_id, job_id_to_anchor[target_id].taxon_label)
        elif target_id in candidate_ids and query_id in job_id_to_anchor:
            key = (target_id, job_id_to_anchor[query_id].taxon_label)
        else:
            continue
        identities[key] = _max_identity_text(identities.get(key, ""), identity)
    return identities


def _uc_percent_identity_to_fraction(percent_identity: str) -> str:
    if not percent_identity or percent_identity == "*":
        return ""
    try:
        value = float(percent_identity)
    except ValueError:
        return ""
    if value > 1:
        value = value / 100
    return f"{value:.3f}"


def _max_identity_text(left: str, right: str) -> str:
    if not left:
        return right
    try:
        return left if float(left) >= float(right) else right
    except ValueError:
        return left or right


def _best_anchor_evidence(
    seq_id: str,
    anchor_labels: list[str],
    candidate_anchor_identities: dict[tuple[str, str], str],
    threshold: float,
) -> tuple[str, str]:
    if not anchor_labels:
        return "", ""
    best_label = anchor_labels[0]
    best_identity = ""
    for label in anchor_labels:
        identity = candidate_anchor_identities.get((seq_id, label), "")
        if identity and _max_identity_text(best_identity, identity) == identity:
            best_label = label
            best_identity = identity
    return best_label, best_identity or f">={threshold:.3f}"


def _apply_candidate_group_result(
    rank: str,
    parent_key: tuple[str, ...],
    candidates: list[ResolveRecordState],
    anchors: list[ResolveAnchor],
    threshold: float,
    cluster_members: dict[str, list[str]],
    job_id_to_anchor: dict[str, ResolveAnchor],
    candidate_anchor_identities: dict[tuple[str, str], str],
    job_threads: int,
    allocator: PlaceholderAllocator,
    active_cluster_to_taxon: dict[str, dict[str, str]],
    taxa_by_cluster_key: dict[str, dict[str, str]],
    evidence_rows: list[dict[str, str]],
) -> None:
    candidate_by_job_id = {state.record.seq_id: state for state in candidates}

    for cluster_id, member_ids in sorted(cluster_members.items()):
        cluster_anchors = [
            job_id_to_anchor[job_id]
            for job_id in member_ids
            if job_id in job_id_to_anchor
        ]
        cluster_candidates = [
            candidate_by_job_id[job_id]
            for job_id in member_ids
            if job_id in candidate_by_job_id
        ]
        if not cluster_candidates:
            continue
        distinct_anchor_labels = sorted({anchor.taxon_label for anchor in cluster_anchors})

        if len(distinct_anchor_labels) == 1:
            label = distinct_anchor_labels[0]
            for state in cluster_candidates:
                best_anchor_taxon, best_anchor_identity = _best_anchor_evidence(
                    seq_id=state.record.seq_id,
                    anchor_labels=distinct_anchor_labels,
                    candidate_anchor_identities=candidate_anchor_identities,
                    threshold=threshold,
                )
                state.resolved_labels[rank] = label
                state.rank_cluster_keys[rank] = _anchor_cluster_key(rank, parent_key, cluster_id)
                evidence_rows.append(
                    _resolve_evidence_row(
                        state=state,
                        rank=rank,
                        threshold=threshold,
                        parent_key=parent_key,
                        input_rank_status="unresolved",
                        anchor_count=len(anchors),
                        best_anchor_taxon=best_anchor_taxon,
                        best_anchor_identity=best_anchor_identity,
                        passing_anchor_count=1,
                        residual_cluster_key="",
                        residual_cluster_size=0,
                        decision="assign_existing",
                        output_taxon=label,
                        reason="unique_same_parent_anchor_cluster",
                        job_size=len(candidates),
                        job_threads=job_threads,
                    )
                )
        elif len(distinct_anchor_labels) > 1:
            for state in cluster_candidates:
                best_anchor_taxon, best_anchor_identity = _best_anchor_evidence(
                    seq_id=state.record.seq_id,
                    anchor_labels=distinct_anchor_labels,
                    candidate_anchor_identities=candidate_anchor_identities,
                    threshold=threshold,
                )
                state.blocked = True
                label = _ambiguous_rank_label(rank)
                state.resolved_labels[rank] = label
                state.warnings.add(f"ambiguous_{rank}_same_parent_anchor")
                evidence_rows.append(
                    _resolve_evidence_row(
                        state=state,
                        rank=rank,
                        threshold=threshold,
                        parent_key=parent_key,
                        input_rank_status="unresolved",
                        anchor_count=len(anchors),
                        best_anchor_taxon=best_anchor_taxon,
                        best_anchor_identity=best_anchor_identity,
                        passing_anchor_count=len(distinct_anchor_labels),
                        competing_anchor_taxa=";".join(distinct_anchor_labels),
                        decision="ambiguous",
                        output_taxon=label,
                        reason="multiple_same_parent_anchor_clusters",
                        job_size=len(candidates),
                        job_threads=job_threads,
                    )
                )
        else:
            cluster_key = _parent_scoped_cluster_key(rank, parent_key, cluster_id)
            representative_seq_id = cluster_candidates[0].record.seq_id
            taxon_row = _get_or_create_placeholder_taxon(
                rank=rank,
                cluster_key=cluster_key,
                parent_taxon_id=_parent_taxon_label(parent_key),
                representative_seq_id=representative_seq_id,
                member_count=len(cluster_candidates),
                warning="",
                allocator=allocator,
                active_cluster_to_taxon=active_cluster_to_taxon,
                taxa_by_cluster_key=taxa_by_cluster_key,
            )
            label = taxon_row["name"]
            for state in cluster_candidates:
                state.resolved_labels[rank] = label
                state.placeholder_taxa.append(label)
                state.rank_placeholders[rank] = label
                state.rank_cluster_keys[rank] = cluster_key
                evidence_rows.append(
                    _resolve_evidence_row(
                        state=state,
                        rank=rank,
                        threshold=threshold,
                        parent_key=parent_key,
                        input_rank_status="unresolved",
                        anchor_count=len(anchors),
                        residual_cluster_key=cluster_key,
                        residual_cluster_size=len(cluster_candidates),
                        decision="create_placeholder",
                        output_taxon=label,
                        reason="no_same_parent_anchor_above_threshold",
                        job_size=len(candidates),
                        job_threads=job_threads,
                    )
                )


def _write_parent_job_fasta(
    rank: str,
    parent_key: tuple[str, ...],
    candidates: list[ResolveRecordState],
    anchors: list[ResolveAnchor],
    threshold: float,
    cluster_dir: Path,
) -> tuple[Path, Path, dict[str, ResolveAnchor]]:
    rank_dir = cluster_dir / rank
    rank_dir.mkdir(parents=True, exist_ok=True)
    parent_hash = hashlib.sha1(";".join(parent_key).encode("utf-8")).hexdigest()[:12]
    stem = f"{rank}_{threshold:.3f}_{parent_hash}"
    fasta_path = rank_dir / f"{stem}.fa"
    uc_path = rank_dir / f"{stem}.uc"
    job_records = [
        FastaRecord(seq_id=anchor.job_id, header=anchor.job_id, sequence=anchor.sequence)
        for anchor in anchors
    ]
    job_records.extend(
        FastaRecord(seq_id=state.record.seq_id, header=state.record.seq_id, sequence=state.record.sequence)
        for state in candidates
    )
    write_fasta(job_records, fasta_path)
    return fasta_path, uc_path, {anchor.job_id: anchor for anchor in anchors}


def _cluster_members_by_id(
    uc_path: Path,
    job_records: list[FastaRecord],
) -> dict[str, list[str]]:
    assignments = _assignments_by_seq_id(uc_path, job_records, "parent")
    members: dict[str, list[str]] = {}
    for seq_id, assignment in assignments.items():
        members.setdefault(assignment.cluster_id, []).append(seq_id)
    return members


def _allocate_parent_job_threads(total_threads: int, job_sizes: list[int]) -> list[int]:
    if not job_sizes:
        return []
    total_threads = max(1, total_threads)
    total_weight = sum(max(1, size) for size in job_sizes)
    return [
        max(1, min(total_threads, round(total_threads * max(1, size) / total_weight)))
        for size in job_sizes
    ]


def _parent_scoped_cluster_key(rank: str, parent_key: tuple[str, ...], cluster_id: str) -> str:
    return f"{rank}:parent={_parent_key_text(parent_key)}:cluster={cluster_id}"


def _anchor_cluster_key(rank: str, parent_key: tuple[str, ...], cluster_id: str) -> str:
    return f"{rank}:parent={_parent_key_text(parent_key)}:anchor_cluster={cluster_id}"


def _parent_key_text(parent_key: tuple[str, ...]) -> str:
    return "/".join(parent_key)


def _parent_taxon_label(parent_key: tuple[str, ...]) -> str:
    return parent_key[-1] if parent_key else ""


def _ambiguous_rank_label(rank: str) -> str:
    return f"{RANK_PREFIXES[rank]}__ambiguous_{rank}"


def _set_unplaced_rank(state: ResolveRecordState, rank: str) -> None:
    state.resolved_labels.setdefault(rank, f"{RANK_PREFIXES[rank]}__unplaced_{rank}")


def _resolve_evidence_row(
    state: ResolveRecordState,
    rank: str,
    threshold: float,
    parent_key: tuple[str, ...],
    input_rank_status: str,
    anchor_count: int,
    decision: str,
    output_taxon: str,
    reason: str,
    job_size: int,
    job_threads: int,
    best_anchor_taxon: str = "",
    best_anchor_identity: str = "",
    passing_anchor_count: int = 0,
    competing_anchor_taxa: str = "",
    residual_cluster_key: str = "",
    residual_cluster_size: int = 0,
) -> dict[str, str]:
    rank_index = SILVA_RANKS.index(rank)
    parent_rank = SILVA_RANKS[rank_index - 1] if rank_index > 0 else ""
    return {
        "source_stage": "silva_resolve",
        "dataset": SILVA_SOURCE,
        "seq_id": state.record.seq_id,
        "rank": rank,
        "parent_rank": parent_rank,
        "parent_taxon": _parent_taxon_label(parent_key),
        "input_rank_value": _row_taxonomy_value(state.source_row, rank),
        "input_rank_status": input_rank_status,
        "threshold": f"{threshold:.3f}",
        "candidate_scope": "same-parent only",
        "anchor_count": str(anchor_count),
        "best_anchor_taxon": best_anchor_taxon,
        "best_anchor_identity": best_anchor_identity,
        "passing_anchor_count": str(passing_anchor_count),
        "competing_anchor_taxa": competing_anchor_taxa,
        "residual_cluster_key": residual_cluster_key,
        "residual_cluster_size": str(residual_cluster_size),
        "decision": decision,
        "output_taxon": output_taxon,
        "reason": reason,
        "job_size": str(job_size),
        "job_threads": str(job_threads),
    }


def _resolve_member_rows(states: list[ResolveRecordState]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for state in states:
        resolved_taxonomy = _state_resolved_taxonomy(state)
        rows.append(
            {
                "seq_id": state.record.seq_id,
                "original_taxonomy": state.source_row.get("original_taxonomy", ""),
                "resolved_taxonomy": resolved_taxonomy,
                "lowest_reliable_rank": state.source_row.get("lowest_reliable_rank", ""),
                "unresolved_ranks": state.source_row.get("unresolved_ranks", ""),
                "genus_placeholder": state.rank_placeholders.get("genus", ""),
                "species_placeholder": state.rank_placeholders.get("species", ""),
                "cluster_key_genus": state.rank_cluster_keys.get("genus", ""),
                "cluster_key_species": state.rank_cluster_keys.get("species", ""),
                "warning": ";".join(sorted(state.warnings)),
            }
        )
    return rows


def _resolve_mapping_rows(states: list[ResolveRecordState]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for state in states:
        rows.append(
            {
                "seq_id": state.record.seq_id,
                "original_silva_taxonomy": state.source_row.get("original_taxonomy", ""),
                "resolved_autotax2_taxonomy": _state_resolved_taxonomy(state),
                "source": SILVA_SOURCE,
                "source_prefix": SILVA_PLACEHOLDER_PREFIX,
                "placeholder_taxa": ",".join(state.placeholder_taxa),
                "warnings": ";".join(sorted(state.warnings)),
            }
        )
    return rows


def _resolve_fasta_records(states: list[ResolveRecordState]) -> list[FastaRecord]:
    return [
        FastaRecord(
            seq_id=state.record.seq_id,
            header=f"{state.record.seq_id} {_state_resolved_taxonomy(state)}",
            sequence=state.record.sequence,
        )
        for state in states
    ]


def _state_resolved_taxonomy(state: ResolveRecordState) -> str:
    for rank in SILVA_RANKS:
        if rank not in state.resolved_labels:
            state.resolved_labels[rank] = _label_for_taxon(rank, _row_taxonomy_value(state.source_row, rank))
    return ";".join(state.resolved_labels[rank] for rank in SILVA_RANKS)


def _mixed_parent_warnings(
    records: list[FastaRecord],
    rows_by_seq_id: dict[str, dict[str, str]],
    genus_assignments: dict[str, ClusterAssignment],
) -> dict[str, str]:
    parents_by_cluster: dict[str, set[str]] = {}
    for record in records:
        row = rows_by_seq_id[record.seq_id]
        parent = _reliable_parent_path(row)
        cluster_key = genus_assignments[record.seq_id].cluster_key
        parents_by_cluster.setdefault(cluster_key, set()).add(parent)
    return {
        cluster_key: MIXED_PARENT_WARNING
        for cluster_key, parents in parents_by_cluster.items()
        if len(parents) > 1
    }


def _reliable_parent_path(row: dict[str, str]) -> str:
    lowest_rank = row.get("lowest_reliable_rank", "")
    if lowest_rank not in SILVA_RANKS:
        return ""
    values = _row_taxonomy_values(row)
    return ";".join(values[: SILVA_RANKS.index(lowest_rank) + 1])


def _resolve_one_record(
    record: FastaRecord,
    source_row: dict[str, str],
    genus_assignment: ClusterAssignment,
    species_assignment: ClusterAssignment,
    parent_warning: str,
    allocator: PlaceholderAllocator,
    active_cluster_to_taxon: dict[str, dict[str, str]],
    taxa_by_cluster_key: dict[str, dict[str, str]],
    genus_counts: dict[str, int],
    species_counts: dict[str, int],
) -> dict[str, Any]:
    unresolved_ranks = _split_unresolved_ranks(source_row.get("unresolved_ranks", ""))
    original_values = _row_taxonomy_values(source_row)
    resolved_labels: list[str] = []
    placeholder_taxa: list[str] = []
    genus_placeholder = ""
    species_placeholder = ""
    previous_taxon_id = ""

    for rank, value in zip(SILVA_RANKS, original_values, strict=True):
        if rank not in unresolved_ranks or rank not in PLACEHOLDER_RANKS:
            label = _prefixed_taxon_label(rank, value)
            resolved_labels.append(label)
            previous_taxon_id = label
            continue

        cluster_key = _placeholder_cluster_key(
            rank=rank,
            source_row=source_row,
            genus_assignment=genus_assignment,
            species_assignment=species_assignment,
            resolved_labels=resolved_labels,
        )
        count_key = species_assignment.cluster_key if rank == "species" else genus_assignment.cluster_key
        representative_seq_id = (
            species_assignment.representative_seq_id
            if rank == "species"
            else genus_assignment.representative_seq_id
        )
        member_count = (
            species_counts.get(species_assignment.cluster_key, 1)
            if rank == "species"
            else genus_counts.get(genus_assignment.cluster_key, 1)
        )
        taxon_row = _get_or_create_placeholder_taxon(
            rank=rank,
            cluster_key=cluster_key,
            parent_taxon_id=previous_taxon_id,
            representative_seq_id=representative_seq_id,
            member_count=member_count,
            warning=parent_warning,
            allocator=allocator,
            active_cluster_to_taxon=active_cluster_to_taxon,
            taxa_by_cluster_key=taxa_by_cluster_key,
        )
        placeholder_name = taxon_row["name"]
        resolved_labels.append(placeholder_name)
        previous_taxon_id = taxon_row["taxon_id"]
        placeholder_taxa.append(placeholder_name)
        if rank == "genus":
            genus_placeholder = placeholder_name
        elif rank == "species":
            species_placeholder = placeholder_name

    resolved_taxonomy = ";".join(resolved_labels)
    warning = parent_warning
    return {
        "member_row": {
            "seq_id": record.seq_id,
            "original_taxonomy": source_row.get("original_taxonomy", ""),
            "resolved_taxonomy": resolved_taxonomy,
            "lowest_reliable_rank": source_row.get("lowest_reliable_rank", ""),
            "unresolved_ranks": source_row.get("unresolved_ranks", ""),
            "genus_placeholder": genus_placeholder,
            "species_placeholder": species_placeholder,
            "cluster_key_genus": genus_assignment.cluster_key,
            "cluster_key_species": _species_cluster_key(source_row, species_assignment, resolved_labels),
            "warning": warning,
        },
        "mapping_row": {
            "seq_id": record.seq_id,
            "original_silva_taxonomy": source_row.get("original_taxonomy", ""),
            "resolved_autotax2_taxonomy": resolved_taxonomy,
            "source": SILVA_SOURCE,
            "source_prefix": SILVA_PLACEHOLDER_PREFIX,
            "placeholder_taxa": ",".join(placeholder_taxa),
            "warnings": warning,
        },
        "fasta_record": FastaRecord(
            seq_id=record.seq_id,
            header=f"{record.seq_id} {resolved_taxonomy}",
            sequence=record.sequence,
        ),
    }


def _placeholder_cluster_key(
    rank: str,
    source_row: dict[str, str],
    genus_assignment: ClusterAssignment,
    species_assignment: ClusterAssignment,
    resolved_labels: list[str],
) -> str:
    if rank == "species":
        return _species_cluster_key(source_row, species_assignment, resolved_labels)
    return f"{rank}:{genus_assignment.cluster_key}"


def _species_cluster_key(
    source_row: dict[str, str],
    species_assignment: ClusterAssignment,
    resolved_labels: list[str],
) -> str:
    genus_context = next(
        (label for label in reversed(resolved_labels) if label.startswith("g__")),
        _prefixed_taxon_label("genus", source_row.get("genus", "")),
    )
    return f"species:{genus_context}:{species_assignment.cluster_key}"


def _get_or_create_placeholder_taxon(
    rank: str,
    cluster_key: str,
    parent_taxon_id: str,
    representative_seq_id: str,
    member_count: int,
    warning: str,
    allocator: PlaceholderAllocator,
    active_cluster_to_taxon: dict[str, dict[str, str]],
    taxa_by_cluster_key: dict[str, dict[str, str]],
) -> dict[str, str]:
    if cluster_key in taxa_by_cluster_key:
        return taxa_by_cluster_key[cluster_key]

    existing = active_cluster_to_taxon.get(cluster_key)
    if existing is not None:
        name = existing["name"]
        taxon_id = existing.get("taxon_id") or name
        allocator.reserve(name)
    else:
        name = allocator.allocate(PLACEHOLDER_RANK_BY_NAME[rank], SILVA_PLACEHOLDER_PREFIX)
        taxon_id = name

    row = {
        "taxon_id": taxon_id,
        "rank": rank,
        "name": name,
        "parent_taxon_id": parent_taxon_id,
        "source": SILVA_SOURCE,
        "source_prefix": SILVA_PLACEHOLDER_PREFIX,
        "cluster_key": cluster_key,
        "status": "active",
        "representative_seq_id": representative_seq_id,
        "member_count": str(member_count),
        "warning": warning,
    }
    taxa_by_cluster_key[cluster_key] = row
    return row


def _resolve_output_paths(silva_dir: Path, dry_run: bool) -> dict[str, Path]:
    suffix = ".dry_run.tsv" if dry_run else ".tsv"
    fasta_name = "silva_unresolved.resolved.dry_run.fa" if dry_run else "silva_unresolved.resolved.fa"
    return {
        "taxa": silva_dir / f"silva_unresolved_taxa{suffix}",
        "members": silva_dir / f"silva_unresolved_members{suffix}",
        "mapping": silva_dir / f"silva_unresolved_mapping{suffix}",
        "evidence": silva_dir / f"silva_unresolved_evidence{suffix}",
        "fasta": silva_dir / fasta_name,
    }


def _write_resolve_outputs(
    taxa_rows: list[dict[str, str]],
    member_rows: list[dict[str, str]],
    mapping_rows: list[dict[str, str]],
    evidence_rows: list[dict[str, str]],
    resolved_records: list[FastaRecord],
    output_paths: dict[str, Path],
) -> None:
    _write_tsv(
        taxa_rows,
        output_paths["taxa"],
        [
            "taxon_id",
            "rank",
            "name",
            "parent_taxon_id",
            "source",
            "source_prefix",
            "cluster_key",
            "status",
            "representative_seq_id",
            "member_count",
            "warning",
        ],
    )
    _write_tsv(
        member_rows,
        output_paths["members"],
        [
            "seq_id",
            "original_taxonomy",
            "resolved_taxonomy",
            "lowest_reliable_rank",
            "unresolved_ranks",
            "genus_placeholder",
            "species_placeholder",
            "cluster_key_genus",
            "cluster_key_species",
            "warning",
        ],
    )
    _write_tsv(
        mapping_rows,
        output_paths["mapping"],
        [
            "seq_id",
            "original_silva_taxonomy",
            "resolved_autotax2_taxonomy",
            "source",
            "source_prefix",
            "placeholder_taxa",
            "warnings",
        ],
    )
    _write_tsv(
        evidence_rows,
        output_paths["evidence"],
        RESOLVE_EVIDENCE_FIELDS,
    )
    write_fasta(resolved_records, output_paths["fasta"])


def _update_resolve_registry(
    registry_dir: Path,
    taxa_rows: list[dict[str, str]],
    allocator: PlaceholderAllocator,
) -> None:
    taxon_nodes_path = registry_dir / "taxon_nodes.tsv"
    name_index_path = registry_dir / "name_index.tsv"
    cluster_to_taxon_path = registry_dir / "cluster_to_taxon.tsv"

    existing_taxon_rows = _read_tsv(taxon_nodes_path)
    parent_lookup = _parent_taxon_lookup(existing_taxon_rows, taxa_rows)
    taxon_fieldnames = _merge_fieldnames(
        [
            "taxon_id",
            "rank",
            "name",
            "parent_taxon_id",
            "taxonomy_path",
            "protected",
            "is_silva_named",
            "source",
            "source_prefix",
            "is_placeholder",
            "is_silva_unresolved",
            "status",
            "cluster_key",
            "representative_seq_id",
            "member_count",
            "warning",
        ],
        existing_taxon_rows,
    )
    rows_by_taxon_id = {row.get("taxon_id", ""): dict(row) for row in existing_taxon_rows}
    for taxon in taxa_rows:
        row = dict(taxon)
        row["parent_taxon_id"] = parent_lookup.get(row.get("parent_taxon_id", ""), row.get("parent_taxon_id", ""))
        row.update(
            {
                "taxonomy_path": taxon["name"],
                "protected": "true",
                "is_silva_named": "false",
                "is_placeholder": "true",
                "is_silva_unresolved": "true",
            }
        )
        rows_by_taxon_id[row["taxon_id"]] = row
    _write_tsv(list(rows_by_taxon_id.values()), taxon_nodes_path, taxon_fieldnames)

    name_rows = _read_tsv(name_index_path)
    name_fieldnames = _merge_fieldnames(["name", "rank", "taxon_id", "protected", "source"], name_rows)
    name_keys = {(row.get("name", ""), row.get("rank", "")) for row in name_rows}
    for taxon in taxa_rows:
        key = (taxon["name"], taxon["rank"])
        if key in name_keys:
            continue
        name_rows.append(
            {
                "name": taxon["name"],
                "rank": taxon["rank"],
                "taxon_id": taxon["taxon_id"],
                "protected": "true",
                "source": SILVA_SOURCE,
            }
        )
        name_keys.add(key)
    _write_tsv(name_rows, name_index_path, name_fieldnames)

    cluster_rows = _read_tsv(cluster_to_taxon_path)
    cluster_fieldnames = _merge_fieldnames(
        ["cluster_key", "taxon_id", "rank", "name", "status", "source_prefix"],
        cluster_rows,
    )
    cluster_rows_by_key = {row.get("cluster_key", ""): dict(row) for row in cluster_rows}
    for taxon in taxa_rows:
        cluster_rows_by_key[taxon["cluster_key"]] = {
            "cluster_key": taxon["cluster_key"],
            "taxon_id": taxon["taxon_id"],
            "rank": taxon["rank"],
            "name": taxon["name"],
            "status": "active",
            "source_prefix": SILVA_PLACEHOLDER_PREFIX,
        }
    _write_tsv(list(cluster_rows_by_key.values()), cluster_to_taxon_path, cluster_fieldnames)

    existing_counters = _read_placeholder_counters(registry_dir)
    silva_counters = existing_counters.setdefault(SILVA_PLACEHOLDER_PREFIX, {})
    for rank in PLACEHOLDER_RANKS:
        silva_counters[rank] = allocator.next_ordinal(
            PLACEHOLDER_RANK_BY_NAME[rank],
            SILVA_PLACEHOLDER_PREFIX,
        )
    _write_placeholder_counters(registry_dir, existing_counters)


def _validate_resolve_threshold_scope(
    family_id: float,
    order_id: float,
    class_id: float,
    phylum_id: float,
) -> None:
    """Fail loudly for resolver threshold controls that are not implemented yet."""
    unsupported = {
        "family_id": (family_id, DEFAULT_FAMILY_ID),
        "order_id": (order_id, DEFAULT_ORDER_ID),
        "class_id": (class_id, DEFAULT_CLASS_ID),
        "phylum_id": (phylum_id, DEFAULT_PHYLUM_ID),
    }
    changed = [
        f"{name}={value}"
        for name, (value, default) in unsupported.items()
        if value != default
    ]
    if changed:
        raise NotImplementedError(
            "The SILVA resolver currently clusters unresolved SILVA records at genus "
            f"and species thresholds only; unsupported threshold overrides: {', '.join(changed)}."
        )


def _read_placeholder_counters(registry_dir: Path) -> dict[str, dict[str, int]]:
    counters: dict[str, dict[str, int]] = {}
    yaml_path = registry_dir / "placeholder_counters.yaml"
    if yaml_path.exists():
        loaded = yaml.safe_load(yaml_path.read_text(encoding="utf-8")) or {}
        if isinstance(loaded, dict):
            for prefix, rank_values in loaded.items():
                if not isinstance(rank_values, dict):
                    continue
                counters[str(prefix)] = {
                    str(rank): int(value)
                    for rank, value in rank_values.items()
                    if str(value).isdigit() or isinstance(value, int)
                }

    for row in _read_tsv(registry_dir / "placeholder_counters.tsv"):
        prefix = row.get("source_prefix", "")
        rank = row.get("rank", "")
        next_ordinal = row.get("next_ordinal", "")
        if not prefix or rank not in PLACEHOLDER_RANK_BY_NAME or not next_ordinal.isdigit():
            continue
        counters.setdefault(prefix, {})[rank] = int(next_ordinal)
    return counters


def _write_placeholder_counters(registry_dir: Path, counters: dict[str, dict[str, int]]) -> None:
    yaml_path = registry_dir / "placeholder_counters.yaml"
    yaml_path.write_text(yaml.safe_dump(counters, sort_keys=True), encoding="utf-8")
    rows = [
        {
            "source_prefix": source_prefix,
            "rank": rank,
            "next_ordinal": str(next_ordinal),
        }
        for source_prefix in sorted(counters)
        for rank, next_ordinal in sorted(counters[source_prefix].items())
    ]
    _write_tsv(
        rows,
        registry_dir / "placeholder_counters.tsv",
        ["source_prefix", "rank", "next_ordinal"],
    )


def _parent_taxon_lookup(
    existing_taxon_rows: list[dict[str, str]],
    taxa_rows: list[dict[str, str]],
) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for row in existing_taxon_rows:
        taxon_id = row.get("taxon_id", "")
        rank = row.get("rank", "")
        name = row.get("name", "")
        if not taxon_id or not rank or not name:
            continue
        lookup[name] = taxon_id
        if rank in RANK_PREFIXES:
            lookup[f"{RANK_PREFIXES[rank]}__{name}"] = taxon_id
    for row in taxa_rows:
        taxon_id = row.get("taxon_id", "")
        name = row.get("name", "")
        if taxon_id and name:
            lookup[name] = taxon_id
    return lookup


def _load_placeholder_state(registry_dir: Path) -> dict[str, Any]:
    issued: set[str] = set()
    deprecated: set[str] = set()
    active_cluster_to_taxon: dict[str, dict[str, str]] = {}

    for row in _read_tsv(registry_dir / "taxon_nodes.tsv"):
        name = row.get("name", "")
        try:
            parse_placeholder_id(name)
        except ValueError:
            continue
        status = row.get("status", "active")
        if status in {"deprecated", "superseded"}:
            deprecated.add(name)
        else:
            issued.add(name)

    for row in _read_tsv(registry_dir / "cluster_to_taxon.tsv"):
        name = row.get("name", "")
        cluster_key = row.get("cluster_key", "")
        if not name or not cluster_key:
            continue
        try:
            parse_placeholder_id(name)
        except ValueError:
            continue
        if row.get("status", "active") == "active":
            issued.add(name)
            active_cluster_to_taxon[cluster_key] = dict(row)
        else:
            deprecated.add(name)

    for row in _read_tsv(registry_dir / "placeholder_counters.tsv"):
        rank = row.get("rank", "")
        next_ordinal = row.get("next_ordinal", "")
        if rank not in PLACEHOLDER_RANK_BY_NAME or not next_ordinal.isdigit():
            continue
        ordinal = int(next_ordinal)
        if ordinal > 1:
            issued.add(
                make_placeholder_id(
                    PLACEHOLDER_RANK_BY_NAME[rank],
                    SILVA_PLACEHOLDER_PREFIX,
                    ordinal - 1,
                )
            )

    return {
        "issued_placeholders": issued,
        "deprecated_placeholders": deprecated,
        "active_cluster_to_taxon": active_cluster_to_taxon,
    }


def _row_taxonomy_values(row: dict[str, str]) -> tuple[str, str, str, str, str, str, str]:
    return tuple(row.get(rank, "") for rank in SILVA_RANKS)  # type: ignore[return-value]


def _row_taxonomy_value(row: dict[str, str], rank: str) -> str:
    return row.get(rank, "")


def _domain_reject_reason(domain: str) -> str:
    if not domain.strip() or domain.startswith(MISSING_VALUE_PREFIX):
        return "empty_domain"
    if _unresolved_reason("domain", domain):
        return "unresolved_domain"
    return ""


def _silva_record_reject_reason(
    fasta_record: FastaRecord,
    taxonomy: SilvaTaxonomy,
    metadata: dict[str, TypeStrainMetadata],
) -> str:
    domain_reject_reason = _domain_reject_reason(taxonomy.domain)
    if domain_reject_reason:
        return domain_reject_reason
    if _invalid_silva_sequence(fasta_record.sequence) and not _is_type_strain_metadata(
        _metadata_for_seq_id(fasta_record.seq_id, metadata)
    ):
        return SILVA_INVALID_SEQUENCE_REASON
    return ""


def _invalid_silva_sequence(sequence: str) -> bool:
    return SILVA_ATGC_RE.fullmatch(sequence) is None


def _split_unresolved_ranks(value: str) -> tuple[str, ...]:
    return tuple(rank.strip() for rank in value.split(",") if rank.strip())


def _prefixed_taxon_label(rank: str, value: str) -> str:
    prefix = RANK_PREFIXES[rank]
    clean_value = value.strip() or f"{MISSING_VALUE_PREFIX}{rank}"
    return f"{prefix}__{clean_value}"


def _label_for_taxon(rank: str, value: str) -> str:
    clean_value = value.strip() or f"{MISSING_VALUE_PREFIX}{rank}"
    prefix = RANK_PREFIXES[rank]
    if clean_value.startswith(f"{prefix}__"):
        return clean_value
    return f"{prefix}__{clean_value}"


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=chr(9)))


def _merge_fieldnames(preferred: list[str], rows: list[dict[str, str]]) -> list[str]:
    fieldnames = list(preferred)
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    return fieldnames


def _unresolved_reason(rank: str, value: str) -> str:
    normalized = value.strip()
    if not normalized:
        return "missing"
    if normalized.lower().startswith(MISSING_VALUE_PREFIX):
        return "missing"
    normalized_value = _normalize_taxon_value(normalized)
    if normalized_value in STANDALONE_UNRESOLVED_VALUES:
        return f"ambiguous value {normalized_value}"
    if rank == "species" and normalized_value in ENVIRONMENTAL_DESCRIPTOR_VALUES:
        return f"environmental descriptor {normalized_value}"
    if rank == "species" and _species_is_sp_placeholder(normalized):
        return "species sp."
    token = _unresolved_token(normalized)
    if token:
        return f"contains {token}"
    return ""


def _species_is_sp_placeholder(value: str) -> bool:
    return re.search(r"(^|\s)sp\.(\s|$)", value, flags=re.IGNORECASE) is not None


def _canonical_species_name(value: str) -> str:
    clean = value.strip().strip(";")
    if not clean:
        return ""
    clean = re.sub(r"[_\s]+", " ", clean)
    clean = clean.strip("'\" ")
    clean = re.sub(
        r"\s+(str\.?|strain|isolate|clone|culture|subsp\.?|serovar)\b.*$",
        "",
        clean,
        flags=re.IGNORECASE,
    ).strip()
    if _species_is_sp_placeholder(clean):
        return ""

    tokens = re.findall(r"[A-Za-z][A-Za-z0-9.-]*", clean)
    if len(tokens) < 2:
        return ""
    if tokens[0].lower() == "candidatus":
        if len(tokens) < 3:
            return ""
        species_tokens = tokens[:3]
    else:
        species_tokens = tokens[:2]

    candidate = " ".join(species_tokens)
    if _unresolved_token(candidate):
        return ""
    normalized = _normalize_taxon_value(candidate)
    bad_words = {
        "bacterium",
        "archaeon",
        "crenarchaeote",
        "euryarchaeote",
        "group",
        "metagenome",
        "uncultured",
        "unidentified",
        "unknown",
    }
    if any(part in bad_words for part in normalized.split()):
        return ""
    return candidate


def _unresolved_token(value: str) -> str:
    normalized = _normalize_taxon_value(value)
    for token in UNRESOLVED_TOKENS:
        if re.search(rf"(?<![a-z]){re.escape(token)}(?![a-z])", normalized):
            return token
    return ""


def _normalize_taxon_value(value: str) -> str:
    normalized = re.sub(r"[_\-]+", " ", value.strip().lower())
    normalized = re.sub(r"[^a-z0-9\s]+", " ", normalized)
    return re.sub(r"\s+", " ", normalized).strip()


def _legal_name_catalog_rows(catalog: LegalNameCatalog) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for rank in SILVA_RANKS:
        for normalized_name in sorted(catalog.names_by_rank.get(rank, frozenset())):
            key = (rank, normalized_name)
            display_name = catalog.display_names_by_rank_name.get(key, normalized_name)
            canonical_species = ""
            if rank == "species":
                canonical = _canonical_species_name(display_name)
                if canonical:
                    canonical_species = catalog.canonical_species_by_name.get(
                        _normalize_taxon_value(canonical),
                        canonical,
                    )
            rows.append(
                {
                    "rank": rank,
                    "name": display_name,
                    "source": ";".join(sorted(catalog.sources_by_rank_name.get(key, frozenset()))),
                    "canonical_species_name": canonical_species,
                }
            )
    return rows


def _named_taxonomy_rows(records: Iterable[SilvaRecord]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for record in records:
        row = {
            "seq_id": record.seq_id,
            "taxonomy_7rank": record.taxonomy.taxonomy_7rank,
        }
        row.update(record.taxonomy.as_rank_dict())
        rows.append(row)
    return rows


def _unresolved_rows(records: Iterable[SilvaRecord]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for record in records:
        row = {
            "seq_id": record.seq_id,
            "original_taxonomy": record.taxonomy.original_taxonomy,
            "lowest_reliable_rank": record.classification.lowest_reliable_rank,
            "unresolved_ranks": ",".join(record.classification.unresolved_ranks),
            "unresolved_reason": record.classification.unresolved_reason,
        }
        row.update(record.taxonomy.as_rank_dict())
        rows.append(row)
    return rows


def _taxon_node_rows(records: Iterable[SilvaRecord]) -> list[dict[str, str]]:
    rows_by_key: dict[tuple[str, str], dict[str, str]] = {}
    next_id = 1

    for record in records:
        parent_taxon_id = ""
        path_parts: list[str] = []
        for rank, name in zip(SILVA_RANKS, record.taxonomy.values, strict=True):
            path_parts.append(name)
            key = (rank, ";".join(path_parts))
            if key not in rows_by_key:
                taxon_id = f"T{next_id:06d}"
                next_id += 1
                rows_by_key[key] = {
                    "taxon_id": taxon_id,
                    "rank": rank,
                    "name": name,
                    "parent_taxon_id": parent_taxon_id,
                    "taxonomy_path": ";".join(path_parts),
                    "protected": "true",
                    "is_silva_named": "true",
                    "source": SILVA_SOURCE,
                }
            parent_taxon_id = rows_by_key[key]["taxon_id"]

    return list(rows_by_key.values())


def _species_taxon_id_by_taxonomy(taxon_rows: Iterable[dict[str, str]]) -> dict[str, str]:
    return {
        row.get("taxonomy_path", ""): row.get("taxon_id", "")
        for row in taxon_rows
        if row.get("rank") == "species" and row.get("taxonomy_path") and row.get("taxon_id")
    }


def _representative_registry_rows(
    named_records: Iterable[SilvaRecord],
    metadata: dict[str, TypeStrainMetadata],
    species_taxon_by_taxonomy: dict[str, str],
) -> list[dict[str, str]]:
    chosen_by_taxon: dict[str, tuple[int, SilvaRecord]] = {}
    for record in named_records:
        taxon_id = species_taxon_by_taxonomy.get(record.taxonomy.taxonomy_7rank, "")
        if not taxon_id:
            continue
        current = chosen_by_taxon.get(taxon_id)
        candidate_score = 1 if _is_type_strain_metadata(_metadata_for_seq_id(record.seq_id, metadata)) else 0
        if current is None or candidate_score > current[0]:
            chosen_by_taxon[taxon_id] = (candidate_score, record)

    rows: list[dict[str, str]] = []
    for taxon_id in sorted(chosen_by_taxon):
        is_type, record = chosen_by_taxon[taxon_id]
        rows.append(
            {
                "representative_seq_id": record.seq_id,
                "taxon_id": taxon_id,
                "dataset": SILVA_SOURCE,
                "source_category": "named_silva",
                "status": "active",
                "protected": "true",
                "source": SILVA_SOURCE,
                "is_type_strain": _bool_text(bool(is_type)),
                "representative_reason": "type_strain" if is_type else "first_silva_named_for_species",
            }
        )
    return rows


def _is_type_strain_metadata(metadata: TypeStrainMetadata | None) -> bool:
    return metadata is not None and metadata.is_type_strain.strip().lower() in {"true", "1", "yes", "y"}


def _metadata_for_seq_id(
    seq_id: str,
    metadata: dict[str, TypeStrainMetadata],
) -> TypeStrainMetadata | None:
    exact = metadata.get(seq_id)
    if exact is not None:
        return exact
    accession = seq_id.split(".", maxsplit=1)[0]
    if accession != seq_id:
        return metadata.get(accession)
    return None


def _sequence_registry_rows(
    records: Iterable[SilvaRecord],
    metadata: dict[str, TypeStrainMetadata],
    species_taxon_by_taxonomy: dict[str, str],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for record in records:
        type_metadata = _metadata_for_seq_id(record.seq_id, metadata)
        is_named = not record.is_unresolved
        taxon_id = species_taxon_by_taxonomy.get(record.taxonomy.taxonomy_7rank, "") if is_named else ""
        rows.append(
            {
                "seq_id": record.seq_id,
                "taxon_id": taxon_id,
                "source": SILVA_SOURCE,
                "sequence_md5": sequence_md5(record.fasta_record.sequence),
                "sequence_length": str(len(record.fasta_record.sequence)),
                "protected": _bool_text(is_named),
                "is_silva_named": _bool_text(is_named),
                "is_silva_unresolved": _bool_text(record.is_unresolved),
                "original_taxonomy": record.taxonomy.original_taxonomy,
                "taxonomy_7rank": record.taxonomy.taxonomy_7rank,
                "is_type_strain": type_metadata.is_type_strain if type_metadata else "",
                "type_species_name": type_metadata.species_name if type_metadata else "",
                "type_strain_id": type_metadata.strain_id if type_metadata else "",
                "type_source": type_metadata.source if type_metadata else "",
                "type_evidence": type_metadata.evidence if type_metadata else "",
            }
        )
    return rows


def _name_index_rows(taxon_rows: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    return [
        {
            "name": row["name"],
            "rank": row["rank"],
            "taxon_id": row["taxon_id"],
            "protected": row["protected"],
            "source": row["source"],
        }
        for row in taxon_rows
    ]


def _write_tsv(rows: Iterable[dict[str, Any]], path: Path, fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter=chr(9),
            lineterminator=chr(10),
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def _bool_text(value: bool) -> str:
    return "true" if value else "false"
