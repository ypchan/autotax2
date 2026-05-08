"""SILVA backbone initialization helpers."""

from __future__ import annotations

import csv
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from autotax2.placeholders import (
    PlaceholderAllocator,
    PlaceholderRank,
    make_placeholder_id,
    parse_placeholder_id,
)
from autotax2.io import FastaRecord, read_fasta, write_fasta
from autotax2.registry import sequence_md5
from autotax2.vsearch import read_uc_assignments, run_vsearch_cluster, write_singleton_uc


SILVA_PLACEHOLDER_PREFIX = "SILVA"
SILVA_SOURCE = "SILVA138.2_NR99"
SILVA_RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
RANK_PREFIXES = {
    "domain": "d",
    "phylum": "p",
    "class": "c",
    "order": "o",
    "family": "f",
    "genus": "g",
    "species": "s",
}
PLACEHOLDER_RANKS = ("class", "order", "family", "genus", "species")
PLACEHOLDER_RANK_BY_NAME = {
    "class": PlaceholderRank.CLASS,
    "order": PlaceholderRank.ORDER,
    "family": PlaceholderRank.FAMILY,
    "genus": PlaceholderRank.GENUS,
    "species": PlaceholderRank.SPECIES,
}
MIXED_PARENT_WARNING = "mixed_silva_parent_unresolved_cluster"
UNRESOLVED_TOKENS = (
    "unidentified",
    "unclassified",
    "uncultured",
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
class ClusterAssignment:
    """Resolved cluster assignment for one sequence."""

    cluster_key: str
    cluster_id: str
    representative_seq_id: str


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
    taxonomy_parts = [part.strip() for part in original_taxonomy.split(";") if part.strip()]
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


def classify_unresolved(taxonomy: SilvaTaxonomy) -> UnresolvedClassification:
    """Detect unresolved ranks in a seven-rank SILVA taxonomy."""
    direct_unresolved: list[str] = []
    reasons: list[str] = []

    for rank, value in zip(SILVA_RANKS, taxonomy.values, strict=True):
        reason = _unresolved_reason(rank, value)
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


def read_silva_records(path: str | Path, domain: str | None = None) -> list[SilvaRecord]:
    """Read SILVA FASTA records, optionally filtering by taxonomy domain."""
    requested_domain = domain.strip() if domain is not None else None
    records: list[SilvaRecord] = []

    for fasta_record in read_fasta(path):
        seq_id, taxonomy = parse_silva_header(fasta_record.header)
        if seq_id != fasta_record.seq_id:
            raise ValueError(
                f"Header parser ID mismatch for {fasta_record.header!r}: "
                f"{seq_id!r} != {fasta_record.seq_id!r}"
            )
        if requested_domain is not None and taxonomy.domain != requested_domain:
            continue
        records.append(
            SilvaRecord(
                fasta_record=fasta_record,
                taxonomy=taxonomy,
                classification=classify_unresolved(taxonomy),
            )
        )

    return records


def initialize_silva_build(
    silva_fasta: str | Path,
    outdir: str | Path,
    domain: str | None = None,
    type_strain_metadata: str | Path | None = None,
) -> dict[str, int]:
    """Initialize an autotax2 build from SILVA FASTA records."""
    output_dir = Path(outdir)
    registry_dir = output_dir / "registry"
    silva_dir = output_dir / "silva"
    logs_dir = output_dir / "logs"
    for directory in (registry_dir, silva_dir, logs_dir):
        directory.mkdir(parents=True, exist_ok=True)

    records = read_silva_records(silva_fasta, domain=domain)
    named_records = [record for record in records if not record.is_unresolved]
    unresolved_records = [record for record in records if record.is_unresolved]
    metadata = read_type_strain_metadata(type_strain_metadata) if type_strain_metadata else {}

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
    taxon_rows = _taxon_node_rows(named_records)
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
        _sequence_registry_rows(records, metadata),
        registry_dir / "sequence_registry.tsv",
        [
            "seq_id",
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
        _name_index_rows(taxon_rows),
        registry_dir / "name_index.tsv",
        ["name", "rank", "taxon_id", "protected", "source"],
    )
    (logs_dir / "init.log").write_text(
        (
            f"silva_fasta={Path(silva_fasta)}\n"
            f"domain={domain or ''}\n"
            f"records={len(records)}\n"
            f"named={len(named_records)}\n"
            f"unresolved={len(unresolved_records)}\n"
        ),
        encoding="utf-8",
    )

    return {
        "records": len(records),
        "named": len(named_records),
        "unresolved": len(unresolved_records),
    }


def read_type_strain_metadata(path: str | Path) -> dict[str, TypeStrainMetadata]:
    """Read optional type-strain metadata keyed by SILVA sequence ID."""
    metadata: dict[str, TypeStrainMetadata] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=chr(9))
        required = {"seq_id", "is_type_strain", "species_name", "strain_id", "source", "evidence"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                "Type strain metadata must contain fields: "
                "seq_id, is_type_strain, species_name, strain_id, source, evidence"
            )
        for row in reader:
            seq_id = (row.get("seq_id") or "").strip()
            if not seq_id:
                continue
            metadata[seq_id] = TypeStrainMetadata(
                seq_id=seq_id,
                is_type_strain=(row.get("is_type_strain") or "").strip(),
                species_name=(row.get("species_name") or "").strip(),
                strain_id=(row.get("strain_id") or "").strip(),
                source=(row.get("source") or "").strip(),
                evidence=(row.get("evidence") or "").strip(),
            )
    return metadata


def resolve_silva_unresolved_build(
    build: str | Path,
    threads: int = 4,
    species_id: float = 0.987,
    genus_id: float = 0.945,
    family_id: float = 0.865,
    order_id: float = 0.820,
    class_id: float = 0.785,
    floor_id: float = 0.750,
    vsearch_bin: str = "vsearch",
    iddef: int = 2,
    dry_run: bool = False,
) -> ResolveSilvaSummary:
    """Resolve SILVA unresolved records into mutable placeholder taxa."""
    del family_id, order_id, class_id, floor_id

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
            resolved_records=[],
            output_paths=output_paths,
        )
        return ResolveSilvaSummary(unresolved_records=0, placeholder_taxa=0, dry_run=dry_run)

    genus_uc = cluster_dir / f"genus_{genus_id:.3f}.uc"
    species_uc = cluster_dir / f"species_{species_id:.3f}.uc"
    _ensure_uc(
        fasta_path=unresolved_fasta_path,
        uc_path=genus_uc,
        identity=genus_id,
        records=records,
        threads=threads,
        vsearch_bin=vsearch_bin,
        iddef=iddef,
    )
    _ensure_uc(
        fasta_path=unresolved_fasta_path,
        uc_path=species_uc,
        identity=species_id,
        records=records,
        threads=threads,
        vsearch_bin=vsearch_bin,
        iddef=iddef,
    )
    genus_assignments = _assignments_by_seq_id(genus_uc, records, "genus")
    species_assignments = _assignments_by_seq_id(species_uc, records, "species")
    genus_counts = _cluster_counts(genus_assignments)
    species_counts = _cluster_counts(species_assignments)
    state = _load_placeholder_state(registry_dir)
    allocator = PlaceholderAllocator(
        issued_ids=state["issued_placeholders"],
        deprecated_ids=state["deprecated_placeholders"],
    )
    active_cluster_to_taxon = state["active_cluster_to_taxon"]
    parent_warnings = _mixed_parent_warnings(records, rows_by_seq_id, genus_assignments)

    taxa_by_cluster_key: dict[str, dict[str, str]] = {}
    member_rows: list[dict[str, str]] = []
    mapping_rows: list[dict[str, str]] = []
    resolved_records: list[FastaRecord] = []

    for record in records:
        source_row = rows_by_seq_id[record.seq_id]
        resolved = _resolve_one_record(
            record=record,
            source_row=source_row,
            genus_assignment=genus_assignments[record.seq_id],
            species_assignment=species_assignments[record.seq_id],
            parent_warning=parent_warnings.get(genus_assignments[record.seq_id].cluster_key, ""),
            allocator=allocator,
            active_cluster_to_taxon=active_cluster_to_taxon,
            taxa_by_cluster_key=taxa_by_cluster_key,
            genus_counts=genus_counts,
            species_counts=species_counts,
        )
        member_rows.append(resolved["member_row"])
        mapping_rows.append(resolved["mapping_row"])
        resolved_records.append(resolved["fasta_record"])

    taxa_rows = list(taxa_by_cluster_key.values())
    _write_resolve_outputs(
        taxa_rows=taxa_rows,
        member_rows=member_rows,
        mapping_rows=mapping_rows,
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
        "fasta": silva_dir / fasta_name,
    }


def _write_resolve_outputs(
    taxa_rows: list[dict[str, str]],
    member_rows: list[dict[str, str]],
    mapping_rows: list[dict[str, str]],
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
    write_fasta(resolved_records, output_paths["fasta"])


def _update_resolve_registry(
    registry_dir: Path,
    taxa_rows: list[dict[str, str]],
    allocator: PlaceholderAllocator,
) -> None:
    taxon_nodes_path = registry_dir / "taxon_nodes.tsv"
    name_index_path = registry_dir / "name_index.tsv"
    cluster_to_taxon_path = registry_dir / "cluster_to_taxon.tsv"
    counters_path = registry_dir / "placeholder_counters.tsv"

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

    _write_tsv(
        [
            {
                "source_prefix": SILVA_PLACEHOLDER_PREFIX,
                "rank": rank,
                "next_ordinal": str(
                    allocator.next_ordinal(PLACEHOLDER_RANK_BY_NAME[rank], SILVA_PLACEHOLDER_PREFIX)
                ),
            }
            for rank in PLACEHOLDER_RANKS
        ],
        counters_path,
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


def _split_unresolved_ranks(value: str) -> tuple[str, ...]:
    return tuple(rank.strip() for rank in value.split(",") if rank.strip())


def _prefixed_taxon_label(rank: str, value: str) -> str:
    prefix = RANK_PREFIXES[rank]
    clean_value = value.strip() or f"{MISSING_VALUE_PREFIX}{rank}"
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
    if rank == "species" and _species_is_sp_placeholder(normalized):
        return "species sp."
    token = _unresolved_token(normalized)
    if token:
        return f"contains {token}"
    return ""


def _species_is_sp_placeholder(value: str) -> bool:
    return re.search(r"(^|\s)sp\.(\s|$)", value, flags=re.IGNORECASE) is not None


def _unresolved_token(value: str) -> str:
    normalized = value.lower()
    for token in UNRESOLVED_TOKENS:
        if re.search(rf"(?<![a-z]){re.escape(token)}(?![a-z])", normalized):
            return token
    return ""


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


def _sequence_registry_rows(
    records: Iterable[SilvaRecord],
    metadata: dict[str, TypeStrainMetadata],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for record in records:
        type_metadata = metadata.get(record.seq_id)
        is_named = not record.is_unresolved
        rows.append(
            {
                "seq_id": record.seq_id,
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
