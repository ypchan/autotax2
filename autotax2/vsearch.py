"""VSEARCH wrappers, clustering, and registry search."""

from __future__ import annotations

import csv
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from autotax2.io import FastaRecord, read_fasta, write_fasta


DEFAULT_IDDEF = 2
EXPECTED_VSEARCH_VERSION = "not_pinned"
CLUSTER_RANKS = ("species", "genus", "family", "order", "class")
CLUSTER_MEMBER_FIELDS = [
    "cluster_id",
    "centroid_id",
    "member_id",
    "identity_to_centroid",
    "record_type",
    "rank_threshold",
]
REGISTRY_HIT_FIELDS = [
    "query",
    "target",
    "identity",
    "alnlen",
    "qlo",
    "qhi",
    "tlo",
    "thi",
    "ql",
    "tl",
    "bits",
    "query_coverage",
    "target_coverage",
]
SUMMARY_FIELDS = [
    "dataset",
    "prefix",
    "input_sequences",
    "species_centroids",
    "genus_centroids",
    "family_centroids",
    "order_centroids",
    "class_centroids",
    "registry_representatives",
    "registry_hits_raw",
    "registry_hits_filtered",
    "iddef",
    "min_query_cov",
    "min_target_cov",
    "near_best_delta",
    "vsearch_version",
]
TOOL_VERSION_FIELDS = [
    "tool",
    "expected_version",
    "detected_version",
    "status",
    "command",
]


@dataclass(frozen=True)
class UcClusterAssignment:
    """One sequence-to-cluster assignment parsed from a VSEARCH .uc file."""

    seq_id: str
    cluster_id: str
    representative_seq_id: str


@dataclass(frozen=True)
class UcRecord:
    """One raw VSEARCH/USEARCH .uc record."""

    record_type: str
    cluster_number: str
    size_or_length: str
    percent_identity: str
    strand: str
    query_start: str
    seed_start: str
    alignment: str
    query_label: str
    target_label: str


@dataclass(frozen=True)
class ClusterMembership:
    """Cluster membership resolved from .uc records."""

    cluster_id: str
    centroid_id: str
    member_id: str
    identity_to_centroid: str
    record_type: str
    rank_threshold: str


@dataclass(frozen=True)
class RegistryHit:
    """One VSEARCH registry search hit with computed coverage."""

    query: str
    target: str
    identity: float
    alnlen: int
    qlo: int
    qhi: int
    tlo: int
    thi: int
    ql: int
    tl: int
    bits: float
    query_coverage: float
    target_coverage: float


@dataclass(frozen=True)
class ClusterSearchSummary:
    """Summary from VSEARCH clustering and registry search."""

    dataset: str
    dataset_dir: Path
    registry_hits_filtered: int


def default_iddef() -> int:
    """Return the default VSEARCH --iddef value."""
    return DEFAULT_IDDEF


def build_cluster_fast_command(
    input_fasta: str | Path,
    uc_path: str | Path,
    centroids_path: str | Path | None,
    identity: float,
    threads: int = 4,
    vsearch_bin: str = "vsearch",
    iddef: int = DEFAULT_IDDEF,
) -> list[str]:
    """Build a VSEARCH cluster_fast command."""
    command = [
        vsearch_bin,
        "--cluster_fast",
        str(input_fasta),
        "--id",
        f"{identity:.3f}",
        "--iddef",
        str(iddef),
    ]
    if centroids_path is not None:
        command.extend(["--centroids", str(centroids_path)])
    command.extend(["--uc", str(uc_path), "--threads", str(threads)])
    return command


def build_usearch_global_command(
    query_fasta: str | Path,
    db_fasta: str | Path,
    userout_path: str | Path,
    identity: float,
    iddef: int = DEFAULT_IDDEF,
    maxaccepts: int = 50,
    maxrejects: int = 256,
    strand: str = "plus",
    vsearch_bin: str = "vsearch",
) -> list[str]:
    """Build a VSEARCH usearch_global command for registry search."""
    return [
        vsearch_bin,
        "--usearch_global",
        str(query_fasta),
        "--db",
        str(db_fasta),
        "--id",
        f"{identity:.3f}",
        "--iddef",
        str(iddef),
        "--maxaccepts",
        str(maxaccepts),
        "--maxrejects",
        str(maxrejects),
        "--strand",
        strand,
        "--userout",
        str(userout_path),
        "--userfields",
        "query+target+id+alnlen+qlo+qhi+tlo+thi+ql+tl+bits",
    ]


def check_vsearch_version(vsearch_bin: str = "vsearch") -> str:
    """Return the detected VSEARCH version string, or an empty string if unavailable."""
    try:
        completed = subprocess.run(
            [vsearch_bin, "--version"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except (FileNotFoundError, OSError, subprocess.CalledProcessError):
        return ""
    output = f"{completed.stdout}\n{completed.stderr}".strip()
    match = re.search(r"vsearch\s+v?(\d+(?:\.\d+)+)", output, re.I)
    if match:
        return match.group(1)
    match = re.search(r"(\d+(?:\.\d+)+)", output)
    return match.group(1) if match else ""


def run_vsearch_cluster(
    fasta_path: str | Path,
    uc_path: str | Path,
    identity: float,
    threads: int = 4,
    vsearch_bin: str = "vsearch",
    iddef: int = DEFAULT_IDDEF,
    centroids_path: str | Path | None = None,
) -> list[str]:
    """Run VSEARCH clustering and write .uc/centroid outputs."""
    command = build_cluster_fast_command(
        input_fasta=fasta_path,
        uc_path=uc_path,
        centroids_path=centroids_path,
        identity=identity,
        threads=threads,
        vsearch_bin=vsearch_bin,
        iddef=iddef,
    )
    subprocess.run(command, check=True)
    return command


def run_vsearch_search(command: list[str]) -> None:
    """Run a prepared VSEARCH registry search command."""
    subprocess.run(command, check=True)


def parse_uc_records(path: str | Path) -> list[UcRecord]:
    """Parse VSEARCH/USEARCH .uc records."""
    records: list[UcRecord] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            records.append(
                UcRecord(
                    record_type=fields[0],
                    cluster_number=fields[1],
                    size_or_length=fields[2],
                    percent_identity=fields[3],
                    strand=fields[4],
                    query_start=fields[5],
                    seed_start=fields[6],
                    alignment=fields[7],
                    query_label=fields[8],
                    target_label=fields[9] if len(fields) > 9 else "",
                )
            )
    return records


def cluster_memberships_from_uc(
    path: str | Path,
    rank_threshold: str,
) -> list[ClusterMembership]:
    """Resolve centroid/member rows from a VSEARCH .uc file."""
    records = parse_uc_records(path)
    centroid_by_cluster: dict[str, str] = {
        record.cluster_number: record.query_label
        for record in records
        if record.record_type == "S"
    }
    memberships: list[ClusterMembership] = []
    for record in records:
        if record.record_type == "S":
            centroid = record.query_label
            identity = "100.0"
        elif record.record_type == "H":
            centroid = record.target_label if record.target_label not in {"", "*"} else ""
            centroid = centroid or centroid_by_cluster.get(record.cluster_number, record.query_label)
            identity = record.percent_identity
        else:
            continue
        memberships.append(
            ClusterMembership(
                cluster_id=record.cluster_number,
                centroid_id=centroid,
                member_id=record.query_label,
                identity_to_centroid=identity,
                record_type=record.record_type,
                rank_threshold=rank_threshold,
            )
        )
    return memberships


def read_uc_assignments(path: str | Path) -> list[UcClusterAssignment]:
    """Read sequence cluster assignments from a VSEARCH .uc file."""
    return [
        UcClusterAssignment(
            seq_id=membership.member_id,
            cluster_id=membership.cluster_id,
            representative_seq_id=membership.centroid_id,
        )
        for membership in cluster_memberships_from_uc(path, rank_threshold="")
    ]


def parse_registry_hits(path: str | Path) -> list[RegistryHit]:
    """Parse VSEARCH userout hits and compute query/target coverage."""
    hits: list[RegistryHit] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("query\t"):
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                continue
            alnlen = int(float(fields[3]))
            ql = int(float(fields[8]))
            tl = int(float(fields[9]))
            hits.append(
                RegistryHit(
                    query=fields[0],
                    target=fields[1],
                    identity=float(fields[2]),
                    alnlen=alnlen,
                    qlo=int(float(fields[4])),
                    qhi=int(float(fields[5])),
                    tlo=int(float(fields[6])),
                    thi=int(float(fields[7])),
                    ql=ql,
                    tl=tl,
                    bits=float(fields[10]),
                    query_coverage=alnlen / ql if ql else 0.0,
                    target_coverage=alnlen / tl if tl else 0.0,
                )
            )
    return hits


def filter_registry_hits(
    hits: Iterable[RegistryHit],
    floor_id: float,
    min_query_cov: float,
    min_target_cov: float,
) -> list[RegistryHit]:
    """Filter registry hits while preserving all passing near-best candidates."""
    return [
        hit
        for hit in hits
        if hit.identity >= floor_id
        and hit.query_coverage >= min_query_cov
        and hit.target_coverage >= min_target_cov
    ]


def cluster_search_dataset(
    build: str | Path,
    dataset: str,
    threads: int = 4,
    vsearch_bin: str = "vsearch",
    strict_tool_version: bool = False,
    iddef: int = DEFAULT_IDDEF,
    species_id: float = 0.987,
    genus_id: float = 0.945,
    family_id: float = 0.865,
    order_id: float = 0.820,
    class_id: float = 0.785,
    floor_id: float = 0.750,
    min_query_cov: float = 0.80,
    min_target_cov: float = 0.0,
    maxaccepts: int = 50,
    maxrejects: int = 256,
    near_best_delta: float = 0.005,
    strand: str = "plus",
) -> ClusterSearchSummary:
    """Run internal dataset clustering and search centroids against current registry."""
    if strand not in {"plus", "both"}:
        raise ValueError("strand must be plus or both.")

    build_dir = Path(build)
    dataset_dir = _find_dataset_dir(build_dir, dataset)
    oriented_fasta = dataset_dir / "sina.oriented.fa"
    input_records = read_fasta(oriented_fasta)
    clusters_dir = dataset_dir / "internal_clusters"
    clusters_dir.mkdir(parents=True, exist_ok=True)
    detected_version = check_vsearch_version(vsearch_bin)
    if strict_tool_version and not detected_version:
        raise RuntimeError(f"Could not parse VSEARCH version from {vsearch_bin} --version.")

    thresholds = {
        "species": species_id,
        "genus": genus_id,
        "family": family_id,
        "order": order_id,
        "class": class_id,
    }
    commands: list[list[str]] = []
    centroid_counts: dict[str, int] = {}
    for rank, identity in thresholds.items():
        uc_path = clusters_dir / f"{rank}_{identity:.3f}.uc"
        centroids_path = clusters_dir / f"{rank}_{identity:.3f}.centroids.fa"
        if not (uc_path.exists() and centroids_path.exists()):
            try:
                command = run_vsearch_cluster(
                    fasta_path=oriented_fasta,
                    uc_path=uc_path,
                    identity=identity,
                    threads=threads,
                    vsearch_bin=vsearch_bin,
                    iddef=iddef,
                    centroids_path=centroids_path,
                )
                commands.append(command)
            except (FileNotFoundError, OSError, subprocess.CalledProcessError) as exc:
                raise RuntimeError(f"VSEARCH cluster_fast failed for {rank}: {vsearch_bin}") from exc
        memberships = cluster_memberships_from_uc(uc_path, rank_threshold=f"{rank}_{identity:.3f}")
        _write_tsv(
            [membership.__dict__ for membership in memberships],
            clusters_dir / f"{rank}_{identity:.3f}.members.tsv",
            CLUSTER_MEMBER_FIELDS,
        )
        centroid_counts[rank] = len(read_fasta(centroids_path)) if centroids_path.exists() else 0

    registry_reps = build_current_representatives(build_dir, current_dataset=dataset)
    reps_fasta = build_dir / "registry" / "current_representatives.fa"
    species_centroids = clusters_dir / f"species_{species_id:.3f}.centroids.fa"
    raw_userout = dataset_dir / "vs_registry.raw.tsv"
    vs_registry = dataset_dir / "vs_registry.tsv"
    vs_filtered = dataset_dir / "vs_registry.filtered.tsv"
    raw_hits: list[RegistryHit] = []

    if registry_reps > 0 and centroid_counts.get("species", 0) > 0:
        search_command = build_usearch_global_command(
            query_fasta=species_centroids,
            db_fasta=reps_fasta,
            userout_path=raw_userout,
            identity=floor_id,
            iddef=iddef,
            maxaccepts=maxaccepts,
            maxrejects=maxrejects,
            strand=strand,
            vsearch_bin=vsearch_bin,
        )
        try:
            run_vsearch_search(search_command)
            commands.append(search_command)
        except (FileNotFoundError, OSError, subprocess.CalledProcessError) as exc:
            raise RuntimeError(f"VSEARCH registry search failed: {vsearch_bin}") from exc
        raw_hits = parse_registry_hits(raw_userout)
    else:
        raw_userout.write_text("", encoding="utf-8")

    filtered_hits = filter_registry_hits(
        raw_hits,
        floor_id=floor_id,
        min_query_cov=min_query_cov,
        min_target_cov=min_target_cov,
    )
    _write_tsv([_hit_row(hit) for hit in raw_hits], vs_registry, REGISTRY_HIT_FIELDS)
    _write_tsv([_hit_row(hit) for hit in filtered_hits], vs_filtered, REGISTRY_HIT_FIELDS)

    dataset_prefix = _dataset_prefix(build_dir, dataset)
    summary_row = {
        "dataset": dataset,
        "prefix": dataset_prefix,
        "input_sequences": str(len(input_records)),
        "species_centroids": str(centroid_counts.get("species", 0)),
        "genus_centroids": str(centroid_counts.get("genus", 0)),
        "family_centroids": str(centroid_counts.get("family", 0)),
        "order_centroids": str(centroid_counts.get("order", 0)),
        "class_centroids": str(centroid_counts.get("class", 0)),
        "registry_representatives": str(registry_reps),
        "registry_hits_raw": str(len(raw_hits)),
        "registry_hits_filtered": str(len(filtered_hits)),
        "iddef": str(iddef),
        "min_query_cov": str(min_query_cov),
        "min_target_cov": str(min_target_cov),
        "near_best_delta": str(near_best_delta),
        "vsearch_version": detected_version,
    }
    _write_tsv([summary_row], dataset_dir / "cluster_search_summary.tsv", SUMMARY_FIELDS)
    _write_tool_version(
        dataset_dir / "tool_versions.tsv",
        detected_version=detected_version,
        status="ok" if detected_version else "warning_unparsed",
        command=" ; ".join(" ".join(command) for command in commands),
    )
    return ClusterSearchSummary(
        dataset=dataset,
        dataset_dir=dataset_dir,
        registry_hits_filtered=len(filtered_hits),
    )


def build_current_representatives(build_dir: Path, current_dataset: str | None = None) -> int:
    """Build a simple current registry representative FASTA and taxonomy TSV."""
    registry_dir = build_dir / "registry"
    registry_dir.mkdir(parents=True, exist_ok=True)
    reps_fasta = registry_dir / "current_representatives.fa"
    reps_tax = registry_dir / "current_representatives.tax.tsv"
    if reps_fasta.exists():
        records = read_fasta(reps_fasta)
        if not reps_tax.exists():
            _write_tsv(
                [{"seq_id": record.seq_id, "taxonomy": "", "source": "existing"} for record in records],
                reps_tax,
                ["seq_id", "taxonomy", "source"],
            )
        return len(records)

    records: list[FastaRecord] = []
    tax_rows: list[dict[str, str]] = []
    _append_representatives_from_fasta(
        build_dir / "silva" / "silva_named_backbone.fa",
        records,
        tax_rows,
        source="silva_named_backbone",
    )
    _append_representatives_from_fasta(
        build_dir / "silva" / "silva_unresolved.resolved.fa",
        records,
        tax_rows,
        source="silva_unresolved_resolved",
    )
    datasets_dir = build_dir / "datasets"
    if datasets_dir.exists():
        for dataset_dir in sorted(path for path in datasets_dir.iterdir() if path.is_dir()):
            dataset_name = dataset_dir.name.split("_", maxsplit=1)[-1]
            if current_dataset is not None and dataset_name == current_dataset:
                continue
            _append_representatives_from_fasta(
                dataset_dir / "sina.oriented.fa",
                records,
                tax_rows,
                source=f"dataset:{dataset_name}",
            )
    write_fasta(records, reps_fasta)
    _write_tsv(tax_rows, reps_tax, ["seq_id", "taxonomy", "source"])
    return len(records)


def write_singleton_uc(records: list[str], path: str | Path) -> None:
    """Write a minimal .uc file assigning every sequence to its own cluster."""
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        for index, seq_id in enumerate(records):
            handle.write(f"S\t{index}\t*\t*\t*\t*\t*\t*\t{seq_id}\t*\n")


def _append_representatives_from_fasta(
    path: Path,
    records: list[FastaRecord],
    tax_rows: list[dict[str, str]],
    source: str,
) -> None:
    if not path.exists():
        return
    for record in read_fasta(path):
        records.append(FastaRecord(seq_id=record.seq_id, header=record.seq_id, sequence=record.sequence))
        taxonomy = ""
        if " " in record.header:
            taxonomy = record.header.split(maxsplit=1)[1]
        tax_rows.append({"seq_id": record.seq_id, "taxonomy": taxonomy, "source": source})


def _find_dataset_dir(build_dir: Path, dataset: str) -> Path:
    datasets_dir = build_dir / "datasets"
    if not datasets_dir.exists():
        raise FileNotFoundError(f"Dataset directory root not found: {datasets_dir}")
    matches = sorted(
        path
        for path in datasets_dir.iterdir()
        if path.is_dir() and path.name.split("_", maxsplit=1)[-1] == dataset
    )
    if not matches:
        raise FileNotFoundError(f"Prepared dataset not found: {dataset}")
    if len(matches) > 1:
        raise ValueError(f"Multiple prepared dataset directories match {dataset!r}.")
    return matches[0]


def _dataset_prefix(build_dir: Path, dataset: str) -> str:
    for row in _read_tsv(build_dir / "registry" / "dataset_registry.tsv"):
        if row.get("dataset_name") == dataset:
            return row.get("prefix", "")
    return ""


def _hit_row(hit: RegistryHit) -> dict[str, str]:
    return {
        "query": hit.query,
        "target": hit.target,
        "identity": f"{hit.identity:.3f}",
        "alnlen": str(hit.alnlen),
        "qlo": str(hit.qlo),
        "qhi": str(hit.qhi),
        "tlo": str(hit.tlo),
        "thi": str(hit.thi),
        "ql": str(hit.ql),
        "tl": str(hit.tl),
        "bits": f"{hit.bits:.3f}",
        "query_coverage": f"{hit.query_coverage:.6f}",
        "target_coverage": f"{hit.target_coverage:.6f}",
    }


def _write_tool_version(
    tool_versions_path: Path,
    detected_version: str,
    status: str,
    command: str,
) -> None:
    rows = [
        row
        for row in _read_tsv(tool_versions_path)
        if row.get("tool") != "vsearch"
    ]
    rows.append(
        {
            "tool": "vsearch",
            "expected_version": EXPECTED_VSEARCH_VERSION,
            "detected_version": detected_version,
            "status": status,
            "command": command,
        }
    )
    _write_tsv(rows, tool_versions_path, TOOL_VERSION_FIELDS)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=chr(9)))


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
