from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from .config import get_reference_path, get_software
from .external import make_blast_db, run_command
from .logging import step, success, warning
from .threads import validate_threads
from .utils import ensure_dir, ensure_file, read_fasta


BLAST_HSP_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
    "qseq",
    "sseq",
]


CANDIDATE_FIELDS = [
    "candidate_id",
    "cluster_id",
    "query_id",
    "subject_id",
    "reference_status",
    "strand",
    "query_intron_start",
    "query_intron_end",
    "query_intron_length",
    "left_qstart",
    "left_qend",
    "right_qstart",
    "right_qend",
    "left_sstart",
    "left_send",
    "right_sstart",
    "right_send",
    "ref_gap",
    "left_hsp_length",
    "right_hsp_length",
    "left_identity",
    "right_identity",
    "mean_identity",
    "left_bitscore",
    "right_bitscore",
]


SUPPORT_FIELDS = [
    "cluster_id",
    "query_id",
    "query_intron_start_median",
    "query_intron_end_median",
    "query_intron_length_median",
    "raw_ref_support",
    "clean_ref_support",
    "reference_intron_support",
    "clean_supporting_species",
    "clean_supporting_genera",
    "clean_supporting_families",
    "clean_supporting_orders",
    "clean_supporting_classes",
    "clean_supporting_phyla",
    "raw_supporting_species",
    "raw_supporting_genera",
    "raw_supporting_families",
    "raw_supporting_orders",
    "raw_supporting_classes",
    "raw_supporting_phyla",
    "median_ref_gap",
    "median_mean_identity",
    "best_mean_identity",
    "max_clean_family_fraction",
    "max_raw_family_fraction",
]


CONFIRMED_FIELDS = [
    "cluster_id",
    "query_id",
    "source_label",
    "query_intron_start",
    "query_intron_end",
    "query_intron_length",
    "raw_ref_support",
    "clean_ref_support",
    "reference_intron_support",
    "clean_supporting_species",
    "clean_supporting_genera",
    "clean_supporting_families",
    "confirmation_subject",
    "confirmation_identity",
    "confirmation_qcov",
    "confirmation_bitscore",
    "intron_status",
]


FAILED_CONFIRMATION_FIELDS = [
    "cluster_id",
    "query_id",
    "source_label",
    "query_intron_start",
    "query_intron_end",
    "query_intron_length",
    "raw_ref_support",
    "clean_ref_support",
    "reference_intron_support",
    "clean_supporting_species",
    "best_confirmation_subject",
    "best_confirmation_identity",
    "best_confirmation_qcov",
    "reason",
]


VERSION_MAP_FIELDS = [
    "analysis_id",
    "original_id",
    "status",
    "removed_count",
    "removed_regions",
    "source_label",
]


@dataclass(frozen=True)
class HSP:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart_raw: int
    qend_raw: int
    sstart_raw: int
    send_raw: int
    evalue: str
    bitscore: float
    qlen: int
    slen: int
    qseq: str
    sseq: str

    @property
    def qstart(self) -> int:
        return min(self.qstart_raw, self.qend_raw)

    @property
    def qend(self) -> int:
        return max(self.qstart_raw, self.qend_raw)

    @property
    def smin(self) -> int:
        return min(self.sstart_raw, self.send_raw)

    @property
    def smax(self) -> int:
        return max(self.sstart_raw, self.send_raw)

    @property
    def strand(self) -> str:
        return "plus" if self.sstart_raw <= self.send_raw else "minus"

    @property
    def qcov(self) -> float:
        if self.qlen <= 0:
            return 0.0
        return (self.qend - self.qstart + 1) / self.qlen


@dataclass(frozen=True)
class ReferenceStatus:
    is_reference_intron: bool
    species: str = ""
    genus: str = ""
    family: str = ""
    order: str = ""
    class_name: str = ""
    phylum: str = ""


@dataclass(frozen=True)
class IntronCandidate:
    candidate_id: str
    cluster_id: str
    query_id: str
    subject_id: str
    strand: str
    query_intron_start: int
    query_intron_end: int
    query_intron_length: int
    ref_gap: int
    left_hsp: HSP
    right_hsp: HSP

    @property
    def mean_identity(self) -> float:
        return (self.left_hsp.pident + self.right_hsp.pident) / 2.0

    @property
    def mean_bitscore(self) -> float:
        return (self.left_hsp.bitscore + self.right_hsp.bitscore) / 2.0


@dataclass(frozen=True)
class ClusterSupport:
    cluster_id: str
    query_id: str
    candidates: List[IntronCandidate]
    raw_refs: Set[str]
    clean_refs: Set[str]
    reference_intron_refs: Set[str]
    raw_species: Set[str]
    raw_genera: Set[str]
    raw_families: Set[str]
    raw_orders: Set[str]
    raw_classes: Set[str]
    raw_phyla: Set[str]
    clean_species: Set[str]
    clean_genera: Set[str]
    clean_families: Set[str]
    clean_orders: Set[str]
    clean_classes: Set[str]
    clean_phyla: Set[str]
    max_clean_family_fraction: float
    max_raw_family_fraction: float

    @property
    def best_candidate(self) -> IntronCandidate:
        return max(
            self.candidates,
            key=lambda item: (
                item.mean_identity,
                item.mean_bitscore,
                -abs(item.ref_gap),
                item.query_intron_length,
            ),
        )

    @property
    def median_start(self) -> int:
        return int(round(median([c.query_intron_start for c in self.candidates])))

    @property
    def median_end(self) -> int:
        return int(round(median([c.query_intron_end for c in self.candidates])))

    @property
    def median_length(self) -> int:
        return int(round(median([c.query_intron_length for c in self.candidates])))

    @property
    def median_ref_gap(self) -> float:
        return float(median([c.ref_gap for c in self.candidates]))

    @property
    def median_mean_identity(self) -> float:
        return float(median([c.mean_identity for c in self.candidates]))

    @property
    def best_mean_identity(self) -> float:
        return max(c.mean_identity for c in self.candidates)


@dataclass(frozen=True)
class ConfirmationResult:
    cluster_id: str
    query_id: str
    confirmation_query_id: str
    subject_id: str
    identity: float
    qcov: float
    bitscore: float
    ok: bool


def safe_float(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return 0.0


def safe_int(value: str) -> int:
    try:
        return int(float(value))
    except ValueError:
        return 0


def is_truthy(value: object) -> bool:
    return str(value).strip().lower() in {"yes", "y", "true", "1", "putative_intron", "intron"}


def fasta_to_dict(path: str | Path) -> Dict[str, Tuple[str, str]]:
    records: Dict[str, Tuple[str, str]] = {}

    for header, seq in read_fasta(path):
        seq_id = header.split()[0].split(";")[0]
        records[seq_id] = (header, "".join(str(seq).split()).upper())

    return records


def write_fasta_record(handle, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")

    for i in range(0, len(sequence), 80):
        handle.write(sequence[i:i + 80] + "\n")


def write_tsv(path: str | Path, fields: Sequence[str], rows: Iterable[Dict[str, object]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t")
        writer.writeheader()

        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def blastdb_exists(prefix: str | Path) -> bool:
    prefix = Path(prefix)
    candidates = [
        Path(str(prefix) + suffix)
        for suffix in (
            ".nhr", ".nin", ".nsq", ".ndb", ".not", ".ntf", ".nto",
            ".00.nhr", ".00.nin", ".00.nsq", ".00.ndb", ".00.not",
            ".00.ntf", ".00.nto",
        )
    ]
    return any(path.exists() for path in candidates)


def prepare_blast_reference(db: str | Path, outdir: str | Path, *, makeblastdb: str | Path) -> str:
    db = Path(db)

    if blastdb_exists(db):
        return str(db)

    if db.exists() and db.is_file():
        prefix = Path(outdir) / "blast" / "refdb"
        make_blast_db(
            input_fasta=db,
            db_prefix=prefix,
            makeblastdb=makeblastdb,
            title="autotax2_intron_reference",
            dbtype="nucl",
        )
        return str(prefix)

    return str(db)


def run_blastn_hsps(
    query_fasta: str | Path,
    db_prefix: str | Path,
    output_tsv: str | Path,
    *,
    blastn: str | Path,
    threads: int,
    perc_identity: float,
    max_hsps: int,
    max_target_seqs: int,
    strand: str,
    evalue: float = 1e-20,
) -> str:
    threads = validate_threads(threads)

    if strand not in {"both", "plus", "minus"}:
        raise ValueError("strand must be one of: both, plus, minus.")

    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    outfmt = (
        "6 qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore qlen slen qseq sseq"
    )

    command = [
        blastn,
        "-query", query_fasta,
        "-db", db_prefix,
        "-task", "blastn",
        "-perc_identity", str(perc_identity),
        "-max_hsps", str(max_hsps),
        "-max_target_seqs", str(max_target_seqs),
        "-evalue", str(evalue),
        "-dust", "no",
        "-strand", strand,
        "-num_threads", str(threads),
        "-outfmt", outfmt,
        "-out", output_tsv,
    ]

    run_command(command)
    return str(output_tsv)


def parse_blast_hsps(path: str | Path) -> List[HSP]:
    hsps: List[HSP] = []
    path = Path(path)

    if not path.exists() or path.stat().st_size == 0:
        return hsps

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue

            values = line.rstrip("\n").split("\t")

            if len(values) < len(BLAST_HSP_FIELDS):
                continue

            row = dict(zip(BLAST_HSP_FIELDS, values))
            hsps.append(
                HSP(
                    qseqid=row["qseqid"],
                    sseqid=row["sseqid"],
                    pident=safe_float(row["pident"]),
                    length=safe_int(row["length"]),
                    mismatch=safe_int(row["mismatch"]),
                    gapopen=safe_int(row["gapopen"]),
                    qstart_raw=safe_int(row["qstart"]),
                    qend_raw=safe_int(row["qend"]),
                    sstart_raw=safe_int(row["sstart"]),
                    send_raw=safe_int(row["send"]),
                    evalue=row["evalue"],
                    bitscore=safe_float(row["bitscore"]),
                    qlen=safe_int(row["qlen"]),
                    slen=safe_int(row["slen"]),
                    qseq=row["qseq"],
                    sseq=row["sseq"],
                )
            )

    return hsps


def reference_gap(left: HSP, right: HSP) -> int:
    if left.strand == "plus":
        return right.smin - left.smax - 1
    return left.smin - right.smax - 1


def detect_candidates_from_hsps(
    hsps: List[HSP],
    *,
    min_intron_len: int,
    min_flank_len: int,
    hsp_identity: float,
    max_ref_gap: int,
    ref_overlap_tolerance: int,
) -> List[IntronCandidate]:
    grouped: Dict[Tuple[str, str, str], List[HSP]] = defaultdict(list)

    for hsp in hsps:
        if hsp.pident < hsp_identity:
            continue
        if hsp.length < min_flank_len:
            continue
        grouped[(hsp.qseqid, hsp.sseqid, hsp.strand)].append(hsp)

    candidates: List[IntronCandidate] = []
    counter = 0

    for (_query_id, _subject_id, _strand), group in grouped.items():
        group = sorted(group, key=lambda item: (item.qstart, item.qend, -item.bitscore))

        for i, left in enumerate(group):
            for right in group[i + 1:]:
                if right.qstart <= left.qend:
                    continue

                query_gap = right.qstart - left.qend - 1
                if query_gap < min_intron_len:
                    continue

                ref_gap = reference_gap(left, right)
                if ref_gap < -ref_overlap_tolerance:
                    continue
                if ref_gap > max_ref_gap:
                    continue

                counter += 1
                candidates.append(
                    IntronCandidate(
                        candidate_id=f"candidate_{counter:06d}",
                        cluster_id="",
                        query_id=left.qseqid,
                        subject_id=left.sseqid,
                        strand=left.strand,
                        query_intron_start=left.qend + 1,
                        query_intron_end=right.qstart - 1,
                        query_intron_length=query_gap,
                        ref_gap=ref_gap,
                        left_hsp=left,
                        right_hsp=right,
                    )
                )

    return candidates


def cluster_candidates(candidates: List[IntronCandidate], *, coordinate_tolerance: int) -> List[IntronCandidate]:
    by_query: Dict[str, List[IntronCandidate]] = defaultdict(list)

    for candidate in candidates:
        by_query[candidate.query_id].append(candidate)

    clustered: List[IntronCandidate] = []
    cluster_counter = 0

    for query_id, query_candidates in by_query.items():
        groups: List[List[IntronCandidate]] = []

        for candidate in sorted(query_candidates, key=lambda item: (item.query_intron_start, item.query_intron_end, item.subject_id)):
            placed = False

            for group in groups:
                center_start = int(round(median([item.query_intron_start for item in group])))
                center_end = int(round(median([item.query_intron_end for item in group])))

                if (
                    abs(candidate.query_intron_start - center_start) <= coordinate_tolerance
                    and abs(candidate.query_intron_end - center_end) <= coordinate_tolerance
                ):
                    group.append(candidate)
                    placed = True
                    break

            if not placed:
                groups.append([candidate])

        for group in groups:
            cluster_counter += 1
            cluster_id = f"{query_id}|intron_cluster_{cluster_counter:06d}"

            for candidate in group:
                clustered.append(
                    IntronCandidate(
                        candidate_id=candidate.candidate_id,
                        cluster_id=cluster_id,
                        query_id=candidate.query_id,
                        subject_id=candidate.subject_id,
                        strand=candidate.strand,
                        query_intron_start=candidate.query_intron_start,
                        query_intron_end=candidate.query_intron_end,
                        query_intron_length=candidate.query_intron_length,
                        ref_gap=candidate.ref_gap,
                        left_hsp=candidate.left_hsp,
                        right_hsp=candidate.right_hsp,
                    )
                )

    return clustered


def read_reference_blacklist(path: str | Path | None) -> Set[str]:
    if path is None:
        return set()

    path = Path(path)
    if not path.exists():
        return set()

    blacklist: Set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            text = line.strip()
            if text and not text.startswith("#"):
                blacklist.add(text.split("\t")[0])
    return blacklist


def read_species_taxonomy(path: str | Path | None, blacklist: Set[str]) -> Dict[str, ReferenceStatus]:
    if path is None:
        return {seq_id: ReferenceStatus(is_reference_intron=True) for seq_id in blacklist}

    path = Path(path)
    if not path.exists():
        return {seq_id: ReferenceStatus(is_reference_intron=True) for seq_id in blacklist}

    taxonomy: Dict[str, ReferenceStatus] = {}

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            rep_id = row.get("rep_id") or row.get("ID") or row.get("original_id")
            if not rep_id:
                continue

            ref_has_intron = (
                is_truthy(row.get("ref_has_putative_intron", ""))
                or rep_id in blacklist
            )

            status = ReferenceStatus(
                is_reference_intron=ref_has_intron,
                species=row.get("species_key") or row.get("Species", "") or "",
                genus=row.get("Genus", "") or "",
                family=row.get("Family", "") or "",
                order=row.get("Order", "") or "",
                class_name=row.get("Class", "") or "",
                phylum=row.get("Phylum", "") or "",
            )

            taxonomy[rep_id] = status

            original_id = row.get("original_id")
            if original_id:
                taxonomy[original_id] = status

    for seq_id in blacklist:
        taxonomy.setdefault(seq_id, ReferenceStatus(is_reference_intron=True))

    return taxonomy


def ref_status(taxonomy: Dict[str, ReferenceStatus], subject_id: str, blacklist: Set[str]) -> ReferenceStatus:
    if subject_id in blacklist:
        current = taxonomy.get(subject_id)
        if current is None:
            return ReferenceStatus(is_reference_intron=True)
        return ReferenceStatus(
            is_reference_intron=True,
            species=current.species,
            genus=current.genus,
            family=current.family,
            order=current.order,
            class_name=current.class_name,
            phylum=current.phylum,
        )
    return taxonomy.get(subject_id, ReferenceStatus(is_reference_intron=False))


def add_nonempty(target: Set[str], value: str) -> None:
    if value:
        target.add(value)


def max_fraction(values: Iterable[str]) -> float:
    counts: Dict[str, int] = defaultdict(int)
    total = 0

    for value in values:
        if not value:
            continue
        counts[value] += 1
        total += 1

    if total == 0:
        return 0.0
    return max(counts.values()) / total


def build_cluster_support(
    candidates: List[IntronCandidate],
    taxonomy: Dict[str, ReferenceStatus],
    blacklist: Set[str],
) -> List[ClusterSupport]:
    by_cluster: Dict[str, List[IntronCandidate]] = defaultdict(list)

    for candidate in candidates:
        by_cluster[candidate.cluster_id].append(candidate)

    supports: List[ClusterSupport] = []

    for cluster_id, group in by_cluster.items():
        raw_refs: Set[str] = set()
        clean_refs: Set[str] = set()
        reference_intron_refs: Set[str] = set()

        raw_species: Set[str] = set()
        raw_genera: Set[str] = set()
        raw_families: Set[str] = set()
        raw_orders: Set[str] = set()
        raw_classes: Set[str] = set()
        raw_phyla: Set[str] = set()

        clean_species: Set[str] = set()
        clean_genera: Set[str] = set()
        clean_families: Set[str] = set()
        clean_orders: Set[str] = set()
        clean_classes: Set[str] = set()
        clean_phyla: Set[str] = set()

        clean_family_values: List[str] = []
        raw_family_values: List[str] = []

        for item in group:
            status = ref_status(taxonomy, item.subject_id, blacklist)
            raw_refs.add(item.subject_id)
            add_nonempty(raw_species, status.species)
            add_nonempty(raw_genera, status.genus)
            add_nonempty(raw_families, status.family)
            add_nonempty(raw_orders, status.order)
            add_nonempty(raw_classes, status.class_name)
            add_nonempty(raw_phyla, status.phylum)
            raw_family_values.append(status.family)

            if status.is_reference_intron:
                reference_intron_refs.add(item.subject_id)
            else:
                clean_refs.add(item.subject_id)
                add_nonempty(clean_species, status.species)
                add_nonempty(clean_genera, status.genus)
                add_nonempty(clean_families, status.family)
                add_nonempty(clean_orders, status.order)
                add_nonempty(clean_classes, status.class_name)
                add_nonempty(clean_phyla, status.phylum)
                clean_family_values.append(status.family)

        supports.append(
            ClusterSupport(
                cluster_id=cluster_id,
                query_id=group[0].query_id,
                candidates=group,
                raw_refs=raw_refs,
                clean_refs=clean_refs,
                reference_intron_refs=reference_intron_refs,
                raw_species=raw_species,
                raw_genera=raw_genera,
                raw_families=raw_families,
                raw_orders=raw_orders,
                raw_classes=raw_classes,
                raw_phyla=raw_phyla,
                clean_species=clean_species,
                clean_genera=clean_genera,
                clean_families=clean_families,
                clean_orders=clean_orders,
                clean_classes=clean_classes,
                clean_phyla=clean_phyla,
                max_clean_family_fraction=max_fraction(clean_family_values),
                max_raw_family_fraction=max_fraction(raw_family_values),
            )
        )

    supports.sort(key=lambda item: (item.query_id, item.median_start, item.median_end, -len(item.clean_refs), -len(item.raw_refs)))
    return supports


def candidate_to_row(candidate: IntronCandidate, taxonomy: Dict[str, ReferenceStatus], blacklist: Set[str]) -> Dict[str, object]:
    status = ref_status(taxonomy, candidate.subject_id, blacklist)
    return {
        "candidate_id": candidate.candidate_id,
        "cluster_id": candidate.cluster_id,
        "query_id": candidate.query_id,
        "subject_id": candidate.subject_id,
        "reference_status": "reference_intron" if status.is_reference_intron else "clean_reference",
        "strand": candidate.strand,
        "query_intron_start": candidate.query_intron_start,
        "query_intron_end": candidate.query_intron_end,
        "query_intron_length": candidate.query_intron_length,
        "left_qstart": candidate.left_hsp.qstart,
        "left_qend": candidate.left_hsp.qend,
        "right_qstart": candidate.right_hsp.qstart,
        "right_qend": candidate.right_hsp.qend,
        "left_sstart": candidate.left_hsp.sstart_raw,
        "left_send": candidate.left_hsp.send_raw,
        "right_sstart": candidate.right_hsp.sstart_raw,
        "right_send": candidate.right_hsp.send_raw,
        "ref_gap": candidate.ref_gap,
        "left_hsp_length": candidate.left_hsp.length,
        "right_hsp_length": candidate.right_hsp.length,
        "left_identity": f"{candidate.left_hsp.pident:.3f}",
        "right_identity": f"{candidate.right_hsp.pident:.3f}",
        "mean_identity": f"{candidate.mean_identity:.3f}",
        "left_bitscore": f"{candidate.left_hsp.bitscore:.3f}",
        "right_bitscore": f"{candidate.right_hsp.bitscore:.3f}",
    }


def support_to_row(support: ClusterSupport) -> Dict[str, object]:
    return {
        "cluster_id": support.cluster_id,
        "query_id": support.query_id,
        "query_intron_start_median": support.median_start,
        "query_intron_end_median": support.median_end,
        "query_intron_length_median": support.median_length,
        "raw_ref_support": len(support.raw_refs),
        "clean_ref_support": len(support.clean_refs),
        "reference_intron_support": len(support.reference_intron_refs),
        "clean_supporting_species": len(support.clean_species),
        "clean_supporting_genera": len(support.clean_genera),
        "clean_supporting_families": len(support.clean_families),
        "clean_supporting_orders": len(support.clean_orders),
        "clean_supporting_classes": len(support.clean_classes),
        "clean_supporting_phyla": len(support.clean_phyla),
        "raw_supporting_species": len(support.raw_species),
        "raw_supporting_genera": len(support.raw_genera),
        "raw_supporting_families": len(support.raw_families),
        "raw_supporting_orders": len(support.raw_orders),
        "raw_supporting_classes": len(support.raw_classes),
        "raw_supporting_phyla": len(support.raw_phyla),
        "median_ref_gap": f"{support.median_ref_gap:.3f}",
        "median_mean_identity": f"{support.median_mean_identity:.3f}",
        "best_mean_identity": f"{support.best_mean_identity:.3f}",
        "max_clean_family_fraction": f"{support.max_clean_family_fraction:.6f}",
        "max_raw_family_fraction": f"{support.max_raw_family_fraction:.6f}",
    }


def write_confirmation_fasta(
    supports: List[ClusterSupport],
    query_records: Dict[str, Tuple[str, str]],
    output_fasta: str | Path,
) -> Dict[str, Tuple[ClusterSupport, str]]:
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    confirmation_map: Dict[str, Tuple[ClusterSupport, str]] = {}

    with output_fasta.open("w", encoding="utf-8") as handle:
        for support in supports:
            best = support.best_candidate
            if best.query_id not in query_records:
                continue

            _header, seq = query_records[best.query_id]
            start0 = best.query_intron_start - 1
            end0 = best.query_intron_end
            trimmed = seq[:start0] + seq[end0:]
            confirmation_id = f"{best.query_id}|{support.cluster_id}|confirmation"
            write_fasta_record(handle, confirmation_id, trimmed)
            confirmation_map[confirmation_id] = (support, trimmed)

    return confirmation_map


def parse_best_confirmation_hits(
    confirmation_hsps: List[HSP],
    *,
    confirm_id_threshold: float,
    confirm_qcov_threshold: float,
) -> Dict[str, ConfirmationResult]:
    best_by_query: Dict[str, HSP] = {}

    for hsp in confirmation_hsps:
        current = best_by_query.get(hsp.qseqid)
        if current is None or (hsp.pident, hsp.qcov, hsp.bitscore) > (current.pident, current.qcov, current.bitscore):
            best_by_query[hsp.qseqid] = hsp

    results: Dict[str, ConfirmationResult] = {}

    for confirmation_query, hsp in best_by_query.items():
        parts = confirmation_query.split("|")
        original_query = parts[0]
        cluster_id = "|".join(parts[1:-1]) if len(parts) >= 3 else confirmation_query
        identity_fraction = hsp.pident / 100.0
        ok = identity_fraction >= confirm_id_threshold and hsp.qcov >= confirm_qcov_threshold
        results[confirmation_query] = ConfirmationResult(
            cluster_id=cluster_id,
            query_id=original_query,
            confirmation_query_id=confirmation_query,
            subject_id=hsp.sseqid,
            identity=identity_fraction,
            qcov=hsp.qcov,
            bitscore=hsp.bitscore,
            ok=ok,
        )

    return results


def intervals_overlap(a: Tuple[int, int], b: Tuple[int, int]) -> bool:
    return max(a[0], b[0]) <= min(a[1], b[1])


def classify_confirmed_status(
    support: ClusterSupport,
    *,
    min_support_refs: int,
    min_support_species: int,
    min_support_genera: int,
) -> Optional[str]:
    clean_refs = len(support.clean_refs)
    clean_species = len(support.clean_species)
    clean_genera = len(support.clean_genera)
    ref_intron_refs = len(support.reference_intron_refs)

    if clean_refs == 0:
        return None

    if clean_refs >= min_support_refs and clean_species >= min_support_species and clean_genera >= min_support_genera:
        if ref_intron_refs > 0:
            return "shared_reference_supported_intron"
        return "query_specific_intron"

    return "lineage_supported_intron"


def choose_confirmed_clusters(
    supports: List[ClusterSupport],
    confirmation_results: Dict[str, ConfirmationResult],
    *,
    min_support_refs: int,
    min_support_species: int,
    min_support_genera: int,
) -> Tuple[Dict[str, ConfirmationResult], Dict[str, str], List[Dict[str, object]]]:
    confirmed: Dict[str, ConfirmationResult] = {}
    status_by_cluster: Dict[str, str] = {}
    failed_rows: List[Dict[str, object]] = []
    confirmation_by_cluster = {result.cluster_id: result for result in confirmation_results.values()}

    for support in supports:
        result = confirmation_by_cluster.get(support.cluster_id)
        best = support.best_candidate

        base_row = {
            "cluster_id": support.cluster_id,
            "query_id": support.query_id,
            "source_label": "",
            "query_intron_start": best.query_intron_start,
            "query_intron_end": best.query_intron_end,
            "query_intron_length": best.query_intron_length,
            "raw_ref_support": len(support.raw_refs),
            "clean_ref_support": len(support.clean_refs),
            "reference_intron_support": len(support.reference_intron_refs),
            "clean_supporting_species": len(support.clean_species),
        }

        if result is None:
            failed_rows.append({**base_row, "best_confirmation_subject": "", "best_confirmation_identity": "", "best_confirmation_qcov": "", "reason": "no_confirmation_hit"})
            continue

        if not result.ok:
            failed_rows.append({
                **base_row,
                "best_confirmation_subject": result.subject_id,
                "best_confirmation_identity": f"{result.identity:.6f}",
                "best_confirmation_qcov": f"{result.qcov:.6f}",
                "reason": "confirmation_below_threshold",
            })
            continue

        status = classify_confirmed_status(
            support,
            min_support_refs=min_support_refs,
            min_support_species=min_support_species,
            min_support_genera=min_support_genera,
        )

        if status is None:
            reason = "only_reference_intron_support" if support.reference_intron_refs else "no_clean_reference_support"
            failed_rows.append({
                **base_row,
                "best_confirmation_subject": result.subject_id,
                "best_confirmation_identity": f"{result.identity:.6f}",
                "best_confirmation_qcov": f"{result.qcov:.6f}",
                "reason": reason,
            })
            continue

        confirmed[support.cluster_id] = result
        status_by_cluster[support.cluster_id] = status

    support_by_cluster = {support.cluster_id: support for support in supports}
    by_query: Dict[str, List[ClusterSupport]] = defaultdict(list)

    for cluster_id in confirmed:
        support = support_by_cluster.get(cluster_id)
        if support is not None:
            by_query[support.query_id].append(support)

    keep: Dict[str, ConfirmationResult] = {}
    keep_status: Dict[str, str] = {}

    for _query_id, query_supports in by_query.items():
        chosen_intervals: List[Tuple[int, int]] = []
        ranked = sorted(
            query_supports,
            key=lambda support: (
                -len(support.clean_refs),
                -len(support.clean_species),
                -len(support.raw_refs),
                -support.best_mean_identity,
                support.median_start,
            ),
        )

        for support in ranked:
            interval = (support.median_start, support.median_end)
            if any(intervals_overlap(interval, existing) for existing in chosen_intervals):
                continue
            chosen_intervals.append(interval)
            keep[support.cluster_id] = confirmed[support.cluster_id]
            keep_status[support.cluster_id] = status_by_cluster[support.cluster_id]

    return keep, keep_status, failed_rows


def write_final_outputs(
    query_records: Dict[str, Tuple[str, str]],
    supports: List[ClusterSupport],
    confirmed: Dict[str, ConfirmationResult],
    status_by_cluster: Dict[str, str],
    outdir: str | Path,
    *,
    source_label: str,
) -> None:
    confirmed_dir = Path(outdir) / "confirmed"
    confirmed_dir.mkdir(parents=True, exist_ok=True)

    analysis_fasta = confirmed_dir / "analysis_sequences.fasta"
    introns_fasta = confirmed_dir / "introns.fasta"
    removed_fasta = confirmed_dir / "removed_insertions.fasta"
    version_map = confirmed_dir / "sequence_version_map.tsv"

    support_by_cluster = {support.cluster_id: support for support in supports}
    confirmed_by_query: Dict[str, List[ClusterSupport]] = defaultdict(list)

    for cluster_id in confirmed:
        support = support_by_cluster.get(cluster_id)
        if support is not None:
            confirmed_by_query[support.query_id].append(support)

    confirmed_rows: List[Dict[str, object]] = []
    version_rows: List[Dict[str, object]] = []

    with analysis_fasta.open("w", encoding="utf-8") as analysis_handle, introns_fasta.open("w", encoding="utf-8") as intron_handle, removed_fasta.open("w", encoding="utf-8") as removed_handle:
        for query_id, (header, seq) in query_records.items():
            query_supports = sorted(confirmed_by_query.get(query_id, []), key=lambda support: support.median_start)
            removed_regions: List[str] = []
            trimmed_parts: List[str] = []
            cursor = 1

            for support in query_supports:
                result = confirmed[support.cluster_id]
                status = status_by_cluster[support.cluster_id]
                start = support.median_start
                end = support.median_end

                if start < cursor:
                    continue

                trimmed_parts.append(seq[cursor - 1:start - 1])
                intron_seq = seq[start - 1:end]
                cursor = end + 1

                region = f"{start}-{end}"
                removed_regions.append(region)

                intron_header = (
                    f"{query_id}|cluster={support.cluster_id}|region={region}|"
                    f"len={len(intron_seq)}|status={status}|clean_refs={len(support.clean_refs)}|"
                    f"ref_intron_refs={len(support.reference_intron_refs)}"
                )
                write_fasta_record(intron_handle, intron_header, intron_seq)
                write_fasta_record(removed_handle, intron_header, intron_seq)

                confirmed_rows.append({
                    "cluster_id": support.cluster_id,
                    "query_id": query_id,
                    "source_label": source_label,
                    "query_intron_start": start,
                    "query_intron_end": end,
                    "query_intron_length": len(intron_seq),
                    "raw_ref_support": len(support.raw_refs),
                    "clean_ref_support": len(support.clean_refs),
                    "reference_intron_support": len(support.reference_intron_refs),
                    "clean_supporting_species": len(support.clean_species),
                    "clean_supporting_genera": len(support.clean_genera),
                    "clean_supporting_families": len(support.clean_families),
                    "confirmation_subject": result.subject_id,
                    "confirmation_identity": f"{result.identity:.6f}",
                    "confirmation_qcov": f"{result.qcov:.6f}",
                    "confirmation_bitscore": f"{result.bitscore:.3f}",
                    "intron_status": status,
                })

            trimmed_parts.append(seq[cursor - 1:])
            analysis_seq = "".join(trimmed_parts)
            write_fasta_record(analysis_handle, header, analysis_seq)

            version_rows.append({
                "analysis_id": query_id,
                "original_id": query_id,
                "status": "confirmed_intron_removed" if removed_regions else "unchanged",
                "removed_count": len(removed_regions),
                "removed_regions": ",".join(removed_regions),
                "source_label": source_label,
            })

    write_tsv(version_map, VERSION_MAP_FIELDS, version_rows)
    write_tsv(confirmed_dir / "confirmed_introns.tsv", CONFIRMED_FIELDS, confirmed_rows)


def write_no_candidate_outputs(query_records: Dict[str, Tuple[str, str]], outdir: str | Path, source_label: str) -> None:
    confirmed_dir = Path(outdir) / "confirmed"
    candidates_dir = Path(outdir) / "candidates"
    blast_dir = Path(outdir) / "blast"
    confirmed_dir.mkdir(parents=True, exist_ok=True)
    candidates_dir.mkdir(parents=True, exist_ok=True)
    blast_dir.mkdir(parents=True, exist_ok=True)

    with (confirmed_dir / "analysis_sequences.fasta").open("w", encoding="utf-8") as handle:
        for _query_id, (header, seq) in query_records.items():
            write_fasta_record(handle, header, seq)

    version_rows = [
        {
            "analysis_id": query_id,
            "original_id": query_id,
            "status": "unchanged",
            "removed_count": 0,
            "removed_regions": "",
            "source_label": source_label,
        }
        for query_id in query_records
    ]

    write_tsv(confirmed_dir / "sequence_version_map.tsv", VERSION_MAP_FIELDS, version_rows)
    write_tsv(confirmed_dir / "confirmed_introns.tsv", CONFIRMED_FIELDS, [])
    write_tsv(candidates_dir / "failed_confirmation.tsv", FAILED_CONFIRMATION_FIELDS, [])
    (confirmed_dir / "introns.fasta").write_text("", encoding="utf-8")
    (confirmed_dir / "removed_insertions.fasta").write_text("", encoding="utf-8")
    (blast_dir / "confirmation_hits.tsv").write_text("", encoding="utf-8")


def detect_introns(
    input_fasta: str | Path,
    db: str | Path,
    outdir: str | Path,
    source_label: Optional[str] = None,
    *,
    threads: int = 4,
    hsp_identity: float = 90.0,
    min_intron_len: int = 20,
    min_flank_len: int = 100,
    max_ref_gap: int = 10,
    ref_overlap_tolerance: int = 5,
    max_hsps: int = 20,
    maxaccepts: int = 500,
    strand: str = "both",
    confirm_id: float = 0.987,
    confirm_qcov: float = 0.80,
    min_support_refs: int = 3,
    min_support_species: int = 2,
    min_support_genera: int = 1,
    coordinate_tolerance: int = 10,
    species_taxonomy: str | Path | None = None,
    reference_blacklist: str | Path | None = None,
    blastn: str | Path | None = None,
    makeblastdb: str | Path | None = None,
) -> Dict[str, str]:
    """Detect intron-like insertions using local broken-HSP BLASTN evidence.

    Reference sequences can themselves contain introns. AutoTax2 therefore
    separates evidence from clean references and references flagged by the
    reference audit. A query insertion is removed only when confirmation succeeds and
    at least one clean reference supports the broken-HSP event. Reference-intron
    support is still reported, but it is not allowed to be the only evidence for
    automatic removal.
    """

    threads = validate_threads(threads)

    if not 0 < confirm_id <= 1:
        raise ValueError("confirm_id must be > 0 and <= 1.")
    if not 0 < confirm_qcov <= 1:
        raise ValueError("confirm_qcov must be > 0 and <= 1.")
    if not 0 < hsp_identity <= 100:
        raise ValueError("hsp_identity must be > 0 and <= 100.")
    if min_intron_len < 1:
        raise ValueError("min_intron_len must be at least 1.")
    if min_flank_len < 1:
        raise ValueError("min_flank_len must be at least 1.")

    source_label = source_label or ""
    input_fasta = ensure_file(input_fasta, "input FASTA")
    outdir = ensure_dir(outdir)

    blast_dir = ensure_dir(outdir / "blast")
    candidates_dir = ensure_dir(outdir / "candidates")
    confirmed_dir = ensure_dir(outdir / "confirmed")

    blastn_exe = str(blastn) if blastn else get_software("blastn", required=True)
    makeblastdb_exe = str(makeblastdb) if makeblastdb else get_software("makeblastdb", required=True)

    taxonomy_path = species_taxonomy
    if taxonomy_path is None:
        taxonomy_path = get_reference_path("species-rep1-taxonomy", required=False)

    blacklist_path = reference_blacklist
    if blacklist_path is None:
        try:
            blacklist_path = get_reference_path("reference-intron-blacklist", required=False)
        except KeyError:
            blacklist_path = None

    step("Loading input FASTA")
    query_records = fasta_to_dict(input_fasta)
    if not query_records:
        raise ValueError(f"No FASTA records found in {input_fasta}")

    step("Preparing BLAST reference")
    db_prefix = prepare_blast_reference(db, outdir, makeblastdb=makeblastdb_exe)

    local_hsps = blast_dir / "local_hsps.tsv"
    run_blastn_hsps(
        query_fasta=input_fasta,
        db_prefix=db_prefix,
        output_tsv=local_hsps,
        blastn=blastn_exe,
        threads=threads,
        perc_identity=hsp_identity,
        max_hsps=max_hsps,
        max_target_seqs=maxaccepts,
        strand=strand,
    )

    step("Parsing local HSPs")
    hsps = parse_blast_hsps(local_hsps)

    step("Detecting broken-HSP intron candidates")
    raw_candidates = detect_candidates_from_hsps(
        hsps,
        min_intron_len=min_intron_len,
        min_flank_len=min_flank_len,
        hsp_identity=hsp_identity,
        max_ref_gap=max_ref_gap,
        ref_overlap_tolerance=ref_overlap_tolerance,
    )
    clustered_candidates = cluster_candidates(raw_candidates, coordinate_tolerance=coordinate_tolerance)

    blacklist = read_reference_blacklist(blacklist_path)
    if blacklist:
        step(f"Loaded {len(blacklist)} reference-intron blacklist entries")
    elif blacklist_path:
        warning(f"Reference-intron blacklist was configured but no entries were loaded: {blacklist_path}")

    step("Loading species taxonomy and reference-intron status")
    taxonomy = read_species_taxonomy(taxonomy_path, blacklist)

    step("Summarizing candidate support")
    supports = build_cluster_support(clustered_candidates, taxonomy, blacklist)

    candidate_tsv = candidates_dir / "fracture_candidates.tsv"
    support_tsv = candidates_dir / "fracture_support_by_rank.tsv"
    failed_confirmation_tsv = candidates_dir / "failed_confirmation.tsv"

    write_tsv(
        candidate_tsv,
        CANDIDATE_FIELDS,
        [candidate_to_row(candidate, taxonomy, blacklist) for candidate in clustered_candidates],
    )
    write_tsv(support_tsv, SUPPORT_FIELDS, [support_to_row(support) for support in supports])

    if not supports:
        step("No intron-like broken-HSP candidates found")
        write_no_candidate_outputs(query_records, outdir, source_label)
        return {
            "analysis_fasta": str(confirmed_dir / "analysis_sequences.fasta"),
            "sequence_version_map": str(confirmed_dir / "sequence_version_map.tsv"),
            "candidate_introns": str(candidate_tsv),
            "support_by_rank": str(support_tsv),
            "confirmed_introns": str(confirmed_dir / "confirmed_introns.tsv"),
            "failed_confirmation": str(failed_confirmation_tsv),
            "introns_fasta": str(confirmed_dir / "introns.fasta"),
            "local_hsps": str(local_hsps),
            "confirmation_hits": str(blast_dir / "confirmation_hits.tsv"),
        }

    confirmation_fasta = blast_dir / "confirmation_candidates.fasta"
    step("Writing confirmation candidate FASTA")
    write_confirmation_fasta(supports, query_records, confirmation_fasta)

    confirmation_hits_tsv = blast_dir / "confirmation_hits.tsv"
    step("Running confirmation BLASTN")
    run_blastn_hsps(
        query_fasta=confirmation_fasta,
        db_prefix=db_prefix,
        output_tsv=confirmation_hits_tsv,
        blastn=blastn_exe,
        threads=threads,
        perc_identity=confirm_id * 100.0,
        max_hsps=5,
        max_target_seqs=50,
        strand=strand,
    )

    confirmation_hsps = parse_blast_hsps(confirmation_hits_tsv)
    confirmation_results = parse_best_confirmation_hits(
        confirmation_hsps,
        confirm_id_threshold=confirm_id,
        confirm_qcov_threshold=confirm_qcov,
    )

    confirmed, status_by_cluster, failed_rows = choose_confirmed_clusters(
        supports,
        confirmation_results,
        min_support_refs=min_support_refs,
        min_support_species=min_support_species,
        min_support_genera=min_support_genera,
    )

    for row in failed_rows:
        row["source_label"] = source_label

    write_tsv(failed_confirmation_tsv, FAILED_CONFIRMATION_FIELDS, failed_rows)

    step("Writing final analysis FASTA and intron reports")
    write_final_outputs(
        query_records,
        supports,
        confirmed,
        status_by_cluster,
        outdir,
        source_label=source_label,
    )

    success(f"Detected {len(confirmed)} confirmed intron-like insertion cluster(s)")

    return {
        "analysis_fasta": str(confirmed_dir / "analysis_sequences.fasta"),
        "sequence_version_map": str(confirmed_dir / "sequence_version_map.tsv"),
        "candidate_introns": str(candidate_tsv),
        "support_by_rank": str(support_tsv),
        "confirmed_introns": str(confirmed_dir / "confirmed_introns.tsv"),
        "failed_confirmation": str(failed_confirmation_tsv),
        "introns_fasta": str(confirmed_dir / "introns.fasta"),
        "removed_insertions_fasta": str(confirmed_dir / "removed_insertions.fasta"),
        "local_hsps": str(local_hsps),
        "confirmation_hits": str(confirmation_hits_tsv),
    }
