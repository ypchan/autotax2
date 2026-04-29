from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .config import get_software
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

COORDINATE_MAP_FIELDS = [
    "query_id",
    "cluster_id",
    "source_label",
    "intron_status",
    "confidence",
    "query_intron_start",
    "query_intron_end",
    "query_intron_length",
    "cleaned_left_pos",
    "cleaned_right_pos",
    "coordinate_ref_id",
    "coordinate_left",
    "coordinate_right",
    "coordinate_site",
    "coordinate_uncertainty",
    "mapping_identity",
    "mapping_qcov",
    "mapping_bitscore",
    "mapping_strand",
    "mapping_status",
]

POSITION_DISTRIBUTION_FIELDS = [
    "coordinate_site",
    "count",
    "high_confidence_count",
    "lineage_confidence_count",
    "shared_reference_supported_count",
    "query_specific_count",
    "median_intron_length",
    "min_intron_length",
    "max_intron_length",
    "query_ids",
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
    def qcov(self) -> float:
        if self.qlen <= 0:
            return 0.0
        return (self.qend - self.qstart + 1) / self.qlen

    @property
    def strand(self) -> str:
        return "plus" if self.sstart_raw <= self.send_raw else "minus"


@dataclass(frozen=True)
class CoordinateHit:
    query_id: str
    ref_id: str
    left_coordinate: Optional[int]
    right_coordinate: Optional[int]
    identity: float
    qcov: float
    bitscore: float
    strand: str
    status: str

    @property
    def coordinate_site(self) -> str:
        if self.left_coordinate is not None and self.right_coordinate is not None:
            return f"{self.left_coordinate}|{self.right_coordinate}"
        if self.left_coordinate is not None:
            return f"{self.left_coordinate}|NA"
        if self.right_coordinate is not None:
            return f"NA|{self.right_coordinate}"
        return ""

    @property
    def uncertainty(self) -> str:
        if self.left_coordinate is None or self.right_coordinate is None:
            return "unknown"
        return str(abs(self.right_coordinate - self.left_coordinate) - 1)


def safe_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def safe_int(value: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


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
            ".00.nhr", ".00.nin", ".00.nsq", ".00.ndb", ".00.not", ".00.ntf", ".00.nto",
        )
    ]
    return any(path.exists() for path in candidates)


def prepare_coordinate_blastdb(
    coordinate_ref: str | Path,
    outdir: str | Path,
    *,
    makeblastdb: str | Path,
) -> str:
    coordinate_ref = Path(coordinate_ref)

    if blastdb_exists(coordinate_ref):
        return str(coordinate_ref)

    coordinate_ref = ensure_file(coordinate_ref, "coordinate reference FASTA or BLAST DB prefix")
    prefix = Path(outdir) / "coordinate_ref"

    make_blast_db(
        input_fasta=coordinate_ref,
        db_prefix=prefix,
        makeblastdb=makeblastdb,
        title="autotax2_coordinate_reference",
        dbtype="nucl",
    )

    return str(prefix)


def run_coordinate_blast(
    analysis_fasta: str | Path,
    coordinate_db: str | Path,
    output_tsv: str | Path,
    *,
    blastn: str | Path,
    threads: int,
    min_identity: float = 80.0,
    max_target_seqs: int = 5,
) -> str:
    threads = validate_threads(threads)

    if not 0 <= min_identity <= 100:
        raise ValueError("min_identity must be between 0 and 100.")

    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    outfmt = (
        "6 qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore qlen slen qseq sseq"
    )

    command = [
        blastn,
        "-query",
        analysis_fasta,
        "-db",
        coordinate_db,
        "-task",
        "blastn",
        "-perc_identity",
        str(min_identity),
        "-max_hsps",
        "20",
        "-max_target_seqs",
        str(max_target_seqs),
        "-evalue",
        "1e-20",
        "-dust",
        "no",
        "-strand",
        "both",
        "-num_threads",
        str(threads),
        "-outfmt",
        outfmt,
        "-out",
        output_tsv,
    ]

    step("Mapping cleaned sequences to coordinate reference")
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


def best_hsps_by_query(hsps: Sequence[HSP]) -> Dict[str, List[HSP]]:
    by_query: Dict[str, List[HSP]] = defaultdict(list)

    for hsp in hsps:
        by_query[hsp.qseqid].append(hsp)

    for query_id, query_hsps in by_query.items():
        query_hsps.sort(key=lambda item: (item.bitscore, item.pident, item.qcov), reverse=True)

    return by_query


def map_query_position_to_subject(hsp: HSP, query_position: int) -> Optional[int]:
    """Map one cleaned-query coordinate to the coordinate reference.

    Coordinates are 1-based. If the query position falls on a gap in the
    subject row, the function returns the nearest subject coordinate observed
    at that aligned column. That keeps coordinate mapping useful near small
    indels.
    """

    if query_position < hsp.qstart or query_position > hsp.qend:
        return None

    q_step = 1 if hsp.qend_raw >= hsp.qstart_raw else -1
    s_step = 1 if hsp.send_raw >= hsp.sstart_raw else -1

    qpos = hsp.qstart_raw
    spos = hsp.sstart_raw
    last_subject_pos: Optional[int] = None

    for q_char, s_char in zip(hsp.qseq, hsp.sseq):
        current_qpos: Optional[int] = None
        current_spos: Optional[int] = None

        if q_char != "-":
            current_qpos = qpos

        if s_char != "-":
            current_spos = spos
            last_subject_pos = spos

        if current_qpos == query_position:
            if current_spos is not None:
                return current_spos
            return last_subject_pos

        if q_char != "-":
            qpos += q_step

        if s_char != "-":
            spos += s_step

    return None


def map_cleaned_site(
    hsps: Sequence[HSP],
    query_id: str,
    cleaned_left_pos: int,
    cleaned_right_pos: int,
) -> CoordinateHit:
    query_hsps = [hsp for hsp in hsps if hsp.qseqid == query_id]

    if not query_hsps:
        return CoordinateHit(query_id, "", None, None, 0.0, 0.0, 0.0, "", "no_coordinate_hit")

    # Prefer a single HSP containing both sides of the insertion site.
    for hsp in query_hsps:
        if cleaned_left_pos >= 1 and hsp.qstart <= cleaned_left_pos <= hsp.qend and hsp.qstart <= cleaned_right_pos <= hsp.qend:
            left = map_query_position_to_subject(hsp, cleaned_left_pos)
            right = map_query_position_to_subject(hsp, cleaned_right_pos)
            if left is not None or right is not None:
                return CoordinateHit(
                    query_id=query_id,
                    ref_id=hsp.sseqid,
                    left_coordinate=left,
                    right_coordinate=right,
                    identity=hsp.pident,
                    qcov=hsp.qcov,
                    bitscore=hsp.bitscore,
                    strand=hsp.strand,
                    status="mapped_both_flanks_same_hsp",
                )

    # Fallback: map left and right independently, using the best HSP that contains each position.
    left_hit: Optional[HSP] = None
    right_hit: Optional[HSP] = None

    for hsp in query_hsps:
        if cleaned_left_pos >= 1 and hsp.qstart <= cleaned_left_pos <= hsp.qend and left_hit is None:
            left_hit = hsp
        if hsp.qstart <= cleaned_right_pos <= hsp.qend and right_hit is None:
            right_hit = hsp

    left = map_query_position_to_subject(left_hit, cleaned_left_pos) if left_hit else None
    right = map_query_position_to_subject(right_hit, cleaned_right_pos) if right_hit else None

    if left_hit and right_hit and left_hit.sseqid == right_hit.sseqid:
        ref_id = left_hit.sseqid
        identity = max(left_hit.pident, right_hit.pident)
        qcov = max(left_hit.qcov, right_hit.qcov)
        bitscore = max(left_hit.bitscore, right_hit.bitscore)
        strand = left_hit.strand if left_hit.strand == right_hit.strand else "mixed"
        status = "mapped_flanks_separate_hsps"
    elif left_hit or right_hit:
        hit = left_hit or right_hit
        ref_id = hit.sseqid
        identity = hit.pident
        qcov = hit.qcov
        bitscore = hit.bitscore
        strand = hit.strand
        status = "partial_coordinate_mapping"
    else:
        ref_id = ""
        identity = 0.0
        qcov = 0.0
        bitscore = 0.0
        strand = ""
        status = "site_outside_coordinate_alignment"

    return CoordinateHit(
        query_id=query_id,
        ref_id=ref_id,
        left_coordinate=left,
        right_coordinate=right,
        identity=identity,
        qcov=qcov,
        bitscore=bitscore,
        strand=strand,
        status=status,
    )


def read_confirmed_introns(path: str | Path) -> List[Dict[str, str]]:
    path = ensure_file(path, "confirmed intron table")

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def compute_cleaned_site_positions(rows: List[Dict[str, str]]) -> List[Dict[str, object]]:
    by_query: Dict[str, List[Dict[str, str]]] = defaultdict(list)

    for row in rows:
        query_id = row.get("query_id", "")
        if query_id:
            by_query[query_id].append(row)

    enriched: List[Dict[str, object]] = []

    for query_id, query_rows in by_query.items():
        query_rows.sort(key=lambda item: safe_int(item.get("query_intron_start", "0")))
        removed_before = 0

        for row in query_rows:
            start = safe_int(row.get("query_intron_start", "0"))
            end = safe_int(row.get("query_intron_end", "0"))
            length = max(0, end - start + 1)

            cleaned_left = start - 1 - removed_before
            cleaned_right = start - removed_before

            enriched_row: Dict[str, object] = dict(row)
            enriched_row["cleaned_left_pos"] = cleaned_left
            enriched_row["cleaned_right_pos"] = cleaned_right
            enriched.append(enriched_row)

            removed_before += length

    enriched.sort(key=lambda item: (str(item.get("query_id", "")), int(item.get("query_intron_start", 0))))
    return enriched


def distribution_rows(mapped_rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    by_site: Dict[str, List[Dict[str, object]]] = defaultdict(list)

    for row in mapped_rows:
        site = str(row.get("coordinate_site", ""))
        if site:
            by_site[site].append(row)

    rows: List[Dict[str, object]] = []

    for site, site_rows in sorted(by_site.items()):
        lengths = [safe_int(str(row.get("query_intron_length", "0"))) for row in site_rows]
        confidence_counts = Counter(str(row.get("confidence", "")) for row in site_rows)
        status_counts = Counter(str(row.get("intron_status", "")) for row in site_rows)
        query_ids = sorted({str(row.get("query_id", "")) for row in site_rows if row.get("query_id")})

        rows.append(
            {
                "coordinate_site": site,
                "count": len(site_rows),
                "high_confidence_count": confidence_counts.get("high_confidence", 0),
                "lineage_confidence_count": confidence_counts.get("lineage_confidence", 0),
                "shared_reference_supported_count": status_counts.get("shared_reference_supported_intron", 0),
                "query_specific_count": status_counts.get("query_specific_intron", 0),
                "median_intron_length": int(round(median(lengths))) if lengths else "",
                "min_intron_length": min(lengths) if lengths else "",
                "max_intron_length": max(lengths) if lengths else "",
                "query_ids": ",".join(query_ids),
            }
        )

    return rows


def map_intron_coordinates(
    confirmed_introns: str | Path,
    analysis_fasta: str | Path,
    coordinate_ref: str | Path,
    outdir: str | Path,
    *,
    blastn: str | Path | None = None,
    makeblastdb: str | Path | None = None,
    threads: int = 4,
    min_identity: float = 80.0,
    max_target_seqs: int = 5,
) -> Dict[str, str]:
    """Map confirmed intron insertion sites onto a coordinate reference.

    The function uses cleaned analysis sequences, because intron-containing
    original sequences can produce broken alignments around the insertion site.
    For each confirmed intron, original query coordinates are translated into
    cleaned-sequence coordinates before mapping to the coordinate reference.
    """

    threads = validate_threads(threads)
    confirmed_introns = ensure_file(confirmed_introns, "confirmed introns TSV")
    analysis_fasta = ensure_file(analysis_fasta, "analysis FASTA")
    outdir = ensure_dir(outdir)

    if blastn is None:
        blastn = get_software("blastn", required=True)

    if makeblastdb is None:
        makeblastdb = get_software("makeblastdb", required=False)

    coordinate_db = prepare_coordinate_blastdb(
        coordinate_ref,
        outdir,
        makeblastdb=makeblastdb,
    )

    coordinate_hits = outdir / "coordinate_alignment_hsps.tsv"
    run_coordinate_blast(
        analysis_fasta=analysis_fasta,
        coordinate_db=coordinate_db,
        output_tsv=coordinate_hits,
        blastn=blastn,
        threads=threads,
        min_identity=min_identity,
        max_target_seqs=max_target_seqs,
    )

    hsps = parse_blast_hsps(coordinate_hits)
    hsp_by_query = best_hsps_by_query(hsps)

    confirmed_rows = read_confirmed_introns(confirmed_introns)
    enriched_rows = compute_cleaned_site_positions(confirmed_rows)

    mapped_rows: List[Dict[str, object]] = []

    for row in enriched_rows:
        query_id = str(row.get("query_id", ""))
        cleaned_left = safe_int(str(row.get("cleaned_left_pos", "0")))
        cleaned_right = safe_int(str(row.get("cleaned_right_pos", "0")))

        hit = map_cleaned_site(
            hsp_by_query.get(query_id, []),
            query_id,
            cleaned_left,
            cleaned_right,
        )

        mapped_rows.append(
            {
                "query_id": query_id,
                "cluster_id": row.get("cluster_id", ""),
                "source_label": row.get("source_label", ""),
                "intron_status": row.get("intron_status", ""),
                "confidence": row.get("confidence", ""),
                "query_intron_start": row.get("query_intron_start", ""),
                "query_intron_end": row.get("query_intron_end", ""),
                "query_intron_length": row.get("query_intron_length", ""),
                "cleaned_left_pos": cleaned_left,
                "cleaned_right_pos": cleaned_right,
                "coordinate_ref_id": hit.ref_id,
                "coordinate_left": hit.left_coordinate if hit.left_coordinate is not None else "",
                "coordinate_right": hit.right_coordinate if hit.right_coordinate is not None else "",
                "coordinate_site": hit.coordinate_site,
                "coordinate_uncertainty": hit.uncertainty,
                "mapping_identity": f"{hit.identity:.3f}",
                "mapping_qcov": f"{hit.qcov:.6f}",
                "mapping_bitscore": f"{hit.bitscore:.3f}",
                "mapping_strand": hit.strand,
                "mapping_status": hit.status,
            }
        )

    coordinate_map = outdir / "intron_coordinate_map.tsv"
    distribution = outdir / "intron_position_distribution.tsv"

    write_tsv(coordinate_map, COORDINATE_MAP_FIELDS, mapped_rows)
    write_tsv(distribution, POSITION_DISTRIBUTION_FIELDS, distribution_rows(mapped_rows))

    success(f"Mapped {len(mapped_rows)} intron insertion site(s) to coordinate reference")

    return {
        "coordinate_alignment_hsps": str(coordinate_hits),
        "intron_coordinate_map": str(coordinate_map),
        "intron_position_distribution": str(distribution),
    }
