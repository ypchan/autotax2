from __future__ import annotations

import csv
import heapq
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from .external import barrnap_annotate, blastn_local_hsps, make_blast_db
from .logging import ProgressCounter, step, success, warning
from .threads import validate_threads
from .utils import ensure_dir, ensure_file, read_fasta, write_manifest


TAXONOMY_COLUMNS = [
    "Domain",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
]


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


REFERENCE_INTRON_FIELDS = [
    "candidate_id",
    "query_id",
    "subject_id",
    "query_species",
    "subject_species",
    "query_intron_start",
    "query_intron_end",
    "query_intron_length",
    "ref_gap",
    "strand",
    "left_qstart",
    "left_qend",
    "right_qstart",
    "right_qend",
    "left_sstart",
    "left_send",
    "right_sstart",
    "right_send",
    "left_identity",
    "right_identity",
    "mean_identity",
    "left_length",
    "right_length",
    "left_bitscore",
    "right_bitscore",
]


REFERENCE_AUDIT_FIELDS = [
    "seq_id",
    "species_key",
    "n_candidates",
    "supporting_refs",
    "supporting_species",
    "median_intron_start",
    "median_intron_end",
    "median_intron_length",
    "max_mean_identity",
    "ref_has_putative_intron",
]


BARRNAP_AUDIT_FIELDS = [
    "seq_id",
    "n_16s_hits",
    "best_16s_start",
    "best_16s_end",
    "best_16s_length",
    "has_partial",
    "has_multiple_16s",
    "barrnap_status",
]


CLEAN_SELECTION_FIELDS = [
    "rep_id",
    "original_id",
    "species_key",
    "sequence_length",
    "selected_for_clean_rep1",
    "selection_reason",
    "ref_has_putative_intron",
    "barrnap_status",
    "Domain",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
]


@dataclass(frozen=True)
class TaxonomyRecord:
    seq_id: str
    domain: str = ""
    phylum: str = ""
    class_name: str = ""
    order: str = ""
    family: str = ""
    genus: str = ""
    species: str = ""

    def as_rank_dict(self) -> Dict[str, str]:
        return {
            "Domain": self.domain,
            "Phylum": self.phylum,
            "Class": self.class_name,
            "Order": self.order,
            "Family": self.family,
            "Genus": self.genus,
            "Species": self.species,
        }


@dataclass(frozen=True)
class RepresentativeRecord:
    species_key: str
    seq_id: str
    header: str
    sequence: str
    length: int
    taxonomy: TaxonomyRecord


@dataclass(frozen=True)
class BlastHSP:
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


@dataclass(frozen=True)
class ReferenceIntronCandidate:
    candidate_id: str
    query_id: str
    subject_id: str
    query_species: str
    subject_species: str
    query_intron_start: int
    query_intron_end: int
    query_intron_length: int
    ref_gap: int
    strand: str
    left_hsp: BlastHSP
    right_hsp: BlastHSP

    @property
    def mean_identity(self) -> float:
        return (self.left_hsp.pident + self.right_hsp.pident) / 2.0

    @property
    def mean_bitscore(self) -> float:
        return (self.left_hsp.bitscore + self.right_hsp.bitscore) / 2.0


def normalize_taxon(value: Optional[str]) -> str:
    if value is None:
        return ""

    value = value.strip()

    if not value:
        return ""

    return value.replace("\t", " ").replace(";", ".")


def fasta_seq_id(header: str) -> str:
    return header.split()[0].split(";")[0]


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


def species_key_from_taxonomy(record: TaxonomyRecord) -> Optional[str]:
    species = normalize_taxon(record.species)

    if not species:
        return None

    genus = normalize_taxon(record.genus)

    if genus and not species.startswith(genus):
        return f"{genus} {species}"

    return species


def read_silva_taxonomy(path: str | Path) -> Dict[str, TaxonomyRecord]:
    path = ensure_file(path, "SILVA taxonomy table")
    taxonomy: Dict[str, TaxonomyRecord] = {}

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"Empty taxonomy table: {path}")

        if "ID" not in reader.fieldnames:
            raise ValueError("SILVA taxonomy table must contain an ID column.")

        for row in reader:
            seq_id = normalize_taxon(row.get("ID"))

            if not seq_id:
                continue

            taxonomy[seq_id] = TaxonomyRecord(
                seq_id=seq_id,
                domain=normalize_taxon(row.get("Domain")),
                phylum=normalize_taxon(row.get("Phylum")),
                class_name=normalize_taxon(row.get("Class")),
                order=normalize_taxon(row.get("Order")),
                family=normalize_taxon(row.get("Family")),
                genus=normalize_taxon(row.get("Genus")),
                species=normalize_taxon(row.get("Species")),
            )

    return taxonomy


def _heap_item(record: RepresentativeRecord) -> Tuple[int, str, RepresentativeRecord]:
    return (record.length, record.seq_id, record)


def select_species_representatives(
    silva_fasta: str | Path,
    taxonomy: Dict[str, TaxonomyRecord],
    *,
    per_species: int,
) -> List[RepresentativeRecord]:
    if per_species < 1:
        raise ValueError("per_species must be at least 1.")

    silva_fasta = ensure_file(silva_fasta, "SILVA FASTA")

    selected: Dict[str, List[Tuple[int, str, RepresentativeRecord]]] = {}
    progress = ProgressCounter("SILVA FASTA records", interval=100000)

    for header, sequence in read_fasta(silva_fasta):
        progress.update()

        seq_id = fasta_seq_id(header)
        tax = taxonomy.get(seq_id)

        if tax is None:
            continue

        species_key = species_key_from_taxonomy(tax)

        if species_key is None:
            continue

        seq = "".join(str(sequence).split()).upper()

        if not seq:
            continue

        record = RepresentativeRecord(
            species_key=species_key,
            seq_id=seq_id,
            header=header,
            sequence=seq,
            length=len(seq),
            taxonomy=tax,
        )

        heap = selected.setdefault(species_key, [])
        item = _heap_item(record)

        if len(heap) < per_species:
            heapq.heappush(heap, item)
        elif item > heap[0]:
            heapq.heapreplace(heap, item)

    progress.finish()

    representatives: List[RepresentativeRecord] = []

    for heap in selected.values():
        for _length, _seq_id, record in heap:
            representatives.append(record)

    representatives.sort(
        key=lambda item: (
            item.species_key,
            -item.length,
            item.seq_id,
        )
    )

    return representatives


def write_fasta_records(
    records: Iterable[RepresentativeRecord],
    output_fasta: str | Path,
) -> None:
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    with output_fasta.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record.header}\n")

            seq = record.sequence
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")


def write_taxonomy_table(
    records: Iterable[RepresentativeRecord],
    output_tsv: str | Path,
    *,
    selected_ids: Optional[Set[str]] = None,
    selection_reason: Optional[Dict[str, str]] = None,
    reference_blacklist: Optional[Set[str]] = None,
    barrnap_status: Optional[Dict[str, str]] = None,
) -> None:
    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    selected_ids = selected_ids or set()
    selection_reason = selection_reason or {}
    reference_blacklist = reference_blacklist or set()
    barrnap_status = barrnap_status or {}

    with output_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CLEAN_SELECTION_FIELDS, delimiter="\t")
        writer.writeheader()

        for record in records:
            tax = record.taxonomy.as_rank_dict()

            writer.writerow(
                {
                    "rep_id": record.seq_id,
                    "original_id": record.seq_id,
                    "species_key": record.species_key,
                    "sequence_length": str(record.length),
                    "selected_for_clean_rep1": "yes" if record.seq_id in selected_ids else "no",
                    "selection_reason": selection_reason.get(record.seq_id, ""),
                    "ref_has_putative_intron": "yes" if record.seq_id in reference_blacklist else "no",
                    "barrnap_status": barrnap_status.get(record.seq_id, ""),
                    "Domain": tax["Domain"],
                    "Phylum": tax["Phylum"],
                    "Class": tax["Class"],
                    "Order": tax["Order"],
                    "Family": tax["Family"],
                    "Genus": tax["Genus"],
                    "Species": tax["Species"],
                }
            )


def parse_blast_hsps(path: str | Path) -> List[BlastHSP]:
    hsps: List[BlastHSP] = []
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
                BlastHSP(
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


def reference_gap(left: BlastHSP, right: BlastHSP) -> int:
    if left.strand == "plus":
        return right.smin - left.smax - 1

    return left.smin - right.smax - 1


def detect_reference_intron_candidates(
    hsps: List[BlastHSP],
    taxonomy_by_id: Dict[str, RepresentativeRecord],
    *,
    hsp_identity: float = 90.0,
    min_flank_len: int = 100,
    min_intron_len: int = 20,
    max_ref_gap: int = 10,
    ref_overlap_tolerance: int = 5,
) -> List[ReferenceIntronCandidate]:
    grouped: Dict[Tuple[str, str, str], List[BlastHSP]] = defaultdict(list)

    for hsp in hsps:
        if hsp.qseqid == hsp.sseqid:
            continue

        if hsp.pident < hsp_identity:
            continue

        if hsp.length < min_flank_len:
            continue

        grouped[(hsp.qseqid, hsp.sseqid, hsp.strand)].append(hsp)

    candidates: List[ReferenceIntronCandidate] = []
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

                query_record = taxonomy_by_id.get(left.qseqid)
                subject_record = taxonomy_by_id.get(left.sseqid)

                counter += 1

                candidates.append(
                    ReferenceIntronCandidate(
                        candidate_id=f"ref_candidate_{counter:06d}",
                        query_id=left.qseqid,
                        subject_id=left.sseqid,
                        query_species=query_record.species_key if query_record else "",
                        subject_species=subject_record.species_key if subject_record else "",
                        query_intron_start=left.qend + 1,
                        query_intron_end=right.qstart - 1,
                        query_intron_length=query_gap,
                        ref_gap=ref_gap,
                        strand=left.strand,
                        left_hsp=left,
                        right_hsp=right,
                    )
                )

    return candidates


def candidate_to_row(candidate: ReferenceIntronCandidate) -> Dict[str, object]:
    return {
        "candidate_id": candidate.candidate_id,
        "query_id": candidate.query_id,
        "subject_id": candidate.subject_id,
        "query_species": candidate.query_species,
        "subject_species": candidate.subject_species,
        "query_intron_start": candidate.query_intron_start,
        "query_intron_end": candidate.query_intron_end,
        "query_intron_length": candidate.query_intron_length,
        "ref_gap": candidate.ref_gap,
        "strand": candidate.strand,
        "left_qstart": candidate.left_hsp.qstart,
        "left_qend": candidate.left_hsp.qend,
        "right_qstart": candidate.right_hsp.qstart,
        "right_qend": candidate.right_hsp.qend,
        "left_sstart": candidate.left_hsp.sstart_raw,
        "left_send": candidate.left_hsp.send_raw,
        "right_sstart": candidate.right_hsp.sstart_raw,
        "right_send": candidate.right_hsp.send_raw,
        "left_identity": f"{candidate.left_hsp.pident:.3f}",
        "right_identity": f"{candidate.right_hsp.pident:.3f}",
        "mean_identity": f"{candidate.mean_identity:.3f}",
        "left_length": candidate.left_hsp.length,
        "right_length": candidate.right_hsp.length,
        "left_bitscore": f"{candidate.left_hsp.bitscore:.3f}",
        "right_bitscore": f"{candidate.right_hsp.bitscore:.3f}",
    }


def write_tsv(
    path: str | Path,
    fields: Sequence[str],
    rows: Iterable[Dict[str, object]],
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t")
        writer.writeheader()

        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def summarize_reference_audit(
    candidates: List[ReferenceIntronCandidate],
    *,
    min_support_refs: int = 2,
) -> Tuple[List[Dict[str, object]], Set[str]]:
    by_query: Dict[str, List[ReferenceIntronCandidate]] = defaultdict(list)

    for candidate in candidates:
        by_query[candidate.query_id].append(candidate)

    rows: List[Dict[str, object]] = []
    blacklist: Set[str] = set()

    for query_id, group in sorted(by_query.items()):
        supporting_refs = {item.subject_id for item in group}
        supporting_species = {item.subject_species for item in group if item.subject_species}
        starts = [item.query_intron_start for item in group]
        ends = [item.query_intron_end for item in group]
        lengths = [item.query_intron_length for item in group]
        max_identity = max(item.mean_identity for item in group) if group else 0.0

        is_flagged = len(supporting_refs) >= min_support_refs

        if is_flagged:
            blacklist.add(query_id)

        rows.append(
            {
                "seq_id": query_id,
                "species_key": group[0].query_species if group else "",
                "n_candidates": len(group),
                "supporting_refs": len(supporting_refs),
                "supporting_species": len(supporting_species),
                "median_intron_start": int(round(median(starts))) if starts else "",
                "median_intron_end": int(round(median(ends))) if ends else "",
                "median_intron_length": int(round(median(lengths))) if lengths else "",
                "max_mean_identity": f"{max_identity:.3f}",
                "ref_has_putative_intron": "yes" if is_flagged else "no",
            }
        )

    return rows, blacklist


def parse_barrnap_gff(path: str | Path) -> Dict[str, Dict[str, object]]:
    path = Path(path)
    per_seq: Dict[str, List[Dict[str, object]]] = defaultdict(list)

    if not path.exists() or path.stat().st_size == 0:
        return {}

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 9:
                continue

            seq_id, _source, feature_type, start, end, _score, _strand, _phase, attrs = parts

            text = f"{feature_type} {attrs}".lower()

            if "16s" not in text and "ssu" not in text:
                continue

            start_i = safe_int(start)
            end_i = safe_int(end)

            per_seq[seq_id].append(
                {
                    "start": min(start_i, end_i),
                    "end": max(start_i, end_i),
                    "length": abs(end_i - start_i) + 1,
                    "partial": "partial" in text,
                }
            )

    summary: Dict[str, Dict[str, object]] = {}

    for seq_id, hits in per_seq.items():
        best = max(hits, key=lambda item: int(item["length"]))
        has_partial = any(bool(item["partial"]) for item in hits)
        has_multiple = len(hits) > 1

        if has_multiple:
            status = "multiple_16s"
        elif has_partial:
            status = "partial_16s"
        else:
            status = "single_complete_16s"

        summary[seq_id] = {
            "seq_id": seq_id,
            "n_16s_hits": len(hits),
            "best_16s_start": best["start"],
            "best_16s_end": best["end"],
            "best_16s_length": best["length"],
            "has_partial": "yes" if has_partial else "no",
            "has_multiple_16s": "yes" if has_multiple else "no",
            "barrnap_status": status,
        }

    return summary


def write_barrnap_audit_table(
    raw_records: List[RepresentativeRecord],
    barrnap_summary: Dict[str, Dict[str, object]],
    output_tsv: str | Path,
) -> Dict[str, str]:
    rows: List[Dict[str, object]] = []
    status_by_id: Dict[str, str] = {}

    for record in raw_records:
        row = barrnap_summary.get(record.seq_id)

        if row is None:
            row = {
                "seq_id": record.seq_id,
                "n_16s_hits": 0,
                "best_16s_start": "",
                "best_16s_end": "",
                "best_16s_length": "",
                "has_partial": "",
                "has_multiple_16s": "",
                "barrnap_status": "no_16s_hit",
            }

        status_by_id[record.seq_id] = str(row["barrnap_status"])
        rows.append(row)

    write_tsv(output_tsv, BARRNAP_AUDIT_FIELDS, rows)
    return status_by_id


def choose_clean_rep1(
    raw_records: List[RepresentativeRecord],
    reference_blacklist: Set[str],
    barrnap_status: Dict[str, str],
) -> Tuple[List[RepresentativeRecord], Dict[str, str]]:
    by_species: Dict[str, List[RepresentativeRecord]] = defaultdict(list)

    for record in raw_records:
        by_species[record.species_key].append(record)

    clean_records: List[RepresentativeRecord] = []
    reason_by_id: Dict[str, str] = {}

    preferred_barrnap_status = {"single_complete_16s", ""}

    for species_key, records in sorted(by_species.items()):
        records = sorted(records, key=lambda item: (-item.length, item.seq_id))

        clean_candidates = [
            record
            for record in records
            if record.seq_id not in reference_blacklist
            and barrnap_status.get(record.seq_id, "") in preferred_barrnap_status
        ]

        if clean_candidates:
            chosen = clean_candidates[0]
            reason = "not_flagged_by_reference_audit_and_barrnap_ok"

        else:
            non_intron_candidates = [
                record
                for record in records
                if record.seq_id not in reference_blacklist
            ]

            if non_intron_candidates:
                chosen = non_intron_candidates[0]
                reason = "not_flagged_by_reference_audit_barrnap_not_ok_or_missing"

            else:
                chosen = records[0]
                reason = "fallback_all_candidates_flagged_or_missing"

        clean_records.append(chosen)
        reason_by_id[chosen.seq_id] = reason

    clean_records.sort(key=lambda item: (-item.length, item.species_key, item.seq_id))
    return clean_records, reason_by_id


def build_raw_species_reference(
    silva_fasta: str | Path,
    silva_taxonomy: str | Path,
    outdir: str | Path,
    *,
    per_species: int,
    makeblastdb: str | Path = "makeblastdb",
    build_blastdb: bool = True,
) -> Tuple[List[RepresentativeRecord], Dict[str, str]]:
    outdir = ensure_dir(outdir)

    label = f"species_rep{per_species}_raw"
    output_fasta = outdir / f"{label}.fasta"
    output_taxonomy = outdir / f"{label}.taxonomy.tsv"
    blastdb_prefix = outdir / label

    step("Reading SILVA taxonomy")
    taxonomy = read_silva_taxonomy(silva_taxonomy)

    step(f"Selecting longest {per_species} raw representative sequence(s) per species")
    raw_records = select_species_representatives(
        silva_fasta,
        taxonomy,
        per_species=per_species,
    )

    if not raw_records:
        raise ValueError(
            "No species representatives were selected. "
            "Check that silva_taxonomy.tsv contains a non-empty Species column."
        )

    step(f"Writing {label} FASTA")
    write_fasta_records(raw_records, output_fasta)

    step(f"Writing {label} taxonomy table")
    write_taxonomy_table(raw_records, output_taxonomy)

    entries: Dict[str, str] = {
        f"{label}_fasta": str(output_fasta),
        f"{label}_taxonomy": str(output_taxonomy),
        f"{label}_n_sequences": str(len(raw_records)),
    }

    if per_species == 10:
        entries["species_rep10_fasta"] = str(output_fasta)
        entries["species_rep10_taxonomy"] = str(output_taxonomy)

    if build_blastdb:
        make_blast_db(
            output_fasta,
            blastdb_prefix,
            makeblastdb=makeblastdb,
            title=label,
            dbtype="nucl",
        )
        entries[f"{label}_blastdb"] = str(blastdb_prefix)

        if per_species == 10:
            entries["species_rep10_blastdb"] = str(blastdb_prefix)

    return raw_records, entries


def audit_reference_introns(
    raw_fasta: str | Path,
    raw_blastdb: str | Path,
    raw_records: List[RepresentativeRecord],
    outdir: str | Path,
    *,
    blastn: str | Path = "blastn",
    threads: int = 4,
    hsp_identity: float = 90.0,
    min_flank_len: int = 100,
    min_intron_len: int = 20,
    max_ref_gap: int = 10,
    ref_overlap_tolerance: int = 5,
    max_hsps: int = 20,
    max_target_seqs: int = 500,
    min_support_refs: int = 2,
) -> Tuple[Set[str], Dict[str, str]]:
    outdir = ensure_dir(outdir)

    audit_hsps = outdir / "reference_local_hsps.tsv"
    candidates_tsv = outdir / "reference_intron_candidates.tsv"
    audit_tsv = outdir / "reference_intron_audit.tsv"
    blacklist_txt = outdir / "reference_intron_blacklist.txt"

    blastn_local_hsps(
        query_fasta=raw_fasta,
        db=raw_blastdb,
        output_tsv=audit_hsps,
        blastn=blastn,
        threads=threads,
        perc_identity=hsp_identity,
        max_hsps=max_hsps,
        max_target_seqs=max_target_seqs,
    )

    taxonomy_by_id = {record.seq_id: record for record in raw_records}

    step("Parsing reference local HSPs")
    hsps = parse_blast_hsps(audit_hsps)

    step("Detecting putative introns in reference sequences")
    candidates = detect_reference_intron_candidates(
        hsps,
        taxonomy_by_id,
        hsp_identity=hsp_identity,
        min_flank_len=min_flank_len,
        min_intron_len=min_intron_len,
        max_ref_gap=max_ref_gap,
        ref_overlap_tolerance=ref_overlap_tolerance,
    )

    write_tsv(
        candidates_tsv,
        REFERENCE_INTRON_FIELDS,
        [candidate_to_row(candidate) for candidate in candidates],
    )

    audit_rows, blacklist = summarize_reference_audit(
        candidates,
        min_support_refs=min_support_refs,
    )

    write_tsv(audit_tsv, REFERENCE_AUDIT_FIELDS, audit_rows)

    with blacklist_txt.open("w", encoding="utf-8") as handle:
        for seq_id in sorted(blacklist):
            handle.write(seq_id + "\n")

    status = {
        row["seq_id"]: row["ref_has_putative_intron"]
        for row in audit_rows
    }

    return blacklist, {
        "reference_local_hsps": str(audit_hsps),
        "reference_intron_candidates": str(candidates_tsv),
        "reference_intron_audit": str(audit_tsv),
        "reference_intron_blacklist": str(blacklist_txt),
        "reference_intron_blacklist_count": str(len(blacklist)),
        **{f"reference_intron_status:{key}": value for key, value in status.items()},
    }


def run_barrnap_reference_audit(
    raw_fasta: str | Path,
    raw_records: List[RepresentativeRecord],
    outdir: str | Path,
    *,
    barrnap: str | Path,
    threads: int,
) -> Tuple[Dict[str, str], Dict[str, str]]:
    outdir = ensure_dir(outdir)

    gff = outdir / "reference_barrnap.gff"
    audit_tsv = outdir / "reference_barrnap_audit.tsv"

    barrnap_annotate(
        input_fasta=raw_fasta,
        output_gff=gff,
        barrnap=barrnap,
        threads=threads,
    )

    summary = parse_barrnap_gff(gff)
    status_by_id = write_barrnap_audit_table(raw_records, summary, audit_tsv)

    return status_by_id, {
        "reference_barrnap_gff": str(gff),
        "reference_barrnap_audit": str(audit_tsv),
    }


def write_clean_rep1(
    clean_records: List[RepresentativeRecord],
    raw_records: List[RepresentativeRecord],
    outdir: str | Path,
    *,
    makeblastdb: str | Path,
    build_blastdb: bool,
    reference_blacklist: Set[str],
    barrnap_status: Dict[str, str],
    selection_reason: Dict[str, str],
) -> Dict[str, str]:
    outdir = ensure_dir(outdir)

    output_fasta = outdir / "species_rep1_clean.fasta"
    output_taxonomy = outdir / "species_rep1_clean.taxonomy.tsv"
    all_selection_tsv = outdir / "species_rep10_raw.selection.tsv"
    blastdb_prefix = outdir / "species_rep1_clean"

    selected_ids = {record.seq_id for record in clean_records}

    step("Writing species_rep1_clean FASTA")
    write_fasta_records(clean_records, output_fasta)

    step("Writing species_rep1_clean taxonomy")
    write_taxonomy_table(
        clean_records,
        output_taxonomy,
        selected_ids=selected_ids,
        selection_reason=selection_reason,
        reference_blacklist=reference_blacklist,
        barrnap_status=barrnap_status,
    )

    step("Writing raw-reference selection audit table")
    write_taxonomy_table(
        raw_records,
        all_selection_tsv,
        selected_ids=selected_ids,
        selection_reason=selection_reason,
        reference_blacklist=reference_blacklist,
        barrnap_status=barrnap_status,
    )

    entries = {
        "species_rep1_fasta": str(output_fasta),
        "species_rep1_taxonomy": str(output_taxonomy),
        "species_rep1_clean_fasta": str(output_fasta),
        "species_rep1_clean_taxonomy": str(output_taxonomy),
        "species_rep10_raw_selection": str(all_selection_tsv),
        "species_rep1_clean_n_sequences": str(len(clean_records)),
    }

    if build_blastdb:
        make_blast_db(
            output_fasta,
            blastdb_prefix,
            makeblastdb=makeblastdb,
            title="species_rep1_clean",
            dbtype="nucl",
        )
        entries["species_rep1_blastdb"] = str(blastdb_prefix)
        entries["species_rep1_clean_blastdb"] = str(blastdb_prefix)

    return entries


def build_species_representative_reference(
    silva_fasta: str | Path,
    silva_taxonomy: str | Path,
    outdir: str | Path,
    *,
    per_species: int = 1,
    makeblastdb: str | Path = "makeblastdb",
    threads: int = 4,
    build_blastdb: bool = True,
) -> Dict[str, str]:
    """Build a raw species representative reference.

    This function is useful for quick testing or custom references.
    For production intron detection, prefer build_default_intron_references(),
    which builds species_rep10_raw, audits reference introns, and produces
    species_rep1_clean.
    """

    validate_threads(threads)

    raw_records, entries = build_raw_species_reference(
        silva_fasta=silva_fasta,
        silva_taxonomy=silva_taxonomy,
        outdir=outdir,
        per_species=per_species,
        makeblastdb=makeblastdb,
        build_blastdb=build_blastdb,
    )

    label = f"species_rep{per_species}"
    raw_label = f"species_rep{per_species}_raw"

    normalized_entries = {
        f"{label}_fasta": entries[f"{raw_label}_fasta"],
        f"{label}_taxonomy": entries[f"{raw_label}_taxonomy"],
        f"{label}_n_sequences": str(len(raw_records)),
    }

    if f"{raw_label}_blastdb" in entries:
        normalized_entries[f"{label}_blastdb"] = entries[f"{raw_label}_blastdb"]

    normalized_entries.update(entries)

    manifest = Path(outdir) / f"{label}.manifest.tsv"
    write_manifest(manifest, normalized_entries)
    normalized_entries["manifest"] = str(manifest)

    success(f"Built {label}: {len(raw_records)} representative sequences")
    return normalized_entries


def build_default_intron_references(
    silva_fasta: str | Path,
    silva_taxonomy: str | Path,
    outdir: str | Path,
    *,
    makeblastdb: str | Path = "makeblastdb",
    blastn: str | Path = "blastn",
    barrnap: str | Path | None = None,
    threads: int = 4,
    build_blastdb: bool = True,
    audit_reference_introns: bool = True,
    run_barrnap: bool = False,
    raw_per_species: int = 10,
    hsp_identity: float = 90.0,
    min_flank_len: int = 100,
    min_intron_len: int = 20,
    max_ref_gap: int = 10,
    ref_overlap_tolerance: int = 5,
    max_hsps: int = 20,
    max_target_seqs: int = 500,
    min_reference_support_refs: int = 2,
) -> Dict[str, str]:
    """Build the production intron-detection reference.

    Output design:
      species_rep10_raw
        Longest 10 raw references per species. Used for audit and sensitivity.

      reference_intron_audit
        Local broken-HSP audit that flags raw references with putative introns.

      species_rep1_clean
        One clean representative per species, avoiding flagged references when
        possible. This is the default reference for detect-intron.
    """

    threads = validate_threads(threads)
    outdir = ensure_dir(outdir)

    raw_records, entries = build_raw_species_reference(
        silva_fasta=silva_fasta,
        silva_taxonomy=silva_taxonomy,
        outdir=outdir,
        per_species=raw_per_species,
        makeblastdb=makeblastdb,
        build_blastdb=build_blastdb,
    )

    raw_label = f"species_rep{raw_per_species}_raw"
    raw_fasta = entries[f"{raw_label}_fasta"]
    raw_blastdb = entries.get(f"{raw_label}_blastdb")

    reference_blacklist: Set[str] = set()
    reference_audit_entries: Dict[str, str] = {}

    if audit_reference_introns:
        if not build_blastdb or raw_blastdb is None:
            warning("Skipping reference intron audit because BLAST DB was not built.")
        else:
            blacklist, reference_audit_entries = audit_reference_introns_fn(
                raw_fasta=raw_fasta,
                raw_blastdb=raw_blastdb,
                raw_records=raw_records,
                outdir=outdir,
                blastn=blastn,
                threads=threads,
                hsp_identity=hsp_identity,
                min_flank_len=min_flank_len,
                min_intron_len=min_intron_len,
                max_ref_gap=max_ref_gap,
                ref_overlap_tolerance=ref_overlap_tolerance,
                max_hsps=max_hsps,
                max_target_seqs=max_target_seqs,
                min_support_refs=min_reference_support_refs,
            )
            reference_blacklist = blacklist

    barrnap_status: Dict[str, str] = {}
    barrnap_entries: Dict[str, str] = {}

    if run_barrnap:
        if barrnap is None:
            warning("Skipping barrnap audit because no barrnap executable was provided.")
        else:
            barrnap_status, barrnap_entries = run_barrnap_reference_audit(
                raw_fasta=raw_fasta,
                raw_records=raw_records,
                outdir=outdir,
                barrnap=barrnap,
                threads=threads,
            )

    clean_records, selection_reason = choose_clean_rep1(
        raw_records=raw_records,
        reference_blacklist=reference_blacklist,
        barrnap_status=barrnap_status,
    )

    clean_entries = write_clean_rep1(
        clean_records=clean_records,
        raw_records=raw_records,
        outdir=outdir,
        makeblastdb=makeblastdb,
        build_blastdb=build_blastdb,
        reference_blacklist=reference_blacklist,
        barrnap_status=barrnap_status,
        selection_reason=selection_reason,
    )

    entries.update(reference_audit_entries)
    entries.update(barrnap_entries)
    entries.update(clean_entries)

    entries["intron_reference_raw_per_species"] = str(raw_per_species)
    entries["intron_reference_audit_enabled"] = str(audit_reference_introns)
    entries["intron_reference_barrnap_enabled"] = str(run_barrnap)

    manifest = outdir / "intron_reference_manifest.tsv"
    write_manifest(manifest, entries)
    entries["manifest"] = str(manifest)

    success(
        f"Built intron reference: raw={len(raw_records)} records, "
        f"clean={len(clean_records)} records, "
        f"flagged_reference_introns={len(reference_blacklist)}"
    )

    return entries


def audit_reference_introns_fn(
    raw_fasta: str | Path,
    raw_blastdb: str | Path,
    raw_records: List[RepresentativeRecord],
    outdir: str | Path,
    *,
    blastn: str | Path,
    threads: int,
    hsp_identity: float,
    min_flank_len: int,
    min_intron_len: int,
    max_ref_gap: int,
    ref_overlap_tolerance: int,
    max_hsps: int,
    max_target_seqs: int,
    min_support_refs: int,
) -> Tuple[Set[str], Dict[str, str]]:
    return audit_reference_introns(
        raw_fasta=raw_fasta,
        raw_blastdb=raw_blastdb,
        raw_records=raw_records,
        outdir=outdir,
        blastn=blastn,
        threads=threads,
        hsp_identity=hsp_identity,
        min_flank_len=min_flank_len,
        min_intron_len=min_intron_len,
        max_ref_gap=max_ref_gap,
        ref_overlap_tolerance=ref_overlap_tolerance,
        max_hsps=max_hsps,
        max_target_seqs=max_target_seqs,
        min_support_refs=min_support_refs,
    )