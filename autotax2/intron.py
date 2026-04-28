from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List, Tuple

from .logging import success
from .utils import ensure_dir, ensure_file, read_fasta, run_cmd, write_fasta
from .vsearch import vsearch_path


def fasta_records_by_id(path: str | Path) -> Dict[str, Tuple[str, str]]:
    return {h.split()[0].split(";")[0]: (h, s) for h, s in read_fasta(path)}


def run_vsearch_alignment(
    input_fasta: str | Path,
    db: str | Path,
    outdir: str | Path,
    vsearch: str = "vsearch",
    threads: int = 8,
    search_id: float = 0.70,
    strand: str = "both",
    maxaccepts: int = 1,
    dry_run: bool = False,
) -> Dict[str, str]:
    outdir = ensure_dir(outdir)
    hits = outdir / "vsearch_alignment_rows.tsv"
    uc = outdir / "vsearch_alignment.uc"
    cmd = [
        vsearch_path(vsearch),
        "--usearch_global", str(input_fasta),
        "--db", str(db),
        "--id", str(search_id),
        "--strand", strand,
        "--maxaccepts", str(maxaccepts),
        "--maxrejects", "0",
        "--userout", str(hits),
        "--userfields", "query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+qstrand+qrow+trow",
        "--uc", str(uc),
        "--threads", str(threads),
    ]
    run_cmd(cmd, dry_run=dry_run)
    return {"hits": str(hits), "uc": str(uc)}


def read_alignment_rows(path: str | Path) -> Dict[str, Dict[str, str]]:
    fields = [
        "query", "target", "id", "alnlen", "mism", "opens",
        "qlo", "qhi", "tlo", "thi", "qstrand", "qrow", "trow",
    ]
    out: Dict[str, Dict[str, str]] = {}
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return out
    with open(p, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip():
                continue
            row = dict(zip(fields, line.rstrip("\n").split("\t")))
            q = row.get("query", "")
            if q and q not in out:
                out[q] = row
    return out


def find_query_insertions(qrow: str, trow: str, min_intron_len: int = 50) -> List[Dict[str, int]]:
    if len(qrow) != len(trow):
        raise ValueError("qrow and trow have different lengths")
    insertions: List[Dict[str, int]] = []
    qpos = 0
    i = 0
    n = len(qrow)
    while i < n:
        if qrow[i] != "-":
            qpos += 1
        if qrow[i] != "-" and trow[i] == "-":
            start_q = qpos
            length = 1
            j = i + 1
            qpos_j = qpos
            while j < n:
                if qrow[j] != "-":
                    qpos_j += 1
                if qrow[j] != "-" and trow[j] == "-":
                    length += 1
                    j += 1
                    continue
                break
            end_q = qpos_j
            if length >= min_intron_len:
                insertions.append({
                    "query_start": start_q,
                    "query_end": end_q,
                    "length": length,
                    "alignment_start": i + 1,
                    "alignment_end": j,
                })
            qpos = qpos_j
            i = j
            continue
        i += 1
    return insertions


def identity_excluding_insertions(qrow: str, trow: str, insertions: List[Dict[str, int]]) -> float:
    masked = set()
    for ins in insertions:
        masked.update(range(ins["alignment_start"] - 1, ins["alignment_end"]))
    matches = 0
    comparable = 0
    for i, (q, t) in enumerate(zip(qrow, trow)):
        if i in masked or q == "-" or t == "-":
            continue
        comparable += 1
        if q.upper() == t.upper():
            matches += 1
    return matches / comparable if comparable else 0.0


def flank_lengths(insertions: List[Dict[str, int]], seq_len: int) -> Tuple[int, int]:
    if not insertions:
        return seq_len, seq_len
    return min(x["query_start"] for x in insertions) - 1, seq_len - max(x["query_end"] for x in insertions)


def remove_intervals(seq: str, intervals: List[Tuple[int, int]]) -> str:
    if not intervals:
        return seq
    parts = []
    pos = 1
    for start, end in sorted(intervals):
        if start > pos:
            parts.append(seq[pos - 1:start - 1])
        pos = max(pos, end + 1)
    if pos <= len(seq):
        parts.append(seq[pos - 1:])
    return "".join(parts)


def write_bed(rows: List[Dict[str, str]], path: str | Path) -> None:
    with open(path, "w", encoding="utf-8") as out:
        for row in rows:
            start_1 = int(row["intron_start"])
            end_1 = int(row["intron_end"])
            out.write(f"{row['sequence_id']}\t{start_1 - 1}\t{end_1}\t{row['intron_id']}\n")


def detect_introns(
    input_fasta: str | Path,
    db: str | Path,
    outdir: str | Path,
    source_label: str = "query",
    vsearch: str = "vsearch",
    threads: int = 8,
    search_id: float = 0.70,
    rescue_id: float = 0.987,
    min_intron_len: int = 50,
    min_flank_len: int = 150,
    strand: str = "both",
    maxaccepts: int = 1,
    dry_run: bool = False,
) -> Dict[str, str]:
    input_fasta = ensure_file(input_fasta, "input FASTA")
    db = ensure_file(db, "reference DB/FASTA")
    outdir = ensure_dir(outdir)
    align = run_vsearch_alignment(input_fasta, db, outdir, vsearch, threads, search_id, strand, maxaccepts, dry_run)
    if dry_run:
        return {
            "alignment_rows": align["hits"],
            "analysis_sequences": str(outdir / "analysis_sequences.fa"),
            "sequence_version_map": str(outdir / "sequence_version_map.tsv"),
            "intron_summary": str(outdir / "intron_summary.tsv"),
            "intron_sequences": str(outdir / "intron_sequences.fa"),
            "intron_regions_bed": str(outdir / "intron_regions.bed"),
        }
    records = fasta_records_by_id(input_fasta)
    hits = read_alignment_rows(align["hits"])
    analysis_records: List[Tuple[str, str]] = []
    intron_records: List[Tuple[str, str]] = []
    summary_rows: List[Dict[str, str]] = []
    version_rows: List[Dict[str, str]] = []
    bed_rows: List[Dict[str, str]] = []
    for seq_id, (_header, seq) in records.items():
        hit = hits.get(seq_id)
        has_intron = "no"
        status = "no_hit" if hit is None else "no_intron_detected"
        analysis_id = seq_id
        analysis_seq = seq
        intron_positions: List[str] = []
        intron_lengths: List[str] = []
        rescued_identity = ""
        if hit is not None and hit.get("qrow") and hit.get("trow"):
            insertions = find_query_insertions(hit["qrow"], hit["trow"], min_intron_len)
            left_flank, right_flank = flank_lengths(insertions, len(seq))
            rescued = identity_excluding_insertions(hit["qrow"], hit["trow"], insertions)
            if insertions:
                intron_positions = [f'{x["query_start"]}-{x["query_end"]}' for x in insertions]
                intron_lengths = [str(x["length"]) for x in insertions]
                rescued_identity = f"{rescued:.6f}"
                if rescued >= rescue_id and left_flank >= min_flank_len and right_flank >= min_flank_len:
                    has_intron = "yes"
                    status = "intron_rescued"
                    analysis_id = f"{seq_id}|intron_free"
                    intervals = [(x["query_start"], x["query_end"]) for x in insertions]
                    analysis_seq = remove_intervals(seq, intervals)
                    for idx, ins in enumerate(insertions, start=1):
                        intron_id = f"{seq_id}|intron_{idx}|{ins['query_start']}-{ins['query_end']}"
                        intron_records.append((intron_id, seq[ins["query_start"] - 1:ins["query_end"]]))
                        bed_rows.append({
                            "sequence_id": seq_id,
                            "intron_id": f"intron_{idx}",
                            "intron_start": str(ins["query_start"]),
                            "intron_end": str(ins["query_end"]),
                        })
                else:
                    status = "candidate_insertion_not_rescued"
        analysis_records.append((analysis_id, analysis_seq))
        version_rows.append({
            "sequence_id": seq_id,
            "analysis_sequence_id": analysis_id,
            "source_label": source_label,
            "has_intron": has_intron,
            "original_length": str(len(seq)),
            "analysis_length": str(len(analysis_seq)),
            "centroid_sequence_allowed": "original_only",
        })
        summary_rows.append({
            "sequence_id": seq_id,
            "source_label": source_label,
            "subject_id": (hit or {}).get("target", ""),
            "original_identity": (hit or {}).get("id", ""),
            "original_length": str(len(seq)),
            "analysis_sequence_id": analysis_id,
            "analysis_length": str(len(analysis_seq)),
            "has_intron": has_intron,
            "n_introns": str(len(intron_positions)),
            "intron_positions_query": ",".join(intron_positions),
            "intron_lengths": ",".join(intron_lengths),
            "rescued_identity": rescued_identity,
            "status": status,
        })
    analysis_path = outdir / "analysis_sequences.fa"
    intron_seq_path = outdir / "intron_sequences.fa"
    summary_path = outdir / "intron_summary.tsv"
    version_path = outdir / "sequence_version_map.tsv"
    bed_path = outdir / "intron_regions.bed"
    write_fasta(analysis_records, analysis_path)
    write_fasta(intron_records, intron_seq_path)
    summary_fields = ["sequence_id", "source_label", "subject_id", "original_identity", "original_length", "analysis_sequence_id", "analysis_length", "has_intron", "n_introns", "intron_positions_query", "intron_lengths", "rescued_identity", "status"]
    version_fields = ["sequence_id", "analysis_sequence_id", "source_label", "has_intron", "original_length", "analysis_length", "centroid_sequence_allowed"]
    with open(summary_path, "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    with open(version_path, "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=version_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(version_rows)
    write_bed(bed_rows, bed_path)
    success(f"Intron detection finished: {outdir}")
    return {
        "alignment_rows": align["hits"],
        "analysis_sequences": str(analysis_path),
        "sequence_version_map": str(version_path),
        "intron_summary": str(summary_path),
        "intron_sequences": str(intron_seq_path),
        "intron_regions_bed": str(bed_path),
    }
