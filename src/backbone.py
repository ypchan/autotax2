from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

from .assign import parse_vsearch_uc_for_new_clusters, relabel_fasta
from .logging import step, success
from .ranks import parent_taxon_id, parse_rank_thresholds, rank_taxon_id
from .threads import validate_threads
from .utils import ensure_dir, ensure_file, parse_manifest, read_fasta, write_fasta
from .vsearch import cluster as vsearch_cluster
from .vsearch import usearch_global


ASSIGNMENT_FIELDS = [
    "sequence_id",
    "analysis_sequence_id",
    "original_sequence_id",
    "source_label",
    "rank",
    "rank_label",
    "identity_threshold",
    "assignment_type",
    "taxon_id",
    "parent_taxon_id",
    "cluster_id",
    "silva_target",
    "silva_identity",
    "silva_qcov",
    "silva_tcov",
    "is_novel",
    "has_intron",
]

SUMMARY_FIELDS = [
    "rank",
    "rank_label",
    "identity_threshold",
    "taxon_id",
    "source_label",
    "n_sequences",
    "is_novel",
    "assignment_types_json",
]


def read_silva_taxonomy(path: str | Path) -> Dict[str, Dict[str, str]]:
    table = Path(path)

    if not table.exists() or table.stat().st_size == 0:
        raise FileNotFoundError(f"SILVA taxonomy table not found: {table}")

    taxonomy: Dict[str, Dict[str, str]] = {}

    with table.open(encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            seq_id = row.get("ID", "")

            if not seq_id:
                continue

            taxonomy[seq_id] = {
                key: value
                for key, value in row.items()
                if key != "ID" and value
            }

    return taxonomy


def read_version_map(path: str | Path | None) -> Dict[str, Dict[str, str]]:
    if not path:
        return {}

    table = Path(path)

    if not table.exists() or table.stat().st_size == 0:
        return {}

    version_map: Dict[str, Dict[str, str]] = {}

    with table.open(encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            analysis_id = (
                row.get("analysis_sequence_id")
                or row.get("analysis_id")
                or row.get("sequence_id")
                or row.get("id")
            )

            if analysis_id:
                version_map[analysis_id] = row

    return version_map


def fasta_dict(path: str | Path) -> Dict[str, Tuple[str, str]]:
    return {
        header.split()[0].split(";")[0]: (header, seq)
        for header, seq in read_fasta(path)
    }


def fasta_ids(path: str | Path) -> Set[str]:
    return {
        header.split()[0].split(";")[0]
        for header, _seq in read_fasta(path)
    }


def safe_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def coverage(start: str, end: str, length: str) -> str:
    start_f = safe_float(start)
    end_f = safe_float(end)
    length_f = safe_float(length)

    if start_f is None or end_f is None or length_f is None or length_f <= 0:
        return ""

    covered = abs(end_f - start_f) + 1.0
    return f"{covered / length_f:.6f}"


def parse_hit_table(userout_tsv: str | Path) -> Dict[str, Dict[str, str]]:
    """Read first VSEARCH global-search hit per query.

    Supported formats:

    Legacy 12-column format:
      query target id alnlen mism opens qlo qhi tlo thi qcov tcov

    Current AutoTax2 14-column format:
      query target id alnlen mism opens qlo qhi tlo thi evalue bits ql tl
    """

    hits: Dict[str, Dict[str, str]] = {}
    table = Path(userout_tsv)

    if not table.exists() or table.stat().st_size == 0:
        return hits

    with table.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue

            values = line.rstrip("\n").split("\t")

            if len(values) >= 14:
                row = {
                    "query": values[0],
                    "target": values[1],
                    "id": values[2],
                    "alnlen": values[3],
                    "mism": values[4],
                    "opens": values[5],
                    "qlo": values[6],
                    "qhi": values[7],
                    "tlo": values[8],
                    "thi": values[9],
                    "evalue": values[10],
                    "bits": values[11],
                    "ql": values[12],
                    "tl": values[13],
                }
                row["qcov"] = coverage(row["qlo"], row["qhi"], row["ql"])
                row["tcov"] = coverage(row["tlo"], row["thi"], row["tl"])

            elif len(values) >= 12:
                row = {
                    "query": values[0],
                    "target": values[1],
                    "id": values[2],
                    "alnlen": values[3],
                    "mism": values[4],
                    "opens": values[5],
                    "qlo": values[6],
                    "qhi": values[7],
                    "tlo": values[8],
                    "thi": values[9],
                    "qcov": values[10],
                    "tcov": values[11],
                }

            else:
                continue

            query = row.get("query", "")

            if query and query not in hits:
                hits[query] = row

    return hits


def write_rows(
    path: str | Path,
    rows: Iterable[Dict[str, str]],
    fields: List[str],
) -> None:
    with open(path, "w", encoding="utf-8", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()

        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def choose_db(
    manifest: Dict[str, str],
    key_fasta: str,
    key_udb: str,
    db_format: str,
) -> str:
    db_format = db_format.lower()

    if db_format == "udb":
        return manifest.get(key_udb) or manifest[key_fasta]

    if db_format == "fasta":
        return manifest[key_fasta]

    if db_format != "auto":
        raise ValueError("db_format must be one of: auto, udb, fasta.")

    return manifest.get(key_udb) or manifest[key_fasta]


def write_unmatched_fasta(
    input_fasta: str | Path,
    matched_ids: Set[str],
    output_fasta: str | Path,
) -> int:
    records: List[Tuple[str, str]] = []

    for header, seq in read_fasta(input_fasta):
        seq_id = header.split()[0].split(";")[0]

        if seq_id not in matched_ids:
            records.append((header, seq))

    write_fasta(records, output_fasta)
    return len(records)


def make_official_centroids(
    assignment_rows: List[Dict[str, str]],
    original_records: Dict[str, Tuple[str, str]],
    analysis_records: Dict[str, Tuple[str, str]],
    out_fasta: str | Path,
) -> None:
    taxon_to_seq: Dict[str, str] = {}

    for row in assignment_rows:
        taxon = row["taxon_id"]

        if taxon not in taxon_to_seq:
            taxon_to_seq[taxon] = (
                row.get("original_sequence_id")
                or row.get("sequence_id")
                or row.get("analysis_sequence_id")
                or ""
            )

    records: List[Tuple[str, str]] = []

    for taxon, seq_id in taxon_to_seq.items():
        if seq_id in original_records:
            records.append((taxon, original_records[seq_id][1]))
        elif seq_id in analysis_records:
            records.append((taxon, analysis_records[seq_id][1]))

    write_fasta(records, out_fasta)


def original_id_from_version_map(
    seq_id: str,
    version_map: Dict[str, Dict[str, str]],
) -> str:
    row = version_map.get(seq_id, {})

    return (
        row.get("original_sequence_id")
        or row.get("original_id")
        or row.get("sequence_id")
        or seq_id
    )


def intron_status_from_version_map(
    seq_id: str,
    version_map: Dict[str, Dict[str, str]],
) -> str:
    row = version_map.get(seq_id, {})

    return (
        row.get("has_intron")
        or row.get("status")
        or ("unknown" if version_map else "no")
    )


def insert_backbone(
    input_fasta: str | Path,
    manifest_path: str | Path,
    outdir: str | Path,
    source_label: str,
    rank_thresholds: str | None = "default",
    original_fasta: str | Path | None = None,
    version_map_path: str | Path | None = None,
    db_format: str = "auto",
    vsearch: str = "vsearch",
    threads: int = 4,
    strand: str = "both",
    maxaccepts: int = 1,
    centroid_priority: str = "original_seed",
    centroid_policy: str | None = None,
) -> Dict[str, str]:
    """Insert extension sequences into the SILVA backbone.

    For each rank threshold:
      1. Search input sequences against SILVA.
      2. Assign matched sequences to existing SILVA taxon IDs.
      3. Write unmatched sequences.
      4. Cluster unmatched sequences as novel taxa.
      5. Write rank-level assignment and summary tables.

    centroid_priority is currently accepted for CLI clarity and future use.
    """

    if maxaccepts < 1:
        raise ValueError("maxaccepts must be at least 1.")

    threads = validate_threads(threads)

    input_fasta = ensure_file(input_fasta, "input FASTA")
    manifest_path = ensure_file(manifest_path, "SILVA manifest")
    outdir = ensure_dir(outdir)

    manifest = parse_manifest(manifest_path)

    if "silva_taxonomy_tsv" not in manifest:
        raise KeyError("Reference manifest is missing 'silva_taxonomy_tsv'.")

    taxonomy = read_silva_taxonomy(manifest["silva_taxonomy_tsv"])
    ranks = parse_rank_thresholds(rank_thresholds)
    db = choose_db(manifest, "silva_fasta", "silva_udb", db_format)

    rank_uc_dir = ensure_dir(outdir / "rank_uc")
    rank_centroids_core = ensure_dir(outdir / "rank_centroids_core")
    rank_centroids_original = ensure_dir(outdir / "rank_centroids_original")
    rank_hits = ensure_dir(outdir / "rank_hits")
    rank_unmatched = ensure_dir(outdir / "rank_unmatched")
    rank_novel_clusters = ensure_dir(outdir / "rank_novel_clusters")

    analysis_records = fasta_dict(input_fasta)
    original_records = fasta_dict(original_fasta) if original_fasta else analysis_records
    version_map = read_version_map(version_map_path)

    all_assignment_rows: List[Dict[str, str]] = []
    all_summary_rows: List[Dict[str, str]] = []

    for rank_threshold in ranks:
        rank = rank_threshold.rank
        threshold = rank_threshold.threshold

        step(f"Backbone insertion: rank={rank}, threshold={threshold}")

        search_out = usearch_global(
            input_fasta=str(input_fasta),
            db=db,
            outdir=str(rank_hits),
            output_prefix=f"{rank}_vs_silva",
            min_id=threshold,
            maxaccepts=maxaccepts,
            strand=strand,
            vsearch=vsearch,
            threads=threads,
        )

        hits = parse_hit_table(search_out["userout"])
        matched_ids = set(hits)

        rank_notmatched = rank_unmatched / f"{rank}.unmatched.fasta"
        write_unmatched_fasta(input_fasta, matched_ids, rank_notmatched)

        assignment_rows: List[Dict[str, str]] = []

        for seq_id, hit in hits.items():
            target_id = hit.get("target", "").split()[0].split(";")[0]
            tax = taxonomy.get(target_id, {})

            taxon = (
                rank_taxon_id(rank, tax, prefix="silva")
                or f"silva|{rank}|target:{target_id}"
            )
            parent = parent_taxon_id(rank, tax, prefix="silva")

            original_id = original_id_from_version_map(seq_id, version_map)
            has_intron = intron_status_from_version_map(seq_id, version_map)

            assignment_rows.append(
                {
                    "sequence_id": seq_id,
                    "analysis_sequence_id": seq_id,
                    "original_sequence_id": original_id,
                    "source_label": source_label,
                    "rank": rank,
                    "rank_label": rank_threshold.label,
                    "identity_threshold": str(threshold),
                    "assignment_type": "silva_existing",
                    "taxon_id": taxon,
                    "parent_taxon_id": parent,
                    "cluster_id": taxon,
                    "silva_target": target_id,
                    "silva_identity": hit.get("id", ""),
                    "silva_qcov": hit.get("qcov", ""),
                    "silva_tcov": hit.get("tcov", ""),
                    "is_novel": "no",
                    "has_intron": has_intron,
                }
            )

        if rank_notmatched.exists() and rank_notmatched.stat().st_size > 0:
            cluster_dir = ensure_dir(rank_novel_clusters / rank)

            result = vsearch_cluster(
                input_fasta=str(rank_notmatched),
                outdir=str(cluster_dir),
                identity=threshold,
                method="cluster_size",
                vsearch=vsearch,
                threads=threads,
                relabel=None,
            )

            raw_centroids = Path(result["centroids"])
            raw_uc = Path(result["uc"])

            novel_uc = rank_uc_dir / f"{rank}.novel.uc"
            core_centroids = rank_centroids_core / f"{rank}.centroids.core.fasta"

            raw_uc.replace(novel_uc)

            prefix = f"novel|{rank}|{source_label}_{rank}_"
            old_to_new = relabel_fasta(raw_centroids, core_centroids, prefix)

            parsed = parse_vsearch_uc_for_new_clusters(novel_uc, old_to_new)

            for row in parsed:
                seq_id = row["query"]
                taxon = row["cluster_id"]
                original_id = original_id_from_version_map(seq_id, version_map)
                has_intron = intron_status_from_version_map(seq_id, version_map)

                assignment_rows.append(
                    {
                        "sequence_id": seq_id,
                        "analysis_sequence_id": seq_id,
                        "original_sequence_id": original_id,
                        "source_label": source_label,
                        "rank": rank,
                        "rank_label": rank_threshold.label,
                        "identity_threshold": str(threshold),
                        "assignment_type": row["assignment_type"],
                        "taxon_id": taxon,
                        "parent_taxon_id": "",
                        "cluster_id": taxon,
                        "silva_target": "",
                        "silva_identity": "",
                        "silva_qcov": "",
                        "silva_tcov": "",
                        "is_novel": "yes",
                        "has_intron": has_intron,
                    }
                )

        else:
            (rank_uc_dir / f"{rank}.novel.uc").write_text("", encoding="utf-8")
            (rank_centroids_core / f"{rank}.centroids.core.fasta").write_text(
                "",
                encoding="utf-8",
            )

        # New usearch_global does not produce UC output.
        # Keep an empty placeholder for downstream directory consistency.
        (rank_uc_dir / f"{rank}.silva.uc").write_text("", encoding="utf-8")

        make_official_centroids(
            assignment_rows,
            original_records,
            analysis_records,
            rank_centroids_original / f"{rank}.centroids.original.fasta",
        )

        counts = Counter(row["taxon_id"] for row in assignment_rows)
        assignment_type_counts = defaultdict(Counter)
        novel_counts = defaultdict(int)

        for row in assignment_rows:
            taxon = row["taxon_id"]
            assignment_type_counts[taxon][row["assignment_type"]] += 1

            if row["is_novel"] == "yes":
                novel_counts[taxon] += 1

        for taxon, count in sorted(counts.items()):
            all_summary_rows.append(
                {
                    "rank": rank,
                    "rank_label": rank_threshold.label,
                    "identity_threshold": str(threshold),
                    "taxon_id": taxon,
                    "source_label": source_label,
                    "n_sequences": str(count),
                    "is_novel": "yes" if novel_counts[taxon] else "no",
                    "assignment_types_json": str(dict(assignment_type_counts[taxon])),
                }
            )

        all_assignment_rows.extend(assignment_rows)

    assignment_path = outdir / "sequence_rank_assignment.tsv"
    summary_path = outdir / "rank_taxa_summary.tsv"

    write_rows(assignment_path, all_assignment_rows, ASSIGNMENT_FIELDS)
    write_rows(summary_path, all_summary_rows, SUMMARY_FIELDS)

    success(f"Backbone insertion finished: {outdir}")

    return {
        "sequence_rank_assignment": str(assignment_path),
        "rank_taxa_summary": str(summary_path),
        "rank_uc_dir": str(rank_uc_dir),
        "rank_centroids_core": str(rank_centroids_core),
        "rank_centroids_original": str(rank_centroids_original),
    }