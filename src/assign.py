from __future__ import annotations

import csv
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

from .threads import validate_threads
from .utils import ensure_dir, ensure_file, read_fasta, write_fasta
from .vsearch import cluster as vsearch_cluster
from .vsearch import usearch_global


ASSIGNMENT_FIELDS = [
    "query",
    "cluster_id",
    "assignment_type",
    "identity",
    "qcov",
    "tcov",
    "source",
]


def _safe_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _coverage(start: str, end: str, length: str) -> str:
    start_f = _safe_float(start)
    end_f = _safe_float(end)
    length_f = _safe_float(length)

    if start_f is None or end_f is None or length_f is None or length_f <= 0:
        return ""

    covered = abs(end_f - start_f) + 1.0
    return f"{covered / length_f:.6f}"


def uc_assignments_from_old_hits(userout_tsv: str | Path) -> List[Dict[str, str]]:
    """Read VSEARCH global-search userout results.

    Supported formats:

    Legacy 12-column format:
      query target id alnlen mism opens qlo qhi tlo thi qcov tcov

    Current AutoTax2 14-column format:
      query target id alnlen mism opens qlo qhi tlo thi evalue bits ql tl

    The first hit per query is treated as the old-cluster assignment.
    """

    rows: List[Dict[str, str]] = []
    seen: Set[str] = set()

    path = Path(userout_tsv)

    if not path.exists() or path.stat().st_size == 0:
        return rows

    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue

            values = line.rstrip("\n").split("\t")

            if len(values) >= 14:
                query = values[0]
                target = values[1]
                identity = values[2]
                qlo = values[6]
                qhi = values[7]
                tlo = values[8]
                thi = values[9]
                qlen = values[12]
                tlen = values[13]
                qcov = _coverage(qlo, qhi, qlen)
                tcov = _coverage(tlo, thi, tlen)

            elif len(values) >= 12:
                query = values[0]
                target = values[1]
                identity = values[2]
                qcov = values[10]
                tcov = values[11]

            else:
                continue

            if not query or query in seen:
                continue

            seen.add(query)

            rows.append(
                {
                    "query": query,
                    "cluster_id": target,
                    "assignment_type": "old_cluster",
                    "identity": identity,
                    "qcov": qcov,
                    "tcov": tcov,
                    "source": "vsearch_usearch_global",
                }
            )

    return rows


def relabel_fasta(
    input_fasta: str | Path,
    output_fasta: str | Path,
    prefix: str,
) -> Dict[str, str]:
    """Relabel FASTA records as prefix1, prefix2, ... .

    Returns
    -------
    dict
        Mapping from old sequence ID to new centroid ID.
    """

    mapping: Dict[str, str] = {}
    records: List[Tuple[str, str]] = []

    for index, (header, seq) in enumerate(read_fasta(input_fasta), start=1):
        old_id = header.split()[0].split(";")[0]
        new_id = f"{prefix}{index}"
        mapping[old_id] = new_id
        records.append((new_id, seq))

    write_fasta(records, output_fasta)
    return mapping


def parse_vsearch_uc_for_new_clusters(
    uc_path: str | Path,
    old_to_new_centroid: Dict[str, str],
) -> List[Dict[str, str]]:
    """Parse VSEARCH UC output from clustering unmatched sequences.

    UC rows used:
      S row: seed/centroid sequence. Query sequence ID is column 9.
      H row: member sequence. Query sequence ID is column 9, target seed is column 10.
    """

    rows: List[Dict[str, str]] = []
    path = Path(uc_path)

    if not path.exists() or path.stat().st_size == 0:
        return rows

    cluster_num_to_seed: Dict[str, str] = {}

    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 9:
                continue

            record_type = parts[0]
            cluster_num = parts[1]
            identity = parts[3] if len(parts) > 3 else ""
            query = parts[8]

            if record_type == "S":
                cluster_num_to_seed[cluster_num] = query
                rows.append(
                    {
                        "query": query,
                        "cluster_id": old_to_new_centroid.get(query, query),
                        "assignment_type": "new_cluster_seed",
                        "identity": "100.0",
                        "qcov": "",
                        "tcov": "",
                        "source": "vsearch_cluster",
                    }
                )

            elif record_type == "H":
                target_seed = (
                    parts[9]
                    if len(parts) > 9
                    else cluster_num_to_seed.get(cluster_num, "")
                )

                rows.append(
                    {
                        "query": query,
                        "cluster_id": old_to_new_centroid.get(target_seed, target_seed),
                        "assignment_type": "new_cluster_member",
                        "identity": identity,
                        "qcov": "",
                        "tcov": "",
                        "source": "vsearch_cluster",
                    }
                )

    return rows


def write_assignment_table(
    rows: List[Dict[str, str]],
    output_tsv: str | Path,
) -> None:
    with open(output_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ASSIGNMENT_FIELDS, delimiter="\t")
        writer.writeheader()

        for row in rows:
            writer.writerow({key: row.get(key, "") for key in ASSIGNMENT_FIELDS})


def concatenate_fastas(
    paths: Iterable[str | Path],
    output: str | Path,
) -> None:
    with open(output, "w", encoding="utf-8") as out_handle:
        for path in paths:
            path = Path(path)

            if not path.exists() or path.stat().st_size == 0:
                continue

            content = path.read_text(encoding="utf-8")
            out_handle.write(content)

            if content and not content.endswith("\n"):
                out_handle.write("\n")


def fasta_ids(fasta: str | Path) -> Set[str]:
    ids: Set[str] = set()

    for header, _seq in read_fasta(fasta):
        seq_id = header.split()[0].split(";")[0]
        ids.add(seq_id)

    return ids


def write_unmatched_fasta(
    input_fasta: str | Path,
    assigned_ids: Set[str],
    output_fasta: str | Path,
) -> int:
    records: List[Tuple[str, str]] = []

    for header, seq in read_fasta(input_fasta):
        seq_id = header.split()[0].split(";")[0]

        if seq_id not in assigned_ids:
            records.append((header, seq))

    write_fasta(records, output_fasta)
    return len(records)


def assign_or_create(
    new_fasta: str | Path,
    old_centroids: str | Path,
    outdir: str | Path,
    identity: float,
    vsearch: str = "vsearch",
    threads: int = 4,
    strand: str = "both",
    new_cluster_prefix: str = "new_",
    cluster_method: str = "cluster_size",
    maxaccepts: int = 1,
) -> Dict[str, str]:
    """Assign new sequences to old clusters or create new clusters.

    Workflow:
      1. Search new sequences against old centroids.
      2. Assign matched sequences to old centroid IDs.
      3. Write unmatched sequences to unmatched.fasta.
      4. Cluster unmatched sequences de novo at the same identity threshold.
      5. Relabel new centroids with new_cluster_prefix.
      6. Concatenate old and new centroids into updated_centroids.fa.
      7. Write assignments.tsv.
    """

    if not 0 < identity <= 1:
        raise ValueError("identity must be > 0 and <= 1.")

    if maxaccepts < 1:
        raise ValueError("maxaccepts must be at least 1.")

    threads = validate_threads(threads)

    new_fasta = ensure_file(new_fasta, "new FASTA")
    old_centroids = ensure_file(old_centroids, "old centroids FASTA")
    outdir = ensure_dir(outdir)

    old_hit = usearch_global(
        input_fasta=str(new_fasta),
        db=str(old_centroids),
        outdir=str(outdir),
        output_prefix="new_vs_old",
        min_id=identity,
        maxaccepts=maxaccepts,
        strand=strand,
        vsearch=vsearch,
        threads=threads,
    )

    unmatched = outdir / "unmatched.fasta"
    new_centroids_raw = outdir / "new_centroids.raw.fasta"
    new_centroids_renamed = outdir / "new_centroids.fasta"
    new_clusters_uc = outdir / "new_clusters.uc"
    updated_centroids = outdir / "updated_centroids.fasta"
    assignment_tsv = outdir / "assignments.tsv"

    old_assignments = uc_assignments_from_old_hits(old_hit["userout"])
    assigned_ids = {row["query"] for row in old_assignments}

    write_unmatched_fasta(new_fasta, assigned_ids, unmatched)

    new_assignments: List[Dict[str, str]] = []

    if unmatched.exists() and unmatched.stat().st_size > 0:
        result = vsearch_cluster(
            input_fasta=str(unmatched),
            outdir=str(outdir),
            identity=identity,
            method=cluster_method,
            vsearch=vsearch,
            threads=threads,
            relabel=None,
        )

        shutil.move(result["centroids"], new_centroids_raw)
        shutil.move(result["uc"], new_clusters_uc)

        old_to_new = relabel_fasta(
            new_centroids_raw,
            new_centroids_renamed,
            new_cluster_prefix,
        )

        new_assignments = parse_vsearch_uc_for_new_clusters(
            new_clusters_uc,
            old_to_new,
        )

        concatenate_fastas(
            [old_centroids, new_centroids_renamed],
            updated_centroids,
        )

    else:
        new_centroids_raw.write_text("", encoding="utf-8")
        new_centroids_renamed.write_text("", encoding="utf-8")
        new_clusters_uc.write_text("", encoding="utf-8")

        concatenate_fastas([old_centroids], updated_centroids)

    all_rows = old_assignments + new_assignments
    write_assignment_table(all_rows, assignment_tsv)

    return {
        "old_assignments": old_hit["userout"],
        "old_assignments_blast6": old_hit["blast6out"],
        "unmatched": str(unmatched),
        "new_centroids": str(new_centroids_renamed),
        "new_clusters_uc": str(new_clusters_uc),
        "updated_centroids": str(updated_centroids),
        "assignments": str(assignment_tsv),
    }