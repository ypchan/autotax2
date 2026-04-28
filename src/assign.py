from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Iterable, List

from .utils import ensure_dir, ensure_file, read_fasta, write_fasta
from .vsearch import cluster as vsearch_cluster
from .vsearch import usearch_global


def uc_assignments_from_old_hits(userout_tsv: str | Path) -> List[Dict[str, str]]:
    """Read VSEARCH --userout results for new-vs-old centroid assignment.

    Expected fields:
      query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+qcov+tcov

    AutoTax2 assign uses --maxaccepts 1 by default, so the first hit per query
    is treated as the final old-cluster assignment.
    """
    fields = ["query", "target", "id", "alnlen", "mism", "opens", "qlo", "qhi", "tlo", "thi", "qcov", "tcov"]
    seen = set()
    rows: List[Dict[str, str]] = []
    p = Path(userout_tsv)
    if not p.exists() or p.stat().st_size == 0:
        return rows

    with open(p, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            values = line.rstrip("\n").split("\t")
            row = dict(zip(fields, values))
            q = row.get("query", "")
            if not q or q in seen:
                continue
            seen.add(q)
            rows.append({
                "query": q,
                "cluster_id": row.get("target", ""),
                "assignment_type": "old_cluster",
                "identity": row.get("id", ""),
                "qcov": row.get("qcov", ""),
                "tcov": row.get("tcov", ""),
                "source": "vsearch_usearch_global",
            })
    return rows


def relabel_fasta(input_fasta: str | Path, output_fasta: str | Path, prefix: str) -> Dict[str, str]:
    """Relabel FASTA records as prefix1, prefix2... and return old_id -> new_id."""
    mapping: Dict[str, str] = {}
    records = []
    for i, (header, seq) in enumerate(read_fasta(input_fasta), start=1):
        old_id = header.split()[0].split(";")[0]
        new_id = f"{prefix}{i}"
        mapping[old_id] = new_id
        records.append((new_id, seq))
    write_fasta(records, output_fasta)
    return mapping


def parse_vsearch_uc_for_new_clusters(uc_path: str | Path, old_to_new_centroid: Dict[str, str]) -> List[Dict[str, str]]:
    """Parse VSEARCH .uc file from clustering unmatched sequences.

    UC rows used:
    - S row: seed/centroid sequence. Query sequence ID is column 9.
    - H row: member sequence. Query sequence ID is column 9, target seed is column 10.
    - Column 4 contains percent identity for H rows.
    """
    rows: List[Dict[str, str]] = []
    p = Path(uc_path)
    if not p.exists() or p.stat().st_size == 0:
        return rows

    cluster_num_to_seed: Dict[str, str] = {}

    with open(p, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            rec_type = parts[0]
            cluster_num = parts[1]
            identity = parts[3] if len(parts) > 3 else ""
            query = parts[8]

            if rec_type == "S":
                cluster_num_to_seed[cluster_num] = query
                rows.append({
                    "query": query,
                    "cluster_id": old_to_new_centroid.get(query, query),
                    "assignment_type": "new_cluster_seed",
                    "identity": "100.0",
                    "qcov": "",
                    "tcov": "",
                    "source": "vsearch_cluster",
                })
            elif rec_type == "H":
                target_seed = parts[9] if len(parts) > 9 else cluster_num_to_seed.get(cluster_num, "")
                rows.append({
                    "query": query,
                    "cluster_id": old_to_new_centroid.get(target_seed, target_seed),
                    "assignment_type": "new_cluster_member",
                    "identity": identity,
                    "qcov": "",
                    "tcov": "",
                    "source": "vsearch_cluster",
                })
    return rows


def write_assignment_table(rows: List[Dict[str, str]], output_tsv: str | Path) -> None:
    fields = ["query", "cluster_id", "assignment_type", "identity", "qcov", "tcov", "source"]
    with open(output_tsv, "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})


def concatenate_fastas(paths: Iterable[str | Path], output: str | Path) -> None:
    with open(output, "w", encoding="utf-8") as out:
        for p in paths:
            p = Path(p)
            if not p.exists() or p.stat().st_size == 0:
                continue
            with open(p, encoding="utf-8") as inp:
                content = inp.read()
                out.write(content)
                if content and not content.endswith("\n"):
                    out.write("\n")


def assign_or_create(
    new_fasta: str | Path,
    old_centroids: str | Path,
    outdir: str | Path,
    identity: float,
    vsearch: str = "vsearch",
    threads: int = 8,
    strand: str = "both",
    new_cluster_prefix: str = "new_",
    cluster_method: str = "cluster_size",
    maxaccepts: int = 1,
    dry_run: bool = False,
) -> Dict[str, str]:
    """Assign new sequences to old clusters or create new clusters.

    Steps:
    1. Run VSEARCH global search of new sequences against old centroids.
    2. Matched sequences are assigned to old centroid IDs.
    3. Unmatched sequences are clustered de novo at the same identity threshold.
    4. New centroids are relabeled with `new_cluster_prefix`.
    5. Old and new centroids are concatenated into `updated_centroids.fa`.
    6. A combined `assignments.tsv` is written.
    """
    new_fasta = ensure_file(new_fasta, "new FASTA")
    old_centroids = ensure_file(old_centroids, "old centroids FASTA")
    outdir = ensure_dir(outdir)

    old_hit = usearch_global(
        query_fasta=str(new_fasta),
        db_fasta=str(old_centroids),
        outdir=str(outdir),
        prefix="new_vs_old",
        identity=identity,
        maxaccepts=maxaccepts,
        strand=strand,
        vsearch=vsearch,
        threads=threads,
        notmatched=True,
        dry_run=dry_run,
    )

    unmatched = Path(old_hit["notmatched"])
    new_centroids_raw = outdir / "new_centroids.raw.fa"
    new_centroids_renamed = outdir / "new_centroids_renamed.fa"
    new_clusters_uc = outdir / "new_clusters.uc"
    updated_centroids = outdir / "updated_centroids.fa"
    assignment_tsv = outdir / "assignments.tsv"

    if dry_run:
        vsearch_cluster(
            input_fasta=str(unmatched),
            outdir=str(outdir),
            identity=identity,
            method=cluster_method,
            vsearch=vsearch,
            threads=threads,
            dry_run=True,
        )
        return {
            "old_assignments": old_hit["userout"],
            "unmatched": str(unmatched),
            "new_centroids": str(new_centroids_renamed),
            "updated_centroids": str(updated_centroids),
            "assignments": str(assignment_tsv),
        }

    old_assignments = uc_assignments_from_old_hits(old_hit["userout"])
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
            dry_run=False,
        )

        Path(result["centroids"]).replace(new_centroids_raw)
        Path(result["uc"]).replace(new_clusters_uc)

        old_to_new = relabel_fasta(new_centroids_raw, new_centroids_renamed, new_cluster_prefix)
        new_assignments = parse_vsearch_uc_for_new_clusters(new_clusters_uc, old_to_new)

        concatenate_fastas([old_centroids, new_centroids_renamed], updated_centroids)
    else:
        new_centroids_raw.write_text("", encoding="utf-8")
        new_centroids_renamed.write_text("", encoding="utf-8")
        new_clusters_uc.write_text("", encoding="utf-8")
        concatenate_fastas([old_centroids], updated_centroids)

    all_rows = old_assignments + new_assignments
    write_assignment_table(all_rows, assignment_tsv)

    return {
        "old_assignments": old_hit["userout"],
        "old_assignments_uc": old_hit["uc"],
        "unmatched": str(unmatched),
        "new_centroids": str(new_centroids_renamed),
        "new_clusters_uc": str(new_clusters_uc),
        "updated_centroids": str(updated_centroids),
        "assignments": str(assignment_tsv),
    }
