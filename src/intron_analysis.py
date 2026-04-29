from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .external import (
    fimo_scan,
    mafft_align,
    meme_discover,
    rnaalifold_predict,
    rnafold_predict,
    streme_discover,
)
from .logging import step, success, warning
from .threads import validate_threads
from .utils import ensure_dir, ensure_file, read_fasta
from .vsearch import cluster as vsearch_cluster


FEATURE_FIELDS = [
    "intron_id",
    "query_id",
    "cluster_id",
    "region",
    "coordinate_site",
    "length",
    "gc_content",
    "first_20bp",
    "last_20bp",
    "confidence",
    "support_refs",
]

RNAFOLD_FIELDS = [
    "intron_id",
    "length",
    "mfe",
    "mfe_per_nt",
    "dot_bracket",
]

MOTIF_SUMMARY_FIELDS = [
    "tool",
    "output_dir",
    "motif_file",
    "fimo_dir",
    "status",
]

CLUSTER_SUMMARY_FIELDS = [
    "identity",
    "centroids",
    "uc",
]


@dataclass(frozen=True)
class IntronRecord:
    intron_id: str
    header: str
    sequence: str
    attributes: Dict[str, str]

    @property
    def query_id(self) -> str:
        if "query_id" in self.attributes:
            return self.attributes["query_id"]
        return self.intron_id.split("|", 1)[0]

    @property
    def cluster_id(self) -> str:
        return self.attributes.get("cluster", self.attributes.get("cluster_id", ""))

    @property
    def region(self) -> str:
        return self.attributes.get("region", "")

    @property
    def confidence(self) -> str:
        return self.attributes.get("confidence", "")

    @property
    def support_refs(self) -> str:
        return self.attributes.get("support_refs", "")


def parse_header_attributes(header: str) -> Dict[str, str]:
    """Parse key=value attributes from an intron FASTA header."""

    attrs: Dict[str, str] = {}

    for part in header.split("|")[1:]:
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        attrs[key.strip()] = value.strip()

    return attrs


def read_intron_records(path: str | Path) -> List[IntronRecord]:
    introns: List[IntronRecord] = []

    for header, seq in read_fasta(path):
        intron_id = header.split()[0]
        sequence = "".join(str(seq).split()).upper()
        introns.append(
            IntronRecord(
                intron_id=intron_id,
                header=header,
                sequence=sequence,
                attributes=parse_header_attributes(header),
            )
        )

    return introns


def gc_content(sequence: str) -> float:
    sequence = sequence.upper()
    valid = [base for base in sequence if base in {"A", "C", "G", "T", "U"}]

    if not valid:
        return 0.0

    gc = sum(1 for base in valid if base in {"G", "C"})
    return gc / len(valid)


def write_tsv(path: str | Path, fields: Sequence[str], rows: Iterable[Dict[str, object]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def read_optional_table(path: str | Path | None, key_fields: Sequence[str]) -> Dict[str, Dict[str, str]]:
    if path is None:
        return {}

    table = Path(path)
    if not table.exists() or table.stat().st_size == 0:
        return {}

    records: Dict[str, Dict[str, str]] = {}

    with table.open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = ""
            for field in key_fields:
                if row.get(field):
                    key = str(row[field])
                    break
            if key:
                records[key] = row

    return records


def feature_rows(
    introns: Sequence[IntronRecord],
    *,
    metadata: Dict[str, Dict[str, str]],
    coordinate_map: Dict[str, Dict[str, str]],
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []

    for record in introns:
        meta = metadata.get(record.cluster_id) or metadata.get(record.query_id) or {}
        coord = coordinate_map.get(record.cluster_id) or coordinate_map.get(record.query_id) or {}

        coordinate_site = (
            coord.get("coordinate_site")
            or coord.get("coordinate_left")
            or coord.get("coordinate_position")
            or ""
        )

        rows.append(
            {
                "intron_id": record.intron_id,
                "query_id": record.query_id,
                "cluster_id": record.cluster_id,
                "region": record.region or meta.get("query_intron_start", "") + "-" + meta.get("query_intron_end", "")
                if meta.get("query_intron_start") and meta.get("query_intron_end")
                else record.region,
                "coordinate_site": coordinate_site,
                "length": len(record.sequence),
                "gc_content": f"{gc_content(record.sequence):.6f}",
                "first_20bp": record.sequence[:20],
                "last_20bp": record.sequence[-20:],
                "confidence": record.confidence or meta.get("confidence", ""),
                "support_refs": record.support_refs or meta.get("clean_ref_support", meta.get("supporting_refs", "")),
            }
        )

    return rows


def parse_rnafold_output(path: str | Path) -> List[Dict[str, object]]:
    """Parse RNAfold FASTA-style output into a compact summary table."""

    output = Path(path)
    if not output.exists() or output.stat().st_size == 0:
        return []

    lines = [line.rstrip("\n") for line in output.read_text(encoding="utf-8", errors="replace").splitlines()]
    rows: List[Dict[str, object]] = []

    current_id = ""
    current_seq = ""

    mfe_re = re.compile(r"\(([-+]?\d+(?:\.\d+)?)\)\s*$")

    for line in lines:
        if not line:
            continue

        if line.startswith(">"):
            current_id = line[1:].split()[0]
            current_seq = ""
            continue

        if current_id and not current_seq:
            current_seq = line.strip().upper()
            continue

        if current_id and current_seq:
            structure_line = line.strip()
            match = mfe_re.search(structure_line)
            mfe = float(match.group(1)) if match else 0.0
            dot = structure_line.split()[0] if structure_line.split() else ""
            length = len(current_seq)
            rows.append(
                {
                    "intron_id": current_id,
                    "length": length,
                    "mfe": f"{mfe:.6f}",
                    "mfe_per_nt": f"{mfe / length:.6f}" if length else "",
                    "dot_bracket": dot,
                }
            )
            current_seq = ""

    return rows


def choose_existing_file(paths: Sequence[str | Path]) -> Optional[str]:
    for path in paths:
        p = Path(path)
        if p.exists() and p.stat().st_size > 0:
            return str(p)
    return None


def run_intron_clustering(
    introns_fasta: str | Path,
    outdir: str | Path,
    *,
    vsearch: str | Path,
    threads: int,
    identities: Sequence[float],
) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    clusters_dir = ensure_dir(Path(outdir) / "clusters")

    for identity in identities:
        label = f"{int(round(identity * 100)):03d}"
        level_dir = ensure_dir(clusters_dir / f"cluster_{label}")

        step(f"Clustering introns at identity {identity}")
        result = vsearch_cluster(
            input_fasta=introns_fasta,
            outdir=level_dir,
            identity=identity,
            method="cluster_size",
            vsearch=vsearch,
            threads=threads,
            relabel=f"intron_{label}_",
        )

        rows.append(
            {
                "identity": str(identity),
                "centroids": result["centroids"],
                "uc": result["uc"],
            }
        )

    return rows


def run_motif_analysis(
    introns_fasta: str | Path,
    outdir: str | Path,
    *,
    meme: str | Path | None,
    streme: str | Path | None,
    fimo: str | Path | None,
    threads: int,
    minw: int,
    maxw: int,
    nmotifs: int,
) -> List[Dict[str, str]]:
    motif_dir = ensure_dir(Path(outdir) / "motif")
    rows: List[Dict[str, str]] = []

    streme_motif_file = ""
    if streme:
        streme_dir = motif_dir / "streme"
        streme_discover(
            introns_fasta,
            streme_dir,
            streme=streme,
            minw=minw,
            maxw=maxw,
            dna=True,
        )
        streme_motif_file = choose_existing_file([streme_dir / "streme.txt", streme_dir / "streme.xml"]) or ""
        rows.append(
            {
                "tool": "streme",
                "output_dir": str(streme_dir),
                "motif_file": streme_motif_file,
                "fimo_dir": "",
                "status": "finished",
            }
        )

    meme_motif_file = ""
    if meme:
        meme_dir = motif_dir / "meme"
        meme_discover(
            introns_fasta,
            meme_dir,
            meme=meme,
            threads=threads,
            nmotifs=nmotifs,
            minw=minw,
            maxw=maxw,
            dna=True,
        )
        meme_motif_file = choose_existing_file([meme_dir / "meme.txt", meme_dir / "meme.xml"]) or ""
        rows.append(
            {
                "tool": "meme",
                "output_dir": str(meme_dir),
                "motif_file": meme_motif_file,
                "fimo_dir": "",
                "status": "finished",
            }
        )

    if fimo:
        for tool_name, motif_file in [("streme", streme_motif_file), ("meme", meme_motif_file)]:
            if not motif_file:
                continue
            fimo_dir = motif_dir / f"fimo_{tool_name}"
            fimo_scan(
                motif_file,
                introns_fasta,
                fimo_dir,
                fimo=fimo,
            )
            rows.append(
                {
                    "tool": f"fimo_{tool_name}",
                    "output_dir": "",
                    "motif_file": motif_file,
                    "fimo_dir": str(fimo_dir),
                    "status": "finished",
                }
            )

    if not rows:
        rows.append(
            {
                "tool": "motif",
                "output_dir": str(motif_dir),
                "motif_file": "",
                "fimo_dir": "",
                "status": "skipped_no_motif_tools",
            }
        )

    return rows


def analyze_introns(
    introns_fasta: str | Path,
    outdir: str | Path,
    *,
    metadata: str | Path | None = None,
    coordinate_map: str | Path | None = None,
    vsearch: str | Path | None = None,
    mafft: str | Path | None = None,
    meme: str | Path | None = None,
    streme: str | Path | None = None,
    fimo: str | Path | None = None,
    rnafold: str | Path | None = None,
    rnaalifold: str | Path | None = None,
    threads: int = 4,
    run_cluster: bool = True,
    run_mafft: bool = True,
    run_motif: bool = True,
    run_structure: bool = True,
    cluster_identities: Sequence[float] = (0.99, 0.97),
    motif_minw: int = 6,
    motif_maxw: int = 50,
    motif_nmotifs: int = 10,
) -> Dict[str, str]:
    """Characterize confirmed intron-like insertion sequences.

    This command is intentionally downstream of detect-intron. It does not
    decide whether a query sequence contains an intron; it analyzes the
    sequences that were already confirmed and extracted.
    """

    threads = validate_threads(threads)
    introns_fasta = ensure_file(introns_fasta, "intron FASTA")
    outdir = ensure_dir(outdir)

    introns = read_intron_records(introns_fasta)
    if not introns:
        raise ValueError(f"No intron records found in {introns_fasta}")

    metadata_rows = read_optional_table(metadata, key_fields=["cluster_id", "query_id", "intron_id"])
    coordinate_rows = read_optional_table(coordinate_map, key_fields=["cluster_id", "query_id", "intron_id"])

    feature_dir = ensure_dir(outdir / "summary")
    feature_tsv = feature_dir / "intron_features.tsv"

    step("Writing intron feature table")
    write_tsv(
        feature_tsv,
        FEATURE_FIELDS,
        feature_rows(introns, metadata=metadata_rows, coordinate_map=coordinate_rows),
    )

    outputs: Dict[str, str] = {
        "intron_features": str(feature_tsv),
    }

    if run_cluster:
        if not vsearch:
            warning("Skipping intron clustering because vsearch was not provided.")
        else:
            cluster_rows = run_intron_clustering(
                introns_fasta,
                outdir,
                vsearch=vsearch,
                threads=threads,
                identities=cluster_identities,
            )
            cluster_summary = Path(outdir) / "clusters" / "intron_cluster_summary.tsv"
            write_tsv(cluster_summary, CLUSTER_SUMMARY_FIELDS, cluster_rows)
            outputs["intron_cluster_summary"] = str(cluster_summary)

    alignment_fasta = Path(outdir) / "alignment" / "introns.mafft.fasta"

    if run_mafft:
        if not mafft:
            warning("Skipping MAFFT alignment because mafft was not provided.")
        else:
            mafft_align(
                introns_fasta,
                alignment_fasta,
                mafft=mafft,
                threads=threads,
                auto=True,
            )
            outputs["mafft_alignment"] = str(alignment_fasta)

    if run_motif:
        motif_rows = run_motif_analysis(
            introns_fasta,
            outdir,
            meme=meme,
            streme=streme,
            fimo=fimo,
            threads=threads,
            minw=motif_minw,
            maxw=motif_maxw,
            nmotifs=motif_nmotifs,
        )
        motif_summary = Path(outdir) / "motif" / "motif_summary.tsv"
        write_tsv(motif_summary, MOTIF_SUMMARY_FIELDS, motif_rows)
        outputs["motif_summary"] = str(motif_summary)

    if run_structure:
        structure_dir = ensure_dir(outdir / "structure")

        if rnafold:
            rnafold_txt = structure_dir / "introns.RNAfold.txt"
            rnafold_predict(
                introns_fasta,
                rnafold_txt,
                rnafold=rnafold,
            )
            rnafold_summary = structure_dir / "introns.RNAfold_summary.tsv"
            write_tsv(rnafold_summary, RNAFOLD_FIELDS, parse_rnafold_output(rnafold_txt))
            outputs["rnafold_output"] = str(rnafold_txt)
            outputs["rnafold_summary"] = str(rnafold_summary)
        else:
            warning("Skipping RNAfold because RNAfold was not provided.")

        if rnaalifold:
            if alignment_fasta.exists() and alignment_fasta.stat().st_size > 0:
                rnaalifold_txt = structure_dir / "introns.RNAalifold.txt"
                rnaalifold_predict(
                    alignment_fasta,
                    rnaalifold_txt,
                    rnaalifold=rnaalifold,
                )
                outputs["rnaalifold_output"] = str(rnaalifold_txt)
            else:
                warning("Skipping RNAalifold because no MAFFT alignment is available.")
        else:
            warning("Skipping RNAalifold because RNAalifold was not provided.")

    success(f"Intron analysis finished: {outdir}")
    return outputs
