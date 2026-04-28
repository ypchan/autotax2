from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from .vsearch import make_udb
from .utils import copy_or_decompress, ensure_dir, ensure_file, read_fasta, strip_fasta_suffix, write_fasta, write_manifest


BAD_TAXON_RE = re.compile(
    r"uncultured|unknown|unidentified|incertae sedis|metagenome|\bbacterium\b|\bpossible\b",
    flags=re.IGNORECASE,
)

SINTAX_PREFIXES = ["d:", "p:", "c:", "o:", "f:", "g:", "s:"]
TAX_COLS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def clean_taxon(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    value = value.strip()
    if not value:
        return None
    value = re.sub(r"[,;]", ".", value)
    if BAD_TAXON_RE.search(value):
        return None
    return value


def parse_silva_header(header: str) -> Tuple[str, List[Optional[str]]]:
    """Parse SILVA header: '<accession> <tax1;tax2;...>'."""
    parts = header.split(" ", 1)
    seq_id = parts[0]
    if len(parts) == 1:
        taxa: List[Optional[str]] = []
    else:
        taxa = [clean_taxon(x) for x in parts[1].split(";")]
    taxa = (taxa + [None] * len(TAX_COLS))[:len(TAX_COLS)]
    return seq_id, taxa


def make_sintax_header(seq_id: str, taxa: List[Optional[str]]) -> str:
    pairs = []
    for prefix, taxon in zip(SINTAX_PREFIXES, taxa):
        if taxon:
            pairs.append(prefix + taxon)
    return f"{seq_id};tax={','.join(pairs)};"


def detect_metadata_columns(header: List[str]) -> Tuple[int, int]:
    lowered = [h.lower() for h in header]
    if "acc" not in lowered:
        raise ValueError("SILVA metadata must contain an 'acc' column.")
    if "flags" not in lowered:
        raise ValueError("SILVA metadata must contain a 'flags' column.")
    return lowered.index("acc"), lowered.index("flags")


def extract_typestrain_accessions(metadata_path: str | Path) -> List[str]:
    """Extract accession IDs whose flags column contains [T] or [t]."""
    ids: List[str] = []
    with open(metadata_path, encoding="utf-8", errors="replace") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        acc_i, flags_i = detect_metadata_columns(header)
        for row in reader:
            if len(row) <= max(acc_i, flags_i):
                continue
            acc = row[acc_i].strip()
            flags = row[flags_i]
            if acc and re.search(r"\[t\]", flags, flags=re.IGNORECASE):
                ids.append(acc)
    return ids


def prepare_silva(
    silva_fasta: str | Path,
    silva_metadata: str | Path,
    outdir: str | Path,
    prefix: Optional[str] = None,
    keep_decompressed: bool = True,
    make_udb_files: bool = False,
    vsearch: str = "vsearch",
    dry_run: bool = False,
) -> Dict[str, str]:
    """Prepare local SILVA FASTA/metadata for AutoTax2 + VSEARCH.

    No downloading is performed.
    """
    silva_fasta = ensure_file(silva_fasta, "SILVA FASTA")
    silva_metadata = ensure_file(silva_metadata, "SILVA metadata")
    outdir = ensure_dir(outdir)

    if prefix is None:
        prefix = strip_fasta_suffix(silva_fasta.name)

    out_fasta = outdir / f"{prefix}.fasta"
    out_sintax = outdir / f"{prefix}_sintax.fasta"
    out_typestrains = outdir / f"{prefix}_typestrains.fasta"
    out_ids = outdir / "typestrains_accessionIDs.txt"
    out_tax = outdir / "silva_taxonomy.tsv"
    manifest = outdir / "autotax2_ref_manifest.tsv"

    # Normalize input FASTA and metadata into plain text local files when needed.
    copy_or_decompress(silva_fasta, out_fasta)

    metadata_plain = outdir / f"{strip_fasta_suffix(silva_metadata.name)}"
    copy_or_decompress(silva_metadata, metadata_plain)

    type_ids = extract_typestrain_accessions(metadata_plain)
    with open(out_ids, "w", encoding="utf-8") as f:
        for acc in type_ids:
            f.write(acc + "\n")

    type_ids_list = list(type_ids)

    sintax_records = []
    type_records = []
    with open(out_tax, "w", encoding="utf-8") as taxout:
        taxout.write("ID\t" + "\t".join(TAX_COLS) + "\n")
        for header, seq in read_fasta(out_fasta):
            seq_id, taxa = parse_silva_header(header)
            taxout.write(seq_id + "\t" + "\t".join(t if t else "" for t in taxa) + "\n")
            sintax_records.append((make_sintax_header(seq_id, taxa), seq))
            if any(acc in header for acc in type_ids_list):
                type_records.append((header, seq))

    write_fasta(sintax_records, out_sintax)
    write_fasta(type_records, out_typestrains)

    entries = {
        "silva_fasta": str(out_fasta),
        "silva_sintax_fasta": str(out_sintax),
        "silva_typestrains_fasta": str(out_typestrains),
        "silva_metadata": str(metadata_plain),
        "typestrains_accession_ids": str(out_ids),
        "silva_taxonomy_tsv": str(out_tax),
    }

    if make_udb_files:
        out_udb = outdir / f"{prefix}.udb"
        out_sintax_udb = outdir / f"{prefix}_sintax.udb"
        out_typestrains_udb = outdir / f"{prefix}_typestrains.udb"
        make_udb(str(out_fasta), str(out_udb), vsearch=vsearch, dry_run=dry_run)
        make_udb(str(out_sintax), str(out_sintax_udb), vsearch=vsearch, dry_run=dry_run)
        make_udb(str(out_typestrains), str(out_typestrains_udb), vsearch=vsearch, dry_run=dry_run)
        entries.update({
            "silva_udb": str(out_udb),
            "silva_sintax_udb": str(out_sintax_udb),
            "silva_typestrains_udb": str(out_typestrains_udb),
        })
    write_manifest(manifest, entries)
    entries["manifest"] = str(manifest)
    return entries
