from __future__ import annotations

import csv
import re
import shutil
from pathlib import Path
from typing import Callable, Dict, List, Optional, Set, Tuple

from .logging import ProgressCounter, step
from .threads import validate_threads
from .utils import copy_or_decompress, ensure_dir, ensure_file, read_fasta, write_manifest
from .vsearch import make_udb


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
    """Parse a SILVA FASTA header.

    Expected format:
        <accession> <tax1;tax2;tax3;...>
    """

    parts = header.split(" ", 1)
    seq_id = parts[0]

    if len(parts) == 1:
        taxa: List[Optional[str]] = []
    else:
        taxa = [clean_taxon(item) for item in parts[1].split(";")]

    taxa = (taxa + [None] * len(TAX_COLS))[:len(TAX_COLS)]

    return seq_id, taxa


def make_sintax_header(seq_id: str, taxa: List[Optional[str]]) -> str:
    pairs: List[str] = []

    for prefix, taxon in zip(SINTAX_PREFIXES, taxa):
        if taxon:
            pairs.append(prefix + taxon)

    return f"{seq_id};tax={','.join(pairs)};"


def detect_metadata_columns(header: List[str]) -> Tuple[int, int]:
    lowered = [item.lower() for item in header]

    if "acc" not in lowered:
        raise ValueError("SILVA metadata must contain an 'acc' column.")

    if "flags" not in lowered:
        raise ValueError("SILVA metadata must contain a 'flags' column.")

    return lowered.index("acc"), lowered.index("flags")


def extract_typestrain_accessions(metadata_path: str | Path) -> Set[str]:
    """Extract accession IDs whose SILVA metadata flags contain [T] or [t]."""

    ids: Set[str] = set()

    with open(metadata_path, encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)

        acc_i, flags_i = detect_metadata_columns(header)

        for row in reader:
            if len(row) <= max(acc_i, flags_i):
                continue

            acc = row[acc_i].strip()
            flags = row[flags_i]

            if acc and re.search(r"\[t\]", flags, flags=re.IGNORECASE):
                ids.add(acc)

    return ids


def wrap_sequence(seq: str, width: int = 80) -> List[str]:
    return [seq[i:i + width] for i in range(0, len(seq), width)]


def write_fasta_record(handle, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")

    for line in wrap_sequence(sequence):
        handle.write(line + "\n")


def normalize_allowed_bases(allowed_bases: str) -> Set[str]:
    allowed = set(allowed_bases.upper().replace(",", "").replace(" ", ""))

    if not allowed:
        raise ValueError("--allowed-bases cannot be empty.")

    invalid = sorted(base for base in allowed if len(base) != 1 or not base.isalpha())

    if invalid:
        raise ValueError(
            "--allowed-bases must contain only base letters, got: "
            + "".join(invalid)
        )

    return allowed


def clean_silva_fasta(
    input_fasta: str | Path,
    output_fasta: str | Path,
    summary_report: str | Path,
    dropped_report: str | Path,
    *,
    allowed_bases: str = "ACGT",
    convert_u_to_t: bool = True,
    log: Callable[[str], None] = step,
) -> Dict[str, int | str]:
    """Clean SILVA FASTA before reference preparation.

    Rules:
      1. Convert U/u to T when convert_u_to_t=True.
      2. Keep only sequences whose characters are in allowed_bases.
      3. Drop the entire sequence if any other character is found.
      4. Write a cleaned FASTA plus summary and dropped-record reports.
    """

    input_fasta = Path(input_fasta)
    output_fasta = Path(output_fasta)
    summary_report = Path(summary_report)
    dropped_report = Path(dropped_report)

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    summary_report.parent.mkdir(parents=True, exist_ok=True)
    dropped_report.parent.mkdir(parents=True, exist_ok=True)

    allowed = normalize_allowed_bases(allowed_bases)

    total_seq = 0
    kept_seq = 0
    dropped_seq = 0

    total_bp = 0
    kept_bp = 0
    dropped_bp = 0

    total_u_to_t = 0
    seq_with_u_to_t = 0

    log("Cleaning SILVA FASTA sequences")
    progress = ProgressCounter("FASTA records", interval=100000)

    with output_fasta.open("w", encoding="utf-8") as out_handle, \
            dropped_report.open("w", encoding="utf-8") as dropped_handle:

        dropped_handle.write(
            "id\theader\tlength\tU_to_T_count\tinvalid_chars\n"
        )

        for header, seq in read_fasta(input_fasta):
            progress.update()

            total_seq += 1

            raw_seq = "".join(str(seq).split()).upper()
            total_bp += len(raw_seq)

            u_count = raw_seq.count("U")
            converted_seq = raw_seq.replace("U", "T") if convert_u_to_t else raw_seq

            total_u_to_t += u_count if convert_u_to_t else 0

            if convert_u_to_t and u_count > 0:
                seq_with_u_to_t += 1

            invalid_chars = sorted(set(converted_seq) - allowed)
            seq_id = header.split()[0] if header else ""

            if invalid_chars:
                dropped_seq += 1
                dropped_bp += len(raw_seq)

                dropped_handle.write(
                    f"{seq_id}\t{header}\t{len(raw_seq)}\t"
                    f"{u_count if convert_u_to_t else 0}\t"
                    f"{''.join(invalid_chars)}\n"
                )
                continue

            kept_seq += 1
            kept_bp += len(converted_seq)

            write_fasta_record(out_handle, header, converted_seq)

    progress.finish()

    with summary_report.open("w", encoding="utf-8") as summary_handle:
        summary_handle.write("metric\tvalue\n")
        summary_handle.write(f"input_fasta\t{input_fasta}\n")
        summary_handle.write(f"cleaned_fasta\t{output_fasta}\n")
        summary_handle.write(f"dropped_detail_report\t{dropped_report}\n")
        summary_handle.write(f"allowed_bases\t{''.join(sorted(allowed))}\n")
        summary_handle.write(f"convert_u_to_t\t{convert_u_to_t}\n")
        summary_handle.write(f"input_sequences\t{total_seq}\n")
        summary_handle.write(f"kept_sequences\t{kept_seq}\n")
        summary_handle.write(f"dropped_sequences\t{dropped_seq}\n")
        summary_handle.write(f"input_bases\t{total_bp}\n")
        summary_handle.write(f"kept_bases\t{kept_bp}\n")
        summary_handle.write(f"dropped_bases\t{dropped_bp}\n")
        summary_handle.write(f"sequences_with_U_converted_to_T\t{seq_with_u_to_t}\n")
        summary_handle.write(f"total_U_to_T_conversions\t{total_u_to_t}\n")

    return {
        "input_sequences": total_seq,
        "kept_sequences": kept_seq,
        "dropped_sequences": dropped_seq,
        "input_bases": total_bp,
        "kept_bases": kept_bp,
        "dropped_bases": dropped_bp,
        "sequences_with_U_converted_to_T": seq_with_u_to_t,
        "total_U_to_T_conversions": total_u_to_t,
        "cleaned_fasta": str(output_fasta),
        "summary_report": str(summary_report),
        "dropped_report": str(dropped_report),
    }


def accession_matches_typestrain(seq_id: str, typestrain_ids: Set[str]) -> bool:
    """Return True if a SILVA sequence accession matches the type-strain set."""

    if seq_id in typestrain_ids:
        return True

    # SILVA FASTA IDs may contain accession versions while metadata may not.
    # Example: AB123456.1 in FASTA, AB123456 in metadata.
    base_id = seq_id.split(".", 1)[0]

    return base_id in typestrain_ids


def prepare_silva(
    silva_fasta: str | Path,
    silva_metadata: str | Path,
    outdir: str | Path,
    *,
    make_udb_files: bool = False,
    vsearch: str | Path = "vsearch",
    threads: int = 4,
    clean_sequences: bool = True,
    allowed_bases: str = "ACGT",
    convert_u_to_t: bool = True,
    keep_decompressed: bool = True,
    log: Callable[[str], None] = step,
) -> Dict[str, str]:
    """Prepare local SILVA FASTA and metadata files for AutoTax2.

    No downloading is performed.

    Main outputs:
      silva.fasta
      silva_sintax.fasta
      silva_typestrains.fasta
      silva_taxonomy.tsv
      typestrains_accessionIDs.txt
      autotax2_ref_manifest.tsv
    """

    threads = validate_threads(threads)

    silva_fasta = ensure_file(silva_fasta, "SILVA FASTA")
    silva_metadata = ensure_file(silva_metadata, "SILVA metadata")
    outdir = ensure_dir(outdir)

    raw_fasta = outdir / "silva_raw.fasta"
    out_fasta = outdir / "silva.fasta"
    out_sintax = outdir / "silva_sintax.fasta"
    out_typestrains = outdir / "silva_typestrains.fasta"
    out_ids = outdir / "typestrains_accessionIDs.txt"
    out_tax = outdir / "silva_taxonomy.tsv"
    manifest = outdir / "autotax2_ref_manifest.tsv"

    metadata_plain = outdir / "silva.full_metadata"

    cleaning_summary = outdir / "silva_fasta_cleaning_summary.tsv"
    cleaning_dropped = outdir / "silva_fasta_cleaning_dropped.tsv"

    log("Copying or decompressing SILVA FASTA")
    copy_or_decompress(silva_fasta, raw_fasta)

    if clean_sequences:
        clean_silva_fasta(
            raw_fasta,
            out_fasta,
            cleaning_summary,
            cleaning_dropped,
            allowed_bases=allowed_bases,
            convert_u_to_t=convert_u_to_t,
            log=log,
        )
    else:
        log("Skipping FASTA cleaning")
        shutil.copyfile(raw_fasta, out_fasta)

        with cleaning_summary.open("w", encoding="utf-8") as handle:
            handle.write("metric\tvalue\n")
            handle.write("clean_sequences\tFalse\n")
            handle.write(f"input_fasta\t{raw_fasta}\n")
            handle.write(f"output_fasta\t{out_fasta}\n")

        with cleaning_dropped.open("w", encoding="utf-8") as handle:
            handle.write("id\theader\tlength\tU_to_T_count\tinvalid_chars\n")

    log("Copying or decompressing SILVA metadata")
    copy_or_decompress(silva_metadata, metadata_plain)

    log("Extracting type-strain accession IDs")
    type_ids = extract_typestrain_accessions(metadata_plain)

    with out_ids.open("w", encoding="utf-8") as handle:
        for acc in sorted(type_ids):
            handle.write(acc + "\n")

    log("Writing taxonomy, SINTAX FASTA, and type-strain FASTA")
    progress = ProgressCounter("cleaned FASTA records", interval=100000)

    with out_tax.open("w", encoding="utf-8") as tax_handle, \
            out_sintax.open("w", encoding="utf-8") as sintax_handle, \
            out_typestrains.open("w", encoding="utf-8") as typestrain_handle:

        tax_handle.write("ID\t" + "\t".join(TAX_COLS) + "\n")

        for header, seq in read_fasta(out_fasta):
            progress.update()

            seq_id, taxa = parse_silva_header(header)

            tax_handle.write(
                seq_id
                + "\t"
                + "\t".join(taxon if taxon else "" for taxon in taxa)
                + "\n"
            )

            write_fasta_record(
                sintax_handle,
                make_sintax_header(seq_id, taxa),
                seq,
            )

            if accession_matches_typestrain(seq_id, type_ids):
                write_fasta_record(typestrain_handle, header, seq)

    progress.finish()

    entries: Dict[str, str] = {
        "silva_raw_fasta": str(raw_fasta),
        "silva_fasta": str(out_fasta),
        "silva_sintax_fasta": str(out_sintax),
        "silva_typestrains_fasta": str(out_typestrains),
        "silva_metadata": str(metadata_plain),
        "typestrains_accession_ids": str(out_ids),
        "silva_taxonomy_tsv": str(out_tax),
        "silva_fasta_cleaning_summary": str(cleaning_summary),
        "silva_fasta_cleaning_dropped": str(cleaning_dropped),
    }

    if make_udb_files:
        log("Building VSEARCH UDB files")

        out_udb = outdir / "silva.udb"
        out_sintax_udb = outdir / "silva_sintax.udb"
        out_typestrains_udb = outdir / "silva_typestrains.udb"

        make_udb(out_fasta, out_udb, vsearch=vsearch, threads=threads)
        make_udb(out_sintax, out_sintax_udb, vsearch=vsearch, threads=threads)
        make_udb(out_typestrains, out_typestrains_udb, vsearch=vsearch, threads=threads)

        entries.update(
            {
                "silva_udb": str(out_udb),
                "silva_sintax_udb": str(out_sintax_udb),
                "silva_typestrains_udb": str(out_typestrains_udb),
            }
        )

    log("Writing reference manifest")
    write_manifest(manifest, entries)

    entries["manifest"] = str(manifest)

    if not keep_decompressed and raw_fasta.exists():
        log("Removing decompressed raw FASTA")
        raw_fasta.unlink()

    return entries