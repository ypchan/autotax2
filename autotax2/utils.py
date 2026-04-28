from __future__ import annotations

import gzip
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple


def open_text(path: str | Path, mode: str = "rt"):
    """Open plain text or .gz files transparently."""
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding=None if "b" in mode else "utf-8")


def ensure_file(path: str | Path, label: str = "file") -> Path:
    p = Path(path)
    if not p.exists() or not p.is_file() or p.stat().st_size == 0:
        raise FileNotFoundError(f"{label} not found or empty: {p}")
    return p


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def run_cmd(cmd: List[str], dry_run: bool = False) -> None:
    """Run an external command and print it for reproducibility."""
    from .logging import info
    printable = " ".join(str(x) for x in cmd)
    info(f"RUN: {printable}")
    if dry_run:
        return
    subprocess.run(cmd, check=True)


def which_or_error(program: str) -> str:
    path = shutil.which(program)
    if path is None:
        raise RuntimeError(
            f"Required executable not found in PATH: {program}. "
            f"Install it or pass the correct --vsearch path."
        )
    return path


def read_fasta(path: str | Path) -> Iterator[Tuple[str, str]]:
    """Yield (header_without_>, sequence) from FASTA or FASTA.gz."""
    header = None
    seq_parts: List[str] = []
    with open_text(path, "rt") as handle:
        for raw in handle:
            line = raw.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            yield header, "".join(seq_parts)


def write_fasta(records: Iterable[Tuple[str, str]], path: str | Path, wrap: int = 80) -> None:
    with open(path, "w", encoding="utf-8") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            if wrap and wrap > 0:
                for i in range(0, len(seq), wrap):
                    out.write(seq[i:i + wrap] + "\n")
            else:
                out.write(seq + "\n")


def copy_or_decompress(src: str | Path, dest: str | Path) -> Path:
    src = Path(src)
    dest = Path(dest)
    if str(src).endswith(".gz"):
        with gzip.open(src, "rb") as inp, open(dest, "wb") as out:
            shutil.copyfileobj(inp, out)
    else:
        shutil.copyfile(src, dest)
    return dest


def parse_manifest(path: str | Path) -> Dict[str, str]:
    data: Dict[str, str] = {}
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        if header[:2] != ["key", "path"]:
            raise ValueError(f"Manifest must start with columns: key<TAB>path. Found: {header}")
        for line in f:
            if not line.strip():
                continue
            key, value = line.rstrip("\n").split("\t", 1)
            data[key] = value
    return data


def write_manifest(path: str | Path, entries: Dict[str, str]) -> None:
    with open(path, "w", encoding="utf-8") as out:
        out.write("key\tpath\n")
        for key, value in entries.items():
            out.write(f"{key}\t{value}\n")


def strip_fasta_suffix(name: str) -> str:
    for suffix in [".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna", ".gz"]:
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return name
