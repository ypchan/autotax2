from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence

from .logging import print_command, step
from .threads import validate_threads


class ExternalCommandError(RuntimeError):
    """Raised when an external command fails."""


@dataclass(frozen=True)
class CommandResult:
    command: Sequence[str | Path]
    returncode: int
    stdout: str = ""
    stderr: str = ""


EXECUTABLE_NAMES: Dict[str, str] = {
    "vsearch": "vsearch",
    "blastn": "blastn",
    "makeblastdb": "makeblastdb",
    "mafft": "mafft",
    "meme": "meme",
    "streme": "streme",
    "fimo": "fimo",
    "rnafold": "RNAfold",
    "rnaalifold": "RNAalifold",
    "barrnap": "barrnap",
}


def command_to_string(command: Sequence[str | Path]) -> str:
    return " ".join(shlex.quote(str(part)) for part in command)


def resolve_executable(executable: str | Path) -> str:
    """Resolve an executable from an absolute path or PATH."""

    text = str(executable).strip()

    if not text:
        raise FileNotFoundError("Executable path is empty.")

    path = Path(text).expanduser()

    if path.exists() and os.access(path, os.X_OK):
        return str(path)

    detected = shutil.which(text)

    if detected:
        return detected

    raise FileNotFoundError(f"Executable not found: {executable}")


def resolve_named_executable(name: str, configured: str | Path | None = None) -> str:
    """Resolve a known tool name, using a configured path when provided."""

    if configured:
        return resolve_executable(configured)

    executable_name = EXECUTABLE_NAMES.get(name.lower(), name)
    return resolve_executable(executable_name)


def ensure_parent(path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def run_command(
    command: Sequence[str | Path],
    *,
    cwd: str | Path | None = None,
    env: Optional[Mapping[str, str]] = None,
    stdout_path: str | Path | None = None,
    stderr_path: str | Path | None = None,
    input_text: str | None = None,
    capture_output: bool = False,
) -> CommandResult:
    """Run an external command.

    Parameters
    ----------
    command
        Command and arguments.
    cwd
        Working directory.
    env
        Extra environment variables.
    stdout_path
        Optional file path for stdout redirection.
    stderr_path
        Optional file path for stderr redirection.
    input_text
        Optional text sent to stdin.
    capture_output
        Capture stdout/stderr in memory. Do not use together with stdout_path
        or stderr_path.

    Returns
    -------
    CommandResult
        Command, return code, and captured output when requested.
    """

    if not command:
        raise ValueError("Command is empty.")

    if capture_output and (stdout_path is not None or stderr_path is not None):
        raise ValueError("capture_output cannot be combined with stdout_path or stderr_path.")

    print_command(command)

    merged_env = os.environ.copy()

    if env:
        merged_env.update({str(key): str(value) for key, value in env.items()})

    stdout_handle = None
    stderr_handle = None

    try:
        stdout_target = subprocess.PIPE if capture_output else None
        stderr_target = subprocess.PIPE if capture_output else None

        if stdout_path is not None:
            stdout_file = ensure_parent(stdout_path)
            stdout_handle = stdout_file.open("w", encoding="utf-8")
            stdout_target = stdout_handle

        if stderr_path is not None:
            stderr_file = ensure_parent(stderr_path)
            stderr_handle = stderr_file.open("w", encoding="utf-8")
            stderr_target = stderr_handle

        completed = subprocess.run(
            [str(part) for part in command],
            cwd=str(cwd) if cwd is not None else None,
            env=merged_env,
            input=input_text,
            text=True,
            stdout=stdout_target,
            stderr=stderr_target,
            check=False,
        )

    finally:
        if stdout_handle is not None:
            stdout_handle.close()

        if stderr_handle is not None:
            stderr_handle.close()

    stdout = completed.stdout or ""
    stderr = completed.stderr or ""

    if completed.returncode != 0:
        message = [
            f"External command failed with exit code {completed.returncode}:",
            command_to_string(command),
        ]

        if stdout:
            message.append("")
            message.append("stdout:")
            message.append(stdout[-4000:])

        if stderr:
            message.append("")
            message.append("stderr:")
            message.append(stderr[-4000:])

        raise ExternalCommandError("\n".join(message))

    return CommandResult(
        command=command,
        returncode=completed.returncode,
        stdout=stdout,
        stderr=stderr,
    )


def capture_command(
    command: Sequence[str | Path],
    *,
    cwd: str | Path | None = None,
    env: Optional[Mapping[str, str]] = None,
) -> CommandResult:
    return run_command(
        command,
        cwd=cwd,
        env=env,
        capture_output=True,
    )


def executable_version(executable: str | Path) -> str:
    """Return a short version string for an executable."""

    resolved = resolve_executable(executable)

    for args in (["--version"], ["-version"], ["-h"]):
        try:
            result = capture_command([resolved, *args])
        except Exception:
            continue

        output = (result.stdout or result.stderr or "").strip()

        if output:
            return output.splitlines()[0]

    return ""


def make_blast_db(
    input_fasta: str | Path,
    db_prefix: str | Path,
    *,
    makeblastdb: str | Path,
    title: str | None = None,
    dbtype: str = "nucl",
) -> str:
    """Build a nucleotide BLAST database."""

    input_fasta = Path(input_fasta)
    db_prefix = Path(db_prefix)
    db_prefix.parent.mkdir(parents=True, exist_ok=True)

    step(f"Building BLAST database: {db_prefix}")

    command: list[str | Path] = [
        makeblastdb,
        "-in",
        input_fasta,
        "-dbtype",
        dbtype,
        "-out",
        db_prefix,
    ]

    if title:
        command.extend(["-title", title])

    run_command(command)

    return str(db_prefix)


def blastn_local_hsps(
    query_fasta: str | Path,
    db: str | Path,
    output_tsv: str | Path,
    *,
    blastn: str | Path,
    threads: int,
    perc_identity: float = 90.0,
    max_hsps: int = 20,
    max_target_seqs: int = 500,
    task: str = "blastn",
    dust: str = "no",
    evalue: float = 1e-20,
) -> str:
    """Run BLASTN for local multi-HSP intron discovery.

    Output fields:
      qseqid sseqid pident length mismatch gapopen qstart qend sstart send
      evalue bitscore qlen slen qseq sseq
    """

    threads = validate_threads(threads)

    if not 0 <= perc_identity <= 100:
        raise ValueError("perc_identity must be between 0 and 100.")

    if max_hsps < 1:
        raise ValueError("max_hsps must be at least 1.")

    if max_target_seqs < 1:
        raise ValueError("max_target_seqs must be at least 1.")

    output_tsv = ensure_parent(output_tsv)

    outfmt = (
        "6 qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore qlen slen qseq sseq"
    )

    step("Running BLASTN local-HSP search")

    command = [
        blastn,
        "-query",
        query_fasta,
        "-db",
        db,
        "-task",
        task,
        "-perc_identity",
        str(perc_identity),
        "-max_hsps",
        str(max_hsps),
        "-max_target_seqs",
        str(max_target_seqs),
        "-evalue",
        str(evalue),
        "-dust",
        dust,
        "-num_threads",
        str(threads),
        "-outfmt",
        outfmt,
        "-out",
        output_tsv,
    ]

    run_command(command)
    return str(output_tsv)


def mafft_align(
    input_fasta: str | Path,
    output_fasta: str | Path,
    *,
    mafft: str | Path,
    threads: int,
    auto: bool = True,
) -> str:
    """Run MAFFT and write aligned FASTA."""

    threads = validate_threads(threads)
    output_fasta = ensure_parent(output_fasta)

    step(f"Running MAFFT alignment: {output_fasta}")

    command: list[str | Path] = [
        mafft,
        "--thread",
        str(threads),
    ]

    if auto:
        command.append("--auto")

    command.append(input_fasta)

    run_command(command, stdout_path=output_fasta)

    return str(output_fasta)


def meme_discover(
    input_fasta: str | Path,
    output_dir: str | Path,
    *,
    meme: str | Path,
    threads: int,
    nmotifs: int = 10,
    minw: int = 6,
    maxw: int = 50,
    dna: bool = True,
) -> str:
    """Run MEME motif discovery."""

    threads = validate_threads(threads)

    if nmotifs < 1:
        raise ValueError("nmotifs must be at least 1.")

    if minw < 1 or maxw < minw:
        raise ValueError("motif widths must satisfy 1 <= minw <= maxw.")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    step(f"Running MEME motif discovery: {output_dir}")

    command: list[str | Path] = [
        meme,
        input_fasta,
        "-oc",
        output_dir,
        "-nmotifs",
        str(nmotifs),
        "-minw",
        str(minw),
        "-maxw",
        str(maxw),
        "-p",
        str(threads),
    ]

    if dna:
        command.append("-dna")

    run_command(command)

    return str(output_dir)


def streme_discover(
    input_fasta: str | Path,
    output_dir: str | Path,
    *,
    streme: str | Path,
    minw: int = 6,
    maxw: int = 30,
    dna: bool = True,
) -> str:
    """Run STREME motif discovery."""

    if minw < 1 or maxw < minw:
        raise ValueError("motif widths must satisfy 1 <= minw <= maxw.")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    step(f"Running STREME motif discovery: {output_dir}")

    command: list[str | Path] = [
        streme,
        "--p",
        input_fasta,
        "--oc",
        output_dir,
        "--minw",
        str(minw),
        "--maxw",
        str(maxw),
    ]

    if dna:
        command.append("--dna")

    run_command(command)

    return str(output_dir)


def fimo_scan(
    motif_file: str | Path,
    sequence_fasta: str | Path,
    output_dir: str | Path,
    *,
    fimo: str | Path,
    thresh: float = 1e-4,
) -> str:
    """Run FIMO to scan motifs against sequences."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    step(f"Running FIMO motif scanning: {output_dir}")

    command = [
        fimo,
        "--oc",
        output_dir,
        "--thresh",
        str(thresh),
        motif_file,
        sequence_fasta,
    ]

    run_command(command)

    return str(output_dir)


def rnafold_predict(
    input_fasta: str | Path,
    output_txt: str | Path,
    *,
    rnafold: str | Path,
) -> str:
    """Run RNAfold for per-sequence secondary-structure prediction."""

    output_txt = ensure_parent(output_txt)

    step(f"Running RNAfold: {output_txt}")

    with Path(input_fasta).open("r", encoding="utf-8", errors="replace") as input_handle:
        input_text = input_handle.read()

    result = run_command(
        [rnafold, "--noPS"],
        input_text=input_text,
        capture_output=True,
    )

    output_txt.write_text(result.stdout, encoding="utf-8")

    return str(output_txt)


def rnaalifold_predict(
    alignment_fasta: str | Path,
    output_txt: str | Path,
    *,
    rnaalifold: str | Path,
) -> str:
    """Run RNAalifold for consensus secondary structure from an alignment."""

    output_txt = ensure_parent(output_txt)

    step(f"Running RNAalifold: {output_txt}")

    with Path(alignment_fasta).open("r", encoding="utf-8", errors="replace") as input_handle:
        input_text = input_handle.read()

    result = run_command(
        [rnaalifold, "--noPS"],
        input_text=input_text,
        capture_output=True,
    )

    output_txt.write_text(result.stdout, encoding="utf-8")

    return str(output_txt)


def barrnap_annotate(
    input_fasta: str | Path,
    output_gff: str | Path,
    *,
    barrnap: str | Path,
    threads: int,
    kingdom: str = "bac",
) -> str:
    """Run barrnap rRNA annotation."""

    threads = validate_threads(threads)
    output_gff = ensure_parent(output_gff)

    step(f"Running barrnap annotation: {output_gff}")

    command = [
        barrnap,
        "--kingdom",
        kingdom,
        "--threads",
        str(threads),
        input_fasta,
    ]

    run_command(command, stdout_path=output_gff)

    return str(output_gff)