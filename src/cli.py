from __future__ import annotations

import inspect
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional

import typer

from .assign import assign_or_create
from .backbone import insert_backbone
from .config import (
    config_path_message,
    get_reference_path,
    get_ref_manifest,
    get_software,
    get_vsearch,
    init_config,
    iter_config_rows,
    parse_manifest,
    set_config_value,
    update_reference_from_manifest,
)
from .dependencies import check_all, check_from_config, print_check_report
from .intron import detect_introns
from .intron_analysis import analyze_introns
from .intron_coordinates import map_intron_coordinates
from .intron_reference import (
    build_default_intron_references,
    build_species_representative_reference,
)
from .orientation import orient_with_sina
from .logging import print_outputs, step, success
from .overlap import overlap_backbone
from .prepare import prepare_silva
from .provenance import annotate_assignments_with_source, provenance_from_uc_levels
from .summarize import summarize
from .threads import DEFAULT_THREADS
from .utils import ensure_dir
from .vsearch import cluster as run_cluster
from .vsearch import dereplicate, sintax, sort_by_size, usearch_global


CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}


app = typer.Typer(
    name="autotax2",
    help="AutoTax2 workflow.",
    no_args_is_help=True,
    add_completion=False,
    context_settings=CONTEXT_SETTINGS,
    rich_markup_mode=None,
)

config_app = typer.Typer(
    name="config",
    help="Configure AutoTax2 software and reference paths.",
    no_args_is_help=True,
    add_completion=False,
    context_settings=CONTEXT_SETTINGS,
    rich_markup_mode=None,
)


def require_options(**values: object) -> None:
    missing = [
        key.replace("_", "-")
        for key, value in values.items()
        if value in (None, "", [])
    ]

    if missing:
        raise typer.BadParameter(
            "Missing required option(s): "
            + ", ".join(f"--{item}" for item in missing)
        )


def parse_id_list(value: str) -> List[float]:
    try:
        ids = [float(item.strip()) for item in value.split(",") if item.strip()]
    except ValueError as exc:
        raise typer.BadParameter(
            "--ids must be a comma-separated list of numbers, for example: 0.99,0.97"
        ) from exc

    if not ids:
        raise typer.BadParameter("--ids cannot be empty.")

    for identity in ids:
        if not 0 < identity <= 1:
            raise typer.BadParameter("--ids values must be between 0 and 1.")

    return ids


def call_supported(func: Any, **kwargs: Any) -> Any:
    signature = inspect.signature(func)
    accepted = {
        key: value
        for key, value in kwargs.items()
        if key in signature.parameters
    }
    return func(**accepted)


def load_manifest_or_config(ref_manifest: Optional[Path]) -> tuple[Path, Dict[str, str]]:
    manifest_path = ref_manifest or get_ref_manifest(required=True)

    if manifest_path is None:
        raise typer.BadParameter(
            "Reference manifest is not configured. "
            "Run: autotax2 config set ref-manifest refdatabases/autotax2_ref_manifest.tsv"
        )

    return manifest_path, parse_manifest(manifest_path)


def resolve_reference_for_detect_intron(db: Optional[Path]) -> Path:
    if db is not None:
        return db

    keys = [
        "species-rep1-clean-blastdb",
        "species-rep1-blastdb",
        "species-rep1-clean-fasta",
        "species-rep1-fasta",
    ]

    for key in keys:
        configured = get_reference_path(key, required=False)
        if configured is not None:
            return configured

    raise typer.BadParameter(
        "Missing intron reference. Use --db or build/configure one with:\n"
        "  autotax2 build-intron-ref\n"
        "or:\n"
        "  autotax2 config set species-rep1-clean-blastdb <blastdb-prefix>"
    )


def resolve_ref_file(
    explicit_path: Optional[Path],
    manifest: Dict[str, str],
    manifest_key: str,
    config_key: str,
    label: str,
) -> Path:
    if explicit_path is not None:
        return explicit_path

    if manifest.get(manifest_key):
        return Path(manifest[manifest_key])

    configured = get_reference_path(config_key, required=False)
    if configured is not None:
        return configured

    raise typer.BadParameter(
        f"Missing {label}. Provide the command-line option, use --ref-manifest, "
        f"or configure it with: autotax2 config set {config_key} <path>"
    )



def optional_software(key: str) -> Optional[str]:
    """Return an optional executable only when it is actually available."""

    try:
        value = get_software(key, required=False)
    except Exception:
        return None

    if not value:
        return None

    path = Path(value).expanduser()

    if path.exists():
        return str(path)

    detected = shutil.which(str(value))

    if detected:
        return detected

    return None


def blastdb_prefix_exists(prefix: str | Path) -> bool:
    prefix = Path(prefix).expanduser()
    suffixes = (
        ".nhr", ".nin", ".nsq", ".ndb", ".not", ".ntf", ".nto",
        ".00.nhr", ".00.nin", ".00.nsq", ".00.ndb", ".00.not", ".00.ntf", ".00.nto",
    )
    return any(Path(str(prefix) + suffix).exists() for suffix in suffixes)


def resolve_sina_reference(sina_ref: Optional[Path]) -> Path:
    if sina_ref is not None:
        return sina_ref

    configured = get_reference_path("sina-ref", required=False)
    if configured is not None:
        return configured

    raise typer.BadParameter(
        "Missing SINA reference database. Provide --sina-ref or configure it with:\n"
        "  autotax2 config set sina-ref <path-to-sina-reference>"
    )


def run_sina_orientation(
    input_fasta: Path,
    outdir: Path,
    *,
    sina_ref: Optional[Path],
    threads: int,
    score_margin: float,
    min_score: float,
    kmer_size: int,
) -> Dict[str, str]:
    sina = get_software("sina", required=True)
    ref = resolve_sina_reference(sina_ref)

    return orient_with_sina(
        input_fasta=input_fasta,
        outdir=outdir,
        sina=sina,
        sina_ref=ref,
        threads=threads,
        score_margin=score_margin,
        min_score=min_score,
        kmer_size=kmer_size,
    )


@config_app.command("init", help="Create the AutoTax2 configuration file.")
def config_init(
    overwrite: bool = typer.Option(
        False,
        "--overwrite",
        help="Overwrite the existing configuration file.",
    ),
) -> None:
    """
    Create the AutoTax2 configuration file.

    Default path:
      ~/.config/autotax2/config.ini

    The command tries to detect configured software from PATH:
      vsearch, blastn, makeblastdb, sina, mafft, meme, streme, fimo,
      RNAfold, RNAalifold, barrnap
    """

    path = init_config(overwrite=overwrite)
    success(f"Configuration file: {path}")


@config_app.command("set", help="Set one persistent configuration value.")
def config_set(
    key: str = typer.Argument(
        ...,
        help=(
            "Configuration key, for example: vsearch, blastn, makeblastdb, "
            "sina, mafft, meme, streme, fimo, RNAfold, RNAalifold, barrnap, "
            "ref-manifest, silva-fasta, species-rep1-blastdb, coordinate-ref, sina-ref, threads."
        ),
    ),
    value: str = typer.Argument(
        ...,
        help="Value to store.",
    ),
) -> None:
    """
    Set a persistent software path, reference path, or runtime value.

    Examples:
      autotax2 config set vsearch /home/user/bin/vsearch
      autotax2 config set blastn /home/user/bin/blastn
      autotax2 config set makeblastdb /home/user/bin/makeblastdb
      autotax2 config set sina /home/user/bin/sina
      autotax2 config set sina-ref /path/to/sina_reference.arb
      autotax2 config set mafft /home/user/bin/mafft
      autotax2 config set RNAfold /home/user/bin/RNAfold
      autotax2 config set barrnap /home/user/bin/barrnap
      autotax2 config set ref-manifest refdatabases/autotax2_ref_manifest.tsv
      autotax2 config set species-rep1-blastdb refdatabases/intron_ref/species_rep1
      autotax2 config set coordinate-ref refdatabases/coordinate_ref.fasta
      autotax2 config set threads 4
    """

    try:
        path = set_config_value(key, value)
    except (KeyError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from exc

    success(f"Updated {key} in {path}")


@config_app.command("show", help="Show the current AutoTax2 configuration.")
def config_show() -> None:
    """
    Print the current configuration file and all configured values.
    """

    typer.echo(f"Config file: {config_path_message()}")
    typer.echo()

    current_section = None

    for section_name, key, value in iter_config_rows():
        if section_name != current_section:
            if current_section is not None:
                typer.echo()
            typer.echo(section_name)
            current_section = section_name

        typer.echo(f"  {key:<30} {value}")


@config_app.command("check", help="Check configured software and reference paths.")
def config_check() -> None:
    """
    Check the current AutoTax2 configuration file.

    Required software:
      vsearch
      blastn
      makeblastdb

    Optional software:
      sina

    Optional intron-analysis software:
      mafft
      meme
      streme
      fimo
      RNAfold
      RNAalifold
      barrnap
    """

    ok = print_check_report(check_from_config(), title="AutoTax2 config check")

    if not ok:
        raise typer.Exit(1)


app.add_typer(config_app, name="config")


@app.command("check", help="Check software, Python packages, and reference data.")
def check(
    ref_manifest: Optional[Path] = typer.Option(
        None,
        "--ref-manifest",
        help="Additional reference manifest to validate. Default: configured manifest.",
    ),
    source_map: Optional[Path] = typer.Option(
        None,
        "--source-map",
        help="Optional source map file to validate.",
    ),
) -> None:
    """
    Check the AutoTax2 runtime environment.

    This command checks:
      1. Config file
      2. Required software: vsearch, blastn, makeblastdb
      3. Optional software: mafft, meme, streme, fimo, RNAfold, RNAalifold, barrnap
      4. Python packages
      5. Reference data paths
      6. Optional input files such as source maps
    """

    ok = check_all(
        ref_manifest=ref_manifest,
        source_map=source_map,
        use_config=True,
    )

    if not ok:
        raise typer.Exit(1)


@app.command("prepare-silva", help="Prepare local SILVA FASTA and metadata files.")
def prepare_silva_cmd(
    silva_fasta: Optional[Path] = typer.Option(
        None,
        "--silva-fasta",
        help="Local SILVA FASTA file, optionally .gz.",
    ),
    silva_metadata: Optional[Path] = typer.Option(
        None,
        "--silva-metadata",
        help="Local SILVA full_metadata file, optionally .gz.",
    ),
    out: Path = typer.Option(
        Path("refdatabases"),
        "--out",
        help="Output reference-data directory.",
    ),
    make_udb: bool = typer.Option(
        False,
        "--make-udb",
        help="Build VSEARCH UDB files for faster downstream search.",
    ),
    clean: bool = typer.Option(
        True,
        "--clean/--no-clean",
        help="Clean FASTA before processing.",
    ),
    allowed_bases: str = typer.Option(
        "ACGT",
        "--allowed-bases",
        help="Allowed bases after optional U-to-T conversion. Default drops N and ambiguity codes.",
    ),
    u_to_t: bool = typer.Option(
        True,
        "--u-to-t/--no-u-to-t",
        help="Convert U/u to T before filtering.",
    ),
    threads: int = typer.Option(
        DEFAULT_THREADS,
        "--threads",
        min=1,
        help="Number of threads for external tools. Python parsing is streaming single-process.",
    ),
) -> None:
    """
    Prepare SILVA reference files for AutoTax2.

    Data source:
      https://www.arb-silva.de/current-release/Exports

    Workflow:
      1. Copy or decompress SILVA FASTA and metadata
      2. Clean FASTA sequences
      3. Convert SILVA taxonomy into SINTAX format
      4. Extract type-strain records
      5. Optionally build VSEARCH UDB databases
      6. Write autotax2_ref_manifest.tsv
      7. Update the persistent AutoTax2 configuration
    """

    require_options(silva_fasta=silva_fasta, silva_metadata=silva_metadata)

    vsearch = get_vsearch(required=False)

    step("Preparing SILVA reference data")
    entries = prepare_silva(
        silva_fasta=silva_fasta,
        silva_metadata=silva_metadata,
        outdir=out,
        keep_decompressed=True,
        make_udb_files=make_udb,
        vsearch=vsearch,
        threads=threads,
        clean_sequences=clean,
        allowed_bases=allowed_bases,
        convert_u_to_t=u_to_t,
        log=step,
    )

    manifest = entries.get("manifest")

    if manifest:
        update_reference_from_manifest(manifest)
        step(f"Configuration updated from manifest: {manifest}")

    print_outputs(entries)
    success(f"SILVA preparation finished: {out}")


@app.command("orient", help="Orient near-full-length or full-length 16S sequences with SINA.")
def orient_cmd(
    input: Optional[Path] = typer.Option(
        None,
        "--input",
        help="Input FASTA to orient.",
    ),
    out: Path = typer.Option(
        Path("oriented"),
        "--out",
        help="Output directory for oriented FASTA and orientation report.",
    ),
    sina_ref: Optional[Path] = typer.Option(
        None,
        "--sina-ref",
        help="SINA reference database. Default: configured sina_ref.",
    ),
    score_margin: float = typer.Option(
        0.05,
        "--score-margin",
        min=0.0,
        help="Minimum k-mer score difference required to call one orientation over the other.",
    ),
    min_score: float = typer.Option(
        0.20,
        "--min-score",
        min=0.0,
        help="Minimum k-mer similarity score required for a non-ambiguous SINA orientation decision.",
    ),
    kmer_size: int = typer.Option(
        15,
        "--kmer-size",
        min=3,
        help="K-mer size used to compare SINA-oriented output with original and reverse-complement sequences.",
    ),
    threads: int = typer.Option(
        DEFAULT_THREADS,
        "--threads",
        min=1,
        help="Number of SINA threads.",
    ),
) -> None:
    """
    Orient 16S sequences with SINA.

    SINA is used only to determine orientation. AutoTax2 writes ungapped
    original sequences to oriented_input.fasta, reverse-complementing only
    sequences with clear SINA evidence.

    Output:
      oriented_input.fasta
      orientation_report.tsv
      reverse_complemented_ids.txt
      ambiguous_orientation_ids.txt
      sina_oriented_aligned.fasta
    """

    require_options(input=input)

    step("Orienting input sequences with SINA")

    outputs = run_sina_orientation(
        input_fasta=input,
        outdir=out,
        sina_ref=sina_ref,
        threads=threads,
        score_margin=score_margin,
        min_score=min_score,
        kmer_size=kmer_size,
    )

    print_outputs(outputs)
    success(f"Orientation finished: {out}")


@app.command("build-intron-ref", help="Build species-level references for intron detection.")
def build_intron_ref_cmd(
    silva_fasta: Optional[Path] = typer.Option(
        None,
        "--silva-fasta",
        help="Cleaned SILVA FASTA. Default: configured silva_fasta.",
    ),
    silva_taxonomy: Optional[Path] = typer.Option(
        None,
        "--silva-taxonomy",
        help="SILVA taxonomy TSV from prepare-silva. Default: configured silva_taxonomy_tsv.",
    ),
    out: Path = typer.Option(
        Path("refdatabases/intron_ref"),
        "--out",
        help="Output directory for intron-detection references.",
    ),
    per_species: Optional[int] = typer.Option(
        None,
        "--per-species",
        min=1,
        help=(
            "Build a raw representative set with this many longest sequences "
            "per species. If omitted, AutoTax2 builds the production reference: "
            "species_rep10_raw plus audited species_rep1_clean."
        ),
    ),
    build_blastdb: bool = typer.Option(
        True,
        "--blastdb/--no-blastdb",
        help="Build BLAST nucleotide databases with makeblastdb.",
    ),
    audit_reference_introns: bool = typer.Option(
        True,
        "--audit-reference-introns/--no-audit-reference-introns",
        help=(
            "Audit raw species representatives for putative reference introns "
            "before building species_rep1_clean."
        ),
    ),
    run_barrnap: bool = typer.Option(
        False,
        "--run-barrnap/--no-run-barrnap",
        help="Run barrnap on raw reference representatives and use the result in clean-reference selection.",
    ),
    raw_per_species: int = typer.Option(
        10,
        "--raw-per-species",
        min=1,
        help="Number of longest raw representatives per species used for reference audit.",
    ),
    hsp_id: float = typer.Option(
        90.0,
        "--hsp-id",
        min=0.0,
        max=100.0,
        help="Minimum identity percentage for local HSPs in the reference-intron audit.",
    ),
    min_flank_len: int = typer.Option(
        100,
        "--min-flank-len",
        min=1,
        help="Minimum HSP length on both sides of a candidate reference intron.",
    ),
    min_intron_len: int = typer.Option(
        20,
        "--min-intron-len",
        min=1,
        help="Minimum candidate reference-intron length.",
    ),
    max_ref_gap: int = typer.Option(
        10,
        "--max-ref-gap",
        min=0,
        help="Maximum allowed gap between the two HSPs on the reference coordinate.",
    ),
    ref_overlap_tolerance: int = typer.Option(
        5,
        "--ref-overlap-tolerance",
        min=0,
        help="Allowed small overlap between the two HSPs on the reference coordinate.",
    ),
    max_hsps: int = typer.Option(
        20,
        "--max-hsps",
        min=1,
        help="Maximum HSPs retained per BLAST target during reference audit.",
    ),
    max_targets: int = typer.Option(
        500,
        "--max-targets",
        min=1,
        help="Maximum BLAST targets retained per raw reference sequence during reference audit.",
    ),
    min_reference_support_refs: int = typer.Option(
        2,
        "--min-reference-support-refs",
        min=1,
        help="Minimum number of supporting references required to flag a raw reference as putatively intron-containing.",
    ),
    threads: int = typer.Option(
        DEFAULT_THREADS,
        "--threads",
        min=1,
        help="Number of threads for BLAST and barrnap.",
    ),
) -> None:
    """
    Build intron-detection references from prepared SILVA files.

    Default production mode:
      1. Build species_rep10_raw from SILVA.
      2. Audit species_rep10_raw for putative reference introns.
      3. Optionally audit raw references with barrnap.
      4. Select species_rep1_clean, avoiding flagged references when possible.
      5. Build species_rep1_clean BLAST DB.
      6. Update AutoTax2 config with species_rep1_clean and species_rep10_raw paths.

    Why reference audit matters:
      SILVA references may themselves contain intron-like insertions. If those
      references are used as clean anchors, query introns may be missed or the
      evidence may be diluted. AutoTax2 therefore audits raw reference sequences
      before selecting the default species_rep1_clean detection database.
    """

    sf = silva_fasta or get_reference_path("silva-fasta", required=True)
    st = silva_taxonomy or get_reference_path("silva-taxonomy-tsv", required=True)

    makeblastdb = None
    blastn = None
    barrnap = None

    if build_blastdb:
        makeblastdb = get_software("makeblastdb", required=True)

    if audit_reference_introns:
        if not build_blastdb:
            raise typer.BadParameter(
                "--audit-reference-introns requires --blastdb. "
                "Use --no-audit-reference-introns when running with --no-blastdb."
            )
        blastn = get_software("blastn", required=True)

    if run_barrnap:
        barrnap = get_software("barrnap", required=True)

    step("Building intron-detection reference")

    if per_species is not None:
        step(f"Building raw species_rep{per_species} reference")
        entries = build_species_representative_reference(
            silva_fasta=sf,
            silva_taxonomy=st,
            outdir=out,
            per_species=per_species,
            makeblastdb=makeblastdb or "makeblastdb",
            threads=threads,
            build_blastdb=build_blastdb,
        )

        if entries.get("manifest"):
            update_reference_from_manifest(entries["manifest"])
            step(f"Configuration updated from manifest: {entries['manifest']}")

        print_outputs(entries)
        success(f"Raw intron reference built: {out}")
        return

    entries = build_default_intron_references(
        silva_fasta=sf,
        silva_taxonomy=st,
        outdir=out,
        makeblastdb=makeblastdb or "makeblastdb",
        blastn=blastn or "blastn",
        barrnap=barrnap,
        threads=threads,
        build_blastdb=build_blastdb,
        audit_reference_introns=audit_reference_introns,
        run_barrnap=run_barrnap,
        raw_per_species=raw_per_species,
        hsp_identity=hsp_id,
        min_flank_len=min_flank_len,
        min_intron_len=min_intron_len,
        max_ref_gap=max_ref_gap,
        ref_overlap_tolerance=ref_overlap_tolerance,
        max_hsps=max_hsps,
        max_target_seqs=max_targets,
        min_reference_support_refs=min_reference_support_refs,
    )

    if entries.get("manifest"):
        update_reference_from_manifest(entries["manifest"])
        step(f"Configuration updated from manifest: {entries['manifest']}")

    print_outputs(entries)
    success(f"Audited intron reference built: {out}")

@app.command("detect-intron", help="Detect intron-like insertions with local broken-HSP evidence.")
def detect_intron_cmd(
    input: Optional[Path] = typer.Option(
        None,
        "--input",
        help="Input FASTA to screen for intron-like insertions.",
    ),
    db: Optional[Path] = typer.Option(
        None,
        "--db",
        help=(
            "BLAST database prefix or FASTA reference. "
            "Default: configured species_rep1_clean BLAST DB or FASTA."
        ),
    ),
    species_taxonomy: Optional[Path] = typer.Option(
        None,
        "--species-taxonomy",
        help=(
            "Taxonomy table for species representatives. "
            "Default: configured species_rep1_clean taxonomy or species_rep1 taxonomy."
        ),
    ),
    reference_blacklist: Optional[Path] = typer.Option(
        None,
        "--reference-blacklist",
        help=(
            "Reference IDs flagged as putatively intron-containing during "
            "build-intron-ref. Default: configured reference_intron_blacklist."
        ),
    ),
    out: Optional[Path] = typer.Option(
        None,
        "--out",
        help="Output directory.",
    ),
    source_label: Optional[str] = typer.Option(
        None,
        "--source-label",
        help="Dataset label written to output reports. It does not affect the algorithm.",
    ),

    orient: bool = typer.Option(
        False,
        "--orient/--no-orient",
        help=(
            "Run SINA orientation normalization before intron detection. "
            "Default: skip orientation normalization."
        ),
    ),
    orientation_method: str = typer.Option(
        "sina",
        "--orientation-method",
        help="Orientation method. Currently supported: sina. Used only with --orient.",
    ),
    sina_ref: Optional[Path] = typer.Option(
        None,
        "--sina-ref",
        help="SINA reference database used for --orient. Default: configured sina_ref.",
    ),
    orientation_score_margin: float = typer.Option(
        0.05,
        "--orientation-score-margin",
        min=0.0,
        help="Minimum SINA k-mer score margin required to call one orientation over the other.",
    ),
    orientation_min_score: float = typer.Option(
        0.20,
        "--orientation-min-score",
        min=0.0,
        help="Minimum SINA k-mer similarity score for a non-ambiguous orientation decision.",
    ),
    orientation_kmer_size: int = typer.Option(
        15,
        "--orientation-kmer-size",
        min=3,
        help="K-mer size used to compare SINA-oriented output with original and reverse-complement sequences.",
    ),

    hsp_id: float = typer.Option(
        90.0,
        "--hsp-id",
        min=0.0,
        max=100.0,
        help="Minimum identity percentage for each local HSP flank.",
    ),
    min_intron_len: int = typer.Option(
        20,
        "--min-intron-len",
        min=1,
        help="Minimum candidate intron-like insertion length in the query.",
    ),
    min_flank_len: int = typer.Option(
        100,
        "--min-flank-len",
        min=1,
        help="Minimum length for both left and right HSP flanks.",
    ),
    max_ref_gap: int = typer.Option(
        10,
        "--max-ref-gap",
        min=0,
        help=(
            "Maximum allowed gap between the two HSPs on the reference. "
            "Small gaps are allowed because local HSP boundaries can shift slightly."
        ),
    ),
    ref_overlap_tolerance: int = typer.Option(
        5,
        "--ref-overlap-tolerance",
        min=0,
        help="Allowed small overlap between the two HSPs on the reference coordinate.",
    ),
    max_hsps: int = typer.Option(
        20,
        "--max-hsps",
        min=1,
        help="Maximum HSPs retained per BLAST target.",
    ),
    max_targets: int = typer.Option(
        500,
        "--max-targets",
        min=1,
        help="Maximum BLAST targets retained per query.",
    ),
    strand: str = typer.Option(
        "both",
        "--strand",
        help="BLAST strand mode: both, plus, or minus.",
    ),

    confirm_id: float = typer.Option(
        0.987,
        "--confirm-id",
        min=0.0,
        max=1.0,
        help=(
            "Minimum identity required after removing a candidate intron-like "
            "insertion and re-aligning the cleaned sequence."
        ),
    ),
    confirm_qcov: float = typer.Option(
        0.80,
        "--confirm-qcov",
        min=0.0,
        max=1.0,
        help=(
            "Minimum query coverage required after removing a candidate "
            "intron-like insertion and re-aligning the cleaned sequence."
        ),
    ),

    min_support_refs: int = typer.Option(
        3,
        "--min-support-refs",
        min=1,
        help="Minimum clean reference support for high-confidence classification.",
    ),
    min_support_species: int = typer.Option(
        2,
        "--min-support-species",
        min=1,
        help="Minimum clean species support for high-confidence classification.",
    ),
    min_support_genera: int = typer.Option(
        1,
        "--min-support-genera",
        min=1,
        help="Minimum clean genus support for high-confidence classification.",
    ),
    coordinate_tolerance: int = typer.Option(
        10,
        "--coordinate-tolerance",
        min=0,
        help="Tolerance in query coordinates for clustering nearby breakpoint candidates.",
    ),

    threads: int = typer.Option(
        DEFAULT_THREADS,
        "--threads",
        min=1,
        help="Number of BLAST or SINA threads.",
    ),
) -> None:
    """
    Detect intron-like insertions using local broken-HSP evidence.

    Orientation behavior:
      By default, this command does not change input sequence orientation.
      Use --orient to run SINA-based orientation normalization first. The
      oriented FASTA is then used for intron detection.

    Output:
      preprocess/oriented_input.fasta, when --orient is used
      preprocess/orientation_report.tsv, when --orient is used
      blast/local_hsps.tsv
      blast/confirmation_candidates.fasta
      blast/confirmation_hits.tsv
      candidates/fracture_candidates.tsv
      candidates/fracture_support_by_rank.tsv
      candidates/failed_confirmation.tsv
      confirmed/analysis_sequences.fasta
      confirmed/sequence_version_map.tsv
      confirmed/confirmed_introns.tsv
      confirmed/introns.fasta
      confirmed/removed_insertions.fasta
    """

    require_options(input=input, out=out)

    detection_input = input
    orientation_outputs: Dict[str, str] = {}

    if orient:
        if orientation_method.lower() != "sina":
            raise typer.BadParameter("Only --orientation-method sina is currently supported.")

        preprocess_dir = Path(out) / "preprocess"
        step("Orienting input sequences with SINA before intron detection")

        orientation_outputs = run_sina_orientation(
            input_fasta=input,
            outdir=preprocess_dir,
            sina_ref=sina_ref,
            threads=threads,
            score_margin=orientation_score_margin,
            min_score=orientation_min_score,
            kmer_size=orientation_kmer_size,
        )
        detection_input = Path(orientation_outputs["oriented_fasta"])

    reference_db = resolve_reference_for_detect_intron(db)

    taxonomy_path = species_taxonomy
    if taxonomy_path is None:
        taxonomy_path = get_reference_path("species-rep1-clean-taxonomy", required=False)

    if taxonomy_path is None:
        taxonomy_path = get_reference_path("species-rep1-taxonomy", required=False)

    blacklist_path = reference_blacklist
    if blacklist_path is None:
        blacklist_path = get_reference_path("reference-intron-blacklist", required=False)

    blastn = get_software("blastn", required=True)
    makeblastdb = get_software("makeblastdb", required=False)

    step("Detecting intron-like insertions with BLAST local HSPs")

    outputs = detect_introns(
        input_fasta=detection_input,
        db=reference_db,
        outdir=out,
        source_label=source_label,
        threads=threads,
        hsp_identity=hsp_id,
        min_intron_len=min_intron_len,
        min_flank_len=min_flank_len,
        max_ref_gap=max_ref_gap,
        ref_overlap_tolerance=ref_overlap_tolerance,
        max_hsps=max_hsps,
        maxaccepts=max_targets,
        strand=strand,
        confirm_id=confirm_id,
        confirm_qcov=confirm_qcov,
        min_support_refs=min_support_refs,
        min_support_species=min_support_species,
        min_support_genera=min_support_genera,
        coordinate_tolerance=coordinate_tolerance,
        species_taxonomy=taxonomy_path,
        reference_blacklist=blacklist_path,
        blastn=blastn,
        makeblastdb=makeblastdb,
    )

    if orientation_outputs:
        outputs = {
            **{f"orientation_{key}": value for key, value in orientation_outputs.items()},
            **outputs,
        }

    print_outputs(outputs)
    success(f"Intron detection finished: {out}")


@app.command("insert-backbone", help="Insert extension sequences into the SILVA backbone.")
def insert_backbone_cmd(
    input: Optional[Path] = typer.Option(None, "--input", help="Input FASTA."),
    out: Optional[Path] = typer.Option(None, "--out", help="Output directory."),
    source_label: Optional[str] = typer.Option(None, "--source-label", help="Dataset label."),
    silva_manifest: Optional[Path] = typer.Option(
        None,
        "--silva-manifest",
        help="Reference manifest from prepare-silva. Default: configured ref-manifest.",
    ),
    original_fasta: Optional[Path] = typer.Option(
        None,
        "--original-fasta",
        help="Original FASTA used when analysis sequences were modified.",
    ),
    version_map: Optional[Path] = typer.Option(
        None,
        "--version-map",
        help="Map from analysis sequence IDs back to original sequence IDs.",
    ),
    rank_thresholds: str = typer.Option(
        "default",
        "--rank-thresholds",
        help="Rank threshold preset or custom rank:identity specification.",
    ),
    db_format: str = typer.Option(
        "auto",
        "--db-format",
        help="Reference database format: auto, udb, or fasta.",
    ),
    maxaccepts: int = typer.Option(1, "--maxaccepts", min=1, help="Maximum hits per query."),
    strand: str = typer.Option("both", "--strand", help="Search strand mode: both or plus."),
    centroid_priority: str = typer.Option(
        "original_seed",
        "--centroid-priority",
        help="Centroid choice rule: original_seed or prefer_non_intron.",
    ),
    threads: int = typer.Option(DEFAULT_THREADS, "--threads", min=1, help="Number of threads."),
) -> None:
    """
    Insert new or extension sequences into the SILVA backbone.
    """

    require_options(input=input, out=out, source_label=source_label)

    manifest_path = silva_manifest or get_ref_manifest(required=True)
    vsearch = get_vsearch(required=True)

    step("Inserting sequences into SILVA backbone")
    outputs = call_supported(
        insert_backbone,
        input_fasta=input,
        manifest_path=manifest_path,
        outdir=out,
        source_label=source_label,
        rank_thresholds=rank_thresholds,
        original_fasta=original_fasta,
        version_map_path=version_map,
        db_format=db_format,
        vsearch=vsearch,
        threads=threads,
        strand=strand,
        maxaccepts=maxaccepts,
        centroid_policy=centroid_priority,
        centroid_priority=centroid_priority,
    )

    print_outputs(outputs)
    success(f"Backbone insertion finished: {out}")


@app.command("overlap-backbone", help="Compare taxon overlap across backbone assignments.")
def overlap_backbone_cmd(
    assignments: Optional[List[Path]] = typer.Option(
        None,
        "--assignments",
        help="One or more backbone assignment TSV files.",
    ),
    labels: Optional[str] = typer.Option(
        None,
        "--labels",
        help="Comma-separated labels matching --assignments.",
    ),
    out: Optional[Path] = typer.Option(None, "--out", help="Output directory."),
) -> None:
    """
    Compare multiple backbone assignment files.
    """

    require_options(assignments=assignments, labels=labels, out=out)

    if assignments is None:
        raise typer.BadParameter("--assignments is required.")

    parsed_labels = [item.strip() for item in labels.split(",") if item.strip()]

    if len(parsed_labels) != len(assignments):
        raise typer.BadParameter("--labels must contain the same number of labels as --assignments.")

    outputs = overlap_backbone(assignments, parsed_labels, out)
    print_outputs(outputs)


@app.command("derep", help="Run standardized VSEARCH dereplication.")
def derep_cmd(
    input: Optional[Path] = typer.Option(None, "--input", help="Input FASTA."),
    out: Optional[Path] = typer.Option(None, "--out", help="Output directory."),
    sort: bool = typer.Option(False, "--sort", help="Also run VSEARCH sortbysize."),
    minsize: Optional[int] = typer.Option(None, "--minsize", min=1, help="Minimum abundance."),
    threads: int = typer.Option(DEFAULT_THREADS, "--threads", min=1, help="Number of threads."),
) -> None:
    """
    Dereplicate identical sequences with VSEARCH.
    """

    require_options(input=input, out=out)

    vsearch = get_vsearch(required=True)

    step("Running dereplication")
    derep_fasta = dereplicate(input, out, vsearch, threads)

    if sort:
        step("Sorting dereplicated sequences by size")
        sorted_fasta = sort_by_size(derep_fasta, out, vsearch, minsize)
        typer.echo(sorted_fasta)
    else:
        typer.echo(derep_fasta)


@app.command("cluster", help="Cluster sequences at one or more identity levels.")
def cluster_cmd(
    input: Optional[Path] = typer.Option(None, "--input", help="Input FASTA."),
    out: Optional[Path] = typer.Option(None, "--out", help="Output directory."),
    ids: str = typer.Option("0.99,0.97,0.95,0.90", "--ids", help="Comma-separated identity thresholds."),
    method: str = typer.Option("cluster_size", "--method", help="cluster_size, cluster_fast, or cluster_smallmem."),
    relabel: Optional[str] = typer.Option(None, "--relabel", help="Optional centroid relabel prefix."),
    threads: int = typer.Option(DEFAULT_THREADS, "--threads", min=1, help="Number of threads."),
) -> None:
    """
    Cluster representative sequences at one or more identity thresholds.
    """

    require_options(input=input, out=out)

    vsearch = get_vsearch(required=True)

    for identity in parse_id_list(ids):
        step(f"Clustering at identity {identity}")
        result = run_cluster(
            input_fasta=input,
            outdir=out,
            identity=identity,
            method=method,
            vsearch=vsearch,
            threads=threads,
            relabel=relabel,
        )
        typer.echo(f"{identity}\t{result['centroids']}\t{result['uc']}")


@app.command("classify", help="Classify representative sequences.")
def classify_cmd(
    input: Optional[Path] = typer.Option(None, "--input", help="Representative FASTA."),
    out: Optional[Path] = typer.Option(None, "--out", help="Output directory."),
    ref_manifest: Optional[Path] = typer.Option(None, "--ref-manifest", help="Reference manifest."),
    silva_fasta: Optional[Path] = typer.Option(None, "--silva-fasta", help="Override SILVA FASTA."),
    silva_sintax: Optional[Path] = typer.Option(None, "--silva-sintax", help="Override SILVA SINTAX FASTA."),
    typestrains: Optional[Path] = typer.Option(None, "--typestrains", help="Override type-strain FASTA."),
    sintax_cutoff: float = typer.Option(0.8, "--sintax-cutoff", min=0.0, max=1.0, help="SINTAX cutoff."),
    min_id: float = typer.Option(0.70, "--min-id", min=0.0, max=1.0, help="Minimum global-search identity."),
    maxaccepts: int = typer.Option(10, "--maxaccepts", min=1, help="Maximum hits retained per query."),
    strand: str = typer.Option("both", "--strand", help="Search strand mode: both or plus."),
    threads: int = typer.Option(DEFAULT_THREADS, "--threads", min=1, help="Number of threads."),
) -> None:
    """
    Classify representative sequences using SINTAX, SILVA hits, and type-strain hits.
    """

    require_options(input=input, out=out)

    manifest: Dict[str, str] = {}
    configured_manifest = get_ref_manifest(required=False)

    if ref_manifest or configured_manifest:
        _, manifest = load_manifest_or_config(ref_manifest)

    sf = resolve_ref_file(silva_fasta, manifest, "silva_fasta", "silva-fasta", "SILVA FASTA")
    ss = resolve_ref_file(silva_sintax, manifest, "silva_sintax_fasta", "silva-sintax", "SILVA SINTAX FASTA")
    ts = resolve_ref_file(typestrains, manifest, "silva_typestrains_fasta", "typestrains", "type-strain FASTA")

    vsearch = get_vsearch(required=True)
    ensure_dir(out)

    step("Running SINTAX classification")
    sintax_tsv = sintax(input, ss, out, sintax_cutoff, vsearch, threads)

    step("Searching against SILVA FASTA")
    silva = usearch_global(input, sf, out, "silva_hits", min_id, maxaccepts, strand, vsearch, threads)

    step("Searching against type-strain FASTA")
    typehit = usearch_global(input, ts, out, "typestrain_hits", min_id, maxaccepts, strand, vsearch, threads)

    summary = str(Path(out) / "taxonomy_summary.tsv")

    step("Writing taxonomy summary")
    summarize(sintax_tsv, silva["userout"], typehit["userout"], summary)

    typer.echo(f"summary\t{summary}")
    success(f"Classification finished: {out}")


@app.command("assign", help="Assign new sequences to old centroids or create new clusters.")
def assign_cmd(
    new: Optional[Path] = typer.Option(None, "--new", help="New FASTA."),
    old_centroids: Optional[Path] = typer.Option(None, "--old-centroids", help="Existing centroid FASTA."),
    out: Optional[Path] = typer.Option(None, "--out", help="Output directory."),
    identity: Optional[float] = typer.Option(None, "--id", min=0.0, max=1.0, help="Assignment identity threshold."),
    method: str = typer.Option("cluster_size", "--method", help="VSEARCH clustering method for unassigned sequences."),
    new_cluster_prefix: str = typer.Option("new_", "--new-cluster-prefix", help="Prefix for newly created cluster IDs."),
    source_map: Optional[Path] = typer.Option(None, "--source-map", help="Optional source map."),
    maxaccepts: int = typer.Option(1, "--maxaccepts", min=1, help="Maximum hits retained per query."),
    strand: str = typer.Option("both", "--strand", help="Search strand mode: both or plus."),
    threads: int = typer.Option(DEFAULT_THREADS, "--threads", min=1, help="Number of threads."),
) -> None:
    """
    Assign new sequences to existing centroids.
    """

    require_options(new=new, old_centroids=old_centroids, out=out, identity=identity)

    vsearch = get_vsearch(required=True)

    step("Assigning new sequences to old centroids")
    outputs = call_supported(
        assign_or_create,
        new_fasta=new,
        old_centroids=old_centroids,
        outdir=out,
        identity=identity,
        vsearch=vsearch,
        threads=threads,
        strand=strand,
        new_cluster_prefix=new_cluster_prefix,
        cluster_method=method,
        method=method,
        maxaccepts=maxaccepts,
    )

    if source_map:
        annotated = str(Path(out) / "assignments_with_source.tsv")
        annotate_assignments_with_source(outputs["assignments"], source_map, annotated)
        outputs["assignments_with_source"] = annotated

    print_outputs(outputs)
    success(f"Assignment finished: {out}")


@app.command("provenance", help="Summarize source composition from UC clustering levels.")
def provenance_cmd(
    source_map: Optional[Path] = typer.Option(None, "--source-map", help="TSV mapping sequence IDs to source labels."),
    derep_uc: Optional[Path] = typer.Option(None, "--derep-uc", help="Dereplication UC file."),
    level_uc: Optional[List[Path]] = typer.Option(None, "--level-uc", help="One or more cluster-level UC files."),
    level_labels: Optional[str] = typer.Option(None, "--level-labels", help="Comma-separated labels for --level-uc."),
    out: Optional[Path] = typer.Option(None, "--out", help="Output directory."),
) -> None:
    """
    Summarize the source composition of dereplicated and clustered sequences.
    """

    require_options(source_map=source_map, derep_uc=derep_uc, level_uc=level_uc, out=out)

    labels = [item.strip() for item in level_labels.split(",")] if level_labels else None

    outputs = provenance_from_uc_levels(source_map, derep_uc, level_uc, out, labels)
    print_outputs(outputs)


@app.command("summarize", help="Summarize existing classify outputs.")
def summarize_cmd(
    sintax_file: Optional[Path] = typer.Option(None, "--sintax", help="SINTAX output TSV."),
    silva_hits: Optional[Path] = typer.Option(None, "--silva-hits", help="SILVA global-search hits table."),
    typestrain_hits: Optional[Path] = typer.Option(None, "--typestrain-hits", help="Type-strain global-search hits table."),
    out: Optional[Path] = typer.Option(None, "--out", help="Output summary TSV."),
) -> None:
    """
    Rebuild taxonomy_summary.tsv from existing classify intermediate files.
    """

    require_options(
        sintax=sintax_file,
        silva_hits=silva_hits,
        typestrain_hits=typestrain_hits,
        out=out,
    )

    summarize(sintax_file, silva_hits, typestrain_hits, out)
    typer.echo(out)


@app.command("map-intron-coordinates", help="Map confirmed intron insertion sites to a coordinate reference.")
def map_intron_coordinates_cmd(
    confirmed_introns: Optional[Path] = typer.Option(
        None,
        "--confirmed-introns",
        help="Confirmed intron table from detect-intron, usually confirmed/confirmed_introns.tsv.",
    ),
    analysis_fasta: Optional[Path] = typer.Option(
        None,
        "--analysis-fasta",
        help="Cleaned analysis FASTA from detect-intron, usually confirmed/analysis_sequences.fasta.",
    ),
    coordinate_ref: Optional[Path] = typer.Option(
        None,
        "--coordinate-ref",
        help="Coordinate reference FASTA or BLAST DB prefix. Default: configured coordinate_ref.",
    ),
    out: Path = typer.Option(
        Path("intron_coordinates"),
        "--out",
        help="Output directory for coordinate mapping results.",
    ),
    min_identity: float = typer.Option(
        80.0,
        "--min-identity",
        min=0.0,
        max=100.0,
        help="Minimum identity percentage for mapping cleaned sequences to the coordinate reference.",
    ),
    max_targets: int = typer.Option(
        5,
        "--max-targets",
        min=1,
        help="Maximum coordinate-reference hits retained per cleaned sequence.",
    ),
    threads: int = typer.Option(
        DEFAULT_THREADS,
        "--threads",
        min=1,
        help="Number of BLAST threads.",
    ),
) -> None:
    """
    Map confirmed intron insertion sites to a common coordinate reference.

    This command uses cleaned analysis sequences, not original intron-containing
    sequences, because the original sequences can produce broken alignments
    around the insertion site.

    Output:
      coordinate_alignment_hsps.tsv
      intron_coordinate_map.tsv
      intron_position_distribution.tsv
    """

    require_options(confirmed_introns=confirmed_introns, analysis_fasta=analysis_fasta)

    cref = coordinate_ref or get_reference_path("coordinate-ref", required=True)

    blastn = get_software("blastn", required=True)

    if blastdb_prefix_exists(cref):
        makeblastdb = get_software("makeblastdb", required=False)
    else:
        makeblastdb = get_software("makeblastdb", required=True)

    step("Mapping confirmed intron insertion sites to coordinate reference")

    outputs = map_intron_coordinates(
        confirmed_introns=confirmed_introns,
        analysis_fasta=analysis_fasta,
        coordinate_ref=cref,
        outdir=out,
        blastn=blastn,
        makeblastdb=makeblastdb,
        threads=threads,
        min_identity=min_identity,
        max_target_seqs=max_targets,
    )

    print_outputs(outputs)
    success(f"Intron coordinate mapping finished: {out}")


@app.command("analyze-intron", help="Analyze confirmed intron-like insertion sequences.")
def analyze_intron_cmd(
    introns: Optional[Path] = typer.Option(
        None,
        "--introns",
        help="Confirmed intron FASTA, usually detect_intron_out/confirmed/introns.fasta.",
    ),
    metadata: Optional[Path] = typer.Option(
        None,
        "--metadata",
        help="Confirmed intron table. Default: confirmed_introns.tsv beside --introns when present.",
    ),
    coordinate_map: Optional[Path] = typer.Option(
        None,
        "--coordinate-map",
        help="Optional intron coordinate map from coordinate-reference mapping.",
    ),
    out: Path = typer.Option(
        Path("intron_analysis"),
        "--out",
        help="Output directory for intron characterization results.",
    ),
    cluster: bool = typer.Option(
        True,
        "--cluster/--no-cluster",
        help="Cluster intron sequences with VSEARCH.",
    ),
    cluster_ids: str = typer.Option(
        "0.99,0.97",
        "--cluster-ids",
        help="Comma-separated intron clustering identity thresholds.",
    ),
    mafft_alignments: bool = typer.Option(
        True,
        "--mafft/--no-mafft",
        help="Run MAFFT multiple sequence alignment on confirmed introns.",
    ),
    motif: bool = typer.Option(
        True,
        "--motif/--no-motif",
        help="Run motif discovery with available MEME/STREME tools and optional FIMO scanning.",
    ),
    structure: bool = typer.Option(
        True,
        "--structure/--no-structure",
        help="Run RNAfold and RNAalifold when available.",
    ),
    motif_minw: int = typer.Option(
        6,
        "--motif-minw",
        min=1,
        help="Minimum motif width for MEME/STREME.",
    ),
    motif_maxw: int = typer.Option(
        50,
        "--motif-maxw",
        min=1,
        help="Maximum motif width for MEME/STREME.",
    ),
    motif_nmotifs: int = typer.Option(
        10,
        "--motif-nmotifs",
        min=1,
        help="Maximum number of motifs for MEME.",
    ),
    threads: int = typer.Option(
        DEFAULT_THREADS,
        "--threads",
        min=1,
        help="Number of threads for VSEARCH, MAFFT, and MEME.",
    ),
) -> None:
    """
    Analyze confirmed intron-like insertion sequences.

    This command is downstream of detect-intron. It does not decide whether a
    query sequence contains an intron. It characterizes already confirmed and
    extracted intron sequences.
    """

    require_options(introns=introns)

    if metadata is None and introns is not None:
        candidate_metadata = introns.parent / "confirmed_introns.tsv"
        if candidate_metadata.exists():
            metadata = candidate_metadata

    if motif_maxw < motif_minw:
        raise typer.BadParameter("--motif-maxw must be greater than or equal to --motif-minw.")

    vsearch = get_vsearch(required=True) if cluster else None
    mafft = optional_software("mafft") if mafft_alignments else None
    meme = optional_software("meme") if motif else None
    streme = optional_software("streme") if motif else None
    fimo = optional_software("fimo") if motif else None
    rnafold = optional_software("rnafold") if structure else None
    rnaalifold = optional_software("rnaalifold") if structure else None

    step("Analyzing confirmed intron-like insertions")

    outputs = analyze_introns(
        introns_fasta=introns,
        metadata=metadata,
        coordinate_map=coordinate_map,
        outdir=out,
        vsearch=vsearch,
        mafft=mafft,
        meme=meme,
        streme=streme,
        fimo=fimo,
        rnafold=rnafold,
        rnaalifold=rnaalifold,
        threads=threads,
        run_cluster=cluster,
        run_mafft=mafft_alignments,
        run_motif=motif,
        run_structure=structure,
        cluster_identities=parse_id_list(cluster_ids),
        motif_minw=motif_minw,
        motif_maxw=motif_maxw,
        motif_nmotifs=motif_nmotifs,
    )

    print_outputs(outputs)
    success(f"Intron analysis finished: {out}")


@app.command("run", help="Run dereplication, clustering, and classification.")
def run_cmd(
    input: Optional[Path] = typer.Option(None, "--input", help="Input FASTA."),
    out: Optional[Path] = typer.Option(None, "--out", help="Output project directory."),
    ref_manifest: Optional[Path] = typer.Option(None, "--ref-manifest", help="Reference manifest."),
    source_map: Optional[Path] = typer.Option(None, "--source-map", help="Optional source map for provenance."),
    ids: str = typer.Option("0.99,0.97,0.95,0.90", "--ids", help="Comma-separated clustering identity thresholds."),
    classify_id: float = typer.Option(0.99, "--classify-id", min=0.0, max=1.0, help="Identity level to classify."),
    method: str = typer.Option("cluster_size", "--method", help="VSEARCH clustering method."),
    minsize: Optional[int] = typer.Option(None, "--minsize", min=1, help="Minimum abundance for sortbysize."),
    sintax_cutoff: float = typer.Option(0.8, "--sintax-cutoff", min=0.0, max=1.0, help="SINTAX cutoff."),
    min_id: float = typer.Option(0.70, "--min-id", min=0.0, max=1.0, help="Minimum global-search identity."),
    maxaccepts: int = typer.Option(10, "--maxaccepts", min=1, help="Maximum hits retained per query."),
    strand: str = typer.Option("both", "--strand", help="Search strand mode: both or plus."),
    threads: int = typer.Option(DEFAULT_THREADS, "--threads", min=1, help="Number of threads."),
) -> None:
    """
    Run the standard AutoTax2 workflow.
    """

    require_options(input=input, out=out)

    manifest_path, manifest = load_manifest_or_config(ref_manifest)

    sf = resolve_ref_file(None, manifest, "silva_fasta", "silva-fasta", "SILVA FASTA")
    ss = resolve_ref_file(None, manifest, "silva_sintax_fasta", "silva-sintax", "SILVA SINTAX FASTA")
    ts = resolve_ref_file(None, manifest, "silva_typestrains_fasta", "typestrains", "type-strain FASTA")

    vsearch = get_vsearch(required=True)

    out_path = Path(out)
    work = out_path / "work"
    clusters = out_path / "clusters"
    classify_dir = out_path / "classify"

    ensure_dir(work)
    ensure_dir(clusters)
    ensure_dir(classify_dir)

    step(f"Using reference manifest: {manifest_path}")

    step("Dereplicating input FASTA")
    derep_fasta = dereplicate(input, str(work), vsearch, threads)

    step("Sorting dereplicated FASTA by size")
    sorted_fasta = sort_by_size(derep_fasta, str(work), vsearch, minsize)

    parsed_ids = parse_id_list(ids)
    chosen_centroids = None

    for identity in parsed_ids:
        step(f"Clustering at identity {identity}")
        result = run_cluster(
            input_fasta=sorted_fasta,
            outdir=str(clusters),
            identity=identity,
            method=method,
            vsearch=vsearch,
            threads=threads,
            relabel=None,
        )

        if abs(identity - classify_id) < 1e-9:
            chosen_centroids = result["centroids"]

    if chosen_centroids is None:
        first = parsed_ids[0]
        chosen_centroids = str(clusters / f"otu{int(round(first * 100)):03d}_centroids.fasta")
        step(f"--classify-id {classify_id} not found in --ids; classifying first level {first}")

    step("Running SINTAX classification")
    sintax_tsv = sintax(chosen_centroids, ss, str(classify_dir), sintax_cutoff, vsearch, threads)

    step("Searching selected centroids against SILVA FASTA")
    silva = usearch_global(chosen_centroids, sf, str(classify_dir), "silva_hits", min_id, maxaccepts, strand, vsearch, threads)

    step("Searching selected centroids against type-strain FASTA")
    typehit = usearch_global(chosen_centroids, ts, str(classify_dir), "typestrain_hits", min_id, maxaccepts, strand, vsearch, threads)

    summary = str(classify_dir / "taxonomy_summary.tsv")

    step("Writing taxonomy summary")
    summarize(sintax_tsv, silva["userout"], typehit["userout"], summary)

    if source_map:
        step("Writing provenance report")
        level_ucs = [
            str(clusters / f"otu{int(round(identity * 100)):03d}.uc")
            for identity in parsed_ids
        ]
        level_labels = [str(identity) for identity in parsed_ids]
        provenance_from_uc_levels(
            source_map,
            str(work / "derep.uc"),
            level_ucs,
            str(out_path / "provenance"),
            level_labels,
        )

    typer.echo(f"summary\t{summary}")
    success(f"AutoTax2 run finished: {out_path}")


def main() -> None:
    app()


if __name__ == "__main__":
    main()