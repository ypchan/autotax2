from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import typer

from .assign import assign_or_create
from .backbone import insert_backbone
from .dependencies import check_core_dependencies, check_paths, summarize_checks
from .intron import detect_introns
from .logging import print_example, success
from .overlap import overlap_backbone
from .prepare import prepare_silva
from .provenance import annotate_assignments_with_source, provenance_from_uc_levels
from .summarize import summarize
from .threads import resolve_threads
from .utils import ensure_dir, parse_manifest
from .vsearch import cluster as run_cluster
from .vsearch import dereplicate, sintax, sort_by_size, usearch_global

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}

app = typer.Typer(
    name="autotax2",
    help="AutoTax2 workflow.",
    no_args_is_help=True,
    add_completion=False,
    context_settings=CONTEXT_SETTINGS,
    rich_markup_mode="rich",
)

EXAMPLES: Dict[str, str] = {
    "check": "autotax2 check --ref-manifest refdatabases/autotax2_ref_manifest.tsv",
    "prepare-silva": """autotax2 prepare-silva \\
  --silva-fasta SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \\
  --silva-metadata SILVA_138.2_SSURef.full_metadata.gz \\
  --out refdatabases \\
  --make-udb \\
  --threads auto""",
    "detect-intron": """autotax2 detect-intron \\
  --input hifimeta.original.fa \\
  --db refdatabases/SILVA_138.2_SSURef_NR99_tax_silva.udb \\
  --source-label hifimeta \\
  --out hifimeta_intron \\
  --rescue-id 0.987 \\
  --min-intron-len 50 \\
  --threads auto""",
    "insert-backbone": """autotax2 insert-backbone \\
  --input hifimeta_intron/analysis_sequences.fa \\
  --original-fasta hifimeta.original.fa \\
  --version-map hifimeta_intron/sequence_version_map.tsv \\
  --source-label hifimeta \\
  --silva-manifest refdatabases/autotax2_ref_manifest.tsv \\
  --rank-thresholds default \\
  --db-format auto \\
  --out hifimeta_inserted \\
  --threads auto""",
    "overlap-backbone": """autotax2 overlap-backbone \\
  --assignments ref2/sequence_rank_assignment.tsv \\
  --assignments ref3/sequence_rank_assignment.tsv \\
  --assignments hifimeta/sequence_rank_assignment.tsv \\
  --labels ref2,ref3,hifimeta \\
  --out ref_overlap""",
    "derep": "autotax2 derep --input input.fa --out work --sort --threads auto",
    "cluster": "autotax2 cluster --input work/derep_sorted.fa --out clusters --ids 0.99,0.97 --threads auto",
    "classify": "autotax2 classify --input clusters/otu099_centroids.fa --ref-manifest refdatabases/autotax2_ref_manifest.tsv --out classify",
    "assign": "autotax2 assign --new new.fa --old-centroids centroids_97.fa --id 0.97 --out assign_97 --threads auto",
    "provenance": "autotax2 provenance --source-map source_map.tsv --derep-uc work/derep.uc --level-uc clusters/otu099.uc --level-labels 0.99 --out provenance",
    "summarize": "autotax2 summarize --sintax classify/sintax.tsv --silva-hits classify/silva_hits.tsv --typestrain-hits classify/typestrain_hits.tsv --out classify/taxonomy_summary.tsv",
    "run": "autotax2 run --input input.fa --ref-manifest refdatabases/autotax2_ref_manifest.tsv --out autotax2_out --threads auto",
}


def show_example(name: str) -> None:
    print_example(f"Example: autotax2 {name}", EXAMPLES[name])


def require_options(example: bool, **values: object) -> None:
    if example:
        return
    missing = [k.replace("_", "-") for k, v in values.items() if v in (None, "", [])]
    if missing:
        raise typer.BadParameter("Missing required option(s): " + ", ".join(f"--{x}" for x in missing))


def print_outputs(outputs: Dict[str, object]) -> None:
    for key, value in outputs.items():
        typer.echo(f"{key}\t{value}")


def parse_id_list(value: str) -> List[float]:
    try:
        return [float(x.strip()) for x in value.split(",") if x.strip()]
    except ValueError as exc:
        raise typer.BadParameter("--ids must be a comma-separated list of numbers, e.g. 0.99,0.97") from exc


@app.command("check", help="Check external dependencies and optional reference files.")
def check(
    vsearch: str = typer.Option("vsearch", help="Path to VSEARCH."),
    ref_manifest: Optional[Path] = typer.Option(None, "--ref-manifest", help="Reference manifest from prepare-silva."),
    source_map: Optional[Path] = typer.Option(None, "--source-map", help="Optional source map file."),
    example: bool = typer.Option(False, "--example", help="Print an example command and exit."),
) -> None:
    if example:
        show_example("check")
        return
    ok = check_core_dependencies(vsearch)
    extra = []
    if ref_manifest:
        manifest = parse_manifest(ref_manifest)
        extra.extend([
            ("silva_fasta", manifest.get("silva_fasta", ""), True),
            ("silva_sintax_fasta", manifest.get("silva_sintax_fasta", ""), True),
            ("silva_typestrains_fasta", manifest.get("silva_typestrains_fasta", ""), True),
        ])
    if source_map:
        extra.append(("source_map", source_map, True))
    if extra:
        ok = summarize_checks(check_paths(extra)) and ok
    if not ok:
        raise typer.Exit(1)
    success("All required dependency checks passed.")


@app.command("prepare-silva", help="Prepare local SILVA FASTA and metadata files.")
def prepare_silva_cmd(
    silva_fasta: Optional[Path] = typer.Option(None, "--silva-fasta", help="Local SILVA FASTA, optionally .gz."),
    silva_metadata: Optional[Path] = typer.Option(None, "--silva-metadata", help="Local SILVA full_metadata, optionally .gz."),
    out: Path = typer.Option(Path("refdatabases"), "--out", help="Output directory."),
    prefix: Optional[str] = typer.Option(None, "--prefix", help="Output prefix."),
    make_udb: bool = typer.Option(False, "--make-udb", help="Build VSEARCH UDB files."),
    clean_sequences: bool = typer.Option(True, "--clean-sequences/--no-clean-sequences", help="Clean FASTA before processing."),
    allowed_bases: str = typer.Option("ACGT", "--allowed-bases", help="Allowed bases after U->T; default drops N and ambiguity codes."),
    u_to_t: bool = typer.Option(True, "--u-to-t/--no-u-to-t", help="Convert U/u to T before filtering."),
    vsearch: str = typer.Option("vsearch", help="Path to VSEARCH."),
    threads: str = typer.Option("auto", help="Accepted for CLI consistency."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Print external commands without running them."),
    example: bool = typer.Option(False, "--example", help="Print an example command and exit."),
) -> None:
    if example:
        show_example("prepare-silva")
        return
    require_options(False, silva_fasta=silva_fasta, silva_metadata=silva_metadata)
    resolve_threads(threads)
    entries = prepare_silva(
        silva_fasta=silva_fasta,
        silva_metadata=silva_metadata,
        outdir=out,
        prefix=prefix,
        make_udb_files=make_udb,
        vsearch=vsearch,
        dry_run=dry_run,
        clean_sequences=clean_sequences,
        allowed_bases=allowed_bases,
        convert_u_to_t=u_to_t,
    )
    print_outputs(entries)


@app.command("detect-intron", help="Detect intron-like insertions and write analysis FASTA files.")
def detect_intron(
    input: Optional[Path] = typer.Option(None, "--input", help="Input FASTA."),
    db: Optional[Path] = typer.Option(None, "--db", help="SILVA FASTA or UDB database."),
    out: Optional[Path] = typer.Option(None, "--out", help="Output directory."),
    source_label: Optional[str] = typer.Option(None, "--source-label", help="Dataset label."),
    search_id: float = typer.Option(0.70, "--search-id", help="Initial low-stringency identity."),
    rescue_id: float = typer.Option(0.987, "--rescue-id", help="Rescue identity after intron removal."),
    min_intron_len: int = typer.Option(50, "--min-intron-len", help="Minimum intron length."),
    min_flank_len: int = typer.Option(150, "--min-flank-len", help="Minimum aligned flank length."),
    maxaccepts: int = typer.Option(1, "--maxaccepts", help="VSEARCH maxaccepts."),
    strand: str = typer.Option("both", "--strand", help="VSEARCH strand: both or plus."),
    vsearch: str = typer.Option("vsearch", help="Path to VSEARCH."),
    threads: str = typer.Option("auto", help="Thread count or auto."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Print external commands without running them."),
    example: bool = typer.Option(False, "--example", help="Print an example command and exit."),
) -> None:
    if example:
        show_example("detect-intron")
        return
    require_options(False, input=input, db=db, out=out)
    outputs = detect_introns(
        input_fasta=input,
        db=db,
        outdir=out,
        source_label=source_label,
        vsearch=vsearch,
        threads=resolve_threads(threads),
        search_id=search_id,
        rescue_id=rescue_id,
        min_intron_len=min_intron_len,
        min_flank_len=min_flank_len,
        strand=strand,
        maxaccepts=maxaccepts,
        dry_run=dry_run,
    )
    print_outputs(outputs)


@app.command("insert-backbone", help="Insert extension sequences into the SILVA backbone.")
def insert_backbone_cmd(
    input: Optional[Path] = typer.Option(None, "--input", help="Input FASTA."),
    silva_manifest: Optional[Path] = typer.Option(None, "--silva-manifest", help="Reference manifest from prepare-silva."),
    out: Optional[Path] = typer.Option(None, "--out", help="Output directory."),
    source_label: Optional[str] = typer.Option(None, "--source-label", help="Dataset label."),
    original_fasta: Optional[Path] = typer.Option(None, "--original-fasta", help="Original FASTA for centroid output."),
    version_map: Optional[Path] = typer.Option(None, "--version-map", help="Analysis-to-original map."),
    rank_thresholds: str = typer.Option("default", "--rank-thresholds", help="Rank threshold preset or custom spec."),
    db_format: str = typer.Option("auto", "--db-format", help="auto, udb, or fasta."),
    maxaccepts: int = typer.Option(1, "--maxaccepts", help="VSEARCH maxaccepts."),
    strand: str = typer.Option("both", "--strand", help="VSEARCH strand: both or plus."),
    centroid_policy: str = typer.Option("original_seed", "--centroid-policy", help="original_seed or prefer_non_intron."),
    vsearch: str = typer.Option("vsearch", help="Path to VSEARCH."),
    threads: str = typer.Option("auto", help="Thread count or auto."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Print external commands without running them."),
    example: bool = typer.Option(False, "--example", help="Print an example command and exit."),
) -> None:
    if example:
        show_example("insert-backbone")
        return
    require_options(False, input=input, silva_manifest=silva_manifest, out=out, source_label=source_label)
    outputs = insert_backbone(
        input_fasta=input,
        manifest_path=silva_manifest,
        outdir=out,
        source_label=source_label,
        rank_thresholds=rank_thresholds,
        original_fasta=original_fasta,
        version_map_path=version_map,
        db_format=db_format,
        vsearch=vsearch,
        threads=resolve_threads(threads),
        strand=strand,
        maxaccepts=maxaccepts,
        centroid_policy=centroid_policy,
        dry_run=dry_run,
    )
    print_outputs(outputs)


@app.command("overlap-backbone", help="Summarize taxon overlap across backbone assignment files.")
def overlap_backbone_cmd(
    assignments: Optional[List[Path]] = typer.Option(None, "--assignments", help="One or more assignment files."),
    labels: Optional[str] = typer.Option(None, "--labels", help="Comma-separated labels matching --assignments."),
    out: Optional[Path] = typer.Option(None, "--out", help="Output directory."),
    example: bool = typer.Option(False, "--example", help="Print an example command and exit."),
) -> None:
    if example:
        show_example("overlap-backbone")
        return
    require_options(False, assignments=assignments, labels=labels, out=out)
    parsed_labels = [x.strip() for x in labels.split(",") if x.strip()]
    outputs = overlap_backbone(assignments, parsed_labels, out)
    print_outputs(outputs)


@app.command("derep", help="Run VSEARCH dereplication.")
def derep_cmd(
    input: Optional[Path] = typer.Option(None, "--input", help="Input FASTA."),
    out: Optional[Path] = typer.Option(None, "--out", help="Output directory."),
    sort: bool = typer.Option(False, "--sort", help="Also run sortbysize."),
    minsize: Optional[int] = typer.Option(None, "--minsize", help="Minimum size for sortbysize."),
    vsearch: str = typer.Option("vsearch", help="Path to VSEARCH."),
    threads: str = typer.Option("auto", help="Thread count or auto."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Print external commands without running them."),
    example: bool = typer.Option(False, "--example", help="Print an example command and exit."),
) -> None:
    if example:
        show_example("derep")
        return
    require_options(False, input=input, out=out)
    t = resolve_threads(threads)
    derep = dereplicate(input, out, vsearch, t, dry_run)
    typer.echo(sort_by_size(derep, out, vsearch, minsize, dry_run) if sort else derep)


@app.command("cluster", help="Run VSEARCH clustering at one or more identity levels.")
def cluster_cmd(
    input: Optional[Path] = typer.Option(None, "--input", help="Input FASTA."),
    out: Optional[Path] = typer.Option(None, "--out", help="Output directory."),
    ids: str = typer.Option("0.99,0.97,0.95,0.90", "--ids", help="Comma-separated identity thresholds."),
    method: str = typer.Option("cluster_size", "--method", help="cluster_size, cluster_fast, or cluster_smallmem."),
    relabel: Optional[str] = typer.Option(None, "--relabel", help="Optional centroid relabel prefix."),
    vsearch: str = typer.Option("vsearch", help="Path to VSEARCH."),
    threads: str = typer.Option("auto", help="Thread count or auto."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Print external commands without running them."),
    example: bool = typer.Option(False, "--example", help="Print an example command and exit."),
) -> None:
    if example:
        show_example("cluster")
        return
    require_options(False, input=input, out=out)
    t = resolve_threads(threads)
    for identity in parse_id_list(ids):
        result = run_cluster(input, out, identity, method, vsearch, t, relabel, dry_run)
        typer.echo(f"{identity}\t{result['centroids']}\t{result['uc']}")


@app.command("classify", help="Classify representative sequences using SINTAX plus SILVA/type-strain hits.")
def classify_cmd(
    input: Optional[Path] = typer.Option(None, "--input", help="Representative FASTA."),
    out: Optional[Path] = typer.Option(None, "--out", help="Output directory."),
    ref_manifest: Optional[Path] = typer.Option(None, "--ref-manifest", help="Reference manifest."),
    silva_fasta: Optional[Path] = typer.Option(None, "--silva-fasta", help="SILVA FASTA override."),
    silva_sintax: Optional[Path] = typer.Option(None, "--silva-sintax", help="SILVA SINTAX FASTA override."),
    typestrains: Optional[Path] = typer.Option(None, "--typestrains", help="Type strain FASTA override."),
    sintax_cutoff: float = typer.Option(0.8, "--sintax-cutoff", help="VSEARCH SINTAX cutoff."),
    min_id: float = typer.Option(0.70, "--min-id", help="Minimum identity for global searches."),
    maxaccepts: int = typer.Option(10, "--maxaccepts", help="VSEARCH maxaccepts."),
    strand: str = typer.Option("both", "--strand", help="VSEARCH strand: both or plus."),
    vsearch: str = typer.Option("vsearch", help="Path to VSEARCH."),
    threads: str = typer.Option("auto", help="Thread count or auto."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Print external commands without running them."),
    example: bool = typer.Option(False, "--example", help="Print an example command and exit."),
) -> None:
    if example:
        show_example("classify")
        return
    require_options(False, input=input, out=out)
    t = resolve_threads(threads)
    manifest = parse_manifest(ref_manifest) if ref_manifest else {}
    sf = silva_fasta or manifest.get("silva_fasta")
    ss = silva_sintax or manifest.get("silva_sintax_fasta")
    ts = typestrains or manifest.get("silva_typestrains_fasta")
    if not sf or not ss or not ts:
        raise typer.BadParameter("Missing reference files. Use --ref-manifest or explicit reference paths.")
    ensure_dir(out)
    sintax_tsv = sintax(input, ss, out, sintax_cutoff, vsearch, t, dry_run)
    silva = usearch_global(input, sf, out, "silva_hits", min_id, maxaccepts, strand, vsearch, t, False, dry_run)
    typehit = usearch_global(input, ts, out, "typestrain_hits", min_id, maxaccepts, strand, vsearch, t, False, dry_run)
    summary = str(Path(out) / "taxonomy_summary.tsv")
    if not dry_run:
        summarize(sintax_tsv, silva["userout"], typehit["userout"], summary)
    typer.echo(f"summary\t{summary}")


@app.command("assign", help="Assign new sequences to old centroids or create new clusters.")
def assign_cmd(
    new: Optional[Path] = typer.Option(None, "--new", help="New FASTA."),
    old_centroids: Optional[Path] = typer.Option(None, "--old-centroids", help="Existing centroid FASTA."),
    out: Optional[Path] = typer.Option(None, "--out", help="Output directory."),
    id: Optional[float] = typer.Option(None, "--id", help="Assignment identity threshold."),
    method: str = typer.Option("cluster_size", "--method", help="cluster_size, cluster_fast, or cluster_smallmem."),
    new_cluster_prefix: str = typer.Option("new_", "--new-cluster-prefix", help="Prefix for new clusters."),
    source_map: Optional[Path] = typer.Option(None, "--source-map", help="Optional source map."),
    maxaccepts: int = typer.Option(1, "--maxaccepts", help="VSEARCH maxaccepts."),
    strand: str = typer.Option("both", "--strand", help="VSEARCH strand: both or plus."),
    vsearch: str = typer.Option("vsearch", help="Path to VSEARCH."),
    threads: str = typer.Option("auto", help="Thread count or auto."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Print external commands without running them."),
    example: bool = typer.Option(False, "--example", help="Print an example command and exit."),
) -> None:
    if example:
        show_example("assign")
        return
    require_options(False, new=new, old_centroids=old_centroids, out=out, id=id)
    outputs = assign_or_create(new, old_centroids, out, id, vsearch, resolve_threads(threads), strand, new_cluster_prefix, method, maxaccepts, dry_run)
    if source_map and not dry_run:
        annotated = str(Path(out) / "assignments_with_source.tsv")
        annotate_assignments_with_source(outputs["assignments"], source_map, annotated)
        outputs["assignments_with_source"] = annotated
    print_outputs(outputs)


@app.command("provenance", help="Summarize source composition from UC clustering levels.")
def provenance_cmd(
    source_map: Optional[Path] = typer.Option(None, "--source-map", help="TSV mapping sequence IDs to sources."),
    derep_uc: Optional[Path] = typer.Option(None, "--derep-uc", help="Dereplication UC file."),
    level_uc: Optional[List[Path]] = typer.Option(None, "--level-uc", help="One or more cluster-level UC files."),
    level_labels: Optional[str] = typer.Option(None, "--level-labels", help="Comma-separated labels for --level-uc."),
    out: Optional[Path] = typer.Option(None, "--out", help="Output directory."),
    example: bool = typer.Option(False, "--example", help="Print an example command and exit."),
) -> None:
    if example:
        show_example("provenance")
        return
    require_options(False, source_map=source_map, derep_uc=derep_uc, level_uc=level_uc, out=out)
    labels = [x.strip() for x in level_labels.split(",")] if level_labels else None
    outputs = provenance_from_uc_levels(source_map, derep_uc, level_uc, out, labels)
    print_outputs(outputs)


@app.command("summarize", help="Summarize existing classify outputs.")
def summarize_cmd(
    sintax_file: Optional[Path] = typer.Option(None, "--sintax", help="SINTAX output TSV."),
    silva_hits: Optional[Path] = typer.Option(None, "--silva-hits", help="SILVA hits table."),
    typestrain_hits: Optional[Path] = typer.Option(None, "--typestrain-hits", help="Type strain hits table."),
    out: Optional[Path] = typer.Option(None, "--out", help="Output TSV."),
    example: bool = typer.Option(False, "--example", help="Print an example command and exit."),
) -> None:
    if example:
        show_example("summarize")
        return
    require_options(False, sintax=sintax_file, silva_hits=silva_hits, typestrain_hits=typestrain_hits, out=out)
    summarize(sintax_file, silva_hits, typestrain_hits, out)
    typer.echo(out)


@app.command("run", help="Run dereplication, clustering and classification in one workflow.")
def run_cmd(
    input: Optional[Path] = typer.Option(None, "--input", help="Input FASTA."),
    ref_manifest: Optional[Path] = typer.Option(None, "--ref-manifest", help="Reference manifest."),
    out: Optional[Path] = typer.Option(None, "--out", help="Output project directory."),
    source_map: Optional[Path] = typer.Option(None, "--source-map", help="Optional source map."),
    ids: str = typer.Option("0.99,0.97,0.95,0.90", "--ids", help="Comma-separated identity thresholds."),
    classify_id: float = typer.Option(0.99, "--classify-id", help="Identity level to classify."),
    method: str = typer.Option("cluster_size", "--method", help="cluster_size, cluster_fast, or cluster_smallmem."),
    minsize: Optional[int] = typer.Option(None, "--minsize", help="Minimum abundance for sortbysize."),
    sintax_cutoff: float = typer.Option(0.8, "--sintax-cutoff", help="VSEARCH SINTAX cutoff."),
    min_id: float = typer.Option(0.70, "--min-id", help="Minimum identity for global searches."),
    maxaccepts: int = typer.Option(10, "--maxaccepts", help="VSEARCH maxaccepts."),
    strand: str = typer.Option("both", "--strand", help="VSEARCH strand: both or plus."),
    vsearch: str = typer.Option("vsearch", help="Path to VSEARCH."),
    threads: str = typer.Option("auto", help="Thread count or auto."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Print external commands without running them."),
    example: bool = typer.Option(False, "--example", help="Print an example command and exit."),
) -> None:
    if example:
        show_example("run")
        return
    require_options(False, input=input, ref_manifest=ref_manifest, out=out)
    t = resolve_threads(threads)
    out_path = Path(out)
    work = out_path / "work"
    clusters = out_path / "clusters"
    classify_dir = out_path / "classify"
    ensure_dir(work); ensure_dir(clusters); ensure_dir(classify_dir)

    derep = dereplicate(input, str(work), vsearch, t, dry_run)
    sorted_fa = sort_by_size(derep, str(work), vsearch, minsize, dry_run)
    parsed_ids = parse_id_list(ids)
    chosen = None
    for identity in parsed_ids:
        result = run_cluster(sorted_fa, str(clusters), identity, method, vsearch, t, None, dry_run)
        if abs(identity - classify_id) < 1e-9:
            chosen = result["centroids"]
    if chosen is None:
        first = parsed_ids[0]
        chosen = str(clusters / f"otu{str(first).replace('.', '')}_centroids.fa")
        typer.echo(f"[autotax2] --classify-id not in --ids; classifying first level {first}")

    manifest = parse_manifest(ref_manifest)
    sf = manifest.get("silva_fasta")
    ss = manifest.get("silva_sintax_fasta")
    ts = manifest.get("silva_typestrains_fasta")
    if not sf or not ss or not ts:
        raise typer.BadParameter("Reference manifest is missing silva_fasta, silva_sintax_fasta, or silva_typestrains_fasta.")

    sintax_tsv = sintax(chosen, ss, str(classify_dir), sintax_cutoff, vsearch, t, dry_run)
    silva = usearch_global(chosen, sf, str(classify_dir), "silva_hits", min_id, maxaccepts, strand, vsearch, t, False, dry_run)
    typehit = usearch_global(chosen, ts, str(classify_dir), "typestrain_hits", min_id, maxaccepts, strand, vsearch, t, False, dry_run)
    summary = str(classify_dir / "taxonomy_summary.tsv")
    if not dry_run:
        summarize(sintax_tsv, silva["userout"], typehit["userout"], summary)
        if source_map:
            level_ucs = [str(clusters / f"otu{str(x).replace('.', '')}.uc") for x in parsed_ids]
            level_labels = [str(x) for x in parsed_ids]
            provenance_from_uc_levels(source_map, str(work / "derep.uc"), level_ucs, str(out_path / "provenance"), level_labels)
    typer.echo(f"summary\t{summary}")
    success(f"AutoTax2 run finished: {out_path}")


def main() -> None:
    app()


if __name__ == "__main__":
    main()
