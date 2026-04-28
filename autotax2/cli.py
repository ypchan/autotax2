from __future__ import annotations

import argparse
from pathlib import Path

try:
    from rich_argparse import RichHelpFormatter
except Exception:  # pragma: no cover
    RichHelpFormatter = argparse.HelpFormatter  # type: ignore

from .assign import assign_or_create
from .backbone import insert_backbone
from .dependencies import check_core_dependencies, check_paths, summarize_checks
from .intron import detect_introns
from .logging import info, print_example, progress_context, progress_step, success
from .overlap import overlap_backbone
from .prepare import prepare_silva
from .provenance import annotate_assignments_with_source, provenance_from_uc_levels
from .summarize import summarize
from .threads import resolve_threads
from .utils import ensure_dir, parse_manifest
from .vsearch import cluster as run_cluster
from .vsearch import dereplicate, sintax, sort_by_size, usearch_global


EXAMPLES = {
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
  --assignments ref2/sequence_rank_assignment.tsv ref3/sequence_rank_assignment.tsv hifimeta/sequence_rank_assignment.tsv \\
  --labels ref2,ref3,hifimeta \\
  --out ref_overlap""",
    "check": "autotax2 check --ref-manifest refdatabases/autotax2_ref_manifest.tsv",
    "prepare-silva": """autotax2 prepare-silva \\
  --silva-fasta SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \\
  --silva-metadata SILVA_138.2_SSURef.full_metadata.gz \\
  --out refdatabases \\
  --make-udb \\
  --threads auto""",
    "derep": "autotax2 derep --input input.fa --out work --sort --threads auto",
    "cluster": "autotax2 cluster --input work/derep_sorted.fa --out clusters --ids 0.99,0.97 --threads auto",
    "classify": "autotax2 classify --input clusters/otu099_centroids.fa --ref-manifest refdatabases/autotax2_ref_manifest.tsv --out classify",
    "assign": "autotax2 assign --new new.fa --old-centroids centroids_97.fa --id 0.97 --out assign_97 --threads auto",
    "provenance": "autotax2 provenance --source-map source_map.tsv --derep-uc work/derep.uc --level-uc clusters/otu099.uc --level-labels 0.99 --out provenance",
    "summarize": "autotax2 summarize --sintax classify/sintax.tsv --silva-hits classify/silva_hits.tsv --typestrain-hits classify/typestrain_hits.tsv --out classify/taxonomy_summary.tsv",
    "run": "autotax2 run --input input.fa --ref-manifest refdatabases/autotax2_ref_manifest.tsv --out autotax2_out --threads auto",
}


def add_common_vsearch_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--vsearch", default="vsearch")
    p.add_argument("--threads", default="auto")
    p.add_argument("--dry-run", action="store_true")


def add_example_arg(p: argparse.ArgumentParser) -> None:
    p.add_argument("--example", action="store_true")


def maybe_print_example(args: argparse.Namespace, name: str) -> bool:
    if getattr(args, "example", False):
        print_example(f"Example: autotax2 {name}", EXAMPLES[name])
        return True
    return False


def cmd_check(args: argparse.Namespace) -> None:
    if maybe_print_example(args, "check"):
        return
    ok = check_core_dependencies(args.vsearch)
    extra = []
    if args.ref_manifest:
        manifest = parse_manifest(args.ref_manifest)
        extra.extend([
            ("silva_fasta", manifest.get("silva_fasta", ""), True),
            ("silva_sintax_fasta", manifest.get("silva_sintax_fasta", ""), True),
            ("silva_typestrains_fasta", manifest.get("silva_typestrains_fasta", ""), True),
        ])
    if args.source_map:
        extra.append(("source_map", args.source_map, True))
    if extra:
        ok = summarize_checks(check_paths(extra)) and ok
    if not ok:
        raise SystemExit(1)
    success("All required dependency checks passed.")


def cmd_detect_intron(args: argparse.Namespace) -> None:
    if maybe_print_example(args, "detect-intron"):
        return
    outputs = detect_introns(
        input_fasta=args.input,
        db=args.db,
        outdir=args.out,
        source_label=args.source_label,
        vsearch=args.vsearch,
        threads=resolve_threads(args.threads),
        search_id=args.search_id,
        rescue_id=args.rescue_id,
        min_intron_len=args.min_intron_len,
        min_flank_len=args.min_flank_len,
        strand=args.strand,
        maxaccepts=args.maxaccepts,
        dry_run=args.dry_run,
    )
    for k, v in outputs.items():
        print(f"{k}\t{v}")


def cmd_insert_backbone(args: argparse.Namespace) -> None:
    if maybe_print_example(args, "insert-backbone"):
        return
    outputs = insert_backbone(
        input_fasta=args.input,
        manifest_path=args.silva_manifest,
        outdir=args.out,
        source_label=args.source_label,
        rank_thresholds=args.rank_thresholds,
        original_fasta=args.original_fasta,
        version_map_path=args.version_map,
        db_format=args.db_format,
        vsearch=args.vsearch,
        threads=resolve_threads(args.threads),
        strand=args.strand,
        maxaccepts=args.maxaccepts,
        centroid_policy=args.centroid_policy,
        dry_run=args.dry_run,
    )
    for k, v in outputs.items():
        print(f"{k}\t{v}")


def cmd_overlap_backbone(args: argparse.Namespace) -> None:
    if maybe_print_example(args, "overlap-backbone"):
        return
    labels = [x.strip() for x in args.labels.split(",") if x.strip()]
    outputs = overlap_backbone(args.assignments, labels, args.out)
    for k, v in outputs.items():
        print(f"{k}\t{v}")


def cmd_prepare_silva(args: argparse.Namespace) -> None:
    if maybe_print_example(args, "prepare-silva"):
        return
    entries = prepare_silva(
        silva_fasta=args.silva_fasta,
        silva_metadata=args.silva_metadata,
        outdir=args.out,
        prefix=args.prefix,
        make_udb_files=args.make_udb,
        vsearch=args.vsearch,
        dry_run=args.dry_run,
    )
    for k, v in entries.items():
        print(f"{k}\t{v}")


def cmd_derep(args: argparse.Namespace) -> None:
    if maybe_print_example(args, "derep"):
        return
    threads = resolve_threads(args.threads)
    derep = dereplicate(args.input, args.out, args.vsearch, threads, args.dry_run)
    if args.sort:
        print(sort_by_size(derep, args.out, args.vsearch, args.minsize, args.dry_run))
    else:
        print(derep)


def parse_id_list(value: str) -> list[float]:
    return [float(x.strip()) for x in value.split(",") if x.strip()]


def cmd_cluster(args: argparse.Namespace) -> None:
    if maybe_print_example(args, "cluster"):
        return
    threads = resolve_threads(args.threads)
    for identity in parse_id_list(args.ids):
        result = run_cluster(args.input, args.out, identity, args.method, args.vsearch, threads, args.relabel, args.dry_run)
        print(f"{identity}\t{result['centroids']}\t{result['uc']}")


def cmd_classify(args: argparse.Namespace) -> None:
    if maybe_print_example(args, "classify"):
        return
    threads = resolve_threads(args.threads)
    manifest = parse_manifest(args.ref_manifest) if args.ref_manifest else {}
    silva_fasta = args.silva_fasta or manifest.get("silva_fasta")
    silva_sintax = args.silva_sintax or manifest.get("silva_sintax_fasta")
    typestrains = args.typestrains or manifest.get("silva_typestrains_fasta")
    if not silva_fasta or not silva_sintax or not typestrains:
        raise SystemExit("Missing reference files. Use --ref-manifest or explicit reference paths.")
    ensure_dir(args.out)
    sintax_tsv = sintax(args.input, silva_sintax, args.out, args.sintax_cutoff, args.vsearch, threads, args.dry_run)
    silva = usearch_global(args.input, silva_fasta, args.out, "silva_hits", args.min_id, args.maxaccepts, args.strand, args.vsearch, threads, False, args.dry_run)
    typehit = usearch_global(args.input, typestrains, args.out, "typestrain_hits", args.min_id, args.maxaccepts, args.strand, args.vsearch, threads, False, args.dry_run)
    summary = str(Path(args.out) / "taxonomy_summary.tsv")
    if not args.dry_run:
        summarize(sintax_tsv, silva["userout"], typehit["userout"], summary)
    print(f"summary\t{summary}")


def cmd_assign(args: argparse.Namespace) -> None:
    if maybe_print_example(args, "assign"):
        return
    outputs = assign_or_create(args.new, args.old_centroids, args.out, args.id, args.vsearch, resolve_threads(args.threads), args.strand, args.new_cluster_prefix, args.method, args.maxaccepts, args.dry_run)
    if args.source_map and not args.dry_run:
        annotated = str(Path(args.out) / "assignments_with_source.tsv")
        annotate_assignments_with_source(outputs["assignments"], args.source_map, annotated)
        outputs["assignments_with_source"] = annotated
    for k, v in outputs.items():
        print(f"{k}\t{v}")


def cmd_provenance(args: argparse.Namespace) -> None:
    if maybe_print_example(args, "provenance"):
        return
    labels = [x.strip() for x in args.level_labels.split(",")] if args.level_labels else None
    outputs = provenance_from_uc_levels(args.source_map, args.derep_uc, args.level_uc, args.out, labels)
    for k, v in outputs.items():
        print(f"{k}\t{v}")


def cmd_summarize(args: argparse.Namespace) -> None:
    if maybe_print_example(args, "summarize"):
        return
    summarize(args.sintax, args.silva_hits, args.typestrain_hits, args.out)
    print(args.out)


def cmd_run(args: argparse.Namespace) -> None:
    if maybe_print_example(args, "run"):
        return
    threads = resolve_threads(args.threads)
    out = Path(args.out)
    work = out / "work"
    clusters = out / "clusters"
    classify_dir = out / "classify"
    ensure_dir(work); ensure_dir(clusters); ensure_dir(classify_dir)
    derep = dereplicate(args.input, str(work), args.vsearch, threads, args.dry_run)
    sorted_fa = sort_by_size(derep, str(work), args.vsearch, args.minsize, args.dry_run)
    chosen = None
    for identity in parse_id_list(args.ids):
        result = run_cluster(sorted_fa, str(clusters), identity, args.method, args.vsearch, threads, None, args.dry_run)
        if abs(identity - args.classify_id) < 1e-9:
            chosen = result["centroids"]
    chosen = chosen or str(clusters / f"otu{str(parse_id_list(args.ids)[0]).replace('.', '')}_centroids.fa")
    class_args = argparse.Namespace(input=chosen, out=str(classify_dir), ref_manifest=args.ref_manifest, silva_fasta=None, silva_sintax=None, typestrains=None, sintax_cutoff=args.sintax_cutoff, min_id=args.min_id, maxaccepts=args.maxaccepts, strand=args.strand, vsearch=args.vsearch, threads=str(threads), dry_run=args.dry_run, example=False)
    cmd_classify(class_args)
    success(f"AutoTax2 run finished: {out}")


def require_args(args: argparse.Namespace, names: list[str]) -> None:
    if getattr(args, "example", False):
        return
    missing = [n for n in names if getattr(args, n.replace("-", "_"), None) in (None, "")]
    if missing:
        raise SystemExit(f"Missing required arguments: {', '.join('--' + x for x in missing)}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="autotax2", description="AutoTax2 workflow.", formatter_class=RichHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    pcx = sub.add_parser("check", formatter_class=RichHelpFormatter)
    pcx.add_argument("--vsearch", default="vsearch"); pcx.add_argument("--ref-manifest", default=None); pcx.add_argument("--source-map", default=None); add_example_arg(pcx); pcx.set_defaults(func=cmd_check)

    pdi = sub.add_parser("detect-intron", formatter_class=RichHelpFormatter)
    pdi.add_argument("--input", required=False); pdi.add_argument("--db", required=False); pdi.add_argument("--out", required=False)
    pdi.add_argument("--source-label", default="query"); pdi.add_argument("--search-id", type=float, default=0.70); pdi.add_argument("--rescue-id", type=float, default=0.987)
    pdi.add_argument("--min-intron-len", type=int, default=50); pdi.add_argument("--min-flank-len", type=int, default=150); pdi.add_argument("--maxaccepts", type=int, default=1); pdi.add_argument("--strand", default="both", choices=["both", "plus"])
    add_common_vsearch_args(pdi); add_example_arg(pdi); pdi.set_defaults(func=cmd_detect_intron)

    pib = sub.add_parser("insert-backbone", formatter_class=RichHelpFormatter)
    pib.add_argument("--input", required=False); pib.add_argument("--original-fasta", default=None); pib.add_argument("--version-map", default=None); pib.add_argument("--source-label", required=False); pib.add_argument("--silva-manifest", required=False); pib.add_argument("--rank-thresholds", default="default"); pib.add_argument("--db-format", default="auto", choices=["auto", "udb", "fasta"]); pib.add_argument("--maxaccepts", type=int, default=1); pib.add_argument("--strand", default="both", choices=["both", "plus"]); pib.add_argument("--centroid-policy", default="original_seed", choices=["original_seed", "prefer_non_intron"]); pib.add_argument("--out", required=False)
    add_common_vsearch_args(pib); add_example_arg(pib); pib.set_defaults(func=cmd_insert_backbone)

    pob = sub.add_parser("overlap-backbone", formatter_class=RichHelpFormatter)
    pob.add_argument("--assignments", nargs="+", required=False); pob.add_argument("--labels", required=False); pob.add_argument("--out", required=False); add_example_arg(pob); pob.set_defaults(func=cmd_overlap_backbone)

    ps = sub.add_parser("prepare-silva", formatter_class=RichHelpFormatter)
    ps.add_argument("--silva-fasta", required=False); ps.add_argument("--silva-metadata", required=False); ps.add_argument("--out", default="refdatabases"); ps.add_argument("--prefix", default=None); ps.add_argument("--make-udb", action="store_true"); add_common_vsearch_args(ps); add_example_arg(ps); ps.set_defaults(func=cmd_prepare_silva)

    pd = sub.add_parser("derep", formatter_class=RichHelpFormatter)
    pd.add_argument("--input", required=False); pd.add_argument("--out", required=False); pd.add_argument("--sort", action="store_true"); pd.add_argument("--minsize", type=int, default=None); add_common_vsearch_args(pd); add_example_arg(pd); pd.set_defaults(func=cmd_derep)

    pc = sub.add_parser("cluster", formatter_class=RichHelpFormatter)
    pc.add_argument("--input", required=False); pc.add_argument("--out", required=False); pc.add_argument("--ids", default="0.99,0.97,0.95,0.90"); pc.add_argument("--method", default="cluster_size", choices=["cluster_size", "cluster_fast", "cluster_smallmem"]); pc.add_argument("--relabel", default=None); add_common_vsearch_args(pc); add_example_arg(pc); pc.set_defaults(func=cmd_cluster)

    pcl = sub.add_parser("classify", formatter_class=RichHelpFormatter)
    pcl.add_argument("--input", required=False); pcl.add_argument("--out", required=False); pcl.add_argument("--ref-manifest", default=None); pcl.add_argument("--silva-fasta", default=None); pcl.add_argument("--silva-sintax", default=None); pcl.add_argument("--typestrains", default=None); pcl.add_argument("--sintax-cutoff", type=float, default=0.8); pcl.add_argument("--min-id", type=float, default=0.70); pcl.add_argument("--maxaccepts", type=int, default=10); pcl.add_argument("--strand", default="both", choices=["both", "plus"]); add_common_vsearch_args(pcl); add_example_arg(pcl); pcl.set_defaults(func=cmd_classify)

    pa = sub.add_parser("assign", formatter_class=RichHelpFormatter)
    pa.add_argument("--new", required=False); pa.add_argument("--old-centroids", required=False); pa.add_argument("--out", required=False); pa.add_argument("--id", type=float, required=False); pa.add_argument("--method", default="cluster_size", choices=["cluster_size", "cluster_fast", "cluster_smallmem"]); pa.add_argument("--new-cluster-prefix", default="new_"); pa.add_argument("--source-map", default=None); pa.add_argument("--maxaccepts", type=int, default=1); pa.add_argument("--strand", default="both", choices=["both", "plus"]); add_common_vsearch_args(pa); add_example_arg(pa); pa.set_defaults(func=cmd_assign)

    pp = sub.add_parser("provenance", formatter_class=RichHelpFormatter)
    pp.add_argument("--source-map", required=False); pp.add_argument("--derep-uc", required=False); pp.add_argument("--level-uc", required=False, nargs="+"); pp.add_argument("--level-labels", default=None); pp.add_argument("--out", required=False); add_example_arg(pp); pp.set_defaults(func=cmd_provenance)

    psu = sub.add_parser("summarize", formatter_class=RichHelpFormatter)
    psu.add_argument("--sintax", required=False); psu.add_argument("--silva-hits", required=False); psu.add_argument("--typestrain-hits", required=False); psu.add_argument("--out", required=False); add_example_arg(psu); psu.set_defaults(func=cmd_summarize)

    pr = sub.add_parser("run", formatter_class=RichHelpFormatter)
    pr.add_argument("--input", required=False); pr.add_argument("--ref-manifest", required=False); pr.add_argument("--source-map", default=None); pr.add_argument("--out", required=False); pr.add_argument("--ids", default="0.99,0.97,0.95,0.90"); pr.add_argument("--classify-id", type=float, default=0.99); pr.add_argument("--method", default="cluster_size", choices=["cluster_size", "cluster_fast", "cluster_smallmem"]); pr.add_argument("--minsize", type=int, default=None); pr.add_argument("--sintax-cutoff", type=float, default=0.8); pr.add_argument("--min-id", type=float, default=0.70); pr.add_argument("--maxaccepts", type=int, default=10); pr.add_argument("--strand", default="both", choices=["both", "plus"]); add_common_vsearch_args(pr); add_example_arg(pr); pr.set_defaults(func=cmd_run)

    return p


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    requirements = {
        "detect-intron": ["input", "db", "out"],
        "insert-backbone": ["input", "source-label", "silva-manifest", "out"],
        "overlap-backbone": ["assignments", "labels", "out"],
        "prepare-silva": ["silva-fasta", "silva-metadata"],
        "derep": ["input", "out"],
        "cluster": ["input", "out"],
        "classify": ["input", "out"],
        "assign": ["new", "old-centroids", "out", "id"],
        "provenance": ["source-map", "derep-uc", "level-uc", "out"],
        "summarize": ["sintax", "silva-hits", "typestrain-hits", "out"],
        "run": ["input", "ref-manifest", "out"],
    }
    if args.cmd in requirements:
        require_args(args, requirements[args.cmd])
    args.func(args)


if __name__ == "__main__":
    main()
