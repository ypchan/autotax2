from __future__ import annotations
import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

from .assign import parse_vsearch_uc_for_new_clusters, relabel_fasta
from .logging import info, success
from .ranks import parent_taxon_id, parse_rank_thresholds, rank_taxon_id
from .utils import ensure_dir, ensure_file, parse_manifest, read_fasta, write_fasta
from .vsearch import cluster as vsearch_cluster
from .vsearch import usearch_global


def read_silva_taxonomy(path: str | Path) -> Dict[str, Dict[str, str]]:
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        raise FileNotFoundError(f'SILVA taxonomy table not found: {p}')
    out: Dict[str, Dict[str, str]] = {}
    with open(p, encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sid = row.get('ID', '')
            if sid:
                out[sid] = {k: v for k, v in row.items() if k != 'ID' and v}
    return out


def read_version_map(path: str | Path | None) -> Dict[str, Dict[str, str]]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return {}
    out: Dict[str, Dict[str, str]] = {}
    with open(p, encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            aid = row.get('analysis_sequence_id') or row.get('sequence_id')
            if aid:
                out[aid] = row
    return out


def fasta_dict(path: str | Path) -> Dict[str, Tuple[str, str]]:
    return {h.split()[0].split(';')[0]: (h, s) for h, s in read_fasta(path)}


def parse_hit_table(userout_tsv: str | Path) -> Dict[str, Dict[str, str]]:
    fields = ['query', 'target', 'id', 'alnlen', 'mism', 'opens', 'qlo', 'qhi', 'tlo', 'thi', 'qcov', 'tcov']
    hits: Dict[str, Dict[str, str]] = {}
    p = Path(userout_tsv)
    if not p.exists() or p.stat().st_size == 0:
        return hits
    with open(p, encoding='utf-8') as f:
        for line in f:
            if not line.strip():
                continue
            row = dict(zip(fields, line.rstrip('\n').split('\t')))
            q = row.get('query', '')
            if q and q not in hits:
                hits[q] = row
    return hits


def write_rows(path: str | Path, rows: List[Dict[str, str]], fields: List[str]) -> None:
    with open(path, 'w', encoding='utf-8', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, '') for k in fields})


def choose_db(manifest: Dict[str, str], key_fasta: str, key_udb: str, db_format: str) -> str:
    if db_format == 'udb':
        return manifest.get(key_udb) or manifest[key_fasta]
    if db_format == 'fasta':
        return manifest[key_fasta]
    return manifest.get(key_udb) or manifest[key_fasta]


def make_official_centroids(assignment_rows: List[Dict[str, str]], original_records: Dict[str, Tuple[str, str]], analysis_records: Dict[str, Tuple[str, str]], out_fasta: str | Path) -> None:
    taxon_to_seq: Dict[str, str] = {}
    for row in assignment_rows:
        taxon = row['taxon_id']
        if taxon not in taxon_to_seq:
            taxon_to_seq[taxon] = row.get('original_sequence_id') or row['sequence_id']
    records = []
    for taxon, oid in taxon_to_seq.items():
        if oid in original_records:
            records.append((taxon, original_records[oid][1]))
        elif oid in analysis_records:
            records.append((taxon, analysis_records[oid][1]))
    write_fasta(records, out_fasta)


def insert_backbone(input_fasta: str | Path, manifest_path: str | Path, outdir: str | Path, source_label: str, rank_thresholds: str | None = 'default', original_fasta: str | Path | None = None, version_map_path: str | Path | None = None, db_format: str = 'auto', vsearch: str = 'vsearch', threads: int = 8, strand: str = 'both', maxaccepts: int = 1, dry_run: bool = False) -> Dict[str, str]:
    input_fasta = ensure_file(input_fasta, 'input FASTA')
    outdir = ensure_dir(outdir)
    manifest = parse_manifest(manifest_path)
    taxonomy = read_silva_taxonomy(manifest['silva_taxonomy_tsv'])
    ranks = parse_rank_thresholds(rank_thresholds)
    db = choose_db(manifest, 'silva_fasta', 'silva_udb', db_format)

    rank_uc_dir = ensure_dir(outdir / 'rank_uc')
    rank_centroids_core = ensure_dir(outdir / 'rank_centroids_core')
    rank_centroids_original = ensure_dir(outdir / 'rank_centroids_original')
    rank_hits = ensure_dir(outdir / 'rank_hits')
    rank_unmatched = ensure_dir(outdir / 'rank_unmatched')

    analysis_records = fasta_dict(input_fasta)
    original_records = fasta_dict(original_fasta) if original_fasta else analysis_records
    version_map = read_version_map(version_map_path)

    all_assignment_rows: List[Dict[str, str]] = []
    all_summary_rows: List[Dict[str, str]] = []

    for rt in ranks:
        rank = rt.rank
        threshold = rt.threshold
        info(f'Backbone insertion for rank={rank}, threshold={threshold}')
        search_out = usearch_global(str(input_fasta), db, str(rank_hits), f'{rank}_vs_silva', threshold, maxaccepts, strand, vsearch, threads, True, dry_run)
        if dry_run:
            continue
        notmatched = Path(search_out['notmatched'])
        rank_notmatched = rank_unmatched / f'{rank}.unmatched.fa'
        if notmatched.exists():
            notmatched.replace(rank_notmatched)
        else:
            rank_notmatched.write_text('', encoding='utf-8')

        hits = parse_hit_table(search_out['userout'])
        assignment_rows: List[Dict[str, str]] = []
        for seq_id, hit in hits.items():
            target_id = hit.get('target', '').split()[0].split(';')[0]
            tax = taxonomy.get(target_id, {})
            taxon = rank_taxon_id(rank, tax, prefix='silva') or f'silva|{rank}|target:{target_id}'
            parent = parent_taxon_id(rank, tax, prefix='silva')
            vm = version_map.get(seq_id, {})
            original_id = vm.get('sequence_id') or seq_id
            assignment_rows.append({'sequence_id': seq_id, 'analysis_sequence_id': seq_id, 'original_sequence_id': original_id, 'source_label': source_label, 'rank': rank, 'rank_label': rt.label, 'identity_threshold': str(threshold), 'assignment_type': 'silva_existing', 'taxon_id': taxon, 'parent_taxon_id': parent, 'cluster_id': taxon, 'silva_target': target_id, 'silva_identity': hit.get('id', ''), 'silva_qcov': hit.get('qcov', ''), 'silva_tcov': hit.get('tcov', ''), 'is_novel': 'no', 'has_intron': vm.get('has_intron', 'unknown' if version_map else 'no')})

        if rank_notmatched.exists() and rank_notmatched.stat().st_size > 0:
            cluster_dir = ensure_dir(outdir / 'rank_novel_clusters' / rank)
            result = vsearch_cluster(str(rank_notmatched), str(cluster_dir), threshold, 'cluster_size', vsearch, threads, None, False)
            raw_centroids = Path(result['centroids'])
            raw_uc = Path(result['uc'])
            novel_uc = rank_uc_dir / f'{rank}.novel.uc'
            raw_uc.replace(novel_uc)
            core_centroids = rank_centroids_core / f'{rank}.centroids.core.fa'
            prefix = f'novel|{rank}|{source_label}_{rank}_'
            old_to_new = relabel_fasta(raw_centroids, core_centroids, prefix)
            parsed = parse_vsearch_uc_for_new_clusters(novel_uc, old_to_new)
            for r in parsed:
                seq_id = r['query']
                taxon = r['cluster_id']
                vm = version_map.get(seq_id, {})
                original_id = vm.get('sequence_id') or seq_id
                assignment_rows.append({'sequence_id': seq_id, 'analysis_sequence_id': seq_id, 'original_sequence_id': original_id, 'source_label': source_label, 'rank': rank, 'rank_label': rt.label, 'identity_threshold': str(threshold), 'assignment_type': r['assignment_type'], 'taxon_id': taxon, 'parent_taxon_id': '', 'cluster_id': taxon, 'silva_target': '', 'silva_identity': '', 'silva_qcov': '', 'silva_tcov': '', 'is_novel': 'yes', 'has_intron': vm.get('has_intron', 'unknown' if version_map else 'no')})
        else:
            (rank_uc_dir / f'{rank}.novel.uc').write_text('', encoding='utf-8')
            (rank_centroids_core / f'{rank}.centroids.core.fa').write_text('', encoding='utf-8')

        Path(search_out['uc']).replace(rank_uc_dir / f'{rank}.silva.uc')
        make_official_centroids(assignment_rows, original_records, analysis_records, rank_centroids_original / f'{rank}.centroids.original.fa')

        counts = Counter(row['taxon_id'] for row in assignment_rows)
        types = defaultdict(Counter)
        novel = defaultdict(int)
        for row in assignment_rows:
            types[row['taxon_id']][row['assignment_type']] += 1
            if row['is_novel'] == 'yes':
                novel[row['taxon_id']] += 1
        for taxon, n in sorted(counts.items()):
            all_summary_rows.append({'rank': rank, 'rank_label': rt.label, 'identity_threshold': str(threshold), 'taxon_id': taxon, 'source_label': source_label, 'n_sequences': str(n), 'is_novel': 'yes' if novel[taxon] else 'no', 'assignment_types_json': str(dict(types[taxon]))})
        all_assignment_rows.extend(assignment_rows)

    assignment_path = outdir / 'sequence_rank_assignment.tsv'
    summary_path = outdir / 'rank_taxa_summary.tsv'
    assignment_fields = ['sequence_id', 'analysis_sequence_id', 'original_sequence_id', 'source_label', 'rank', 'rank_label', 'identity_threshold', 'assignment_type', 'taxon_id', 'parent_taxon_id', 'cluster_id', 'silva_target', 'silva_identity', 'silva_qcov', 'silva_tcov', 'is_novel', 'has_intron']
    summary_fields = ['rank', 'rank_label', 'identity_threshold', 'taxon_id', 'source_label', 'n_sequences', 'is_novel', 'assignment_types_json']
    if not dry_run:
        write_rows(assignment_path, all_assignment_rows, assignment_fields)
        write_rows(summary_path, all_summary_rows, summary_fields)
    success(f'Backbone insertion finished: {outdir}')
    return {'sequence_rank_assignment': str(assignment_path), 'rank_taxa_summary': str(summary_path), 'rank_uc_dir': str(rank_uc_dir), 'rank_centroids_core': str(rank_centroids_core), 'rank_centroids_original': str(rank_centroids_original)}
