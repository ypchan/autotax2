from __future__ import annotations
import csv
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Set
from .utils import ensure_dir

def read_assignment(path: str | Path, label: str | None = None) -> List[Dict[str, str]]:
    rows = []
    with open(path, encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if label:
                row['source_label'] = label
            rows.append(row)
    return rows

def write_tsv(path: str | Path, rows: List[Dict[str, str]], fields: List[str]) -> None:
    with open(path, 'w', encoding='utf-8', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, '') for k in fields})

def overlap_backbone(assignments: List[str | Path], labels: List[str], outdir: str | Path) -> Dict[str, str]:
    if len(assignments) != len(labels):
        raise ValueError('--assignments and --labels must have the same length')
    outdir = ensure_dir(outdir)
    rows = []
    for path, label in zip(assignments, labels):
        rows.extend(read_assignment(path, label=label))
    taxon_source_counts: Dict[tuple[str, str], Counter] = defaultdict(Counter)
    for row in rows:
        taxon_source_counts[(row['rank'], row['taxon_id'])][row.get('source_label', '')] += 1
    presence_rows = []
    count_rows = []
    for (rank, taxon), counts in sorted(taxon_source_counts.items()):
        source_set = sorted([s for s, n in counts.items() if n > 0])
        pres = {'rank': rank, 'taxon_id': taxon, 'sources': ','.join(source_set), 'n_sources': str(len(source_set)), 'is_source_specific': 'yes' if len(source_set) == 1 else 'no', 'source_specific_to': source_set[0] if len(source_set) == 1 else ''}
        for label in labels:
            pres[label] = 'yes' if counts.get(label, 0) > 0 else 'no'
            count_rows.append({'rank': rank, 'taxon_id': taxon, 'source_label': label, 'n_sequences': str(counts.get(label, 0)), 'present': 'yes' if counts.get(label, 0) > 0 else 'no'})
        presence_rows.append(pres)
    pairwise_rows = []
    for rank in sorted(set(r['rank'] for r in rows)):
        source_taxa: Dict[str, Set[str]] = defaultdict(set)
        for (r, taxon), counts in taxon_source_counts.items():
            if r != rank:
                continue
            for source, n in counts.items():
                if n > 0:
                    source_taxa[source].add(taxon)
        for a, b in combinations(labels, 2):
            A, B = source_taxa[a], source_taxa[b]
            shared, union = A & B, A | B
            pairwise_rows.append({'rank': rank, 'source_a': a, 'source_b': b, 'shared_taxon_count': str(len(shared)), 'source_a_taxon_count': str(len(A)), 'source_b_taxon_count': str(len(B)), 'union_taxon_count': str(len(union)), 'jaccard_taxon_overlap': f'{len(shared)/len(union):.6f}' if union else '0'})
    unique_rows = [{'rank': r['rank'], 'taxon_id': r['taxon_id'], 'source_label': r['source_specific_to']} for r in presence_rows if r['is_source_specific'] == 'yes']
    presence_path = outdir/'taxon_presence_by_source.tsv'
    count_path = outdir/'taxon_count_by_source.tsv'
    pairwise_path = outdir/'source_pairwise_overlap_by_rank.tsv'
    unique_path = outdir/'source_unique_taxa_by_rank.tsv'
    write_tsv(presence_path, presence_rows, ['rank','taxon_id','sources','n_sources','is_source_specific','source_specific_to'] + labels)
    write_tsv(count_path, count_rows, ['rank','taxon_id','source_label','n_sequences','present'])
    write_tsv(pairwise_path, pairwise_rows, ['rank','source_a','source_b','shared_taxon_count','source_a_taxon_count','source_b_taxon_count','union_taxon_count','jaccard_taxon_overlap'])
    write_tsv(unique_path, unique_rows, ['rank','taxon_id','source_label'])
    return {'taxon_presence_by_source': str(presence_path), 'taxon_count_by_source': str(count_path), 'source_pairwise_overlap_by_rank': str(pairwise_path), 'source_unique_taxa_by_rank': str(unique_path)}
