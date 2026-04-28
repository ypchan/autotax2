from __future__ import annotations
from dataclasses import dataclass
from typing import List

@dataclass(frozen=True)
class RankThreshold:
    rank: str
    threshold: float
    label: str

DEFAULT_RANK_THRESHOLDS = [
    RankThreshold('species', 0.987, 'species-like'),
    RankThreshold('genus', 0.945, 'genus-like'),
    RankThreshold('family', 0.865, 'family-like'),
    RankThreshold('order', 0.820, 'order-like'),
    RankThreshold('class', 0.785, 'class-like'),
    RankThreshold('phylum', 0.750, 'phylum-like'),
]

RANK_ORDER = ['phylum', 'class', 'order', 'family', 'genus', 'species']
RANK_TO_TAX_COL = {'phylum': 'Phylum', 'class': 'Class', 'order': 'Order', 'family': 'Family', 'genus': 'Genus', 'species': 'Species'}
RANK_TO_PARENT = {'species': 'genus', 'genus': 'family', 'family': 'order', 'order': 'class', 'class': 'phylum', 'phylum': 'domain'}

def parse_rank_thresholds(value: str | None) -> List[RankThreshold]:
    if value is None or value.strip().lower() == 'default':
        return list(DEFAULT_RANK_THRESHOLDS)
    defaults = {x.rank: x.threshold for x in DEFAULT_RANK_THRESHOLDS}
    out = []
    for item in value.split(','):
        item = item.strip()
        if not item:
            continue
        if ':' in item:
            rank, threshold = item.split(':', 1)
            rank = rank.strip().lower()
            threshold_f = float(threshold)
        else:
            rank = item.lower()
            if rank not in defaults:
                raise ValueError(f'Unknown rank: {rank}. Use one of {sorted(defaults)} or rank:threshold.')
            threshold_f = defaults[rank]
        out.append(RankThreshold(rank, threshold_f, f'{rank}-like'))
    return out

def rank_taxon_id(rank: str, taxonomy: dict[str, str], prefix: str = 'silva') -> str:
    rank = rank.lower()
    if rank not in RANK_TO_TAX_COL:
        raise ValueError(f'Unsupported rank: {rank}')
    parts = []
    for r in RANK_ORDER:
        col = RANK_TO_TAX_COL[r]
        if taxonomy.get(col):
            parts.append(f'{r}:{taxonomy[col]}')
        if r == rank:
            break
    return f'{prefix}|{rank}|' + ';'.join(parts) if parts else ''

def parent_taxon_id(rank: str, taxonomy: dict[str, str], prefix: str = 'silva') -> str:
    parent = RANK_TO_PARENT.get(rank.lower())
    if not parent or parent == 'domain':
        domain = taxonomy.get('Domain') or taxonomy.get('Kingdom') or ''
        return f'{prefix}|domain|d:{domain}' if domain else ''
    return rank_taxon_id(parent, taxonomy, prefix=prefix)
