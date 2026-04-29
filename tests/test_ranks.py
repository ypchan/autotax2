from src.ranks import parse_rank_thresholds, rank_taxon_id

def test_default_thresholds():
    r = parse_rank_thresholds('default')
    assert r[0].rank == 'species'
    assert abs(r[0].threshold - 0.987) < 1e-9

def test_rank_taxon_id():
    tax = {'Phylum':'Firmicutes','Class':'Bacilli','Order':'Lactobacillales','Family':'Lactobacillaceae','Genus':'Lactobacillus'}
    tid = rank_taxon_id('genus', tax)
    assert tid.startswith('silva|genus|')
    assert 'genus:Lactobacillus' in tid
