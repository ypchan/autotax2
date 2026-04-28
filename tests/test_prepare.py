from autotax2.prepare import clean_taxon, make_sintax_header, parse_silva_header


def test_clean_taxon():
    assert clean_taxon("uncultured bacterium") is None
    assert clean_taxon("Firmicutes") == "Firmicutes"
    assert clean_taxon("A,B") == "A.B"


def test_parse_header():
    seq_id, taxa = parse_silva_header("AB1 Bacteria;Firmicutes;Bacilli;")
    assert seq_id == "AB1"
    assert taxa[0] == "Bacteria"
    assert taxa[1] == "Firmicutes"


def test_sintax_header():
    h = make_sintax_header("AB1", ["Bacteria", "Firmicutes", None, None, None, None, None])
    assert h == "AB1;tax=d:Bacteria,p:Firmicutes;"
