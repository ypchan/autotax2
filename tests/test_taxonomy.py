from autotax2.taxonomy import infer_anchor_rank, make_placeholder, parse_tax_string


def test_parse_tax_string():
    tax = parse_tax_string("d__Bacteria;p__P;c__C;o__O;f__F;g__G;s__S;")
    assert tax["domain"] == "d__Bacteria"
    assert tax["genus"] == "g__G"


def test_infer_anchor_rank():
    tax = parse_tax_string("d__Bacteria;p__P;c__C;o__O;f__F;g__G;s__S;")
    thresholds = {"species": 97.2, "genus": 90.1, "family": 80.1, "order": 72.9, "class": 72.2, "phylum": 69.6}
    assert infer_anchor_rank(98.0, tax, thresholds) == "species"
    assert infer_anchor_rank(88.0, tax, thresholds) == "family"
    assert infer_anchor_rank(60.0, tax, thresholds) is None


def test_placeholder():
    assert make_placeholder("species", "midas", 101) == "s__midas_s000101"
    assert make_placeholder("genus", "mfd", 2) == "g__mfd_g000002"
