from autotax2.sina import parse_sina_header


def test_parse_sina_header():
    h = "seq1 [align_ident_slv=88.5] [align_quality_slv=90] [lca_tax_gtdb=d__Bacteria;p__P;] [turn=reversed and complemented]"
    ann = parse_sina_header(h)
    assert ann.seq_id == "seq1"
    assert ann.align_ident == 88.5
    assert ann.align_quality == 90
    assert ann.turn == "reversed and complemented"
    assert ann.lca_tax_ref == "d__Bacteria;p__P;"
