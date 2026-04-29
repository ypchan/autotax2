from src.intron import find_query_insertions, identity_excluding_insertions, remove_intervals


def test_find_query_insertions():
    qrow = "AAAAACCCCCGGGGG"
    trow = "AAAAA-----GGGGG"
    insertions = find_query_insertions(qrow, trow, min_intron_len=5)
    assert len(insertions) == 1
    assert insertions[0]["query_start"] == 6
    assert insertions[0]["query_end"] == 10


def test_remove_intervals():
    assert remove_intervals("AAAAACCCCCGGGGG", [(6, 10)]) == "AAAAAGGGGG"


def test_identity_excluding_insertions():
    qrow = "AAAAACCCCCGGGGG"
    trow = "AAAAA-----GGGGA"
    insertions = find_query_insertions(qrow, trow, min_intron_len=5)
    ident = identity_excluding_insertions(qrow, trow, insertions)
    assert 0.8 < ident < 1.0
