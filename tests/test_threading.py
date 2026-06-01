import pandas as pd

from autotax2.core import build_cluster_groups, resolve_group_parallelism


def test_resolve_group_parallelism_respects_thread_budget():
    assert resolve_group_parallelism(group_count=10, threads=8, group_jobs=4) == (4, 2)
    assert resolve_group_parallelism(group_count=10, threads=8, group_jobs=99) == (8, 1)
    assert resolve_group_parallelism(group_count=2, threads=8, group_jobs=4) == (2, 4)
    assert resolve_group_parallelism(group_count=0, threads=8, group_jobs=4) == (0, 0)


def test_build_cluster_groups_skips_fully_anchored_species():
    table = pd.DataFrame(
        [
            {
                "seq_id": "seq1",
                "anchor_rank": "species",
                "phylum": "p__P",
                "class": "c__C",
                "order": "o__O",
                "family": "f__F",
                "genus": "g__G",
                "species": "s__S",
            },
            {
                "seq_id": "seq2",
                "anchor_rank": "family",
                "phylum": "p__P",
                "class": "c__C",
                "order": "o__O",
                "family": "f__F",
                "genus": "",
                "species": "",
            },
        ]
    )

    groups = build_cluster_groups(table)

    assert len(groups) == 1
    assert groups[0]["seq_ids"] == {"seq2"}
    assert groups[0]["novel_ranks_fine"] == ["species", "genus"]
