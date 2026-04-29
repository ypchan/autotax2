from src.threads import resolve_threads


def test_resolve_threads_auto():
    assert resolve_threads("auto") >= 1


def test_resolve_threads_int():
    assert resolve_threads(1) == 1
