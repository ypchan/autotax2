from __future__ import annotations

import re
import shutil
import uuid
from collections.abc import Iterator
from pathlib import Path

import pytest

from autotax2.audit import write_event_log


@pytest.fixture
def audit_tmp_dir() -> Iterator[Path]:
    tmp_dir = Path("tests") / f"tmp_audit_{uuid.uuid4().hex}"
    tmp_dir.mkdir(parents=True)
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def test_write_event_log_uses_dated_filename_and_key_values(audit_tmp_dir: Path) -> None:
    path = write_event_log(
        audit_tmp_dir / "build",
        "prepare",
        {"dataset": "digester2020", "records": 3},
    )

    assert re.fullmatch(r"prepare_date[0-9]{14}\.log", path.name)
    content = path.read_text(encoding="utf-8")
    assert "command=prepare\n" in content
    assert "dataset=digester2020\n" in content
    assert "records=3\n" in content
