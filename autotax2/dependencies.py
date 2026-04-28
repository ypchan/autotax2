from __future__ import annotations

import importlib.util
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional

from .logging import print_kv


@dataclass
class CheckResult:
    name: str
    ok: bool
    detail: str
    required: bool = True


def check_executable(name: str, path: Optional[str] = None, version_args: Optional[list[str]] = None) -> CheckResult:
    exe = path or name
    found = shutil.which(exe) if "/" not in exe else (exe if Path(exe).exists() else None)
    if not found:
        return CheckResult(name, False, f"not found: {exe}")

    detail = str(found)
    if version_args:
        try:
            proc = subprocess.run([found] + version_args, capture_output=True, text=True, timeout=20)
            out = (proc.stdout or proc.stderr).strip().splitlines()
            if out:
                detail += f" | {out[0]}"
        except Exception as exc:
            detail += f" | version check failed: {exc}"
    return CheckResult(name, True, detail)


def check_python_module(module: str, required: bool = True) -> CheckResult:
    spec = importlib.util.find_spec(module)
    if spec is None:
        return CheckResult(module, False, "not importable", required=required)
    return CheckResult(module, True, "importable", required=required)


def check_paths(paths: Iterable[tuple[str, str, bool]]) -> List[CheckResult]:
    results: List[CheckResult] = []
    for name, path, required in paths:
        p = Path(path)
        ok = p.exists() and p.stat().st_size > 0
        detail = str(p) if ok else f"missing or empty: {p}"
        results.append(CheckResult(name, ok, detail, required=required))
    return results


def summarize_checks(results: list[CheckResult]) -> bool:
    rows = {}
    all_ok = True
    for r in results:
        status = "OK" if r.ok else ("MISSING" if r.required else "OPTIONAL_MISSING")
        rows[r.name] = f"{status} | {r.detail}"
        if r.required and not r.ok:
            all_ok = False
    print_kv("Dependency check", rows)
    return all_ok


def check_core_dependencies(vsearch: str = "vsearch") -> bool:
    results = [
        check_executable("vsearch", vsearch, ["--version"]),
        check_python_module("rich", required=False),
        check_python_module("rich_argparse", required=False),
    ]
    return summarize_checks(results)
