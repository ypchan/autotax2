from __future__ import annotations

from typer.testing import CliRunner

from autotax2.cli import app


runner = CliRunner()


def test_root_help_works() -> None:
    result = runner.invoke(app, ["--help"])

    assert result.exit_code == 0
    assert "autotax2" in result.output


def test_placeholder_command_help_works() -> None:
    for command in [
        "init",
        "resolve-silva",
        "prepare-dataset",
        "orient-sina",
        "cluster-search",
        "place",
        "add",
        "export",
        "summarize",
        "validate",
    ]:
        result = runner.invoke(app, [command, "--help"])

        assert result.exit_code == 0, result.output
        assert command in result.output
