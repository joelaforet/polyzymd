"""Stored records of timeseries analyses: what keys reuse and what is hashed."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

import polyzymd as pz
from polyzymd.cli.main import cli
from tests._support.analysis_testkit import write_committed_study

pytest.importorskip("MDAnalysis")
pytestmark = pytest.mark.filterwarnings("ignore")


def _analyze(*arguments: str):
    return CliRunner().invoke(cli, ["analyze", *arguments], catch_exceptions=False)


def _json(result) -> dict:
    return json.loads(result.output[result.output.index("{") :])


@pytest.mark.usefixtures("git_identity")
def test_reordered_parts_recompute(tmp_path: Path) -> None:
    """The order of parts names the stored columns, so it keys reuse."""
    root = write_committed_study(
        tmp_path,
        "  pair: {function: analyses/two.py:two, kind: timeseries, universe: u, parts: [a, b]}\n",
    )
    (root / "analyses" / "two.py").write_text("def two(u):\n    return {'a': 1.0, 'b': 100.0}\n")
    options = ["--study", str(root), "--no-plots", "--no-eq-check", "--run", "a"]
    first = _json(_analyze("pair", *options, "--format", "json"))
    study_yaml = root / "study.yaml"
    study_yaml.write_text(study_yaml.read_text().replace("parts: [a, b]", "parts: [b, a]"))
    second = _json(_analyze("pair", *options, "--format", "json"))
    assert [c["mean"] for c in first["conditions"]] == [1.0, 1.0]
    assert [c["mean"] for c in second["conditions"]] == [1.0, 1.0]
    means = pz.Study(root).results("pair").table.groupby("part")["value"].mean()
    assert means["a"] == 1.0 and means["b"] == 100.0


def test_shipped_analyses_hash_every_analysis_file(tmp_path: Path) -> None:
    """A fix in any analysis module, even one imported two levels down, changes the hash."""
    import shutil

    from polyzymd.analyses import functions, timeseries

    assert timeseries._function_record(functions.rms_decomposition)["hash_of"] == "polyzymd_modules"
    package = Path(timeseries.__file__).parent
    before, after = tmp_path / "before", tmp_path / "after"
    for copy in (before, after):
        shutil.copytree(package, copy, ignore=shutil.ignore_patterns("__pycache__"))
    with open(after / "shared" / "centroid.py", "a") as handle:
        handle.write("# fixed\n")
    assert timeseries._shipped_code_hash(before) != timeseries._shipped_code_hash(after)


def test_stray_files_do_not_change_the_code_hash(tmp_path: Path, caplog) -> None:
    """A function's hash covers the Python files beside it and data/, not other files.

    The analysis log names the files it leaves out.
    """
    import logging

    from polyzymd.analyses.timeseries import folder_hash

    analyses = tmp_path / "analyses"
    (analyses / "data").mkdir(parents=True)
    (analyses / "f.py").write_text("def f(u):\n    return 1.0\n")
    (analyses / "data" / "table.csv").write_text("1\n")
    before = folder_hash(analyses / "f.py")
    (analyses / "notes.txt").write_text("my notes")
    (analyses / "copy.xtc").write_bytes(b"\0" * 1000)
    with caplog.at_level(logging.INFO):
        assert folder_hash(analyses / "f.py") == before
    assert "copy.xtc, notes.txt" in caplog.text and "data/" in caplog.text
    (analyses / "data" / "table.csv").write_text("2\n")
    assert folder_hash(analyses / "f.py") != before
