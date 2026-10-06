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


def test_shipped_analyses_hash_their_modules(monkeypatch) -> None:
    """A fix in a module a shipped analysis imports changes its record."""
    from polyzymd.analyses import functions, timeseries

    record = timeseries._function_record(functions.rms_decomposition)
    assert record["hash_of"] == "polyzymd_modules"
    timeseries._shipped_code_hash.cache_clear()
    reference = Path(__import__("polyzymd.analyses.reference").analyses.reference.__file__)
    real = Path.read_bytes

    def edited(self):
        data = real(self)
        return data + b"# fixed\n" if self == reference else data

    monkeypatch.setattr(Path, "read_bytes", edited)
    try:
        assert timeseries._function_record(functions.rms_decomposition)["hash"] != record["hash"]
    finally:
        monkeypatch.undo()
        timeseries._shipped_code_hash.cache_clear()
