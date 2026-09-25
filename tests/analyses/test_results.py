"""The tidy tables ``load_results`` builds from a study's comparison results."""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path
from typing import Any, ClassVar

import numpy as np
import pytest
from click.testing import CliRunner
from pydantic import BaseModel

from polyzymd.analyses import load_results
from polyzymd.analyses.contract import Observable, contract_analysis, iter_frames
from polyzymd.analyses.orchestrator import run_comparison
from polyzymd.analyses.testing import synthetic_universe
from tests.analyses.conftest import make_comparison

VALUES = {"No Polymer": 1.0, "SBMA": 2.0}


class _Settings(BaseModel):
    pass


class _Probe:
    """A scalar that differs by condition and a three-residue profile."""

    name: ClassVar[str] = "results_probe"
    Settings: ClassVar[type[BaseModel]] = _Settings
    references: ClassVar[tuple[str, ...]] = ()

    def compute(self, universe: Any, frames: Any, settings: _Settings) -> list[Observable]:
        base = VALUES[universe.probe_label] + 0.1 * universe.probe_replicate
        series = [base for _ in iter_frames(universe, frames)]
        return [
            Observable(name="value", kind="mean_of_timeseries", unit="A", values=series),
            Observable(
                name="per_residue",
                kind="profile",
                unit="A",
                values=np.asarray([base, 2 * base, 3 * base]),
                index=np.asarray([10, 11, 12]),
                index_label="Residue",
            ),
        ]


ProbeAnalysis = contract_analysis(_Probe)


@pytest.fixture
def study(tmp_path: Path, serve_replicates: Any) -> Path:
    """A study with two comparisons of the probe."""

    def universe_for(label: str, replicate: int) -> Any:
        universe = synthetic_universe()
        universe.probe_label, universe.probe_replicate = label, replicate
        return universe

    serve_replicates(universe_for)
    root = tmp_path / "thermal"
    (root / "comparisons").mkdir(parents=True)
    (root / "study.yaml").write_text("name: Lipase thermal stability\n")
    for folder, title in (("calb_343K", "CALB 343 K"), ("rml_333K", "RML 333 K")):
        directory = root / "comparisons" / folder
        directory.mkdir()
        (directory / "comparison.yaml").write_text(f"name: {title}\n")
        run_comparison(
            ProbeAnalysis(),
            make_comparison(directory, labels=tuple(VALUES), control="No Polymer"),
        )
    return root


def test_conditions_hold_one_row_per_condition_and_scalar(study: Path) -> None:
    results = load_results(study)

    assert len(results.conditions) == 4
    row = next(
        r
        for r in results.conditions
        if r["comparison"] == "CALB 343 K" and r["condition"] == "SBMA"
    )
    assert row["study"] == "Lipase thermal stability"
    assert row["analysis"] == "results_probe"
    assert row["observable"] == "value"
    assert row["n_replicates"] == 3
    assert row["mean"] == pytest.approx(2.2)
    assert row["replicate_values"] == pytest.approx([2.1, 2.2, 2.3])
    assert row["is_control"] is False
    assert row["ci95_low"] < row["mean"] < row["ci95_high"]


def test_comparisons_hold_one_row_per_test(study: Path) -> None:
    rows = load_results(study).comparisons

    assert {row["comparison"] for row in rows} == {"CALB 343 K", "RML 333 K"}
    row = next(r for r in rows if r["comparison"] == "CALB 343 K" and r["observable"] == "value")
    assert (row["control"], row["condition"]) == ("No Polymer", "SBMA")
    assert row["delta"] == pytest.approx(1.0)
    assert row["p_adjusted"] is not None


def test_profiles_hold_one_row_per_index(study: Path) -> None:
    rows = [r for r in load_results(study).profiles if r["comparison"] == "CALB 343 K"]

    assert len(rows) == 6
    sbma = [r for r in rows if r["condition"] == "SBMA"]
    assert [r["index"] for r in sbma] == [10, 11, 12]
    assert [r["mean"] for r in sbma] == pytest.approx([2.2, 4.4, 6.6])
    assert sbma[0]["index_label"] == "Residue"


def test_csv_and_dataframe_exports(study: Path, tmp_path: Path) -> None:
    results = load_results(study)

    written = results.to_csv(tmp_path / "out")
    frame = results.to_dataframe("comparisons")

    assert [path.name for path in written] == [
        "conditions.csv",
        "comparisons.csv",
        "profiles.csv",
    ]
    with written[0].open() as stream:
        first = next(csv.DictReader(stream))
    assert json.loads(first["replicate_values"])
    assert set(frame["comparison"]) == {"CALB 343 K", "RML 333 K"}


def test_only_the_comparison_results_are_needed(study: Path) -> None:
    """A published study without replicate artifacts or trajectories still reads."""
    for analysis_dir in study.glob("comparisons/*/analysis"):
        shutil.rmtree(analysis_dir)

    assert len(load_results(study).conditions) == 4


def test_sources_can_be_a_comparison_file_or_a_glob(study: Path) -> None:
    one = load_results(study / "comparisons" / "calb_343K" / "comparison.yaml")
    both = load_results(str(study / "comparisons" / "*"))

    assert {row["comparison"] for row in one.conditions} == {"CALB 343 K"}
    assert len(both.conditions) == 4


def test_analysis_filter_and_missing_results(study: Path) -> None:
    (study / "comparisons" / "empty").mkdir()
    (study / "comparisons" / "empty" / "comparison.yaml").write_text("name: empty\n")

    results = load_results(study, analyses=["not_run"])

    assert results.conditions == []
    assert len(results.warnings) == 3


def test_the_study_results_command_writes_the_tables(
    study: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from polyzymd.cli.main import cli

    monkeypatch.chdir(study / "comparisons")
    outcome = CliRunner().invoke(cli, ["study", "results"])

    assert outcome.exit_code == 0, outcome.output
    assert "4 condition rows" in outcome.output
    assert (study / "results" / "comparisons.csv").is_file()
