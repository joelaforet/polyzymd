"""Incomplete data is computed, and every result says how incomplete it is.

A study is looked at long before it finishes: some replicates have not
started, others are still being written. Nothing here is refused. Each test
checks that a result computed from partial data records what it used and says
so wherever the result is read: the aggregate, the comparison, the results
table, the agent report, the figure footnote and the published package.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any, ClassVar

import pytest
from pydantic import BaseModel

from polyzymd.analyses import loading, orchestrator
from polyzymd.analyses.completeness import summary
from polyzymd.analyses.contract import Observable, contract_analysis, iter_frames
from polyzymd.analyses.testing import synthetic_universe
from tests.analyses.conftest import make_comparison, make_condition

TIMESTEP_PS = 40.0
N_FRAMES = 5  # 5 frames at 40 ps cover 0.2 ns


class _Settings(BaseModel):
    pass


class _Probe:
    name: ClassVar[str] = "partial_probe"
    Settings: ClassVar[type[BaseModel]] = _Settings
    references: ClassVar[tuple[str, ...]] = ()

    def compute(self, universe: Any, frames: Any, settings: _Settings) -> list[Observable]:
        values = [float(universe.probe_replicate) for _ in iter_frames(universe, frames)]
        return [Observable(name="value", kind="mean_of_timeseries", unit="A", values=values)]


ProbeAnalysis = contract_analysis(_Probe)


@pytest.fixture
def serve(serve_replicates: Any, monkeypatch: pytest.MonkeyPatch) -> Any:
    """Serve replicates; returns a setter for planned length, missing and running ones."""
    state = SimpleNamespace(planned_ns=0.2, missing=set(), running=set(), require_complete=[])

    def universe_for(label: str, replicate: int) -> Any:
        if (label, replicate) in state.missing:
            raise FileNotFoundError(f"no trajectory for {label} replicate {replicate}")
        universe = synthetic_universe(n_frames=N_FRAMES)
        universe.probe_replicate = replicate
        return universe

    def status(replicate: int) -> dict[str, str]:
        return {"0": "completed", "1": "running" if replicate in state.running else "completed"}

    serve_replicates(universe_for, (), timestep_ps=TIMESTEP_PS, segment_status=status)
    served_open = loading.open_replicate

    def recording_open(config: Any, replicate: int, equilibration: str, **kwargs: Any) -> Any:
        state.require_complete.append(kwargs.get("require_complete", True))
        return served_open(config, replicate, equilibration, **kwargs)

    monkeypatch.setattr(loading, "open_replicate", recording_open)
    original = orchestrator.Condition.from_condition_config

    def with_plan(cfg: Any) -> Any:
        condition = original(cfg)
        _plan(condition, state)
        return condition

    monkeypatch.setattr(orchestrator.Condition, "from_condition_config", staticmethod(with_plan))
    return state


def _plan(condition: Any, state: Any) -> Any:
    condition.sim_config.simulation_phases = SimpleNamespace(
        production=SimpleNamespace(duration=state.planned_ns)
    )
    return condition


def _condition(state: Any, tmp_path: Path, label: str = "A", replicates=(1, 2, 3)) -> Any:
    return _plan(make_condition(label, tmp_path, replicates), state)


def _run(state: Any, tmp_path: Path, **kwargs: Any) -> Any:
    return orchestrator.run_analysis(
        ProbeAnalysis(),
        _condition(state, tmp_path),
        _Settings(),
        "0ns",
        tmp_path / "analysis" / "A",
        **kwargs,
    )


def test_a_replicate_records_the_time_it_covers_against_the_plan(
    serve: Any, tmp_path: Path
) -> None:
    serve.planned_ns = 0.4

    artifact = orchestrator.run_replicate_once(
        ProbeAnalysis(), _condition(serve, tmp_path), _Settings(), "0ns", tmp_path / "r", 1, True
    )

    record = artifact.metadata["completeness"]
    assert record["frames_used"] == N_FRAMES
    assert record["covered_ns"] == pytest.approx(0.2)
    assert record["planned_ns"] == pytest.approx(0.4)
    assert record["production_fraction"] == pytest.approx(0.5)
    assert record["complete"] is False


def test_full_data_is_complete(serve: Any, tmp_path: Path) -> None:
    record = _run(serve, tmp_path).metadata["completeness"]

    assert record["complete"] is True
    assert record["replicates_used"] == [1, 2, 3]
    assert summary("A", record) is None


def test_a_missing_replicate_is_left_out_and_recorded(serve: Any, tmp_path: Path) -> None:
    serve.missing.add(("A", 2))

    artifact = _run(serve, tmp_path)

    record = artifact.metadata["completeness"]
    assert artifact.replicates == [1, 3]
    assert record["replicates_listed"] == [1, 2, 3]
    assert record["replicates_used"] == [1, 3]
    assert record["dropped"][0]["replicate"] == 2
    assert "no trajectory" in record["dropped"][0]["reason"]
    assert record["complete"] is False
    assert summary("A", record) == "A: replicates 1, 3 of 1-3"


def test_a_running_segment_marks_the_replicate_incomplete(serve: Any, tmp_path: Path) -> None:
    serve.running.add(2)
    serve.planned_ns = 0.4

    record = _run(serve, tmp_path).metadata["completeness"]

    assert record["replicates"]["2"]["unfinished_segments"] == {"1": "running"}
    assert record["complete"] is False
    assert "replicate 2 at 50% of 0.4 ns (still running)" in summary("A", record)


def test_include_running_reads_unfinished_segments(serve: Any, tmp_path: Path) -> None:
    _run(serve, tmp_path, include_running=True)

    assert serve.require_complete and not any(serve.require_complete)


class TestWhereItShows:
    """A partial comparison says so wherever its numbers are read."""

    @pytest.fixture
    def partial(self, serve: Any, tmp_path: Path) -> dict[str, Any]:
        serve.missing.add(("B", 3))
        return orchestrator.run_comparison(ProbeAnalysis(), make_comparison(tmp_path))

    def test_the_comparison_is_not_complete(self, partial: dict[str, Any]) -> None:
        completeness = partial["completeness"]

        assert completeness["complete"] is False
        assert completeness["conditions"]["A"]["complete"] is True
        assert partial["comparison"].metadata["completeness"] == completeness

    def test_the_text_report_leads_with_it(self, partial: dict[str, Any]) -> None:
        text = ProbeAnalysis().format(partial["comparison"])

        assert text.splitlines()[1] == "PARTIAL: B: replicates 1-2 of 1-3"

    def test_the_results_table_marks_every_row(
        self, partial: dict[str, Any], tmp_path: Path
    ) -> None:
        from polyzymd.analyses import load_results

        results = load_results(tmp_path)
        by_condition = {row["condition"]: row for row in results.conditions}

        assert by_condition["A"]["complete"] is True
        assert by_condition["B"]["complete"] is False
        assert by_condition["B"]["replicates_listed"] == [1, 2, 3]
        assert all(row["complete"] is False for row in results.comparisons)

    def test_the_agent_report_leads_with_it(self, partial: dict[str, Any], tmp_path: Path) -> None:
        from polyzymd.analyses.protocols import build_report

        report = build_report(ProbeAnalysis(), make_comparison(tmp_path), partial)

        assert report.complete is False
        assert report.partial == ["B: replicates 1-2 of 1-3"]
        assert report.to_agent_text().splitlines()[1] == "PARTIAL: B: replicates 1-2 of 1-3"

    def test_the_figure_footnote_says_it(
        self, serve: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from polyzymd.analyses.shared import plotting

        written: list[str] = []
        original = plotting.add_uncertainty_footnote

        def recording(fig: Any, **kwargs: Any) -> str:
            written.append(original(fig, **kwargs))
            return written[-1]

        monkeypatch.setattr(plotting, "add_uncertainty_footnote", recording)
        serve.missing.add(("B", 3))
        orchestrator.run_comparison(ProbeAnalysis(), make_comparison(tmp_path))

        assert written and all(
            text.endswith("Partial: B: replicates 1-2 of 1-3.") for text in written
        )

    def test_export_names_the_partial_result(self, partial: dict[str, Any], tmp_path: Path) -> None:
        from polyzymd.study_bundle import plan_export

        (tmp_path / "study.yaml").write_text("name: s\n")
        (tmp_path / "comparisons").mkdir()
        (tmp_path / "comparison.yaml").write_text("name: project\nconditions: []\n")

        plan = plan_export(tmp_path)

        assert plan.partial == ["comparison/partial_probe: B: replicates 1-2 of 1-3"]
