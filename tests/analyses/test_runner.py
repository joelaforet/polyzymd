"""What the analysis runner does, pinned through its public entry points.

Every test here runs one contract plugin on in-memory replicates served by
``serve_replicates`` and calls the functions ``cli/compare.py`` and the SLURM
workers call. None of them reaches into how the runner is built, so they hold
for any implementation: they are the specification of the runner.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, ClassVar

import pytest
from pydantic import BaseModel

from polyzymd.analyses import discovery, orchestrator
from polyzymd.analyses.contract import Observable, contract_analysis, iter_frames
from polyzymd.analyses.exceptions import (
    PluginContractError,
    ReplicateError,
    ReplicateSkippedError,
)
from polyzymd.analyses.mda.universe import FileIdentity
from polyzymd.analyses.testing import synthetic_universe
from tests.analyses.conftest import make_comparison, make_condition

NAME = "runner_probe"


class _ProbeSettings(BaseModel):
    """Settings for the probe."""

    scale: float = 1.0


class _Probe:
    """Contract plugin whose value is ``scale`` times the replicate number.

    ``failures`` maps a (condition label, replicate) pair to the exception
    ``compute()`` raises for it, and ``calls`` records every pair computed.
    """

    name: ClassVar[str] = NAME
    Settings: ClassVar[type[BaseModel]] = _ProbeSettings
    references: ClassVar[tuple[str, ...]] = ()
    calls: ClassVar[list[tuple[str, int]]] = []
    failures: ClassVar[dict[tuple[str, int], BaseException]] = {}

    def compute(self, universe: Any, frames: Any, settings: _ProbeSettings) -> list[Observable]:
        """Report ``scale * replicate`` on every frame."""
        key = (universe.probe_label, universe.probe_replicate)
        type(self).calls.append(key)
        if key in type(self).failures:
            raise type(self).failures[key]
        value = settings.scale * universe.probe_replicate
        values = [value for _ in iter_frames(universe, frames)]
        return [Observable(name="value", kind="mean_of_timeseries", unit="A", values=values)]


ProbeAnalysis = contract_analysis(_Probe)


@pytest.fixture
def served(serve_replicates: Any, tmp_path: Path) -> Any:
    """Serve one universe per replicate and register the probe by name.

    Returns the mapping of (label, replicate) to an exception raised while that
    replicate loads, standing in for missing trajectory files.
    """
    _Probe.calls = []
    _Probe.failures = {}
    load_failures: dict[tuple[str, int], BaseException] = {}

    def universe_for(label: str, replicate: int) -> Any:
        if (label, replicate) in load_failures:
            raise load_failures[(label, replicate)]
        universe = synthetic_universe()
        universe.probe_label = label
        universe.probe_replicate = replicate
        return universe

    trajectory = tmp_path / "prod.dcd"
    trajectory.write_bytes(b"DCD")
    serve_replicates(universe_for, [FileIdentity.from_path(trajectory).as_dict()])
    discovery.register_analysis(ProbeAnalysis)
    yield load_failures
    discovery.clear_cache()


@pytest.fixture(autouse=True)
def _serve(served: Any) -> None:
    """Every test in this module runs against the served probe."""


def _condition_dir(root: Path, label: str = "A") -> Path:
    return root / "analysis" / label / NAME


def _observable(artifact: Any) -> dict[str, Any]:
    return artifact.payload["observables"][0]


class TestReplicates:
    """One replicate at a time."""

    def test_run_replicate_once_writes_the_canonical_result(self, tmp_path: Path) -> None:
        run_dir = _condition_dir(tmp_path) / "run_1"

        artifact = orchestrator.run_replicate_once(
            ProbeAnalysis(),
            make_condition("A", tmp_path),
            _ProbeSettings(),
            "0ns",
            run_dir,
            1,
            True,
        )

        assert (run_dir / "result.json").is_file()
        assert artifact.replicate == 1
        assert _observable(artifact)["value"] == pytest.approx(1.0)

    def test_run_analysis_equals_replicates_then_aggregate_from_disk(self, tmp_path: Path) -> None:
        condition = make_condition("A", tmp_path, (1, 2))
        settings = _ProbeSettings(scale=2.0)
        direct = orchestrator.run_analysis(
            ProbeAnalysis(), condition, settings, "0ns", tmp_path / "direct"
        )
        stepwise_dir = tmp_path / "stepwise"
        for replicate in (1, 2):
            orchestrator.run_replicate_once(
                ProbeAnalysis(),
                condition,
                settings,
                "0ns",
                stepwise_dir / f"run_{replicate}",
                replicate,
                False,
            )
        stepwise = orchestrator.aggregate_condition_from_disk(
            ProbeAnalysis(), condition, settings, "0ns", stepwise_dir, (1, 2)
        )

        assert _observable(direct)["mean"] == _observable(stepwise)["mean"] == pytest.approx(3.0)
        assert direct.replicates == stepwise.replicates == [1, 2]

    def test_run_analysis_writes_replicate_and_aggregate_files(self, tmp_path: Path) -> None:
        output = _condition_dir(tmp_path)

        orchestrator.run_analysis(
            ProbeAnalysis(), make_condition("A", tmp_path, (1, 2)), _ProbeSettings(), "0ns", output
        )

        assert (output / "run_1" / "result.json").is_file()
        assert (output / "run_2" / "result.json").is_file()
        aggregate = json.loads((output / "aggregated" / "result.json").read_text())
        assert aggregate["replicates"] == [1, 2]
        assert aggregate["metadata"]["settings_fingerprint"]


def test_loader_warnings_reach_the_replicate_artifact(
    tmp_path: Path, serve_replicates: Any
) -> None:
    """What the loader noticed about the trajectory is kept with the numbers."""
    universe = synthetic_universe()
    universe.probe_label, universe.probe_replicate = "A", 1
    serve_replicates(universe, warnings=["segment 3 is still running"])

    artifact = orchestrator.run_replicate_once(
        ProbeAnalysis(),
        make_condition("A", tmp_path),
        _ProbeSettings(),
        "0ns",
        tmp_path / "run_1",
        1,
        True,
    )

    assert "segment 3 is still running" in artifact.warnings


class TestReplicateFailures:
    """How a failing replicate affects its condition."""

    def test_missing_trajectory_skips_the_replicate(
        self, tmp_path: Path, served: dict[tuple[str, int], BaseException]
    ) -> None:
        served[("A", 2)] = FileNotFoundError("prod.dcd not found")

        artifact = orchestrator.run_analysis(
            ProbeAnalysis(), make_condition("A", tmp_path), _ProbeSettings(), "0ns", tmp_path / "o"
        )

        assert artifact.replicates == [1, 3]

    def test_skipped_replicate_is_logged(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        _Probe.failures[("A", 2)] = ReplicateSkippedError("segment still running")

        with caplog.at_level(logging.WARNING):
            artifact = orchestrator.run_analysis(
                ProbeAnalysis(),
                make_condition("A", tmp_path),
                _ProbeSettings(),
                "0ns",
                tmp_path / "o",
            )

        assert artifact.replicates == [1, 3]
        assert "segment still running" in caplog.text

    def test_every_replicate_skipped_fails_the_minimum(self, tmp_path: Path) -> None:
        for replicate in (1, 2):
            _Probe.failures[("A", replicate)] = ReplicateSkippedError("not ready")

        with pytest.raises(ValueError, match="need at least 1"):
            orchestrator.run_analysis(
                ProbeAnalysis(),
                make_condition("A", tmp_path, (1, 2)),
                _ProbeSettings(),
                "0ns",
                tmp_path / "o",
            )

    def test_unexpected_error_is_a_replicate_error(self, tmp_path: Path) -> None:
        _Probe.failures[("X", 1)] = RuntimeError("boom")

        with pytest.raises(ReplicateError, match="condition='X' replicate=1"):
            orchestrator.run_analysis(
                ProbeAnalysis(),
                make_condition("X", tmp_path),
                _ProbeSettings(),
                "0ns",
                tmp_path / "o",
            )

    def test_contract_error_propagates_unwrapped(self, tmp_path: Path) -> None:
        _Probe.failures[("A", 1)] = PluginContractError("bad observable")

        with pytest.raises(PluginContractError, match="bad observable"):
            orchestrator.run_analysis(
                ProbeAnalysis(),
                make_condition("A", tmp_path),
                _ProbeSettings(),
                "0ns",
                tmp_path / "o",
            )


class TestAggregateFromDisk:
    """Aggregating replicate results a worker already wrote."""

    def _write(self, root: Path, replicates: tuple[int, ...]) -> None:
        for replicate in replicates:
            orchestrator.run_replicate_once(
                ProbeAnalysis(),
                make_condition("A", root, (1, 2)),
                _ProbeSettings(),
                "0ns",
                _condition_dir(root) / f"run_{replicate}",
                replicate,
                False,
            )

    def test_loads_replicates_and_writes_the_aggregate(self, tmp_path: Path) -> None:
        self._write(tmp_path, (1, 2))

        artifact = orchestrator.aggregate_condition_from_disk(
            ProbeAnalysis(),
            make_condition("A", tmp_path, (1, 2)),
            _ProbeSettings(),
            "0ns",
            _condition_dir(tmp_path),
            (1, 2),
        )

        assert artifact.replicates == [1, 2]
        assert (_condition_dir(tmp_path) / "aggregated" / "result.json").is_file()

    def test_a_missing_replicate_is_tolerated_above_the_minimum(self, tmp_path: Path) -> None:
        self._write(tmp_path, (1,))

        artifact = orchestrator.aggregate_condition_from_disk(
            ProbeAnalysis(),
            make_condition("A", tmp_path, (1, 2)),
            _ProbeSettings(),
            "0ns",
            _condition_dir(tmp_path),
            (1, 2),
        )

        assert artifact.replicates == [1]

    def test_no_replicates_reports_the_expected_paths(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError) as excinfo:
            orchestrator.aggregate_condition_from_disk(
                ProbeAnalysis(),
                make_condition("A", tmp_path, (1, 2)),
                _ProbeSettings(),
                "0ns",
                _condition_dir(tmp_path),
                (1, 2),
            )

        message = str(excinfo.value)
        assert "replicate result(s) on disk" in message
        assert str(_condition_dir(tmp_path) / "run_1" / "result.json") in message

    def test_recompute_clears_the_aggregate_directory(self, tmp_path: Path) -> None:
        self._write(tmp_path, (1, 2))
        stale = _condition_dir(tmp_path) / "aggregated" / "stale_sidecar.txt"
        stale.parent.mkdir(parents=True)
        stale.write_text("stale")

        orchestrator.aggregate_condition_from_disk(
            ProbeAnalysis(),
            make_condition("A", tmp_path, (1, 2)),
            _ProbeSettings(),
            "0ns",
            _condition_dir(tmp_path),
            (1, 2),
            recompute=True,
        )

        assert not stale.exists()

    def test_a_stale_aggregate_is_replaced_without_recompute(self, tmp_path: Path) -> None:
        output = _condition_dir(tmp_path)
        (output / "aggregated").mkdir(parents=True)
        (output / "aggregated" / "result.json").write_text(json.dumps({"stale": True}))

        orchestrator.run_analysis(
            ProbeAnalysis(), make_condition("A", tmp_path, (1, 2)), _ProbeSettings(), "0ns", output
        )

        assert "stale" not in json.loads((output / "aggregated" / "result.json").read_text())


class TestComparison:
    """The whole pipeline for one analysis."""

    def test_prepare_resolves_conditions_settings_and_window(self, tmp_path: Path) -> None:
        config = make_comparison(tmp_path, settings={NAME: _ProbeSettings(scale=2.0)})

        prepared = orchestrator.prepare_comparison_run(ProbeAnalysis(), config, "5ns")

        assert [condition.label for condition in prepared["valid_conditions"]] == ["A", "B"]
        assert prepared["settings"].scale == 2.0
        assert prepared["equilibration"] == "5ns"
        assert prepared["analysis_root"] == tmp_path / "analysis"

    def test_run_comparison_returns_aggregates_comparison_and_plots(self, tmp_path: Path) -> None:
        result = orchestrator.run_comparison(ProbeAnalysis(), make_comparison(tmp_path))

        assert set(result["aggregated"]) == {"A", "B"}
        assert result["comparison"].control_label == "A"
        assert result["comparison_path"] == tmp_path / "comparison" / NAME / "result.json"
        assert result["comparison_path"].is_file()
        assert result["plots"] and all(path.exists() for path in result["plots"])
        assert all(path.parent == tmp_path / "figures" / NAME for path in result["plots"])

    def test_recompute_replaces_outputs_and_keeps_sibling_figures(self, tmp_path: Path) -> None:
        config = make_comparison(tmp_path)
        orchestrator.run_comparison(ProbeAnalysis(), config)
        stale_figure = tmp_path / "figures" / NAME / "stale.png"
        stale_figure.write_text("stale")
        sibling = tmp_path / "figures" / "other_analysis" / "keep.png"
        sibling.parent.mkdir(parents=True)
        sibling.write_text("keep")
        stale_aggregate = _condition_dir(tmp_path) / "aggregated" / "stale_sidecar.txt"
        stale_aggregate.write_text("stale")

        orchestrator.run_comparison(ProbeAnalysis(), config, recompute=True)

        assert not stale_figure.exists()
        assert not stale_aggregate.exists()
        assert sibling.exists()

    def test_contract_error_stops_before_later_conditions(self, tmp_path: Path) -> None:
        _Probe.failures[("A", 1)] = PluginContractError("bad observable")

        with pytest.raises(PluginContractError):
            orchestrator.run_comparison(ProbeAnalysis(), make_comparison(tmp_path))

        assert all(label == "A" for label, _ in _Probe.calls)

    def test_a_failed_condition_makes_the_comparison_fail(
        self, tmp_path: Path, served: dict[tuple[str, int], BaseException]
    ) -> None:
        for replicate in (1, 2, 3):
            served[("B", replicate)] = FileNotFoundError("no trajectory")

        with pytest.raises(ValueError, match="missing aggregated results"):
            orchestrator.run_comparison(ProbeAnalysis(), make_comparison(tmp_path))

    def test_execution_summary_is_logged(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        with caplog.at_level(logging.INFO, logger="polyzymd.analyses"):
            orchestrator.run_comparison(ProbeAnalysis(), make_comparison(tmp_path))

        assert "Mode: sequential (local)" in caplog.text
        assert "total replicate tasks" in caplog.text

    @pytest.mark.parametrize(
        ("sbatch", "expected"),
        [
            ("/usr/bin/sbatch", "Consider submitting to SLURM"),
            (None, "If you have access to an HPC cluster with SLURM"),
        ],
    )
    def test_expensive_analyses_suggest_slurm(
        self,
        tmp_path: Path,
        caplog: pytest.LogCaptureFixture,
        monkeypatch: pytest.MonkeyPatch,
        sbatch: str | None,
        expected: str,
    ) -> None:
        monkeypatch.setattr(orchestrator.shutil, "which", lambda command: sbatch)
        monkeypatch.setattr(ProbeAnalysis, "execution_cost_hint", "high")

        with caplog.at_level(logging.WARNING, logger="polyzymd.analyses"):
            orchestrator.run_comparison(ProbeAnalysis(), make_comparison(tmp_path))

        assert expected in caplog.text
        assert f"polyzymd compare submit {NAME}" in caplog.text


class TestFinalize:
    """Compare and plot from aggregates the aggregate workers wrote."""

    def _aggregate(self, root: Path, labels: tuple[str, ...]) -> dict[str, Path]:
        dirs = {}
        for label in labels:
            dirs[label] = _condition_dir(root, label)
            orchestrator.run_analysis(
                ProbeAnalysis(),
                make_condition(label, root, (1, 2)),
                _ProbeSettings(),
                "0ns",
                dirs[label],
            )
        return dirs

    def _finalize(self, root: Path, dirs: dict[str, Path], **kwargs: Any) -> dict[str, Any]:
        config = make_comparison(root, replicates=(1, 2))
        return orchestrator.finalize_comparison_from_disk(
            analysis=ProbeAnalysis(),
            config=config,
            analysis_dirs=dirs,
            aggregated_results={},
            results_dir=root / "comparison" / NAME,
            figures_dir=root / "figures" / NAME,
            settings=_ProbeSettings(),
            effective_control=config.control,
            **kwargs,
        )

    def test_loads_aggregates_from_disk_and_writes_outputs(self, tmp_path: Path) -> None:
        result = self._finalize(tmp_path, self._aggregate(tmp_path, ("A", "B")))

        assert result["comparison_path"].is_file()
        assert result["plots"]
        assert result["comparison"].conditions == ["A", "B"]

    def test_strict_finalize_names_the_missing_aggregate(self, tmp_path: Path) -> None:
        dirs = self._aggregate(tmp_path, ("A",))
        dirs["B"] = _condition_dir(tmp_path, "B")

        with pytest.raises(ValueError) as excinfo:
            self._finalize(tmp_path, dirs)

        message = str(excinfo.value)
        assert "missing aggregated results for condition(s)" in message
        assert str(dirs["B"] / "aggregated" / "result.json") in message
        assert "--allow-partial" in message

    def test_strict_finalize_refuses_a_missing_control(self, tmp_path: Path) -> None:
        dirs = self._aggregate(tmp_path, ("B",))
        dirs["A"] = _condition_dir(tmp_path, "A")

        with pytest.raises(ValueError, match="missing aggregated results"):
            self._finalize(tmp_path, dirs)

    def test_partial_finalize_without_the_control_compares_all_pairs(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        dirs = self._aggregate(tmp_path, ("B",))
        dirs["A"] = _condition_dir(tmp_path, "A")

        with caplog.at_level(logging.WARNING):
            result = self._finalize(tmp_path, dirs, allow_partial=True)

        assert result["comparison"].control_label is None
        assert result["comparison"].conditions == ["B"]
        assert "without a designated control (all-vs-all)" in caplog.text

    def test_partial_finalize_with_no_survivor_fails(self, tmp_path: Path) -> None:
        dirs = {label: _condition_dir(tmp_path, label) for label in ("A", "B")}

        with pytest.raises(ValueError, match="no aggregate files were found"):
            self._finalize(tmp_path, dirs, allow_partial=True)

    def test_an_aggregate_under_other_settings_counts_as_missing(self, tmp_path: Path) -> None:
        dirs = self._aggregate(tmp_path, ("A", "B"))
        config = make_comparison(tmp_path, replicates=(1, 2))

        with pytest.raises(ValueError, match="missing aggregated results"):
            orchestrator.finalize_comparison_from_disk(
                analysis=ProbeAnalysis(),
                config=config,
                analysis_dirs=dirs,
                aggregated_results={},
                results_dir=tmp_path / "comparison" / NAME,
                figures_dir=tmp_path / "figures" / NAME,
                settings=_ProbeSettings(scale=5.0),
                effective_control="A",
            )


class TestPlotOnly:
    """Redrawing figures from results already on disk."""

    def test_returns_the_figures(self, tmp_path: Path) -> None:
        config = make_comparison(tmp_path)
        orchestrator.run_comparison(ProbeAnalysis(), config)

        paths, failures = orchestrator.run_plot_only(ProbeAnalysis(), config)

        assert paths and failures == []

    def test_an_expected_plot_error_is_reported_not_raised(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from polyzymd.analyses import contract_plots

        config = make_comparison(tmp_path)
        orchestrator.run_comparison(ProbeAnalysis(), config)

        def broken(*args: Any, **kwargs: Any) -> Any:
            raise ValueError("no data")

        monkeypatch.setattr(contract_plots, "plot_observables", broken)
        paths, failures = orchestrator.run_plot_only(ProbeAnalysis(), config)

        assert paths == []
        assert failures == [
            (NAME, f"{NAME}: plot failed for comparison='project': ValueError: no data")
        ]

    def test_an_unexpected_plot_error_propagates(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from polyzymd.analyses import contract_plots

        config = make_comparison(tmp_path)
        orchestrator.run_comparison(ProbeAnalysis(), config)

        def broken(*args: Any, **kwargs: Any) -> Any:
            raise RuntimeError("bug")

        monkeypatch.setattr(contract_plots, "plot_observables", broken)
        with pytest.raises(RuntimeError, match="bug"):
            orchestrator.run_plot_only(ProbeAnalysis(), config)

    def test_run_all_plots_passes_the_equilibration_through(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        seen: list[Any] = []
        monkeypatch.setattr(
            orchestrator,
            "run_plot_only",
            lambda analysis, config, equilibration=None: seen.append(equilibration) or ([], []),
        )

        orchestrator.run_all_plots(make_comparison(tmp_path), [NAME], equilibration="7ns")

        assert seen == ["7ns"]


class TestSeveralAnalyses:
    """Entry points that take analysis names."""

    def test_order_keeps_the_requested_order_without_duplicates(self) -> None:
        assert orchestrator.order_analyses_for_execution([NAME, "rg", NAME]) == [NAME, "rg"]

    def test_run_all_comparisons_reports_failures_and_skips_unknown_names(
        self, tmp_path: Path, served: dict[tuple[str, int], BaseException]
    ) -> None:
        served[("B", 1)] = FileNotFoundError("no trajectory")
        served[("B", 2)] = FileNotFoundError("no trajectory")
        served[("B", 3)] = FileNotFoundError("no trajectory")

        results = orchestrator.run_all_comparisons(
            make_comparison(tmp_path), [NAME, "not_an_analysis"]
        )

        assert set(results) == {NAME}
        assert "error" in results[NAME]

    def test_run_all_comparisons_reraises_contract_errors(self, tmp_path: Path) -> None:
        _Probe.failures[("A", 1)] = PluginContractError("bad observable")

        with pytest.raises(PluginContractError):
            orchestrator.run_all_comparisons(make_comparison(tmp_path), [NAME])
