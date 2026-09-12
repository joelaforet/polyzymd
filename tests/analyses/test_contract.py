"""Known-answer tests for the observable contract and its lifecycle adapter."""

from __future__ import annotations

import importlib.util
import math
import types
from pathlib import Path
from typing import Any

import MDAnalysis as mda
import numpy as np
import pytest
from MDAnalysis.coordinates.memory import MemoryReader

from polyzymd.analyses._framework.contexts import ComparisonContext, Condition
from polyzymd.analyses._framework.lifecycle import AnalysisLifecycle
from polyzymd.analyses.contract import (
    Observable,
    ObservableAggregate,
    aggregate_observables,
    compare_observables,
    reduce_observable,
)
from polyzymd.analyses.exceptions import PluginContractError
from polyzymd.analyses.mda.artifacts import ComparisonArtifact, ConditionArtifact
from polyzymd.analyses.rg_contract import Rg2Analysis, RgSettings
from polyzymd.analyses.shared.window import TrajectoryWindow

CROSS = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]], dtype=np.float32)


def _series(values: Any, kind: str = "mean_of_timeseries", name: str = "x") -> Observable:
    """Build one observable of the requested kind."""
    return Observable(name=name, kind=kind, unit="A", values=values)


def _universe(scale: float = 1.0, n_frames: int = 8) -> mda.Universe:
    """Four unit-mass atoms on a cross, held still, at the requested scale."""
    universe = mda.Universe.empty(4, n_residues=1, atom_resindex=[0] * 4, trajectory=True)
    universe.add_TopologyAttr("masses", [1.0] * 4)
    universe.load_new(np.stack([CROSS * scale] * n_frames), format=MemoryReader)
    return universe


def _sim_config(label: str) -> types.SimpleNamespace:
    """Minimal stand-in for the fields the config hash reads."""
    return types.SimpleNamespace(
        name=label,
        enzyme=types.SimpleNamespace(name="enzyme", pdb_path="/tmp/enzyme.pdb"),
        thermodynamics=types.SimpleNamespace(temperature=300.0, pressure=1.0),
        output=types.SimpleNamespace(
            projects_directory="/tmp/projects",
            effective_scratch_directory="/tmp/scratch",
            naming_template="{name}",
        ),
        substrate=None,
        polymers=None,
    )


def _stub_analysis_class(scale: float) -> type:
    """Rg2 with the trajectory loader replaced by an in-memory universe."""

    class _StubProvider:
        def __init__(self, config: Any, loader: Any = None) -> None:
            self.config = config

        def load_universe(self, replicate: int) -> mda.Universe:
            return _universe(scale=scale + 0.5 * replicate)

        def provenance_for(self, replicate: int) -> dict[str, Any]:
            return {
                "topology": {
                    "path": "/tmp/topology.pdb",
                    "format": "pdb",
                    "size_bytes": 1,
                    "mtime_ns": 2,
                },
                "trajectories": [],
                "warnings": [],
            }

    class _StubLoader:
        def __init__(self, config: Any) -> None:
            self.config = config

    class _StubRg2(Rg2Analysis):
        def _trajectory_loader_factory(self) -> type:
            return _StubLoader

        def _mda_universe_provider_factory(self) -> type:
            return _StubProvider

        def get_trajectory_window(
            self, ctx: Any, replicate: int, loader: Any, universe: Any
        ) -> TrajectoryWindow:
            n_frames = len(universe.trajectory)
            return TrajectoryWindow(
                start=0,
                stop=n_frames,
                step=1,
                equilibration_start=0,
                n_frames_total=n_frames,
                n_frames_selected=n_frames,
                timestep_ps=1.0,
                equilibration_ps=0.0,
                equilibration="0ns",
            )

    return _StubRg2


def _run_condition(label: str, scale: float, root: Path) -> tuple[Condition, ConditionArtifact]:
    """Run rg2 for one condition with three replicates and return its aggregate."""
    condition = Condition(
        label=label,
        config_path=root / f"{label}.yaml",
        replicates=(1, 2, 3),
        sim_config=_sim_config(label),
    )
    lifecycle = AnalysisLifecycle(_stub_analysis_class(scale)())
    aggregate = lifecycle.run_analysis(
        condition,
        RgSettings(runs=[{"label": "protein", "selection": "all"}]),
        "0ns",
        root / "analysis" / label / "rg2",
    )
    return condition, aggregate


def test_mean_of_timeseries_gives_sigma_over_root_n() -> None:
    """The SEM is the replicate standard deviation over the square root of n."""
    replicates = [[_series([value - 0.5, value + 0.5])] for value in (1.0, 2.0, 3.0, 4.0, 5.0)]

    aggregate = aggregate_observables(replicates)[0]

    assert aggregate.mean == pytest.approx(3.0)
    assert aggregate.sem == pytest.approx(math.sqrt(2.5) / math.sqrt(5))
    assert aggregate.n_replicates == 5
    assert aggregate.ci_method == "student_t"
    assert aggregate.coverage == pytest.approx(0.95)
    half_width = aggregate.ci95_high - aggregate.mean
    assert half_width == pytest.approx(2.776445 * aggregate.sem, rel=1e-5)


def test_fraction_reduces_to_the_bernoulli_rate() -> None:
    """A fraction observable reduces to the mean of its indicator series."""
    replicates = [
        [_series([1, 1, 0, 0, 0, 0, 0, 0, 0, 0], kind="fraction")],
        [_series([1, 1, 1, 1, 0, 0, 0, 0, 0, 0], kind="fraction")],
    ]

    aggregate = aggregate_observables(replicates)[0]

    assert aggregate.replicate_values == pytest.approx([0.2, 0.4])
    assert aggregate.mean == pytest.approx(0.3)
    assert aggregate.sem == pytest.approx(0.1)


def test_fraction_outside_the_unit_interval_is_rejected() -> None:
    """A fraction cannot carry values outside zero to one."""
    with pytest.raises(ValueError, match="fraction outside"):
        _series([0.0, 1.5], kind="fraction")


def test_fluctuation_reduces_to_the_sample_standard_deviation() -> None:
    """A fluctuation observable reduces to the standard deviation of its series."""
    estimate = reduce_observable(_series([1.0, 2.0, 3.0, 4.0, 5.0], kind="fluctuation"))

    assert estimate.value == pytest.approx(math.sqrt(2.5))
    assert estimate.n_frames == 5


def test_profile_averages_each_index_across_replicates() -> None:
    """A profile keeps its index and reports a mean and a SEM per index."""
    replicates = [
        [Observable(name="rmsf", kind="profile", unit="A", values=[1.0, 3.0], index=[10, 11])],
        [Observable(name="rmsf", kind="profile", unit="A", values=[3.0, 5.0], index=[10, 11])],
    ]

    aggregate = aggregate_observables(replicates)[0]

    assert aggregate.index == [10.0, 11.0]
    assert aggregate.profile_mean == pytest.approx([2.0, 4.0])
    assert aggregate.profile_sem == pytest.approx([1.0, 1.0])
    assert aggregate.mean is None


def test_statistical_inefficiency_is_a_diagnostic_not_a_correction() -> None:
    """N_eff is reported per condition and does not enter the SEM."""
    rng = np.random.default_rng(3)
    replicates = [[_series(rng.normal(0.0, 1.0, 400))] for _ in range(4)]

    aggregate = aggregate_observables(replicates)[0]

    assert aggregate.n_eff_min is not None
    assert aggregate.n_eff_min <= 400
    assert aggregate.sem == pytest.approx(
        float(np.std(aggregate.replicate_values, ddof=1)) / math.sqrt(4)
    )


def test_missing_observable_in_one_replicate_is_rejected() -> None:
    """Replicates must report the same observables."""
    replicates = [[_series([1.0, 2.0], name="a")], [_series([1.0, 2.0], name="b")]]

    with pytest.raises(PluginContractError, match="missing from some"):
        aggregate_observables(replicates)


def test_comparison_uses_one_benjamini_hochberg_family() -> None:
    """Every test in the run is adjusted together, and the largest is unchanged."""
    conditions = {
        "control": aggregate_observables(
            [[_series([v], name="a"), _series([v], name="b")] for v in (1.0, 1.1, 0.9)]
        ),
        "low": aggregate_observables(
            [[_series([v], name="a"), _series([v], name="b")] for v in (1.02, 1.11, 0.93)]
        ),
        "high": aggregate_observables(
            [[_series([v], name="a"), _series([v], name="b")] for v in (5.0, 5.1, 4.9)]
        ),
    }

    comparisons = compare_observables(conditions, control_label="control", fdr_alpha=0.05)

    assert len(comparisons) == 4
    assert {c.correction for c in comparisons} == {"benjamini_hochberg"}
    assert {c.test for c in comparisons} == {"student_t"}
    raw = [c.p_value for c in comparisons]
    adjusted = [c.p_adjusted for c in comparisons]
    assert all(a >= r for a, r in zip(adjusted, raw, strict=True))
    assert max(adjusted) == pytest.approx(max(raw))
    high = [c for c in comparisons if c.condition == "high"]
    assert all(c.significant for c in high)
    assert all(c.delta == pytest.approx(4.0, abs=1e-6) for c in high)


def test_profiles_are_not_tested_pairwise() -> None:
    """A profile carries no replicate sample, so it produces no comparison."""
    profile = [
        [Observable(name="p", kind="profile", unit="A", values=[1.0], index=[1])] for _ in range(2)
    ]
    aggregates = aggregate_observables(profile)

    assert compare_observables({"a": aggregates, "b": aggregates}) == []


def test_runner_executes_rg2_end_to_end(tmp_path: Path) -> None:
    """run_analysis computes, persists and aggregates a contract plugin."""
    _, aggregate = _run_condition("A", 1.0, tmp_path)

    assert isinstance(aggregate, ConditionArtifact)
    assert aggregate.replicates == [1, 2, 3]
    observable = ObservableAggregate.model_validate(aggregate.payload["observables"][0])
    assert observable.name == "protein"
    assert observable.unit == "A"
    assert observable.n_replicates == 3
    # Rg of a unit cross scaled by s is exactly s.
    assert observable.replicate_values == pytest.approx([1.5, 2.0, 2.5], abs=1e-5)
    assert observable.mean == pytest.approx(2.0, abs=1e-5)

    replicate_dir = tmp_path / "analysis" / "A" / "rg2" / "run_1"
    assert (replicate_dir / "result.json").exists()
    assert (replicate_dir / "observables.npz").exists()
    identity = aggregate.provenance["identity"]
    assert identity["plugin"] == "rg2"
    assert identity["polyzymd_version"]
    assert identity["config_hash"]


def test_runner_reuses_a_replicate_whose_identity_matches(tmp_path: Path) -> None:
    """A second run reuses the cached replicate instead of recomputing it."""
    _run_condition("A", 1.0, tmp_path)
    marker = tmp_path / "analysis" / "A" / "rg2" / "run_1" / "observables.npz"
    stamp = marker.stat().st_mtime_ns

    _run_condition("A", 1.0, tmp_path)

    assert marker.stat().st_mtime_ns == stamp


def test_runner_recomputes_when_settings_change(tmp_path: Path) -> None:
    """A different settings fingerprint invalidates the cached replicate."""
    condition, _ = _run_condition("A", 1.0, tmp_path)
    analysis = _stub_analysis_class(1.0)()
    lifecycle = AnalysisLifecycle(analysis)

    lifecycle.run_analysis(
        condition,
        RgSettings(runs=[{"label": "everything", "selection": "all"}]),
        "0ns",
        tmp_path / "analysis" / "A" / "rg2",
    )

    artifact = analysis._load_aggregated_result(tmp_path / "analysis" / "A" / "rg2" / "aggregated")
    assert artifact.payload["observables"][0]["name"] == "everything"


def test_runner_compares_two_conditions(tmp_path: Path) -> None:
    """The adapter produces a comparison artifact and a readable report."""
    control, control_aggregate = _run_condition("A", 1.0, tmp_path)
    treatment, treatment_aggregate = _run_condition("B", 3.0, tmp_path)
    analysis = _stub_analysis_class(1.0)()

    comparison = analysis.compare(
        ComparisonContext(
            name="project",
            conditions=[control, treatment],
            excluded_conditions=[],
            control_label="A",
            analysis_dirs={
                "A": tmp_path / "analysis" / "A" / "rg2",
                "B": tmp_path / "analysis" / "B" / "rg2",
            },
            results_dir=tmp_path / "results",
            equilibration="0ns",
            settings=RgSettings(runs=[{"label": "protein", "selection": "all"}]),
            aggregated_results={"A": control_aggregate, "B": treatment_aggregate},
        )
    )

    assert isinstance(comparison, ComparisonArtifact)
    entry = comparison.payload["comparisons"][0]
    assert entry["control"] == "A"
    assert entry["condition"] == "B"
    assert entry["delta"] == pytest.approx(2.0, abs=1e-5)
    assert entry["correction"] == "benjamini_hochberg"
    assert entry["significant"] is True

    report = analysis.format(comparison)
    assert "p_adj" in report
    assert "mean 2 A" in report


def test_contract_scaffold_renders_and_imports(tmp_path: Path) -> None:
    """The contract scaffold produces an importable plugin and a test file."""
    from polyzymd.cli.scaffold import generate_scaffold

    created = generate_scaffold("probe_contract", tmp_path, style="contract")
    plugin_path = tmp_path / "src" / "polyzymd" / "analyses" / "probe_contract.py"
    assert set(created) == {
        plugin_path,
        tmp_path / "tests" / "analyses" / "plugins" / "test_probe_contract.py",
    }
    assert len(plugin_path.read_text().splitlines()) < 60

    spec = importlib.util.spec_from_file_location("probe_contract_scaffold", plugin_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    assert module.ProbeContractAnalysis.name == "probe_contract"
    assert module.ProbeContractAnalysis.Settings is module.ProbeContractSettings
    observables = module.ProbeContract().compute(
        _universe(scale=2.0, n_frames=3),
        types.SimpleNamespace(start=0, stop=3, step=1, frames=None),
        module.ProbeContractSettings(),
    )
    assert reduce_observable(observables[0]).value == pytest.approx(2.0, abs=1e-5)
