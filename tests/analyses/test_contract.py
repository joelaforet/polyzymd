"""Known-answer tests for the observable contract and its lifecycle adapter."""

from __future__ import annotations

import importlib.util
import math
import types
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pytest

from polyzymd.analyses._framework.contexts import ComparisonContext, Condition
from polyzymd.analyses.contract import (
    AnalysisProtocol,
    Observable,
    ObservableAggregate,
    aggregate_observables,
    compare_observables,
    reduce_observable,
    reduce_replicate,
)
from polyzymd.analyses.contract_runner import contract_analysis
from polyzymd.analyses.exceptions import PluginContractError
from polyzymd.analyses.mda.artifacts import ComparisonArtifact, ConditionArtifact
from polyzymd.analyses.rg_contract import Rg2Analysis, RgContract, RgSettings
from tests.analyses.conftest import make_simulation_config, make_synthetic_universe

RG_SETTINGS = RgSettings(runs=[{"label": "protein", "selection": "all"}])
OTHER_INPUTS = ({"path": "/tmp/run.xtc", "format": "xtc", "size_bytes": 99, "mtime_ns": 7},)


def _series(values: Any, kind: str = "mean_of_timeseries", name: str = "x") -> Observable:
    """Build one observable of the requested kind."""
    return Observable(name=name, kind=kind, unit="A", values=values)


def _universes_from(base: float) -> Callable[[int], Any]:
    """Replicate factory whose radius of gyration rises by 0.5 per replicate."""
    return lambda replicate: make_synthetic_universe(scale=base + 0.5 * replicate)


_scaled_universes = _universes_from(1.0)


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


def test_the_scaffold_unit_placeholder_is_rejected() -> None:
    """An author must replace the placeholder unit with a real one."""
    with pytest.raises(ValueError, match="has not stated its unit"):
        Observable(name="x", kind="mean_of_timeseries", unit="TODO", values=[1.0])


def test_fluctuation_reduces_to_the_sample_standard_deviation() -> None:
    """A fluctuation observable reduces to the standard deviation of its series."""
    estimate = reduce_observable(_series([1.0, 2.0, 3.0, 4.0, 5.0], kind="fluctuation"))

    assert estimate.value == pytest.approx(math.sqrt(2.5))
    assert estimate.n_frames == 5


def test_single_frame_fluctuation_has_no_estimate() -> None:
    """One frame cannot fluctuate, so the replicate contributes no value."""
    assert reduce_observable(_series([2.0], kind="fluctuation")).value is None


def test_too_few_estimable_replicates_is_rejected() -> None:
    """Aggregation refuses to report a mean over fewer than two estimates."""
    replicates = [
        [_series([1.0, 2.0], kind="fluctuation")],
        [_series([3.0], kind="fluctuation")],
    ]

    with pytest.raises(PluginContractError, match="estimable replicate"):
        aggregate_observables(replicates)


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
    assert aggregate.ci_method == "student_t"
    assert aggregate.coverage == pytest.approx(0.95)
    assert aggregate.mean is None


def _profile(name: str, values: Any, **fields: Any) -> Observable:
    """Build one profile observable indexed 1..n."""
    return Observable(
        name=name,
        kind="profile",
        unit="A",
        values=values,
        index=list(range(1, len(values) + 1)),
        **fields,
    )


def test_a_profile_can_declare_a_scalar_mean_over_its_index() -> None:
    """The declared reduction adds one scalar observable beside the profile."""
    replicates = [
        [_profile("rmsf", [1.0, 3.0], reduce="mean_over_index", reduced_kind="fluctuation")],
        [_profile("rmsf", [3.0, 7.0], reduce="mean_over_index", reduced_kind="fluctuation")],
    ]

    aggregates = {aggregate.name: aggregate for aggregate in aggregate_observables(replicates)}

    assert sorted(aggregates) == ["rmsf", "rmsf_mean"]
    assert aggregates["rmsf_mean"].kind == "fluctuation"
    assert aggregates["rmsf_mean"].unit == "A"
    assert aggregates["rmsf_mean"].replicate_values == pytest.approx([2.0, 5.0])
    assert aggregates["rmsf_mean"].mean == pytest.approx(3.5)


def test_a_profile_can_declare_a_total_over_its_index() -> None:
    """A sum is named for the operation and is not called a mean."""
    replicates = [
        [_profile("area", [1.0, 3.0], reduce="sum_over_index")],
        [_profile("area", [2.0, 6.0], reduce="sum_over_index")],
    ]

    aggregates = {aggregate.name: aggregate for aggregate in aggregate_observables(replicates)}

    assert sorted(aggregates) == ["area", "area_total"]
    assert aggregates["area_total"].kind == "mean_of_timeseries"
    assert aggregates["area_total"].replicate_values == pytest.approx([4.0, 8.0])


def test_a_declared_scalar_keeps_the_frame_count_and_drops_the_diagnostics() -> None:
    """The scalar counts the frames the profile came from, not its indices."""
    profile = _profile(
        "rmsf",
        [1.0, 2.0, 3.0, 4.0],
        n_frames=500,
        reduce="mean_over_index",
        reduced_kind="fluctuation",
    )

    estimates = reduce_replicate([profile])

    assert [estimate.name for estimate in estimates] == ["rmsf", "rmsf_mean"]
    assert estimates[1].value == pytest.approx(2.5)
    assert estimates[1].n_frames == 500
    assert estimates[1].statistical_inefficiency is None
    assert estimates[1].n_eff is None


def test_a_reduction_on_a_non_profile_is_rejected() -> None:
    """Only a profile has an index to reduce over."""
    with pytest.raises(ValueError, match="kind is not 'profile'"):
        Observable(
            name="x", kind="mean_of_timeseries", unit="A", values=[1.0], reduce="mean_over_index"
        )


def test_a_profile_cannot_reduce_to_a_profile() -> None:
    """The scalar a profile declares is a scalar."""
    with pytest.raises(ValueError, match="reduce a profile to a profile"):
        _profile("x", [1.0], reduce="mean_over_index", reduced_kind="profile")


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
        label: aggregate_observables(
            [[_series([v], name="a"), _series([v], name="b")] for v in values]
        )
        for label, values in (
            ("control", (1.0, 1.1, 0.9)),
            ("low", (1.02, 1.11, 0.93)),
            ("high", (5.0, 5.1, 4.9)),
        )
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
    assert all(c.significant and c.testable for c in high)
    assert all(c.delta == pytest.approx(4.0, abs=1e-6) for c in high)


def test_an_untested_observable_is_reported_but_never_tested() -> None:
    """tested=False keeps a dependent quantity out of the tests and the family.

    The untested observable is still aggregated with its mean and SEM, it
    produces no pairwise row, and the adjusted p-values of the tested
    observables are what they would be if it had never been reported.
    """

    def condition(values: tuple[float, ...], *, with_dependent: bool) -> list[Any]:
        replicates = []
        for value in values:
            observables = [_series([value], name="a"), _series([value], name="b")]
            if with_dependent:
                observables.append(
                    Observable(
                        name="dependent",
                        kind="mean_of_timeseries",
                        unit="A",
                        values=[10.0 - value],
                        tested=False,
                    )
                )
            replicates.append(observables)
        return aggregate_observables(replicates)

    samples = {"control": (1.0, 1.1, 0.9), "high": (5.0, 5.1, 4.9)}
    with_dependent = {
        label: condition(values, with_dependent=True) for label, values in samples.items()
    }
    without = {label: condition(values, with_dependent=False) for label, values in samples.items()}

    dependent = next(agg for agg in with_dependent["control"] if agg.name == "dependent")
    assert dependent.n_replicates == 3
    assert dependent.mean == pytest.approx(9.0, abs=1e-6)
    assert dependent.sem == pytest.approx(0.1 / 3**0.5, abs=1e-6)
    assert dependent.tested is False

    comparisons = compare_observables(with_dependent, control_label="control")

    assert [c.name for c in comparisons] == ["a", "b"]
    baseline = compare_observables(without, control_label="control")
    assert [c.p_adjusted for c in comparisons] == pytest.approx(
        [c.p_adjusted for c in baseline], abs=1e-12
    )


def test_unknown_control_label_is_rejected() -> None:
    """A control that names no compared condition is an error, not a fallback."""
    aggregates = aggregate_observables([[_series([1.0])], [_series([2.0])]])

    with pytest.raises(PluginContractError, match="is not among the compared conditions"):
        compare_observables({"a": aggregates, "b": aggregates}, control_label="A")


def test_a_single_replicate_condition_is_not_testable() -> None:
    """One replicate gives no test, and the pair says so instead of no difference."""
    control = aggregate_observables([[_series([1.0])], [_series([1.2])], [_series([0.9])]])
    thin = aggregate_observables([[_series([5.0])]])

    comparison = compare_observables({"control": control, "thin": thin})[0]

    assert comparison.testable is False
    assert comparison.note == "single replicate"
    assert comparison.significant is False


def test_profiles_are_not_tested_pairwise() -> None:
    """A profile carries no replicate sample, so it produces no comparison."""
    profile = [
        [Observable(name="p", kind="profile", unit="A", values=[1.0], index=[1])] for _ in range(2)
    ]
    aggregates = aggregate_observables(profile)

    assert compare_observables({"a": aggregates, "b": aggregates}) == []


def test_a_plugin_that_misses_the_protocol_is_named() -> None:
    """contract_analysis enforces the protocol and says what is missing."""

    class Incomplete:
        name = "incomplete"

    with pytest.raises(PluginContractError, match=r"Settings.*compute.*references"):
        contract_analysis(Incomplete)


def test_a_plugin_keeps_its_slurm_resource_hint() -> None:
    """A plugin that needs more memory than the default says so on the class."""
    from polyzymd.analyses.base import SlurmResourceHint
    from polyzymd.analyses.rg_contract import RgContract

    class Hungry(RgContract):
        name = "hungry"
        slurm_resource_hint = SlurmResourceHint(mem="16G")

    assert contract_analysis(Hungry).slurm_resource_hint == SlurmResourceHint(mem="16G")
    assert contract_analysis(RgContract).slurm_resource_hint is None


def test_rg2_satisfies_the_protocol() -> None:
    """The prototype port is an instance of the runtime-checkable protocol."""
    from polyzymd.analyses.rg_contract import RgContract

    assert isinstance(RgContract(), AnalysisProtocol)


def test_runner_executes_rg2_end_to_end(tmp_path: Path, run_contract_analysis: Any) -> None:
    """run_analysis computes, persists and aggregates a contract plugin."""
    aggregate = run_contract_analysis(Rg2Analysis, RG_SETTINGS, _scaled_universes, root=tmp_path)

    assert isinstance(aggregate, ConditionArtifact)
    assert aggregate.replicates == [1, 2, 3]
    observable = ObservableAggregate.model_validate(aggregate.payload["observables"][0])
    assert observable.name == "protein"
    assert observable.unit == "A"
    assert observable.n_replicates == 3
    assert observable.replicate_values == pytest.approx([1.5, 2.0, 2.5], abs=1e-5)
    assert observable.mean == pytest.approx(2.0, abs=1e-5)

    replicate_dir = tmp_path / "analysis" / "A" / "rg2" / "run_1"
    assert (replicate_dir / "result.json").exists()
    assert (replicate_dir / "observables.npz").exists()
    identity = aggregate.provenance["identity"]
    assert identity["plugin"] == "rg2"
    assert identity["polyzymd_version"]
    assert identity["config_hash"]


def test_runner_reuses_a_replicate_whose_identity_matches(
    tmp_path: Path, run_contract_analysis: Any
) -> None:
    """A second run reuses the cached replicate instead of recomputing it."""
    run_contract_analysis(Rg2Analysis, RG_SETTINGS, _scaled_universes, root=tmp_path)
    marker = tmp_path / "analysis" / "A" / "rg2" / "run_1" / "observables.npz"
    stamp = marker.stat().st_mtime_ns

    run_contract_analysis(Rg2Analysis, RG_SETTINGS, _scaled_universes, root=tmp_path)

    assert marker.stat().st_mtime_ns == stamp


def test_runner_recomputes_when_a_shared_version_is_bumped(
    tmp_path: Path, run_contract_analysis: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Bumping a shared-module version invalidates every cached replicate.

    ``plugin_code_hash`` only covers the plugin module, so a fix in a shared
    module has to announce itself through its version constant.
    """
    from polyzymd.analyses import contract_runner

    run_contract_analysis(Rg2Analysis, RG_SETTINGS, _scaled_universes, root=tmp_path)
    marker = tmp_path / "analysis" / "A" / "rg2" / "run_1" / "observables.npz"
    stamp = marker.stat().st_mtime_ns

    monkeypatch.setattr(contract_runner, "ALIGNMENT_VERSION", "99")
    run_contract_analysis(Rg2Analysis, RG_SETTINGS, _scaled_universes, root=tmp_path)

    assert marker.stat().st_mtime_ns != stamp


def test_shared_versions_are_recorded_in_the_identity_block(
    tmp_path: Path, run_contract_analysis: Any
) -> None:
    """The identity block names the shared modules a cached result depends on."""
    from polyzymd.analyses.contract_runner import _shared_versions

    aggregate = run_contract_analysis(Rg2Analysis, RG_SETTINGS, _scaled_universes, root=tmp_path)

    assert aggregate.provenance["identity"]["shared_versions"] == _shared_versions()


def test_runner_recomputes_when_an_input_file_changes(
    tmp_path: Path, run_contract_analysis: Any
) -> None:
    """Extending a trajectory changes its file identity and invalidates the cache."""
    run_contract_analysis(Rg2Analysis, RG_SETTINGS, _scaled_universes, root=tmp_path)
    marker = tmp_path / "analysis" / "A" / "rg2" / "run_1" / "observables.npz"
    stamp = marker.stat().st_mtime_ns

    run_contract_analysis(
        Rg2Analysis, RG_SETTINGS, _scaled_universes, root=tmp_path, inputs=OTHER_INPUTS
    )

    assert marker.stat().st_mtime_ns != stamp


def test_runner_recomputes_when_a_declared_settings_file_changes(
    tmp_path: Path, run_contract_analysis: Any
) -> None:
    """A file the plugin names through identity_files is part of the identity.

    The path does not change, only the contents, which is what happens when a
    reference structure is regenerated in place.
    """
    reference = tmp_path / "reference.pdb"
    reference.write_text("first", encoding="utf-8")

    class _WithFile(RgContract):
        name = "rg_with_file"

        @staticmethod
        def identity_files(settings: Any) -> tuple[Path, ...]:
            del settings
            return (reference,)

    analysis_cls = contract_analysis(_WithFile)
    run_contract_analysis(analysis_cls, RG_SETTINGS, _scaled_universes, root=tmp_path)
    marker = tmp_path / "analysis" / "A" / "rg_with_file" / "run_1" / "observables.npz"
    stamp = marker.stat().st_mtime_ns

    run_contract_analysis(analysis_cls, RG_SETTINGS, _scaled_universes, root=tmp_path)
    assert marker.stat().st_mtime_ns == stamp

    reference.write_text("second, longer contents", encoding="utf-8")
    run_contract_analysis(analysis_cls, RG_SETTINGS, _scaled_universes, root=tmp_path)

    assert marker.stat().st_mtime_ns != stamp


def test_runner_records_a_declared_settings_file_that_is_missing(
    tmp_path: Path, run_contract_analysis: Any
) -> None:
    """A file that is absent is recorded as absent, so creating it invalidates."""

    class _WithMissingFile(RgContract):
        name = "rg_missing_file"

        @staticmethod
        def identity_files(settings: Any) -> tuple[Path, ...]:
            del settings
            return (tmp_path / "not_here.pdb",)

    aggregate = run_contract_analysis(
        contract_analysis(_WithMissingFile), RG_SETTINGS, _scaled_universes, root=tmp_path
    )

    identity = aggregate.provenance["identity"]
    assert identity["settings_files"] == [{"path": str(tmp_path / "not_here.pdb"), "missing": True}]


def test_runner_recomputes_when_the_plugin_source_changes(
    tmp_path: Path, run_contract_analysis: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A fix inside the plugin invalidates every replicate it already wrote.

    This is what keeps a corrected estimator from being averaged with numbers
    the superseded one produced, without a hand-maintained version key in the
    plugin.
    """
    import polyzymd.analyses.contract_runner as runner

    run_contract_analysis(Rg2Analysis, RG_SETTINGS, _scaled_universes, root=tmp_path)
    marker = tmp_path / "analysis" / "A" / "rg2" / "run_1" / "observables.npz"
    stamp = marker.stat().st_mtime_ns

    monkeypatch.setattr(runner, "_code_hash", lambda plugin: "a different plugin")
    run_contract_analysis(Rg2Analysis, RG_SETTINGS, _scaled_universes, root=tmp_path)

    assert marker.stat().st_mtime_ns != stamp


def test_runner_recomputes_when_settings_change(tmp_path: Path, run_contract_analysis: Any) -> None:
    """A different settings fingerprint invalidates the cached replicate."""
    run_contract_analysis(Rg2Analysis, RG_SETTINGS, _scaled_universes, root=tmp_path)

    aggregate = run_contract_analysis(
        Rg2Analysis,
        RgSettings(runs=[{"label": "everything", "selection": "all"}]),
        _scaled_universes,
        root=tmp_path,
    )

    assert aggregate.payload["observables"][0]["name"] == "everything"


def test_runner_compares_two_conditions(tmp_path: Path, run_contract_analysis: Any) -> None:
    """The adapter produces a comparison artifact and a readable report."""
    aggregates = {
        label: run_contract_analysis(
            Rg2Analysis,
            RG_SETTINGS,
            _universes_from(base),
            label=label,
            root=tmp_path,
        )
        for label, base in (("A", 1.0), ("B", 3.0))
    }
    conditions = [
        Condition(
            label=label,
            config_path=tmp_path / f"{label}.yaml",
            replicates=(1, 2, 3),
            sim_config=make_simulation_config(label),
        )
        for label in aggregates
    ]

    comparison = Rg2Analysis().compare(
        ComparisonContext(
            name="project",
            conditions=conditions,
            excluded_conditions=[],
            control_label="A",
            analysis_dirs={label: tmp_path / "analysis" / label / "rg2" for label in aggregates},
            results_dir=tmp_path / "results",
            equilibration="0ns",
            settings=RG_SETTINGS,
            aggregated_results=aggregates,
        )
    )

    assert isinstance(comparison, ComparisonArtifact)
    entry = comparison.payload["comparisons"][0]
    assert entry["control"] == "A"
    assert entry["condition"] == "B"
    assert entry["delta"] == pytest.approx(2.0, abs=1e-5)
    assert entry["correction"] == "benjamini_hochberg"
    assert entry["significant"] is True

    report = Rg2Analysis().format(comparison)
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
    with pytest.raises(ValueError, match="has not stated its unit"):
        module.ProbeContract().compute(
            make_synthetic_universe(scale=2.0, n_frames=3),
            types.SimpleNamespace(start=0, stop=3, step=1, frames=None),
            module.ProbeContractSettings(),
        )


def test_metadata_and_index_label_reach_the_aggregate() -> None:
    """Per-observable provenance survives the reduction and the aggregation."""
    replicates = [
        [
            Observable(
                name="p",
                kind="profile",
                unit="A",
                values=[1.0, 2.0],
                index=[7.0, 8.0],
                index_label="residue index",
                metadata={"chunk_size": 50},
            )
        ]
        for _ in range(2)
    ]

    aggregate = aggregate_observables(replicates)[0]

    assert aggregate.index_label == "residue index"
    assert aggregate.metadata == {"chunk_size": 50}
    assert reduce_observable(replicates[0][0]).metadata == {"chunk_size": 50}


def test_index_label_needs_a_profile() -> None:
    """Labelling an index that does not exist is a contract error."""
    with pytest.raises(ValueError, match="has an index but kind is not 'profile'"):
        Observable(
            name="x",
            kind="mean_of_timeseries",
            unit="A",
            values=[1.0],
            index_label="residue index",
        )
