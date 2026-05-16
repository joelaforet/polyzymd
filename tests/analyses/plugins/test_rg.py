"""Tests for the MDAnalysis-native Rg analysis plugin."""

from __future__ import annotations

import sys
from datetime import datetime
from types import ModuleType, SimpleNamespace
from typing import Any

import numpy as np
import pytest

from polyzymd.analyses._analysis_lifecycle import AnalysisLifecycle
from polyzymd.analyses.base import Analysis, PlotContext
from polyzymd.analyses.discovery import get_analysis
from polyzymd.analyses.mda import ArtifactStore, ConditionArtifact, ReplicateArtifact
from polyzymd.analyses.rg import RgAnalysis, RgRunSettings, RgSettings
from polyzymd.analyses.rg._comparison_results import (
    RgComparisonResult,
    RgConditionSummary,
    RgRunSummary,
)
from polyzymd.analyses.rg._mda import (
    RG_FRAGMENT_POLICY_WARNING,
    RG_PBC_POLICY_WARNING,
    compute_rg_run,
)
from polyzymd.analyses.rg._plotters import _load_condition_aggregated, _load_replicate_timeseries
from polyzymd.analyses.rg._results import (
    RgAggregatedResult,
    RgResult,
    RgRunAggregatedResult,
)
from polyzymd.analyses.shared.config_hash import settings_fingerprint
from polyzymd.config.comparison import PlotSettings
from tests._support.analysis_testkit import (
    make_aggregate_context,
    make_comparison_context,
    make_condition,
    make_replicate_context,
)


class _FakeTrajectory:
    """Small trajectory object that records the current frame."""

    def __init__(self, universe: _FakeUniverse, n_frames: int = 5, dt: float = 10.0) -> None:
        self.universe = universe
        self.n_frames = n_frames
        self.dt = dt

    def __len__(self) -> int:
        """Return the number of frames."""

        return self.n_frames

    def __getitem__(self, index: int) -> SimpleNamespace:
        """Select one frame and expose frame/time metadata."""

        self.universe.current_frame = int(index)
        return SimpleNamespace(frame=int(index), time=float(index) * self.dt)


class _FakeAtomGroup:
    """Atom group with deterministic frame-dependent Rg values."""

    def __init__(
        self,
        universe: _FakeUniverse,
        values: list[float],
        *,
        masses: float = 1.0,
        fragments: list[_FakeAtomGroup] | None = None,
    ) -> None:
        self.universe = universe
        self.values = list(values)
        self._mass = float(masses)
        self._fragments = fragments
        self.indices = np.arange(len(values) or 1, dtype=np.int64)

    def __len__(self) -> int:
        """Return a nonzero atom count when values exist."""

        return len(self.values)

    @property
    def fragments(self) -> list[_FakeAtomGroup]:
        """Return static topology fragments."""

        return list(self._fragments or [])

    def radius_of_gyration(self) -> float:
        """Return the value assigned to the current trajectory frame."""

        index = min(self.universe.current_frame, len(self.values) - 1)
        return float(self.values[index])

    def total_mass(self) -> float:
        """Return a static fragment mass."""

        return self._mass


class _FakeUniverse:
    """Universe with configurable selections for Rg tests."""

    def __init__(self, n_frames: int = 5, dt: float = 10.0) -> None:
        self.current_frame = 0
        self.trajectory = _FakeTrajectory(self, n_frames=n_frames, dt=dt)
        self.selection_map: dict[str, _FakeAtomGroup] = {}

    def select_atoms(self, selection: str) -> _FakeAtomGroup:
        """Return the configured atom group for a selection."""

        return self.selection_map[selection]


class _FakeAnalysisBase:
    """Minimal AnalysisBase replacement that calls ``_single_frame`` per frame."""

    def __init__(self, trajectory: _FakeTrajectory) -> None:
        self._trajectory = trajectory
        self.results = SimpleNamespace()

    def run(self, start: int = 0, stop: int | None = None, step: int = 1, **kwargs: Any) -> Any:
        """Run frame iteration using the AnalysisBase lifecycle shape."""

        if kwargs:
            raise ValueError(f"Unexpected backend kwargs: {sorted(kwargs)}")
        if stop is None:
            stop = len(self._trajectory)
        self.frames = np.asarray(list(range(start, stop, step)), dtype=np.int64)
        self.times = self.frames.astype(np.float64) * float(self._trajectory.dt)
        self._prepare()
        for frame in self.frames:
            self._trajectory[int(frame)]
            self._single_frame()
        self._conclude()
        return self


@pytest.fixture
def fake_mdanalysis(monkeypatch: pytest.MonkeyPatch) -> None:
    """Install an import-light fake MDAnalysis AnalysisBase module."""

    mda_module = ModuleType("MDAnalysis")
    mda_module.__version__ = "test-mda"
    analysis_module = ModuleType("MDAnalysis.analysis")
    base_module = ModuleType("MDAnalysis.analysis.base")
    base_module.AnalysisBase = _FakeAnalysisBase
    monkeypatch.setitem(sys.modules, "MDAnalysis", mda_module)
    monkeypatch.setitem(sys.modules, "MDAnalysis.analysis", analysis_module)
    monkeypatch.setitem(sys.modules, "MDAnalysis.analysis.base", base_module)


def _settings() -> RgSettings:
    """Return two-run Rg settings for artifact tests."""

    return RgSettings(
        runs=[
            RgRunSettings(label="protein_rg", selection="protein"),
            RgRunSettings(label="polymer/frags", selection="polymer", calculation_mode="fragments"),
        ]
    )


def _make_comparison_result() -> RgComparisonResult:
    """Create a compact comparison result for plotting and formatting tests."""

    return RgComparisonResult(
        metric="mean_rg",
        name="rg_compare",
        n_runs=1,
        run_labels=["protein_rg"],
        control_label="Control",
        conditions=[
            RgConditionSummary(
                label="Control",
                config_path="/fake/control.yaml",
                n_replicates=2,
                run_summaries=[
                    RgRunSummary(
                        label="protein_rg",
                        selection="protein",
                        mean_rg=10.5,
                        sem_rg=0.5,
                        per_replicate_means=[10.0, 11.0],
                        replicates=[1, 2],
                        n_replicates=2,
                    )
                ],
            )
        ],
        pairwise_comparisons=[],
        anova_by_run=None,
        ranking_by_run={"protein_rg": ["Control"]},
        equilibration_time="0ns",
        created_at=datetime.now(),
        polyzymd_version="test",
    )


def _replicate_artifact(
    root: Any,
    *,
    condition_label: str,
    replicate: int,
    settings: RgSettings,
    run_label: str = "protein_rg",
    rg_values: list[float] | None = None,
    skipped_runs: list[dict[str, Any]] | None = None,
) -> ReplicateArtifact:
    """Write one canonical Rg replicate artifact with a sidecar."""

    rg_array = np.asarray(rg_values or [10.0 + replicate, 11.0 + replicate], dtype=np.float64)
    run_dir = root / f"run_{replicate}"
    store = ArtifactStore(run_dir)
    sidecar = store.write_npz_sidecar(
        f"sidecars/rg_{replicate}_{run_label}.npz",
        rg_values=rg_array,
        reduced_rg_values=rg_array,
        time_ns=np.arange(rg_array.size, dtype=np.float64),
        frames=np.arange(rg_array.size, dtype=np.int64),
        raw_timestep_ps=np.asarray(10.0, dtype=np.float64),
        frame_stride=np.asarray(1, dtype=np.int64),
        effective_timestep_ps=np.asarray(10.0, dtype=np.float64),
        metadata={"run_label": run_label},
    )
    mean_rg = float(np.mean(rg_array))
    run_payload = {
        "config_hash": "hash123",
        "polyzymd_version": "test",
        "replicate": replicate,
        "equilibration_time": 0.0,
        "equilibration_unit": "ns",
        "selection_string": "protein",
        "run_label": run_label,
        "selection": "protein",
        "calculation_mode": "selection",
        "fragment_weighting": None,
        "mean_rg": mean_rg,
        "std_rg": float(np.std(rg_array, ddof=0)),
        "median_rg": float(np.median(rg_array)),
        "min_rg": float(np.min(rg_array)),
        "max_rg": float(np.max(rg_array)),
        "final_rg": float(rg_array[-1]),
        "sem_rg": None,
        "correlation_time_unit": None,
        "statistical_inefficiency": None,
        "autocorrelation_warning": None,
        "n_frames_total": int(rg_array.size),
        "n_frames_used": int(rg_array.size),
        "npz_path": sidecar.path,
        "sidecar": sidecar.model_dump(mode="json"),
        "time_unit": "ns",
        "timestep_ps": 10.0,
        "raw_timestep_ps": 10.0,
        "frame_stride": 1,
    }
    metrics = {f"{run_label}.mean_rg": mean_rg}
    artifact = ReplicateArtifact(
        analysis_name="rg",
        condition_label=condition_label,
        replicate=replicate,
        payload={
            "runs": [run_payload],
            "skipped_runs": skipped_runs or [],
            "metrics": metrics,
            "replicate_metrics": metrics,
            "n_frames_total": int(rg_array.size),
            "n_frames_used": int(rg_array.size),
        },
        sidecars=[sidecar],
        provenance={"pbc_policy": {"coordinates": "as_loaded"}},
        metadata={
            "result_kind": "rg_mda_replicate",
            "settings_fingerprint": settings_fingerprint(settings),
            "config_hash": "hash123",
            "polyzymd_version": "test",
            "equilibration_time": 0.0,
            "equilibration_unit": "ns",
            "selection_string": "protein",
        },
        warnings=[RG_PBC_POLICY_WARNING],
    )
    store.write_replicate_result(artifact)
    return artifact


def test_rg_plugin_discovered() -> None:
    """Rg plugin should be auto-discovered and MDA-backed."""

    cls = get_analysis("rg")

    assert cls is RgAnalysis
    assert RgAnalysis.name == "rg"
    assert RgAnalysis.min_replicates == 1
    assert RgAnalysis.build_runner is Analysis.build_runner
    assert RgAnalysis.ReplicateResultClass is None


def test_rg_settings_validate_fragment_options() -> None:
    """Rg settings should preserve fragment-mode validation."""

    settings = RgRunSettings(label="frags", selection="segid C", calculation_mode="fragments")

    assert settings.fragment_weighting == "equal"
    with pytest.raises(ValueError, match="fragment_weighting is only meaningful"):
        RgRunSettings(
            label="bad",
            selection="protein",
            calculation_mode="selection",
            fragment_weighting="mass",
        )
    with pytest.raises(ValueError, match="At least one Rg run"):
        RgSettings(runs=[])


def test_compute_rg_run_selection_mode_uses_analysisbase(fake_mdanalysis: None) -> None:
    """Custom Rg AnalysisBase should collect frame/time arrays in selection mode."""

    universe = _FakeUniverse(n_frames=4, dt=20.0)
    universe.selection_map["protein"] = _FakeAtomGroup(universe, [1.0, 2.0, 3.0, 4.0])

    payload = compute_rg_run(
        universe=universe,
        run=RgRunSettings(label="protein_rg", selection="protein"),
        replicate=1,
        start=1,
        stop=4,
        step=2,
        timestep_ps=20.0,
    )

    assert payload.rg_values.tolist() == [2.0, 4.0]
    assert payload.frames.tolist() == [1, 3]
    assert payload.time_ns.tolist() == [0.02, 0.06]
    assert payload.effective_timestep_ps == pytest.approx(40.0)


def test_compute_rg_run_fragment_mode_records_distributions(fake_mdanalysis: None) -> None:
    """Fragment mode should reduce fragment Rg values and retain masses/counts."""

    universe = _FakeUniverse(n_frames=3, dt=10.0)
    frag_a = _FakeAtomGroup(universe, [1.0, 2.0, 3.0], masses=2.0)
    frag_b = _FakeAtomGroup(universe, [3.0, 4.0, 5.0], masses=1.0)
    universe.selection_map["polymer"] = _FakeAtomGroup(
        universe,
        [0.0, 0.0, 0.0],
        fragments=[frag_a, frag_b],
    )

    payload = compute_rg_run(
        universe=universe,
        run=RgRunSettings(
            label="polymer_frags",
            selection="polymer",
            calculation_mode="fragments",
            fragment_weighting="mass",
        ),
        replicate=1,
        start=0,
        stop=3,
        step=1,
        timestep_ps=10.0,
    )

    assert payload.rg_values.tolist() == pytest.approx([5.0 / 3.0, 8.0 / 3.0, 11.0 / 3.0])
    assert payload.fragment_counts_per_frame.tolist() == [2, 2, 2]
    assert payload.fragment_masses.tolist() == [2.0, 1.0]
    assert payload.fragment_rg_values.tolist() == [1.0, 3.0, 2.0, 4.0, 3.0, 5.0]


def test_run_replicate_persists_artifact_and_empty_selection_skip(
    fake_mdanalysis: None,
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """MDA lifecycle should persist a ReplicateArtifact and explicit skip provenance."""

    analysis = RgAnalysis()
    settings = RgSettings(
        runs=[
            RgRunSettings(label="protein_rg", selection="protein"),
            RgRunSettings(label="missing_rg", selection="missing"),
        ]
    )
    condition = make_condition(label="Control", replicates=(1,))
    ctx = make_replicate_context(
        condition=condition,
        replicate=1,
        output_dir=tmp_path / "run_1",
        settings=settings,
        equilibration="0ns",
    )
    universe = _FakeUniverse(n_frames=3, dt=10.0)
    universe.selection_map["protein"] = _FakeAtomGroup(universe, [10.0, 11.0, 12.0])
    universe.selection_map["missing"] = _FakeAtomGroup(universe, [])

    class _Provider:
        def __init__(self, *_args: Any, **_kwargs: Any) -> None:
            pass

        @classmethod
        def from_config(cls, *_args: Any, **_kwargs: Any) -> _Provider:
            return cls()

        def load_universe(self, _replicate: int) -> _FakeUniverse:
            return universe

        def provenance_for(self, _replicate: int) -> dict[str, Any]:
            return {"warnings": []}

    class _Loader:
        def __init__(self, *_args: Any, **_kwargs: Any) -> None:
            pass

        def get_timestep(self, _replicate: int, unit: str = "ps") -> float:
            """Return a fixed timestep for frame-window resolution."""

            del unit
            return 10.0

    monkeypatch.setattr(analysis, "_mda_universe_provider_factory", lambda: _Provider)
    monkeypatch.setattr(analysis, "_trajectory_loader_factory", lambda: _Loader)
    monkeypatch.setattr("polyzymd.analyses.rg._mda.compute_config_hash", lambda _config: "hash123")

    artifact = AnalysisLifecycle(analysis).run_replicate_once(
        condition,
        settings,
        "0ns",
        ctx.output_dir,
        1,
        recompute=True,
    )

    assert isinstance(artifact, ReplicateArtifact)
    assert [run["run_label"] for run in artifact.payload["runs"]] == ["protein_rg"]
    assert artifact.payload["skipped_runs"][0]["run_label"] == "missing_rg"
    assert RG_PBC_POLICY_WARNING in artifact.warnings
    persisted = ArtifactStore(ctx.output_dir).read_replicate_result("result.json")
    assert persisted.replicate == 1
    assert persisted.sidecars[0].path.startswith("sidecars/rg_000_protein_rg_")


def test_aggregate_rejects_legacy_inputs(tmp_path) -> None:
    """Rg aggregation should reject legacy RgResult inputs with recompute guidance."""

    analysis = RgAnalysis()
    settings = RgSettings(runs=[RgRunSettings(label="protein_rg", selection="protein")])
    ctx = make_aggregate_context(
        condition=make_condition(label="Control", replicates=(1,)),
        replicates=(1,),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="0ns",
    )
    legacy = RgResult(
        config_hash="hash",
        polyzymd_version="test",
        replicate=1,
        equilibration_time=0.0,
        equilibration_unit="ns",
        selection_string="protein",
        run_results=[],
        n_frames_total=1,
        n_frames_used=1,
    )

    with pytest.raises(TypeError, match="Legacy Rg replicate caches are incompatible"):
        analysis.aggregate(ctx, [legacy])


def test_aggregate_artifacts_and_compare_via_adapter(tmp_path) -> None:
    """Rg should aggregate ReplicateArtifact inputs and compare condition artifacts."""

    analysis = RgAnalysis()
    settings = RgSettings(runs=[RgRunSettings(label="protein_rg", selection="protein")])
    control = make_condition(label="Control", replicates=(1, 2))
    treated = make_condition(label="Treated", replicates=(1, 2))
    control_dir = tmp_path / "control"
    treated_dir = tmp_path / "treated"
    control_artifacts = [
        _replicate_artifact(control_dir, condition_label="Control", replicate=1, settings=settings),
        _replicate_artifact(control_dir, condition_label="Control", replicate=2, settings=settings),
    ]
    treated_artifacts = [
        _replicate_artifact(
            treated_dir,
            condition_label="Treated",
            replicate=1,
            settings=settings,
            rg_values=[8.0, 9.0],
        ),
        _replicate_artifact(
            treated_dir,
            condition_label="Treated",
            replicate=2,
            settings=settings,
            rg_values=[9.0, 10.0],
        ),
    ]
    control_ctx = make_aggregate_context(
        condition=control,
        replicates=(1, 2),
        output_dir=control_dir / "aggregated",
        settings=settings,
        equilibration="0ns",
    )
    treated_ctx = make_aggregate_context(
        condition=treated,
        replicates=(1, 2),
        output_dir=treated_dir / "aggregated",
        settings=settings,
        equilibration="0ns",
    )

    control_agg = analysis.aggregate(control_ctx, control_artifacts)
    treated_agg = analysis.aggregate(treated_ctx, treated_artifacts)

    assert isinstance(control_agg, ConditionArtifact)
    assert control_agg.payload["runs"][0]["overall_mean"] == pytest.approx(12.0)
    assert control_agg.payload["runs"][0]["reduced_histogram_edges"] is not None
    compare_ctx = make_comparison_context(
        name="rg_compare",
        conditions=[control, treated],
        analysis_dirs={"Control": control_dir, "Treated": treated_dir},
        results_dir=tmp_path / "comparison",
        settings=settings,
        control_label="Control",
        equilibration="0ns",
        aggregated_results={"Control": control_agg, "Treated": treated_agg},
    )
    comparison = analysis.compare(compare_ctx)
    assert comparison is not None
    assert comparison.ranking_by_run["protein_rg"] == ["Treated", "Control"]


def test_aggregate_honors_skip_only_run(tmp_path) -> None:
    """Aggregation should omit runs only when every missing entry has skip provenance."""

    analysis = RgAnalysis()
    settings = RgSettings(
        runs=[
            RgRunSettings(label="protein_rg", selection="protein"),
            RgRunSettings(label="polymer_rg", selection="polymer"),
        ]
    )
    condition = make_condition(label="No Polymer", replicates=(1, 2))

    def skipped(rep: int) -> list[dict[str, Any]]:
        """Return skip provenance for one replicate."""

        return [
            {
                "run_label": "polymer_rg",
                "selection": "polymer",
                "replicate": rep,
                "reason": "selection matched no atoms",
                "reason_code": "empty_selection",
            }
        ]

    artifacts = [
        _replicate_artifact(
            tmp_path,
            condition_label="No Polymer",
            replicate=rep,
            settings=settings,
            skipped_runs=skipped(rep),
        )
        for rep in (1, 2)
    ]
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="0ns",
    )

    aggregate = analysis.aggregate(ctx, artifacts)

    assert [run["run_label"] for run in aggregate.payload["runs"]] == ["protein_rg"]
    assert len(aggregate.payload["skipped_runs"]) == 2


def test_aggregate_rejects_missing_artifact_sidecar(tmp_path) -> None:
    """Aggregation should validate canonical NPZ sidecars."""

    analysis = RgAnalysis()
    settings = RgSettings(runs=[RgRunSettings(label="protein_rg", selection="protein")])
    artifact = _replicate_artifact(
        tmp_path,
        condition_label="Control",
        replicate=1,
        settings=settings,
    )
    sidecar_path = ArtifactStore(tmp_path / "run_1").resolve_sidecar(artifact.sidecars[0])
    sidecar_path.unlink()
    ctx = make_aggregate_context(
        condition=make_condition(label="Control", replicates=(1,)),
        replicates=(1,),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="0ns",
    )

    with pytest.raises(Exception, match="Missing sidecar"):
        analysis.aggregate(ctx, [artifact])


def test_plotters_load_canonical_artifacts_only(tmp_path) -> None:
    """Rg plot helpers should read replicate and condition artifacts, not legacy names."""

    settings = RgSettings(runs=[RgRunSettings(label="protein_rg", selection="protein")])
    analysis_dir = tmp_path / "control"
    artifacts = [
        _replicate_artifact(
            analysis_dir, condition_label="Control", replicate=1, settings=settings
        ),
        _replicate_artifact(
            analysis_dir, condition_label="Control", replicate=2, settings=settings
        ),
    ]
    analysis = RgAnalysis()
    ctx = make_aggregate_context(
        condition=make_condition(label="Control", replicates=(1, 2)),
        replicates=(1, 2),
        output_dir=analysis_dir / "aggregated",
        settings=settings,
        equilibration="0ns",
    )
    analysis.aggregate(ctx, artifacts)

    time_ns, matrix = _load_replicate_timeseries(analysis_dir, "protein_rg", [1, 2])
    aggregate_payload = _load_condition_aggregated(analysis_dir)

    assert time_ns.tolist() == [0.0, 1.0]
    assert matrix.shape == (2, 2)
    assert aggregate_payload is not None
    assert aggregate_payload["runs"][0]["run_label"] == "protein_rg"


def test_plot_rejects_legacy_aggregated_result_from_disk(tmp_path) -> None:
    """plot should fail loudly when aggregated Rg cache lacks artifact identity."""

    analysis = RgAnalysis()
    settings = RgSettings(runs=[RgRunSettings(label="protein_rg", selection="protein")])
    analysis_dir = tmp_path / "analysis" / "control" / "rg"
    aggregated_dir = analysis_dir / "aggregated"
    aggregated_dir.mkdir(parents=True)
    legacy_run = RgRunAggregatedResult(
        config_hash="hash",
        polyzymd_version="test",
        replicate=None,
        equilibration_time=0.0,
        equilibration_unit="ns",
        selection_string="protein",
        replicates=[1, 2],
        n_replicates=2,
        run_label="protein_rg",
        selection="protein",
        overall_mean=1.0,
        overall_sem=0.0,
        overall_median=1.0,
        per_replicate_means=[1.0, 1.0],
        per_replicate_stds=[0.0, 0.0],
        per_replicate_medians=[1.0, 1.0],
    )
    RgAggregatedResult(
        config_hash="hash",
        polyzymd_version="test",
        replicate=None,
        equilibration_time=0.0,
        equilibration_unit="ns",
        selection_string="protein",
        replicates=[1, 2],
        n_replicates=2,
        run_results=[legacy_run],
    ).save(aggregated_dir / "result.json")
    results_dir = tmp_path / "comparison"
    results_dir.mkdir()
    _make_comparison_result().save(results_dir / "result.json")
    plot_ctx = PlotContext(
        conditions=[make_condition(label="Control", replicates=(1, 2))],
        analysis_dirs={"Control": analysis_dir},
        results_dir=results_dir,
        output_dir=tmp_path / "figures",
        settings=settings,
        plot_settings=PlotSettings(),
        equilibration="0ns",
    )

    with pytest.raises(ValueError, match="missing a settings fingerprint"):
        analysis.plot(plot_ctx)


def test_condition_artifact_records_pbc_and_fragment_assumptions(tmp_path) -> None:
    """Aggregates should preserve explicit PBC and fragment-topology provenance."""

    settings = _settings()
    artifact = _replicate_artifact(
        tmp_path,
        condition_label="Control",
        replicate=1,
        settings=settings,
        skipped_runs=[
            {
                "run_label": "polymer/frags",
                "selection": "polymer",
                "replicate": 1,
                "reason": "selection matched no atoms",
                "reason_code": "empty_selection",
            }
        ],
    )
    artifact.warnings.append(RG_FRAGMENT_POLICY_WARNING)
    ctx = make_aggregate_context(
        condition=make_condition(label="Control", replicates=(1,)),
        replicates=(1,),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="0ns",
    )

    aggregate = RgAnalysis().aggregate(ctx, [artifact])

    assert aggregate.provenance["pbc_policy"]["coordinates"] == "as_loaded"
    assert RG_PBC_POLICY_WARNING in aggregate.warnings
    assert RG_FRAGMENT_POLICY_WARNING in aggregate.warnings
