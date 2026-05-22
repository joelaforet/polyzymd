"""Tests for the hydrogen-bonds analysis plugin scaffold."""

from __future__ import annotations

import sys
import types
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from pydantic import ValidationError

from polyzymd.analyses import get_analysis, list_analyses
from polyzymd.analyses.base import (
    AggregateContext,
    ANOVAResult,
    ComparisonContext,
    ComparisonResult,
    Condition,
    ConditionSummary,
    PairwiseResult,
    PlotContext,
    ReplicateContext,
    SlurmResourceHint,
)
from polyzymd.analyses.hydrogen_bonds import (
    HydrogenBondCompositionSettings,
    HydrogenBondsAnalysis,
    HydrogenBondSettings,
    HydrogenBondSummarySettings,
    _settings_hash,
)
from polyzymd.analyses.hydrogen_bonds._mda import (
    HBOND_EVENT_COLUMNS,
    HydrogenBondMDAAnalysis,
    artifact_to_hydrogen_bond_result,
)
from polyzymd.analyses.hydrogen_bonds._results import (
    AggregatedCompositionEntry,
    CompositionEntry,
    DirectedResiduePairResult,
    HydrogenBondAggregatedResult,
    HydrogenBondAggregatedSummary,
    HydrogenBondReplicateSummary,
    HydrogenBondResult,
    ResidueRef,
    UndirectedPairAggregate,
    UndirectedResiduePairResult,
)
from polyzymd.analyses.mda import ArtifactStore, ConditionArtifact, ReplicateArtifact
from polyzymd.analyses.stats import (
    default_scalar_comparison,
    format_pct,
    interpret_direction,
    rank_conditions,
)
from polyzymd.config.comparison import PlotSettings, PlotTheme


class _MockAtom:
    def __init__(self, index: int, chain_id: str, resid: int, resname: str, resindex: int) -> None:
        self.index = index
        self.segid = chain_id
        self.chainID = chain_id
        self.resid = resid
        self.resname = resname
        self.resindex = resindex


class _MockAtomCollection:
    def __init__(self, atoms_by_index: dict[int, _MockAtom]) -> None:
        self._atoms_by_index = atoms_by_index

    def __getitem__(self, item: int) -> _MockAtom:
        return self._atoms_by_index[item]


class _MockAtomGroup:
    def __init__(self, indices: list[int]) -> None:
        self.indices = np.array(indices, dtype=int)

    def __len__(self) -> int:
        return int(self.indices.size)

    def __or__(self, other: "_MockAtomGroup") -> "_MockAtomGroup":
        merged = sorted(set(self.indices.tolist()) | set(other.indices.tolist()))
        return _MockAtomGroup(merged)

    def __and__(self, other: "_MockAtomGroup") -> "_MockAtomGroup":
        overlap = sorted(set(self.indices.tolist()) & set(other.indices.tolist()))
        return _MockAtomGroup(overlap)


class _MockUniverseProvider:
    """Minimal MDA universe provider for plugin unit tests."""

    def __init__(self, loader: Any) -> None:
        self._loader = loader

    @classmethod
    def from_config(cls, config: Any, *, loader: Any) -> "_MockUniverseProvider":
        """Return a provider using the injected test loader."""

        del config
        return cls(loader)

    def load_universe(self, replicate: int) -> Any:
        """Load the universe from the injected test loader."""

        return self._loader.load_universe(replicate)

    def provenance_for(self, replicate: int) -> None:
        """Return no filesystem provenance for synthetic unit tests."""

        del replicate
        return None


@pytest.fixture(autouse=True)
def _use_mock_universe_provider(monkeypatch: pytest.MonkeyPatch) -> None:
    """Avoid filesystem provenance discovery in synthetic plugin tests."""

    monkeypatch.setattr(
        HydrogenBondsAnalysis,
        "_mda_universe_provider_factory",
        lambda _self: _MockUniverseProvider,
    )


def _make_mdanalysis_module(mock_hbond_cls: type) -> dict[str, types.ModuleType]:
    mda_module = types.ModuleType("MDAnalysis")
    mda_module.__version__ = "mock"

    analysis_module = types.ModuleType("MDAnalysis.analysis")
    hbonds_module = types.ModuleType("MDAnalysis.analysis.hydrogenbonds")
    hbond_analysis_module = types.ModuleType("MDAnalysis.analysis.hydrogenbonds.hbond_analysis")
    hbond_analysis_module.HydrogenBondAnalysis = mock_hbond_cls

    return {
        "MDAnalysis": mda_module,
        "MDAnalysis.analysis": analysis_module,
        "MDAnalysis.analysis.hydrogenbonds": hbonds_module,
        "MDAnalysis.analysis.hydrogenbonds.hbond_analysis": hbond_analysis_module,
    }


def _as_hydrogen_bond_result(value: Any) -> HydrogenBondResult:
    """Adapt MDA replicate artifacts to the legacy assertion model."""

    if isinstance(value, ReplicateArtifact):
        return artifact_to_hydrogen_bond_result(value)
    if isinstance(value, HydrogenBondResult):
        return value
    raise TypeError(f"Expected hydrogen-bond result or artifact, got {type(value).__name__}")


def _as_hydrogen_bond_aggregate(value: Any) -> HydrogenBondAggregatedResult:
    """Adapt MDA condition artifacts to the legacy assertion model."""

    if isinstance(value, ConditionArtifact):
        from polyzymd.analyses.hydrogen_bonds._mda import condition_artifact_to_legacy_result

        return condition_artifact_to_legacy_result(value)
    if isinstance(value, HydrogenBondAggregatedResult):
        return value
    raise TypeError(f"Expected hydrogen-bond aggregate or artifact, got {type(value).__name__}")


def test_discovered() -> None:
    """Hydrogen-bonds plugin should be auto-discovered."""
    analyses = list_analyses()
    assert "hydrogen_bonds" in analyses


def test_aliases() -> None:
    """Alias names should resolve to the plugin class."""
    assert get_analysis("hbonds") is HydrogenBondsAnalysis
    assert get_analysis("hbond") is HydrogenBondsAnalysis


def test_class_variables() -> None:
    """Plugin class variables should follow analysis contract."""
    cls = get_analysis("hydrogen_bonds")
    assert cls.name == "hydrogen_bonds"
    assert cls.Settings is HydrogenBondSettings
    assert cls.min_replicates == 1
    assert cls.has_compute_stage is True
    assert cls.has_aggregate_stage is True
    assert cls.slurm_resource_hint == SlurmResourceHint(mem="16G")


def test_settings_defaults() -> None:
    """Default settings should include protein/polymer and one summary."""
    settings = HydrogenBondSettings()
    assert settings.groups == {"protein": "chainid A", "polymer": "chainid C"}
    assert len(settings.summaries) == 1
    assert settings.summaries[0].name == "protein_polymer"
    assert settings.summaries[0].between == ("protein", "polymer")
    assert settings.allow_empty_groups is True


def test_settings_custom() -> None:
    """Custom groups and summaries should parse correctly."""
    settings = HydrogenBondSettings(
        groups={"enzyme": "chainid A", "ligand": "chainid B", "polymer": "chainid C"},
        summaries=[
            HydrogenBondSummarySettings(name="enzyme_polymer", between=("enzyme", "polymer")),
            HydrogenBondSummarySettings(name="ligand_internal", within="ligand"),
        ],
    )
    assert set(settings.groups) == {"enzyme", "ligand", "polymer"}
    assert [summary.name for summary in settings.summaries] == ["enzyme_polymer", "ligand_internal"]


def test_settings_accepts_mapping_summaries_form() -> None:
    """Settings should accept summary mappings keyed by summary name."""
    settings = HydrogenBondSettings(
        groups={"protein_all": "protein", "polymer_all": "chainid C"},
        summaries={
            "protein_polymer": {"between": ("protein_all", "polymer_all")},
        },
    )

    assert len(settings.summaries) == 1
    assert settings.summaries[0].name == "protein_polymer"
    assert settings.summaries[0].between == ("protein_all", "polymer_all")


def test_settings_mapping_summary_name_conflict_rejected() -> None:
    """Mapping key and explicit summary name must match when both are set."""
    with pytest.raises(ValidationError, match="Summary mapping key must match"):
        HydrogenBondSettings(
            groups={"protein_all": "protein", "polymer_all": "chainid C"},
            summaries={
                "protein_polymer": {
                    "name": "other_name",
                    "between": ("protein_all", "polymer_all"),
                },
            },
        )


def test_between_summary() -> None:
    """Between-mode summary should validate."""
    summary = HydrogenBondSummarySettings(name="protein_polymer", between=("protein", "polymer"))
    assert summary.between == ("protein", "polymer")
    assert summary.within is None


def test_within_summary() -> None:
    """Within-mode summary should validate."""
    summary = HydrogenBondSummarySettings(name="protein_internal", within="protein")
    assert summary.within == "protein"
    assert summary.between is None


def test_summary_both_modes_rejected() -> None:
    """Summary with both between and within should fail validation."""
    with pytest.raises(ValidationError, match="Exactly one of 'between' or 'within' must be set"):
        HydrogenBondSummarySettings(
            name="invalid",
            between=("protein", "polymer"),
            within="protein",
        )


def test_summary_neither_mode_rejected() -> None:
    """Summary with neither between nor within should fail validation."""
    with pytest.raises(ValidationError, match="Exactly one of 'between' or 'within' must be set"):
        HydrogenBondSummarySettings(name="invalid")


def test_summary_invalid_group_ref() -> None:
    """Summary group names must exist in the group mapping."""
    with pytest.raises(ValidationError, match="references unknown group"):
        HydrogenBondSettings(
            groups={"protein": "chainid A"},
            summaries=[HydrogenBondSummarySettings(name="bad", between=("protein", "polymer"))],
        )


def test_summary_duplicate_names_rejected() -> None:
    """Summary names should be unique."""
    with pytest.raises(ValidationError, match="Duplicate summary name"):
        HydrogenBondSettings(
            summaries=[
                HydrogenBondSummarySettings(name="dup", between=("protein", "polymer")),
                HydrogenBondSummarySettings(name="dup", within="protein"),
            ]
        )


def test_composition_settings() -> None:
    """Composition partitions should parse correctly."""
    settings = HydrogenBondSettings(
        composition=HydrogenBondCompositionSettings(
            partitions={"protein": "chainid A", "polymer": "chainid C"}
        )
    )
    assert settings.composition is not None
    assert settings.composition.partitions["protein"] == "chainid A"


def test_distance_cutoff_validation() -> None:
    """Distance cutoff must be strictly positive."""
    with pytest.raises(ValidationError):
        HydrogenBondSettings(distance_cutoff=0)


def test_settings_hash_changes_with_cutoff() -> None:
    """Settings hash should change when distance cutoff changes."""
    settings_a = HydrogenBondSettings(distance_cutoff=3.0)
    settings_b = HydrogenBondSettings(distance_cutoff=3.5)
    assert _settings_hash(settings_a) != _settings_hash(settings_b)


def test_settings_hash_changes_with_groups() -> None:
    """Settings hash should change when group selections change."""
    settings_a = HydrogenBondSettings(groups={"protein": "chainid A", "polymer": "chainid C"})
    settings_b = HydrogenBondSettings(groups={"protein": "chainid A", "polymer": "chainid B"})
    assert _settings_hash(settings_a) != _settings_hash(settings_b)


def test_settings_hash_stable_and_cache_name_stable() -> None:
    """Settings hash and cache key should be stable for identical settings."""
    settings_a = HydrogenBondSettings()
    settings_b = HydrogenBondSettings()

    hash_a = _settings_hash(settings_a)
    hash_b = _settings_hash(settings_b)
    assert hash_a == hash_b
    assert f"hbonds_eq0ns_{hash_a}.json" == f"hbonds_eq0ns_{hash_b}.json"


def test_get_trajectory_window_uses_timestep_override(tmp_path: Path) -> None:
    """Runner window resolution should honor settings.timestep_ps overrides."""
    analysis = HydrogenBondsAnalysis()
    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    settings = HydrogenBondSettings(timestep_ps=50.0)
    ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="1ns",
        recompute=True,
        settings=settings,
    )
    loader = MagicMock()
    loader.get_timestep.return_value = 10.0
    universe = MagicMock()
    universe.trajectory = [object()] * 100

    window = analysis.get_trajectory_window(ctx, 1, loader, universe)

    assert window.start == 20
    assert window.timestep_ps == pytest.approx(50.0)


def test_mda_analysis_records_effective_timestep_metadata() -> None:
    """Hydrogen-bond MDA job should record raw timestep, stride, and effective spacing."""

    class MockHydrogenBondAnalysis:
        def __init__(self, **kwargs) -> None:
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            del start, stop, step, verbose
            self.results.hbonds = np.empty((0, 6), dtype=float)

    universe = MagicMock()
    universe.trajectory = [object()] * 9
    universe.atoms = _MockAtomCollection(
        {
            0: _MockAtom(0, "A", 10, "SER", 0),
            1: _MockAtom(1, "C", 100, "OEG", 1),
        }
    )
    selections = {
        "chainid A": _MockAtomGroup([0]),
        "chainid C": _MockAtomGroup([1]),
    }
    universe.select_atoms.side_effect = lambda selection, updating: selections[selection]
    analysis = HydrogenBondMDAAnalysis(
        universe=universe,
        settings=HydrogenBondSettings(),
        condition_label="test",
        replicate=1,
        raw_timestep_ps=10.0,
    )

    with patch.dict(sys.modules, _make_mdanalysis_module(MockHydrogenBondAnalysis)):
        analysis.run(start=0, stop=9, step=3)

    assert analysis.plan is not None
    assert analysis.plan.timestep_ps == pytest.approx(30.0)
    assert analysis.plan.raw_timestep_ps == pytest.approx(10.0)
    assert analysis.plan.frame_stride == 3


def test_replicate_artifact_preserves_timing_metadata(tmp_path: Path) -> None:
    """Hydrogen-bond artifacts should preserve MDA timing metadata."""
    analysis = HydrogenBondsAnalysis()
    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    settings = HydrogenBondSettings()
    ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="10ns",
        recompute=True,
        settings=settings,
    )
    summary = HydrogenBondReplicateSummary(
        name="protein_polymer",
        mode="between",
        group_names=["protein", "polymer"],
        n_frames_used=3,
        mean_hbonds_per_frame=0.0,
        fraction_frames_with_any_hbond=0.0,
        counts_per_frame=[0, 0, 0],
    )
    artifact = ReplicateArtifact(
        analysis_name="hydrogen_bonds",
        condition_label="test",
        replicate=1,
        payload={"summaries": [summary.model_dump(mode="json")], "composition_entries": []},
        metadata={
            "settings_fingerprint": _settings_hash(settings),
            "config_hash": "hash123",
            "equilibration_time": 10.0,
            "equilibration_unit": "ns",
            "selection_string": "(chainid A) or (chainid C)",
            "timestep_ps": 30.0,
            "raw_timestep_ps": 10.0,
            "frame_stride": 3,
        },
    )

    del analysis, ctx
    result = artifact_to_hydrogen_bond_result(
        artifact, settings_fingerprint=_settings_hash(settings)
    )

    assert result.timestep_ps == pytest.approx(30.0)
    assert result.raw_timestep_ps == pytest.approx(10.0)
    assert result.frame_stride == 3


def test_angle_cutoff_validation() -> None:
    """Angle cutoff must be in the range (0, 180]."""
    with pytest.raises(ValidationError):
        HydrogenBondSettings(angle_cutoff=0)
    with pytest.raises(ValidationError):
        HydrogenBondSettings(angle_cutoff=181)


def test_replicate_result_save_load(tmp_path: Path) -> None:
    """Replicate result models should round-trip with save/load."""
    settings = HydrogenBondSettings()
    result = HydrogenBondResult(
        replicate=1,
        settings_fingerprint=_settings_hash(settings),
        summaries=[
            HydrogenBondReplicateSummary(
                name="protein_polymer",
                mode="between",
                group_names=["protein", "polymer"],
                n_frames_used=100,
                mean_hbonds_per_frame=3.2,
                fraction_frames_with_any_hbond=0.82,
                counts_per_frame=[3, 2, 4],
            )
        ],
    )

    path = result.save(tmp_path / "hbonds_replicate.json")
    loaded = HydrogenBondResult.load(path)

    assert loaded.replicate == 1
    assert loaded.settings_fingerprint == _settings_hash(settings)
    assert loaded.summaries[0].name == "protein_polymer"
    assert loaded.summaries[0].mean_hbonds_per_frame == pytest.approx(3.2)


def test_aggregated_result_save_load(tmp_path: Path) -> None:
    """Aggregated result models should round-trip with save/load."""
    settings = HydrogenBondSettings()
    result = HydrogenBondAggregatedResult(
        settings_fingerprint=_settings_hash(settings),
        replicates=[1, 2],
        n_replicates=2,
        summaries=[
            HydrogenBondAggregatedSummary(
                name="protein_polymer",
                mode="between",
                group_names=["protein", "polymer"],
                n_replicates=2,
                mean_hbonds_per_frame=3.0,
                sem_hbonds_per_frame=0.2,
                per_replicate_mean_hbonds=[2.8, 3.2],
                mean_fraction_with_any=0.8,
                sem_fraction_with_any=0.05,
                per_replicate_fraction_with_any=[0.75, 0.85],
            )
        ],
    )

    path = result.save(tmp_path / "hbonds_aggregated.json")
    loaded = HydrogenBondAggregatedResult.load(path)

    assert loaded.n_replicates == 2
    assert loaded.replicates == [1, 2]
    assert loaded.settings_fingerprint == _settings_hash(settings)
    assert loaded.summaries[0].name == "protein_polymer"


def test_load_aggregated_result_rejects_stale_settings_fingerprint(tmp_path: Path) -> None:
    """Aggregated cache loads should reject results from different settings."""
    analysis = HydrogenBondsAnalysis()
    current_settings = HydrogenBondSettings()
    stale_settings = HydrogenBondSettings(distance_cutoff=3.5)

    aggregated_dir = tmp_path / "aggregated"
    aggregated_dir.mkdir()
    HydrogenBondAggregatedResult(
        settings_fingerprint=_settings_hash(stale_settings),
        replicates=[1, 2],
        n_replicates=2,
        summaries=[_make_aggregated_summary("protein_polymer", 3.0, 0.2, [2.8, 3.2])],
    ).save(aggregated_dir / "result.json")

    with pytest.raises(ValueError, match="current settings require"):
        analysis._load_aggregated_result(
            aggregated_dir,
            settings=current_settings,
            condition_label="CondA",
        )


def test_load_aggregated_result_ignores_noncanonical_json(tmp_path: Path) -> None:
    """Aggregated loads should not fall back to arbitrary JSON files."""
    analysis = HydrogenBondsAnalysis()
    aggregated_dir = tmp_path / "aggregated"
    aggregated_dir.mkdir()
    HydrogenBondAggregatedResult(
        settings_fingerprint=_settings_hash(HydrogenBondSettings()),
        replicates=[1, 2],
        n_replicates=2,
        summaries=[_make_aggregated_summary("protein_polymer", 3.0, 0.2, [2.8, 3.2])],
    ).save(aggregated_dir / "legacy_summary.json")

    assert analysis._load_aggregated_result(aggregated_dir) is None


def test_aggregated_summary_deserializes_legacy_std_unique_pairs_field() -> None:
    """Legacy std field should populate sem_unique_pairs_per_frame."""
    legacy_payload = {
        "name": "protein_polymer",
        "mode": "between",
        "group_names": ["protein", "polymer"],
        "n_replicates": 2,
        "mean_hbonds_per_frame": 3.0,
        "sem_hbonds_per_frame": 0.2,
        "per_replicate_mean_hbonds": [2.8, 3.2],
        "mean_unique_pairs_per_frame": 1.5,
        "std_unique_pairs_per_frame": 0.4,
        "mean_fraction_with_any": 0.8,
        "sem_fraction_with_any": 0.05,
        "per_replicate_fraction_with_any": [0.75, 0.85],
    }
    summary = HydrogenBondAggregatedSummary.model_validate(legacy_payload)
    assert summary.sem_unique_pairs_per_frame == pytest.approx(0.4)


def test_aggregated_summary_deserializes_new_sem_unique_pairs_field() -> None:
    """New sem field name should deserialize directly."""
    payload = {
        "name": "protein_polymer",
        "mode": "between",
        "group_names": ["protein", "polymer"],
        "n_replicates": 2,
        "mean_hbonds_per_frame": 3.0,
        "sem_hbonds_per_frame": 0.2,
        "per_replicate_mean_hbonds": [2.8, 3.2],
        "mean_unique_pairs_per_frame": 1.5,
        "sem_unique_pairs_per_frame": 0.33,
        "mean_fraction_with_any": 0.8,
        "sem_fraction_with_any": 0.05,
        "per_replicate_fraction_with_any": [0.75, 0.85],
    }
    summary = HydrogenBondAggregatedSummary.model_validate(payload)
    assert summary.sem_unique_pairs_per_frame == pytest.approx(0.33)


def test_result_summary() -> None:
    """Result summary should return a non-empty string."""
    result = HydrogenBondResult(
        summaries=[
            HydrogenBondReplicateSummary(
                name="protein_polymer",
                mode="between",
                group_names=["protein", "polymer"],
                n_frames_used=10,
                mean_hbonds_per_frame=1.2,
                fraction_frames_with_any_hbond=0.5,
                counts_per_frame=[1, 2],
            )
        ]
    )
    text = result.summary()
    assert isinstance(text, str)
    assert text


def test_residue_ref_label() -> None:
    """ResidueRef.label should follow resname-resid-chain format."""
    residue = ResidueRef(chain_id="A", resid=123, resname="SER")
    assert residue.label == "SER123:A"


def test_run_replicate_basic(tmp_path: Path) -> None:
    """run_replicate should return populated summary metrics."""

    instances: list[MockHydrogenBondAnalysis] = []

    class MockHydrogenBondAnalysis:
        def __init__(self, **kwargs) -> None:
            self.kwargs = kwargs
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))
            self.run_args: dict[str, int | bool | None] | None = None
            instances.append(self)

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            self.run_args = {"start": start, "stop": stop, "step": step, "verbose": verbose}
            self.results.hbonds = np.array(
                [
                    [0, 0, 10, 1, 2.8, 160.0],
                    [1, 0, 10, 1, 2.9, 158.0],
                ],
                dtype=float,
            )

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "C", 100, "OEG", 1),
    }
    universe = MagicMock()
    universe.trajectory = [object(), object(), object(), object(), object()]
    universe.atoms = _MockAtomCollection(atoms)
    selections = {
        "chainid A": _MockAtomGroup([0]),
        "chainid C": _MockAtomGroup([1]),
    }
    universe.select_atoms.side_effect = lambda selection, updating: selections[selection]

    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    settings = HydrogenBondSettings()
    ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=settings,
    )

    analysis = HydrogenBondsAnalysis()
    mock_modules = _make_mdanalysis_module(MockHydrogenBondAnalysis)

    with (
        patch.dict(sys.modules, mock_modules),
        patch("polyzymd.analyses.hydrogen_bonds.TrajectoryLoader") as mock_loader_cls,
        patch("polyzymd.analyses.hydrogen_bonds._mda.compute_config_hash", return_value="abc123"),
    ):
        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader
        mock_loader.load_universe.return_value = universe
        mock_loader.get_timestep.return_value = 10.0

        artifact = analysis.run_replicate(ctx, 1)
        result = _as_hydrogen_bond_result(artifact)

    assert isinstance(artifact, ReplicateArtifact)
    assert artifact.payload["n_events"] == 2
    assert "hbonds" not in artifact.payload
    assert artifact.sidecars[0].path == "sidecars/hydrogen_bonds_events.npz"
    with np.load(tmp_path / "run_1" / artifact.sidecars[0].path) as sidecar:
        assert sidecar["hbonds"].shape == (2, 6)
    assert isinstance(result, HydrogenBondResult)
    assert result.selection_string == "(chainid A) or (chainid C)"
    assert len(result.summaries) == 1
    summary = result.summaries[0]
    assert summary.name == "protein_polymer"
    assert summary.mean_hbonds_per_frame == pytest.approx(0.4)
    assert summary.fraction_frames_with_any_hbond == pytest.approx(0.4)
    assert summary.counts_per_frame == [1, 1, 0, 0, 0]
    assert len(summary.directed_residue_pairs) == 1
    assert len(summary.undirected_residue_pairs) == 1
    assert len(instances) == 1
    assert instances[0].kwargs["d_a_cutoff"] == pytest.approx(settings.distance_cutoff)
    assert instances[0].kwargs["d_h_a_angle_cutoff"] == pytest.approx(settings.angle_cutoff)
    assert instances[0].run_args is not None
    assert instances[0].run_args["start"] == 0


def test_run_replicate_warns_once_for_default_groups_and_summaries(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Default groups/summaries should emit one warning per analysis run."""

    class MockHydrogenBondAnalysis:
        def __init__(self, **kwargs) -> None:
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            self.results.hbonds = np.empty((0, 6), dtype=float)

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "C", 100, "OEG", 1),
    }
    universe = MagicMock()
    universe.trajectory = [object(), object(), object()]
    universe.atoms = _MockAtomCollection(atoms)
    universe.select_atoms.side_effect = lambda selection, updating: {
        "chainid A": _MockAtomGroup([0]),
        "chainid C": _MockAtomGroup([1]),
    }[selection]

    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1, 2),
        sim_config=MagicMock(),
    )
    settings = HydrogenBondSettings()
    ctx_1 = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=settings,
    )
    ctx_2 = ReplicateContext(
        condition=condition,
        replicate=2,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_2",
        equilibration="0ns",
        recompute=True,
        settings=settings,
    )

    analysis = HydrogenBondsAnalysis()
    HydrogenBondsAnalysis._defaults_warned = False
    mock_modules = _make_mdanalysis_module(MockHydrogenBondAnalysis)

    with (
        patch.dict(sys.modules, mock_modules),
        patch("polyzymd.analyses.hydrogen_bonds.TrajectoryLoader") as mock_loader_cls,
        patch("polyzymd.analyses.hydrogen_bonds._mda.compute_config_hash", return_value="abc123"),
        caplog.at_level("WARNING", logger="polyzymd.analyses.hydrogen_bonds"),
    ):
        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader
        mock_loader.load_universe.return_value = universe
        mock_loader.get_timestep.return_value = 10.0

        analysis.run_replicate(ctx_1, 1)
        analysis.run_replicate(ctx_2, 2)

    assert caplog.text.count("No explicit groups/summaries in YAML config — using defaults") == 1


def test_run_replicate_empty_selection(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """run_replicate should handle empty groups with zeroed summaries."""

    class MockHydrogenBondAnalysis:
        def __init__(self, **kwargs) -> None:
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            return None

    universe = MagicMock()
    universe.trajectory = [object(), object(), object()]
    universe.atoms = _MockAtomCollection({})
    universe.select_atoms.return_value = _MockAtomGroup([])

    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    settings = HydrogenBondSettings()
    ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=settings,
    )

    analysis = HydrogenBondsAnalysis()
    mock_modules = _make_mdanalysis_module(MockHydrogenBondAnalysis)

    with (
        patch.dict(sys.modules, mock_modules),
        patch("polyzymd.analyses.hydrogen_bonds.TrajectoryLoader") as mock_loader_cls,
        patch("polyzymd.analyses.hydrogen_bonds._mda.compute_config_hash", return_value="abc123"),
        caplog.at_level("WARNING"),
    ):
        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader
        mock_loader.load_universe.return_value = universe
        mock_loader.get_timestep.return_value = 10.0

        result = _as_hydrogen_bond_result(analysis.run_replicate(ctx, 1))

    assert isinstance(result, HydrogenBondResult)
    assert len(result.summaries) == 1
    summary = result.summaries[0]
    assert summary.mean_hbonds_per_frame == 0.0
    assert summary.counts_per_frame == [0, 0, 0]
    assert "matched 0 atoms" in caplog.text
    assert "No atoms selected for any summary" in caplog.text


def test_run_replicate_skips_only_empty_summary_and_keeps_other_summaries(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Empty-group summaries should be skipped without failing other summaries."""

    class MockHydrogenBondAnalysis:
        def __init__(self, **kwargs) -> None:
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            self.results.hbonds = np.empty((0, 6), dtype=float)

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "A", 11, "GLY", 1),
    }
    universe = MagicMock()
    universe.trajectory = [object(), object(), object()]
    universe.atoms = _MockAtomCollection(atoms)
    selections = {
        "chainid A": _MockAtomGroup([0, 1]),
        "chainid C": _MockAtomGroup([]),
    }
    universe.select_atoms.side_effect = lambda selection, updating: selections[selection]

    condition = Condition(
        label="No Polymer (Control)",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    settings = HydrogenBondSettings(
        groups={"protein": "chainid A", "polymer": "chainid C"},
        summaries=[
            HydrogenBondSummarySettings(name="protein_polymer", between=("protein", "polymer")),
            HydrogenBondSummarySettings(name="protein_internal", within="protein"),
        ],
    )
    ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=settings,
    )

    analysis = HydrogenBondsAnalysis()
    mock_modules = _make_mdanalysis_module(MockHydrogenBondAnalysis)

    with (
        patch.dict(sys.modules, mock_modules),
        patch("polyzymd.analyses.hydrogen_bonds.TrajectoryLoader") as mock_loader_cls,
        patch("polyzymd.analyses.hydrogen_bonds._mda.compute_config_hash", return_value="abc123"),
        caplog.at_level("WARNING", logger="polyzymd.analyses.hydrogen_bonds"),
    ):
        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader
        mock_loader.load_universe.return_value = universe
        mock_loader.get_timestep.return_value = 10.0

        result = _as_hydrogen_bond_result(analysis.run_replicate(ctx, 1))

    assert isinstance(result, HydrogenBondResult)
    summaries = {summary.name: summary for summary in result.summaries}
    assert set(summaries) == {"protein_polymer", "protein_internal"}

    skipped_summary = summaries["protein_polymer"]
    assert skipped_summary.mean_hbonds_per_frame == 0.0
    assert skipped_summary.counts_per_frame == [0, 0, 0]

    computed_summary = summaries["protein_internal"]
    assert computed_summary.mean_hbonds_per_frame == 0.0
    assert computed_summary.counts_per_frame == [0, 0, 0]

    assert "skipping summary 'protein_polymer'" in caplog.text
    assert "condition='No Polymer (Control)' replicate=1" in caplog.text


def test_run_replicate_empty_group_raises_by_default(tmp_path: Path) -> None:
    """Empty group selections should raise when strict mode is enabled."""

    class MockHydrogenBondAnalysis:
        def __init__(self, **kwargs) -> None:
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            return None

    universe = MagicMock()
    universe.trajectory = [object(), object(), object()]
    universe.atoms = _MockAtomCollection({})
    universe.select_atoms.return_value = _MockAtomGroup([])

    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=HydrogenBondSettings(allow_empty_groups=False),
    )

    analysis = HydrogenBondsAnalysis()
    mock_modules = _make_mdanalysis_module(MockHydrogenBondAnalysis)

    with (
        patch.dict(sys.modules, mock_modules),
        patch("polyzymd.analyses.hydrogen_bonds.TrajectoryLoader") as mock_loader_cls,
        patch("polyzymd.analyses.hydrogen_bonds._mda.compute_config_hash", return_value="abc123"),
    ):
        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader
        mock_loader.load_universe.return_value = universe
        mock_loader.get_timestep.return_value = 10.0

        with pytest.raises(Exception, match="allow_empty_groups: true"):
            analysis.run_replicate(ctx, 1)


def test_load_replicate_result_ignores_legacy_custom_cache(tmp_path: Path) -> None:
    """Legacy custom hbonds_eq caches should not be active after MDA migration."""
    run_dir = tmp_path / "run_1"
    run_dir.mkdir()
    HydrogenBondResult(replicate=1, config_hash="abc123", summaries=[]).save(
        run_dir / "hbonds_eq0ns_deadbeef.json"
    )
    analysis = HydrogenBondsAnalysis()

    assert analysis._load_replicate_result(run_dir) is None


def test_equilibration_exceeds_trajectory_raises(tmp_path: Path) -> None:
    """Equilibration beyond trajectory length should raise a clear error."""

    class MockHydrogenBondAnalysis:
        def __init__(self, **kwargs) -> None:
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            return None

    universe = MagicMock()
    universe.trajectory = [object(), object(), object()]
    universe.atoms = _MockAtomCollection({})
    universe.select_atoms.side_effect = lambda selection, updating: _MockAtomGroup([0])

    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="100ns",
        recompute=True,
        settings=HydrogenBondSettings(),
    )

    analysis = HydrogenBondsAnalysis()
    mock_modules = _make_mdanalysis_module(MockHydrogenBondAnalysis)

    with (
        patch.dict(sys.modules, mock_modules),
        patch("polyzymd.analyses.hydrogen_bonds.TrajectoryLoader") as mock_loader_cls,
        patch("polyzymd.analyses.hydrogen_bonds._mda.compute_config_hash", return_value="abc123"),
    ):
        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader
        mock_loader.load_universe.return_value = universe
        mock_loader.get_timestep.return_value = 10.0

        with pytest.raises(ValueError) as exc_info:
            analysis.run_replicate(ctx, 1)

    message = str(exc_info.value)
    assert "Equilibration time" in message
    assert "trajectory length" in message
    assert "Reduce --eq-time or check your trajectory" in message


def test_equilibration_leaves_one_frame_warns(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """One remaining frame after equilibration should warn and still run."""

    class MockHydrogenBondAnalysis:
        def __init__(self, **kwargs) -> None:
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))
            self.run_args: dict[str, int | bool | None] | None = None

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            self.run_args = {"start": start, "stop": stop, "step": step, "verbose": verbose}
            self.results.hbonds = np.array([[4, 0, 10, 1, 2.8, 160.0]], dtype=float)

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "C", 100, "OEG", 1),
    }
    universe = MagicMock()
    universe.trajectory = [object(), object(), object(), object(), object()]
    universe.atoms = _MockAtomCollection(atoms)
    selections = {
        "chainid A": _MockAtomGroup([0]),
        "chainid C": _MockAtomGroup([1]),
    }
    universe.select_atoms.side_effect = lambda selection, updating: selections[selection]

    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    settings = HydrogenBondSettings()
    ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="40ps",
        recompute=True,
        settings=settings,
    )

    analysis = HydrogenBondsAnalysis()
    mock_modules = _make_mdanalysis_module(MockHydrogenBondAnalysis)

    with (
        patch.dict(sys.modules, mock_modules),
        patch("polyzymd.analyses.hydrogen_bonds.TrajectoryLoader") as mock_loader_cls,
        patch("polyzymd.analyses.hydrogen_bonds._mda.compute_config_hash", return_value="abc123"),
        caplog.at_level("WARNING"),
    ):
        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader
        mock_loader.load_universe.return_value = universe
        mock_loader.get_timestep.return_value = 10.0

        result = _as_hydrogen_bond_result(analysis.run_replicate(ctx, 1))

    assert isinstance(result, HydrogenBondResult)
    assert result.summaries[0].counts_per_frame == [1]

    warning_messages = [record.getMessage() for record in caplog.records]
    assert any(
        "Warning: Skipping" in message and "trajectory for equilibration" in message
        for message in warning_messages
    )
    assert any(
        "Only 1 frame(s) remain after equilibration window [4:5:1]" in message
        for message in warning_messages
    )


def test_intra_residue_exclusion(tmp_path: Path) -> None:
    """H-bond events within the same residue should be excluded."""

    class MockHydrogenBondAnalysis:
        def __init__(self, **kwargs) -> None:
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            self.results.hbonds = np.array(
                [
                    [0, 0, 10, 1, 2.8, 160.0],
                    [1, 0, 10, 1, 2.8, 160.0],
                ],
                dtype=float,
            )

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "A", 10, "SER", 0),
    }
    universe = MagicMock()
    universe.trajectory = [object(), object(), object()]
    universe.atoms = _MockAtomCollection(atoms)
    universe.select_atoms.side_effect = lambda selection, updating: _MockAtomGroup([0, 1])

    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    settings = HydrogenBondSettings(
        groups={"protein": "chainid A"},
        summaries=[HydrogenBondSummarySettings(name="protein_internal", within="protein")],
    )
    ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=settings,
    )

    analysis = HydrogenBondsAnalysis()
    mock_modules = _make_mdanalysis_module(MockHydrogenBondAnalysis)

    with (
        patch.dict(sys.modules, mock_modules),
        patch("polyzymd.analyses.hydrogen_bonds.TrajectoryLoader") as mock_loader_cls,
        patch("polyzymd.analyses.hydrogen_bonds._mda.compute_config_hash", return_value="abc123"),
    ):
        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader
        mock_loader.load_universe.return_value = universe
        mock_loader.get_timestep.return_value = 10.0

        result = _as_hydrogen_bond_result(analysis.run_replicate(ctx, 1))

    summary = result.summaries[0]
    assert summary.mean_hbonds_per_frame == 0.0
    assert summary.fraction_frames_with_any_hbond == 0.0
    assert summary.counts_per_frame == [0, 0, 0]


def test_compute_composition_basic() -> None:
    """Composition should classify partition pairs with correct fractions."""
    analysis = HydrogenBondsAnalysis()
    composition_settings = HydrogenBondCompositionSettings(
        partitions={"protein": "chainid A", "polymer": "chainid C", "substrate": "chainid B"}
    )

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "C", 100, "OEG", 1),
        2: _MockAtom(2, "B", 200, "LIG", 2),
    }
    universe = MagicMock()
    universe.atoms = _MockAtomCollection(atoms)
    universe.select_atoms.side_effect = lambda selection: {
        "chainid A": _MockAtomGroup([0]),
        "chainid C": _MockAtomGroup([1]),
        "chainid B": _MockAtomGroup([2]),
    }[selection]

    hbond_array = np.array(
        [
            [0, 0, 10, 1, 2.8, 160.0],
            [1, 0, 10, 1, 2.9, 158.0],
            [2, 2, 11, 0, 2.9, 155.0],
        ],
        dtype=float,
    )

    entries = analysis._compute_composition(
        composition_settings=composition_settings,
        hbond_array=hbond_array,
        universe=universe,
        start_frame=0,
        n_frames=5,
        allow_overlapping=True,
    )

    by_pair = {(entry.donor_partition, entry.acceptor_partition): entry for entry in entries}
    assert by_pair[("protein", "polymer")].mean_hbonds_per_frame == pytest.approx(2 / 5)
    assert by_pair[("protein", "polymer")].fraction_of_total == pytest.approx(2 / 3)
    assert by_pair[("substrate", "protein")].mean_hbonds_per_frame == pytest.approx(1 / 5)
    assert by_pair[("substrate", "protein")].fraction_of_total == pytest.approx(1 / 3)


def test_composition_disjoint_warning(caplog: pytest.LogCaptureFixture) -> None:
    """Overlapping partitions should emit a warning."""
    analysis = HydrogenBondsAnalysis()
    composition_settings = HydrogenBondCompositionSettings(
        partitions={"protein": "chainid A", "polymer": "chainid C"}
    )

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "C", 100, "OEG", 1),
    }
    universe = MagicMock()
    universe.atoms = _MockAtomCollection(atoms)
    universe.select_atoms.side_effect = lambda selection: {
        "chainid A": _MockAtomGroup([0, 1]),
        "chainid C": _MockAtomGroup([1]),
    }[selection]

    with caplog.at_level("WARNING"):
        _ = analysis._compute_composition(
            composition_settings=composition_settings,
            hbond_array=np.empty((0, 6), dtype=float),
            universe=universe,
            start_frame=0,
            n_frames=3,
            allow_overlapping=True,
        )

    assert "Overlapping atoms will be counted in BOTH partitions" in caplog.text
    assert "composition fractions may exceed 1.0" in caplog.text


def test_composition_warns_dynamic_selection(caplog: pytest.LogCaptureFixture) -> None:
    """Coordinate-dependent partition selections should emit a warning."""
    analysis = HydrogenBondsAnalysis()
    composition_settings = HydrogenBondCompositionSettings(
        partitions={"dynamic": "around 5.0 protein", "protein": "chainid A"}
    )

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "C", 100, "OEG", 1),
    }
    universe = MagicMock()
    universe.atoms = _MockAtomCollection(atoms)
    universe.select_atoms.side_effect = lambda selection: {
        "around 5.0 protein": _MockAtomGroup([0]),
        "chainid A": _MockAtomGroup([0]),
    }[selection]

    with caplog.at_level("WARNING"):
        _ = analysis._compute_composition(
            composition_settings=composition_settings,
            hbond_array=np.empty((0, 6), dtype=float),
            universe=universe,
            start_frame=0,
            n_frames=3,
            allow_overlapping=True,
        )

    assert "uses coordinate-dependent selection" in caplog.text


def test_settings_reject_dynamic_group_selection_when_updating() -> None:
    """Dynamic group selections should be rejected with update_selections=True."""
    with pytest.raises(ValidationError, match="coordinate-dependent"):
        HydrogenBondSettings(
            groups={"protein": "chainid A", "polymer": "around 5.0 chainid A"},
            update_selections=True,
        )


def test_settings_allow_dynamic_group_selection_when_not_updating() -> None:
    """Dynamic group selections should validate with update_selections=False."""
    settings = HydrogenBondSettings(
        groups={"protein": "chainid A", "polymer": "around 5.0 chainid A"},
        update_selections=False,
    )
    assert settings.update_selections is False


def test_dynamic_selection_policy_recorded_in_artifact(tmp_path: Path) -> None:
    """Dynamic selection/update policy should be explicit in MDA artifacts."""

    class MockHydrogenBondAnalysis:
        def __init__(self, **kwargs) -> None:
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            del start, stop, step, verbose
            self.results.hbonds = np.empty((0, 6), dtype=float)

    atoms = {0: _MockAtom(0, "A", 10, "SER", 0)}
    universe = MagicMock()
    universe.trajectory = [object(), object()]
    universe.atoms = _MockAtomCollection(atoms)
    universe.select_atoms.side_effect = lambda selection, updating: {
        "around 5.0 chainid A": _MockAtomGroup([0]),
    }[selection]
    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    settings = HydrogenBondSettings(
        groups={"dynamic": "around 5.0 chainid A"},
        summaries=[HydrogenBondSummarySettings(name="dynamic_internal", within="dynamic")],
        update_selections=False,
    )
    ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=settings,
    )
    analysis = HydrogenBondsAnalysis()

    with (
        patch.dict(sys.modules, _make_mdanalysis_module(MockHydrogenBondAnalysis)),
        patch("polyzymd.analyses.hydrogen_bonds.TrajectoryLoader") as mock_loader_cls,
        patch("polyzymd.analyses.hydrogen_bonds._mda.compute_config_hash", return_value="abc123"),
    ):
        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader
        mock_loader.load_universe.return_value = universe
        mock_loader.get_timestep.return_value = 10.0

        artifact = analysis.run_replicate(ctx, 1)

    assert isinstance(artifact, ReplicateArtifact)
    policy = artifact.provenance["dynamic_selection_policy"]
    assert policy["update_selections"] is False
    assert policy["dynamic_group_selections"] == {"dynamic": "around 5.0 chainid A"}
    assert any("membership is evaluated once" in warning for warning in artifact.warnings)


def test_dynamic_composition_warning_recorded_in_artifact(tmp_path: Path) -> None:
    """Coordinate-dependent composition warnings should be artifact provenance."""

    class MockHydrogenBondAnalysis:
        def __init__(self, **kwargs) -> None:
            del kwargs
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            del start, stop, step, verbose
            self.results.hbonds = np.array([[0, 0, 10, 1, 2.8, 160.0]], dtype=float)

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "C", 100, "OEG", 1),
    }
    universe = MagicMock()
    universe.trajectory = [object(), object()]
    universe.atoms = _MockAtomCollection(atoms)

    def select_atoms(selection: str, updating: bool | None = None) -> _MockAtomGroup:
        del updating
        return {
            "chainid A": _MockAtomGroup([0]),
            "chainid C": _MockAtomGroup([1]),
            "around 5.0 chainid A": _MockAtomGroup([0]),
        }[selection]

    universe.select_atoms.side_effect = select_atoms
    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    settings = HydrogenBondSettings(
        composition=HydrogenBondCompositionSettings(
            partitions={"dynamic": "around 5.0 chainid A", "polymer": "chainid C"}
        )
    )
    ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=settings,
    )
    analysis = HydrogenBondsAnalysis()

    with (
        patch.dict(sys.modules, _make_mdanalysis_module(MockHydrogenBondAnalysis)),
        patch("polyzymd.analyses.hydrogen_bonds.TrajectoryLoader") as mock_loader_cls,
        patch("polyzymd.analyses.hydrogen_bonds._mda.compute_config_hash", return_value="abc123"),
    ):
        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader
        mock_loader.load_universe.return_value = universe
        mock_loader.get_timestep.return_value = 10.0

        artifact = analysis.run_replicate(ctx, 1)

    assert isinstance(artifact, ReplicateArtifact)
    expected = "Composition partition 'dynamic' uses coordinate-dependent selection"
    assert any(expected in warning for warning in artifact.warnings)
    assert any(expected in warning for warning in artifact.provenance["composition_warnings"])


def test_composition_skips_unpartitioned_atoms() -> None:
    """Events with donor or acceptor outside partitions should be skipped."""
    analysis = HydrogenBondsAnalysis()
    composition_settings = HydrogenBondCompositionSettings(
        partitions={"protein": "chainid A", "polymer": "chainid C"}
    )

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "C", 100, "OEG", 1),
        2: _MockAtom(2, "D", 200, "HOH", 2),
    }
    universe = MagicMock()
    universe.atoms = _MockAtomCollection(atoms)
    universe.select_atoms.side_effect = lambda selection: {
        "chainid A": _MockAtomGroup([0]),
        "chainid C": _MockAtomGroup([1]),
    }[selection]

    hbond_array = np.array(
        [
            [0, 0, 10, 1, 2.8, 160.0],
            [1, 0, 10, 2, 2.9, 158.0],
            [2, 2, 10, 1, 2.9, 158.0],
        ],
        dtype=float,
    )

    entries = analysis._compute_composition(
        composition_settings=composition_settings,
        hbond_array=hbond_array,
        universe=universe,
        start_frame=0,
        n_frames=4,
        allow_overlapping=True,
    )

    assert len(entries) == 1
    entry = entries[0]
    assert (entry.donor_partition, entry.acceptor_partition) == ("protein", "polymer")
    assert entry.mean_hbonds_per_frame == pytest.approx(0.25)
    assert entry.fraction_of_total == pytest.approx(1.0)


def test_load_replicate_timeseries_corrupt_json(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Corrupt replicate cache JSON should be skipped without crashing."""
    analysis = HydrogenBondsAnalysis()
    analysis_dir = tmp_path / "analysis" / "conda" / "hydrogen_bonds"
    run_dir = analysis_dir / "run_1"
    run_dir.mkdir(parents=True)
    (run_dir / "result.json").write_text("{not valid json", encoding="utf-8")

    data = {"CondA": {"analysis_dir": analysis_dir, "replicates": [1]}}

    with caplog.at_level("DEBUG"):
        loaded = analysis._load_replicate_timeseries(data, ["CondA"])

    assert loaded == {}
    assert "Could not load replicate result" in caplog.text


def test_load_replicate_timeseries_rejects_stale_settings_fingerprint(tmp_path: Path) -> None:
    """Timeseries loading should reject stale replicate caches during plotting."""
    analysis = HydrogenBondsAnalysis()
    current_settings = HydrogenBondSettings()
    stale_settings = HydrogenBondSettings(distance_cutoff=3.5)
    analysis_dir = tmp_path / "analysis" / "conda" / "hydrogen_bonds"
    run_dir = analysis_dir / "run_1"
    run_dir.mkdir(parents=True)

    HydrogenBondResult(
        config_hash="cfg123",
        settings_fingerprint=_settings_hash(stale_settings),
        replicate=1,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="chainid A",
        summaries=[
            HydrogenBondReplicateSummary(
                name="protein_polymer",
                mode="between",
                group_names=["protein", "polymer"],
                n_frames_used=3,
                mean_hbonds_per_frame=2.0,
                fraction_frames_with_any_hbond=0.5,
                counts_per_frame=[1, 2, 3],
            )
        ],
    ).save(run_dir / "result.json")

    data = {"CondA": {"analysis_dir": analysis_dir, "replicates": [1]}}

    with pytest.raises(ValueError, match="current settings require"):
        analysis._load_replicate_timeseries(data, ["CondA"], settings=current_settings)


def test_aggregate_composition() -> None:
    """Composition aggregation should compute mean and SEM with zero-fill."""
    analysis = HydrogenBondsAnalysis()
    results = [
        HydrogenBondResult(
            replicate=1,
            summaries=[],
            composition_entries=[
                CompositionEntry(
                    donor_partition="protein",
                    acceptor_partition="polymer",
                    mean_hbonds_per_frame=0.40,
                    fraction_of_total=0.80,
                )
            ],
        ),
        HydrogenBondResult(
            replicate=2,
            summaries=[],
            composition_entries=[
                CompositionEntry(
                    donor_partition="protein",
                    acceptor_partition="polymer",
                    mean_hbonds_per_frame=0.20,
                    fraction_of_total=0.50,
                ),
                CompositionEntry(
                    donor_partition="substrate",
                    acceptor_partition="protein",
                    mean_hbonds_per_frame=0.10,
                    fraction_of_total=0.50,
                ),
            ],
        ),
    ]

    aggregated = analysis._aggregate_composition(results)
    by_pair = {(entry.donor_partition, entry.acceptor_partition): entry for entry in aggregated}

    pp = by_pair[("protein", "polymer")]
    assert isinstance(pp, AggregatedCompositionEntry)
    assert pp.per_replicate_hbonds == pytest.approx([0.40, 0.20])
    assert pp.mean_hbonds_per_frame == pytest.approx(0.30)
    assert pp.sem_hbonds_per_frame == pytest.approx(np.std([0.40, 0.20], ddof=1) / np.sqrt(2))
    assert pp.per_replicate_fraction == pytest.approx([0.80, 0.50])

    sp = by_pair[("substrate", "protein")]
    assert sp.per_replicate_hbonds == pytest.approx([0.0, 0.10])
    assert sp.per_replicate_fraction == pytest.approx([0.0, 0.50])


def test_composition_not_configured(tmp_path: Path) -> None:
    """run_replicate should emit no composition entries when unset."""

    class MockHydrogenBondAnalysis:
        def __init__(self, **kwargs) -> None:
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            self.results.hbonds = np.array([[0, 0, 10, 1, 2.8, 160.0]], dtype=float)

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "C", 100, "OEG", 1),
    }
    universe = MagicMock()
    universe.trajectory = [object(), object(), object()]
    universe.atoms = _MockAtomCollection(atoms)
    selections = {
        "chainid A": _MockAtomGroup([0]),
        "chainid C": _MockAtomGroup([1]),
    }
    universe.select_atoms.side_effect = lambda selection, updating: selections[selection]

    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=HydrogenBondSettings(composition=None),
    )

    analysis = HydrogenBondsAnalysis()
    mock_modules = _make_mdanalysis_module(MockHydrogenBondAnalysis)

    with (
        patch.dict(sys.modules, mock_modules),
        patch("polyzymd.analyses.hydrogen_bonds.TrajectoryLoader") as mock_loader_cls,
        patch("polyzymd.analyses.hydrogen_bonds._mda.compute_config_hash", return_value="abc123"),
    ):
        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader
        mock_loader.load_universe.return_value = universe
        mock_loader.get_timestep.return_value = 10.0

        result = _as_hydrogen_bond_result(analysis.run_replicate(ctx, 1))

    assert result.composition_entries == []


def _make_directed_pair(
    donor: tuple[str, int, str],
    acceptor: tuple[str, int, str],
    occupancy: float,
    events_per_frame: float,
) -> DirectedResiduePairResult:
    return DirectedResiduePairResult(
        donor=ResidueRef(chain_id=donor[0], resid=donor[1], resname=donor[2]),
        acceptor=ResidueRef(chain_id=acceptor[0], resid=acceptor[1], resname=acceptor[2]),
        frames_present=1,
        occupancy=occupancy,
        event_count=1,
        mean_events_per_frame=events_per_frame,
    )


def _make_undirected_pair(
    residue_a: tuple[str, int, str],
    residue_b: tuple[str, int, str],
    occupancy: float,
    events_per_frame: float,
) -> UndirectedResiduePairResult:
    return UndirectedResiduePairResult(
        residue_a=ResidueRef(chain_id=residue_a[0], resid=residue_a[1], resname=residue_a[2]),
        residue_b=ResidueRef(chain_id=residue_b[0], resid=residue_b[1], resname=residue_b[2]),
        frames_present=1,
        occupancy=occupancy,
        event_count=1,
        mean_events_per_frame=events_per_frame,
    )


def _make_replicate_result(
    replicate: int,
    mean_hbonds: float,
    fraction_with_any: float,
    directed_pairs: list[DirectedResiduePairResult] | None = None,
    undirected_pairs: list[UndirectedResiduePairResult] | None = None,
    summaries: list[HydrogenBondReplicateSummary] | None = None,
    settings: HydrogenBondSettings | None = None,
) -> HydrogenBondResult:
    if summaries is None:
        summaries = [
            HydrogenBondReplicateSummary(
                name="protein_polymer",
                mode="between",
                group_names=["protein", "polymer"],
                n_frames_used=10,
                mean_hbonds_per_frame=mean_hbonds,
                fraction_frames_with_any_hbond=fraction_with_any,
                counts_per_frame=[0] * 10,
                directed_residue_pairs=directed_pairs or [],
                undirected_residue_pairs=undirected_pairs or [],
            )
        ]

    return HydrogenBondResult(
        config_hash="cfg123",
        settings_fingerprint=_settings_hash(settings or HydrogenBondSettings()),
        replicate=replicate,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="(chainid A) or (chainid C)",
        summaries=summaries,
    )


def _make_replicate_artifact(
    result: HydrogenBondResult,
    ctx: AggregateContext,
) -> ReplicateArtifact:
    """Wrap a legacy assertion model in a valid hydrogen-bond artifact."""

    store = ArtifactStore(ctx.output_dir.parent / f"run_{result.replicate}")
    sidecar = store.write_npz_sidecar(
        "sidecars/hydrogen_bond_events.npz",
        hbonds=np.empty((0, 6), dtype=float),
        metadata={"kind": "hydrogen_bond_events", "columns": list(HBOND_EVENT_COLUMNS)},
    )
    return ReplicateArtifact(
        analysis_name="hydrogen_bonds",
        condition_label=ctx.condition.label,
        replicate=result.replicate,
        payload={
            "summaries": [summary.model_dump(mode="json") for summary in result.summaries],
            "composition_entries": [
                entry.model_dump(mode="json") for entry in result.composition_entries
            ],
            "metrics": {},
            "replicate_metrics": {},
            "n_frames_used": result.summaries[0].n_frames_used if result.summaries else 0,
            "n_events": 0,
            "event_sidecar": sidecar.path,
        },
        sidecars=[sidecar],
        provenance={"source": "test_legacy_adapter"},
        metadata={
            "result_kind": "hydrogen_bonds_mda_replicate",
            "settings_fingerprint": result.settings_fingerprint,
            "config_hash": result.config_hash,
            "polyzymd_version": result.polyzymd_version,
            "equilibration_time": result.equilibration_time,
            "equilibration_unit": result.equilibration_unit,
            "selection_string": result.selection_string,
            "timestep_ps": result.timestep_ps,
            "raw_timestep_ps": result.raw_timestep_ps,
            "frame_stride": result.frame_stride,
            "event_columns": list(HBOND_EVENT_COLUMNS),
        },
    )


def _write_hbond_condition_artifact(
    aggregated_dir: Path,
    *,
    condition_label: str,
    result: HydrogenBondAggregatedResult,
) -> ConditionArtifact:
    """Write a canonical hydrogen-bond condition artifact for plot tests."""

    artifact = ConditionArtifact(
        analysis_name="hydrogen_bonds",
        condition_label=condition_label,
        replicates=[int(replicate) for replicate in result.replicates],
        payload={
            "summaries": [summary.model_dump(mode="json") for summary in result.summaries],
            "composition_entries": [
                entry.model_dump(mode="json") for entry in result.composition_entries
            ],
            "metrics": {},
            "replicate_metrics": {},
            "n_replicates": result.n_replicates,
        },
        metadata={
            "result_kind": "hydrogen_bonds_mda_condition",
            "settings_fingerprint": result.settings_fingerprint,
            "config_hash": result.config_hash,
            "polyzymd_version": result.polyzymd_version,
            "equilibration_time": result.equilibration_time,
            "equilibration_unit": result.equilibration_unit,
            "selection_string": result.selection_string,
            "timestep_ps": result.timestep_ps,
            "raw_timestep_ps": result.raw_timestep_ps,
            "frame_stride": result.frame_stride,
            "n_replicates": result.n_replicates,
        },
    )
    ArtifactStore(aggregated_dir).write_condition_result(artifact)
    return artifact


def _aggregate_replicate_results(
    analysis: HydrogenBondsAnalysis,
    ctx: AggregateContext,
    results: list[HydrogenBondResult],
) -> HydrogenBondAggregatedResult:
    """Aggregate legacy assertion models through the artifact-only lifecycle."""

    artifacts = [_make_replicate_artifact(result, ctx) for result in results]
    return _as_hydrogen_bond_aggregate(analysis.aggregate(ctx, artifacts))


def _make_aggregate_context(
    tmp_path: Path,
    top_n_pairs: int = 15,
    replicates: tuple[int, ...] = (1, 2, 3),
) -> AggregateContext:
    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=replicates,
        sim_config=MagicMock(),
    )
    settings = HydrogenBondSettings(top_n_pairs=top_n_pairs)
    return AggregateContext(
        condition=condition,
        replicates=replicates,
        output_dir=tmp_path / "aggregated",
        equilibration="10ns",
        settings=settings,
    )


def test_aggregate_basic(tmp_path: Path) -> None:
    """aggregate should compute scalar mean and SEM across replicates."""
    analysis = HydrogenBondsAnalysis()
    ctx = _make_aggregate_context(tmp_path)

    results = [
        _make_replicate_result(1, 2.0, 0.40),
        _make_replicate_result(2, 4.0, 0.50),
        _make_replicate_result(3, 6.0, 0.60),
    ]

    aggregated = _aggregate_replicate_results(analysis, ctx, results)

    assert isinstance(aggregated, HydrogenBondAggregatedResult)
    assert aggregated.n_replicates == 3
    assert len(aggregated.summaries) == 1

    summary = aggregated.summaries[0]
    assert summary.per_replicate_mean_hbonds == [2.0, 4.0, 6.0]
    assert summary.mean_hbonds_per_frame == pytest.approx(4.0)
    assert summary.sem_hbonds_per_frame == pytest.approx(
        np.std([2.0, 4.0, 6.0], ddof=1) / np.sqrt(3)
    )
    assert summary.per_replicate_fraction_with_any == [0.40, 0.50, 0.60]
    assert summary.mean_fraction_with_any == pytest.approx(0.50)
    assert summary.sem_fraction_with_any == pytest.approx(
        np.std([0.40, 0.50, 0.60], ddof=1) / np.sqrt(3)
    )
    assert not (ctx.output_dir / "result.json").exists()


def test_aggregate_rejects_legacy_replicate_results(tmp_path: Path) -> None:
    """aggregate should reject stale legacy replicate models after MDA migration."""
    analysis = HydrogenBondsAnalysis()
    ctx = _make_aggregate_context(tmp_path, replicates=(1,))

    with pytest.raises(TypeError, match="ReplicateArtifact.*recompute=True"):
        _ = analysis.aggregate(ctx, [_make_replicate_result(1, 2.0, 0.40)])


def test_aggregate_propagates_timing_metadata(tmp_path: Path) -> None:
    """Hydrogen-bond aggregate should preserve replicate timing metadata."""
    analysis = HydrogenBondsAnalysis()
    ctx = _make_aggregate_context(tmp_path, replicates=(1, 2))
    results = [
        _make_replicate_result(rep, 1.0 + rep, 0.2).model_copy(
            update={"timestep_ps": 30.0, "raw_timestep_ps": 10.0, "frame_stride": 3}
        )
        for rep in (1, 2)
    ]

    aggregated = _aggregate_replicate_results(analysis, ctx, results)

    assert aggregated.timestep_ps == pytest.approx(30.0)
    assert aggregated.raw_timestep_ps == pytest.approx(10.0)
    assert aggregated.frame_stride == 3


def test_aggregate_pair_merging(tmp_path: Path) -> None:
    """aggregate should merge residue pairs and zero-fill missing replicate values."""
    analysis = HydrogenBondsAnalysis()
    ctx = _make_aggregate_context(tmp_path, replicates=(1, 2))

    donor_a = ("A", 10, "SER")
    acceptor_b = ("C", 100, "OEG")
    acceptor_c = ("C", 101, "OEG")
    acceptor_d = ("C", 102, "OEG")

    rep1 = _make_replicate_result(
        1,
        3.0,
        0.7,
        directed_pairs=[
            _make_directed_pair(donor_a, acceptor_b, occupancy=0.5, events_per_frame=0.5),
            _make_directed_pair(donor_a, acceptor_c, occupancy=0.3, events_per_frame=0.3),
        ],
        undirected_pairs=[
            _make_undirected_pair(donor_a, acceptor_b, occupancy=0.5, events_per_frame=0.5),
            _make_undirected_pair(donor_a, acceptor_c, occupancy=0.3, events_per_frame=0.3),
        ],
    )
    rep2 = _make_replicate_result(
        2,
        2.5,
        0.6,
        directed_pairs=[
            _make_directed_pair(donor_a, acceptor_b, occupancy=0.6, events_per_frame=0.6),
            _make_directed_pair(donor_a, acceptor_d, occupancy=0.2, events_per_frame=0.2),
        ],
        undirected_pairs=[
            _make_undirected_pair(donor_a, acceptor_b, occupancy=0.6, events_per_frame=0.6),
            _make_undirected_pair(donor_a, acceptor_d, occupancy=0.2, events_per_frame=0.2),
        ],
    )
    aggregated = _aggregate_replicate_results(analysis, ctx, [rep1, rep2])
    summary = aggregated.summaries[0]

    directed_by_pair = {
        (pair.donor.label, pair.acceptor.label): pair for pair in summary.directed_pairs
    }
    ab = directed_by_pair[("SER10:A", "OEG100:C")]
    ac = directed_by_pair[("SER10:A", "OEG101:C")]
    ad = directed_by_pair[("SER10:A", "OEG102:C")]

    assert ab.mean_occupancy == pytest.approx(0.55)
    assert ac.per_replicate_occupancy == pytest.approx([0.3, 0.0])
    assert ad.per_replicate_occupancy == pytest.approx([0.0, 0.2])

    undirected_by_pair = {
        tuple(sorted([pair.residue_a.label, pair.residue_b.label])): pair
        for pair in summary.undirected_pairs
    }
    ac_u = undirected_by_pair[tuple(sorted(["SER10:A", "OEG101:C"]))]
    ad_u = undirected_by_pair[tuple(sorted(["SER10:A", "OEG102:C"]))]
    assert ac_u.per_replicate_occupancy == pytest.approx([0.3, 0.0])
    assert ad_u.per_replicate_occupancy == pytest.approx([0.0, 0.2])


def test_aggregate_top_n_truncation(tmp_path: Path) -> None:
    """aggregate should keep only top_n_pairs by mean occupancy."""
    analysis = HydrogenBondsAnalysis()
    ctx = _make_aggregate_context(tmp_path, top_n_pairs=2)

    donor = ("A", 10, "SER")
    pairs = [
        ("C", 100, "OEG", 0.9),
        ("C", 101, "OEG", 0.7),
        ("C", 102, "OEG", 0.5),
    ]

    directed = [
        _make_directed_pair(donor, (chain, resid, resname), occupancy=occ, events_per_frame=occ)
        for chain, resid, resname, occ in pairs
    ]
    undirected = [
        _make_undirected_pair(donor, (chain, resid, resname), occupancy=occ, events_per_frame=occ)
        for chain, resid, resname, occ in pairs
    ]

    results = [
        _make_replicate_result(
            1,
            2.0,
            0.4,
            directed_pairs=directed,
            undirected_pairs=undirected,
            settings=ctx.settings,
        ),
        _make_replicate_result(
            2,
            2.0,
            0.4,
            directed_pairs=directed,
            undirected_pairs=undirected,
            settings=ctx.settings,
        ),
        _make_replicate_result(
            3,
            2.0,
            0.4,
            directed_pairs=directed,
            undirected_pairs=undirected,
            settings=ctx.settings,
        ),
    ]

    aggregated = _aggregate_replicate_results(analysis, ctx, results)
    summary = aggregated.summaries[0]

    assert len(summary.directed_pairs) == 2
    assert len(summary.undirected_pairs) == 2
    assert summary.directed_pairs[0].mean_occupancy >= summary.directed_pairs[1].mean_occupancy


def test_aggregate_error_message_includes_replicate_details(tmp_path: Path) -> None:
    """Aggregate mismatch errors should include expected replicate IDs."""
    analysis = HydrogenBondsAnalysis()
    ctx = _make_aggregate_context(tmp_path, replicates=(1, 2, 3))
    with pytest.raises(ValueError, match=r"missing \[2, 3\]"):
        _ = analysis.aggregate(
            ctx,
            [_make_replicate_artifact(_make_replicate_result(1, 1.0, 0.2), ctx)],
        )


def test_aggregate_rejects_duplicate_replicate_ids(tmp_path: Path) -> None:
    """Aggregate should reject duplicate replicate IDs even when result count matches."""
    analysis = HydrogenBondsAnalysis()
    ctx = _make_aggregate_context(tmp_path, replicates=(1, 2, 3))

    results = [
        _make_replicate_result(1, 1.0, 0.2),
        _make_replicate_result(2, 1.2, 0.3),
        _make_replicate_result(2, 1.4, 0.4),
    ]

    with pytest.raises(
        ValueError,
        match=r"Duplicate hydrogen-bond artifact for replicate 2",
    ):
        _ = analysis.aggregate(
            ctx,
            [_make_replicate_artifact(result, ctx) for result in results],
        )


def test_aggregate_rejects_unexpected_replicate_ids(tmp_path: Path) -> None:
    """Aggregate should reject unexpected replicate IDs instead of trusting list order."""
    analysis = HydrogenBondsAnalysis()
    ctx = _make_aggregate_context(tmp_path, replicates=(1, 2, 3))

    results = [
        _make_replicate_result(1, 1.0, 0.2),
        _make_replicate_result(2, 1.2, 0.3),
        _make_replicate_result(4, 1.4, 0.4),
    ]

    with pytest.raises(
        ValueError,
        match=r"missing \[3\], unexpected \[4\]",
    ):
        _ = analysis.aggregate(
            ctx,
            [_make_replicate_artifact(result, ctx) for result in results],
        )


def test_framework_can_persist_returned_aggregate(tmp_path: Path) -> None:
    """Framework-owned persistence should write returned aggregate artifacts."""
    analysis = HydrogenBondsAnalysis()
    ctx = _make_aggregate_context(tmp_path)
    results = [
        _make_replicate_result(1, 1.0, 0.2),
        _make_replicate_result(2, 1.2, 0.3),
        _make_replicate_result(3, 1.4, 0.4),
    ]

    artifacts = [_make_replicate_artifact(result, ctx) for result in results]
    aggregate = analysis.aggregate(ctx, artifacts)
    ArtifactStore(ctx.output_dir).write_condition_result(aggregate)

    assert (ctx.output_dir / "result.json").exists()


def test_aggregate_rejects_missing_summary_in_one_replicate(tmp_path: Path) -> None:
    """Aggregate should fail loudly on stale replicate outputs missing summaries."""
    analysis = HydrogenBondsAnalysis()
    ctx = _make_aggregate_context(tmp_path, replicates=(1, 2, 3))

    rep1 = _make_replicate_result(1, 2.0, 0.4)
    rep2 = _make_replicate_result(2, 0.0, 0.0, summaries=[])
    rep3 = _make_replicate_result(3, 4.0, 0.6)

    with pytest.raises(
        ValueError,
        match=r"condition 'test'.*summary 'protein_polymer'.*replicate results \[2\]",
    ):
        _ = _aggregate_replicate_results(analysis, ctx, [rep1, rep2, rep3])


def test_aggregate_single_replicate(tmp_path: Path) -> None:
    """Single replicate aggregate should report SEM as 0."""
    analysis = HydrogenBondsAnalysis()
    ctx = _make_aggregate_context(tmp_path, replicates=(1,))

    aggregated = _aggregate_replicate_results(
        analysis,
        ctx,
        [_make_replicate_result(1, 3.5, 0.9)],
    )
    summary = aggregated.summaries[0]

    assert summary.mean_hbonds_per_frame == pytest.approx(3.5)
    assert summary.sem_hbonds_per_frame == pytest.approx(0.0)
    assert summary.mean_fraction_with_any == pytest.approx(0.9)
    assert summary.sem_fraction_with_any == pytest.approx(0.0)


def test_aggregate_empty_results(tmp_path: Path) -> None:
    """Aggregating empty results should raise ValueError."""
    analysis = HydrogenBondsAnalysis()
    ctx = _make_aggregate_context(tmp_path, replicates=())

    with pytest.raises(ValueError, match="no replicate results provided"):
        _ = analysis.aggregate(ctx, [])


def test_aggregate_zero_hbonds(tmp_path: Path) -> None:
    """All-zero replicate values should produce all-zero aggregated summary."""
    analysis = HydrogenBondsAnalysis()
    ctx = _make_aggregate_context(tmp_path, replicates=(1, 2, 3))

    results = [
        _make_replicate_result(1, 0.0, 0.0),
        _make_replicate_result(2, 0.0, 0.0),
        _make_replicate_result(3, 0.0, 0.0),
    ]
    aggregated = _aggregate_replicate_results(analysis, ctx, results)
    summary = aggregated.summaries[0]

    assert summary.mean_hbonds_per_frame == pytest.approx(0.0)
    assert summary.sem_hbonds_per_frame == pytest.approx(0.0)
    assert summary.mean_fraction_with_any == pytest.approx(0.0)
    assert summary.sem_fraction_with_any == pytest.approx(0.0)
    assert summary.per_replicate_mean_hbonds == pytest.approx([0.0, 0.0, 0.0])
    assert summary.per_replicate_fraction_with_any == pytest.approx([0.0, 0.0, 0.0])


def test_aggregate_pair_alignment(tmp_path: Path) -> None:
    """Pair occupancy arrays should align to full replicate order."""
    analysis = HydrogenBondsAnalysis()
    ctx = _make_aggregate_context(tmp_path, replicates=(1, 2, 3))

    donor_a = ("A", 10, "SER")
    acceptor_b = ("C", 100, "OEG")
    acceptor_c = ("C", 101, "OEG")

    rep1 = _make_replicate_result(
        1,
        1.0,
        0.4,
        directed_pairs=[
            _make_directed_pair(donor_a, acceptor_b, occupancy=0.5, events_per_frame=0.5)
        ],
    )
    rep2 = _make_replicate_result(2, 0.0, 0.0, directed_pairs=[])
    rep3 = _make_replicate_result(
        3,
        1.0,
        0.4,
        directed_pairs=[
            _make_directed_pair(donor_a, acceptor_c, occupancy=0.2, events_per_frame=0.2)
        ],
    )

    aggregated = _aggregate_replicate_results(analysis, ctx, [rep1, rep2, rep3])
    summary = aggregated.summaries[0]
    directed_by_pair = {
        (pair.donor.label, pair.acceptor.label): pair for pair in summary.directed_pairs
    }

    assert directed_by_pair[("SER10:A", "OEG100:C")].per_replicate_occupancy == pytest.approx(
        [0.5, 0.0, 0.0]
    )
    assert directed_by_pair[("SER10:A", "OEG101:C")].per_replicate_occupancy == pytest.approx(
        [0.0, 0.0, 0.2]
    )


def test_composition_overlap_counts_all_partition_combinations() -> None:
    """Overlapping atoms should contribute to all matching partition pairs."""
    analysis = HydrogenBondsAnalysis()
    composition_settings = HydrogenBondCompositionSettings(
        partitions={
            "z_partition": "chainid A",
            "a_partition": "chainid C",
        }
    )

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "C", 100, "OEG", 1),
    }
    universe = MagicMock()
    universe.atoms = _MockAtomCollection(atoms)
    universe.select_atoms.side_effect = lambda selection: {
        "chainid A": _MockAtomGroup([0]),
        "chainid C": _MockAtomGroup([0, 1]),
    }[selection]

    entries = analysis._compute_composition(
        composition_settings=composition_settings,
        hbond_array=np.array([[0, 0, 10, 1, 2.8, 160.0]], dtype=float),
        universe=universe,
        start_frame=0,
        n_frames=2,
        allow_overlapping=True,
    )

    by_pair = {(entry.donor_partition, entry.acceptor_partition): entry for entry in entries}
    assert len(by_pair) == 2
    assert ("z_partition", "a_partition") in by_pair
    assert ("a_partition", "a_partition") in by_pair


def _make_aggregated_summary(name: str, mean: float, sem: float, reps: list[float]) -> Any:
    return HydrogenBondAggregatedSummary(
        name=name,
        mode="between",
        group_names=["protein", "polymer"],
        n_replicates=len(reps),
        mean_hbonds_per_frame=mean,
        sem_hbonds_per_frame=sem,
        per_replicate_mean_hbonds=reps,
        mean_fraction_with_any=0.0,
        sem_fraction_with_any=0.0,
        per_replicate_fraction_with_any=[0.0 for _ in reps],
    )


def _make_aggregated_result_for_metrics() -> HydrogenBondAggregatedResult:
    return HydrogenBondAggregatedResult(
        replicates=[1, 2, 3],
        n_replicates=3,
        summaries=[
            _make_aggregated_summary("protein_polymer", 2.5, 0.2, [2.2, 2.5, 2.8]),
            _make_aggregated_summary("protein_internal", 1.2, 0.1, [1.0, 1.3, 1.3]),
        ],
    )


def test_extract_metrics_basic() -> None:
    """extract_metrics should return one MetricValue with expected fields."""
    analysis = HydrogenBondsAnalysis()
    summary = HydrogenBondAggregatedResult(
        replicates=[1, 2],
        n_replicates=2,
        summaries=[_make_aggregated_summary("protein_polymer", 3.0, 0.3, [2.8, 3.2])],
    )

    metrics = analysis.extract_metrics(summary)

    assert set(metrics) == {"mean_hbonds_protein_polymer"}
    metric = metrics["mean_hbonds_protein_polymer"]
    assert metric.mean == pytest.approx(3.0)
    assert metric.sem == pytest.approx(0.3)
    assert metric.replicate_values == pytest.approx([2.8, 3.2])
    assert metric.higher_is_better is None
    assert metric.direction_labels == ("fewer H-bonds", "similar", "more H-bonds")


def test_direction_labels_correct() -> None:
    """direction_labels should map negative and positive changes correctly."""
    analysis = HydrogenBondsAnalysis()
    metric = analysis.extract_metrics(_make_aggregated_result_for_metrics())[
        "mean_hbonds_protein_polymer"
    ]

    assert interpret_direction(-5.0, metric.direction_labels) == "fewer H-bonds"
    assert interpret_direction(0.0, metric.direction_labels) == "similar"
    assert interpret_direction(5.0, metric.direction_labels) == "more H-bonds"


def test_extract_metrics_multiple_summaries() -> None:
    """extract_metrics should return one metric per aggregated summary."""
    analysis = HydrogenBondsAnalysis()
    metrics = analysis.extract_metrics(_make_aggregated_result_for_metrics())

    assert set(metrics) == {
        "mean_hbonds_protein_polymer",
        "mean_hbonds_protein_internal",
    }


def test_extract_metrics_dict_input() -> None:
    """extract_metrics should support dict-style aggregated summaries."""
    analysis = HydrogenBondsAnalysis()
    metrics = analysis.extract_metrics(
        {
            "summaries": [
                {
                    "name": "protein_polymer",
                    "mean_hbonds_per_frame": 2.0,
                    "sem_hbonds_per_frame": 0.25,
                    "per_replicate_mean_hbonds": [1.8, 2.2],
                }
            ]
        }
    )

    metric = metrics["mean_hbonds_protein_polymer"]
    assert metric.mean == pytest.approx(2.0)
    assert metric.sem == pytest.approx(0.25)
    assert metric.replicate_values == pytest.approx([1.8, 2.2])


def test_format_comparison_result() -> None:
    """format should return non-empty text for ComparisonResult inputs."""
    analysis = HydrogenBondsAnalysis()
    result = ComparisonResult(
        analysis_type="hydrogen_bonds",
        name="hbonds_cmp",
        conditions=[
            ConditionSummary(
                label="Control",
                n_replicates=3,
                mean_hbonds_protein_polymer_mean=2.0,
                mean_hbonds_protein_polymer_sem=0.1,
            ),
            ConditionSummary(
                label="Treatment",
                n_replicates=3,
                mean_hbonds_protein_polymer_mean=2.5,
                mean_hbonds_protein_polymer_sem=0.1,
            ),
        ],
        pairwise_comparisons=[
            PairwiseResult(
                condition_a="Control",
                condition_b="Treatment",
                metric="mean_hbonds_protein_polymer",
                t_statistic=1.0,
                p_value=0.2,
                cohens_d=0.4,
                effect_size_interpretation="small",
                direction="more H-bonds",
                significant=False,
                percent_change=25.0,
            )
        ],
        ranking=["Treatment", "Control"],
        rankings_by_metric={"mean_hbonds_protein_polymer": ["Treatment", "Control"]},
        equilibration_time="10ns",
        created_at="2026-01-01T00:00:00",
        polyzymd_version="test",
    )

    text = analysis.format(result, output_format="text")

    assert isinstance(text, str)
    assert text.strip()
    assert "Hydrogen Bond Analysis" in text
    assert "no direction preference" in text


def test_format_multi_summary() -> None:
    """format should isolate per-metric ranking and statistics by section."""
    analysis = HydrogenBondsAnalysis()
    result = ComparisonResult(
        analysis_type="hydrogen_bonds",
        name="hbonds_cmp",
        conditions=[
            ConditionSummary(
                label="Control",
                n_replicates=3,
                mean_hbonds_protein_polymer_mean=2.0,
                mean_hbonds_protein_polymer_sem=0.1,
                mean_hbonds_protein_internal_mean=1.0,
                mean_hbonds_protein_internal_sem=0.1,
            ),
            ConditionSummary(
                label="Treatment",
                n_replicates=3,
                mean_hbonds_protein_polymer_mean=2.5,
                mean_hbonds_protein_polymer_sem=0.1,
                mean_hbonds_protein_internal_mean=1.4,
                mean_hbonds_protein_internal_sem=0.1,
            ),
        ],
        pairwise_comparisons=[
            PairwiseResult(
                condition_a="Control",
                condition_b="Treatment",
                metric="mean_hbonds_protein_polymer",
                t_statistic=1.0,
                p_value=0.22,
                cohens_d=0.4,
                effect_size_interpretation="small",
                direction="more H-bonds",
                significant=False,
                percent_change=25.0,
            ),
            PairwiseResult(
                condition_a="Control",
                condition_b="Treatment",
                metric="mean_hbonds_protein_internal",
                t_statistic=1.2,
                p_value=0.01,
                cohens_d=0.5,
                effect_size_interpretation="small",
                direction="more H-bonds",
                significant=True,
                percent_change=40.0,
            ),
        ],
        ranking=["Treatment", "Control"],
        anova=[
            ANOVAResult(
                metric="mean_hbonds_protein_polymer",
                f_statistic=2.1,
                p_value=0.20,
                significant=False,
            ),
            ANOVAResult(
                metric="mean_hbonds_protein_internal",
                f_statistic=7.7,
                p_value=0.01,
                significant=True,
            ),
        ],
        rankings_by_metric={
            "mean_hbonds_protein_polymer": ["Treatment", "Control"],
            "mean_hbonds_protein_internal": ["Control", "Treatment"],
        },
        equilibration_time="10ns",
        created_at="2026-01-01T00:00:00",
        polyzymd_version="test",
    )

    text = analysis.format(result, output_format="text")

    assert "H-bonds: protein_polymer" in text
    assert "H-bonds: protein_internal" in text

    polymer_section = text.split("H-bonds: protein_polymer", maxsplit=1)[1].split(
        "H-bonds: protein_internal", maxsplit=1
    )[0]
    internal_section = text.split("H-bonds: protein_internal", maxsplit=1)[1]

    assert "Highest value: Treatment" in polymer_section
    assert "Highest value: Control" not in polymer_section
    assert "0.2200" in polymer_section
    assert "0.0100" not in polymer_section
    assert "Metric: mean_hbonds_protein_polymer" in polymer_section
    assert "Metric: mean_hbonds_protein_internal" not in polymer_section

    assert "Highest value: Control" in internal_section
    assert "Highest value: Treatment" not in internal_section
    assert "0.0100" in internal_section
    assert "0.2200" not in internal_section
    assert "Metric: mean_hbonds_protein_internal" in internal_section
    assert "Metric: mean_hbonds_protein_polymer" not in internal_section


def test_format_multi_summary_markdown() -> None:
    """Markdown path should also isolate per-metric ranking and statistics."""
    analysis = HydrogenBondsAnalysis()
    result = ComparisonResult(
        analysis_type="hydrogen_bonds",
        name="hbonds_cmp",
        conditions=[
            ConditionSummary(
                label="Control",
                n_replicates=3,
                mean_hbonds_protein_polymer_mean=2.0,
                mean_hbonds_protein_polymer_sem=0.1,
                mean_hbonds_protein_internal_mean=1.0,
                mean_hbonds_protein_internal_sem=0.1,
            ),
            ConditionSummary(
                label="Treatment",
                n_replicates=3,
                mean_hbonds_protein_polymer_mean=2.5,
                mean_hbonds_protein_polymer_sem=0.1,
                mean_hbonds_protein_internal_mean=1.4,
                mean_hbonds_protein_internal_sem=0.1,
            ),
        ],
        pairwise_comparisons=[
            PairwiseResult(
                condition_a="Control",
                condition_b="Treatment",
                metric="mean_hbonds_protein_polymer",
                t_statistic=1.0,
                p_value=0.22,
                cohens_d=0.4,
                effect_size_interpretation="small",
                direction="more H-bonds",
                significant=False,
                percent_change=25.0,
            ),
            PairwiseResult(
                condition_a="Control",
                condition_b="Treatment",
                metric="mean_hbonds_protein_internal",
                t_statistic=1.2,
                p_value=0.01,
                cohens_d=0.5,
                effect_size_interpretation="small",
                direction="more H-bonds",
                significant=True,
                percent_change=40.0,
            ),
        ],
        ranking=["Treatment", "Control"],
        anova=[
            ANOVAResult(
                metric="mean_hbonds_protein_polymer",
                f_statistic=2.1,
                p_value=0.20,
                significant=False,
            ),
            ANOVAResult(
                metric="mean_hbonds_protein_internal",
                f_statistic=7.7,
                p_value=0.01,
                significant=True,
            ),
        ],
        rankings_by_metric={
            "mean_hbonds_protein_polymer": ["Treatment", "Control"],
            "mean_hbonds_protein_internal": ["Control", "Treatment"],
        },
        equilibration_time="10ns",
        created_at="2026-01-01T00:00:00",
        polyzymd_version="test",
    )

    md = analysis.format(result, output_format="markdown")

    assert "H-bonds: protein_polymer" in md
    assert "H-bonds: protein_internal" in md

    polymer_section = md.split("H-bonds: protein_polymer", maxsplit=1)[1].split(
        "H-bonds: protein_internal", maxsplit=1
    )[0]
    internal_section = md.split("H-bonds: protein_internal", maxsplit=1)[1]

    assert "**Highest value:** Treatment" in polymer_section
    assert "**Highest value:** Control" not in polymer_section

    assert "**Highest value:** Control" in internal_section
    assert "**Highest value:** Treatment" not in internal_section

    # P-values should be isolated to their respective sections
    assert "0.22" in polymer_section
    assert "0.01" not in polymer_section

    assert "0.01" in internal_section
    assert "0.22" not in internal_section


def test_format_neutral_metric_keeps_condition_rows_and_neutral_language() -> None:
    """Neutral metrics should render condition rows without best-language."""
    analysis = HydrogenBondsAnalysis()

    cond_a = HydrogenBondAggregatedResult(
        replicates=[1, 2, 3],
        n_replicates=3,
        summaries=[_make_aggregated_summary("protein_polymer", 2.0, 0.1, [1.8, 2.0, 2.2])],
    )
    cond_b = HydrogenBondAggregatedResult(
        replicates=[1, 2, 3],
        n_replicates=3,
        summaries=[_make_aggregated_summary("protein_polymer", 3.0, 0.1, [2.8, 3.0, 3.2])],
    )

    comparison = default_scalar_comparison(
        analysis_name="hydrogen_bonds",
        project_name="hbonds_cmp",
        metrics_by_condition={
            "Alpha": analysis.extract_metrics(cond_a),
            "Beta": analysis.extract_metrics(cond_b),
        },
        control_label="Alpha",
        equilibration="10ns",
    )

    text = analysis.format(comparison, output_format="text")

    assert "Condition Summary" in text
    assert "Alpha" in text
    assert "Beta" in text
    assert "no direction preference" in text
    assert "Best" not in text
    assert "best" not in text


def test_singleton_hydrogen_bond_comparison_not_testable() -> None:
    """Hydrogen-bond default comparison should mark singleton runs not testable."""
    analysis = HydrogenBondsAnalysis()
    cond_a = HydrogenBondAggregatedResult(
        replicates=[1],
        n_replicates=1,
        summaries=[_make_aggregated_summary("protein_polymer", 2.0, 0.0, [2.0])],
    )
    cond_b = HydrogenBondAggregatedResult(
        replicates=[1],
        n_replicates=1,
        summaries=[_make_aggregated_summary("protein_polymer", 3.0, 0.0, [3.0])],
    )

    comparison = default_scalar_comparison(
        analysis_name="hydrogen_bonds",
        project_name="hbonds_singleton",
        metrics_by_condition={
            "Alpha": analysis.extract_metrics(cond_a),
            "Beta": analysis.extract_metrics(cond_b),
        },
        control_label="Alpha",
        equilibration="10ns",
    )

    pairwise = comparison.pairwise_comparisons[0]
    assert pairwise.testable is False
    assert pairwise.significant is False
    assert pairwise.note is not None
    assert "not testable" in analysis.format(comparison, output_format="text")


def test_format_non_comparison_result() -> None:
    """format should fall back to base formatter for non-ComparisonResult objects."""
    analysis = HydrogenBondsAnalysis()
    payload = {"hello": "world"}

    text = analysis.format(payload, output_format="text")

    assert text == str(payload)


def test_compare_rejects_stale_preloaded_aggregated_result(tmp_path: Path) -> None:
    """compare should fail loudly when aggregated results use stale settings."""
    analysis = HydrogenBondsAnalysis()
    current_settings = HydrogenBondSettings()
    stale_settings = HydrogenBondSettings(angle_cutoff=140.0)
    condition = Condition(
        label="CondA",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1, 2),
        sim_config=MagicMock(),
    )

    ctx = ComparisonContext(
        name="hbonds_compare",
        conditions=[condition],
        excluded_conditions=[],
        control_label="CondA",
        analysis_dirs={"CondA": tmp_path / "analysis" / "conda" / "hydrogen_bonds"},
        results_dir=tmp_path / "comparison",
        equilibration="10ns",
        settings=current_settings,
        recompute=False,
        aggregated_results={
            "CondA": HydrogenBondAggregatedResult(
                settings_fingerprint=_settings_hash(stale_settings),
                replicates=[1, 2],
                n_replicates=2,
                summaries=[_make_aggregated_summary("protein_polymer", 2.0, 0.2, [1.8, 2.2])],
            )
        },
    )

    with pytest.raises(ValueError, match="current settings require"):
        analysis.compare(ctx)


def test_plot_rejects_stale_aggregated_result_from_disk(tmp_path: Path) -> None:
    """plot should fail loudly instead of silently plotting stale aggregated caches."""
    analysis = HydrogenBondsAnalysis()
    current_settings = HydrogenBondSettings()
    stale_settings = HydrogenBondSettings(groups={"protein": "chainid A", "polymer": "chainid B"})
    condition = Condition(
        label="CondA",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1, 2),
        sim_config=MagicMock(),
    )
    analysis_dir = tmp_path / "analysis" / "conda" / "hydrogen_bonds"
    aggregated_dir = analysis_dir / "aggregated"
    aggregated_dir.mkdir(parents=True)
    _write_hbond_condition_artifact(
        aggregated_dir,
        condition_label="CondA",
        result=HydrogenBondAggregatedResult(
            settings_fingerprint=_settings_hash(stale_settings),
            replicates=[1, 2],
            n_replicates=2,
            summaries=[_make_aggregated_summary("protein_polymer", 2.0, 0.2, [1.8, 2.2])],
        ),
    )

    ctx = PlotContext(
        conditions=[condition],
        analysis_dirs={"CondA": analysis_dir},
        results_dir=tmp_path / "comparison",
        output_dir=tmp_path / "figures",
        settings=current_settings,
        plot_settings=PlotSettings(),
    )

    with pytest.raises(ValueError, match="current settings require"):
        analysis.plot(ctx)


def test_plot_rejects_legacy_replicate_cache_from_disk(tmp_path: Path) -> None:
    """plot should fail loudly instead of silently plotting legacy replicate caches."""
    analysis = HydrogenBondsAnalysis()
    current_settings = HydrogenBondSettings()
    condition = Condition(
        label="CondA",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    analysis_dir = tmp_path / "analysis" / "conda" / "hydrogen_bonds"
    aggregated_dir = analysis_dir / "aggregated"
    run_dir = analysis_dir / "run_1"
    aggregated_dir.mkdir(parents=True)
    run_dir.mkdir(parents=True)

    _write_hbond_condition_artifact(
        aggregated_dir,
        condition_label="CondA",
        result=HydrogenBondAggregatedResult(
            settings_fingerprint=_settings_hash(current_settings),
            replicates=[1],
            n_replicates=1,
            summaries=[_make_aggregated_summary("protein_polymer", 2.0, 0.2, [2.0])],
        ),
    )

    HydrogenBondResult(
        config_hash="cfg123",
        replicate=1,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="chainid A",
        summaries=[
            HydrogenBondReplicateSummary(
                name="protein_polymer",
                mode="between",
                group_names=["protein", "polymer"],
                n_frames_used=3,
                mean_hbonds_per_frame=2.0,
                fraction_frames_with_any_hbond=0.5,
                counts_per_frame=[1, 2, 3],
            )
        ],
    ).save(run_dir / "result.json")

    ctx = PlotContext(
        conditions=[condition],
        analysis_dirs={"CondA": analysis_dir},
        results_dir=tmp_path / "comparison",
        output_dir=tmp_path / "figures",
        settings=current_settings,
        plot_settings=PlotSettings(),
    )

    with pytest.raises(ValueError, match="missing a settings fingerprint"):
        analysis.plot(ctx)


def test_plot_returns_paths(tmp_path: Path) -> None:
    """plot should return paths emitted by plotter functions."""
    analysis = HydrogenBondsAnalysis()
    condition = Condition(
        label="CondA",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1, 2),
        sim_config=MagicMock(),
    )
    analysis_dir = tmp_path / "analysis" / "conda" / "hydrogen_bonds"
    analysis_dir.mkdir(parents=True)
    settings = HydrogenBondSettings()

    ctx = PlotContext(
        conditions=[condition],
        analysis_dirs={"CondA": analysis_dir},
        results_dir=tmp_path / "comparison",
        output_dir=tmp_path / "figures",
        settings=settings,
        plot_settings=PlotSettings(),
    )

    loaded_result = HydrogenBondAggregatedResult(
        settings_fingerprint=_settings_hash(settings),
        replicates=[1, 2],
        n_replicates=2,
        summaries=[_make_aggregated_summary("protein_polymer", 2.0, 0.2, [1.8, 2.2])],
        composition_entries=[
            AggregatedCompositionEntry(
                donor_partition="protein",
                acceptor_partition="polymer",
                mean_hbonds_per_frame=1.1,
                sem_hbonds_per_frame=0.1,
                per_replicate_hbonds=[1.0, 1.2],
                mean_fraction_of_total=0.6,
                sem_fraction_of_total=0.05,
                per_replicate_fraction=[0.55, 0.65],
            )
        ],
    )
    _write_hbond_condition_artifact(
        analysis_dir / "aggregated",
        condition_label="CondA",
        result=loaded_result,
    )

    with (
        patch.object(analysis, "_load_aggregated_result", return_value=loaded_result),
        patch.object(
            analysis,
            "_load_replicate_timeseries",
            return_value={"CondA": {"protein_polymer": [[1, 2, 3]]}},
        ),
        patch(
            "polyzymd.analyses.hydrogen_bonds._plotters.plot_summary_comparison",
            return_value=tmp_path / "figures" / "hbond_summary_comparison.png",
        ),
        patch(
            "polyzymd.analyses.hydrogen_bonds._plotters.plot_timeseries",
            return_value=tmp_path / "figures" / "hbond_timeseries_protein_polymer.png",
        ),
        patch(
            "polyzymd.analyses.hydrogen_bonds._plotters.plot_top_pairs",
            return_value=tmp_path / "figures" / "hbond_top_pairs_protein_polymer.png",
        ),
        patch(
            "polyzymd.analyses.hydrogen_bonds._plotters.plot_composition_absolute",
            autospec=True,
            return_value=tmp_path / "figures" / "hbond_composition_absolute.png",
        ),
        patch(
            "polyzymd.analyses.hydrogen_bonds._plotters.plot_composition_fraction",
            autospec=True,
            return_value=tmp_path / "figures" / "hbond_composition_fraction.png",
        ),
    ):
        plots = analysis.plot(ctx)

    assert all(isinstance(path, Path) for path in plots)
    assert len(plots) == 5


def test_plot_empty_data(tmp_path: Path) -> None:
    """plot should return empty list when no aggregated results are loaded."""
    analysis = HydrogenBondsAnalysis()
    condition = Condition(
        label="CondA",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    analysis_dir = tmp_path / "analysis" / "conda" / "hydrogen_bonds"
    analysis_dir.mkdir(parents=True)

    ctx = PlotContext(
        conditions=[condition],
        analysis_dirs={"CondA": analysis_dir},
        results_dir=tmp_path / "comparison",
        output_dir=tmp_path / "figures",
        settings=HydrogenBondSettings(),
        plot_settings=PlotSettings(),
    )

    with patch.object(analysis, "_load_aggregated_result", return_value=None):
        plots = analysis.plot(ctx)

    assert plots == []


def test_plot_filters_none_paths(tmp_path: Path) -> None:
    """plot should drop None values returned by plotter helpers."""
    analysis = HydrogenBondsAnalysis()
    condition = Condition(
        label="CondA",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    analysis_dir = tmp_path / "analysis" / "conda" / "hydrogen_bonds"
    analysis_dir.mkdir(parents=True)

    ctx = PlotContext(
        conditions=[condition],
        analysis_dirs={"CondA": analysis_dir},
        results_dir=tmp_path / "comparison",
        output_dir=tmp_path / "figures",
        settings=HydrogenBondSettings(),
        plot_settings=PlotSettings(),
    )

    loaded_result = HydrogenBondAggregatedResult(
        replicates=[1],
        n_replicates=1,
        summaries=[_make_aggregated_summary("protein_polymer", 2.0, 0.2, [2.0])],
    )
    expected_path = tmp_path / "figures" / "hbond_timeseries_protein_polymer.png"

    with (
        patch.object(analysis, "_load_aggregated_result", return_value=loaded_result),
        patch.object(
            analysis,
            "_load_replicate_timeseries",
            return_value={"CondA": {"protein_polymer": [[1, 2, 3]]}},
        ),
        patch(
            "polyzymd.analyses.hydrogen_bonds._plotters.plot_summary_comparison",
            return_value=None,
        ),
        patch(
            "polyzymd.analyses.hydrogen_bonds._plotters.plot_timeseries",
            return_value=expected_path,
        ),
        patch(
            "polyzymd.analyses.hydrogen_bonds._plotters.plot_top_pairs",
            return_value=None,
        ),
    ):
        plots = analysis.plot(ctx)

    assert plots == [expected_path]


def test_plot_composition_gating(tmp_path: Path) -> None:
    """Composition plots should run only when composition entries exist."""
    analysis = HydrogenBondsAnalysis()
    condition = Condition(
        label="CondA",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    analysis_dir = tmp_path / "analysis" / "conda" / "hydrogen_bonds"
    analysis_dir.mkdir(parents=True)

    ctx = PlotContext(
        conditions=[condition],
        analysis_dirs={"CondA": analysis_dir},
        results_dir=tmp_path / "comparison",
        output_dir=tmp_path / "figures",
        settings=HydrogenBondSettings(),
        plot_settings=PlotSettings(),
    )

    no_comp_result = HydrogenBondAggregatedResult(
        replicates=[1],
        n_replicates=1,
        summaries=[_make_aggregated_summary("protein_polymer", 2.0, 0.2, [2.0])],
    )
    with (
        patch.object(analysis, "_load_aggregated_result", return_value=no_comp_result),
        patch.object(
            analysis,
            "_load_replicate_timeseries",
            return_value={"CondA": {"protein_polymer": [[1, 2, 3]]}},
        ),
        patch(
            "polyzymd.analyses.hydrogen_bonds._plotters.plot_summary_comparison", return_value=None
        ),
        patch("polyzymd.analyses.hydrogen_bonds._plotters.plot_timeseries", return_value=None),
        patch("polyzymd.analyses.hydrogen_bonds._plotters.plot_top_pairs", return_value=None),
        patch(
            "polyzymd.analyses.hydrogen_bonds._plotters.plot_composition_absolute",
            autospec=True,
        ) as mock_comp_abs,
        patch(
            "polyzymd.analyses.hydrogen_bonds._plotters.plot_composition_fraction",
            autospec=True,
        ) as mock_comp_frac,
    ):
        _ = analysis.plot(ctx)

    mock_comp_abs.assert_not_called()
    mock_comp_frac.assert_not_called()

    with_comp_result = HydrogenBondAggregatedResult(
        replicates=[1],
        n_replicates=1,
        summaries=[_make_aggregated_summary("protein_polymer", 2.0, 0.2, [2.0])],
        composition_entries=[
            AggregatedCompositionEntry(
                donor_partition="protein",
                acceptor_partition="polymer",
                mean_hbonds_per_frame=1.0,
                sem_hbonds_per_frame=0.1,
                per_replicate_hbonds=[1.0],
                mean_fraction_of_total=1.0,
                sem_fraction_of_total=0.0,
                per_replicate_fraction=[1.0],
            )
        ],
    )
    with (
        patch.object(analysis, "_load_aggregated_result", return_value=with_comp_result),
        patch.object(
            analysis,
            "_load_replicate_timeseries",
            return_value={"CondA": {"protein_polymer": [[1, 2, 3]]}},
        ),
        patch(
            "polyzymd.analyses.hydrogen_bonds._plotters.plot_summary_comparison", return_value=None
        ),
        patch("polyzymd.analyses.hydrogen_bonds._plotters.plot_timeseries", return_value=None),
        patch("polyzymd.analyses.hydrogen_bonds._plotters.plot_top_pairs", return_value=None),
        patch(
            "polyzymd.analyses.hydrogen_bonds._plotters.plot_composition_absolute",
            autospec=True,
            return_value=None,
        ) as mock_comp_abs,
        patch(
            "polyzymd.analyses.hydrogen_bonds._plotters.plot_composition_fraction",
            autospec=True,
            return_value=None,
        ) as mock_comp_frac,
    ):
        _ = analysis.plot(ctx)

    mock_comp_abs.assert_called_once()
    mock_comp_frac.assert_called_once()


def test_plot_top_pairs_receives_top_n(tmp_path: Path) -> None:
    """plot should pass settings.top_n_pairs to plot_top_pairs."""
    analysis = HydrogenBondsAnalysis()
    condition = Condition(
        label="CondA",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    analysis_dir = tmp_path / "analysis" / "conda" / "hydrogen_bonds"
    analysis_dir.mkdir(parents=True)

    settings = HydrogenBondSettings(top_n_pairs=7)
    ctx = PlotContext(
        conditions=[condition],
        analysis_dirs={"CondA": analysis_dir},
        results_dir=tmp_path / "comparison",
        output_dir=tmp_path / "figures",
        settings=settings,
        plot_settings=PlotSettings(),
    )

    loaded_result = HydrogenBondAggregatedResult(
        replicates=[1],
        n_replicates=1,
        summaries=[_make_aggregated_summary("protein_polymer", 2.0, 0.2, [2.0])],
    )

    with (
        patch.object(analysis, "_load_aggregated_result", return_value=loaded_result),
        patch.object(
            analysis,
            "_load_replicate_timeseries",
            return_value={"CondA": {"protein_polymer": [[1, 2, 3]]}},
        ),
        patch(
            "polyzymd.analyses.hydrogen_bonds._plotters.plot_summary_comparison", return_value=None
        ),
        patch("polyzymd.analyses.hydrogen_bonds._plotters.plot_timeseries", return_value=None),
        patch(
            "polyzymd.analyses.hydrogen_bonds._plotters.plot_top_pairs",
            return_value=None,
        ) as mock_top_pairs,
    ):
        _ = analysis.plot(ctx)

    mock_top_pairs.assert_called_once()
    assert mock_top_pairs.call_args.kwargs["top_n"] == 7


def test_settings_within_mode() -> None:
    """Settings should validate with a within-mode summary."""
    settings = HydrogenBondSettings(
        groups={"protein": "chainid A"},
        summaries=[HydrogenBondSummarySettings(name="protein_internal", within="protein")],
    )
    assert settings.summaries[0].within == "protein"
    assert settings.summaries[0].between is None


def test_settings_empty_groups() -> None:
    """Settings should allow empty groups when summaries do not reference them."""
    settings = HydrogenBondSettings(groups={}, summaries=[])
    assert settings.groups == {}
    assert settings.summaries == []


def test_settings_custom_cutoffs() -> None:
    """Settings should accept non-default geometric cutoffs."""
    settings = HydrogenBondSettings(distance_cutoff=3.7, angle_cutoff=175.0)
    assert settings.distance_cutoff == pytest.approx(3.7)
    assert settings.angle_cutoff == pytest.approx(175.0)


def test_composition_settings_empty_partitions() -> None:
    """Composition settings should allow empty partitions."""
    settings = HydrogenBondCompositionSettings(partitions={})
    assert settings.partitions == {}


def test_result_empty_summaries() -> None:
    """Replicate result should support an empty summaries list."""
    result = HydrogenBondResult(summaries=[])
    assert result.summaries == []
    assert "0 summaries" in result.summary()


def test_aggregated_result_zero_replicates_field() -> None:
    """Aggregated result should allow a zero-replicate metadata edge case."""
    result = HydrogenBondAggregatedResult(replicates=[], n_replicates=0, summaries=[])
    assert result.n_replicates == 0
    assert result.replicates == []


def test_residue_ref_equality() -> None:
    """ResidueRef models with identical fields should compare equal."""
    residue_a = ResidueRef(chain_id="A", resid=10, resname="SER")
    residue_b = ResidueRef(chain_id="A", resid=10, resname="SER")
    assert residue_a == residue_b


def test_directed_pair_occupancy_bounds() -> None:
    """Directed pair occupancy should be constrained to [0, 1]."""
    with pytest.raises(ValidationError):
        DirectedResiduePairResult(
            donor=ResidueRef(chain_id="A", resid=10, resname="SER"),
            acceptor=ResidueRef(chain_id="C", resid=100, resname="OEG"),
            frames_present=1,
            occupancy=-0.1,
            event_count=1,
            mean_events_per_frame=0.1,
        )

    with pytest.raises(ValidationError):
        DirectedResiduePairResult(
            donor=ResidueRef(chain_id="A", resid=10, resname="SER"),
            acceptor=ResidueRef(chain_id="C", resid=100, resname="OEG"),
            frames_present=1,
            occupancy=1.1,
            event_count=1,
            mean_events_per_frame=0.1,
        )


def test_undirected_pair_result_roundtrip(tmp_path: Path) -> None:
    """UndirectedResiduePairResult should round-trip through JSON file I/O."""
    pair = UndirectedResiduePairResult(
        residue_a=ResidueRef(chain_id="A", resid=10, resname="SER"),
        residue_b=ResidueRef(chain_id="C", resid=100, resname="OEG"),
        frames_present=4,
        occupancy=0.4,
        event_count=4,
        mean_events_per_frame=0.4,
    )
    file_path = tmp_path / "undirected_pair.json"
    file_path.write_text(pair.model_dump_json(indent=2), encoding="utf-8")
    loaded = UndirectedResiduePairResult.model_validate_json(file_path.read_text(encoding="utf-8"))
    assert loaded == pair


def test_compute_no_hbonds_found(tmp_path: Path) -> None:
    """run_replicate should return zero summaries when no events are found."""

    class MockHydrogenBondAnalysis:
        def __init__(self, **kwargs) -> None:
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            self.results.hbonds = np.empty((0, 6), dtype=float)

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "C", 100, "OEG", 1),
    }
    universe = MagicMock()
    universe.trajectory = [object(), object(), object(), object()]
    universe.atoms = _MockAtomCollection(atoms)
    universe.select_atoms.side_effect = lambda selection, updating: {
        "chainid A": _MockAtomGroup([0]),
        "chainid C": _MockAtomGroup([1]),
    }[selection]

    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=HydrogenBondSettings(),
    )

    analysis = HydrogenBondsAnalysis()
    mock_modules = _make_mdanalysis_module(MockHydrogenBondAnalysis)

    with (
        patch.dict(sys.modules, mock_modules),
        patch("polyzymd.analyses.hydrogen_bonds.TrajectoryLoader") as mock_loader_cls,
        patch("polyzymd.analyses.hydrogen_bonds._mda.compute_config_hash", return_value="abc123"),
    ):
        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader
        mock_loader.load_universe.return_value = universe
        mock_loader.get_timestep.return_value = 10.0

        result = _as_hydrogen_bond_result(analysis.run_replicate(ctx, 1))

    summary = result.summaries[0]
    assert summary.mean_hbonds_per_frame == pytest.approx(0.0)
    assert summary.fraction_frames_with_any_hbond == pytest.approx(0.0)
    assert summary.counts_per_frame == [0, 0, 0, 0]


def test_compute_between_mode_classification(tmp_path: Path) -> None:
    """Between-mode should classify donor→acceptor events in both directions."""

    class MockHydrogenBondAnalysis:
        def __init__(self, **kwargs) -> None:
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            self.results.hbonds = np.array(
                [
                    [0, 0, 10, 1, 2.8, 160.0],
                    [1, 1, 10, 0, 2.9, 158.0],
                    [2, 0, 10, 2, 2.9, 158.0],
                ],
                dtype=float,
            )

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "C", 100, "OEG", 1),
        2: _MockAtom(2, "A", 11, "THR", 2),
    }
    universe = MagicMock()
    universe.trajectory = [object(), object(), object()]
    universe.atoms = _MockAtomCollection(atoms)
    universe.select_atoms.side_effect = lambda selection, updating: {
        "chainid A": _MockAtomGroup([0, 2]),
        "chainid C": _MockAtomGroup([1]),
    }[selection]

    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=HydrogenBondSettings(),
    )

    analysis = HydrogenBondsAnalysis()
    mock_modules = _make_mdanalysis_module(MockHydrogenBondAnalysis)

    with (
        patch.dict(sys.modules, mock_modules),
        patch("polyzymd.analyses.hydrogen_bonds.TrajectoryLoader") as mock_loader_cls,
        patch("polyzymd.analyses.hydrogen_bonds._mda.compute_config_hash", return_value="abc123"),
    ):
        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader
        mock_loader.load_universe.return_value = universe
        mock_loader.get_timestep.return_value = 10.0

        result = _as_hydrogen_bond_result(analysis.run_replicate(ctx, 1))

    summary = result.summaries[0]
    assert summary.counts_per_frame == [1, 1, 0]
    assert summary.mean_hbonds_per_frame == pytest.approx(2.0 / 3.0)
    assert summary.fraction_frames_with_any_hbond == pytest.approx(2.0 / 3.0)


def test_compute_within_mode_classification(tmp_path: Path) -> None:
    """Within-mode should only count events where both atoms are in one group."""

    class MockHydrogenBondAnalysis:
        def __init__(self, **kwargs) -> None:
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            self.results.hbonds = np.array(
                [
                    [0, 0, 10, 1, 2.8, 160.0],
                    [1, 0, 10, 2, 2.9, 158.0],
                    [2, 2, 10, 3, 2.9, 158.0],
                ],
                dtype=float,
            )

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "A", 11, "THR", 1),
        2: _MockAtom(2, "C", 100, "OEG", 2),
        3: _MockAtom(3, "C", 101, "OEG", 3),
    }
    universe = MagicMock()
    universe.trajectory = [object(), object(), object()]
    universe.atoms = _MockAtomCollection(atoms)
    universe.select_atoms.side_effect = lambda selection, updating: {
        "chainid A": _MockAtomGroup([0, 1]),
        "chainid C": _MockAtomGroup([2, 3]),
    }[selection]

    settings = HydrogenBondSettings(
        groups={"protein": "chainid A", "polymer": "chainid C"},
        summaries=[HydrogenBondSummarySettings(name="protein_internal", within="protein")],
    )

    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=settings,
    )

    analysis = HydrogenBondsAnalysis()
    mock_modules = _make_mdanalysis_module(MockHydrogenBondAnalysis)

    with (
        patch.dict(sys.modules, mock_modules),
        patch("polyzymd.analyses.hydrogen_bonds.TrajectoryLoader") as mock_loader_cls,
        patch("polyzymd.analyses.hydrogen_bonds._mda.compute_config_hash", return_value="abc123"),
    ):
        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader
        mock_loader.load_universe.return_value = universe
        mock_loader.get_timestep.return_value = 10.0

        result = _as_hydrogen_bond_result(analysis.run_replicate(ctx, 1))

    summary = result.summaries[0]
    assert summary.counts_per_frame == [1, 0, 0]
    assert summary.mean_hbonds_per_frame == pytest.approx(1.0 / 3.0)
    assert summary.fraction_frames_with_any_hbond == pytest.approx(1.0 / 3.0)


def test_compute_multiple_summaries(tmp_path: Path) -> None:
    """Overlapping groups should allow one event to populate multiple summaries."""

    class MockHydrogenBondAnalysis:
        def __init__(self, **kwargs) -> None:
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            self.results.hbonds = np.array([[0, 0, 10, 1, 2.8, 160.0]], dtype=float)

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "C", 100, "OEG", 1),
    }
    universe = MagicMock()
    universe.trajectory = [object(), object(), object()]
    universe.atoms = _MockAtomCollection(atoms)
    universe.select_atoms.side_effect = lambda selection, updating: {
        "chainid A": _MockAtomGroup([0]),
        "chainid C": _MockAtomGroup([1]),
        "chainid A or chainid C": _MockAtomGroup([0, 1]),
    }[selection]

    settings = HydrogenBondSettings(
        groups={
            "protein": "chainid A",
            "polymer": "chainid C",
            "all_atoms": "chainid A or chainid C",
        },
        allow_overlapping_composition=True,
        summaries=[
            HydrogenBondSummarySettings(name="protein_polymer", between=("protein", "polymer")),
            HydrogenBondSummarySettings(name="within_all", within="all_atoms"),
        ],
    )

    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=settings,
    )

    analysis = HydrogenBondsAnalysis()
    mock_modules = _make_mdanalysis_module(MockHydrogenBondAnalysis)

    with (
        patch.dict(sys.modules, mock_modules),
        patch("polyzymd.analyses.hydrogen_bonds.TrajectoryLoader") as mock_loader_cls,
        patch("polyzymd.analyses.hydrogen_bonds._mda.compute_config_hash", return_value="abc123"),
    ):
        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader
        mock_loader.load_universe.return_value = universe
        mock_loader.get_timestep.return_value = 10.0

        result = _as_hydrogen_bond_result(analysis.run_replicate(ctx, 1))

    by_name = {summary.name: summary for summary in result.summaries}
    assert by_name["protein_polymer"].counts_per_frame == [1, 0, 0]
    assert by_name["within_all"].counts_per_frame == [1, 0, 0]


def test_aggregate_undirected_pair_merging(tmp_path: Path) -> None:
    """Undirected pair keys should merge residue order reversals across replicates."""
    analysis = HydrogenBondsAnalysis()
    ctx = _make_aggregate_context(tmp_path, replicates=(1, 2))

    rep1 = _make_replicate_result(
        1,
        2.0,
        0.5,
        undirected_pairs=[_make_undirected_pair(("A", 10, "SER"), ("C", 100, "OEG"), 0.4, 0.4)],
    )
    rep2 = _make_replicate_result(
        2,
        2.0,
        0.5,
        undirected_pairs=[_make_undirected_pair(("C", 100, "OEG"), ("A", 10, "SER"), 0.6, 0.6)],
    )

    aggregated = _aggregate_replicate_results(analysis, ctx, [rep1, rep2])
    summary = aggregated.summaries[0]
    assert len(summary.undirected_pairs) == 1
    pair = summary.undirected_pairs[0]
    assert pair.per_replicate_occupancy == pytest.approx([0.4, 0.6])
    assert pair.mean_occupancy == pytest.approx(0.5)


def test_aggregate_all_summaries_present(tmp_path: Path) -> None:
    """Aggregation should preserve values when all summaries are present."""
    analysis = HydrogenBondsAnalysis()

    settings = HydrogenBondSettings(
        groups={"protein": "chainid A", "polymer": "chainid C"},
        summaries=[
            HydrogenBondSummarySettings(name="protein_polymer", between=("protein", "polymer")),
            HydrogenBondSummarySettings(name="protein_internal", within="protein"),
        ],
    )
    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1, 2),
        sim_config=MagicMock(),
    )
    ctx = AggregateContext(
        condition=condition,
        replicates=(1, 2),
        output_dir=tmp_path / "aggregated",
        equilibration="10ns",
        settings=settings,
    )

    rep1 = HydrogenBondResult(
        replicate=1,
        settings_fingerprint=_settings_hash(settings),
        summaries=[
            HydrogenBondReplicateSummary(
                name="protein_polymer",
                mode="between",
                group_names=["protein", "polymer"],
                n_frames_used=10,
                mean_hbonds_per_frame=2.0,
                fraction_frames_with_any_hbond=0.5,
                counts_per_frame=[0] * 10,
            ),
            HydrogenBondReplicateSummary(
                name="protein_internal",
                mode="within",
                group_names=["protein"],
                n_frames_used=10,
                mean_hbonds_per_frame=1.0,
                fraction_frames_with_any_hbond=0.4,
                counts_per_frame=[0] * 10,
            ),
        ],
    )
    rep2 = HydrogenBondResult(
        replicate=2,
        settings_fingerprint=_settings_hash(settings),
        summaries=[
            HydrogenBondReplicateSummary(
                name="protein_polymer",
                mode="between",
                group_names=["protein", "polymer"],
                n_frames_used=10,
                mean_hbonds_per_frame=4.0,
                fraction_frames_with_any_hbond=0.7,
                counts_per_frame=[0] * 10,
            ),
            HydrogenBondReplicateSummary(
                name="protein_internal",
                mode="within",
                group_names=["protein"],
                n_frames_used=10,
                mean_hbonds_per_frame=2.0,
                fraction_frames_with_any_hbond=0.6,
                counts_per_frame=[0] * 10,
            ),
        ],
    )

    aggregated = _aggregate_replicate_results(analysis, ctx, [rep1, rep2])
    by_name = {summary.name: summary for summary in aggregated.summaries}

    assert by_name["protein_polymer"].per_replicate_mean_hbonds == pytest.approx([2.0, 4.0])
    assert by_name["protein_internal"].per_replicate_mean_hbonds == pytest.approx([1.0, 2.0])


def test_extract_metrics_empty_summaries() -> None:
    """extract_metrics should return an empty dict when summaries are empty."""
    analysis = HydrogenBondsAnalysis()
    summary = HydrogenBondAggregatedResult(replicates=[1, 2], n_replicates=2, summaries=[])
    assert analysis.extract_metrics(summary) == {}


def test_extract_metrics_preserves_replicate_values() -> None:
    """extract_metrics should preserve per-replicate values in MetricValue."""
    analysis = HydrogenBondsAnalysis()
    summary = HydrogenBondAggregatedResult(
        replicates=[1, 2, 3],
        n_replicates=3,
        summaries=[_make_aggregated_summary("protein_polymer", 3.0, 0.2, [2.5, 3.0, 3.5])],
    )
    metric = analysis.extract_metrics(summary)["mean_hbonds_protein_polymer"]
    assert metric.replicate_values == pytest.approx([2.5, 3.0, 3.5])


def _make_aggregated_result_with_pair(
    mean_occupancy: float,
    sem_occupancy: float,
    per_replicate_occupancy: list[float],
) -> HydrogenBondAggregatedResult:
    """Build an aggregated result containing one top-pair candidate."""

    summary = _make_aggregated_summary(
        "protein_polymer",
        mean_occupancy,
        sem_occupancy,
        per_replicate_occupancy,
    )
    pair = UndirectedPairAggregate(
        residue_a=ResidueRef(chain_id="A", resid=10, resname="SER"),
        residue_b=ResidueRef(chain_id="C", resid=100, resname="OEG"),
        mean_occupancy=mean_occupancy,
        sem_occupancy=sem_occupancy,
        mean_events_per_frame=mean_occupancy,
        sem_events_per_frame=sem_occupancy,
        per_replicate_occupancy=per_replicate_occupancy,
    )
    return HydrogenBondAggregatedResult(
        replicates=list(range(1, len(per_replicate_occupancy) + 1)),
        n_replicates=len(per_replicate_occupancy),
        summaries=[summary.model_copy(update={"undirected_pairs": [pair]})],
    )


def test_plot_top_pairs_uses_sem_occupancy_for_xerr(tmp_path: Path) -> None:
    """Top-pairs plot should pass aggregate SEM values as horizontal xerr."""
    pytest.importorskip("matplotlib")
    import matplotlib.axes

    from polyzymd.analyses.hydrogen_bonds._plotters import plot_top_pairs

    cond_a_sem = 0.123
    cond_b_sem = 0.217
    # Sentinels intentionally differ from SEMs recomputed from replicate values
    results = {
        "CondA": _make_aggregated_result_with_pair(0.35, cond_a_sem, [0.31, 0.39]),
        "CondB": _make_aggregated_result_with_pair(0.55, cond_b_sem, [0.49, 0.61]),
    }
    captured_kwargs: list[dict[str, Any]] = []
    original_barh = matplotlib.axes.Axes.barh

    def _capture_barh(self: Any, *args: Any, **kwargs: Any) -> Any:
        captured_kwargs.append(kwargs.copy())
        return original_barh(self, *args, **kwargs)

    with patch("matplotlib.axes.Axes.barh", autospec=True, side_effect=_capture_barh):
        path = plot_top_pairs(
            results,
            ["CondA", "CondB"],
            "protein_polymer",
            tmp_path / "plots",
            PlotSettings(),
        )

    assert path is not None
    assert len(captured_kwargs) == 2
    assert captured_kwargs[0]["xerr"] == pytest.approx([cond_a_sem])
    assert captured_kwargs[1]["xerr"] == pytest.approx([cond_b_sem])
    assert all(kwargs["capsize"] == PlotSettings().theme.bar_capsize for kwargs in captured_kwargs)


def test_plot_top_pairs_suppresses_singleton_errorbars(tmp_path: Path) -> None:
    """Top-pairs plot should omit xerr when only singleton replicate data exist."""
    pytest.importorskip("matplotlib")
    import matplotlib.axes

    from polyzymd.analyses.hydrogen_bonds._plotters import plot_top_pairs

    results = {
        "CondA": _make_aggregated_result_with_pair(0.35, 0.04, [0.35]),
        "CondB": _make_aggregated_result_with_pair(0.55, 0.06, [0.55]),
    }
    captured_kwargs: list[dict[str, Any]] = []
    original_barh = matplotlib.axes.Axes.barh

    def _capture_barh(self: Any, *args: Any, **kwargs: Any) -> Any:
        captured_kwargs.append(kwargs.copy())
        return original_barh(self, *args, **kwargs)

    with patch("matplotlib.axes.Axes.barh", autospec=True, side_effect=_capture_barh):
        path = plot_top_pairs(
            results,
            ["CondA", "CondB"],
            "protein_polymer",
            tmp_path / "plots",
            PlotSettings(),
        )

    assert path is not None
    assert captured_kwargs
    assert all("xerr" not in kwargs for kwargs in captured_kwargs)
    assert all("capsize" not in kwargs for kwargs in captured_kwargs)


def test_plot_summary_comparison_smoke(tmp_path: Path) -> None:
    """plot_summary_comparison should create a figure file with minimal data."""
    pytest.importorskip("matplotlib")
    from polyzymd.analyses.hydrogen_bonds._plotters import plot_summary_comparison

    output_dir = tmp_path / "plots"
    output_dir.mkdir(parents=True)
    results = {
        "CondA": HydrogenBondAggregatedResult(
            replicates=[1],
            n_replicates=1,
            summaries=[_make_aggregated_summary("protein_polymer", 2.0, 0.1, [2.0])],
        )
    }

    path = plot_summary_comparison(results, ["CondA"], output_dir, PlotSettings())
    assert path is not None
    assert path.exists()


def test_plot_summary_comparison_facets_multiple_summaries(tmp_path: Path) -> None:
    """Summary comparison plot should create one subplot per summary."""
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    from polyzymd.analyses.hydrogen_bonds._plotters import plot_summary_comparison

    output_dir = tmp_path / "plots"
    output_dir.mkdir(parents=True)
    results = {
        "CondA": HydrogenBondAggregatedResult(
            replicates=[1],
            n_replicates=1,
            summaries=[
                _make_aggregated_summary("s1", 2.0, 0.1, [2.0]),
                _make_aggregated_summary("s2", 10.0, 0.2, [10.0]),
            ],
        ),
        "CondB": HydrogenBondAggregatedResult(
            replicates=[1],
            n_replicates=1,
            summaries=[
                _make_aggregated_summary("s1", 3.0, 0.1, [3.0]),
                _make_aggregated_summary("s2", 12.0, 0.2, [12.0]),
            ],
        ),
    }

    captured_nrows: list[int] = []
    original_subplots = plt.subplots

    def _capture_subplots(*args, **kwargs):
        captured_nrows.append(args[0] if args else kwargs.get("nrows", 1))
        return original_subplots(*args, **kwargs)

    with patch("matplotlib.pyplot.subplots", side_effect=_capture_subplots):
        path = plot_summary_comparison(results, ["CondA", "CondB"], output_dir, PlotSettings())

    assert path is not None
    assert path.exists()
    assert captured_nrows and captured_nrows[0] == 2


def test_plot_composition_absolute_smoke(tmp_path: Path) -> None:
    """plot_composition_absolute should create a figure file with minimal data."""
    pytest.importorskip("matplotlib")
    from polyzymd.analyses.hydrogen_bonds._plotters import plot_composition_absolute

    output_dir = tmp_path / "plots"
    output_dir.mkdir(parents=True)
    results = {
        "CondA": HydrogenBondAggregatedResult(
            replicates=[1],
            n_replicates=1,
            summaries=[_make_aggregated_summary("protein_polymer", 2.0, 0.1, [2.0])],
            composition_entries=[
                AggregatedCompositionEntry(
                    donor_partition="protein",
                    acceptor_partition="polymer",
                    mean_hbonds_per_frame=1.0,
                    sem_hbonds_per_frame=0.0,
                    per_replicate_hbonds=[1.0],
                    mean_fraction_of_total=1.0,
                    sem_fraction_of_total=0.0,
                    per_replicate_fraction=[1.0],
                )
            ],
        )
    }

    path = plot_composition_absolute(results, ["CondA"], output_dir, PlotSettings())
    assert path is not None
    assert path.exists()


def test_plot_composition_absolute_uses_replicate_specific_bases(tmp_path: Path) -> None:
    """Absolute composition dots should use per-replicate cumulative bases."""
    pytest.importorskip("matplotlib")
    from polyzymd.analyses.hydrogen_bonds._plotters import plot_composition_absolute

    output_dir = tmp_path / "plots"
    output_dir.mkdir(parents=True)
    results = {
        "CondA": HydrogenBondAggregatedResult(
            replicates=[1, 2],
            n_replicates=2,
            summaries=[_make_aggregated_summary("protein_polymer", 2.0, 0.1, [1.0, 3.0])],
            composition_entries=[
                AggregatedCompositionEntry(
                    donor_partition="groupA",
                    acceptor_partition="groupA",
                    mean_hbonds_per_frame=2.0,
                    sem_hbonds_per_frame=1.0,
                    per_replicate_hbonds=[1.0, 3.0],
                    mean_fraction_of_total=0.4,
                    sem_fraction_of_total=0.0,
                    per_replicate_fraction=[0.2, 0.6],
                ),
                AggregatedCompositionEntry(
                    donor_partition="groupA",
                    acceptor_partition="groupB",
                    mean_hbonds_per_frame=3.0,
                    sem_hbonds_per_frame=1.0,
                    per_replicate_hbonds=[4.0, 2.0],
                    mean_fraction_of_total=0.6,
                    sem_fraction_of_total=0.0,
                    per_replicate_fraction=[0.8, 0.4],
                ),
            ],
        )
    }

    with patch(
        "polyzymd.analyses.hydrogen_bonds._plotters.scatter_stacked_segment_replicates",
        autospec=True,
    ) as scatter:
        path = plot_composition_absolute(results, ["CondA"], output_dir, PlotSettings())

    assert path is not None
    assert scatter.call_count == 2
    assert scatter.call_args_list[0].kwargs["replicate_base_values"] == [0.0, 0.0]
    assert scatter.call_args_list[0].kwargs["placement"] == "end"
    assert scatter.call_args_list[1].kwargs["replicate_base_values"] == [1.0, 3.0]
    assert scatter.call_args_list[1].kwargs["placement"] == "end"


def test_plot_composition_fraction_overlap_exceeds_one(tmp_path: Path) -> None:
    """Composition fraction plot should not clip stacked fractions above 1.0."""
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    from polyzymd.analyses.hydrogen_bonds._plotters import plot_composition_fraction

    output_dir = tmp_path / "plots"
    output_dir.mkdir(parents=True)
    # Two partitions with fractions summing to 1.4, simulating overlap
    results = {
        "CondA": HydrogenBondAggregatedResult(
            replicates=[1],
            n_replicates=1,
            summaries=[_make_aggregated_summary("s", 5.0, 0.1, [5.0])],
            composition_entries=[
                AggregatedCompositionEntry(
                    donor_partition="groupA",
                    acceptor_partition="groupA",
                    mean_hbonds_per_frame=3.5,
                    sem_hbonds_per_frame=0.0,
                    per_replicate_hbonds=[3.5],
                    mean_fraction_of_total=0.7,
                    sem_fraction_of_total=0.0,
                    per_replicate_fraction=[0.7],
                ),
                AggregatedCompositionEntry(
                    donor_partition="groupA",
                    acceptor_partition="groupB",
                    mean_hbonds_per_frame=3.5,
                    sem_hbonds_per_frame=0.0,
                    per_replicate_hbonds=[3.5],
                    mean_fraction_of_total=0.7,
                    sem_fraction_of_total=0.0,
                    per_replicate_fraction=[0.7],
                ),
            ],
        )
    }
    captured_upper_limits: list[float] = []

    def _capture_save_figure(fig: Any, output_path: Path, plot_settings: PlotSettings) -> Path:
        """Capture the y-axis limit before closing the test figure."""
        del plot_settings
        captured_upper_limits.append(float(fig.axes[0].get_ylim()[1]))
        plt.close(fig)
        return output_path

    with patch(
        "polyzymd.analyses.hydrogen_bonds._plotters.save_figure",
        autospec=True,
        side_effect=_capture_save_figure,
    ):
        path = plot_composition_fraction(results, ["CondA"], output_dir, PlotSettings())

    assert path == output_dir / "hbond_composition_fraction.png"
    assert captured_upper_limits
    assert captured_upper_limits[0] > 1.4


def test_plot_composition_fraction_uses_replicate_specific_bases(tmp_path: Path) -> None:
    """Fraction composition dots should use per-replicate cumulative bases."""
    pytest.importorskip("matplotlib")
    from polyzymd.analyses.hydrogen_bonds._plotters import plot_composition_fraction

    output_dir = tmp_path / "plots"
    output_dir.mkdir(parents=True)
    results = {
        "CondA": HydrogenBondAggregatedResult(
            replicates=[1, 2],
            n_replicates=2,
            summaries=[_make_aggregated_summary("protein_polymer", 2.0, 0.1, [2.0, 2.0])],
            composition_entries=[
                AggregatedCompositionEntry(
                    donor_partition="groupA",
                    acceptor_partition="groupA",
                    mean_hbonds_per_frame=2.0,
                    sem_hbonds_per_frame=0.0,
                    per_replicate_hbonds=[2.0, 2.0],
                    mean_fraction_of_total=0.5,
                    sem_fraction_of_total=0.1,
                    per_replicate_fraction=[0.4, 0.6],
                ),
                AggregatedCompositionEntry(
                    donor_partition="groupA",
                    acceptor_partition="groupB",
                    mean_hbonds_per_frame=2.0,
                    sem_hbonds_per_frame=0.0,
                    per_replicate_hbonds=[2.0, 2.0],
                    mean_fraction_of_total=0.5,
                    sem_fraction_of_total=0.1,
                    per_replicate_fraction=[0.7, 0.3],
                ),
            ],
        )
    }

    with patch(
        "polyzymd.analyses.hydrogen_bonds._plotters.scatter_stacked_segment_replicates",
        autospec=True,
    ) as scatter:
        path = plot_composition_fraction(results, ["CondA"], output_dir, PlotSettings())

    assert path is not None
    assert scatter.call_count == 2
    assert scatter.call_args_list[0].kwargs["replicate_base_values"] == [0.0, 0.0]
    assert scatter.call_args_list[0].kwargs["placement"] == "end"
    assert scatter.call_args_list[1].kwargs["replicate_base_values"] == [0.4, 0.6]
    assert scatter.call_args_list[1].kwargs["placement"] == "end"


def test_plot_composition_fraction_ylim_covers_skewed_replicate_endpoints(
    tmp_path: Path,
) -> None:
    """Fraction composition y-axis should include end-placed replicate endpoints."""
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    from polyzymd.analyses.hydrogen_bonds._plotters import plot_composition_fraction

    output_dir = tmp_path / "plots"
    output_dir.mkdir(parents=True)
    results = {
        "CondA": HydrogenBondAggregatedResult(
            replicates=[1, 2],
            n_replicates=2,
            summaries=[_make_aggregated_summary("protein_polymer", 2.0, 0.1, [2.0, 2.0])],
            composition_entries=[
                AggregatedCompositionEntry(
                    donor_partition="groupA",
                    acceptor_partition="groupA",
                    mean_hbonds_per_frame=2.0,
                    sem_hbonds_per_frame=0.0,
                    per_replicate_hbonds=[2.0, 2.0],
                    mean_fraction_of_total=0.5,
                    sem_fraction_of_total=0.0,
                    per_replicate_fraction=[0.9, 0.1],
                ),
                AggregatedCompositionEntry(
                    donor_partition="groupA",
                    acceptor_partition="groupB",
                    mean_hbonds_per_frame=2.0,
                    sem_hbonds_per_frame=0.0,
                    per_replicate_hbonds=[2.0, 2.0],
                    mean_fraction_of_total=0.5,
                    sem_fraction_of_total=0.0,
                    per_replicate_fraction=[0.9, 0.1],
                ),
            ],
        )
    }
    captured_upper_limits: list[float] = []

    def _capture_save_figure(fig: Any, output_path: Path, plot_settings: PlotSettings) -> Path:
        """Capture the y-axis limit before closing the test figure."""
        del plot_settings
        captured_upper_limits.append(float(fig.axes[0].get_ylim()[1]))
        plt.close(fig)
        return output_path

    with patch(
        "polyzymd.analyses.hydrogen_bonds._plotters.save_figure",
        autospec=True,
        side_effect=_capture_save_figure,
    ):
        path = plot_composition_fraction(results, ["CondA"], output_dir, PlotSettings())

    assert path == output_dir / "hbond_composition_fraction.png"
    assert captured_upper_limits
    assert captured_upper_limits[0] > 1.8


@pytest.mark.parametrize(
    "theme",
    [PlotTheme(dot_size=0, dot_alpha=0.7), PlotTheme(dot_size=8, dot_alpha=0.0)],
)
def test_plot_composition_fraction_ylim_ignores_hidden_replicate_endpoints(
    tmp_path: Path,
    theme: PlotTheme,
) -> None:
    """Fraction composition y-axis should ignore hidden replicate endpoints."""
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    from polyzymd.analyses.hydrogen_bonds._plotters import plot_composition_fraction

    output_dir = tmp_path / "plots"
    output_dir.mkdir(parents=True)
    results = {
        "CondA": HydrogenBondAggregatedResult(
            replicates=[1, 2],
            n_replicates=2,
            summaries=[_make_aggregated_summary("protein_polymer", 2.0, 0.1, [2.0, 2.0])],
            composition_entries=[
                AggregatedCompositionEntry(
                    donor_partition="groupA",
                    acceptor_partition="groupA",
                    mean_hbonds_per_frame=2.0,
                    sem_hbonds_per_frame=0.0,
                    per_replicate_hbonds=[2.0, 2.0],
                    mean_fraction_of_total=0.5,
                    sem_fraction_of_total=0.0,
                    per_replicate_fraction=[0.9, 0.1],
                ),
                AggregatedCompositionEntry(
                    donor_partition="groupA",
                    acceptor_partition="groupB",
                    mean_hbonds_per_frame=2.0,
                    sem_hbonds_per_frame=0.0,
                    per_replicate_hbonds=[2.0, 2.0],
                    mean_fraction_of_total=0.5,
                    sem_fraction_of_total=0.0,
                    per_replicate_fraction=[0.9, 0.1],
                ),
            ],
        )
    }
    captured_upper_limits: list[float] = []

    def _capture_save_figure(fig: Any, output_path: Path, plot_settings: PlotSettings) -> Path:
        """Capture the y-axis limit before closing the test figure."""
        del plot_settings
        captured_upper_limits.append(float(fig.axes[0].get_ylim()[1]))
        plt.close(fig)
        return output_path

    with patch(
        "polyzymd.analyses.hydrogen_bonds._plotters.save_figure",
        autospec=True,
        side_effect=_capture_save_figure,
    ):
        path = plot_composition_fraction(
            results,
            ["CondA"],
            output_dir,
            PlotSettings(theme=theme),
        )

    assert path == output_dir / "hbond_composition_fraction.png"
    assert captured_upper_limits
    assert captured_upper_limits[0] < 1.8


def test_plot_timeseries_smoke(tmp_path: Path) -> None:
    """plot_timeseries should create a figure file with minimal traces."""
    pytest.importorskip("matplotlib")
    from polyzymd.analyses.hydrogen_bonds._plotters import plot_timeseries

    output_dir = tmp_path / "plots"
    output_dir.mkdir(parents=True)
    results = {
        "CondA": HydrogenBondAggregatedResult(
            replicates=[1],
            n_replicates=1,
            summaries=[_make_aggregated_summary("protein_polymer", 2.0, 0.1, [2.0])],
        )
    }
    replicate_data = {"CondA": {"protein_polymer": [[1, 2, 3, 2]]}}

    path = plot_timeseries(
        results,
        replicate_data,
        ["CondA"],
        "protein_polymer",
        output_dir,
        PlotSettings(),
    )
    assert path is not None
    assert path.exists()


def test_condition_plots_use_semantic_colors_but_composition_remains_categorical(
    tmp_path: Path,
) -> None:
    """H-bond condition plots should not reuse semantic colors for composition partitions."""
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba

    from polyzymd.analyses.hydrogen_bonds import _plotters

    output_dir = tmp_path / "plots"
    output_dir.mkdir(parents=True)
    results = {
        "Treatment": HydrogenBondAggregatedResult(
            replicates=[1],
            n_replicates=1,
            summaries=[_make_aggregated_summary("protein_polymer", 4.0, 0.1, [4.0])],
            composition_entries=[
                AggregatedCompositionEntry(
                    donor_partition="protein",
                    acceptor_partition="polymer",
                    mean_hbonds_per_frame=2.0,
                    sem_hbonds_per_frame=0.0,
                    per_replicate_hbonds=[2.0],
                    mean_fraction_of_total=1.0,
                    sem_fraction_of_total=0.0,
                    per_replicate_fraction=[1.0],
                )
            ],
        ),
        "Control": HydrogenBondAggregatedResult(
            replicates=[1],
            n_replicates=1,
            summaries=[_make_aggregated_summary("protein_polymer", 2.0, 0.1, [2.0])],
            composition_entries=[
                AggregatedCompositionEntry(
                    donor_partition="protein",
                    acceptor_partition="polymer",
                    mean_hbonds_per_frame=1.0,
                    sem_hbonds_per_frame=0.0,
                    per_replicate_hbonds=[1.0],
                    mean_fraction_of_total=1.0,
                    sem_fraction_of_total=0.0,
                    per_replicate_fraction=[1.0],
                )
            ],
        ),
    }
    plot_settings = PlotSettings(
        semantic_colors={
            "enabled": True,
            "order": ["Control", "Treatment"],
            "conditions": {
                "Control": {"role": "control"},
                "Treatment": {"color": "#ff7f0e"},
            },
            "control_color": "#111111",
        }
    )
    captured = []

    def _capture_save_figure(fig, output_path, settings):
        captured.append(fig)
        return output_path

    with patch.object(_plotters, "save_figure", side_effect=_capture_save_figure):
        summary_path = _plotters.plot_summary_comparison(
            results,
            ["Treatment", "Control"],
            output_dir,
            plot_settings,
            control_label="Control",
        )
        composition_path = _plotters.plot_composition_absolute(
            results,
            ["Treatment", "Control"],
            output_dir,
            plot_settings,
        )

    assert summary_path == output_dir / "hbond_summary_comparison.png"
    assert composition_path == output_dir / "hbond_composition_absolute.png"
    summary_ax = captured[0].axes[0]
    composition_ax = captured[1].axes[0]
    assert [tick.get_text() for tick in summary_ax.get_xticklabels()] == [
        "Control",
        "Treatment",
    ]
    assert to_rgba(summary_ax.patches[0].get_facecolor())[:3] == to_rgba("#111111")[:3]
    assert to_rgba(summary_ax.patches[1].get_facecolor())[:3] == to_rgba("#ff7f0e")[:3]
    assert [tick.get_text() for tick in composition_ax.get_xticklabels()] == [
        "Control",
        "Treatment",
    ]
    assert to_rgba(composition_ax.patches[0].get_facecolor())[:3] != to_rgba("#111111")[:3]
    for fig in captured:
        plt.close(fig)


def test_rank_conditions_neutral_sorted_by_value_descending() -> None:
    """Neutral rankings should be ordered by descending metric values."""
    metrics = {
        "Alpha": get_analysis("hydrogen_bonds")().extract_metrics(
            HydrogenBondAggregatedResult(
                replicates=[1],
                n_replicates=1,
                summaries=[_make_aggregated_summary("s", 1.0, 0.1, [1.0])],
            )
        )["mean_hbonds_s"],
        "Beta": get_analysis("hydrogen_bonds")().extract_metrics(
            HydrogenBondAggregatedResult(
                replicates=[1],
                n_replicates=1,
                summaries=[_make_aggregated_summary("s", 3.0, 0.1, [3.0])],
            )
        )["mean_hbonds_s"],
        "Gamma": get_analysis("hydrogen_bonds")().extract_metrics(
            HydrogenBondAggregatedResult(
                replicates=[1],
                n_replicates=1,
                summaries=[_make_aggregated_summary("s", 2.0, 0.1, [2.0])],
            )
        )["mean_hbonds_s"],
    }
    assert rank_conditions(metrics) == ["Beta", "Gamma", "Alpha"]


def test_format_pct_uses_semantic_infinity_labels() -> None:
    """Percent formatter should use semantic labels for infinite values."""
    assert format_pct(np.inf) == "new (baseline=0)"
    assert format_pct(-np.inf) == "gone (current=0)"
    assert format_pct(np.nan) == "undefined"


def test_overlap_composition_raises_by_default() -> None:
    """Overlapping composition partitions should raise by default."""
    analysis = HydrogenBondsAnalysis()
    composition_settings = HydrogenBondCompositionSettings(
        partitions={"protein": "chainid A", "polymer": "chainid C"}
    )
    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "C", 100, "OEG", 1),
    }
    universe = MagicMock()
    universe.atoms = _MockAtomCollection(atoms)
    universe.select_atoms.side_effect = lambda selection: {
        "chainid A": _MockAtomGroup([0, 1]),
        "chainid C": _MockAtomGroup([1]),
    }[selection]

    with pytest.raises(ValueError, match="allow_overlapping_composition: true"):
        _ = analysis._compute_composition(
            composition_settings=composition_settings,
            hbond_array=np.empty((0, 6), dtype=float),
            universe=universe,
            start_frame=0,
            n_frames=3,
        )


def test_unique_pairs_metric_and_aggregation(tmp_path: Path) -> None:
    """Replicate and aggregate results should include unique-pairs statistics."""
    analysis = HydrogenBondsAnalysis()
    results = [
        _make_replicate_result(1, 2.0, 0.4),
        _make_replicate_result(2, 4.0, 0.6),
    ]
    results[0].summaries[0].mean_unique_pairs_per_frame = 1.0
    results[1].summaries[0].mean_unique_pairs_per_frame = 3.0
    ctx = _make_aggregate_context(tmp_path, replicates=(1, 2))

    aggregated = _aggregate_replicate_results(analysis, ctx, results)
    summary = aggregated.summaries[0]
    assert hasattr(summary, "mean_unique_pairs_per_frame")
    assert summary.mean_unique_pairs_per_frame == pytest.approx(2.0)
    assert summary.sem_unique_pairs_per_frame > 0.0


def test_replicate_does_not_truncate_pairs_before_aggregation(tmp_path: Path) -> None:
    """run_replicate should keep all pairs and defer top_n truncation."""

    class MockHydrogenBondAnalysis:
        def __init__(self, **kwargs) -> None:
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            self.results.hbonds = np.array(
                [
                    [0, 0, 10, 1, 2.8, 160.0],
                    [0, 2, 10, 3, 2.8, 160.0],
                ],
                dtype=float,
            )

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "C", 100, "OEG", 1),
        2: _MockAtom(2, "A", 11, "THR", 2),
        3: _MockAtom(3, "C", 101, "OEG", 3),
    }
    universe = MagicMock()
    universe.trajectory = [object(), object()]
    universe.atoms = _MockAtomCollection(atoms)
    universe.select_atoms.side_effect = lambda selection, updating: {
        "chainid A": _MockAtomGroup([0, 2]),
        "chainid C": _MockAtomGroup([1, 3]),
    }[selection]

    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    settings = HydrogenBondSettings(top_n_pairs=1)
    ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=settings,
    )

    analysis = HydrogenBondsAnalysis()
    mock_modules = _make_mdanalysis_module(MockHydrogenBondAnalysis)

    with (
        patch.dict(sys.modules, mock_modules),
        patch("polyzymd.analyses.hydrogen_bonds.TrajectoryLoader") as mock_loader_cls,
        patch("polyzymd.analyses.hydrogen_bonds._mda.compute_config_hash", return_value="abc123"),
    ):
        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader
        mock_loader.load_universe.return_value = universe
        mock_loader.get_timestep.return_value = 10.0

        result = _as_hydrogen_bond_result(analysis.run_replicate(ctx, 1))

    assert len(result.summaries[0].directed_residue_pairs) == 2
    assert len(result.summaries[0].undirected_residue_pairs) == 2


def test_full_lifecycle_mocked(tmp_path: Path) -> None:
    """Mocked pipeline should flow compute to aggregate to metrics and formatting."""

    class MockHydrogenBondAnalysis:
        def __init__(self, universe: MagicMock, **kwargs) -> None:
            self._universe = universe
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            self.results.hbonds = self._universe.hbonds_data

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "C", 100, "OEG", 1),
    }
    universe_1 = MagicMock()
    universe_1.trajectory = [object(), object(), object(), object()]
    universe_1.atoms = _MockAtomCollection(atoms)
    universe_1.hbonds_data = np.array(
        [
            [0, 0, 10, 1, 2.8, 160.0],
            [1, 0, 10, 1, 2.9, 158.0],
        ],
        dtype=float,
    )
    universe_1.select_atoms.side_effect = lambda selection, updating: {
        "chainid A": _MockAtomGroup([0]),
        "chainid C": _MockAtomGroup([1]),
    }[selection]

    universe_2 = MagicMock()
    universe_2.trajectory = [object(), object(), object(), object()]
    universe_2.atoms = _MockAtomCollection(atoms)
    universe_2.hbonds_data = np.array(
        [
            [0, 0, 10, 1, 2.8, 160.0],
            [1, 0, 10, 1, 2.9, 158.0],
            [2, 0, 10, 1, 2.9, 158.0],
        ],
        dtype=float,
    )
    universe_2.select_atoms.side_effect = lambda selection, updating: {
        "chainid A": _MockAtomGroup([0]),
        "chainid C": _MockAtomGroup([1]),
    }[selection]

    analysis = HydrogenBondsAnalysis()
    condition = Condition(
        label="Lifecycle",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1, 2),
        sim_config=MagicMock(),
    )
    settings = HydrogenBondSettings()

    rep_ctx_1 = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=settings,
    )
    rep_ctx_2 = ReplicateContext(
        condition=condition,
        replicate=2,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_2",
        equilibration="0ns",
        recompute=True,
        settings=settings,
    )

    mock_modules = _make_mdanalysis_module(MockHydrogenBondAnalysis)
    with (
        patch.dict(sys.modules, mock_modules),
        patch("polyzymd.analyses.hydrogen_bonds.TrajectoryLoader") as mock_loader_cls,
        patch("polyzymd.analyses.hydrogen_bonds._mda.compute_config_hash", return_value="abc123"),
    ):
        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader
        mock_loader.load_universe.side_effect = [universe_1, universe_2]
        mock_loader.get_timestep.return_value = 10.0

        rep_result_1 = analysis.run_replicate(rep_ctx_1, 1)
        rep_result_2 = analysis.run_replicate(rep_ctx_2, 2)

    aggregate_ctx = AggregateContext(
        condition=condition,
        replicates=(1, 2),
        output_dir=tmp_path / "aggregated",
        equilibration="0ns",
        settings=settings,
    )
    aggregated = analysis.aggregate(aggregate_ctx, [rep_result_1, rep_result_2])
    legacy_aggregated = _as_hydrogen_bond_aggregate(aggregated)
    legacy_rep_result_1 = _as_hydrogen_bond_result(rep_result_1)
    legacy_rep_result_2 = _as_hydrogen_bond_result(rep_result_2)
    metrics = analysis.extract_metrics(aggregated)

    metric_key = "mean_hbonds_protein_polymer"
    metric = metrics[metric_key]
    comparison = ComparisonResult(
        analysis_type="hydrogen_bonds",
        name="hbonds_lifecycle",
        conditions=[
            ConditionSummary(
                label="Lifecycle",
                n_replicates=2,
                **{f"{metric_key}_mean": metric.mean, f"{metric_key}_sem": metric.sem},
            )
        ],
        pairwise_comparisons=[],
        ranking=["Lifecycle"],
        rankings_by_metric={metric_key: ["Lifecycle"]},
        equilibration_time="0ns",
        created_at="2026-01-01T00:00:00",
        polyzymd_version="test",
    )
    formatted = analysis.format(comparison, output_format="text")

    assert isinstance(rep_result_1, ReplicateArtifact)
    assert isinstance(rep_result_2, ReplicateArtifact)
    assert isinstance(aggregated, ConditionArtifact)
    assert isinstance(legacy_aggregated, HydrogenBondAggregatedResult)
    assert metric_key in metrics
    assert metric.replicate_values == pytest.approx(
        [
            legacy_rep_result_1.summaries[0].mean_hbonds_per_frame,
            legacy_rep_result_2.summaries[0].mean_hbonds_per_frame,
        ]
    )
    assert isinstance(formatted, str)
    assert "Hydrogen Bond Analysis" in formatted


def test_run_replicate_hydrogens_sel_explicit(tmp_path: Path) -> None:
    """run_replicate should pass explicit element H selection, not None."""

    captured_kwargs: list[dict] = []

    class MockHydrogenBondAnalysis:
        def __init__(self, **kwargs) -> None:
            self.kwargs = kwargs
            self.results = types.SimpleNamespace(hbonds=np.empty((0, 6), dtype=float))
            captured_kwargs.append(kwargs)

        def run(self, start: int, stop: int | None, step: int, verbose: bool) -> None:
            pass

    atoms = {
        0: _MockAtom(0, "A", 10, "SER", 0),
        1: _MockAtom(1, "C", 100, "OEG", 1),
    }
    universe = MagicMock()
    universe.trajectory = [object(), object(), object()]
    universe.atoms = _MockAtomCollection(atoms)
    selections = {
        "chainid A": _MockAtomGroup([0]),
        "chainid C": _MockAtomGroup([1]),
    }
    universe.select_atoms.side_effect = lambda selection, updating: selections[selection]

    condition = Condition(
        label="test",
        config_path=Path("/tmp/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    settings = HydrogenBondSettings()
    ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=settings,
    )

    analysis = HydrogenBondsAnalysis()
    mock_modules = _make_mdanalysis_module(MockHydrogenBondAnalysis)

    with (
        patch.dict(sys.modules, mock_modules),
        patch("polyzymd.analyses.hydrogen_bonds.TrajectoryLoader") as mock_loader_cls,
        patch("polyzymd.analyses.hydrogen_bonds._mda.compute_config_hash", return_value="abc123"),
    ):
        mock_loader = MagicMock()
        mock_loader_cls.return_value = mock_loader
        mock_loader.load_universe.return_value = universe
        mock_loader.get_timestep.return_value = 10.0

        analysis.run_replicate(ctx, 1)

    assert len(captured_kwargs) == 1
    hydrogens_sel = captured_kwargs[0]["hydrogens_sel"]
    assert hydrogens_sel is not None, "hydrogens_sel must not be None (causes NoDataError on PDB)"
    assert "element H" in hydrogens_sel
    expected = "((chainid A) or (chainid C)) and element H"
    assert hydrogens_sel == expected


class TestLoadReplicateResult:
    """Tests for the custom _load_replicate_result override."""

    def test_loads_canonical_result_json(self, tmp_path: Path) -> None:
        """Should prefer canonical result.json when it exists."""
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()

        rep_result = HydrogenBondResult(
            config_hash="abc123",
            replicate=1,
            equilibration_time=0.0,
            equilibration_unit="ns",
            selection_string="chainid A",
            summaries=[],
            composition_entries=[],
        )
        result_path = run_dir / "result.json"
        result_path.write_text(rep_result.model_dump_json())

        analysis = HydrogenBondsAnalysis()
        loaded = analysis._load_replicate_result(run_dir)
        assert loaded is not None

    def test_ignores_custom_cache(self, tmp_path: Path) -> None:
        """Should ignore legacy hbonds_eq*.json files when result.json is absent."""
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()

        rep_result = HydrogenBondResult(
            config_hash="abc123",
            replicate=1,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="chainid A",
            summaries=[],
            composition_entries=[],
        )
        cache_path = run_dir / "hbonds_eq10ns_deadbeef.json"
        cache_path.write_text(rep_result.model_dump_json())

        analysis = HydrogenBondsAnalysis()
        loaded = analysis._load_replicate_result(run_dir)
        assert loaded is None

    def test_returns_none_for_empty_directory(self, tmp_path: Path) -> None:
        """Should return None when no result files exist."""
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()

        analysis = HydrogenBondsAnalysis()
        loaded = analysis._load_replicate_result(run_dir)
        assert loaded is None

    def test_returns_none_for_missing_directory(self, tmp_path: Path) -> None:
        """Should return None when directory doesn't exist."""
        run_dir = tmp_path / "run_nonexistent"

        analysis = HydrogenBondsAnalysis()
        loaded = analysis._load_replicate_result(run_dir)
        assert loaded is None

    def test_ignores_multiple_legacy_cache_files(self, tmp_path: Path) -> None:
        """Should ignore multiple legacy custom cache files."""
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()

        rep_result = HydrogenBondResult(
            config_hash="abc123",
            replicate=1,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="chainid A",
            summaries=[],
            composition_entries=[],
        )
        json_content = rep_result.model_dump_json()

        (run_dir / "hbonds_eq10ns_aaaa1111.json").write_text(json_content)
        (run_dir / "hbonds_eq10ns_bbbb2222.json").write_text(json_content)

        analysis = HydrogenBondsAnalysis()
        loaded = analysis._load_replicate_result(run_dir)

        assert loaded is None

    def test_ignores_legacy_caches_across_equilibration_keys(self, tmp_path: Path) -> None:
        """Legacy caches with multiple equilibration keys should be ignored."""
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()

        eq10_result = HydrogenBondResult(
            config_hash="eq10_hash",
            replicate=1,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="chainid A",
            summaries=[],
            composition_entries=[],
        )
        eq20_result = HydrogenBondResult(
            config_hash="eq20_hash",
            replicate=1,
            equilibration_time=20.0,
            equilibration_unit="ns",
            selection_string="chainid A",
            summaries=[],
            composition_entries=[],
        )

        (run_dir / "hbonds_eq10ns_aaaa1111.json").write_text(eq10_result.model_dump_json())
        (run_dir / "hbonds_eq20ns_bbbb2222.json").write_text(eq20_result.model_dump_json())

        analysis = HydrogenBondsAnalysis()
        loaded = analysis._load_replicate_result(run_dir)

        assert loaded is None

    def test_handles_stat_failure_gracefully(self, tmp_path: Path) -> None:
        """Should return None when globbing cache files raises OSError."""
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()

        rep_result = HydrogenBondResult(
            config_hash="abc123",
            replicate=1,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="chainid A",
            summaries=[],
            composition_entries=[],
        )
        cache_path = run_dir / "hbonds_eq10ns_deadbeef.json"
        cache_path.write_text(rep_result.model_dump_json())

        analysis = HydrogenBondsAnalysis()
        with patch.object(Path, "glob", autospec=True, side_effect=OSError("Permission denied")):
            loaded = analysis._load_replicate_result(run_dir)

        assert loaded is None

    def test_handles_corrupt_cache_gracefully(self, tmp_path: Path) -> None:
        """Should return None when cache file contains invalid JSON."""
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()

        cache_path = run_dir / "hbonds_eq10ns_deadbeef.json"
        cache_path.write_text("THIS IS NOT JSON")

        analysis = HydrogenBondsAnalysis()
        loaded = analysis._load_replicate_result(run_dir)
        assert loaded is None


def test_execution_cost_hint_is_high() -> None:
    """execution_cost_hint should be 'high' — H-bond analysis is expensive."""
    assert HydrogenBondsAnalysis.execution_cost_hint == "high"
