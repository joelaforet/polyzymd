"""Tests for the polymer bridging analysis plugin."""

from __future__ import annotations

from pathlib import Path
from types import ModuleType, SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from pydantic import BaseModel

from polyzymd.analyses.base import AggregateContext, Condition, PlotContext, ReplicateContext
from polyzymd.analyses.polymer_bridging import (
    PolymerBridgingAggregatedResult,
    PolymerBridgingAnalysis,
    PolymerBridgingObservation,
    PolymerBridgingReplicateResult,
    PolymerBridgingSettings,
    _compute_bridging_statistics_from_frames,
)


def _make_mock_mdanalysis_modules() -> dict[str, ModuleType]:
    """Build minimal MDAnalysis modules for lazy distance imports.

    Returns
    -------
    dict[str, ModuleType]
        Modules keyed by import name for patching ``sys.modules``.
    """

    mdanalysis_module = ModuleType("MDAnalysis")
    exceptions_module = ModuleType("MDAnalysis.exceptions")
    lib_module = ModuleType("MDAnalysis.lib")
    distances_module = ModuleType("MDAnalysis.lib.distances")

    class NoDataError(Exception):
        pass

    def _distance_array(a, b, box=None):
        del box
        return np.linalg.norm(
            np.asarray(a)[:, None, :] - np.asarray(b)[None, :, :],
            axis=-1,
        )

    distances_module.capped_distance = MagicMock(return_value=([], []))
    distances_module.distance_array = MagicMock(side_effect=_distance_array)
    exceptions_module.NoDataError = NoDataError
    mdanalysis_module.exceptions = exceptions_module
    lib_module.distances = distances_module
    mdanalysis_module.lib = lib_module
    return {
        "MDAnalysis": mdanalysis_module,
        "MDAnalysis.exceptions": exceptions_module,
        "MDAnalysis.lib": lib_module,
        "MDAnalysis.lib.distances": distances_module,
    }


def _make_hashable_sim_config(tmp_path: Path):
    """Build a lightweight config object compatible with ``compute_config_hash``."""

    return SimpleNamespace(
        name="bridging-test",
        enzyme=SimpleNamespace(name="enzyme", pdb_path=tmp_path / "enzyme.pdb"),
        substrate=None,
        polymers=None,
        thermodynamics=SimpleNamespace(temperature=300.0, pressure=1.0),
        output=SimpleNamespace(
            projects_directory=tmp_path / "projects",
            effective_scratch_directory=tmp_path / "scratch",
            naming_template="run_{replicate}",
        ),
    )


class TestDiscovery:
    def test_discovery_by_name(self):
        from polyzymd.analyses.discovery import get_analysis

        assert get_analysis("polymer_bridging") is PolymerBridgingAnalysis

    def test_discovery_by_alias(self):
        from polyzymd.analyses.discovery import get_analysis

        assert get_analysis("bridging") is PolymerBridgingAnalysis


class TestSettings:
    def test_defaults(self):
        settings = PolymerBridgingSettings()
        assert settings.polymer_selection == "chainID C"
        assert settings.min_ca_distance_angstrom == pytest.approx(0.0)

    def test_negative_distance_rejected(self):
        with pytest.raises(ValueError):
            PolymerBridgingSettings(min_ca_distance_angstrom=-1.0)


class TestCacheFingerprint:
    """Verify cache filenames include settings fingerprint."""

    def test_cache_tag_changes_with_settings(self):
        """Different settings must produce different cache tags."""
        analysis = PolymerBridgingAnalysis()
        tag_default = analysis._make_settings_cache_tag(PolymerBridgingSettings())
        tag_custom = analysis._make_settings_cache_tag(
            PolymerBridgingSettings(cutoff=5.0, min_ca_distance_angstrom=10.0)
        )
        assert tag_default != tag_custom
        assert len(tag_default) == 8
        assert len(tag_custom) == 8

    def test_cache_paths_include_window_identity(self, tmp_path):
        """Replicate and aggregate sidecars should encode the analysis window."""
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()

        replicate_path = analysis._replicate_cache_path(tmp_path / "run_1", settings, "10ns")
        aggregate_path = analysis._aggregate_cache_path(
            tmp_path / "aggregated",
            settings,
            "10ns",
            (1, 2),
        )

        assert "eq10ns" in replicate_path.name
        assert "eq10ns" in aggregate_path.name
        assert "reps1-2" in aggregate_path.name

    def test_cache_window_comparison_normalizes_units(self):
        """Equivalent windows in different units should match cache metadata."""

        result = SimpleNamespace(equilibration_time=1000.0, equilibration_unit="ps")

        assert PolymerBridgingAnalysis._cache_matches_window(result, "1ns")


class TestCoreComputation:
    def test_bridging_stats_without_distance_threshold(self):
        frame_contacts = [{10}, {10, 35}, {10, 35}, {60}]
        stats = _compute_bridging_statistics_from_frames(
            frame_contacts,
            min_ca_distance_angstrom=0.0,
            ca_distances={},
        )
        assert stats["contacting_observations"] == 4
        assert stats["multisite_observations"] == 2
        assert stats["high_valency_observations"] == 0
        assert stats["multisite_fraction"] == pytest.approx(0.5)
        assert stats["mean_contacts_per_contacting_oligomer"] == pytest.approx(1.5)

    def test_bridging_stats_with_distance_threshold(self):
        frame_contacts = [{10}, {10, 35}, {10, 35}, {60}]
        stats = _compute_bridging_statistics_from_frames(
            frame_contacts,
            min_ca_distance_angstrom=100.0,
            ca_distances={(10, 35): 95.0, (10, 60): 190.0, (35, 60): 95.0},
        )
        assert stats["multisite_fraction"] == pytest.approx(0.0)
        assert stats["mean_contacts_per_contacting_oligomer"] == pytest.approx(1.0)

    def test_bridging_stats_use_dynamic_observation_distances(self):
        frame_contacts = [
            ({10, 35}, {(10, 35): 7.0}),
            ({10, 35}, {(10, 35): 9.0}),
            ({60}, {}),
        ]
        stats = _compute_bridging_statistics_from_frames(
            frame_contacts,
            min_ca_distance_angstrom=8.0,
        )
        assert stats["contacting_observations"] == 3
        assert stats["multisite_observations"] == 1
        assert stats["multisite_fraction"] == pytest.approx(1.0 / 3.0)
        assert stats["mean_contacts_per_contacting_oligomer"] == pytest.approx(4.0 / 3.0)

    def test_bridging_stats_capture_anchor_and_signature_metadata(self):
        observation = PolymerBridgingObservation(
            protein_residues={10, 35, 60},
            protein_resnames={10: "PHE", 35: "SER", 60: "LEU"},
            protein_groups={10: "aromatic", 35: "polar", 60: "nonpolar"},
            contacting_polymer_resids={101, 103},
            polymer_resnames={101: "SBM", 102: "EGM", 103: "SBM", 104: "EGM", 105: "SBM"},
            fragment_signature=("SBM", "EGM", "SBM", "EGM", "SBM"),
            ca_distances={(10, 35): 9.0, (10, 60): 12.0, (35, 60): 7.0},
            pair_min_distances={(101, 10): 3.2, (103, 35): 3.8, (103, 60): 4.2},
        )
        stats = _compute_bridging_statistics_from_frames(
            [observation],
            min_ca_distance_angstrom=8.0,
        )

        assert stats["anchor_protein_group_probabilities"]["aromatic"] == pytest.approx(1.0)
        assert stats["polymer_anchor_type_probabilities"]["SBM"] == pytest.approx(1.0)
        assert stats["multivalent_protein_group_probabilities"]["nonpolar"] > 0.0
        assert "SBM-EGM-SBM-EGM-SBM" in stats["fragment_signature_probabilities"]


class TestCAValidation:
    """Verify loud failure when CA atoms are missing for distance filtering."""

    def test_no_ca_atoms_with_distance_filter_raises(self):
        """Raise ValueError when CA filtering is requested but no CA atoms exist."""
        condition = MagicMock()
        condition.sim_config = MagicMock()

        mock_loader = MagicMock()
        mock_universe = MagicMock()

        atom = MagicMock()
        atom.resid = 1
        atom.resname = "ALA"

        mock_protein = MagicMock()
        mock_protein.__len__.return_value = 1
        mock_protein.atoms = [atom]
        mock_residue = MagicMock()
        mock_residue.resid = 1
        mock_ca_selection = MagicMock()
        mock_ca_selection.__len__.return_value = 0
        mock_residue.atoms.select_atoms.return_value = mock_ca_selection
        mock_protein.residues = [mock_residue]

        mock_polymer = MagicMock()
        mock_polymer.__len__.return_value = 1

        mock_universe.trajectory.__len__.return_value = 1
        mock_universe.select_atoms.side_effect = [mock_protein, mock_polymer]
        mock_loader.load_universe.return_value = mock_universe
        mock_loader.get_timestep.return_value = 10.0

        with (
            patch(
                "polyzymd.analyses.shared.loader.TrajectoryLoader",
                return_value=mock_loader,
            ),
            patch.dict("sys.modules", _make_mock_mdanalysis_modules()),
        ):
            with pytest.raises(ValueError, match="(?i)no CA atoms"):
                PolymerBridgingAnalysis._compute_frame_contacts(
                    condition,
                    1,
                    protein_selection="protein",
                    polymer_selection="chainID C",
                    cutoff=4.5,
                    equilibration="0ns",
                    min_ca_distance_angstrom=8.0,
                )

    def test_no_ca_atoms_without_distance_filter_succeeds(self):
        """Allow execution when CA filtering is disabled."""
        condition = MagicMock()
        condition.sim_config = MagicMock()

        mock_loader = MagicMock()
        mock_universe = MagicMock()

        atom = MagicMock()
        atom.resid = 1
        atom.resname = "ALA"

        mock_protein = MagicMock()
        mock_protein.__len__.return_value = 1
        mock_protein.atoms = [atom]
        mock_residue = MagicMock()
        mock_residue.resid = 1
        mock_ca_selection = MagicMock()
        mock_ca_selection.__len__.return_value = 0
        mock_residue.atoms.select_atoms.return_value = mock_ca_selection
        mock_protein.residues = [mock_residue]

        mock_polymer = MagicMock()
        mock_polymer.__len__.return_value = 1
        mock_polymer.fragments = []

        mock_universe.select_atoms.side_effect = [mock_protein, mock_polymer]
        mock_universe.trajectory.__len__.return_value = 1
        mock_universe.trajectory.__getitem__.return_value = [
            SimpleNamespace(dimensions=None, positions=[])
        ]
        mock_loader.load_universe.return_value = mock_universe
        mock_loader.get_timestep.return_value = 10.0

        with (
            patch(
                "polyzymd.analyses.shared.loader.TrajectoryLoader",
                return_value=mock_loader,
            ),
            patch.dict("sys.modules", _make_mock_mdanalysis_modules()),
        ):
            observations, n_frames, timestep_ps = PolymerBridgingAnalysis._compute_frame_contacts(
                condition,
                1,
                protein_selection="protein",
                polymer_selection="chainID C",
                cutoff=4.5,
                equilibration="0ns",
                min_ca_distance_angstrom=0.0,
            )

        assert observations == []
        assert n_frames == 1
        assert timestep_ps == pytest.approx(10.0)

    def test_fragment_lookup_falls_back_without_bond_topology(self):
        """No-bond topologies should be treated as one polymer fragment."""
        from polyzymd.analyses.polymer_bridging._runner import _fragments_or_single

        modules = _make_mock_mdanalysis_modules()
        no_data_error = modules["MDAnalysis.exceptions"].NoDataError

        class _AtomGroup:
            @property
            def fragments(self):
                raise no_data_error("No bond information")

        atoms = _AtomGroup()

        with patch.dict("sys.modules", modules):
            assert _fragments_or_single(atoms, context="test") == [atoms]

    def test_fragment_lookup_propagates_unrelated_errors(self):
        """Only MDAnalysis NoDataError should trigger the no-bond fallback."""
        from polyzymd.analyses.polymer_bridging._runner import _fragments_or_single

        class CustomFragmentError(RuntimeError):
            pass

        class _AtomGroup:
            @property
            def fragments(self):
                raise CustomFragmentError("unexpected fragment failure")

        with patch.dict("sys.modules", _make_mock_mdanalysis_modules()):
            with pytest.raises(CustomFragmentError, match="unexpected fragment failure"):
                _fragments_or_single(_AtomGroup(), context="test")


class TestPBCDistanceHandling:
    """Regression tests for PBC-aware bridging distance calculations."""

    def test_ca_distances_use_periodic_boundary_conditions(self):
        """CA-distance filtering should receive minimum-image distances."""

        from polyzymd.analyses.polymer_bridging._runner import _observation_ca_distances

        def _minimum_image_distance(a, b, box=None):
            delta = np.asarray(a)[:, None, :] - np.asarray(b)[None, :, :]
            if box is not None:
                lengths = np.asarray(box[:3], dtype=float)
                delta -= lengths * np.round(delta / lengths)
            return np.linalg.norm(delta, axis=-1)

        modules = _make_mock_mdanalysis_modules()
        modules["MDAnalysis.lib.distances"].distance_array = MagicMock(
            side_effect=_minimum_image_distance
        )
        positions = np.asarray([[0.2, 0.0, 0.0], [9.8, 0.0, 0.0]], dtype=float)

        with patch.dict("sys.modules", modules):
            distances = _observation_ca_distances(
                {10, 20},
                positions,
                {10: 0, 20: 1},
                np.asarray([10.0, 10.0, 10.0, 90.0, 90.0, 90.0]),
            )

        assert distances[(10, 20)] == pytest.approx(0.4)

    def test_pair_min_distances_use_pbc_distances_from_contact_search(self):
        """Anchor distances should use PBC-aware capped-distance outputs."""

        from polyzymd.analyses.polymer_bridging._runner import _compute_pair_min_distances

        fragment = SimpleNamespace(atoms=[SimpleNamespace(resid=101)])
        protein = SimpleNamespace(atoms=[SimpleNamespace(resid=10)])

        distances = _compute_pair_min_distances(
            fragment,
            protein,
            np.asarray([0], dtype=np.int64),
            np.asarray([0], dtype=np.int64),
            np.asarray([10], dtype=np.int64),
            np.asarray([0.35], dtype=float),
        )

        assert distances[(101, 10)] == pytest.approx(0.35)


def _make_replicate_result(replicate: int) -> PolymerBridgingReplicateResult:
    """Build a minimal polymer bridging replicate result for aggregation tests."""

    return PolymerBridgingReplicateResult(
        replicate=replicate,
        n_frames=4,
        timestep_ps=10.0,
        min_ca_distance_angstrom=0.0,
        contacting_observations=4,
        multisite_observations=2,
        high_valency_observations=0,
        mean_contacts_per_contacting_oligomer=1.5,
        multisite_fraction=0.5,
        high_valency_fraction=0.0,
        valency_probabilities={"1": 0.5, "2": 0.5, "3+": 0.0},
    )


def _make_aggregated_result(
    replicates: tuple[int, ...] = (1, 2),
    *,
    mean_contacts: float = 1.5,
    equilibration_time: float = 0.0,
    settings_fingerprint: str | None = None,
    config_hash: str = "unknown",
) -> PolymerBridgingAggregatedResult:
    """Build a minimal polymer bridging aggregate result for cache tests."""

    n_replicates = len(replicates)
    replicate_values = [mean_contacts] * n_replicates
    return PolymerBridgingAggregatedResult(
        settings_fingerprint=settings_fingerprint,
        config_hash=config_hash,
        equilibration_time=equilibration_time,
        equilibration_unit="ns",
        n_replicates=n_replicates,
        replicates=list(replicates),
        min_ca_distance_angstrom=0.0,
        mean_contacts_per_contacting_oligomer=mean_contacts,
        mean_contacts_sem=0.1,
        multisite_fraction=0.5,
        multisite_fraction_sem=0.05,
        high_valency_fraction=0.1,
        high_valency_fraction_sem=0.02,
        mean_contacts_per_contacting_oligomer_replicates=replicate_values,
        multisite_fraction_replicates=[0.45] * n_replicates,
        high_valency_fraction_replicates=[0.08] * n_replicates,
        valency_probabilities_mean={"1": 0.5, "2": 0.4, "3+": 0.1},
        valency_probabilities_sem={"1": 0.01, "2": 0.02, "3+": 0.01},
        valency_probabilities_per_replicate={"1": [0.52] * n_replicates},
    )


def test_aggregated_result_loads_legacy_payload_without_valency_replicates() -> None:
    """Legacy aggregate caches should load when valency replicate profiles are absent."""

    payload = _make_aggregated_result().model_dump()
    payload.pop("valency_probabilities_per_replicate")

    result = PolymerBridgingAggregatedResult.model_validate(payload)

    assert result.valency_probabilities_per_replicate == {}


class TestLifecycle:
    def test_compute_replicate_uses_runner_seam(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        condition = Condition(
            label="Cond",
            config_path=tmp_path / "condition" / "config.yaml",
            replicates=(1,),
            sim_config=MagicMock(),
        )
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "out",
            equilibration="0ns",
            recompute=False,
            settings=PolymerBridgingSettings(),
            result_path=tmp_path / "out" / "result.json",
        )

        runner_result = SimpleNamespace(
            observations=[
                ({10}, {}),
                ({10, 35}, {(10, 35): 20.0}),
                ({10, 35}, {(10, 35): 20.0}),
                ({60}, {}),
            ],
            n_frames=4,
            timestep_ps=10.0,
        )
        runner = MagicMock()
        runner.results = runner_result
        runner.run.return_value = runner
        mock_universe = MagicMock()
        mock_universe.trajectory.__len__.return_value = 4
        mock_loader = MagicMock()
        mock_loader.load_universe.return_value = mock_universe
        mock_loader.get_timestep.return_value = 10.0

        with (
            patch("polyzymd.analyses.shared.loader.TrajectoryLoader", return_value=mock_loader),
            patch(
                "polyzymd.analyses.polymer_bridging._runner.PolymerBridgingReplicateRunner",
                return_value=runner,
            ) as mock_runner_cls,
        ):
            result = analysis.run_replicate(ctx, 1)

        assert result.replicate == 1
        assert result.multisite_fraction == pytest.approx(0.5)
        assert result.settings_fingerprint == analysis._make_settings_cache_tag(ctx.settings)
        runner.run.assert_called_once_with(start=0, stop=4, step=1)
        mock_runner_cls.assert_called_once()
        assert ctx.result_path.exists()
        assert list(ctx.output_dir.glob("polymer_bridging_*.json"))

    def test_compute_replicate_uses_cached_result_when_available(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        condition = Condition(
            label="Cond",
            config_path=tmp_path / "condition" / "config.yaml",
            replicates=(1,),
            sim_config=MagicMock(),
        )
        result_path = analysis._replicate_cache_path(tmp_path / "out", settings, "0ns")
        cached = PolymerBridgingReplicateResult(
            settings_fingerprint=analysis._make_settings_cache_tag(settings),
            replicate=1,
            n_frames=4,
            timestep_ps=10.0,
            min_ca_distance_angstrom=0.0,
            contacting_observations=4,
            multisite_observations=2,
            high_valency_observations=0,
            mean_contacts_per_contacting_oligomer=1.5,
            multisite_fraction=0.5,
            high_valency_fraction=0.0,
            valency_probabilities={"1": 0.5, "2": 0.5, "3+": 0.0},
        )
        result_path.parent.mkdir(parents=True, exist_ok=True)
        cached.save(result_path)

        analysis._compute_frame_contacts = MagicMock(side_effect=AssertionError("should not run"))

        loaded = analysis.run_replicate(
            ReplicateContext(
                condition=condition,
                replicate=1,
                sim_config=condition.sim_config,
                output_dir=tmp_path / "out",
                equilibration="0ns",
                recompute=False,
                settings=settings,
            ),
            1,
        )

        assert loaded == cached

    def test_compute_replicate_prefers_sidecar_over_canonical_cache(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        condition = Condition(
            label="Cond",
            config_path=tmp_path / "condition" / "config.yaml",
            replicates=(1,),
            sim_config=MagicMock(),
        )
        output_dir = tmp_path / "out"
        canonical_path = output_dir / "result.json"
        sidecar_path = analysis._replicate_cache_path(output_dir, settings, "0ns")
        canonical = _make_replicate_result(1).model_copy(
            update={
                "settings_fingerprint": analysis._make_settings_cache_tag(settings),
                "mean_contacts_per_contacting_oligomer": 9.0,
            }
        )
        sidecar = _make_replicate_result(1).model_copy(
            update={
                "settings_fingerprint": analysis._make_settings_cache_tag(settings),
                "mean_contacts_per_contacting_oligomer": 1.5,
            }
        )
        canonical.save(canonical_path)
        sidecar.save(sidecar_path)

        loaded = analysis.run_replicate(
            ReplicateContext(
                condition=condition,
                replicate=1,
                sim_config=condition.sim_config,
                output_dir=output_dir,
                equilibration="0ns",
                recompute=False,
                settings=settings,
                result_path=canonical_path,
            ),
            1,
        )

        assert loaded.mean_contacts_per_contacting_oligomer == pytest.approx(1.5)

    def test_compute_replicate_rejects_missing_fingerprint_sidecar(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        condition = Condition(
            label="Cond",
            config_path=tmp_path / "condition" / "config.yaml",
            replicates=(1,),
            sim_config=MagicMock(),
        )
        output_dir = tmp_path / "out"
        sidecar_path = analysis._replicate_cache_path(output_dir, settings, "0ns")
        _make_replicate_result(1).save(sidecar_path)

        with patch(
            "polyzymd.analyses.base.Analysis._run_replicate_default",
            side_effect=RuntimeError("should recompute"),
        ):
            with pytest.raises(RuntimeError, match="should recompute"):
                analysis.run_replicate(
                    ReplicateContext(
                        condition=condition,
                        replicate=1,
                        sim_config=condition.sim_config,
                        output_dir=output_dir,
                        equilibration="0ns",
                        recompute=False,
                        settings=settings,
                    ),
                    1,
                )

    def test_compute_replicate_rejects_and_does_not_copy_missing_fingerprint_canonical(
        self,
        tmp_path,
    ):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        condition = Condition(
            label="Cond",
            config_path=tmp_path / "condition" / "config.yaml",
            replicates=(1,),
            sim_config=MagicMock(),
        )
        output_dir = tmp_path / "out"
        canonical_path = output_dir / "result.json"
        sidecar_path = analysis._replicate_cache_path(output_dir, settings, "0ns")
        _make_replicate_result(1).save(canonical_path)

        with patch(
            "polyzymd.analyses.base.Analysis._run_replicate_default",
            side_effect=RuntimeError("should recompute"),
        ):
            with pytest.raises(RuntimeError, match="should recompute"):
                analysis.run_replicate(
                    ReplicateContext(
                        condition=condition,
                        replicate=1,
                        sim_config=condition.sim_config,
                        output_dir=output_dir,
                        equilibration="0ns",
                        recompute=False,
                        settings=settings,
                        result_path=canonical_path,
                    ),
                    1,
                )

        assert not sidecar_path.exists()

    def test_compute_replicate_ignores_window_mismatched_cache(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        condition = Condition(
            label="Cond",
            config_path=tmp_path / "condition" / "config.yaml",
            replicates=(1,),
            sim_config=MagicMock(),
        )
        result_path = analysis._replicate_cache_path(tmp_path / "out", settings, "0ns")
        cached = _make_replicate_result(1)
        cached.equilibration_time = 0.0
        cached.equilibration_unit = "ns"
        result_path.parent.mkdir(parents=True, exist_ok=True)
        cached.save(result_path)
        fresh = _make_replicate_result(1)
        fresh.equilibration_time = 10.0
        fresh.equilibration_unit = "ns"

        with patch(
            "polyzymd.analyses.base.Analysis._run_replicate_default",
            return_value=fresh,
        ) as mock_compute:
            loaded = analysis.run_replicate(
                ReplicateContext(
                    condition=condition,
                    replicate=1,
                    sim_config=condition.sim_config,
                    output_dir=tmp_path / "out",
                    equilibration="10ns",
                    recompute=False,
                    settings=settings,
                ),
                1,
            )

        assert loaded is fresh
        mock_compute.assert_called_once()

    def test_compute_replicate_ignores_replicate_id_mismatched_cache(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        condition = Condition(
            label="Cond",
            config_path=tmp_path / "condition" / "config.yaml",
            replicates=(1,),
            sim_config=MagicMock(),
        )
        result_path = analysis._replicate_cache_path(tmp_path / "out", settings, "0ns")
        cached = _make_replicate_result(2).model_copy(
            update={"settings_fingerprint": analysis._make_settings_cache_tag(settings)}
        )
        result_path.parent.mkdir(parents=True, exist_ok=True)
        cached.save(result_path)
        fresh = _make_replicate_result(1)

        with patch(
            "polyzymd.analyses.base.Analysis._run_replicate_default",
            return_value=fresh,
        ) as mock_compute:
            loaded = analysis.run_replicate(
                ReplicateContext(
                    condition=condition,
                    replicate=1,
                    sim_config=condition.sim_config,
                    output_dir=tmp_path / "out",
                    equilibration="0ns",
                    recompute=False,
                    settings=settings,
                ),
                1,
            )

        assert loaded is fresh
        mock_compute.assert_called_once()

    def test_aggregate_rejects_replicate_id_mismatched_inputs(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        condition = Condition(
            label="Cond",
            config_path=tmp_path / "condition" / "config.yaml",
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2),
            output_dir=tmp_path / "aggregated",
            equilibration="0ns",
            settings=settings,
        )

        with pytest.raises(ValueError, match="replicate IDs do not match"):
            analysis.aggregate(ctx, [_make_replicate_result(1), _make_replicate_result(3)])


class TestConditionHasPolymer:
    """Regression tests for polymer-selection-aware condition filtering."""

    def test_topology_chain_check_uses_requested_polymer_selection(self):
        """A non-C polymer selection should match the requested topology chain."""

        from polyzymd.analyses.polymer_bridging import _condition_has_polymer

        sim_config = MagicMock()
        sim_config.polymers = None
        sim_config.topology.chains = [SimpleNamespace(chain_id="B")]
        condition = SimpleNamespace(sim_config=sim_config, replicates=())

        assert _condition_has_polymer(condition, polymer_selection="chainID B") is True

    def test_topology_shortcut_ignores_complex_polymer_selection(self):
        """A chain-only topology hit should not satisfy a richer selection."""

        from polyzymd.analyses.polymer_bridging import _condition_has_polymer

        sim_config = MagicMock()
        sim_config.polymers = None
        sim_config.topology.chains = [SimpleNamespace(chain_id="C")]
        condition = SimpleNamespace(sim_config=sim_config, replicates=())

        assert (
            _condition_has_polymer(condition, polymer_selection="chainID C and resname PEG")
            is False
        )

    def test_topology_chain_check_does_not_hardcode_chain_c(self):
        """Chain C should not satisfy a different polymer selection."""

        from polyzymd.analyses.polymer_bridging import _condition_has_polymer

        sim_config = MagicMock()
        sim_config.polymers = None
        sim_config.topology.chains = [SimpleNamespace(chain_id="C")]
        condition = SimpleNamespace(sim_config=sim_config, replicates=())

        assert _condition_has_polymer(condition, polymer_selection="chainID B") is False

    def test_enabled_polymer_config_does_not_bypass_requested_selection(self):
        """Enabled polymer config should not override selection-based detection."""

        from polyzymd.analyses.polymer_bridging import _condition_has_polymer

        sim_config = MagicMock()
        sim_config.polymers = SimpleNamespace(enabled=True)
        sim_config.topology.chains = [SimpleNamespace(chain_id="C")]
        condition = SimpleNamespace(sim_config=sim_config, replicates=())

        assert _condition_has_polymer(condition, polymer_selection="chainID B") is False

    def test_aggregate_builds_typed_result(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        condition = Condition(
            label="Cond",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2),
            output_dir=tmp_path / "agg",
            equilibration="0ns",
            settings=PolymerBridgingSettings(),
            result_path=tmp_path / "agg" / "result.json",
        )
        results = [
            _make_replicate_result(1),
            _make_replicate_result(2),
        ]

        aggregated = analysis.aggregate(ctx, results)

        assert isinstance(aggregated, PolymerBridgingAggregatedResult)
        assert aggregated.n_replicates == 2
        assert ctx.result_path.exists()

    def test_aggregate_uses_cached_result_when_available(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        condition = Condition(
            label="Cond",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        result_path = analysis._aggregate_cache_path(tmp_path / "agg", settings, "0ns", (1, 2))
        cached = PolymerBridgingAggregatedResult(
            settings_fingerprint=analysis._make_settings_cache_tag(settings),
            n_replicates=2,
            replicates=[1, 2],
            min_ca_distance_angstrom=0.0,
            mean_contacts_per_contacting_oligomer=1.5,
            mean_contacts_sem=0.1,
            multisite_fraction=0.5,
            multisite_fraction_sem=0.05,
            high_valency_fraction=0.1,
            high_valency_fraction_sem=0.02,
            mean_contacts_per_contacting_oligomer_replicates=[1.4, 1.6],
            multisite_fraction_replicates=[0.45, 0.55],
            high_valency_fraction_replicates=[0.08, 0.12],
            valency_probabilities_mean={"1": 0.5, "2": 0.4, "3+": 0.1},
            valency_probabilities_sem={"1": 0.01, "2": 0.02, "3+": 0.01},
            valency_probabilities_per_replicate={
                "1": [0.52, 0.48],
                "2": [0.38, 0.42],
                "3+": [0.1, 0.1],
            },
        )
        result_path.parent.mkdir(parents=True, exist_ok=True)
        cached.save(result_path)

        loaded = analysis.aggregate(
            AggregateContext(
                condition=condition,
                replicates=(1, 2),
                output_dir=tmp_path / "agg",
                equilibration="0ns",
                settings=settings,
            ),
            [],
        )

        assert loaded == cached

    def test_aggregate_prefers_sidecar_over_canonical_cache(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        condition = Condition(
            label="Cond",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2),
            output_dir=tmp_path / "agg",
            equilibration="0ns",
            settings=settings,
            result_path=tmp_path / "agg" / "result.json",
        )
        canonical = _make_aggregated_result(
            (1, 2),
            mean_contacts=9.0,
            settings_fingerprint=analysis._make_settings_cache_tag(settings),
        )
        sidecar = _make_aggregated_result(
            (1, 2),
            mean_contacts=1.5,
            settings_fingerprint=analysis._make_settings_cache_tag(settings),
        )
        canonical.save(ctx.result_path)
        sidecar.save(analysis._aggregate_cache_path(ctx.output_dir, settings, "0ns", (1, 2)))

        loaded = analysis.aggregate(ctx, [])

        assert loaded.mean_contacts_per_contacting_oligomer == pytest.approx(1.5)

    def test_aggregate_rejects_replicate_set_mismatch_cache(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        condition = Condition(
            label="Cond",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2, 3),
            output_dir=tmp_path / "agg",
            equilibration="0ns",
            settings=settings,
            result_path=tmp_path / "agg" / "result.json",
        )
        stale = _make_aggregated_result(
            (1, 2),
            mean_contacts=9.0,
            settings_fingerprint=analysis._make_settings_cache_tag(settings),
        )
        current_path = analysis._aggregate_cache_path(ctx.output_dir, settings, "0ns", (1, 2, 3))
        stale.save(current_path)

        loaded = analysis.aggregate(
            ctx,
            [_make_replicate_result(1), _make_replicate_result(2), _make_replicate_result(3)],
        )

        assert loaded.mean_contacts_per_contacting_oligomer == pytest.approx(1.5)
        assert loaded.replicates == [1, 2, 3]

    def test_aggregate_ignores_window_mismatched_cache(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        condition = Condition(
            label="Cond",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        result_path = analysis._aggregate_cache_path(tmp_path / "agg", settings, "0ns", (1, 2))
        cached = PolymerBridgingAggregatedResult(
            n_replicates=2,
            replicates=[1, 2],
            min_ca_distance_angstrom=0.0,
            mean_contacts_per_contacting_oligomer=9.9,
            mean_contacts_sem=0.1,
            multisite_fraction=0.5,
            multisite_fraction_sem=0.05,
            high_valency_fraction=0.1,
            high_valency_fraction_sem=0.02,
            mean_contacts_per_contacting_oligomer_replicates=[9.9, 9.9],
            multisite_fraction_replicates=[0.45, 0.55],
            high_valency_fraction_replicates=[0.08, 0.12],
            valency_probabilities_mean={"1": 0.5, "2": 0.4, "3+": 0.1},
            valency_probabilities_sem={"1": 0.01, "2": 0.02, "3+": 0.01},
            valency_probabilities_per_replicate={"1": [0.52, 0.48]},
            equilibration_time=0.0,
            equilibration_unit="ns",
        )
        result_path.parent.mkdir(parents=True, exist_ok=True)
        cached.save(result_path)

        loaded = analysis.aggregate(
            AggregateContext(
                condition=condition,
                replicates=(1, 2),
                output_dir=tmp_path / "agg",
                equilibration="10ns",
                settings=settings,
            ),
            [_make_replicate_result(1), _make_replicate_result(2)],
        )

        assert loaded is not cached
        assert loaded.equilibration_time == pytest.approx(10.0)

    def test_aggregate_rejects_stale_config_hash_cache(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        sim_config = _make_hashable_sim_config(tmp_path)
        condition = Condition(
            label="Cond",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=sim_config,
        )
        result_path = analysis._aggregate_cache_path(tmp_path / "agg", settings, "0ns", (1, 2))
        cached = PolymerBridgingAggregatedResult(
            settings_fingerprint=analysis._make_settings_cache_tag(settings),
            config_hash="stale-config",
            n_replicates=2,
            replicates=[1, 2],
            min_ca_distance_angstrom=0.0,
            mean_contacts_per_contacting_oligomer=9.9,
            mean_contacts_sem=0.1,
            multisite_fraction=0.5,
            multisite_fraction_sem=0.05,
            high_valency_fraction=0.1,
            high_valency_fraction_sem=0.02,
            mean_contacts_per_contacting_oligomer_replicates=[9.9, 9.9],
            multisite_fraction_replicates=[0.45, 0.55],
            high_valency_fraction_replicates=[0.08, 0.12],
            valency_probabilities_mean={"1": 0.5, "2": 0.4, "3+": 0.1},
            valency_probabilities_sem={"1": 0.01, "2": 0.02, "3+": 0.01},
            valency_probabilities_per_replicate={"1": [0.52, 0.48]},
        )
        result_path.parent.mkdir(parents=True, exist_ok=True)
        cached.save(result_path)

        loaded = analysis.aggregate(
            AggregateContext(
                condition=condition,
                replicates=(1, 2),
                output_dir=tmp_path / "agg",
                equilibration="0ns",
                settings=settings,
            ),
            [_make_replicate_result(1), _make_replicate_result(2)],
        )

        assert loaded is not cached
        assert loaded.mean_contacts_per_contacting_oligomer == pytest.approx(1.5)

    def test_aggregate_cache_checks_include_sim_config(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        sim_config = MagicMock()
        condition = Condition(
            label="Cond",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=sim_config,
        )
        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2),
            output_dir=tmp_path / "agg",
            equilibration="0ns",
            settings=settings,
            result_path=tmp_path / "agg" / "result.json",
        )

        with patch.object(analysis, "_check_cache", return_value=None) as mock_check_cache:
            analysis.aggregate(ctx, [_make_replicate_result(1), _make_replicate_result(2)])

        assert mock_check_cache.call_count == 2
        assert mock_check_cache.call_args_list[0].kwargs["sim_config"] is sim_config
        assert mock_check_cache.call_args_list[1].kwargs["sim_config"] is sim_config

    def test_plot_returns_paths(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        aggregated = PolymerBridgingAggregatedResult(
            n_replicates=2,
            replicates=[1, 2],
            min_ca_distance_angstrom=0.0,
            mean_contacts_per_contacting_oligomer=1.5,
            mean_contacts_sem=0.1,
            multisite_fraction=0.5,
            multisite_fraction_sem=0.05,
            high_valency_fraction=0.1,
            high_valency_fraction_sem=0.02,
            mean_contacts_per_contacting_oligomer_replicates=[1.4, 1.6],
            multisite_fraction_replicates=[0.45, 0.55],
            high_valency_fraction_replicates=[0.08, 0.12],
            valency_probabilities_mean={"1": 0.5, "2": 0.4, "3+": 0.1},
            valency_probabilities_sem={"1": 0.01, "2": 0.02, "3+": 0.01},
            valency_probabilities_per_replicate={
                "1": [0.52, 0.48],
                "2": [0.38, 0.42],
                "3+": [0.1, 0.1],
            },
        )
        comparison_json = tmp_path / "comparison" / "polymer_bridging" / "result.json"
        comparison_json.parent.mkdir(parents=True, exist_ok=True)
        comparison_json.write_text(
            '{"analysis_type":"polymer_bridging","name":"Test","control_label":null,'
            '"conditions":[{"label":"A","n_replicates":2,'
            '"multisite_fraction_mean":0.5,"multisite_fraction_sem":0.05,'
            '"multisite_fraction_replicate_values":[0.45,0.55],'
            '"mean_contacts_per_contacting_oligomer_mean":1.5,'
            '"mean_contacts_per_contacting_oligomer_sem":0.1,'
            '"mean_contacts_per_contacting_oligomer_replicate_values":[1.4,1.6],'
            '"high_valency_fraction_mean":0.1,"high_valency_fraction_sem":0.02,'
            '"high_valency_fraction_replicate_values":[0.08,0.12]},'
            '{"label":"B","n_replicates":2,'
            '"multisite_fraction_mean":0.7,"multisite_fraction_sem":0.04,'
            '"multisite_fraction_replicate_values":[0.66,0.74],'
            '"mean_contacts_per_contacting_oligomer_mean":1.8,'
            '"mean_contacts_per_contacting_oligomer_sem":0.08,'
            '"mean_contacts_per_contacting_oligomer_replicate_values":[1.72,1.88],'
            '"high_valency_fraction_mean":0.2,"high_valency_fraction_sem":0.03,'
            '"high_valency_fraction_replicate_values":[0.17,0.23]}],'
            '"pairwise_comparisons":[],"anova":null,"ranking":["B","A"],'
            '"rankings_by_metric":{"multisite_fraction":["B","A"],'
            '"mean_contacts_per_contacting_oligomer":["B","A"],'
            '"high_valency_fraction":["B","A"]},"equilibration_time":"0ns",'
            '"created_at":"now","polyzymd_version":"test"}'
        )
        cond_a = Condition(
            label="A", config_path=Path("/tmp/a.yaml"), replicates=(1, 2, 3), sim_config=MagicMock()
        )
        cond_b = Condition(
            label="B", config_path=Path("/tmp/b.yaml"), replicates=(1, 2, 3), sim_config=MagicMock()
        )
        analysis_dir_a = tmp_path / "analysis" / "A" / "polymer_bridging" / "aggregated"
        analysis_dir_b = tmp_path / "analysis" / "B" / "polymer_bridging" / "aggregated"
        analysis_dir_a.mkdir(parents=True, exist_ok=True)
        analysis_dir_b.mkdir(parents=True, exist_ok=True)
        aggregated.save(analysis_dir_a / "result.json")
        aggregated.model_copy(
            update={
                "multisite_fraction": 0.7,
                "mean_contacts_per_contacting_oligomer": 1.8,
                "high_valency_fraction": 0.2,
                "valency_probabilities_mean": {"1": 0.3, "2": 0.5, "3+": 0.2},
            }
        ).save(analysis_dir_b / "result.json")
        from polyzymd.config.comparison import PlotSettings

        ctx = PlotContext(
            conditions=[cond_a, cond_b],
            analysis_dirs={"A": analysis_dir_a.parent, "B": analysis_dir_b.parent},
            results_dir=comparison_json.parent,
            output_dir=tmp_path / "figures",
            settings=PolymerBridgingSettings(),
            plot_settings=PlotSettings(),
            comparison_path=comparison_json,
        )

        paths = analysis.plot(ctx)

        assert len(paths) >= 3
        assert all(path.exists() for path in paths)

    def test_plot_stacked_bars_overlay_component_replicates(self, tmp_path):
        """Stacked probability plots should overlay component replicate dots."""
        from polyzymd.analyses.base import ComparisonResult, ConditionSummary
        from polyzymd.config.comparison import PlotSettings

        analysis = PolymerBridgingAnalysis()
        condition = Condition(
            label="A",
            config_path=Path("/tmp/a.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        comparison = ComparisonResult(
            analysis_type="polymer_bridging",
            name="Test",
            conditions=[ConditionSummary(label="A", n_replicates=2)],
            pairwise_comparisons=[],
            ranking=["A"],
        )
        summary = _make_aggregated_result((1, 2)).model_copy(
            update={
                "valency_probabilities_mean": {"1": 0.6, "2": 0.4, "3+": 0.0},
                "valency_probabilities_per_replicate": {
                    "1": [0.58, 0.62],
                    "2": [0.42, 0.38],
                },
                "multivalent_protein_group_probabilities_mean": {
                    "aromatic": 0.25,
                    "polar": 0.75,
                },
                "multivalent_protein_group_probabilities_per_replicate": {
                    "aromatic": [0.2, 0.3],
                    "polar": [0.8, 0.7],
                },
            }
        )
        plot_settings = PlotSettings(
            polymer_bridging={
                "generate_multisite_bars": False,
                "generate_mean_contacts_bars": False,
                "generate_valency_stack": True,
                "generate_anchor_group_bars": False,
                "generate_protein_group_stack": True,
                "generate_anchor_peripheral_heatmap": False,
                "generate_polymer_anchor_heatmap": False,
                "generate_fragment_signature_bars": False,
            }
        )
        ctx = PlotContext(
            conditions=[condition],
            analysis_dirs={"A": tmp_path / "analysis" / "A" / "polymer_bridging"},
            results_dir=tmp_path / "comparison",
            output_dir=tmp_path / "figures",
            settings=PolymerBridgingSettings(),
            plot_settings=plot_settings,
        )

        with (
            patch.object(analysis, "_load_comparison_result", return_value=comparison),
            patch.object(
                analysis,
                "_load_validated_aggregated_result",
                side_effect=[summary, summary, summary, summary],
            ),
            patch(
                "polyzymd.analyses.polymer_bridging.scatter_stacked_segment_replicates"
            ) as scatter,
            patch(
                "polyzymd.analyses.polymer_bridging.save_figure",
                side_effect=lambda _fig, path, *_args, **_kwargs: path,
            ),
        ):
            analysis.plot(ctx)

        assert scatter.call_count == 4
        assert scatter.call_args_list[0].args[2] == pytest.approx(0.0)
        assert scatter.call_args_list[0].args[3] == [0.58, 0.62]
        assert scatter.call_args_list[0].kwargs["replicate_base_values"] == [0.0, 0.0]
        assert scatter.call_args_list[1].args[2] == pytest.approx(0.6)
        assert scatter.call_args_list[1].args[3] == [0.42, 0.38]
        assert scatter.call_args_list[1].kwargs["replicate_base_values"] == [0.58, 0.62]
        assert scatter.call_args_list[2].args[2] == pytest.approx(0.0)
        assert scatter.call_args_list[2].args[3] == [0.2, 0.3]
        assert scatter.call_args_list[2].kwargs["replicate_base_values"] == [0.0, 0.0]
        assert scatter.call_args_list[3].args[2] == pytest.approx(0.25)
        assert scatter.call_args_list[3].args[3] == [0.8, 0.7]
        assert scatter.call_args_list[3].kwargs["replicate_base_values"] == [0.2, 0.3]

    def test_format_emits_sections(self):
        analysis = PolymerBridgingAnalysis()
        from polyzymd.analyses.base import ComparisonResult, ConditionSummary

        result = ComparisonResult(
            analysis_type="polymer_bridging",
            name="Test",
            conditions=[
                ConditionSummary(
                    label="A",
                    n_replicates=2,
                    multisite_fraction_mean=0.5,
                    multisite_fraction_sem=0.05,
                    multisite_fraction_replicate_values=[0.45, 0.55],
                    mean_contacts_per_contacting_oligomer_mean=1.5,
                    mean_contacts_per_contacting_oligomer_sem=0.1,
                    mean_contacts_per_contacting_oligomer_replicate_values=[1.4, 1.6],
                    high_valency_fraction_mean=0.1,
                    high_valency_fraction_sem=0.02,
                    high_valency_fraction_replicate_values=[0.08, 0.12],
                ),
                ConditionSummary(
                    label="B",
                    n_replicates=2,
                    multisite_fraction_mean=0.7,
                    multisite_fraction_sem=0.04,
                    multisite_fraction_replicate_values=[0.66, 0.74],
                    mean_contacts_per_contacting_oligomer_mean=1.8,
                    mean_contacts_per_contacting_oligomer_sem=0.08,
                    mean_contacts_per_contacting_oligomer_replicate_values=[1.72, 1.88],
                    high_valency_fraction_mean=0.2,
                    high_valency_fraction_sem=0.03,
                    high_valency_fraction_replicate_values=[0.17, 0.23],
                ),
            ],
            pairwise_comparisons=[],
            anova=None,
            ranking=["B", "A"],
            rankings_by_metric={
                "multisite_fraction": ["B", "A"],
                "mean_contacts_per_contacting_oligomer": ["B", "A"],
                "high_valency_fraction": ["B", "A"],
            },
            equilibration_time="0ns",
            created_at="now",
            polyzymd_version="test",
        )

        text = analysis.format(result, output_format="text")

        assert "WARNING: Experimental analysis" in text
        assert "Polymer Bridging Comparison" in text
        assert "Average Oligomer Valency" in text
        assert "High-Valency Oligomers" in text

    def test_compare_sanitizes_nan_stats(self):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        condition_a = Condition(
            label="A",
            config_path=Path("/tmp/a.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        condition_b = Condition(
            label="B",
            config_path=Path("/tmp/b.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        aggregated_a = PolymerBridgingAggregatedResult(
            settings_fingerprint=analysis._make_settings_cache_tag(settings),
            n_replicates=3,
            replicates=[1, 2, 3],
            min_ca_distance_angstrom=0.0,
            mean_contacts_per_contacting_oligomer=1.0,
            mean_contacts_sem=0.0,
            multisite_fraction=0.0,
            multisite_fraction_sem=0.0,
            high_valency_fraction=0.0,
            high_valency_fraction_sem=0.0,
            mean_contacts_per_contacting_oligomer_replicates=[1.0, 1.0, 1.0],
            multisite_fraction_replicates=[0.0, 0.0, 0.0],
            high_valency_fraction_replicates=[0.0, 0.0, 0.0],
            valency_probabilities_mean={"1": 1.0, "2": 0.0, "3+": 0.0},
            valency_probabilities_sem={"1": 0.0, "2": 0.0, "3+": 0.0},
            valency_probabilities_per_replicate={
                "1": [1.0, 1.0, 1.0],
                "2": [0.0, 0.0, 0.0],
                "3+": [0.0, 0.0, 0.0],
            },
        )
        aggregated_b = aggregated_a.model_copy(
            update={
                "mean_contacts_per_contacting_oligomer": 2.0,
                "mean_contacts_per_contacting_oligomer_replicates": [2.0, 2.0, 2.0],
                "multisite_fraction": 1.0,
                "multisite_fraction_replicates": [1.0, 1.0, 1.0],
                "high_valency_fraction": 1.0,
                "high_valency_fraction_replicates": [1.0, 1.0, 1.0],
                "valency_probabilities_mean": {"1": 0.0, "2": 0.0, "3+": 1.0},
                "valency_probabilities_per_replicate": {
                    "1": [0.0, 0.0, 0.0],
                    "2": [0.0, 0.0, 0.0],
                    "3+": [1.0, 1.0, 1.0],
                },
            }
        )
        from polyzymd.analyses.base import ComparisonContext

        ctx = ComparisonContext(
            name="Test",
            conditions=[condition_a, condition_b],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"A": Path("/tmp/a"), "B": Path("/tmp/b")},
            results_dir=Path("/tmp/results"),
            equilibration="0ns",
            settings=settings,
            recompute=False,
            aggregated_results={"A": aggregated_a, "B": aggregated_b},
        )

        result = analysis.compare(ctx)

        assert result is not None
        assert all(pair.p_value >= 0.0 for pair in result.pairwise_comparisons)
        for pair in result.pairwise_comparisons:
            assert pair.direction == "not_testable"
            assert pair.effect_size_interpretation == "not_testable"

    def test_compare_accepts_successful_replicate_subset_aggregate(self, tmp_path):
        """Finalize should accept aggregates that record successful replicate subsets."""
        from polyzymd.analyses.base import ComparisonContext

        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        condition = Condition(
            label="A",
            config_path=Path("/tmp/a.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        aggregated = _make_aggregated_result(
            (1, 2),
            mean_contacts=1.5,
            settings_fingerprint=analysis._make_settings_cache_tag(settings),
        )
        ctx = ComparisonContext(
            name="Test",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"A": tmp_path / "A" / "polymer_bridging"},
            results_dir=tmp_path / "results",
            equilibration="0ns",
            settings=settings,
            recompute=False,
            aggregated_results={"A": aggregated},
        )

        result = analysis.compare(ctx)

        assert result is not None
        assert result.conditions[0].n_replicates == 2

    def test_compare_rejects_replicate_subset_below_minimum(self, tmp_path):
        """Finalize should not accept successful subsets below plugin minimum."""
        from polyzymd.analyses.base import ComparisonContext

        analysis = PolymerBridgingAnalysis()
        condition = Condition(
            label="A",
            config_path=Path("/tmp/a.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        aggregated = _make_aggregated_result((1,), mean_contacts=1.5)
        ctx = ComparisonContext(
            name="Test",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"A": tmp_path / "A" / "polymer_bridging"},
            results_dir=tmp_path / "results",
            equilibration="0ns",
            settings=PolymerBridgingSettings(),
            recompute=False,
            aggregated_results={"A": aggregated},
        )

        result = analysis.compare(ctx)

        assert result is None

    def test_compare_rejects_preloaded_settings_mismatch(self, tmp_path):
        """Finalize should not use preloaded aggregates from different settings."""
        from polyzymd.analyses.base import ComparisonContext
        from polyzymd.analyses.shared.config_hash import settings_fingerprint

        analysis = PolymerBridgingAnalysis()
        condition = Condition(
            label="A",
            config_path=Path("/tmp/a.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        stale_settings = PolymerBridgingSettings(cutoff=9.0)
        stale_summary = PolymerBridgingAggregatedResult(
            settings_fingerprint=settings_fingerprint(stale_settings),
            n_replicates=2,
            replicates=[1, 2],
            min_ca_distance_angstrom=0.0,
            mean_contacts_per_contacting_oligomer=1.0,
            mean_contacts_sem=0.1,
            multisite_fraction=0.5,
            multisite_fraction_sem=0.05,
            high_valency_fraction=0.1,
            high_valency_fraction_sem=0.02,
            mean_contacts_per_contacting_oligomer_replicates=[0.9, 1.1],
            multisite_fraction_replicates=[0.45, 0.55],
            high_valency_fraction_replicates=[0.08, 0.12],
            valency_probabilities_mean={"1": 0.5, "2": 0.4, "3+": 0.1},
            valency_probabilities_sem={"1": 0.01, "2": 0.02, "3+": 0.01},
            valency_probabilities_per_replicate={"1": [0.5, 0.5]},
        )
        ctx = ComparisonContext(
            name="Test",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"A": tmp_path / "A" / "polymer_bridging"},
            results_dir=tmp_path / "results",
            equilibration="0ns",
            settings=PolymerBridgingSettings(),
            recompute=False,
            aggregated_results={"A": stale_summary},
        )

        assert analysis.compare(ctx) is None

    def test_compare_keeps_normal_direction_for_testable_stats(self):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        condition_a = Condition(
            label="A",
            config_path=Path("/tmp/a.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        condition_b = Condition(
            label="B",
            config_path=Path("/tmp/b.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        aggregated_a = PolymerBridgingAggregatedResult(
            settings_fingerprint=analysis._make_settings_cache_tag(settings),
            n_replicates=3,
            replicates=[1, 2, 3],
            min_ca_distance_angstrom=0.0,
            mean_contacts_per_contacting_oligomer=1.0,
            mean_contacts_sem=0.057735,
            multisite_fraction=0.15,
            multisite_fraction_sem=0.028868,
            high_valency_fraction=0.07,
            high_valency_fraction_sem=0.008819,
            mean_contacts_per_contacting_oligomer_replicates=[0.9, 1.0, 1.1],
            multisite_fraction_replicates=[0.1, 0.15, 0.2],
            high_valency_fraction_replicates=[0.06, 0.07, 0.08],
            valency_probabilities_mean={"1": 0.7, "2": 0.2, "3+": 0.1},
            valency_probabilities_sem={"1": 0.02, "2": 0.01, "3+": 0.01},
            valency_probabilities_per_replicate={
                "1": [0.72, 0.7, 0.68],
                "2": [0.18, 0.2, 0.22],
                "3+": [0.1, 0.1, 0.1],
            },
        )
        aggregated_b = PolymerBridgingAggregatedResult(
            settings_fingerprint=analysis._make_settings_cache_tag(settings),
            n_replicates=3,
            replicates=[1, 2, 3],
            min_ca_distance_angstrom=0.0,
            mean_contacts_per_contacting_oligomer=1.6,
            mean_contacts_sem=0.057735,
            multisite_fraction=0.45,
            multisite_fraction_sem=0.028868,
            high_valency_fraction=0.17,
            high_valency_fraction_sem=0.008819,
            mean_contacts_per_contacting_oligomer_replicates=[1.5, 1.6, 1.7],
            multisite_fraction_replicates=[0.4, 0.45, 0.5],
            high_valency_fraction_replicates=[0.16, 0.17, 0.18],
            valency_probabilities_mean={"1": 0.4, "2": 0.3, "3+": 0.3},
            valency_probabilities_sem={"1": 0.02, "2": 0.01, "3+": 0.01},
            valency_probabilities_per_replicate={
                "1": [0.42, 0.4, 0.38],
                "2": [0.28, 0.3, 0.32],
                "3+": [0.3, 0.3, 0.3],
            },
        )
        from polyzymd.analyses.base import ComparisonContext

        ctx = ComparisonContext(
            name="Test",
            conditions=[condition_a, condition_b],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"A": Path("/tmp/a"), "B": Path("/tmp/b")},
            results_dir=Path("/tmp/results"),
            equilibration="0ns",
            settings=settings,
            recompute=False,
            aggregated_results={"A": aggregated_a, "B": aggregated_b},
        )

        result = analysis.compare(ctx)

        assert result is not None
        assert result.pairwise_comparisons
        for pair in result.pairwise_comparisons:
            assert pair.direction != "not_testable"
            assert pair.effect_size_interpretation != "not_testable"

    def test_compare_passes_welch_method_to_pairwise_helper(self, monkeypatch):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        condition_a = Condition(
            label="A",
            config_path=Path("/tmp/a.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        condition_b = Condition(
            label="B",
            config_path=Path("/tmp/b.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        aggregated_a = PolymerBridgingAggregatedResult(
            settings_fingerprint=analysis._make_settings_cache_tag(settings),
            n_replicates=2,
            replicates=[1, 2],
            min_ca_distance_angstrom=0.0,
            mean_contacts_per_contacting_oligomer=1.0,
            mean_contacts_sem=0.1,
            multisite_fraction=0.4,
            multisite_fraction_sem=0.05,
            high_valency_fraction=0.1,
            high_valency_fraction_sem=0.02,
            mean_contacts_per_contacting_oligomer_replicates=[0.9, 1.1],
            multisite_fraction_replicates=[0.35, 0.45],
            high_valency_fraction_replicates=[0.09, 0.11],
            valency_probabilities_mean={"1": 0.6, "2": 0.3, "3+": 0.1},
            valency_probabilities_sem={"1": 0.01, "2": 0.01, "3+": 0.01},
            valency_probabilities_per_replicate={
                "1": [0.59, 0.61],
                "2": [0.31, 0.29],
                "3+": [0.1, 0.1],
            },
        )
        aggregated_b = aggregated_a.model_copy(
            update={
                "mean_contacts_per_contacting_oligomer": 1.4,
                "mean_contacts_per_contacting_oligomer_replicates": [1.3, 1.5],
                "multisite_fraction": 0.6,
                "multisite_fraction_replicates": [0.55, 0.65],
                "high_valency_fraction": 0.2,
                "high_valency_fraction_replicates": [0.18, 0.22],
            }
        )

        captured: dict[str, object] = {}

        def _fake_pairwise(metrics_by_condition, control_label, **kwargs):
            captured["ttest_method"] = kwargs.get("ttest_method")
            captured["posthoc_method"] = kwargs.get("posthoc_method")
            captured["fdr_alpha"] = kwargs.get("fdr_alpha")
            return []

        monkeypatch.setattr(
            "polyzymd.analyses.polymer_bridging.pairwise_comparisons",
            _fake_pairwise,
        )

        from polyzymd.analyses.base import ComparisonContext

        ctx = ComparisonContext(
            name="Test",
            conditions=[condition_a, condition_b],
            excluded_conditions=[],
            control_label="A",
            analysis_dirs={"A": Path("/tmp/a"), "B": Path("/tmp/b")},
            results_dir=Path("/tmp/results"),
            equilibration="0ns",
            settings=settings,
            recompute=False,
            aggregated_results={"A": aggregated_a, "B": aggregated_b},
            ttest_method="welch",
            posthoc_method="ttest_bh",
            fdr_alpha=0.1,
        )

        result = analysis.compare(ctx)

        assert result is not None
        assert captured["ttest_method"] == "welch"
        assert captured["posthoc_method"] == "ttest_bh"
        assert captured["fdr_alpha"] == pytest.approx(0.1)
        assert result.ttest_method == "welch"

    def test_compare_honors_tukey_hsd_posthoc_method(self):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        condition_a = Condition(
            label="A",
            config_path=Path("/tmp/a.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        condition_b = Condition(
            label="B",
            config_path=Path("/tmp/b.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        condition_c = Condition(
            label="C",
            config_path=Path("/tmp/c.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        aggregated_a = PolymerBridgingAggregatedResult(
            settings_fingerprint=analysis._make_settings_cache_tag(settings),
            n_replicates=3,
            replicates=[1, 2, 3],
            min_ca_distance_angstrom=0.0,
            mean_contacts_per_contacting_oligomer=1.0,
            mean_contacts_sem=0.1,
            multisite_fraction=0.2,
            multisite_fraction_sem=0.02,
            high_valency_fraction=0.1,
            high_valency_fraction_sem=0.01,
            mean_contacts_per_contacting_oligomer_replicates=[0.9, 1.0, 1.1],
            multisite_fraction_replicates=[0.18, 0.2, 0.22],
            high_valency_fraction_replicates=[0.08, 0.1, 0.12],
            valency_probabilities_mean={"1": 0.7, "2": 0.2, "3+": 0.1},
            valency_probabilities_sem={"1": 0.01, "2": 0.01, "3+": 0.01},
            valency_probabilities_per_replicate={
                "1": [0.7, 0.69, 0.71],
                "2": [0.2, 0.21, 0.19],
                "3+": [0.1, 0.1, 0.1],
            },
        )
        aggregated_b = aggregated_a.model_copy(
            update={
                "mean_contacts_per_contacting_oligomer": 1.4,
                "mean_contacts_per_contacting_oligomer_replicates": [1.3, 1.4, 1.5],
                "multisite_fraction": 0.5,
                "multisite_fraction_replicates": [0.48, 0.5, 0.52],
                "high_valency_fraction": 0.2,
                "high_valency_fraction_replicates": [0.18, 0.2, 0.22],
            }
        )
        aggregated_c = aggregated_a.model_copy(
            update={
                "mean_contacts_per_contacting_oligomer": 1.8,
                "mean_contacts_per_contacting_oligomer_replicates": [1.7, 1.8, 1.9],
                "multisite_fraction": 0.8,
                "multisite_fraction_replicates": [0.78, 0.8, 0.82],
                "high_valency_fraction": 0.4,
                "high_valency_fraction_replicates": [0.38, 0.4, 0.42],
            }
        )

        from polyzymd.analyses.base import ComparisonContext

        ctx = ComparisonContext(
            name="Test",
            conditions=[condition_a, condition_b, condition_c],
            excluded_conditions=[],
            control_label="A",
            analysis_dirs={"A": Path("/tmp/a"), "B": Path("/tmp/b"), "C": Path("/tmp/c")},
            results_dir=Path("/tmp/results"),
            equilibration="0ns",
            settings=settings,
            recompute=False,
            aggregated_results={"A": aggregated_a, "B": aggregated_b, "C": aggregated_c},
            ttest_method="welch",
            posthoc_method="tukey_hsd",
            fdr_alpha=0.05,
        )

        result = analysis.compare(ctx)

        assert result is not None
        assert result.posthoc_method == "tukey_hsd"
        assert result.ttest_method == "welch"
        assert result.pairwise_comparisons
        assert all(pair.posthoc_method == "tukey_hsd" for pair in result.pairwise_comparisons)


class TestConfigCompatibility:
    def test_legacy_analysis_settings_promoted_to_plugins(self, tmp_path):
        yaml_path = tmp_path / "comparison.yaml"
        cond_dir = tmp_path / "cond"
        cond_dir.mkdir()
        (cond_dir / "config.yaml").write_text("name: test\n")
        yaml_path.write_text(
            "name: Legacy\n"
            "defaults:\n  equilibration_time: 0ns\n"
            "conditions:\n"
            f"  - label: A\n    config: {cond_dir / 'config.yaml'}\n    replicates: [1, 2]\n"
            f"  - label: B\n    config: {cond_dir / 'config.yaml'}\n    replicates: [1, 2]\n"
            "analysis_settings:\n"
            "  polymer_bridging:\n"
            "    min_ca_distance_angstrom: 12.0\n"
        )

        from polyzymd.config.comparison import ComparisonConfig

        cfg = ComparisonConfig.from_yaml(yaml_path)

        assert "polymer_bridging" in cfg.plugins.get_enabled_plugins()
        settings = cfg.plugins.get("polymer_bridging")
        assert settings.min_ca_distance_angstrom == pytest.approx(12.0)


class TestSharedPlotSettingsReexport:
    def test_plotsettings_can_be_imported_from_shared(self):
        from polyzymd.analyses.shared import PlotSettings

        assert issubclass(PlotSettings, BaseModel)
