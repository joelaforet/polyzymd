"""Tests for the contacts analysis plugin.

Covers discovery, settings, run_replicate, aggregate, compare (full override),
filter_conditions, binding preference sub-pipeline, plot delegation,
AggregatedResultClass, and per-replicate metric helpers.
"""

from __future__ import annotations

from pathlib import Path
from types import ModuleType, SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


class TestDiscovery:
    """The plugin is found by the automatic discovery system."""

    def test_discovery_by_name(self):
        from polyzymd.analyses.contacts import ContactsAnalysis
        from polyzymd.analyses.discovery import get_analysis

        cls = get_analysis("contacts")
        assert cls is ContactsAnalysis
        assert cls.name == "contacts"

    def test_listed(self):
        from polyzymd.analyses.discovery import list_analyses

        analyses = list_analyses()
        names = list(analyses.keys())
        assert "contacts" in names

    def test_all_names(self):
        from polyzymd.analyses.discovery import list_all_names

        names = list_all_names()
        assert "contacts" in names


# ---------------------------------------------------------------------------
# Class Attributes
# ---------------------------------------------------------------------------


class TestClassAttributes:
    """Verify ClassVar declarations on ContactsAnalysis."""

    def test_name(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        assert ContactsAnalysis.name == "contacts"

    def test_settings_type(self):
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        assert ContactsAnalysis.Settings is ContactsSettings

    def test_aliases(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        assert ContactsAnalysis.aliases == ()

    def test_dependencies(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        assert ContactsAnalysis.dependencies == ()

    def test_min_replicates(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        assert ContactsAnalysis.min_replicates == 2


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class TestSettings:
    """Validate ContactsSettings."""

    def test_defaults(self):
        from polyzymd.analyses.contacts import ContactsSettings

        s = ContactsSettings()
        assert s.polymer_selection == "chainID C"
        assert s.protein_selection == "protein"
        assert s.cutoff == 4.5
        assert s.grouping == "aa_class"
        assert s.compute_residence_times is True
        assert s.compute_binding_preference is False
        assert s.fdr_alpha == 0.05

    def test_custom_values(self):
        from polyzymd.analyses.contacts import ContactsSettings

        s = ContactsSettings(
            polymer_selection="chainID D",
            protein_selection="protein and name CA",
            cutoff=5.0,
            grouping="none",
            compute_binding_preference=True,
            surface_exposure_threshold=0.3,
            fdr_alpha=0.01,
        )
        assert s.polymer_selection == "chainID D"
        assert s.cutoff == 5.0
        assert s.grouping == "none"
        assert s.compute_binding_preference is True
        assert s.surface_exposure_threshold == 0.3
        assert s.fdr_alpha == 0.01

    def test_invalid_grouping(self):
        from polyzymd.analyses.contacts import ContactsSettings

        with pytest.raises(ValueError):
            ContactsSettings(grouping="invalid")

    def test_invalid_fdr_alpha(self):
        from polyzymd.analyses.contacts import ContactsSettings

        with pytest.raises(ValueError):
            ContactsSettings(fdr_alpha=1.5)

    def test_protein_partitions_validation(self):
        from polyzymd.analyses.contacts import ContactsSettings

        # Requires protein_groups to be defined
        with pytest.raises(ValueError):
            ContactsSettings(
                protein_partitions={"part1": ["group1"]},
                # protein_groups not defined
            )

    def test_protein_partitions_undefined_group(self):
        from polyzymd.analyses.contacts import ContactsSettings

        with pytest.raises(ValueError):
            ContactsSettings(
                protein_groups={"existing": [1, 2]},
                protein_partitions={"part1": ["nonexistent"]},
            )

    def test_protein_partitions_overlapping_groups(self):
        from polyzymd.analyses.contacts import ContactsSettings

        with pytest.raises(ValueError):
            ContactsSettings(
                protein_groups={"g1": [1, 2, 3], "g2": [3, 4, 5]},
                protein_partitions={"part1": ["g1", "g2"]},  # resid 3 overlaps
            )

    def test_valid_protein_partitions(self):
        from polyzymd.analyses.contacts import ContactsSettings

        s = ContactsSettings(
            protein_groups={"g1": [1, 2, 3], "g2": [4, 5, 6]},
            protein_partitions={"part1": ["g1", "g2"]},
        )
        assert s.protein_partitions is not None

    def test_serialization_roundtrip(self):
        from polyzymd.analyses.contacts import ContactsSettings

        s = ContactsSettings(
            cutoff=5.0,
            polymer_types=["SBM", "EGM"],
            compute_binding_preference=True,
        )
        d = s.model_dump()
        s2 = ContactsSettings.model_validate(d)
        assert s2.cutoff == 5.0
        assert s2.polymer_types == ["SBM", "EGM"]
        assert s2.compute_binding_preference is True

    def test_settings_fingerprint_changes_with_cutoff(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint

        low = ContactsSettings(cutoff=4.0)
        high = ContactsSettings(cutoff=4.5)

        assert contacts_settings_fingerprint(low) != contacts_settings_fingerprint(high)

    def test_contacts_domain_fingerprint_excludes_binding_preference_fields(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint

        base = ContactsSettings(cutoff=4.5)
        downstream_only = ContactsSettings(
            cutoff=4.5,
            compute_binding_preference=True,
            surface_exposure_threshold=0.5,
            protein_groups={"active": [1, 2]},
            protein_partitions={"site": ["active"]},
            top_residues=25,
        )

        assert contacts_settings_fingerprint(base) == contacts_settings_fingerprint(downstream_only)

    def test_contacts_domain_fingerprint_changes_with_contact_fields(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint

        unfiltered = ContactsSettings(polymer_selection="chainID C")
        filtered = ContactsSettings(polymer_selection="chainID C", polymer_types=["PEG"])

        assert contacts_settings_fingerprint(unfiltered) != contacts_settings_fingerprint(filtered)

    def test_contacts_domain_fingerprint_changes_with_residence_time_toggle(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint

        enabled = ContactsSettings(compute_residence_times=True)
        disabled = ContactsSettings(compute_residence_times=False)

        assert contacts_settings_fingerprint(enabled) != contacts_settings_fingerprint(disabled)

    def test_contacts_domain_candidates_include_legacy_only_when_residence_times_enabled(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint_candidates

        enabled = ContactsSettings(compute_residence_times=True)
        disabled = ContactsSettings(compute_residence_times=False)

        enabled_candidates = contacts_settings_fingerprint_candidates(enabled)
        disabled_candidates = contacts_settings_fingerprint_candidates(disabled)

        assert len(enabled_candidates) > 1
        assert len(disabled_candidates) == 1


# ---------------------------------------------------------------------------
# run_replicate
# ---------------------------------------------------------------------------


def _make_mock_contact_result(replicate: int = 1):
    """Create a mock ContactResult."""
    mock = MagicMock()
    mock.replicate = replicate
    mock.save = MagicMock()
    mock.has_per_residue_statistics = MagicMock(return_value=True)
    return mock


def _make_hashable_sim_config(tmp_path: Path):
    """Build a lightweight config object compatible with ``compute_config_hash``."""

    return SimpleNamespace(
        name="contacts-test",
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


def _make_mdanalysis_exception_modules() -> dict[str, ModuleType]:
    """Build minimal MDAnalysis exception modules for fragment tests."""

    mdanalysis_module = ModuleType("MDAnalysis")
    exceptions_module = ModuleType("MDAnalysis.exceptions")

    class NoDataError(Exception):
        pass

    exceptions_module.NoDataError = NoDataError
    mdanalysis_module.exceptions = exceptions_module
    return {
        "MDAnalysis": mdanalysis_module,
        "MDAnalysis.exceptions": exceptions_module,
    }


class TestRunReplicate:
    """run_replicate delegates to ParallelContactAnalyzer."""

    def test_delegates_to_analyzer(self, tmp_path):
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings()

        mock_sim_config = _make_hashable_sim_config(tmp_path)
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2, 3),
            sim_config=mock_sim_config,
        )
        output_dir = tmp_path / "run_1"
        ctx = ReplicateContext(
            condition=cond,
            replicate=1,
            sim_config=mock_sim_config,
            output_dir=output_dir,
            equilibration="10ns",
            recompute=False,
            settings=settings,
        )

        mock_result = _make_mock_contact_result(replicate=1)

        with (
            patch("polyzymd.analyses.contacts.ParallelContactAnalyzer") as MockAnalyzer,
            patch("polyzymd.analyses.shared.loader.TrajectoryLoader") as MockLoader,
        ):
            mock_universe = MagicMock()
            mock_universe.trajectory.dt = 10.0
            mock_universe.trajectory.__len__.return_value = 2000
            MockLoader.return_value.load_universe.return_value = mock_universe
            MockLoader.return_value.get_timestep.return_value = 10.0

            MockAnalyzer.return_value.run.return_value = mock_result

            result = analysis.run_replicate(ctx, 1)

        assert result is mock_result
        MockAnalyzer.assert_called_once()
        call_kwargs = MockAnalyzer.call_args[1]
        assert call_kwargs["cutoff"] == 4.5

    def test_returns_none_on_missing_trajectory(self, tmp_path):
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.exceptions import ReplicateSkippedError

        analysis = ContactsAnalysis()
        settings = ContactsSettings()

        mock_sim_config = _make_hashable_sim_config(tmp_path)
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim_config,
        )
        ctx = ReplicateContext(
            condition=cond,
            replicate=1,
            sim_config=mock_sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=False,
            settings=settings,
        )

        with patch("polyzymd.analyses.shared.loader.TrajectoryLoader") as MockLoader:
            MockLoader.return_value.load_universe.side_effect = FileNotFoundError("No trajectory")
            with pytest.raises(ReplicateSkippedError, match="No trajectory data found"):
                analysis.run_replicate(ctx, 1)

    def test_cache_filename_includes_settings_fingerprint(self, tmp_path):
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint

        analysis = ContactsAnalysis()
        settings = ContactsSettings(cutoff=4.5)

        mock_sim_config = _make_hashable_sim_config(tmp_path)
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim_config,
        )
        ctx = ReplicateContext(
            condition=cond,
            replicate=1,
            sim_config=mock_sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=False,
            settings=settings,
        )

        mock_result = _make_mock_contact_result(replicate=1)
        expected_fp = contacts_settings_fingerprint(settings)

        with (
            patch("polyzymd.analyses.contacts.ParallelContactAnalyzer") as MockAnalyzer,
            patch("polyzymd.analyses.shared.loader.TrajectoryLoader") as MockLoader,
        ):
            mock_universe = MagicMock()
            mock_universe.trajectory.dt = 10.0
            mock_universe.trajectory.__len__.return_value = 2000
            MockLoader.return_value.load_universe.return_value = mock_universe
            MockLoader.return_value.get_timestep.return_value = 10.0
            MockAnalyzer.return_value.run.return_value = mock_result

            analysis.run_replicate(ctx, 1)

        save_path = mock_result.save.call_args[0][0]
        assert f"_s{expected_fp}_" in str(save_path)

    def test_compute_writes_canonical_and_legacy_paths(self, tmp_path):
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._results import ContactResult

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        sim_config = _make_hashable_sim_config(tmp_path)
        condition = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1,),
            sim_config=sim_config,
        )
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=True,
            settings=settings,
            result_path=tmp_path / "run_1" / "result.json",
        )
        result = ContactResult(
            residue_contacts=[],
            n_frames=10,
            timestep_ps=10.0,
            criteria_label="any_atom_4.5A",
            criteria_cutoff=4.5,
            metadata={},
        )

        with (
            patch("polyzymd.analyses.contacts.ParallelContactAnalyzer") as MockAnalyzer,
            patch("polyzymd.analyses.shared.loader.TrajectoryLoader") as MockLoader,
        ):
            mock_universe = MagicMock()
            mock_universe.trajectory.dt = 10.0
            mock_universe.trajectory.__len__.return_value = 2000
            MockLoader.return_value.load_universe.return_value = mock_universe
            MockLoader.return_value.get_timestep.return_value = 10.0
            MockAnalyzer.return_value.run.return_value = result

            analysis.run_replicate(ctx, 1)

        assert ctx.result_path.exists()
        assert list(ctx.output_dir.glob("contacts_eq10ns_cut4.5_s*_rep1.json"))

    def test_aggregate_writes_canonical_and_legacy_paths(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        sim_config = _make_hashable_sim_config(tmp_path)
        condition = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=sim_config,
        )
        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2),
            output_dir=tmp_path / "aggregated",
            equilibration="10ns",
            settings=settings,
            result_path=tmp_path / "aggregated" / "result.json",
        )
        mock_agg = _make_mock_agg_result(n_replicates=2)

        def _write_mock_agg(path):
            path = Path(path)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("{}")
            return path

        mock_agg.save.side_effect = _write_mock_agg

        with patch(
            "polyzymd.analyses.contacts._aggregator.aggregate_contact_results",
            return_value=mock_agg,
        ):
            analysis.aggregate(ctx, [_make_mock_contact_result(1), _make_mock_contact_result(2)])

        assert ctx.result_path.exists()
        assert list(ctx.output_dir.glob("contacts_aggregated_eq10ns_cut4.5_s*_reps1-2.json"))

    def test_cache_miss_when_cutoff_changes(self, tmp_path):
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._results import ContactResult

        analysis = ContactsAnalysis()
        old_settings = ContactsSettings(cutoff=4.0)
        new_settings = ContactsSettings(cutoff=4.5)
        sim_config = _make_hashable_sim_config(tmp_path)
        condition = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1,),
            sim_config=sim_config,
        )
        legacy_cache = analysis._replicate_sidecar_path(tmp_path / "run_1", old_settings, "10ns", 1)
        ContactResult(
            config_hash="stale",
            residue_contacts=[],
            n_frames=10,
            timestep_ps=10.0,
            criteria_label="any_atom_4.0A",
            criteria_cutoff=4.0,
            metadata={},
        ).save(legacy_cache)
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=False,
            settings=new_settings,
        )
        fresh_result = _make_mock_contact_result(replicate=1)

        with (
            patch("polyzymd.analyses.contacts.ParallelContactAnalyzer") as MockAnalyzer,
            patch("polyzymd.analyses.shared.loader.TrajectoryLoader") as MockLoader,
        ):
            mock_universe = MagicMock()
            mock_universe.trajectory.dt = 10.0
            mock_universe.trajectory.__len__.return_value = 2000
            MockLoader.return_value.load_universe.return_value = mock_universe
            MockLoader.return_value.get_timestep.return_value = 10.0
            MockAnalyzer.return_value.run.return_value = fresh_result

            result = analysis.run_replicate(ctx, 1)

        assert result is fresh_result
        MockAnalyzer.assert_called_once()

    def test_canonical_cache_miss_when_cutoff_changes(self, tmp_path):
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint
        from polyzymd.analyses.contacts._results import ContactResult

        analysis = ContactsAnalysis()
        old_settings = ContactsSettings(cutoff=4.0)
        new_settings = ContactsSettings(cutoff=4.5)
        sim_config = _make_hashable_sim_config(tmp_path)
        condition = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1,),
            sim_config=sim_config,
        )
        canonical = tmp_path / "run_1" / "result.json"
        ContactResult(
            config_hash="stale",
            residue_contacts=[],
            n_frames=10,
            timestep_ps=10.0,
            criteria_label="any_atom_4.0A",
            criteria_cutoff=4.0,
            metadata={"settings_fingerprint": contacts_settings_fingerprint(old_settings)},
        ).save(canonical)
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=False,
            settings=new_settings,
            result_path=canonical,
        )
        fresh_result = _make_mock_contact_result(replicate=1)

        with (
            patch("polyzymd.analyses.contacts.ParallelContactAnalyzer") as MockAnalyzer,
            patch("polyzymd.analyses.shared.loader.TrajectoryLoader") as MockLoader,
        ):
            mock_universe = MagicMock()
            mock_universe.trajectory.dt = 10.0
            mock_universe.trajectory.__len__.return_value = 2000
            MockLoader.return_value.load_universe.return_value = mock_universe
            MockLoader.return_value.get_timestep.return_value = 10.0
            MockAnalyzer.return_value.run.return_value = fresh_result

            result = analysis.run_replicate(ctx, 1)

        assert result is fresh_result
        MockAnalyzer.assert_called_once()

    def test_cached_replicate_id_mismatch_is_rejected(self, tmp_path):
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._results import ContactResult

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        sim_config = _make_hashable_sim_config(tmp_path)
        condition = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1,),
            sim_config=sim_config,
        )
        cache = analysis._replicate_sidecar_path(tmp_path / "run_1", settings, "10ns", 1)
        stale = ContactResult(
            residue_contacts=[],
            n_frames=10,
            timestep_ps=10.0,
            criteria_label="any_atom_4.5A",
            criteria_cutoff=4.5,
            replicate=2,
            equilibration_time=10.0,
            equilibration_unit="ns",
        )
        stale.save(cache)
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=False,
            settings=settings,
        )
        fresh_result = _make_mock_contact_result(replicate=1)

        with (
            patch("polyzymd.analyses.contacts.ParallelContactAnalyzer") as MockAnalyzer,
            patch("polyzymd.analyses.shared.loader.TrajectoryLoader") as MockLoader,
        ):
            mock_universe = MagicMock()
            mock_universe.trajectory.dt = 10.0
            mock_universe.trajectory.__len__.return_value = 2000
            MockLoader.return_value.load_universe.return_value = mock_universe
            MockLoader.return_value.get_timestep.return_value = 10.0
            MockAnalyzer.return_value.run.return_value = fresh_result

            result = analysis.run_replicate(ctx, 1)

        assert result is fresh_result
        MockAnalyzer.assert_called_once()

    def test_canonical_cache_miss_when_equilibration_changes(self, tmp_path):
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint
        from polyzymd.analyses.contacts._results import ContactResult

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        sim_config = _make_hashable_sim_config(tmp_path)
        condition = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1,),
            sim_config=sim_config,
        )
        canonical = tmp_path / "run_1" / "result.json"
        ContactResult(
            config_hash="stale",
            residue_contacts=[],
            n_frames=10,
            timestep_ps=10.0,
            criteria_label="any_atom_4.5A",
            criteria_cutoff=4.5,
            equilibration_time=0.0,
            equilibration_unit="ns",
            metadata={"settings_fingerprint": contacts_settings_fingerprint(settings)},
        ).save(canonical)
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=False,
            settings=settings,
            result_path=canonical,
        )
        fresh_result = _make_mock_contact_result(replicate=1)

        with (
            patch("polyzymd.analyses.contacts.ParallelContactAnalyzer") as MockAnalyzer,
            patch("polyzymd.analyses.shared.loader.TrajectoryLoader") as MockLoader,
        ):
            mock_universe = MagicMock()
            mock_universe.trajectory.dt = 10.0
            mock_universe.trajectory.__len__.return_value = 2000
            MockLoader.return_value.load_universe.return_value = mock_universe
            MockLoader.return_value.get_timestep.return_value = 10.0
            MockAnalyzer.return_value.run.return_value = fresh_result

            result = analysis.run_replicate(ctx, 1)

        assert result is fresh_result
        MockAnalyzer.assert_called_once()

    def test_compute_prefers_sidecar_over_canonical_cache(self, tmp_path):
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint
        from polyzymd.analyses.contacts._results import ContactResult

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        sim_config = _make_hashable_sim_config(tmp_path)
        condition = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1,),
            sim_config=sim_config,
        )
        output_dir = tmp_path / "run_1"
        canonical = output_dir / "result.json"
        sidecar = analysis._replicate_sidecar_path(output_dir, settings, "10ns", 1)
        ContactResult(
            replicate=1,
            residue_contacts=[],
            n_frames=10,
            timestep_ps=10.0,
            criteria_label="canonical",
            criteria_cutoff=4.5,
            equilibration_time=10.0,
            equilibration_unit="ns",
        ).save(canonical)
        ContactResult(
            replicate=1,
            residue_contacts=[],
            n_frames=10,
            timestep_ps=10.0,
            criteria_label="sidecar",
            criteria_cutoff=4.5,
            equilibration_time=10.0,
            equilibration_unit="ns",
            metadata={"settings_fingerprint": contacts_settings_fingerprint(settings)},
        ).save(sidecar)
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=sim_config,
            output_dir=output_dir,
            equilibration="10ns",
            recompute=False,
            settings=settings,
            result_path=canonical,
        )

        result = analysis.run_replicate(ctx, 1)

        assert result.criteria_label == "sidecar"

    def test_fingerprintless_canonical_not_copied_to_sidecar(self, tmp_path):
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._results import ContactResult

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        sim_config = _make_hashable_sim_config(tmp_path)
        output_dir = tmp_path / "run_1"
        canonical = output_dir / "result.json"
        sidecar = analysis._replicate_sidecar_path(output_dir, settings, "10ns", 1)
        ContactResult(
            replicate=1,
            residue_contacts=[],
            n_frames=10,
            timestep_ps=10.0,
            criteria_label="canonical",
            criteria_cutoff=4.5,
            equilibration_time=10.0,
            equilibration_unit="ns",
        ).save(canonical)
        condition = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1,),
            sim_config=sim_config,
        )
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=sim_config,
            output_dir=output_dir,
            equilibration="10ns",
            recompute=False,
            settings=settings,
            result_path=canonical,
        )

        with patch(
            "polyzymd.analyses.base.Analysis._run_replicate_default",
            side_effect=RuntimeError("should recompute"),
        ):
            with pytest.raises(RuntimeError, match="should recompute"):
                analysis.run_replicate(ctx, 1)

        assert not sidecar.exists()

    def test_filename_only_sidecar_fingerprint_rejected(self, tmp_path):
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._results import ContactResult

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        sim_config = _make_hashable_sim_config(tmp_path)
        output_dir = tmp_path / "run_1"
        sidecar = analysis._replicate_sidecar_path(output_dir, settings, "10ns", 1)
        ContactResult(
            replicate=1,
            residue_contacts=[],
            n_frames=10,
            timestep_ps=10.0,
            criteria_label="sidecar",
            criteria_cutoff=4.5,
            equilibration_time=10.0,
            equilibration_unit="ns",
        ).save(sidecar)
        condition = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1,),
            sim_config=sim_config,
        )
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=sim_config,
            output_dir=output_dir,
            equilibration="10ns",
            recompute=False,
            settings=settings,
        )

        with patch(
            "polyzymd.analyses.base.Analysis._run_replicate_default",
            side_effect=RuntimeError("should recompute"),
        ):
            with pytest.raises(RuntimeError, match="should recompute"):
                analysis.run_replicate(ctx, 1)

    def test_disabled_residence_times_rejects_fingerprintless_legacy_cache(self):
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._results import ContactResult

        settings = ContactsSettings(compute_residence_times=False)
        legacy = ContactResult(
            replicate=1,
            residue_contacts=[],
            n_frames=10,
            timestep_ps=10.0,
            criteria_label="any_atom_4.5A",
            criteria_cutoff=4.5,
            selection_string="protein : chainID C",
            equilibration_time=10.0,
            equilibration_unit="ns",
            metadata={"grouping": "aa_class"},
        )

        assert not ContactsAnalysis._cache_matches_contacts_settings(legacy, settings)

    def test_enabled_residence_times_accepts_fingerprintless_legacy_cache(self):
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._results import ContactResult

        settings = ContactsSettings(compute_residence_times=True)
        legacy = ContactResult(
            replicate=1,
            residue_contacts=[],
            n_frames=10,
            timestep_ps=10.0,
            criteria_label="any_atom_4.5A",
            criteria_cutoff=4.5,
            selection_string="protein : chainID C",
            equilibration_time=10.0,
            equilibration_unit="ns",
            metadata={"grouping": "aa_class"},
        )

        assert ContactsAnalysis._cache_matches_contacts_settings(legacy, settings)
        assert legacy.metadata["compute_residence_times"] is True


class TestParallelContactAnalyzerResidueIdentity:
    """Ensure residue identity preserves chain-separated duplicate resid values."""

    def test_build_atom_to_residue_map_distinguishes_duplicate_resids_by_chain(self):
        from polyzymd.analyses.contacts import ParallelContactAnalyzer

        class _Residue:
            def __init__(self, resid: int, resname: str, chain_id: str):
                self.resid = resid
                self.resname = resname
                self.chainID = chain_id

        class _Atom:
            def __init__(self, residue):
                self.residue = residue

        residue_a = _Residue(1, "SBM", "C")
        residue_b = _Residue(1, "SBM", "D")
        residues = [residue_a, residue_b]
        atoms = [_Atom(residue_a), _Atom(residue_b)]

        atom_to_res = ParallelContactAnalyzer._build_atom_to_residue_map(atoms, residues)

        assert atom_to_res.tolist() == [0, 1]

    def test_run_keeps_duplicate_resids_distinct_across_chains(self):
        import numpy as np

        from polyzymd.analyses.contacts import ParallelContactAnalyzer

        class _Residue:
            def __init__(self, resid: int, resname: str, chain_id: str, ix: int):
                self.resid = resid
                self.resname = resname
                self.chainID = chain_id
                self.ix = ix

        class _AtomGroup(list):
            @property
            def atoms(self):
                return self

            @property
            def fragments(self):
                return [self]

        class _Fragment:
            def __init__(self, residues):
                self.residues = residues

        class _QueryAtomGroup(_AtomGroup):
            @property
            def fragments(self):
                return [_Fragment([query_residues[0]]), _Fragment([query_residues[1]])]

        class _ResidueGroup(list):
            def __init__(self, residues, atoms):
                super().__init__(residues)
                self.atoms = atoms

        class _TargetAtom:
            def __init__(self, residue):
                self.residue = residue

        class _QueryAtom:
            def __init__(self, residue):
                self.residue = residue

        protein_residue = _Residue(10, "ALA", "A", ix=100)
        query_residues = [
            _Residue(1, "SBM", "C", ix=200),
            _Residue(1, "SBM", "D", ix=201),
        ]

        target_atoms = _AtomGroup([_TargetAtom(protein_residue)])
        query_atoms = _QueryAtomGroup(
            [_QueryAtom(query_residues[0]), _QueryAtom(query_residues[1])]
        )

        selector_result_target = MagicMock()
        selector_result_target.atoms = target_atoms
        selector_result_target.residues = [protein_residue]
        selector_result_query = MagicMock()
        selector_result_query.atoms = query_atoms
        selector_result_query.residues = _ResidueGroup(query_residues, query_atoms)

        target_selector = MagicMock()
        target_selector.select.return_value = selector_result_target
        target_selector.label = "protein"
        query_selector = MagicMock()
        query_selector.select.return_value = selector_result_query
        query_selector.label = "polymer"

        analyzer = ParallelContactAnalyzer(
            target_selector=target_selector,
            query_selector=query_selector,
            cutoff=4.5,
        )

        contact_matrix = np.array([[[1, 1]]], dtype=np.uint8)
        fake_analysis = MagicMock()
        fake_analysis.results.contact_matrix = contact_matrix
        fake_analysis.run.return_value = None

        universe = MagicMock()
        universe.trajectory.dt = 10.0

        with patch(
            "polyzymd.analyses.contacts._runner._get_contact_analysis_base_cls"
        ) as mock_factory:
            mock_factory.return_value.return_value = fake_analysis
            result = analyzer.run(universe, start=0, stop=1, step=1)

        assert result.n_frames == 1
        residue_contacts = result.residue_contacts[0].segment_contacts
        assert len(residue_contacts) == 2
        assert sorted(sc.polymer_chain_idx for sc in residue_contacts) == [0, 1]
        assert [sc.polymer_resid for sc in residue_contacts] == [1, 1]

    def test_fragment_lookup_falls_back_without_bond_topology(self):
        """No-bond topologies should be treated as one polymer fragment."""
        from polyzymd.analyses.contacts._runner import _fragments_or_single

        modules = _make_mdanalysis_exception_modules()
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
        from polyzymd.analyses.contacts._runner import _fragments_or_single

        class CustomFragmentError(RuntimeError):
            pass

        class _AtomGroup:
            @property
            def fragments(self):
                raise CustomFragmentError("unexpected fragment failure")

        with patch.dict("sys.modules", _make_mdanalysis_exception_modules()):
            with pytest.raises(CustomFragmentError, match="unexpected fragment failure"):
                _fragments_or_single(_AtomGroup(), context="test")


class TestContactsRunnerGrouping:
    """Regression tests for settings grouping across the runner seam."""

    def _make_context(self, tmp_path, settings):
        """Build a replicate context for runner construction tests."""

        from polyzymd.analyses.base import Condition, ReplicateContext

        condition = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1,),
            sim_config=MagicMock(),
        )
        return ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="0ns",
            recompute=False,
            settings=settings,
        )

    def test_build_runner_threads_none_grouping(self, tmp_path):
        """The ``none`` grouping should preserve residue-name labels."""

        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        ctx = self._make_context(tmp_path, ContactsSettings(grouping="none"))

        runner = analysis.build_runner(ctx, 1, MagicMock(), MagicMock())

        assert runner.grouping.classify("ALA") == "ALA"

    def test_build_runner_threads_secondary_structure_grouping(self, tmp_path):
        """Secondary-structure grouping should use residue annotations."""

        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        ctx = self._make_context(
            tmp_path,
            ContactsSettings(grouping="secondary_structure"),
        )

        runner = analysis.build_runner(ctx, 1, MagicMock(), MagicMock())

        residue = SimpleNamespace(resname="ALA", secondary_structure="E")
        assert runner.grouping.classify_residue(residue) == "sheet"

    def test_build_runner_threads_aa_class_grouping(self, tmp_path):
        """AA-class grouping should retain the default biochemical classes."""

        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        ctx = self._make_context(tmp_path, ContactsSettings(grouping="aa_class"))

        runner = analysis.build_runner(ctx, 1, MagicMock(), MagicMock())

        assert runner.grouping.classify("ALA") == "nonpolar"

    def test_build_runner_applies_polymer_types_to_query_selection(self, tmp_path):
        """polymer_types should constrain the runner-backed polymer selection."""

        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings(polymer_selection="chainID C", polymer_types=["SBM", "EGM"])
        ctx = self._make_context(tmp_path, settings)

        runner = analysis.build_runner(ctx, 1, MagicMock(), MagicMock())

        assert runner.query_selector.selection == "(chainID C) and (resname SBM EGM)"

    def test_settings_fingerprint_changes_with_polymer_types(self):
        """polymer_types should participate in contacts cache identity."""

        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.shared.config_hash import settings_fingerprint

        unfiltered = ContactsSettings(polymer_selection="chainID C")
        filtered = ContactsSettings(polymer_selection="chainID C", polymer_types=["SBM"])

        assert settings_fingerprint(unfiltered) != settings_fingerprint(filtered)


# ---------------------------------------------------------------------------
# aggregate
# ---------------------------------------------------------------------------


def _make_mock_agg_result(n_replicates: int = 3, n_residues: int = 5):
    """Create a mock AggregatedContactResult."""
    mock = MagicMock()
    mock.n_replicates = n_replicates
    mock.replicates = list(range(1, n_replicates + 1))
    mock.n_residues = n_residues
    mock.coverage_mean = 0.8
    mock.coverage_sem = 0.05
    mock.mean_contact_fraction = 0.3
    mock.mean_contact_fraction_sem = 0.02
    mock.residence_time_by_polymer_type = {"SBM": (10.0, 1.0)}
    mock.metadata = {}
    mock.save = MagicMock()

    residue_stats = []
    for i in range(n_residues):
        rs = MagicMock()
        rs.protein_resid = i + 1
        # Per-replicate contact fractions
        fracs = [0.2 + i * 0.05 + rep * 0.01 for rep in range(n_replicates)]
        rs.contact_fraction_per_replicate = fracs
        residue_stats.append(rs)
    mock.residue_stats = residue_stats

    return mock


def _make_valid_agg_result(settings, n_replicates: int = 2, n_residues: int = 2):
    """Create a valid aggregated contacts result for compare regression tests.

    Parameters
    ----------
    settings : ContactsSettings
        Contacts settings used to stamp cache identity metadata.
    n_replicates : int, optional
        Number of replicate vectors to include, by default 2.
    n_residues : int, optional
        Number of residue summary rows to include, by default 2.

    Returns
    -------
    AggregatedContactResult
        Aggregate result with window, replicate, and settings metadata.
    """
    from polyzymd.analyses.contacts._aggregator import (
        AggregatedContactResult,
        AggregatedResidueStats,
    )
    from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint

    residue_stats = []
    for residue_index in range(n_residues):
        fractions = [0.25 + 0.05 * residue_index + 0.05 * rep for rep in range(n_replicates)]
        residue_stats.append(
            AggregatedResidueStats(
                protein_resid=residue_index + 1,
                protein_resname="ALA",
                contact_fraction_mean=sum(fractions) / len(fractions),
                contact_fraction_sem=0.01,
                contact_fraction_per_replicate=fractions,
            )
        )

    return AggregatedContactResult(
        n_replicates=n_replicates,
        replicates=list(range(1, n_replicates + 1)),
        residue_stats=residue_stats,
        total_frames_per_replicate=[10] * n_replicates,
        criteria_cutoff=float(settings.cutoff),
        coverage_mean=1.0,
        coverage_sem=0.0,
        mean_contact_fraction=0.325,
        mean_contact_fraction_sem=0.01,
        equilibration_time=10.0,
        equilibration_unit="ns",
        metadata={
            "settings_fingerprint": contacts_settings_fingerprint(settings),
            "compute_residence_times": bool(settings.compute_residence_times),
            "residence_times_computed": bool(settings.compute_residence_times),
            "replicates": list(range(1, n_replicates + 1)),
        },
    )


class TestAggregate:
    """aggregate delegates to aggregate_contact_results."""

    def test_aggregate_delegates(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        output_dir = tmp_path / "aggregated"

        mock_sim_config = _make_hashable_sim_config(tmp_path)
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2, 3),
            sim_config=mock_sim_config,
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1, 2, 3),
            output_dir=output_dir,
            equilibration="10ns",
            settings=settings,
        )

        mock_results = [_make_mock_contact_result(i) for i in range(1, 4)]
        mock_agg = _make_mock_agg_result(n_replicates=3)

        with patch(
            "polyzymd.analyses.contacts._aggregator.aggregate_contact_results",
            return_value=mock_agg,
        ) as mock_fn:
            result = analysis.aggregate(ctx, mock_results)

        mock_fn.assert_called_once()
        assert mock_fn.call_args.kwargs["compute_residence_times"] is True
        assert result is mock_agg
        assert mock_agg.metadata["compute_residence_times"] is True
        assert mock_agg.metadata["residence_times_computed"] is True
        assert mock_agg.metadata["settings_fingerprint_domain"] == "contacts-v2"
        mock_agg.save.assert_called_once()

    def test_aggregate_passes_disabled_residence_time_setting(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings(compute_residence_times=False)
        mock_sim_config = _make_hashable_sim_config(tmp_path)
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=mock_sim_config,
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1, 2),
            output_dir=tmp_path / "aggregated",
            equilibration="10ns",
            settings=settings,
        )
        mock_agg = _make_mock_agg_result(n_replicates=2)

        with patch(
            "polyzymd.analyses.contacts._aggregator.aggregate_contact_results",
            return_value=mock_agg,
        ) as mock_fn:
            analysis.aggregate(ctx, [_make_mock_contact_result(1), _make_mock_contact_result(2)])

        assert mock_fn.call_args.kwargs["compute_residence_times"] is False
        assert mock_agg.metadata["compute_residence_times"] is False
        assert mock_agg.metadata["residence_times_computed"] is False

    def test_aggregate_rejects_replicate_id_mismatch(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        mock_sim_config = _make_hashable_sim_config(tmp_path)
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=mock_sim_config,
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1, 2),
            output_dir=tmp_path / "aggregated",
            equilibration="10ns",
            settings=settings,
        )

        with pytest.raises(ValueError, match="replicate IDs do not match"):
            analysis.aggregate(ctx, [_make_mock_contact_result(1), _make_mock_contact_result(3)])

    def test_aggregate_saves_file(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        output_dir = tmp_path / "aggregated"

        mock_sim_config = _make_hashable_sim_config(tmp_path)
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 3),
            sim_config=mock_sim_config,
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1, 3),
            output_dir=output_dir,
            equilibration="10ns",
            settings=settings,
        )

        mock_results = [_make_mock_contact_result(1), _make_mock_contact_result(3)]
        mock_agg = _make_mock_agg_result(n_replicates=2)

        with patch(
            "polyzymd.analyses.contacts._aggregator.aggregate_contact_results",
            return_value=mock_agg,
        ):
            analysis.aggregate(ctx, mock_results)

        # Confirm the sidecar path encodes the exact replicate set
        save_call = mock_agg.save.call_args[0][0]
        assert "reps1_3" in str(save_call)

    def test_aggregate_rejects_stale_config_hash_cache(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._aggregator import AggregatedContactResult
        from polyzymd.analyses.shared.config_hash import settings_fingerprint

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        sim_config = _make_hashable_sim_config(tmp_path)
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=sim_config,
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1, 2),
            output_dir=tmp_path / "aggregated",
            equilibration="10ns",
            settings=settings,
        )
        stale = AggregatedContactResult(
            config_hash="stale-config",
            equilibration_time=10.0,
            equilibration_unit="ns",
            metadata={"settings_fingerprint": settings_fingerprint(settings)},
            n_replicates=2,
            criteria_cutoff=4.5,
            coverage_mean=9.9,
            mean_contact_fraction=9.9,
        )
        cache_path = analysis._aggregated_sidecar_path(
            ctx.output_dir,
            settings,
            ctx.equilibration,
            ctx.replicates,
        )
        stale.save(cache_path)
        fresh = _make_mock_agg_result(n_replicates=2)

        with patch(
            "polyzymd.analyses.contacts._aggregator.aggregate_contact_results",
            return_value=fresh,
        ):
            result = analysis.aggregate(
                ctx, [_make_mock_contact_result(1), _make_mock_contact_result(2)]
            )

        assert result is fresh

    def test_aggregate_prefers_sidecar_over_canonical_cache(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._aggregator import AggregatedContactResult
        from polyzymd.analyses.shared.config_hash import settings_fingerprint

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        sim_config = _make_hashable_sim_config(tmp_path)
        condition = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=sim_config,
        )
        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2),
            output_dir=tmp_path / "aggregated",
            equilibration="10ns",
            settings=settings,
            result_path=tmp_path / "aggregated" / "result.json",
        )
        metadata = {"settings_fingerprint": settings_fingerprint(settings), "replicates": [1, 2]}
        canonical = AggregatedContactResult(
            n_replicates=2,
            replicates=[1, 2],
            total_frames_per_replicate=[10, 10],
            criteria_cutoff=4.5,
            coverage_mean=9.0,
            mean_contact_fraction=9.0,
            equilibration_time=10.0,
            equilibration_unit="ns",
            metadata=metadata,
        )
        sidecar = canonical.model_copy(update={"coverage_mean": 0.2})
        canonical.save(ctx.result_path)
        sidecar.save(
            analysis._aggregated_sidecar_path(
                ctx.output_dir,
                settings,
                ctx.equilibration,
                ctx.replicates,
            )
        )

        result = analysis.aggregate(ctx, [])

        assert result.coverage_mean == pytest.approx(0.2)

    def test_aggregate_rejects_replicate_set_mismatch_cache(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._aggregator import AggregatedContactResult
        from polyzymd.analyses.shared.config_hash import settings_fingerprint

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        sim_config = _make_hashable_sim_config(tmp_path)
        condition = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2, 3),
            sim_config=sim_config,
        )
        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2, 3),
            output_dir=tmp_path / "aggregated",
            equilibration="10ns",
            settings=settings,
            result_path=tmp_path / "aggregated" / "result.json",
        )
        AggregatedContactResult(
            n_replicates=2,
            replicates=[1, 2],
            criteria_cutoff=4.5,
            coverage_mean=9.0,
            mean_contact_fraction=9.0,
            equilibration_time=10.0,
            equilibration_unit="ns",
            metadata={"settings_fingerprint": settings_fingerprint(settings), "replicates": [1, 2]},
        ).save(ctx.result_path)
        fresh = _make_mock_agg_result(n_replicates=3)

        with patch(
            "polyzymd.analyses.contacts._aggregator.aggregate_contact_results",
            return_value=fresh,
        ):
            result = analysis.aggregate(
                ctx,
                [
                    _make_mock_contact_result(1),
                    _make_mock_contact_result(2),
                    _make_mock_contact_result(3),
                ],
            )

        assert result is fresh
        assert fresh.replicates == [1, 2, 3]


class TestAggregatePolymerTypeCoverage:
    """Regression tests for polymer-type per-replicate aggregation."""

    @staticmethod
    def _make_residence_time_results():
        """Create minimal contact results with PEG events in two replicates."""

        from polyzymd.analyses.contacts._results import (
            ContactEvent,
            ContactResult,
            PolymerSegmentContacts,
            ResidueContactData,
        )

        return [
            ContactResult(
                residue_contacts=[
                    ResidueContactData(
                        protein_resid=1,
                        protein_resname="ALA",
                        protein_group="nonpolar",
                        segment_contacts=[
                            PolymerSegmentContacts(
                                polymer_resname="PEG",
                                polymer_resid=1,
                                polymer_chain_idx=0,
                                events=[ContactEvent(start_frame=0, duration=5)],
                            )
                        ],
                        statistical_inefficiency=1.0,
                        n_effective=10.0,
                    )
                ],
                n_frames=10,
                criteria_label="distance",
                criteria_cutoff=4.5,
                replicate=1,
            ),
            ContactResult(
                residue_contacts=[
                    ResidueContactData(
                        protein_resid=1,
                        protein_resname="ALA",
                        protein_group="nonpolar",
                        segment_contacts=[
                            PolymerSegmentContacts(
                                polymer_resname="PEG",
                                polymer_resid=1,
                                polymer_chain_idx=0,
                                events=[ContactEvent(start_frame=0, duration=3)],
                            )
                        ],
                        statistical_inefficiency=1.0,
                        n_effective=10.0,
                    )
                ],
                n_frames=10,
                criteria_label="distance",
                criteria_cutoff=4.5,
                replicate=2,
            ),
        ]

    def test_includes_zero_contact_replicates_in_polymer_type_vectors(self):
        from polyzymd.analyses.contacts._aggregator import aggregate_contact_results
        from polyzymd.analyses.contacts._results import (
            ContactEvent,
            ContactResult,
            PolymerSegmentContacts,
            ResidueContactData,
        )

        rep1 = ContactResult(
            residue_contacts=[
                ResidueContactData(
                    protein_resid=1,
                    protein_resname="ALA",
                    protein_group="nonpolar",
                    segment_contacts=[
                        PolymerSegmentContacts(
                            polymer_resname="PEG",
                            polymer_resid=1,
                            polymer_chain_idx=0,
                            events=[ContactEvent(start_frame=0, duration=5)],
                        )
                    ],
                    statistical_inefficiency=1.0,
                    n_effective=10.0,
                )
            ],
            n_frames=10,
            criteria_label="distance",
            criteria_cutoff=4.5,
            replicate=1,
        )
        rep2 = ContactResult(
            residue_contacts=[
                ResidueContactData(
                    protein_resid=1,
                    protein_resname="ALA",
                    protein_group="nonpolar",
                    segment_contacts=[
                        PolymerSegmentContacts(
                            polymer_resname="PEG",
                            polymer_resid=1,
                            polymer_chain_idx=0,
                            events=[ContactEvent(start_frame=0, duration=3)],
                        )
                    ],
                    statistical_inefficiency=1.0,
                    n_effective=10.0,
                )
            ],
            n_frames=10,
            criteria_label="distance",
            criteria_cutoff=4.5,
            replicate=2,
        )
        rep3 = ContactResult(
            residue_contacts=[
                ResidueContactData(
                    protein_resid=1,
                    protein_resname="ALA",
                    protein_group="nonpolar",
                    segment_contacts=[],
                    statistical_inefficiency=1.0,
                    n_effective=10.0,
                )
            ],
            n_frames=10,
            criteria_label="distance",
            criteria_cutoff=4.5,
            replicate=3,
        )

        aggregated = aggregate_contact_results([rep1, rep2, rep3])
        residue = aggregated.residue_stats[0]

        assert residue.by_polymer_type_per_replicate["PEG"] == [0.5, 0.3, 0.0]
        assert residue.by_polymer_type["PEG"][0] == pytest.approx((0.5 + 0.3 + 0.0) / 3.0)

    def test_sparse_residence_time_vectors_record_replicate_identity(self):
        from polyzymd.analyses.contacts import ContactsAnalysis
        from polyzymd.analyses.contacts._aggregator import aggregate_contact_results
        from polyzymd.analyses.contacts._results import (
            ContactEvent,
            ContactResult,
            PolymerSegmentContacts,
            ResidueContactData,
        )

        rep1 = ContactResult(
            residue_contacts=[
                ResidueContactData(
                    protein_resid=1,
                    protein_resname="ALA",
                    protein_group="nonpolar",
                    segment_contacts=[
                        PolymerSegmentContacts(
                            polymer_resname="PEG",
                            polymer_resid=1,
                            polymer_chain_idx=0,
                            events=[ContactEvent(start_frame=0, duration=5)],
                        )
                    ],
                    statistical_inefficiency=1.0,
                    n_effective=10.0,
                )
            ],
            n_frames=10,
            criteria_label="distance",
            criteria_cutoff=4.5,
            replicate=1,
        )
        rep2 = ContactResult(
            residue_contacts=[
                ResidueContactData(
                    protein_resid=1,
                    protein_resname="ALA",
                    protein_group="nonpolar",
                    segment_contacts=[],
                    statistical_inefficiency=1.0,
                    n_effective=10.0,
                )
            ],
            n_frames=10,
            criteria_label="distance",
            criteria_cutoff=4.5,
            replicate=2,
        )

        aggregated = aggregate_contact_results([rep1, rep2])
        residue = aggregated.residue_stats[0]

        assert residue.residence_time_by_polymer_type_per_replicate["PEG"] == [5.0]
        assert residue.residence_time_by_polymer_type_replicates["PEG"] == [1]
        assert ContactsAnalysis()._cache_matches_replicates(aggregated, (1, 2))

    def test_residence_time_summaries_are_computed_by_default(self):
        from polyzymd.analyses.contacts._aggregator import aggregate_contact_results

        aggregated = aggregate_contact_results(self._make_residence_time_results())
        residue = aggregated.residue_stats[0]

        assert aggregated.residence_time_by_polymer_type["PEG"][0] == pytest.approx(4.0)
        assert aggregated.residence_time_by_polymer_type_replicates["PEG"] == [1, 2]
        assert residue.residence_time_by_polymer_type["PEG"][0] == pytest.approx(4.0)
        assert residue.residence_time_by_polymer_type_replicates["PEG"] == [1, 2]

    def test_residence_time_summaries_empty_when_disabled(self):
        from polyzymd.analyses.contacts._aggregator import aggregate_contact_results

        aggregated = aggregate_contact_results(
            self._make_residence_time_results(),
            compute_residence_times=False,
        )
        residue = aggregated.residue_stats[0]

        assert aggregated.residence_time_by_polymer_type == {}
        assert aggregated.residence_time_by_polymer_type_replicates == {}
        assert residue.residence_time_by_polymer_type == {}
        assert residue.residence_time_by_polymer_type_per_replicate == {}
        assert residue.residence_time_by_polymer_type_replicates == {}

    def test_sparse_residence_time_mismatched_identity_is_rejected(self):
        from polyzymd.analyses.contacts import ContactsAnalysis
        from polyzymd.analyses.contacts._aggregator import (
            AggregatedContactResult,
            AggregatedResidueStats,
        )

        aggregated = AggregatedContactResult(
            n_replicates=2,
            replicates=[1, 2],
            total_frames_per_replicate=[10, 10],
            residue_stats=[
                AggregatedResidueStats(
                    protein_resid=1,
                    protein_resname="ALA",
                    contact_fraction_per_replicate=[0.5, 0.0],
                    residence_time_by_polymer_type_per_replicate={"PEG": [5.0, 7.0]},
                    residence_time_by_polymer_type_replicates={"PEG": [1]},
                )
            ],
        )

        assert not ContactsAnalysis()._cache_matches_replicates(aggregated, (1, 2))

    def test_sparse_residence_time_expands_to_aggregate_replicate_order(self):
        from polyzymd.analyses.contacts._aggregator import (
            AggregatedContactResult,
            AggregatedResidueStats,
        )

        aggregated = AggregatedContactResult(
            n_replicates=2,
            replicates=[1, 3],
            total_frames_per_replicate=[10, 10],
            timestep_ps=1000.0,
            residue_stats=[
                AggregatedResidueStats(
                    protein_resid=1,
                    protein_resname="ALA",
                    protein_group="nonpolar",
                    contact_fraction_per_replicate=[0.0, 0.5],
                    residence_time_by_polymer_type_per_replicate={"PEG": [7.0]},
                    residence_time_by_polymer_type_replicates={"PEG": [3]},
                )
            ],
        )

        assert aggregated.subset_residence_time_per_replicate(
            [1], polymer_type="PEG", units="frames"
        ) == [0.0, 7.0]
        assert aggregated.group_residence_time_per_replicate(
            polymer_type="PEG", units="frames"
        ) == {"nonpolar": [0.0, 7.0]}


# ---------------------------------------------------------------------------
# Per-replicate metric computation
# ---------------------------------------------------------------------------


class TestPerReplicateMetrics:
    """Test _compute_coverage_per_replicate and _compute_contact_fraction_per_replicate."""

    def test_coverage_per_replicate(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        agg = _make_mock_agg_result(n_replicates=3, n_residues=5)
        # By default all residues have >0 fractions, so coverage should be 1.0
        coverages = ContactsAnalysis._compute_coverage_per_replicate(agg)
        assert len(coverages) == 3
        for c in coverages:
            assert c == 1.0  # all residues have >0 contact fraction

    def test_coverage_with_zero_fractions(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        agg = MagicMock()
        agg.n_replicates = 2
        agg.n_residues = 3

        # Residue 1: contacted in both reps
        # Residue 2: contacted only in rep 0
        # Residue 3: never contacted
        rs1 = MagicMock()
        rs1.contact_fraction_per_replicate = [0.5, 0.3]
        rs2 = MagicMock()
        rs2.contact_fraction_per_replicate = [0.1, 0.0]
        rs3 = MagicMock()
        rs3.contact_fraction_per_replicate = [0.0, 0.0]
        agg.residue_stats = [rs1, rs2, rs3]

        coverages = ContactsAnalysis._compute_coverage_per_replicate(agg)
        assert len(coverages) == 2
        assert coverages[0] == pytest.approx(2 / 3)  # rs1 + rs2 contacted in rep 0
        assert coverages[1] == pytest.approx(1 / 3)  # only rs1 contacted in rep 1

    def test_contact_fraction_per_replicate(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        agg = MagicMock()
        agg.n_replicates = 2

        rs1 = MagicMock()
        rs1.contact_fraction_per_replicate = [0.4, 0.6]
        rs2 = MagicMock()
        rs2.contact_fraction_per_replicate = [0.2, 0.8]
        agg.residue_stats = [rs1, rs2]

        fractions = ContactsAnalysis._compute_contact_fraction_per_replicate(agg)
        assert len(fractions) == 2
        assert fractions[0] == pytest.approx(0.3)  # (0.4 + 0.2) / 2
        assert fractions[1] == pytest.approx(0.7)  # (0.6 + 0.8) / 2


# ---------------------------------------------------------------------------
# Residue set validation
# ---------------------------------------------------------------------------


class TestResidueSetValidation:
    """Test _validate_residue_sets."""

    def test_matching_residue_sets(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        cond_a = MagicMock()
        cond_a.label = "A"
        data_a = {"agg_result": MagicMock()}
        rs_a = [MagicMock(protein_resid=i) for i in [1, 2, 3]]
        data_a["agg_result"].residue_stats = rs_a

        cond_b = MagicMock()
        cond_b.label = "B"
        data_b = {"agg_result": MagicMock()}
        rs_b = [MagicMock(protein_resid=i) for i in [1, 2, 3]]
        data_b["agg_result"].residue_stats = rs_b

        # Should not raise
        ContactsAnalysis._validate_residue_sets([(cond_a, data_a), (cond_b, data_b)])

    def test_mismatched_residue_sets(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        cond_a = MagicMock()
        cond_a.label = "A"
        data_a = {"agg_result": MagicMock()}
        data_a["agg_result"].residue_stats = [MagicMock(protein_resid=i) for i in [1, 2, 3]]

        cond_b = MagicMock()
        cond_b.label = "B"
        data_b = {"agg_result": MagicMock()}
        data_b["agg_result"].residue_stats = [MagicMock(protein_resid=i) for i in [1, 2, 4]]

        with pytest.raises(ValueError, match="Residue set mismatch"):
            ContactsAnalysis._validate_residue_sets([(cond_a, data_a), (cond_b, data_b)])


# ---------------------------------------------------------------------------
# filter_conditions
# ---------------------------------------------------------------------------


class TestFilterConditions:
    """Test filter_conditions polymer detection logic."""

    def test_stale_cached_results_do_not_bypass_topology_selection(self, tmp_path):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()

        # A stale cache alone should not prove polymer presence
        projects_dir = tmp_path / "projects"
        contacts_dir = projects_dir / "analysis" / "contacts"
        contacts_dir.mkdir(parents=True)
        (contacts_dir / "contacts_rep1.json").write_text("{}")
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        (run_dir / "solvated_system.pdb").write_text("ATOM ...")

        mock_sim_config = MagicMock()
        mock_sim_config.output.projects_directory = projects_dir
        mock_sim_config.get_working_directory.return_value = run_dir

        cond = Condition(
            label="With Polymer",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim_config,
        )

        with patch("MDAnalysis.Universe") as MockUniverse:
            mock_universe = MagicMock()
            mock_universe.select_atoms.return_value = MagicMock(__len__=lambda s: 0)
            MockUniverse.return_value = mock_universe

            result = analysis.filter_conditions(
                [cond], settings=ContactsSettings(polymer_selection="chainID C and resname PEG")
            )

        assert result == []
        mock_universe.select_atoms.assert_called_with("chainID C and resname PEG")

    def test_condition_without_polymer_excluded(self, tmp_path):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()

        # No cached results, topology exists but no polymer atoms
        projects_dir = tmp_path / "projects"
        projects_dir.mkdir()
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()

        mock_sim_config = MagicMock()
        mock_sim_config.output.projects_directory = projects_dir
        mock_sim_config.get_working_directory.return_value = run_dir

        # Create a dummy topology file
        (run_dir / "solvated_system.pdb").write_text("ATOM ...")

        cond = Condition(
            label="No Polymer",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim_config,
        )

        # Mock MDAnalysis to return empty selection
        with patch("MDAnalysis.Universe") as MockUniverse:
            mock_universe = MagicMock()
            mock_universe.select_atoms.return_value = MagicMock(__len__=lambda s: 0)
            MockUniverse.return_value = mock_universe

            result = analysis.filter_conditions([cond])

        assert len(result) == 0

    def test_error_during_check_includes_condition(self, tmp_path):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()

        mock_sim_config = MagicMock()
        mock_sim_config.output.projects_directory = tmp_path / "projects"
        mock_sim_config.get_working_directory.side_effect = OSError("boom")

        cond = Condition(
            label="ErrorCond",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim_config,
        )

        result = analysis.filter_conditions([cond])
        assert len(result) == 1  # Included because of error (fail-open)

    def test_missing_topology_includes_condition(self, tmp_path):
        """Conditions are fail-open when no replicate topology can be inspected."""
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()

        mock_sim_config = MagicMock()
        mock_sim_config.output.projects_directory = tmp_path / "projects"
        mock_sim_config.get_working_directory.return_value = run_dir
        cond = Condition(
            label="Unchecked",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim_config,
        )

        result = analysis.filter_conditions([cond])

        assert result == [cond]


# ---------------------------------------------------------------------------
# compare — full override
# ---------------------------------------------------------------------------


def _make_condition_and_data(label: str, n_reps: int = 3, n_residues: int = 5):
    """Build a mock condition + data dict for compare() tests."""
    from polyzymd.analyses.base import Condition

    mock_sim_config = MagicMock()
    cond = Condition(
        label=label,
        config_path=Path(f"/tmp/{label}/config.yaml"),
        replicates=tuple(range(1, n_reps + 1)),
        sim_config=mock_sim_config,
    )

    agg = _make_mock_agg_result(n_replicates=n_reps, n_residues=n_residues)
    data = {
        "agg_result": agg,
        "coverage_per_replicate": [0.8 + i * 0.01 for i in range(n_reps)],
        "contact_fraction_per_replicate": [0.3 + i * 0.005 for i in range(n_reps)],
    }
    return cond, data


class TestCompare:
    """Test compare() full override."""

    def _make_ctx(
        self,
        tmp_path,
        n_conditions: int = 3,
        control: str | None = "Control",
    ):
        """Create a ComparisonContext with mock data."""
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsSettings

        labels = ["Control", "Treatment A", "Treatment B"][:n_conditions]
        conditions = []
        for label in labels:
            mock_sim = MagicMock()
            conditions.append(
                Condition(
                    label=label,
                    config_path=Path(f"/tmp/{label}/config.yaml"),
                    replicates=(1, 2, 3),
                    sim_config=mock_sim,
                )
            )

        analysis_dirs = {}
        for cond in conditions:
            agg_dir = tmp_path / cond.label / "contacts" / "aggregated"
            agg_dir.mkdir(parents=True, exist_ok=True)
            analysis_dirs[cond.label] = tmp_path / cond.label / "contacts"

        return ComparisonContext(
            name="test_comparison",
            conditions=conditions,
            excluded_conditions=[],
            control_label=control,
            analysis_dirs=analysis_dirs,
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=ContactsSettings(),
            recompute=False,
        )

    def test_compare_returns_result(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        ctx = self._make_ctx(tmp_path, n_conditions=3, control="Control")

        # Mock _load_aggregated_result to return proper aggregated results
        mock_agg_results = {}
        for cond in ctx.conditions:
            mock_agg = _make_mock_agg_result(n_replicates=3, n_residues=5)
            mock_agg_results[cond.label] = mock_agg

        def side_effect(agg_dir):
            # Extract label from path
            label = agg_dir.parent.parent.name
            return mock_agg_results.get(label)

        with (
            patch.object(analysis, "_load_aggregated_result", side_effect=side_effect),
            patch.object(analysis, "_load_or_compute_binding_preference", return_value=None),
        ):
            result = analysis.compare(ctx)

        assert result is not None
        assert result.name == "test_comparison"
        assert len(result.conditions) == 3
        assert result.control_label == "Control"
        assert len(result.ranking_by_coverage) == 3
        assert len(result.ranking_by_contact_fraction) == 3

    def test_compare_pairwise_with_control(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        ctx = self._make_ctx(tmp_path, n_conditions=3, control="Control")

        mock_agg_results = {}
        for cond in ctx.conditions:
            mock_agg_results[cond.label] = _make_mock_agg_result(3, 5)

        def side_effect(agg_dir):
            label = agg_dir.parent.parent.name
            return mock_agg_results.get(label)

        with (
            patch.object(analysis, "_load_aggregated_result", side_effect=side_effect),
            patch.object(analysis, "_load_or_compute_binding_preference", return_value=None),
        ):
            result = analysis.compare(ctx)

        # With 3 conditions and a control, expect 2 pairwise comparisons
        assert len(result.pairwise_comparisons) == 2
        for comp in result.pairwise_comparisons:
            assert comp.condition_a == "Control"
            assert len(comp.aggregate_comparisons) == 2  # coverage + contact_fraction

    def test_compare_pairwise_without_control(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        ctx = self._make_ctx(tmp_path, n_conditions=3, control=None)

        mock_agg_results = {}
        for cond in ctx.conditions:
            mock_agg_results[cond.label] = _make_mock_agg_result(3, 5)

        def side_effect(agg_dir):
            label = agg_dir.parent.parent.name
            return mock_agg_results.get(label)

        with (
            patch.object(analysis, "_load_aggregated_result", side_effect=side_effect),
            patch.object(analysis, "_load_or_compute_binding_preference", return_value=None),
        ):
            result = analysis.compare(ctx)

        # Without control: all pairs = C(3,2) = 3
        assert len(result.pairwise_comparisons) == 3

    def test_compare_falls_back_when_control_aggregate_missing(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        ctx = self._make_ctx(tmp_path, n_conditions=3, control="Control")
        mock_agg_results = {
            "Treatment A": _make_mock_agg_result(3, 5),
            "Treatment B": _make_mock_agg_result(3, 5),
        }

        def side_effect(agg_dir, **_kwargs):
            label = agg_dir.parent.parent.name
            return mock_agg_results.get(label)

        with (
            patch.object(analysis, "_load_validated_aggregated_result", side_effect=side_effect),
            patch.object(analysis, "_load_or_compute_binding_preference", return_value=None),
        ):
            result = analysis.compare(ctx)

        assert result is not None
        assert result.control_label is None
        assert len(result.pairwise_comparisons) == 1
        assert result.pairwise_comparisons[0].condition_a == "Treatment A"
        assert result.pairwise_comparisons[0].condition_b == "Treatment B"

    def test_compare_anova_with_three_conditions(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        ctx = self._make_ctx(tmp_path, n_conditions=3)

        mock_agg_results = {}
        for cond in ctx.conditions:
            mock_agg_results[cond.label] = _make_mock_agg_result(3, 5)

        def side_effect(agg_dir):
            label = agg_dir.parent.parent.name
            return mock_agg_results.get(label)

        with (
            patch.object(analysis, "_load_aggregated_result", side_effect=side_effect),
            patch.object(analysis, "_load_or_compute_binding_preference", return_value=None),
        ):
            result = analysis.compare(ctx)

        # ANOVA with 3+ conditions: 2 summaries (coverage + contact_fraction)
        assert len(result.anova) == 2
        metrics = {a.metric for a in result.anova}
        assert metrics == {"coverage", "mean_contact_fraction"}

    def test_compare_no_anova_with_two_conditions(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        ctx = self._make_ctx(tmp_path, n_conditions=2, control="Control")

        mock_agg_results = {}
        for cond in ctx.conditions:
            mock_agg_results[cond.label] = _make_mock_agg_result(3, 5)

        def side_effect(agg_dir):
            label = agg_dir.parent.parent.name
            return mock_agg_results.get(label)

        with (
            patch.object(analysis, "_load_aggregated_result", side_effect=side_effect),
            patch.object(analysis, "_load_or_compute_binding_preference", return_value=None),
        ):
            result = analysis.compare(ctx)

        # No ANOVA with < 3 conditions
        assert len(result.anova) == 0

    def test_compare_returns_result_with_single_condition(self, tmp_path):
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()

        mock_sim = MagicMock()
        cond = Condition(
            label="Only",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=mock_sim,
        )
        ctx = ComparisonContext(
            name="test",
            conditions=[cond],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"Only": tmp_path / "Only" / "contacts"},
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=ContactsSettings(),
            recompute=False,
        )

        with (
            patch.object(
                analysis, "_load_aggregated_result", return_value=_make_mock_agg_result(3, 5)
            ),
            patch.object(analysis, "_load_or_compute_binding_preference", return_value=None),
        ):
            result = analysis.compare(ctx)

        assert result is not None
        assert len(result.conditions) == 1
        assert result.pairwise_comparisons == []

    def test_compare_uses_aggregated_results_when_recompute_true(self, tmp_path):
        """Finalize should consume in-memory aggregates even with recompute enabled."""
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        condition = Condition(
            label="Finalized",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        analysis_dir = tmp_path / "Finalized" / "contacts"
        analysis_dir.mkdir(parents=True)
        aggregate = _make_valid_agg_result(settings)
        ctx = ComparisonContext(
            name="test",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"Finalized": analysis_dir},
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=settings,
            aggregated_results={"Finalized": aggregate},
            recompute=True,
        )

        with (
            patch.object(analysis, "_load_validated_aggregated_result") as load_mock,
            patch.object(analysis, "_load_or_compute_binding_preference", return_value=None),
        ):
            result = analysis.compare(ctx)

        assert result is not None
        assert result.conditions[0].label == "Finalized"
        load_mock.assert_not_called()

    def test_compare_loads_existing_aggregate_when_recompute_true(self, tmp_path):
        """Comparison-stage reads should not suppress valid aggregate JSON files."""
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        condition = Condition(
            label="Cached",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        analysis_dir = tmp_path / "Cached" / "contacts"
        agg_dir = analysis_dir / "aggregated"
        agg_dir.mkdir(parents=True)
        _make_valid_agg_result(settings).save(agg_dir / "result.json")
        ctx = ComparisonContext(
            name="test",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"Cached": analysis_dir},
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=settings,
            recompute=True,
        )

        with patch.object(analysis, "_load_or_compute_binding_preference", return_value=None):
            result = analysis.compare(ctx)

        assert result is not None
        assert result.conditions[0].label == "Cached"

    def test_compare_excluded_conditions_recorded(self, tmp_path):
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()

        mock_sim = MagicMock()
        conditions = [
            Condition(
                label=label,
                config_path=Path(f"/tmp/{label}/config.yaml"),
                replicates=(1, 2, 3),
                sim_config=mock_sim,
            )
            for label in ["A", "B"]
        ]
        excluded = [
            Condition(
                label="No Polymer",
                config_path=Path("/tmp/np/config.yaml"),
                replicates=(1,),
                sim_config=mock_sim,
            )
        ]

        analysis_dirs = {}
        for cond in conditions:
            agg_dir = tmp_path / cond.label / "contacts" / "aggregated"
            agg_dir.mkdir(parents=True, exist_ok=True)
            analysis_dirs[cond.label] = tmp_path / cond.label / "contacts"

        ctx = ComparisonContext(
            name="test",
            conditions=conditions,
            excluded_conditions=excluded,
            control_label=None,
            analysis_dirs=analysis_dirs,
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=ContactsSettings(),
            recompute=False,
        )

        mock_agg = _make_mock_agg_result(3, 5)

        with (
            patch.object(analysis, "_load_aggregated_result", return_value=mock_agg),
            patch.object(analysis, "_load_or_compute_binding_preference", return_value=None),
        ):
            result = analysis.compare(ctx)

        assert "No Polymer" in result.excluded_conditions

    def test_aggregate_context_validation_rejects_settings_mismatch(self, tmp_path):
        """Compare/finalize should not accept aggregates from different settings."""
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.shared.config_hash import settings_fingerprint

        analysis = ContactsAnalysis()
        current_settings = ContactsSettings(cutoff=4.5)
        stale_settings = ContactsSettings(cutoff=4.0)
        summary = SimpleNamespace(
            equilibration_time=10.0,
            equilibration_unit="ns",
            metadata={"settings_fingerprint": settings_fingerprint(stale_settings)},
            config_hash="unknown",
        )

        assert not analysis._cache_matches_context(
            summary,
            settings=current_settings,
            equilibration="10ns",
            sim_config=_make_hashable_sim_config(tmp_path),
        )

    def test_compare_accepts_successful_replicate_subset_aggregate(self, tmp_path):
        """Finalize should accept aggregates that record successful replicate subsets."""
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._aggregator import (
            AggregatedContactResult,
            AggregatedResidueStats,
        )
        from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        condition = Condition(
            label="Subset",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        analysis_dir = tmp_path / "Subset" / "contacts"
        agg_dir = analysis_dir / "aggregated"
        agg_dir.mkdir(parents=True)
        AggregatedContactResult(
            n_replicates=2,
            replicates=[1, 2],
            residue_stats=[
                AggregatedResidueStats(
                    protein_resid=1,
                    protein_resname="ALA",
                    contact_fraction_mean=0.5,
                    contact_fraction_per_replicate=[0.4, 0.6],
                )
            ],
            total_frames_per_replicate=[10, 10],
            criteria_cutoff=4.5,
            coverage_mean=1.0,
            coverage_sem=0.0,
            mean_contact_fraction=0.5,
            mean_contact_fraction_sem=0.1,
            equilibration_time=10.0,
            equilibration_unit="ns",
            metadata={
                "settings_fingerprint": contacts_settings_fingerprint(settings),
                "replicates": [1, 2],
            },
        ).save(agg_dir / "result.json")
        ctx = ComparisonContext(
            name="test",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"Subset": analysis_dir},
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=settings,
            recompute=False,
        )

        with patch.object(analysis, "_load_or_compute_binding_preference", return_value=None):
            result = analysis.compare(ctx)

        assert result is not None
        assert result.conditions[0].n_replicates == 2

    def test_compare_rejects_replicate_subset_below_minimum(self, tmp_path):
        """Finalize should not accept successful subsets below plugin minimum."""
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._aggregator import (
            AggregatedContactResult,
            AggregatedResidueStats,
        )
        from polyzymd.analyses.shared.config_hash import settings_fingerprint

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        condition = Condition(
            label="Subset",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        analysis_dir = tmp_path / "Subset" / "contacts"
        agg_dir = analysis_dir / "aggregated"
        agg_dir.mkdir(parents=True)
        AggregatedContactResult(
            n_replicates=1,
            replicates=[1],
            residue_stats=[
                AggregatedResidueStats(
                    protein_resid=1,
                    protein_resname="ALA",
                    contact_fraction_mean=0.5,
                    contact_fraction_per_replicate=[0.5],
                )
            ],
            total_frames_per_replicate=[10],
            criteria_cutoff=4.5,
            coverage_mean=1.0,
            coverage_sem=0.0,
            mean_contact_fraction=0.5,
            mean_contact_fraction_sem=0.0,
            equilibration_time=10.0,
            equilibration_unit="ns",
            metadata={"settings_fingerprint": settings_fingerprint(settings), "replicates": [1]},
        ).save(agg_dir / "result.json")
        ctx = ComparisonContext(
            name="test",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"Subset": analysis_dir},
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=settings,
            recompute=False,
        )

        with patch.object(analysis, "_load_or_compute_binding_preference", return_value=None):
            result = analysis.compare(ctx)

        assert result is None


# ---------------------------------------------------------------------------
# Pairwise comparison details
# ---------------------------------------------------------------------------


class TestPairwiseComparison:
    """Test _compare_contacts_pair static method."""

    def test_pair_produces_two_metrics(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        summary_a = MagicMock()
        summary_a.coverage_mean = 0.8
        summary_a.coverage_sem = 0.05
        summary_a.mean_contact_fraction = 0.3
        summary_a.mean_contact_fraction_sem = 0.02

        summary_b = MagicMock()
        summary_b.coverage_mean = 0.9
        summary_b.coverage_sem = 0.04
        summary_b.mean_contact_fraction = 0.4
        summary_b.mean_contact_fraction_sem = 0.03

        data_a = {
            "coverage_per_replicate": [0.78, 0.80, 0.82],
            "contact_fraction_per_replicate": [0.28, 0.30, 0.32],
        }
        data_b = {
            "coverage_per_replicate": [0.88, 0.90, 0.92],
            "contact_fraction_per_replicate": [0.38, 0.40, 0.42],
        }

        comp = ContactsAnalysis._compare_contacts_pair(
            "A", summary_a, data_a, "B", summary_b, data_b
        )

        assert comp.condition_a == "A"
        assert comp.condition_b == "B"
        assert len(comp.aggregate_comparisons) == 2

        metrics = {ac.metric for ac in comp.aggregate_comparisons}
        assert metrics == {"coverage", "mean_contact_fraction"}

    def test_pair_direction_labels(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        summary_a = MagicMock(
            coverage_mean=0.5,
            coverage_sem=0.05,
            mean_contact_fraction=0.2,
            mean_contact_fraction_sem=0.02,
        )
        summary_b = MagicMock(
            coverage_mean=0.9,
            coverage_sem=0.04,
            mean_contact_fraction=0.5,
            mean_contact_fraction_sem=0.03,
        )

        data_a = {
            "coverage_per_replicate": [0.48, 0.50, 0.52],
            "contact_fraction_per_replicate": [0.18, 0.20, 0.22],
        }
        data_b = {
            "coverage_per_replicate": [0.88, 0.90, 0.92],
            "contact_fraction_per_replicate": [0.48, 0.50, 0.52],
        }

        comp = ContactsAnalysis._compare_contacts_pair(
            "Control", summary_a, data_a, "Treatment", summary_b, data_b
        )

        for ac in comp.aggregate_comparisons:
            assert ac.direction in ("increased", "decreased", "unchanged")


# ---------------------------------------------------------------------------
# ANOVA
# ---------------------------------------------------------------------------


class TestANOVA:
    """Test _compute_contacts_anova."""

    def test_anova_returns_two_summaries(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        condition_data = []
        for i in range(3):
            cond = MagicMock()
            cond.label = f"Cond{i}"
            data = {
                "coverage_per_replicate": [0.7 + i * 0.05 + j * 0.01 for j in range(3)],
                "contact_fraction_per_replicate": [0.2 + i * 0.03 + j * 0.005 for j in range(3)],
            }
            condition_data.append((cond, data))

        results = ContactsAnalysis._compute_contacts_anova(condition_data)
        assert len(results) == 2
        metrics = {r.metric for r in results}
        assert metrics == {"coverage", "mean_contact_fraction"}

        for r in results:
            assert hasattr(r, "f_statistic")
            assert hasattr(r, "p_value")
            assert hasattr(r, "significant")


# ---------------------------------------------------------------------------
# AggregatedResultClass and _deserialize_result
# ---------------------------------------------------------------------------


class TestDeserializeResult:
    """Test AggregatedResultClass and _deserialize_result hook."""

    def test_aggregated_result_class_set(self):
        from polyzymd.analyses.contacts import ContactsAnalysis
        from polyzymd.analyses.contacts._aggregator import AggregatedContactResult

        assert ContactsAnalysis.AggregatedResultClass is AggregatedContactResult

    def test_deserialize_delegates_to_load(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()

        mock_result = MagicMock()
        with patch.object(analysis.AggregatedResultClass, "load", return_value=mock_result):
            result = analysis._deserialize_result(tmp_path / "test.json")

        assert result is mock_result


# ---------------------------------------------------------------------------
# plot delegation
# ---------------------------------------------------------------------------


_PLOT_FUNCTIONS = [
    "_plot_contact_fraction_profile",
    "_plot_residence_time_profile",
    "_plot_cf_by_aa_class_bars",
    "_plot_cf_by_partition_bars",
    "_plot_rt_by_aa_class_bars",
    "_plot_rt_by_partition_bars",
    "_plot_user_partition_bars",
    "_plot_system_coverage_bars",
    "_plot_system_coverage_heatmap",
    "_plot_binding_preference_bars",
    "_plot_binding_preference_heatmap",
]


def _write_contacts_plot_aggregate(
    analysis_dir: Path,
    settings=None,
    replicates: tuple[int, ...] = (1, 2),
) -> Path:
    """Write a valid contacts aggregate for plot orchestration tests."""

    from polyzymd.analyses.contacts import ContactsSettings
    from polyzymd.analyses.contacts._aggregator import (
        AggregatedContactResult,
        AggregatedResidueStats,
    )
    from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint

    resolved_settings = settings or ContactsSettings()
    aggregated_dir = analysis_dir / "aggregated"
    aggregated_dir.mkdir(parents=True, exist_ok=True)
    aggregate = AggregatedContactResult(
        n_replicates=len(replicates),
        replicates=list(replicates),
        residue_stats=[
            AggregatedResidueStats(
                protein_resid=1,
                protein_resname="ALA",
                protein_group="nonpolar",
                contact_fraction_mean=0.5,
                contact_fraction_sem=0.1,
                contact_fraction_per_replicate=[0.4, 0.6][: len(replicates)],
                by_polymer_type={"PEG": (0.5, 0.1)},
                by_polymer_type_per_replicate={"PEG": [0.4, 0.6][: len(replicates)]},
            ),
            AggregatedResidueStats(
                protein_resid=2,
                protein_resname="ASP",
                protein_group="charged",
                contact_fraction_mean=0.2,
                contact_fraction_sem=0.05,
                contact_fraction_per_replicate=[0.1, 0.3][: len(replicates)],
                by_polymer_type={"PEG": (0.2, 0.05)},
                by_polymer_type_per_replicate={"PEG": [0.1, 0.3][: len(replicates)]},
            ),
        ],
        total_frames_per_replicate=[10] * len(replicates),
        criteria_cutoff=resolved_settings.cutoff,
        coverage_mean=1.0,
        mean_contact_fraction=0.35,
        metadata={
            "settings_fingerprint": contacts_settings_fingerprint(resolved_settings),
            "replicates": list(replicates),
        },
    )
    return aggregate.save(aggregated_dir / "result.json")


class TestPlot:
    """Test plot() delegates to private module-level plotting functions."""

    def test_plot_creates_output_dir(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.config.comparison import PlotSettings

        analysis = ContactsAnalysis()

        mock_sim = MagicMock()
        conditions = [
            Condition(
                label="A",
                config_path=Path("/tmp/a/config.yaml"),
                replicates=(1, 2),
                sim_config=mock_sim,
            )
        ]

        output_dir = tmp_path / "plots"
        ctx = PlotContext(
            conditions=conditions,
            analysis_dirs={"A": tmp_path / "A" / "contacts"},
            results_dir=tmp_path / "results",
            output_dir=output_dir,
            settings=ContactsSettings(),
            plot_settings=PlotSettings(),
        )
        _write_contacts_plot_aggregate(ctx.analysis_dirs["A"], ctx.settings)

        # Mock all 11 private plot functions to return empty lists
        patches = {
            fn: patch(f"polyzymd.analyses.contacts.{fn}", return_value=[]) for fn in _PLOT_FUNCTIONS
        }
        mocks = {name: p.start() for name, p in patches.items()}
        try:
            analysis.plot(ctx)
        finally:
            for p in patches.values():
                p.stop()

        assert output_dir.exists()
        # Verify all 11 functions were called
        for name, mock_fn in mocks.items():
            assert mock_fn.called, f"{name} was not called"

    def test_plot_returns_combined_paths(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.config.comparison import PlotSettings

        analysis = ContactsAnalysis()

        mock_sim = MagicMock()
        conditions = [
            Condition(
                label="A",
                config_path=Path("/tmp/a/config.yaml"),
                replicates=(1, 2),
                sim_config=mock_sim,
            )
        ]

        ctx = PlotContext(
            conditions=conditions,
            analysis_dirs={"A": tmp_path / "A" / "contacts"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=ContactsSettings(),
            plot_settings=PlotSettings(),
        )
        _write_contacts_plot_aggregate(ctx.analysis_dirs["A"], ctx.settings)

        # Make first two functions return paths, rest return empty
        path_a = tmp_path / "plot_a.png"
        path_b = tmp_path / "plot_b.png"

        patches = {}
        for fn in _PLOT_FUNCTIONS:
            patches[fn] = patch(f"polyzymd.analyses.contacts.{fn}", return_value=[])

        patches[_PLOT_FUNCTIONS[0]] = patch(
            f"polyzymd.analyses.contacts.{_PLOT_FUNCTIONS[0]}",
            return_value=[path_a],
        )
        patches[_PLOT_FUNCTIONS[1]] = patch(
            f"polyzymd.analyses.contacts.{_PLOT_FUNCTIONS[1]}",
            return_value=[path_b],
        )

        for p in patches.values():
            p.start()
        try:
            result = analysis.plot(ctx)
        finally:
            for p in patches.values():
                p.stop()

        assert path_a in result
        assert path_b in result
        assert len(result) == 2

    def test_plot_catches_plotter_exceptions(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.config.comparison import PlotSettings

        analysis = ContactsAnalysis()

        mock_sim = MagicMock()
        conditions = [
            Condition(
                label="A",
                config_path=Path("/tmp/a/config.yaml"),
                replicates=(1, 2),
                sim_config=mock_sim,
            )
        ]

        ctx = PlotContext(
            conditions=conditions,
            analysis_dirs={"A": tmp_path / "A" / "contacts"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=ContactsSettings(),
            plot_settings=PlotSettings(),
        )
        _write_contacts_plot_aggregate(ctx.analysis_dirs["A"], ctx.settings)

        # Make all plot functions raise
        patches = {
            fn: patch(
                f"polyzymd.analyses.contacts.{fn}",
                side_effect=RuntimeError("plot failed"),
            )
            for fn in _PLOT_FUNCTIONS
        }
        for p in patches.values():
            p.start()
        try:
            result = analysis.plot(ctx)
        finally:
            for p in patches.values():
                p.stop()

        # Should return empty list, not raise
        assert result == []

    def test_plot_passes_plot_settings(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.config.comparison import PlotSettings

        analysis = ContactsAnalysis()
        ps = PlotSettings()

        mock_sim = MagicMock()
        conditions = [
            Condition(
                label="A",
                config_path=Path("/tmp/a/config.yaml"),
                replicates=(1, 2),
                sim_config=mock_sim,
            )
        ]

        ctx = PlotContext(
            conditions=conditions,
            analysis_dirs={"A": tmp_path / "A" / "contacts"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=ContactsSettings(),
            plot_settings=ps,
        )
        _write_contacts_plot_aggregate(ctx.analysis_dirs["A"], ctx.settings)

        patches = {
            fn: patch(f"polyzymd.analyses.contacts.{fn}", return_value=[]) for fn in _PLOT_FUNCTIONS
        }
        mocks = {name: p.start() for name, p in patches.items()}
        try:
            analysis.plot(ctx)
        finally:
            for p in patches.values():
                p.stop()

        # Verify PlotSettings was passed as 4th positional arg to each function
        for name, mock_fn in mocks.items():
            assert mock_fn.call_args[0][3] is ps, f"{name} did not receive PlotSettings as 4th arg"

    def test_plot_skips_residence_time_plotters_when_disabled(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.config.comparison import PlotSettings

        analysis = ContactsAnalysis()
        condition = Condition(
            label="A",
            config_path=Path("/tmp/a/config.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        ctx = PlotContext(
            conditions=[condition],
            analysis_dirs={"A": tmp_path / "A" / "contacts"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=ContactsSettings(compute_residence_times=False),
            plot_settings=PlotSettings(),
        )
        _write_contacts_plot_aggregate(ctx.analysis_dirs["A"], ctx.settings)
        residence_plotters = {
            "_plot_residence_time_profile",
            "_plot_rt_by_aa_class_bars",
            "_plot_rt_by_partition_bars",
        }

        patches = {
            fn: patch(f"polyzymd.analyses.contacts.{fn}", return_value=[]) for fn in _PLOT_FUNCTIONS
        }
        mocks = {name: patcher.start() for name, patcher in patches.items()}
        try:
            analysis.plot(ctx)
        finally:
            for patcher in patches.values():
                patcher.stop()

        for name in residence_plotters:
            assert not mocks[name].called
        for name in set(_PLOT_FUNCTIONS) - residence_plotters:
            assert mocks[name].called

    def test_plot_accepts_successful_replicate_subset_aggregate(self, tmp_path):
        """Plot validation should accept aggregates with successful replicate subsets."""
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._aggregator import (
            AggregatedContactResult,
            AggregatedResidueStats,
        )
        from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint
        from polyzymd.config.comparison import PlotSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        condition = Condition(
            label="A",
            config_path=Path("/tmp/a/config.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        analysis_dir = tmp_path / "A" / "contacts"
        agg_dir = analysis_dir / "aggregated"
        agg_dir.mkdir(parents=True)
        AggregatedContactResult(
            n_replicates=2,
            replicates=[1, 2],
            residue_stats=[
                AggregatedResidueStats(
                    protein_resid=1,
                    protein_resname="ALA",
                    contact_fraction_mean=0.5,
                    contact_fraction_per_replicate=[0.4, 0.6],
                )
            ],
            total_frames_per_replicate=[10, 10],
            criteria_cutoff=4.5,
            coverage_mean=1.0,
            mean_contact_fraction=0.5,
            equilibration_time=10.0,
            equilibration_unit="ns",
            metadata={
                "settings_fingerprint": contacts_settings_fingerprint(settings),
                "replicates": [1, 2],
            },
        ).save(agg_dir / "result.json")
        ctx = PlotContext(
            conditions=[condition],
            analysis_dirs={"A": analysis_dir},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=settings,
            plot_settings=PlotSettings(),
        )

        with patch(f"polyzymd.analyses.contacts.{_PLOT_FUNCTIONS[0]}", return_value=[]) as plot_fn:
            patches = [
                patch(f"polyzymd.analyses.contacts.{fn}", return_value=[])
                for fn in _PLOT_FUNCTIONS[1:]
            ]
            for patcher in patches:
                patcher.start()
            try:
                analysis.plot(ctx)
            finally:
                for patcher in patches:
                    patcher.stop()

        assert plot_fn.called

    def test_plot_skips_condition_without_aggregate_data(self, tmp_path, caplog):
        """Plot should not invoke plotters when aggregate artifacts are absent."""
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.config.comparison import PlotSettings

        analysis = ContactsAnalysis()
        condition = Condition(
            label="A",
            config_path=Path("/tmp/a/config.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        ctx = PlotContext(
            conditions=[condition],
            analysis_dirs={"A": tmp_path / "A" / "contacts"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=ContactsSettings(),
            plot_settings=PlotSettings(),
        )

        with patch(f"polyzymd.analyses.contacts.{_PLOT_FUNCTIONS[0]}", return_value=[]) as plot_fn:
            with caplog.at_level("INFO", logger="polyzymd.analyses.contacts"):
                result = analysis.plot(ctx)

        assert result == []
        assert not plot_fn.called
        assert "no aggregated contacts JSON found" in caplog.text

    def test_contact_plotter_loads_current_aggregate_layout(self, tmp_path):
        """Contact plotters should load canonical aggregates below aggregated/."""
        from polyzymd.analyses.contacts._plotters import _plot_contact_fraction_profile
        from polyzymd.config.comparison import PlotSettings

        analysis_dir = tmp_path / "A" / "contacts"
        _write_contacts_plot_aggregate(analysis_dir)
        data = {
            "A": {
                "analysis_dir": analysis_dir,
                "aggregated_dir": analysis_dir / "aggregated",
            },
            "__meta__": {"settings": None},
        }

        paths = _plot_contact_fraction_profile(data, ["A"], tmp_path / "plots", PlotSettings())

        assert paths
        assert all(path.exists() for path in paths)

    def test_plot_propagates_validated_sidecar_when_canonical_is_malformed(self, tmp_path):
        """Plot gating should pass the exact validated aggregate artifact to plotters."""
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.config.comparison import PlotSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        condition = Condition(
            label="A",
            config_path=Path("/tmp/a/config.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        analysis_dir = tmp_path / "A" / "contacts"
        canonical = _write_contacts_plot_aggregate(analysis_dir, settings)
        sidecar = analysis._aggregated_sidecar_path(
            canonical.parent,
            settings,
            "0ns",
            condition.replicates,
        )
        canonical.rename(sidecar)
        canonical.write_text("{malformed json")
        ctx = PlotContext(
            conditions=[condition],
            analysis_dirs={"A": analysis_dir},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=settings,
            plot_settings=PlotSettings(),
        )
        captured: dict[str, Path | None] = {}

        def capture_artifact(data, labels, output_dir, plot_settings):
            captured["artifact"] = data["A"].get("aggregated_result_path")
            return []

        patches = {
            fn: patch(f"polyzymd.analyses.contacts.{fn}", return_value=[]) for fn in _PLOT_FUNCTIONS
        }
        patches[_PLOT_FUNCTIONS[0]] = patch(
            f"polyzymd.analyses.contacts.{_PLOT_FUNCTIONS[0]}",
            side_effect=capture_artifact,
        )
        for patcher in patches.values():
            patcher.start()
        try:
            analysis.plot(ctx)
        finally:
            for patcher in patches.values():
                patcher.stop()

        assert captured["artifact"] == sidecar

    def test_contact_result_loader_falls_back_after_unloadable_candidate(self, tmp_path):
        """Contact plot loaders should continue after a stale malformed candidate."""
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._plotters import _load_aggregated_contact_results

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        replicates = (1, 2)
        analysis_dir = tmp_path / "A" / "contacts"
        canonical = _write_contacts_plot_aggregate(analysis_dir, settings, replicates=replicates)
        sidecar = analysis._aggregated_sidecar_path(canonical.parent, settings, "0ns", replicates)
        canonical.rename(sidecar)
        canonical.write_text("{malformed json")
        data = {
            "A": {
                "analysis_dir": analysis_dir,
                "aggregated_dir": analysis_dir / "aggregated",
                "aggregated_result_path": canonical,
            },
            "__meta__": {"settings": settings},
        }

        results = _load_aggregated_contact_results(data, ["A"])

        assert "A" in results
        assert results["A"].mean_contact_fraction == pytest.approx(0.35)


class TestPartitionDefinitions:
    """Tests for contacts partition auto-fill behavior."""

    def test_load_partition_definitions_does_not_mutate_settings(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._plotters import _load_partition_definitions

        settings = ContactsSettings(
            protein_groups={"nterm": [1, 2], "cterm": [5, 6]},
            protein_partitions={"left": ["nterm"], "right": ["cterm"]},
        )
        original_groups = {k: list(v) for k, v in settings.protein_groups.items()}
        original_partitions = {k: list(v) for k, v in settings.protein_partitions.items()}

        data = {"__meta__": {"settings": settings}}
        groups, partitions = _load_partition_definitions(data, all_resids={1, 2, 3, 4, 5, 6})

        assert settings.protein_groups == original_groups
        assert settings.protein_partitions == original_partitions
        assert groups is not settings.protein_groups
        assert partitions is not settings.protein_partitions

    def test_auto_fill_creates_unique_remainder_group_per_partition(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._plotters import _load_partition_definitions

        settings = ContactsSettings(
            protein_groups={"nterm": [1, 2], "cterm": [5, 6]},
            protein_partitions={"left": ["nterm"], "right": ["cterm"]},
        )

        data = {"__meta__": {"settings": settings}}
        groups, partitions = _load_partition_definitions(data, all_resids={1, 2, 3, 4, 5, 6})

        left_remainder = partitions["left"][-1]
        right_remainder = partitions["right"][-1]

        assert left_remainder == "_rest_of_left"
        assert right_remainder == "_rest_of_right"
        assert left_remainder != right_remainder
        assert set(groups[left_remainder]) == {3, 4, 5, 6}
        assert set(groups[right_remainder]) == {1, 2, 3, 4}

    def test_repeated_calls_are_idempotent(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._plotters import _load_partition_definitions

        settings = ContactsSettings(
            protein_groups={"nterm": [1, 2], "cterm": [5, 6]},
            protein_partitions={"left": ["nterm"], "right": ["cterm"]},
        )
        data = {"__meta__": {"settings": settings}}

        groups_1, partitions_1 = _load_partition_definitions(data, all_resids={1, 2, 3, 4, 5, 6})
        groups_2, partitions_2 = _load_partition_definitions(data, all_resids={1, 2, 3, 4, 5, 6})

        assert groups_1 == groups_2
        assert partitions_1 == partitions_2


# ---------------------------------------------------------------------------
# Binding preference sub-pipeline
# ---------------------------------------------------------------------------


class TestBindingPreference:
    """Test _load_or_compute_binding_preference and helpers."""

    def test_returns_none_when_disabled_and_no_cached(self, tmp_path):
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings(compute_binding_preference=False)

        mock_sim = MagicMock()
        conditions = [
            Condition(
                label=label,
                config_path=Path(f"/tmp/{label}/config.yaml"),
                replicates=(1, 2),
                sim_config=mock_sim,
            )
            for label in ["A", "B"]
        ]
        analysis_dirs = {
            "A": tmp_path / "A" / "contacts",
            "B": tmp_path / "B" / "contacts",
        }
        for d in analysis_dirs.values():
            d.mkdir(parents=True, exist_ok=True)

        ctx = ComparisonContext(
            name="test",
            conditions=conditions,
            excluded_conditions=[],
            control_label=None,
            analysis_dirs=analysis_dirs,
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=settings,
            recompute=False,
        )

        condition_data = [(cond, {"agg_result": MagicMock()}) for cond in conditions]

        result = analysis._load_or_compute_binding_preference(ctx, condition_data)
        assert result is None

    def test_try_load_cached_bp_not_found(self, tmp_path):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.shared.binding_preference_helpers import (
            try_load_cached_binding_preference,
        )

        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=MagicMock(),
        )

        # Empty directory — no cached results
        result = try_load_cached_binding_preference(cond, tmp_path)
        assert result is None

    def test_optional_bp_uses_contacts_domain_fingerprint_for_contacts_sidecar(self, tmp_path):
        """Contacts BP should resolve upstream contacts by contacts-domain identity."""

        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint
        from polyzymd.analyses.shared.config_hash import settings_fingerprint

        analysis = ContactsAnalysis()
        settings = ContactsSettings(
            compute_binding_preference=True,
            cutoff=6.0,
            grouping="none",
            surface_exposure_threshold=0.4,
            protein_groups={"site": [1, 2]},
        )
        contacts_fp = contacts_settings_fingerprint(settings)
        assert contacts_fp != settings_fingerprint(settings)

        condition = Condition(
            label="A",
            config_path=tmp_path / "A" / "config.yaml",
            replicates=(1,),
            sim_config=MagicMock(),
        )
        analysis_dir = tmp_path / "A" / "contacts"
        analysis_dir.mkdir(parents=True)
        sidecar = analysis_dir / f"contacts_eq10ns_cut6.0_s{contacts_fp}_rep1.json"
        sidecar.write_text(f'{{"contacts_settings_fingerprint": "{contacts_fp}"}}')
        enzyme_pdb = tmp_path / "enzyme.pdb"
        enzyme_pdb.write_text("ATOM\n")
        computed = SimpleNamespace(surface_exposure_threshold=settings.surface_exposure_threshold)
        ctx = ComparisonContext(
            name="test",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={condition.label: analysis_dir},
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=settings,
            recompute=True,
        )

        with (
            patch(
                "polyzymd.analyses.shared.binding_preference_helpers.resolve_enzyme_pdb",
                return_value=enzyme_pdb,
            ),
            patch(
                "polyzymd.analyses.shared.binding_preference.compute_condition_binding_preference",
                return_value=computed,
            ) as mock_compute,
        ):
            result = analysis._load_or_compute_binding_preference(
                ctx,
                [(condition, {"agg_result": MagicMock()})],
            )

        assert result is not None
        assert mock_compute.call_args.kwargs["contact_results_by_replicate"] == {1: sidecar}
        assert mock_compute.call_args.kwargs["contact_settings_fp"] == contacts_fp


# ---------------------------------------------------------------------------
# Binding preference pairwise p-values
# ---------------------------------------------------------------------------


class TestBPPairwisePValues:
    """Test _compute_bp_pairwise_pvalues."""

    def test_returns_empty_with_one_condition(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        result = ContactsAnalysis._compute_bp_pairwise_pvalues({"A": [1.0, 2.0, 3.0]})
        assert result == {}

    def test_returns_pvalues_with_two_conditions(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        data = {
            "A": [1.0, 1.1, 1.2, 0.9, 1.05],
            "B": [2.0, 2.1, 2.2, 1.9, 2.05],
        }
        result = ContactsAnalysis._compute_bp_pairwise_pvalues(data)
        assert "A_vs_B" in result
        assert 0.0 <= result["A_vs_B"] <= 1.0

    def test_skips_insufficient_data(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        data = {
            "A": [1.0],  # Only 1 value — not enough for t-test
            "B": [2.0, 2.1, 2.2],
        }
        result = ContactsAnalysis._compute_bp_pairwise_pvalues(data)
        assert result == {}  # Skipped due to insufficient data


# ---------------------------------------------------------------------------
# Integration: full lifecycle
# ---------------------------------------------------------------------------


class TestLifecycle:
    """Verify the plugin can be instantiated and key methods are callable."""

    def test_instantiation(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        assert analysis.name == "contacts"
        assert repr(analysis) == "<ContactsAnalysis(name='contacts')>"

    def test_extract_metrics_returns_empty(self):
        """contacts overrides compare() entirely, so extract_metrics() returns {}."""
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        result = analysis.extract_metrics(MagicMock())
        assert result == {}

    def test_is_subclass_of_analysis(self):
        from polyzymd.analyses.base import Analysis
        from polyzymd.analyses.contacts import ContactsAnalysis

        assert issubclass(ContactsAnalysis, Analysis)

    def test_has_all_required_methods(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        assert callable(analysis.run_replicate)
        assert callable(analysis.aggregate)
        assert callable(analysis.compare)
        assert callable(analysis.plot)
        assert callable(analysis.filter_conditions)


class TestContactsCacheAmbiguity:
    """Tests for ambiguous cache file detection in contacts path resolution."""

    def test_canonical_run_result_precedes_legacy_cache(self, tmp_path):
        """Canonical ``run_N/result.json`` should win over legacy sidecars."""
        from polyzymd.analyses.contacts._paths import find_contact_result_for_replicate

        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        canonical = run_dir / "result.json"
        canonical.write_text("{}")
        (tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json").write_text("{}")

        assert find_contact_result_for_replicate(tmp_path, 1) == canonical

    def test_find_results_for_replicates_includes_canonical_paths(self, tmp_path):
        """Replicate mapping should expose canonical result paths for consumers."""
        from polyzymd.analyses.contacts._paths import find_contact_results_for_replicates

        run_1 = tmp_path / "run_1"
        run_2 = tmp_path / "run_2"
        run_1.mkdir()
        run_2.mkdir()
        (run_1 / "result.json").write_text("{}")
        (run_2 / "result.json").write_text("{}")

        found = find_contact_results_for_replicates(tmp_path, (1, 2))

        assert found == {1: run_1 / "result.json", 2: run_2 / "result.json"}

    def test_raises_on_ambiguous_top_level_glob(self, tmp_path):
        """Multiple contacts_eq*_cut*_rep1.json files should raise ValueError."""
        from polyzymd.analyses.contacts._paths import find_contact_result_for_replicate

        (tmp_path / "contacts_eq10ns_cut4.0_rep1.json").write_text("{}")
        (tmp_path / "contacts_eq10ns_cut4.5_rep1.json").write_text("{}")

        with pytest.raises(ValueError, match="Ambiguous contacts cache"):
            find_contact_result_for_replicate(tmp_path, 1)

    def test_raises_on_ambiguous_run_subdir_glob(self, tmp_path):
        """Multiple matches in run_{rep}/ subdir should raise ValueError."""
        from polyzymd.analyses.contacts._paths import find_contact_result_for_replicate

        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        (run_dir / "contacts_eq5ns_cut4.0_rep1.json").write_text("{}")
        (run_dir / "contacts_eq10ns_cut4.5_rep1.json").write_text("{}")

        with pytest.raises(ValueError, match="Ambiguous contacts cache"):
            find_contact_result_for_replicate(tmp_path, 1)

    def test_single_glob_match_succeeds(self, tmp_path):
        """A single glob match should return that file without error."""
        from polyzymd.analyses.contacts._paths import find_contact_result_for_replicate

        (tmp_path / "contacts_eq10ns_cut4.5_rep1.json").write_text("{}")
        result = find_contact_result_for_replicate(tmp_path, 1)

        assert result is not None
        assert result.name == "contacts_eq10ns_cut4.5_rep1.json"

    def test_fingerprinted_cache_precedes_legacy_when_metadata_proves_settings(
        self,
        tmp_path,
    ):
        """Fingerprinted cache should be preferred only with embedded settings proof."""
        from polyzymd.analyses.contacts._paths import find_contact_result_for_replicate

        legacy = tmp_path / "contacts_eq10ns_cut4.5_rep1.json"
        fingerprinted = tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json"
        legacy.write_text("{}")
        fingerprinted.write_text('{"contacts_settings_fingerprint": "deadbeef"}')

        result = find_contact_result_for_replicate(tmp_path, 1, settings_fp="deadbeef")

        assert result == fingerprinted

    def test_strict_lookup_rejects_filename_only_fingerprinted_sidecar(self, tmp_path):
        """Filename-only ``_s`` tokens are candidate locators, not identity proof."""
        from polyzymd.analyses.contacts._paths import find_contact_result_for_replicate

        sidecar = tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json"
        sidecar.write_text("{}")

        result = find_contact_result_for_replicate(tmp_path, 1, settings_fp="deadbeef")

        assert result is None

    def test_strict_lookup_accepts_sidecar_with_embedded_contacts_metadata(self, tmp_path):
        """Strict sidecar lookup should accept explicit contacts-domain metadata."""
        from polyzymd.analyses.contacts._paths import find_contact_result_for_replicate

        sidecar = tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json"
        sidecar.write_text('{"contacts_settings_fingerprint": "deadbeef"}')

        result = find_contact_result_for_replicate(tmp_path, 1, settings_fp="deadbeef")

        assert result == sidecar

    @pytest.mark.parametrize(
        ("relative_name", "expected_relative_name"),
        [
            ("contacts_eq10ns_cut4.5_rep1.json", "contacts_eq10ns_cut4.5_rep1.json"),
            ("run_1/contacts_eq10ns_cut4.5_rep1.json", "run_1/contacts_eq10ns_cut4.5_rep1.json"),
            ("contacts_rep1.json", "contacts_rep1.json"),
            ("run_1/contacts_rep1.json", "run_1/contacts_rep1.json"),
        ],
    )
    def test_strict_lookup_accepts_metadata_proven_legacy_sidecar(
        self,
        tmp_path,
        relative_name,
        expected_relative_name,
    ):
        """Explicit lookups may use legacy contacts sidecars proven by metadata."""
        from polyzymd.analyses.contacts._paths import find_contact_result_for_replicate

        legacy = tmp_path / relative_name
        legacy.parent.mkdir(parents=True, exist_ok=True)
        legacy.write_text(
            '{"metadata": {"contacts_settings_fingerprint": "deadbeef", '
            '"equilibration": "10ns"}}'
        )

        result = find_contact_result_for_replicate(
            tmp_path,
            1,
            settings_fp="deadbeef",
            equilibration="10ns",
        )

        assert result == tmp_path / expected_relative_name

    @pytest.mark.parametrize(
        "fingerprinted_payload",
        [
            "{}",
            '{"metadata": {"contacts_settings_fingerprint": "badcafe0"}}',
        ],
    )
    def test_strict_lookup_blocks_legacy_when_fingerprinted_sidecar_unproven_or_mismatched(
        self,
        tmp_path,
        fingerprinted_payload,
    ):
        """Any unproven or mismatched ``_s*`` sidecar blocks legacy fallback."""
        from polyzymd.analyses.contacts._paths import find_contact_result_for_replicate

        fingerprinted = tmp_path / "contacts_eq10ns_cut4.5_sbadcafe0_rep1.json"
        legacy = tmp_path / "contacts_eq10ns_cut4.5_rep1.json"
        fingerprinted.write_text(fingerprinted_payload)
        legacy.write_text(
            '{"metadata": {"contacts_settings_fingerprint": "deadbeef", '
            '"equilibration": "10ns"}}'
        )

        result = find_contact_result_for_replicate(
            tmp_path,
            1,
            settings_fp="deadbeef",
            equilibration="10ns",
        )

        assert result is None

    def test_lookup_filters_contacts_by_equilibration_window(self, tmp_path):
        """Contacts artifact lookup should not cross equilibration windows."""
        from polyzymd.analyses.contacts._paths import find_contact_result_for_replicate

        wrong = tmp_path / "contacts_eq0ns_cut4.5_sdeadbeef_rep1.json"
        expected = tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json"
        wrong.write_text("{}")
        expected.write_text('{"contacts_settings_fingerprint": "deadbeef"}')

        result = find_contact_result_for_replicate(
            tmp_path,
            1,
            settings_fp="deadbeef",
            equilibration="10ns",
        )

        assert result == expected

    def test_lookup_rejects_canonical_with_wrong_window_metadata(self, tmp_path):
        """Canonical contacts outputs must prove the requested window."""
        from polyzymd.analyses.contacts._paths import find_contact_result_for_replicate

        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        canonical = run_dir / "result.json"
        canonical.write_text('{"equilibration_time": 0.0, "equilibration_unit": "ns"}')

        assert find_contact_result_for_replicate(tmp_path, 1, equilibration="10ns") is None

    def test_fingerprinted_cache_precedes_canonical_when_settings_fp_provided(self, tmp_path):
        """Matching fingerprinted cache should win over potentially stale canonical output."""
        from polyzymd.analyses.contacts._paths import find_contact_result_for_replicate

        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        canonical = run_dir / "result.json"
        fingerprinted = tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json"
        canonical.write_text("{}")
        fingerprinted.write_text('{"contacts_settings_fingerprint": "deadbeef"}')

        result = find_contact_result_for_replicate(tmp_path, 1, settings_fp="deadbeef")

        assert result == fingerprinted

    def test_strict_lookup_accepts_canonical_only_with_matching_fingerprint(self, tmp_path):
        """Explicit settings lookups may use canonical only when metadata proves identity."""
        from polyzymd.analyses.contacts._paths import find_contact_result_for_replicate

        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        canonical = run_dir / "result.json"
        canonical.write_text(
            '{"metadata": {"settings_fingerprint": "deadbeef", "equilibration": "10ns"}}'
        )

        result = find_contact_result_for_replicate(
            tmp_path,
            1,
            settings_fp="deadbeef",
            equilibration="10ns",
        )

        assert result == canonical

    def test_strict_lookup_rejects_window_only_canonical(self, tmp_path):
        """Explicit settings lookups should not fall back to unproven canonical files."""
        from polyzymd.analyses.contacts._paths import find_contact_result_for_replicate

        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        canonical = run_dir / "result.json"
        canonical.write_text('{"metadata": {"equilibration": "10ns"}}')

        result = find_contact_result_for_replicate(
            tmp_path,
            1,
            settings_fp="deadbeef",
            equilibration="10ns",
        )

        assert result is None

    def test_mismatched_fingerprinted_cache_rejected_when_settings_fp_provided(self, tmp_path):
        """Non-matching fingerprinted sidecars should not satisfy explicit lookups."""
        from polyzymd.analyses.contacts._paths import find_contact_result_for_replicate

        mismatched = tmp_path / "contacts_eq10ns_cut4.5_sfeedface_rep1.json"
        mismatched.write_text("{}")

        result = find_contact_result_for_replicate(tmp_path, 1, settings_fp="deadbeef")

        assert result is None

    def test_infers_settings_fingerprint_from_contact_metadata(self, tmp_path):
        """Downstream consumers should use actual contacts artifact metadata."""
        from polyzymd.analyses.contacts._paths import infer_contact_results_settings_fingerprint

        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        canonical = run_dir / "result.json"
        canonical.write_text('{"metadata": {"settings_fingerprint": "cafebabe"}}')

        assert infer_contact_results_settings_fingerprint({1: canonical}) == "cafebabe"

    def test_infer_settings_fingerprint_rejects_conflicting_artifacts(self, tmp_path):
        """Mixed contacts contracts should not be silently accepted."""
        from polyzymd.analyses.contacts._paths import infer_contact_results_settings_fingerprint

        first = tmp_path / "contacts_eq10ns_cut4.5_s11111111_rep1.json"
        second = tmp_path / "contacts_eq10ns_cut4.5_s22222222_rep2.json"
        first.write_text('{"metadata": {"settings_fingerprint": "11111111"}}')
        second.write_text('{"metadata": {"settings_fingerprint": "22222222"}}')

        with pytest.raises(ValueError, match="inconsistent settings fingerprints"):
            infer_contact_results_settings_fingerprint({1: first, 2: second})

    def test_artifact_inference_prefers_fingerprinted_sidecar_over_stale_canonical(
        self,
        tmp_path,
    ):
        """Current fingerprinted contacts sidecars should define cache identity."""
        from polyzymd.analyses.contacts._paths import infer_contacts_artifact_settings_fingerprint

        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        (run_dir / "result.json").write_text(
            "{"
            '"criteria_cutoff": 4.0,'
            '"criteria_label": "Distance <= 4.0 A",'
            '"n_frames": 10,'
            '"residue_contacts": [],'
            '"metadata": {"settings_fingerprint": "badcafe0", "equilibration": "10ns"}'
            "}"
        )
        (tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json").write_text(
            '{"metadata": {"contacts_settings_fingerprint": "deadbeef"}}'
        )

        result = infer_contacts_artifact_settings_fingerprint(
            tmp_path,
            (1,),
            equilibration="10ns",
        )

        assert result == "deadbeef"

    def test_contact_artifact_match_requires_embedded_metadata(self, tmp_path):
        """Strict artifact matching should ignore filename-only settings tokens."""
        from polyzymd.analyses.contacts._paths import (
            contact_artifact_matches_settings_fingerprint,
        )

        filename_only = tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json"
        embedded = tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep2.json"
        filename_only.write_text("{}")
        embedded.write_text('{"metadata": {"contacts_settings_fingerprint": "deadbeef"}}')

        assert not contact_artifact_matches_settings_fingerprint(filename_only, "deadbeef")
        assert contact_artifact_matches_settings_fingerprint(embedded, "deadbeef")

    def test_infer_contacts_artifact_settings_ignores_filename_only_sidecars(
        self,
        tmp_path,
    ):
        """Contacts identity inference should fail closed for filename-only sidecars."""
        from polyzymd.analyses.contacts._paths import infer_contacts_artifact_settings_fingerprint

        (tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json").write_text("{}")

        result = infer_contacts_artifact_settings_fingerprint(
            tmp_path,
            (1,),
            equilibration="10ns",
        )

        assert result is None

    def test_infer_contacts_artifact_settings_uses_metadata_proven_legacy_sidecars(
        self,
        tmp_path,
    ):
        """Legacy contacts sidecars with metadata should define contacts identity."""
        from polyzymd.analyses.contacts._paths import infer_contacts_artifact_settings_fingerprint

        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        (run_dir / "result.json").write_text('{"metadata": {"equilibration": "10ns"}}')
        (tmp_path / "contacts_eq10ns_cut6.0_rep1.json").write_text(
            '{"metadata": {"contacts_settings_fingerprint": "deadbeef", '
            '"equilibration": "10ns"}}'
        )

        result = infer_contacts_artifact_settings_fingerprint(
            tmp_path,
            (1,),
            equilibration="10ns",
        )

        assert result == "deadbeef"

    def test_filename_only_sidecar_blocks_legacy_artifact_inference(self, tmp_path):
        """Unproven fingerprinted contacts sidecars should block legacy inference."""
        from polyzymd.analyses.contacts._paths import infer_contacts_artifact_settings_fingerprint

        (tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json").write_text("{}")
        (tmp_path / "contacts_eq10ns_cut6.0_rep1.json").write_text(
            '{"metadata": {"contacts_settings_fingerprint": "feedface", '
            '"equilibration": "10ns"}}'
        )

        result = infer_contacts_artifact_settings_fingerprint(
            tmp_path,
            (1,),
            equilibration="10ns",
        )

        assert result is None

    def test_strict_identity_blocks_canonical_when_fingerprinted_sidecar_unproven(
        self,
        tmp_path,
    ):
        """Strict downstream lookup must fail closed on filename-only sidecars."""
        from polyzymd.analyses.contacts._paths import find_contact_result_for_replicate

        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        canonical = run_dir / "result.json"
        canonical.write_text('{"metadata": {"equilibration": "10ns"}}')
        (tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json").write_text("{}")

        result = find_contact_result_for_replicate(
            tmp_path,
            1,
            equilibration="10ns",
            strict_identity=True,
        )

        assert result is None

    def test_strict_identity_accepts_proven_fingerprinted_sidecar_without_requested_fp(
        self,
        tmp_path,
    ):
        """Strict downstream lookup can use sidecars with embedded contacts identity."""
        from polyzymd.analyses.contacts._paths import find_contact_result_for_replicate

        sidecar = tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json"
        sidecar.write_text('{"metadata": {"contacts_settings_fingerprint": "deadbeef"}}')

        result = find_contact_result_for_replicate(
            tmp_path,
            1,
            equilibration="10ns",
            strict_identity=True,
        )

        assert result == sidecar

    def test_infer_contacts_artifact_settings_stops_before_bp_when_sidecar_unproven(
        self,
        tmp_path,
    ):
        """Filename-only contacts sidecars are an identity barrier for BP fallback."""
        from polyzymd.analyses.contacts._paths import infer_contacts_artifact_settings_fingerprint

        (tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json").write_text("{}")
        (tmp_path / "binding_preference_sfeedface_rep1.json").write_text(
            '{"metadata": {"contacts_settings_fingerprint": "feedface", "equilibration": "10ns"}}'
        )

        result = infer_contacts_artifact_settings_fingerprint(
            tmp_path,
            (1,),
            equilibration="10ns",
        )

        assert result is None

    def test_binding_preference_generic_settings_do_not_infer_contacts_identity(
        self,
        tmp_path,
    ):
        """BP settings fingerprints must not masquerade as contacts fingerprints."""
        from polyzymd.analyses.contacts._paths import infer_contacts_artifact_settings_fingerprint

        generic = tmp_path / "binding_preference_sdeadbeef_rep1.json"
        generic.write_text(
            "{"
            '"settings_fingerprint": "badcafe0",'
            '"settings_fp": "feedface",'
            '"metadata": {"settings_fingerprint": "cafebabe"}'
            "}"
        )

        result = infer_contacts_artifact_settings_fingerprint(tmp_path, (1,))

        assert result is None

    def test_binding_preference_contacts_metadata_infers_contacts_identity(self, tmp_path):
        """BP artifacts may declare contacts identity with explicit contacts keys."""
        from polyzymd.analyses.contacts._paths import infer_contacts_artifact_settings_fingerprint

        artifact = tmp_path / "binding_preference_sbadcafe_rep1.json"
        artifact.write_text('{"metadata": {"contacts_settings_fingerprint": "deadbeef"}}')

        result = infer_contacts_artifact_settings_fingerprint(tmp_path, (1,))

        assert result == "deadbeef"

    def test_aggregated_sidecar_names_noncontiguous_replicates_uniquely(self, tmp_path):
        """Non-contiguous replicate IDs should not collide with contiguous ranges."""
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        settings = ContactsSettings()

        sparse = ContactsAnalysis._aggregated_sidecar_path(tmp_path, settings, "10ns", (1, 3))
        dense = ContactsAnalysis._aggregated_sidecar_path(tmp_path, settings, "10ns", (1, 2, 3))

        assert sparse.name != dense.name
        assert "_reps1_3.json" in sparse.name
        assert "_reps1-3.json" in dense.name
