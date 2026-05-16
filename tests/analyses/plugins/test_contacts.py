"""Tests for the contacts analysis plugin.

Covers discovery, settings, run_replicate, aggregate, compare (full override),
filter_conditions, plot delegation, AggregatedResultClass, and per-replicate
metric helpers.
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

        assert ContactsAnalysis.min_replicates == 1


def test_singleton_contacts_pairwise_not_testable() -> None:
    """Contacts singleton pairwise metrics should be marked not testable."""
    from polyzymd.analyses.contacts._comparison import compare_contacts_pair

    summary_a = SimpleNamespace(
        coverage_mean=0.4,
        coverage_sem=0.0,
        mean_contact_fraction=0.1,
        mean_contact_fraction_sem=0.0,
    )
    summary_b = SimpleNamespace(
        coverage_mean=0.6,
        coverage_sem=0.0,
        mean_contact_fraction=0.2,
        mean_contact_fraction_sem=0.0,
    )

    result = compare_contacts_pair(
        "A",
        summary_a,
        {"coverage_per_replicate": [0.4], "contact_fraction_per_replicate": [0.1]},
        "B",
        summary_b,
        {"coverage_per_replicate": [0.6], "contact_fraction_per_replicate": [0.2]},
    )

    assert all(comparison.testable is False for comparison in result.aggregate_comparisons)
    assert all(comparison.significant is False for comparison in result.aggregate_comparisons)
    assert all(comparison.p_value_adjusted is None for comparison in result.aggregate_comparisons)


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class TestSettings:
    """Validate ContactsSettings."""

    def test_defaults(self):
        from polyzymd.analyses.contacts import ContactsSettings

        s = ContactsSettings()
        assert s.polymer_selection == "chainid C"
        assert s.protein_selection == "chainid A"
        assert s.cutoff == 4.5
        assert s.grouping == "aa_class"
        assert s.compute_residence_times is True
        assert s.fdr_alpha == 0.05

    def test_custom_values(self):
        from polyzymd.analyses.contacts import ContactsSettings

        s = ContactsSettings(
            polymer_selection="chainid D",
            protein_selection="chainid A and name CA",
            cutoff=5.0,
            grouping="none",
            protein_groups={"active": [1, 2]},
            protein_partitions={"site": ["active"]},
            fdr_alpha=0.01,
        )
        assert s.polymer_selection == "chainid D"
        assert s.cutoff == 5.0
        assert s.grouping == "none"
        assert s.protein_groups == {"active": [1, 2]}
        assert s.protein_partitions == {"site": ["active"]}
        assert s.fdr_alpha == 0.01

    @pytest.mark.parametrize(
        "key,value",
        [
            ("compute_binding_preference", True),
            ("surface_exposure_threshold", 0.3),
            ("enzyme_pdb_for_sasa", "enzyme.pdb"),
            ("include_default_aa_groups", False),
            ("polymer_type_selections", {"PEG": "resname PEG"}),
            ("polymer_chain", "C"),
            ("enrichment_normalization", "residue"),
        ],
    )
    def test_archived_binding_preference_settings_rejected(self, key, value):
        from polyzymd.analyses.contacts import ContactsSettings

        with pytest.raises(ValueError, match="archive_experimental_analysis") as excinfo:
            ContactsSettings.model_validate({key: value})
        assert "feature/mda-analysis-migration" in str(excinfo.value)

    @pytest.mark.parametrize(
        "key,value",
        [
            ("generate_enrichment_heatmap", True),
            ("generate_enrichment_bars", True),
            ("figsize_enrichment_heatmap", (10, 6)),
            ("figsize_enrichment_bars", (10, 6)),
            ("enrichment_colormap", "RdBu_r"),
            ("show_enrichment_error", True),
            ("generate_system_coverage_heatmap", True),
            ("generate_system_coverage_bars", True),
            ("figsize_system_coverage_heatmap", (10, 6)),
            ("figsize_system_coverage_bars", (10, 6)),
            ("show_system_coverage_error", True),
            ("generate_user_partition_bars", True),
            ("figsize_user_partition_bars", (10, 6)),
            ("show_user_partition_error", True),
        ],
    )
    def test_archived_binding_preference_plot_settings_rejected(self, key, value):
        from polyzymd.analyses.contacts._plot_settings import ContactsPlotSettings

        with pytest.raises(ValueError, match="archive_experimental_analysis") as excinfo:
            ContactsPlotSettings.model_validate({key: value})
        assert "feature/mda-analysis-migration" in str(excinfo.value)

    def test_retained_contacts_plot_settings_still_validate(self):
        from polyzymd.analyses.contacts._plot_settings import ContactsPlotSettings

        settings = ContactsPlotSettings(
            generate_contact_fraction_profile=True,
            generate_residence_time_profile=True,
            generate_cf_by_aa_class_bars=True,
            generate_cf_by_partition_bars=True,
            generate_rt_by_aa_class_bars=True,
            generate_rt_by_partition_bars=True,
            highlight_residues=[1, 2, 3],
        )

        assert settings.highlight_residues == [1, 2, 3]

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
            protein_groups={"active": [1, 2]},
            protein_partitions={"site": ["active"]},
        )
        d = s.model_dump()
        s2 = ContactsSettings.model_validate(d)
        assert s2.cutoff == 5.0
        assert s2.polymer_types == ["SBM", "EGM"]
        assert s2.protein_groups == {"active": [1, 2]}
        assert s2.protein_partitions == {"site": ["active"]}

    def test_settings_fingerprint_changes_with_cutoff(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint

        low = ContactsSettings(cutoff=4.0)
        high = ContactsSettings(cutoff=4.5)

        assert contacts_settings_fingerprint(low) != contacts_settings_fingerprint(high)

    def test_contacts_domain_fingerprint_excludes_comparison_partition_fields(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint

        base = ContactsSettings(cutoff=4.5)
        downstream_only = ContactsSettings(
            cutoff=4.5,
            protein_groups={"active": [1, 2]},
            protein_partitions={"site": ["active"]},
            top_residues=25,
        )

        assert contacts_settings_fingerprint(base) == contacts_settings_fingerprint(downstream_only)

    def test_contacts_domain_fingerprint_changes_with_contact_fields(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint

        unfiltered = ContactsSettings(polymer_selection="chainid C")
        filtered = ContactsSettings(polymer_selection="chainid C", polymer_types=["PEG"])

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


# ---------------------------------------------------------------------------
# MDAnalysis sparse contact detection
# ---------------------------------------------------------------------------


def _make_contact_universe(coords=None, dimensions=None):
    """Build a tiny MDAnalysis universe for sparse contact tests."""

    import MDAnalysis as mda
    import numpy as np

    if coords is None:
        coords = np.asarray(
            [
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
                [[0.0, 0.0, 0.0], [8.0, 0.0, 0.0]],
            ],
            dtype=np.float32,
        )
    if dimensions is None:
        dimensions = np.tile(np.asarray([20.0, 20.0, 20.0, 90.0, 90.0, 90.0]), (len(coords), 1))
    universe = mda.Universe.empty(
        2,
        n_residues=2,
        atom_resindex=[0, 1],
        n_segments=2,
        residue_segindex=[0, 1],
        trajectory=True,
    )
    universe.add_TopologyAttr("names", ["CA", "C1"])
    universe.add_TopologyAttr("types", ["C", "C"])
    universe.add_TopologyAttr("resids", [1, 10])
    universe.add_TopologyAttr("resnames", ["ALA", "PEG"])
    universe.add_TopologyAttr("chainIDs", ["A", "C"])
    universe.add_TopologyAttr("segids", ["A", "C"])
    universe.load_new(coords, order="fac", dimensions=dimensions)
    return universe


def _make_contacts_collector_context(tmp_path, settings):
    """Build a collector context for direct contact artifact tests."""

    from polyzymd.analyses.base import Condition, ReplicateContext
    from polyzymd.analyses.mda import ArtifactStore, FrameSelection, MDACollectorContext
    from polyzymd.analyses.mda.job import MDAUniversePolicy

    sim_config = _make_hashable_sim_config(tmp_path)
    condition = Condition(
        label="test",
        config_path=tmp_path / "config.yaml",
        replicates=(1,),
        sim_config=sim_config,
    )
    replicate_ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=sim_config,
        output_dir=tmp_path,
        equilibration="0ns",
        recompute=True,
        settings=settings,
        result_path=tmp_path / "result.json",
    )
    return MDACollectorContext(
        analysis_name="contacts",
        replicate_context=replicate_ctx,
        frame_selection=FrameSelection(start=0, stop=2, step=1, timestep_ps=1.0, n_frames_total=2),
        universe_policy=MDAUniversePolicy(condition_label="test", replicate=1),
        artifact_store=ArtifactStore(tmp_path),
        settings_fingerprint="contacts-test-fp",
    )


class TestContactsMDAArtifacts:
    """Contacts detection uses sparse MDAnalysis event artifacts."""

    def test_plugin_uses_mda_artifacts_without_runner_hook(self):
        from polyzymd.analyses.base import Analysis
        from polyzymd.analyses.contacts import ContactsAnalysis

        assert ContactsAnalysis.ReplicateResultClass is None
        assert ContactsAnalysis.build_runner is Analysis.build_runner
        assert ContactsAnalysis.summarize_replicate is Analysis.summarize_replicate

    def test_sparse_analysis_streams_contact_events(self):
        from polyzymd.analyses.contacts._mda import build_contact_event_analysis

        analysis = build_contact_event_analysis(
            universe=_make_contact_universe(),
            protein_selection="chainid A",
            polymer_selection="chainid C",
            cutoff=4.5,
            grouping_mode="aa_class",
            raw_timestep_ps=1.0,
        )
        analysis.run(start=0, stop=2, step=1)

        assert analysis.results.event_start_sample_index.tolist() == [0]
        assert analysis.results.event_start_frame.tolist() == [0]
        assert analysis.results.event_duration_samples.tolist() == [1]
        assert analysis.results.event_duration_ps.tolist() == [1.0]
        assert analysis.results.protein_residue_index.tolist() == [0]
        assert analysis.results.polymer_residue_index.tolist() == [0]
        assert analysis.results.contact_sample_counts.tolist() == [1]
        assert not hasattr(analysis.results, "contact_matrix")

    def test_capped_distance_receives_timestep_box(self):
        import numpy as np

        from polyzymd.analyses.contacts._mda import build_contact_event_analysis

        coords = np.asarray([[[0.2, 0.0, 0.0], [9.8, 0.0, 0.0]]], dtype=np.float32)
        dims = np.asarray([[10.0, 10.0, 10.0, 90.0, 90.0, 90.0]])
        captured_boxes = []

        def _fake_capped_distance(*args, **kwargs):
            captured_boxes.append(np.asarray(kwargs["box"], dtype=float).copy())
            return np.asarray([[0, 0]], dtype=np.int64), np.asarray([0.4], dtype=float)

        analysis = build_contact_event_analysis(
            universe=_make_contact_universe(coords=coords, dimensions=dims),
            protein_selection="chainid A",
            polymer_selection="chainid C",
            cutoff=0.5,
            grouping_mode="none",
            raw_timestep_ps=1.0,
        )
        with patch("MDAnalysis.lib.distances.capped_distance", side_effect=_fake_capped_distance):
            analysis.run(start=0, stop=1, step=1)

        assert captured_boxes
        np.testing.assert_allclose(captured_boxes[0], dims[0])

    def test_zero_contact_output_keeps_identity_arrays(self):
        from polyzymd.analyses.contacts._mda import build_contact_event_analysis

        analysis = build_contact_event_analysis(
            universe=_make_contact_universe(),
            protein_selection="chainid A",
            polymer_selection="chainid C",
            cutoff=0.1,
            grouping_mode="none",
            raw_timestep_ps=1.0,
        )
        analysis.run(start=0, stop=2, step=1)

        assert analysis.results.event_duration_samples.size == 0
        assert analysis.results.n_frames_used == 2
        assert [row.resid for row in analysis.results.protein_identity] == [1]
        assert [row.resid for row in analysis.results.polymer_identity] == [10]

    def test_collector_writes_summary_json_and_npz_sidecar(self, tmp_path):
        import numpy as np

        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._mda import (
            ContactsArtifactCollector,
            build_contact_event_analysis,
        )
        from polyzymd.analyses.mda import FrameSelection, MDAAnalysisJob

        settings = ContactsSettings()
        analysis = build_contact_event_analysis(
            universe=_make_contact_universe(),
            protein_selection="chainid A",
            polymer_selection="chainid C",
            cutoff=4.5,
            grouping_mode="aa_class",
            raw_timestep_ps=1.0,
        )
        job = MDAAnalysisJob(
            name="contacts",
            analysis=analysis,
            frame_selection=FrameSelection(
                start=0, stop=2, step=1, timestep_ps=1.0, n_frames_total=2
            ),
        )
        completed = job.run()
        artifact = ContactsArtifactCollector()(
            _make_contacts_collector_context(tmp_path, settings), [completed]
        )

        assert artifact.payload["n_contact_events"] == 1
        assert artifact.payload["n_protein_residues"] == 1
        assert artifact.payload["metrics"] == {"coverage": 1.0, "mean_contact_fraction": 0.5}
        assert artifact.payload["event_sidecar"] == "sidecars/contact_events.npz"
        assert "event_start_sample_index" not in artifact.model_dump(mode="json")["payload"]
        with np.load(tmp_path / artifact.payload["event_sidecar"]) as data:
            assert data["event_duration_samples"].tolist() == [1]
            assert data["event_duration_ns"].tolist() == [0.001]
            assert data["protein_resids"].tolist() == [1]
            assert data["polymer_resids"].tolist() == [10]

    def test_artifact_adapter_preserves_legacy_contact_semantics(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._mda import (
            ContactsArtifactCollector,
            artifact_to_contact_result,
            build_contact_event_analysis,
        )
        from polyzymd.analyses.mda import FrameSelection, MDAAnalysisJob

        settings = ContactsSettings()
        analysis = build_contact_event_analysis(
            universe=_make_contact_universe(),
            protein_selection="chainid A",
            polymer_selection="chainid C",
            cutoff=4.5,
            grouping_mode="aa_class",
            raw_timestep_ps=1.0,
        )
        job = MDAAnalysisJob(
            name="contacts",
            analysis=analysis,
            frame_selection=FrameSelection(
                start=0, stop=2, step=1, timestep_ps=1.0, n_frames_total=2
            ),
        )
        artifact = ContactsArtifactCollector()(
            _make_contacts_collector_context(tmp_path, settings), [job.run()]
        )

        result = artifact_to_contact_result(
            artifact, run_dir=tmp_path, settings_fingerprint="contacts-test-fp"
        )

        assert result.n_frames == 2
        assert result.coverage_fraction() == 1.0
        assert result.mean_contact_fraction() == 0.5
        assert result.residue_contacts[0].segment_contacts[0].events[0].start_frame == 0


class TestContactsSparseEventUtilities:
    """Contacts-local helpers preserve legacy scientific semantics."""

    def test_build_atom_to_residue_map_distinguishes_duplicate_resids_by_chain(self):
        from polyzymd.analyses.contacts._events import build_atom_to_residue_map

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
        atoms = [_Atom(residue_a), _Atom(residue_b)]

        atom_to_res = build_atom_to_residue_map(atoms, [residue_a, residue_b])

        assert atom_to_res.tolist() == [0, 1]

    def test_fragment_lookup_falls_back_without_bond_topology(self):
        from polyzymd.analyses.contacts._events import fragments_or_single

        mdanalysis_module = ModuleType("MDAnalysis")
        exceptions_module = ModuleType("MDAnalysis.exceptions")

        class NoDataError(Exception):
            pass

        exceptions_module.NoDataError = NoDataError
        mdanalysis_module.exceptions = exceptions_module

        class _AtomGroup:
            @property
            def fragments(self):
                raise NoDataError("No bond information")

        warnings = []
        with patch.dict(
            "sys.modules",
            {"MDAnalysis": mdanalysis_module, "MDAnalysis.exceptions": exceptions_module},
        ):
            atoms = _AtomGroup()
            assert fragments_or_single(atoms, context="test", warnings=warnings) == [atoms]

        assert warnings == [
            "test: topology has no bond information; treating the selected polymer as one fragment"
        ]

    def test_detection_fingerprint_records_cutoff_and_detection_policy(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._identity import (
            CONTACTS_PBC_POLICY,
            contacts_detection_fingerprint,
            contacts_detection_identity_payload,
        )

        low = ContactsSettings(cutoff=4.0)
        high = ContactsSettings(cutoff=4.5)
        payload = contacts_detection_identity_payload(high)

        assert contacts_detection_fingerprint(low) != contacts_detection_fingerprint(high)
        assert payload["cutoff"] == {"value": 4.5, "units": "angstrom"}
        assert payload["pbc_policy"] == CONTACTS_PBC_POLICY
        assert payload["contact_semantics"] == "any_atom_residue_pair-v1"


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


def _write_contacts_replicate_artifact(
    analysis_dir: Path,
    settings,
    *,
    condition_label: str = "test",
    replicate: int = 1,
    equilibration: str = "10ns",
    contact_fractions: list[float] | None = None,
    event_durations_ns: list[tuple[int, str, float]] | None = None,
):
    """Write a synthetic contacts replicate artifact and event sidecar."""

    import numpy as np

    from polyzymd.analyses.contacts._identity import (
        CONTACT_SEMANTICS_VERSION,
        CONTACTS_PBC_POLICY,
        contacts_detection_fingerprint,
        contacts_detection_identity_payload,
    )
    from polyzymd.analyses.mda import ArtifactStore, ReplicateArtifact

    fractions = contact_fractions or [0.5, 0.0]
    events = event_durations_ns or []
    run_dir = analysis_dir / f"run_{replicate}"
    store = ArtifactStore(run_dir)
    protein_resids = np.arange(1, len(fractions) + 1, dtype=np.int64)
    protein_resnames = np.asarray(["ALA", "ASP"][: len(fractions)], dtype="U16")
    protein_groups = np.asarray(["nonpolar", "charged"][: len(fractions)], dtype="U32")
    polymer_resnames = np.asarray(["PEG"], dtype="U16")
    start_indices = np.arange(len(events), dtype=np.int64)
    durations_ns = np.asarray([event[2] for event in events], dtype=np.float64)
    protein_indices = np.asarray([event[0] for event in events], dtype=np.int64)
    sidecar = store.write_npz_sidecar(
        "sidecars/contact_events.npz",
        event_start_sample_index=start_indices,
        event_start_frame=start_indices,
        event_duration_samples=np.ones(len(events), dtype=np.int64),
        event_duration_ps=durations_ns * 1000.0,
        event_duration_ns=durations_ns,
        protein_residue_index=protein_indices,
        polymer_residue_index=np.zeros(len(events), dtype=np.int64),
        polymer_chain_index=np.zeros(len(events), dtype=np.int64),
        frame_indices=np.arange(10, dtype=np.int64),
        time_ps=np.arange(10, dtype=np.float64) * 1000.0,
        protein_resids=protein_resids,
        protein_resnames=protein_resnames,
        protein_groups=protein_groups,
        protein_chainids=np.asarray(["A"] * len(fractions), dtype="U16"),
        polymer_resids=np.asarray([10], dtype=np.int64),
        polymer_resnames=polymer_resnames,
        polymer_chain_indices=np.asarray([0], dtype=np.int64),
        polymer_chainids=np.asarray(["C"], dtype="U16"),
        metadata={"kind": "contact_events"},
    )
    detection = contacts_detection_identity_payload(settings)
    coverage = sum(value > 0 for value in fractions) / len(fractions)
    artifact = ReplicateArtifact(
        analysis_name="contacts",
        condition_label=condition_label,
        replicate=replicate,
        payload={
            "metrics": {
                "coverage": coverage,
                "mean_contact_fraction": sum(fractions) / len(fractions),
            },
            "replicate_metrics": {
                "coverage": coverage,
                "mean_contact_fraction": sum(fractions) / len(fractions),
            },
            "n_frames_total": 10,
            "n_frames_used": 10,
            "n_contact_events": len(events),
            "n_contacted_residues": sum(value > 0 for value in fractions),
            "n_protein_residues": len(fractions),
            "n_polymer_residues": 1,
            "event_sidecar": sidecar.path,
            "criteria_cutoff": settings.cutoff,
            "contact_semantics": CONTACT_SEMANTICS_VERSION,
            "pbc_policy": CONTACTS_PBC_POLICY,
            "polymer_types": ["PEG"],
            "protein_residues": [
                {
                    "protein_residue_index": index,
                    "protein_resid": index + 1,
                    "protein_resname": str(protein_resnames[index]),
                    "protein_chain_id": "A",
                    "protein_group": str(protein_groups[index]),
                    "contact_fraction": fraction,
                    "event_count": sum(event[0] == index for event in events),
                    "polymer_type_contact_fractions": {"PEG": fraction},
                }
                for index, fraction in enumerate(fractions)
            ],
        },
        sidecars=[sidecar],
        provenance={
            "frame_selection": {
                "start": 0,
                "stop": 10,
                "step": 1,
                "equilibration": equilibration,
                "n_frames_total": 10,
                "n_frames_selected": 10,
                "timestep_ps": 1000.0,
            },
            "detection_identity": detection,
            "protein_selection": detection["protein_selection"],
            "polymer_selection": detection["polymer_selection"],
            "effective_polymer_selection": detection["effective_polymer_selection"],
            "grouping": detection["grouping"],
        },
        metadata={
            "settings_fingerprint": contacts_detection_fingerprint(settings),
            "contacts_detection_fingerprint": contacts_detection_fingerprint(settings),
            "time_axis_policy": "regular_selected_time_axis",
            "equilibration": equilibration,
        },
    )
    store.write_replicate_result(artifact)
    return artifact


def _make_condition_artifact(
    label: str,
    settings,
    replicates=(1, 2, 3),
    n_residues: int = 2,
    equilibration: str = "10ns",
    compute_residence_times: bool | None = None,
    include_residence_metadata: bool = True,
):
    """Create a condition artifact for contacts comparison tests."""

    import numpy as np

    from polyzymd.analyses.contacts._identity import contacts_detection_fingerprint
    from polyzymd.analyses.mda import ConditionArtifact

    residence_enabled = (
        bool(settings.compute_residence_times)
        if compute_residence_times is None
        else bool(compute_residence_times)
    )
    coverage_values = [1.0 for _ in replicates]
    contact_values = [0.25 + 0.05 * index for index, _rep in enumerate(replicates)]

    def metric(name, values):
        mean = float(np.mean(values))
        std = 0.0 if len(values) == 1 else float(np.std(values, ddof=1))
        sem = 0.0 if len(values) == 1 else std / (len(values) ** 0.5)
        return {
            "name": name,
            "values": list(values),
            "mean": mean,
            "sem": sem,
            "std": std,
            "n": len(values),
        }

    payload = {
        "metrics": {
            "coverage": metric("coverage", coverage_values),
            "mean_contact_fraction": metric("mean_contact_fraction", contact_values),
        },
        "replicate_metrics": {
            str(rep): {"coverage": 1.0, "mean_contact_fraction": contact_values[index]}
            for index, rep in enumerate(replicates)
        },
        "n_residues": n_residues,
        "residue_stats": [
            {
                "protein_resid": index + 1,
                "protein_resname": "ALA",
                "protein_chain_id": "A",
                "contact_fraction_mean": 0.5 - index * 0.1,
            }
            for index in range(n_residues)
        ],
    }
    if residence_enabled:
        payload["residence_time_by_polymer_type"] = {
            "PEG": {
                "mean_ns": 4.0,
                "sem_ns": 1.0,
                "n_events": 2,
                "replicates_with_events": [replicates[0]],
                "replicate_means_ns": [4.0],
            }
        }
    metadata = {
        "contacts_detection_fingerprint": contacts_detection_fingerprint(settings),
        "equilibration": equilibration,
    }
    if include_residence_metadata:
        metadata["compute_residence_times"] = residence_enabled
    return ConditionArtifact(
        analysis_name="contacts",
        condition_label=label,
        replicates=list(replicates),
        payload=payload,
        provenance={
            "frame_selection": {
                "start": 0,
                "stop": 10,
                "step": 1,
                "equilibration": equilibration,
                "n_frames_total": 10,
                "n_frames_selected": 10,
                "timestep_ps": 1000.0,
            }
        },
        metadata=metadata,
    )


class TestAggregate:
    """aggregate consumes contacts replicate artifacts only."""

    def test_aggregate_artifacts_writes_condition_artifact(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.mda import ArtifactStore, ConditionArtifact
        from polyzymd.analyses.mda.aggregation import MDAAggregationError

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        analysis_dir = tmp_path / "contacts"
        output_dir = analysis_dir / "aggregated"
        sim_config = _make_hashable_sim_config(tmp_path)
        cond = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2, 3),
            sim_config=sim_config,
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1, 2, 3),
            output_dir=output_dir,
            equilibration="10ns",
            settings=settings,
        )
        artifacts = [
            _write_contacts_replicate_artifact(
                analysis_dir,
                settings,
                replicate=replicate,
                contact_fractions=[0.5, 0.0 if replicate == 3 else 0.2],
            )
            for replicate in (1, 2, 3)
        ]

        result = analysis.aggregate(ctx, artifacts)

        assert isinstance(result, ConditionArtifact)
        assert (output_dir / "result.json").exists()
        assert ArtifactStore(output_dir).validate_sidecar(result.sidecars[0]).exists()
        assert result.payload["metrics"]["coverage"]["values"] == [1.0, 1.0, 0.5]
        assert result.payload["metrics"]["mean_contact_fraction"]["n"] == 3

    def test_aggregate_passes_disabled_residence_time_setting(self, tmp_path):
        import numpy as np

        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.mda import ArtifactStore

        analysis = ContactsAnalysis()
        settings = ContactsSettings(compute_residence_times=False)
        analysis_dir = tmp_path / "contacts"
        sim_config = _make_hashable_sim_config(tmp_path)
        cond = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=sim_config,
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1, 2),
            output_dir=analysis_dir / "aggregated",
            equilibration="10ns",
            settings=settings,
        )
        artifacts = [
            _write_contacts_replicate_artifact(
                analysis_dir,
                settings,
                replicate=replicate,
                event_durations_ns=[(0, "PEG", 5.0)],
            )
            for replicate in (1, 2)
        ]

        result = analysis.aggregate(ctx, artifacts)

        assert result.payload["residence_time_by_polymer_type"] == {}
        assert "residence_time_by_polymer_type" not in result.payload["residue_stats"][0]
        assert result.metadata["compute_residence_times"] is False
        profile_path = ArtifactStore(ctx.output_dir).validate_sidecar(result.sidecars[0])

        with np.load(profile_path) as profile:
            assert "residence_time_mean_ns" not in profile.files
            assert "residence_time_sem_ns" not in profile.files
            assert "residence_time_event_counts" not in profile.files

    def test_aggregate_rejects_stale_replicate_equilibration(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.mda.aggregation import MDAAggregationError

        settings = ContactsSettings()
        analysis_dir = tmp_path / "contacts"
        artifact = _write_contacts_replicate_artifact(
            analysis_dir,
            settings,
            replicate=1,
            equilibration="0ns",
        )
        condition = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1,),
            sim_config=_make_hashable_sim_config(tmp_path),
        )
        ctx = AggregateContext(
            condition=condition,
            replicates=(1,),
            output_dir=analysis_dir / "aggregated",
            equilibration="10ns",
            settings=settings,
        )

        with pytest.raises(MDAAggregationError, match="equilibration mismatch"):
            ContactsAnalysis().aggregate(ctx, [artifact])

    def test_aggregate_rejects_legacy_contact_results(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        cond = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=_make_hashable_sim_config(tmp_path),
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1, 2),
            output_dir=tmp_path / "aggregated",
            equilibration="10ns",
            settings=settings,
        )

        with pytest.raises(ValueError, match="ReplicateArtifact inputs"):
            analysis.aggregate(ctx, [_make_mock_contact_result(1), _make_mock_contact_result(2)])

    def test_aggregate_rejects_sidecar_hash_mismatch(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.mda.aggregation import MDAAggregationError

        settings = ContactsSettings()
        analysis_dir = tmp_path / "contacts"
        artifact = _write_contacts_replicate_artifact(analysis_dir, settings, replicate=1)
        (analysis_dir / "run_1" / artifact.payload["event_sidecar"]).write_bytes(b"stale")
        condition = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1,),
            sim_config=_make_hashable_sim_config(tmp_path),
        )
        ctx = AggregateContext(
            condition=condition,
            replicates=(1,),
            output_dir=analysis_dir / "aggregated",
            equilibration="10ns",
            settings=settings,
        )

        with pytest.raises(MDAAggregationError, match="invalid event sidecar|stale sidecar"):
            ContactsAnalysis().aggregate(ctx, [artifact])


class TestAggregatePolymerTypeCoverage:
    """Regression tests for polymer-type and RT artifact aggregation."""

    def test_includes_zero_contact_replicates_in_polymer_type_vectors(self, tmp_path):
        result = self._aggregate(tmp_path=tmp_path)
        residue = result.payload["residue_stats"][0]

        assert residue["by_polymer_type_per_replicate"]["PEG"] == [0.5, 0.3, 0.0]
        assert residue["by_polymer_type"]["PEG"]["mean"] == pytest.approx((0.5 + 0.3) / 3.0)

    def test_sparse_residence_time_vectors_record_replicate_identity(self, tmp_path):
        result = self._aggregate(tmp_path=tmp_path, events=[[(0, "PEG", 5.0)], [], []])
        residue = result.payload["residue_stats"][0]

        assert residue["residence_time_by_polymer_type"]["PEG"]["replicate_means_ns"] == [5.0]
        assert residue["residence_time_by_polymer_type"]["PEG"]["replicates_with_events"] == [1]

    def test_residence_time_summaries_are_computed_by_default(self, tmp_path):
        result = self._aggregate(
            tmp_path=tmp_path, events=[[(0, "PEG", 5.0)], [(0, "PEG", 3.0)], []]
        )

        summary = result.payload["residence_time_by_polymer_type"]["PEG"]
        assert summary["mean_ns"] == pytest.approx(4.0)
        assert summary["replicates_with_events"] == [1, 2]
        assert summary["n_events"] == 2

    def test_residence_time_summaries_empty_when_disabled(self, tmp_path):
        result = self._aggregate(
            tmp_path=tmp_path,
            settings_kwargs={"compute_residence_times": False},
            events=[[(0, "PEG", 5.0)], [(0, "PEG", 3.0)], []],
        )

        assert result.payload["residence_time_by_polymer_type"] == {}

    @staticmethod
    def _aggregate(tmp_path, settings_kwargs=None, events=None):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        settings = ContactsSettings(**(settings_kwargs or {}))
        analysis_dir = tmp_path / "contacts"
        event_sets = events or [[(0, "PEG", 5.0)], [(0, "PEG", 3.0)], []]
        fractions = [[0.5, 0.0], [0.3, 0.0], [0.0, 0.0]]
        artifacts = [
            _write_contacts_replicate_artifact(
                analysis_dir,
                settings,
                replicate=replicate,
                contact_fractions=fractions[index],
                event_durations_ns=event_sets[index],
            )
            for index, replicate in enumerate((1, 2, 3))
        ]
        condition = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2, 3),
            sim_config=_make_hashable_sim_config(tmp_path),
        )
        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2, 3),
            output_dir=analysis_dir / "aggregated",
            equilibration="10ns",
            settings=settings,
        )
        return ContactsAnalysis().aggregate(ctx, artifacts)


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
                [cond], settings=ContactsSettings(polymer_selection="chainid C and resname PEG")
            )

        assert result == []
        mock_universe.select_atoms.assert_called_with("chainid C and resname PEG")

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
        settings = ContactsSettings()
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
            settings=settings,
            aggregated_results={
                label: _make_condition_artifact(label, settings, replicates=(1, 2, 3))
                for label in labels
            },
            recompute=False,
        )

    def test_compare_returns_result(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        ctx = self._make_ctx(tmp_path, n_conditions=3, control="Control")

        result = analysis.compare(ctx)

        assert result is not None
        assert result.name == "test_comparison"
        assert len(result.conditions) == 3
        assert result.control_label == "Control"
        assert len(result.ranking_by_coverage) == 3
        assert len(result.ranking_by_contact_fraction) == 3
        assert "binding_preference" not in result.model_dump()
        assert "binding_preference" not in result.model_dump_json()

    def test_compare_pairwise_with_control(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        ctx = self._make_ctx(tmp_path, n_conditions=3, control="Control")

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

        result = analysis.compare(ctx)

        # Without control: all pairs = C(3,2) = 3
        assert len(result.pairwise_comparisons) == 3

    def test_compare_falls_back_when_control_aggregate_missing(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        ctx = self._make_ctx(tmp_path, n_conditions=3, control="Control")
        ctx.aggregated_results.pop("Control")

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

        result = analysis.compare(ctx)

        # ANOVA with 3+ conditions: 2 summaries (coverage + contact_fraction)
        assert len(result.anova) == 2
        metrics = {a.metric for a in result.anova}
        assert metrics == {"coverage", "mean_contact_fraction"}

    def test_compare_no_anova_with_two_conditions(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        ctx = self._make_ctx(tmp_path, n_conditions=2, control="Control")

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
            aggregated_results={
                "Only": _make_condition_artifact("Only", ContactsSettings(), replicates=(1, 2))
            },
            recompute=False,
        )

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
        aggregate = _make_condition_artifact("Finalized", settings, replicates=(1, 2))
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

        with patch.object(analysis, "_load_validated_aggregated_result") as load_mock:
            result = analysis.compare(ctx)

        assert result is not None
        assert result.conditions[0].label == "Finalized"
        load_mock.assert_not_called()

    def test_compare_rejects_stale_condition_equilibration(self, tmp_path):
        """Comparison should reject condition artifacts from a different window."""
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.mda import ArtifactStore

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        condition = Condition(
            label="Stale",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        analysis_dir = tmp_path / "Stale" / "contacts"
        agg_dir = analysis_dir / "aggregated"
        agg_dir.mkdir(parents=True)
        ArtifactStore(agg_dir).write_condition_result(
            _make_condition_artifact("Stale", settings, replicates=(1, 2), equilibration="0ns")
        )
        ctx = ComparisonContext(
            name="test",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"Stale": analysis_dir},
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=settings,
            recompute=False,
        )

        with pytest.raises(ValueError, match="equilibration mismatch"):
            analysis.compare(ctx)

    @pytest.mark.parametrize(
        ("case", "match"),
        [
            ("length", "length mismatch"),
            ("nonfinite", "non-finite values"),
            ("mean", "mean mismatch"),
            ("n_nonintegral", "invalid n"),
            ("n_bool", "invalid n"),
            ("replicate_keys", "replicate_metrics keys"),
        ],
    )
    def test_compare_rejects_invalid_condition_metric_integrity(self, tmp_path, case, match):
        """Comparison should reject stale or corrupted metric summaries."""
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        aggregate = _make_condition_artifact("Corrupt", settings, replicates=(1, 2))
        if case == "length":
            aggregate.payload["metrics"]["coverage"]["values"].append(0.5)
        elif case == "nonfinite":
            aggregate.payload["metrics"]["coverage"]["values"][0] = float("nan")
        elif case == "mean":
            aggregate.payload["metrics"]["mean_contact_fraction"]["mean"] = 42.0
        elif case == "n_nonintegral":
            aggregate.payload["metrics"]["coverage"]["n"] = 2.9
        elif case == "n_bool":
            aggregate.payload["metrics"]["coverage"]["n"] = True
        elif case == "replicate_keys":
            aggregate.payload["replicate_metrics"].pop("2")
        condition = Condition(
            label="Corrupt",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        ctx = ComparisonContext(
            name="test",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"Corrupt": tmp_path / "Corrupt" / "contacts"},
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=settings,
            aggregated_results={"Corrupt": aggregate},
            recompute=False,
        )

        with pytest.raises(ValueError, match=match):
            analysis.compare(ctx)

    def test_compare_rejects_residence_time_setting_mismatch(self, tmp_path):
        """Comparison should reject non-RT aggregates when RT is requested."""
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        current_settings = ContactsSettings(compute_residence_times=True)
        aggregate = _make_condition_artifact(
            "NoRT",
            current_settings,
            replicates=(1, 2),
            compute_residence_times=False,
        )
        condition = Condition(
            label="NoRT",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        ctx = ComparisonContext(
            name="test",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"NoRT": tmp_path / "NoRT" / "contacts"},
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=current_settings,
            aggregated_results={"NoRT": aggregate},
            recompute=False,
        )

        with pytest.raises(ValueError, match="compute_residence_times mismatch"):
            analysis.compare(ctx)

    def test_compare_rejects_missing_residence_time_setting_metadata(self, tmp_path):
        """Comparison should reject aggregates missing RT identity metadata."""
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        aggregate = _make_condition_artifact(
            "MissingRTMeta",
            settings,
            replicates=(1, 2),
            include_residence_metadata=False,
        )
        condition = Condition(
            label="MissingRTMeta",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        ctx = ComparisonContext(
            name="test",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"MissingRTMeta": tmp_path / "MissingRTMeta" / "contacts"},
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=settings,
            aggregated_results={"MissingRTMeta": aggregate},
            recompute=False,
        )

        with pytest.raises(ValueError, match="lacks compute_residence_times metadata"):
            analysis.compare(ctx)

    def test_compare_rejects_ambiguous_residence_time_setting_metadata(self, tmp_path):
        """Comparison should reject stringified RT identity metadata."""
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings(compute_residence_times=True)
        aggregate = _make_condition_artifact("StringRTMeta", settings, replicates=(1, 2))
        aggregate.metadata["compute_residence_times"] = "False"
        condition = Condition(
            label="StringRTMeta",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        ctx = ComparisonContext(
            name="test",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"StringRTMeta": tmp_path / "StringRTMeta" / "contacts"},
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=settings,
            aggregated_results={"StringRTMeta": aggregate},
            recompute=False,
        )

        with pytest.raises(ValueError, match="must be a boolean"):
            analysis.compare(ctx)

    def test_compare_rejects_legacy_aggregate_when_recompute_true(self, tmp_path):
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

        with pytest.raises(ValueError, match="canonical ConditionArtifact"):
            analysis.compare(ctx)

    def test_compare_ignores_root_level_legacy_aggregate(self, tmp_path):
        """Comparison loading should not scan root-level legacy contacts aggregates."""
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        condition = Condition(
            label="LegacyRoot",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        analysis_dir = tmp_path / "LegacyRoot" / "contacts"
        (analysis_dir / "aggregated").mkdir(parents=True)
        legacy_path = analysis_dir / "contacts_aggregated_eq10ns_cut4.5_reps1-2.json"
        _make_valid_agg_result(settings).save(legacy_path)
        ctx = ComparisonContext(
            name="test",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"LegacyRoot": analysis_dir},
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=settings,
            recompute=False,
        )

        result = analysis.compare(ctx)

        assert result is None

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
            aggregated_results={
                label: _make_condition_artifact(label, ContactsSettings(), replicates=(1, 2, 3))
                for label in ["A", "B"]
            },
            recompute=False,
        )

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
        from polyzymd.analyses.mda import ArtifactStore

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
        ArtifactStore(agg_dir).write_condition_result(
            _make_condition_artifact("Subset", settings, replicates=(1, 2))
        )
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

        result = analysis.compare(ctx)

        assert result is not None
        assert result.conditions[0].n_replicates == 2

    def test_compare_accepts_singleton_replicate_subset(self, tmp_path):
        """Finalize should accept one successful replicate for smoke-test comparisons."""
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.mda import ArtifactStore

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
        ArtifactStore(agg_dir).write_condition_result(
            _make_condition_artifact("Subset", settings, replicates=(1,))
        )
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

        result = analysis.compare(ctx)

        assert result is not None
        assert result.conditions[0].n_replicates == 1
        assert result.pairwise_comparisons == []


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
            '{"metadata": {"contacts_settings_fingerprint": "deadbeef", "equilibration": "10ns"}}'
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
            '{"metadata": {"contacts_settings_fingerprint": "deadbeef", "equilibration": "10ns"}}'
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
            '{"metadata": {"contacts_settings_fingerprint": "deadbeef", "equilibration": "10ns"}}'
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
            '{"metadata": {"contacts_settings_fingerprint": "feedface", "equilibration": "10ns"}}'
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

    def test_aggregated_sidecar_names_noncontiguous_replicates_uniquely(self, tmp_path):
        """Non-contiguous replicate IDs should not collide with contiguous ranges."""
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        settings = ContactsSettings()

        sparse = ContactsAnalysis._aggregated_sidecar_path(tmp_path, settings, "10ns", (1, 3))
        dense = ContactsAnalysis._aggregated_sidecar_path(tmp_path, settings, "10ns", (1, 2, 3))

        assert sparse.name != dense.name
        assert "_reps1_3.json" in sparse.name
        assert "_reps1-3.json" in dense.name
