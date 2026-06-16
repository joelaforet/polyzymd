"""Tests for the contacts analysis plugin.

Covers discovery, settings, compute-stage dispatch, aggregate, compare (full override),
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

        with pytest.raises(ValueError, match="not shipped as an active PolyzyMD") as excinfo:
            ContactsSettings.model_validate({key: value})
        assert "feature/mda-analysis-migration" not in str(excinfo.value)

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

        with pytest.raises(ValueError, match="not shipped as active PolyzyMD") as excinfo:
            ContactsPlotSettings.model_validate({key: value})
        assert "feature/mda-analysis-migration" not in str(excinfo.value)

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

    def test_detection_fingerprint_changes_with_cutoff(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_detection_fingerprint

        low = ContactsSettings(cutoff=4.0)
        high = ContactsSettings(cutoff=4.5)

        assert contacts_detection_fingerprint(low) != contacts_detection_fingerprint(high)

    def test_detection_fingerprint_excludes_comparison_partition_fields(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_detection_fingerprint

        base = ContactsSettings(cutoff=4.5)
        downstream_only = ContactsSettings(
            cutoff=4.5,
            protein_groups={"active": [1, 2]},
            protein_partitions={"site": ["active"]},
            top_residues=25,
        )

        assert contacts_detection_fingerprint(base) == contacts_detection_fingerprint(
            downstream_only
        )

    def test_detection_fingerprint_changes_with_contact_fields(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_detection_fingerprint

        unfiltered = ContactsSettings(polymer_selection="chainid C")
        filtered = ContactsSettings(polymer_selection="chainid C", polymer_types=["PEG"])

        assert contacts_detection_fingerprint(unfiltered) != contacts_detection_fingerprint(
            filtered
        )


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
        assert not hasattr(Analysis, "build_runner")
        assert not hasattr(Analysis, "summarize_replicate")
        assert "build_runner" not in ContactsAnalysis.__dict__
        assert "summarize_replicate" not in ContactsAnalysis.__dict__

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


class TestContactsSparseEventUtilities:
    """Contacts-local helpers preserve existing scientific semantics."""

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


def _write_contacts_replicate_artifact(
    analysis_dir: Path,
    settings,
    *,
    condition_label: str = "test",
    replicate: int = 1,
    equilibration: str = "10ns",
    contact_fractions: list[float] | None = None,
    event_durations_ns: list[tuple[int, str, float]] | None = None,
    frame_selection_overrides: dict[str, object] | None = None,
    sidecar_frame_indices: object | None = None,
    sidecar_time_ps: list[float] | None = None,
):
    """Write a synthetic contacts replicate artifact and event sidecar.

    Parameters
    ----------
    analysis_dir : Path
        Contacts analysis directory that will receive ``run_<replicate>`` artifacts.
    settings : Any
        Contacts settings used to build detection identity metadata.
    condition_label : str, optional
        Condition label for the artifact, by default ``"test"``.
    replicate : int, optional
        Replicate identifier, by default ``1``.
    equilibration : str, optional
        Equilibration window string, by default ``"10ns"``.
    contact_fractions : list[float] | None, optional
        Synthetic per-residue contact fractions, by default ``None``.
    event_durations_ns : list[tuple[int, str, float]] | None, optional
        Synthetic residence-time events, by default ``None``.
    frame_selection_overrides : dict[str, object] | None, optional
        Frame-selection provenance fields to override, by default ``None``.
    sidecar_frame_indices : object, optional
        Explicit sidecar frame-index array to write, by default ``None``.
    sidecar_time_ps : list of float or None, optional
        Explicit sidecar time axis to write, by default ``None``.

    Returns
    -------
    ReplicateArtifact
        Written contacts replicate artifact.
    """

    import numpy as np

    from polyzymd.analyses.contacts._identity import (
        CONTACT_SEMANTICS_VERSION,
        CONTACTS_PBC_POLICY,
        contacts_detection_fingerprint,
        contacts_detection_identity_payload,
    )
    from polyzymd.analyses.mda import ArtifactStore, ReplicateArtifact
    from polyzymd.analyses.shared.loader import convert_time, parse_time_string

    fractions = contact_fractions or [0.5, 0.0]
    events = event_durations_ns or []
    run_dir = analysis_dir / f"run_{replicate}"
    store = ArtifactStore(run_dir)
    equilibration_value, equilibration_unit = parse_time_string(equilibration)
    equilibration_ps = float(convert_time(equilibration_value, equilibration_unit, "ps"))
    frame_selection = {
        "start": 0,
        "stop": 10,
        "step": 1,
        "frames": None,
        "equilibration": equilibration,
        "equilibration_start": 0,
        "equilibration_ps": equilibration_ps,
        "timestep_ps": 1000.0,
        "first_frame_time_ps": equilibration_ps,
        "selected_start_time_ps": equilibration_ps,
        "equilibration_time_reference": "trajectory_timestamp",
        "n_frames_total": 10,
        "n_frames_selected": 10,
        "warning_message": None,
    }
    if frame_selection_overrides is not None:
        frame_selection.update(frame_selection_overrides)
    n_frames_selected = int(frame_selection["n_frames_selected"])
    start = int(frame_selection["start"])
    stop = int(frame_selection["stop"])
    step = int(frame_selection.get("step") or 1)
    timestep_ps = float(frame_selection["timestep_ps"])
    first_frame_time_ps = frame_selection.get("first_frame_time_ps")
    time_reference = str(frame_selection.get("equilibration_time_reference"))
    time_origin_ps = (
        float(first_frame_time_ps)
        if time_reference == "trajectory_timestamp" and first_frame_time_ps is not None
        else 0.0
    )
    frame_indices = np.arange(start, stop, step, dtype=np.int64)
    if frame_indices.size != n_frames_selected:
        frame_indices = frame_indices[:n_frames_selected]
    sidecar_frame_indices_array = (
        np.asarray(sidecar_frame_indices) if sidecar_frame_indices is not None else frame_indices
    )
    time_ps = time_origin_ps + frame_indices.astype(np.float64) * timestep_ps
    if sidecar_time_ps is not None:
        time_ps = np.asarray(sidecar_time_ps, dtype=np.float64)
    protein_resids = np.arange(1, len(fractions) + 1, dtype=np.int64)
    protein_resnames = np.asarray(["ALA", "ASP"][: len(fractions)], dtype="U16")
    protein_groups = np.asarray(["nonpolar", "charged"][: len(fractions)], dtype="U32")
    polymer_resnames = np.asarray(["PEG"], dtype="U16")
    start_indices = np.arange(len(events), dtype=np.int64)
    event_start_frames = (
        frame_indices[start_indices] if frame_indices.size >= start_indices.size else start_indices
    )
    durations_ns = np.asarray([event[2] for event in events], dtype=np.float64)
    protein_indices = np.asarray([event[0] for event in events], dtype=np.int64)
    sidecar = store.write_npz_sidecar(
        "sidecars/contact_events.npz",
        event_start_sample_index=start_indices,
        event_start_frame=event_start_frames,
        event_duration_samples=np.ones(len(events), dtype=np.int64),
        event_duration_ps=durations_ns * 1000.0,
        event_duration_ns=durations_ns,
        protein_residue_index=protein_indices,
        polymer_residue_index=np.zeros(len(events), dtype=np.int64),
        polymer_chain_index=np.zeros(len(events), dtype=np.int64),
        frame_indices=sidecar_frame_indices_array,
        time_ps=time_ps,
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
            "n_frames_used": n_frames_selected,
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
            "frame_selection": frame_selection,
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

    def test_aggregate_artifacts_returns_condition_artifact_without_persistence(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.mda import ArtifactStore, ConditionArtifact

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
        assert not (output_dir / "result.json").exists()
        assert ArtifactStore(output_dir).validate_sidecar(result.sidecars[0]).exists()
        assert result.payload["metrics"]["coverage"]["values"] == [1.0, 1.0, 0.5]
        assert result.payload["metrics"]["mean_contact_fraction"]["n"] == 3

    def test_aggregate_accepts_timestamp_provenance_differences(self, tmp_path):
        """Aggregation should allow per-replicate timestamp provenance values."""

        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.mda import ConditionArtifact

        settings = ContactsSettings()
        analysis_dir = tmp_path / "contacts"
        artifacts = [
            _write_contacts_replicate_artifact(
                analysis_dir,
                settings,
                replicate=1,
                equilibration="200ns",
                frame_selection_overrides={
                    "start": 4,
                    "stop": 14,
                    "equilibration_start": 4,
                    "timestep_ps": 400.0,
                    "first_frame_time_ps": 198_400.0,
                    "selected_start_time_ps": 200_000.0,
                    "n_frames_total": 20,
                    "warning_message": "replicate 1 warning",
                },
            ),
            _write_contacts_replicate_artifact(
                analysis_dir,
                settings,
                replicate=2,
                equilibration="200ns",
                frame_selection_overrides={
                    "start": 0,
                    "stop": 10,
                    "equilibration_start": 0,
                    "timestep_ps": 400.0,
                    "first_frame_time_ps": 202_400.0,
                    "selected_start_time_ps": 202_400.0,
                    "warning_message": "replicate 2 warning",
                },
            ),
        ]
        cond = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=_make_hashable_sim_config(tmp_path),
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1, 2),
            output_dir=analysis_dir / "aggregated",
            equilibration="200ns",
            settings=settings,
        )

        result = ContactsAnalysis().aggregate(ctx, artifacts)

        assert isinstance(result, ConditionArtifact)
        assert result.payload["metrics"]["mean_contact_fraction"]["n"] == 2

    def test_aggregate_accepts_mixed_timestamp_derived_starts(self, tmp_path):
        """Aggregation should allow timestamp-derived start-frame differences."""

        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        settings = ContactsSettings()
        analysis_dir = tmp_path / "contacts"
        artifacts = [
            _write_contacts_replicate_artifact(
                analysis_dir,
                settings,
                replicate=1,
                frame_selection_overrides={
                    "start": 4,
                    "stop": 14,
                    "equilibration_start": 4,
                    "first_frame_time_ps": 6_000.0,
                    "selected_start_time_ps": 10_000.0,
                    "n_frames_total": 20,
                },
            ),
            _write_contacts_replicate_artifact(
                analysis_dir,
                settings,
                replicate=2,
                frame_selection_overrides={
                    "start": 0,
                    "stop": 10,
                    "equilibration_start": 0,
                    "first_frame_time_ps": 10_000.0,
                    "selected_start_time_ps": 10_000.0,
                    "n_frames_total": 10,
                },
            ),
        ]
        cond = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=_make_hashable_sim_config(tmp_path),
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1, 2),
            output_dir=analysis_dir / "aggregated",
            equilibration="10ns",
            settings=settings,
        )

        result = ContactsAnalysis().aggregate(ctx, artifacts)

        assert result.payload["total_frames_per_replicate"] == [10, 10]

    def test_aggregate_accepts_rounded_timestamp_sidecar_axis(self, tmp_path):
        """Aggregation should tolerate harmless timestamp rounding drift."""

        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        settings = ContactsSettings()
        analysis_dir = tmp_path / "contacts"
        artifact = _write_contacts_replicate_artifact(
            analysis_dir,
            settings,
            equilibration="0ns",
            frame_selection_overrides={
                "start": 0,
                "stop": 3,
                "equilibration_start": 0,
                "equilibration_ps": 0.0,
                "timestep_ps": 33.959999,
                "first_frame_time_ps": 0.0,
                "selected_start_time_ps": 0.0,
                "n_frames_total": 3,
                "n_frames_selected": 3,
            },
            sidecar_time_ps=[0.0, 33.96, 67.92],
        )
        cond = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1,),
            sim_config=_make_hashable_sim_config(tmp_path),
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1,),
            output_dir=analysis_dir / "aggregated",
            equilibration="0ns",
            settings=settings,
        )

        result = ContactsAnalysis().aggregate(ctx, [artifact])

        assert result.payload["total_frames_per_replicate"] == [3]

    def test_aggregate_rejects_shifted_timestamp_sidecar_axis(self, tmp_path):
        """Aggregation should still reject clearly shifted timestamp sidecars."""

        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.mda.aggregation import MDAAggregationError

        settings = ContactsSettings()
        analysis_dir = tmp_path / "contacts"
        artifact = _write_contacts_replicate_artifact(
            analysis_dir,
            settings,
            equilibration="0ns",
            frame_selection_overrides={
                "start": 0,
                "stop": 3,
                "equilibration_start": 0,
                "equilibration_ps": 0.0,
                "timestep_ps": 33.959999,
                "first_frame_time_ps": 0.0,
                "selected_start_time_ps": 0.0,
                "n_frames_total": 3,
                "n_frames_selected": 3,
            },
            sidecar_time_ps=[1.0, 34.96, 68.92],
        )
        cond = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1,),
            sim_config=_make_hashable_sim_config(tmp_path),
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1,),
            output_dir=analysis_dir / "aggregated",
            equilibration="0ns",
            settings=settings,
        )

        with pytest.raises(MDAAggregationError, match="sidecar first time mismatch"):
            ContactsAnalysis().aggregate(ctx, [artifact])

    @pytest.mark.parametrize(
        "sidecar_frame_indices",
        [[0.1, 1.1], [0.0, float("nan")], [0.0, float("inf")]],
    )
    def test_aggregate_rejects_non_integer_timestamp_sidecar_frame_indices(
        self, tmp_path, sidecar_frame_indices
    ):
        """Aggregation should reject fractional or non-finite sidecar frame indices."""

        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.mda.aggregation import MDAAggregationError

        settings = ContactsSettings()
        analysis_dir = tmp_path / "contacts"
        artifact = _write_contacts_replicate_artifact(
            analysis_dir,
            settings,
            equilibration="0ns",
            frame_selection_overrides={
                "start": 0,
                "stop": 2,
                "equilibration_start": 0,
                "equilibration_ps": 0.0,
                "timestep_ps": 1000.0,
                "first_frame_time_ps": 0.0,
                "selected_start_time_ps": 0.0,
                "n_frames_total": 2,
                "n_frames_selected": 2,
            },
            sidecar_frame_indices=sidecar_frame_indices,
        )
        cond = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1,),
            sim_config=_make_hashable_sim_config(tmp_path),
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1,),
            output_dir=analysis_dir / "aggregated",
            equilibration="0ns",
            settings=settings,
        )

        with pytest.raises(MDAAggregationError, match="sidecar frame_indices"):
            ContactsAnalysis().aggregate(ctx, [artifact])

    def test_aggregate_rejects_stale_timestamp_window_for_requested_equilibration(self, tmp_path):
        """Aggregation should validate timestamp-derived starts per artifact."""

        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.mda.aggregation import MDAAggregationError

        settings = ContactsSettings()
        analysis_dir = tmp_path / "contacts"
        artifact = _write_contacts_replicate_artifact(
            analysis_dir,
            settings,
            replicate=1,
            equilibration="200ns",
            frame_selection_overrides={
                "start": 0,
                "stop": 10,
                "equilibration_start": 0,
                "timestep_ps": 400.0,
                "first_frame_time_ps": 198_400.0,
                "selected_start_time_ps": 198_400.0,
                "n_frames_total": 20,
            },
        )
        cond = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1,),
            sim_config=_make_hashable_sim_config(tmp_path),
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1,),
            output_dir=analysis_dir / "aggregated",
            equilibration="200ns",
            settings=settings,
        )

        with pytest.raises(MDAAggregationError, match="frame-selection start mismatch"):
            ContactsAnalysis().aggregate(ctx, [artifact])

    @pytest.mark.parametrize(
        ("candidate_overrides", "message"),
        [
            (
                {"start": 11, "stop": 21, "selected_start_time_ps": 11_000.0},
                "frame-selection start mismatch",
            ),
            ({"equilibration_start": 11}, "frame-selection equilibration_start mismatch"),
        ],
    )
    def test_aggregate_rejects_noncanonical_frame_selection_window_mismatch(
        self,
        tmp_path,
        candidate_overrides,
        message,
    ):
        """Aggregation should reject mismatched loaded-frame-relative windows."""

        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.mda.aggregation import MDAAggregationError

        settings = ContactsSettings()
        analysis_dir = tmp_path / "contacts"
        loaded_frame_zero_window = {
            "start": 10,
            "stop": 20,
            "equilibration_start": 10,
            "first_frame_time_ps": None,
            "selected_start_time_ps": 10_000.0,
            "equilibration_time_reference": "loaded_frame_zero",
            "n_frames_total": 20,
            "n_frames_selected": 10,
        }
        artifacts = [
            _write_contacts_replicate_artifact(
                analysis_dir,
                settings,
                replicate=1,
                frame_selection_overrides=loaded_frame_zero_window,
            ),
            _write_contacts_replicate_artifact(
                analysis_dir,
                settings,
                replicate=2,
                frame_selection_overrides=loaded_frame_zero_window | candidate_overrides,
            ),
        ]
        cond = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=_make_hashable_sim_config(tmp_path),
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1, 2),
            output_dir=analysis_dir / "aggregated",
            equilibration="10ns",
            settings=settings,
        )

        with pytest.raises(MDAAggregationError, match=message):
            ContactsAnalysis().aggregate(ctx, artifacts)

    def test_aggregate_rejects_single_stale_noncanonical_loaded_frame_zero_window(self, tmp_path):
        """Aggregation should validate the first artifact against the requested window."""

        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.mda.aggregation import MDAAggregationError

        settings = ContactsSettings()
        analysis_dir = tmp_path / "contacts"
        stale_loaded_frame_zero_window = {
            "start": 0,
            "equilibration_start": 0,
            "first_frame_time_ps": None,
            "selected_start_time_ps": 0.0,
            "equilibration_time_reference": "loaded_frame_zero",
        }
        artifact = _write_contacts_replicate_artifact(
            analysis_dir,
            settings,
            replicate=1,
            frame_selection_overrides=stale_loaded_frame_zero_window,
        )
        cond = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1,),
            sim_config=_make_hashable_sim_config(tmp_path),
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1,),
            output_dir=analysis_dir / "aggregated",
            equilibration="10ns",
            settings=settings,
        )

        with pytest.raises(
            MDAAggregationError,
            match="frame-selection equilibration_start mismatch for replicate 1",
        ):
            ContactsAnalysis().aggregate(ctx, [artifact])

    def test_aggregate_rejects_loaded_frame_zero_stale_equilibration_ps(self, tmp_path):
        """Aggregation should validate stored offsets against requested equilibration."""

        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.mda.aggregation import MDAAggregationError

        settings = ContactsSettings()
        analysis_dir = tmp_path / "contacts"
        stale_loaded_frame_zero_window = {
            "start": 0,
            "equilibration_start": 0,
            "equilibration_ps": 0.0,
            "first_frame_time_ps": None,
            "selected_start_time_ps": 0.0,
            "equilibration_time_reference": "loaded_frame_zero",
        }
        artifact = _write_contacts_replicate_artifact(
            analysis_dir,
            settings,
            replicate=1,
            equilibration="10ns",
            frame_selection_overrides=stale_loaded_frame_zero_window,
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

        with pytest.raises(
            MDAAggregationError,
            match="frame-selection equilibration_ps mismatch for replicate 1",
        ):
            ContactsAnalysis().aggregate(ctx, [artifact])

    def test_aggregate_rejects_matching_stale_noncanonical_loaded_frame_zero_windows(
        self, tmp_path
    ):
        """Aggregation should not trust the first stale artifact as the baseline."""

        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.mda.aggregation import MDAAggregationError

        settings = ContactsSettings()
        analysis_dir = tmp_path / "contacts"
        stale_loaded_frame_zero_window = {
            "start": 0,
            "equilibration_start": 0,
            "first_frame_time_ps": None,
            "selected_start_time_ps": 0.0,
            "equilibration_time_reference": "loaded_frame_zero",
        }
        artifacts = [
            _write_contacts_replicate_artifact(
                analysis_dir,
                settings,
                replicate=replicate,
                frame_selection_overrides=stale_loaded_frame_zero_window,
            )
            for replicate in (1, 2)
        ]
        cond = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=_make_hashable_sim_config(tmp_path),
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1, 2),
            output_dir=analysis_dir / "aggregated",
            equilibration="10ns",
            settings=settings,
        )

        with pytest.raises(
            MDAAggregationError,
            match="frame-selection equilibration_start mismatch for replicate 1",
        ):
            ContactsAnalysis().aggregate(ctx, artifacts)

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

    def test_aggregate_rejects_noncanonical_contact_results(self, tmp_path):
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
            analysis.aggregate(ctx, [object(), object()])

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
    """Test contacts replicate metric helpers."""

    def test_coverage_per_replicate(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._comparison import compute_coverage_per_replicate

        artifact = _make_condition_artifact("A", ContactsSettings(), replicates=(1, 2, 3))

        assert compute_coverage_per_replicate(artifact) == [1.0, 1.0, 1.0]

    def test_coverage_with_zero_fractions(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._comparison import compute_coverage_per_replicate

        artifact = _make_condition_artifact("A", ContactsSettings(), replicates=(1, 2))
        artifact.payload["metrics"]["coverage"]["values"] = [2 / 3, 1 / 3]

        coverages = compute_coverage_per_replicate(artifact)
        assert coverages[0] == pytest.approx(2 / 3)
        assert coverages[1] == pytest.approx(1 / 3)

    def test_contact_fraction_per_replicate(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._comparison import (
            compute_contact_fraction_per_replicate,
        )

        artifact = _make_condition_artifact("A", ContactsSettings(), replicates=(1, 2))
        artifact.payload["metrics"]["mean_contact_fraction"]["values"] = [0.3, 0.7]

        assert compute_contact_fraction_per_replicate(artifact) == [0.3, 0.7]


# ---------------------------------------------------------------------------
# Residue set validation
# ---------------------------------------------------------------------------


class TestResidueSetValidation:
    """Test contacts residue-set validation."""

    def test_matching_residue_sets(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._comparison import validate_residue_sets

        cond_a = MagicMock()
        cond_a.label = "A"
        data_a = {"agg_result": _make_condition_artifact("A", ContactsSettings(), n_residues=3)}

        cond_b = MagicMock()
        cond_b.label = "B"
        data_b = {"agg_result": _make_condition_artifact("B", ContactsSettings(), n_residues=3)}

        # Should not raise
        validate_residue_sets([(cond_a, data_a), (cond_b, data_b)])

    def test_mismatched_residue_sets(self):
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._comparison import validate_residue_sets

        cond_a = MagicMock()
        cond_a.label = "A"
        data_a = {"agg_result": _make_condition_artifact("A", ContactsSettings(), n_residues=3)}

        cond_b = MagicMock()
        cond_b.label = "B"
        data_b = {"agg_result": _make_condition_artifact("B", ContactsSettings(), n_residues=3)}
        data_b["agg_result"].payload["residue_stats"][2]["protein_resid"] = 4

        with pytest.raises(ValueError, match="Residue set mismatch"):
            validate_residue_sets([(cond_a, data_a), (cond_b, data_b)])


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
        mock_sim_config.engine = "openmm"
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
        mock_sim_config.engine = "openmm"
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

    def test_missing_config_engine_raises_value_error(self, tmp_path):
        """Missing engine configuration should not be swallowed by fail-open filtering."""
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()

        mock_sim_config = MagicMock()
        del mock_sim_config.engine
        mock_sim_config.output.projects_directory = tmp_path / "projects"
        mock_sim_config.get_working_directory.return_value = run_dir
        cond = Condition(
            label="MissingEngine",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim_config,
        )

        with pytest.raises(ValueError, match="non-empty string engine"):
            analysis.filter_conditions([cond])

    def test_invalid_config_engine_raises_value_error(self, tmp_path):
        """Unknown engine configuration should not be swallowed by fail-open filtering."""
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()

        mock_sim_config = MagicMock()
        mock_sim_config.engine = "namd"
        mock_sim_config.output.projects_directory = tmp_path / "projects"
        mock_sim_config.get_working_directory.return_value = run_dir
        cond = Condition(
            label="InvalidEngine",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim_config,
        )

        with pytest.raises(ValueError, match="Unknown engine"):
            analysis.filter_conditions([cond])

    def test_missing_topology_includes_condition(self, tmp_path):
        """Conditions are fail-open when no replicate topology can be inspected."""
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()

        mock_sim_config = MagicMock()
        mock_sim_config.engine = "openmm"
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

        result = analysis.compare(ctx)

        assert result is not None
        assert result.conditions[0].label == "Finalized"

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
        # Simulate framework-owned condition artifact persistence for compare
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

    def test_compare_rejects_noncanonical_aggregate_when_recompute_true(self, tmp_path):
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
        (agg_dir / "result.json").write_text(
            '{"analysis_type":"contacts_aggregated","n_replicates":2}', encoding="utf-8"
        )
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

    def test_compare_ignores_root_level_noncanonical_aggregate(self, tmp_path):
        """Comparison loading should not scan root-level non-canonical contacts aggregates."""
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        condition = Condition(
            label="Non-canonicalRoot",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        analysis_dir = tmp_path / "Non-canonicalRoot" / "contacts"
        (analysis_dir / "aggregated").mkdir(parents=True)
        noncanonical_path = analysis_dir / "contacts_aggregated_eq10ns_cut4.5_reps1-2.json"
        noncanonical_path.write_text(
            '{"analysis_type":"contacts_aggregated","n_replicates":2}', encoding="utf-8"
        )
        ctx = ComparisonContext(
            name="test",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"Non-canonicalRoot": analysis_dir},
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

    def test_compare_rejects_detection_fingerprint_mismatch(self, tmp_path):
        """Compare should not accept condition artifacts from different settings."""
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        current_settings = ContactsSettings(cutoff=4.5)
        stale_settings = ContactsSettings(cutoff=4.0)
        condition = Condition(
            label="A",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        ctx = ComparisonContext(
            name="test",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"A": tmp_path / "A" / "contacts"},
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=current_settings,
            aggregated_results={
                "A": _make_condition_artifact("A", stale_settings, replicates=(1, 2, 3))
            },
        )

        with pytest.raises(ValueError, match="detection fingerprint mismatch"):
            analysis.compare(ctx)

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
        # Simulate framework-owned condition artifact persistence for compare
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
        # Simulate framework-owned condition artifact persistence for compare
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
    """Test contacts pairwise comparison helper."""

    def test_pair_produces_two_metrics(self):
        from polyzymd.analyses.contacts._comparison import compare_contacts_pair

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

        comp = compare_contacts_pair("A", summary_a, data_a, "B", summary_b, data_b)

        assert comp.condition_a == "A"
        assert comp.condition_b == "B"
        assert len(comp.aggregate_comparisons) == 2

        metrics = {ac.metric for ac in comp.aggregate_comparisons}
        assert metrics == {"coverage", "mean_contact_fraction"}

    def test_pair_direction_labels(self):
        from polyzymd.analyses.contacts._comparison import compare_contacts_pair

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

        comp = compare_contacts_pair("Control", summary_a, data_a, "Treatment", summary_b, data_b)

        for ac in comp.aggregate_comparisons:
            assert ac.direction in ("increased", "decreased", "unchanged")


# ---------------------------------------------------------------------------
# ANOVA
# ---------------------------------------------------------------------------


class TestANOVA:
    """Test contacts ANOVA helper."""

    def test_anova_returns_two_summaries(self):
        from polyzymd.analyses.contacts._comparison import compute_contacts_anova

        condition_data = []
        for i in range(3):
            cond = MagicMock()
            cond.label = f"Cond{i}"
            data = {
                "coverage_per_replicate": [0.7 + i * 0.05 + j * 0.01 for j in range(3)],
                "contact_fraction_per_replicate": [0.2 + i * 0.03 + j * 0.005 for j in range(3)],
            }
            condition_data.append((cond, data))

        results = compute_contacts_anova(condition_data)
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
        from polyzymd.analyses.mda import ConditionArtifact

        assert ContactsAnalysis.AggregatedResultClass is ConditionArtifact

    def test_deserialize_loads_condition_artifact(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.mda import ArtifactStore, ConditionArtifact

        analysis = ContactsAnalysis()
        artifact = _make_condition_artifact("Canonical", ContactsSettings(), replicates=(1, 2))
        ArtifactStore(tmp_path).write_condition_result(artifact)

        result = analysis._deserialize_result(tmp_path / "result.json")

        assert isinstance(result, ConditionArtifact)
        assert result.condition_label == "Canonical"


# ---------------------------------------------------------------------------
# artifact-only contacts plotting
# ---------------------------------------------------------------------------


def _write_contacts_plot_artifact(
    analysis_dir: Path,
    settings,
    *,
    label: str = "A",
    replicates: tuple[int, ...] = (1, 2),
    compute_residence_times: bool | None = None,
) -> Path:
    """Write a synthetic condition artifact and profile sidecar for plot tests."""

    import numpy as np

    from polyzymd.analyses.contacts._identity import contacts_detection_fingerprint
    from polyzymd.analyses.mda import ArtifactStore, ConditionArtifact

    residence_enabled = (
        bool(settings.compute_residence_times)
        if compute_residence_times is None
        else bool(compute_residence_times)
    )
    aggregated_dir = analysis_dir / "aggregated"
    store = ArtifactStore(aggregated_dir)
    contact_by_replicate = np.asarray(
        [[0.40, 0.20, 0.10], [0.60, 0.30, 0.20]][: len(replicates)], dtype=np.float64
    )
    profile_arrays = {
        "replicates": np.asarray(replicates, dtype=np.int64),
        "protein_resids": np.asarray([1, 2, 3], dtype=np.int64),
        "protein_resnames": np.asarray(["ALA", "ASP", "SER"], dtype="U16"),
        "protein_groups": np.asarray(["nonpolar", "charged_negative", "polar"], dtype="U32"),
        "contact_fraction_by_replicate": contact_by_replicate,
        "contact_fraction_mean": np.mean(contact_by_replicate, axis=0),
        "contact_fraction_sem": np.std(contact_by_replicate, axis=0, ddof=1)
        / (len(replicates) ** 0.5),
        "polymer_types": np.asarray(["PEG"], dtype="U16"),
        "contact_fraction_by_polymer_type": contact_by_replicate.reshape(1, len(replicates), 3),
    }
    if residence_enabled:
        profile_arrays.update(
            residence_time_mean_ns=np.asarray([[4.0, 2.0, 1.0]], dtype=np.float64),
            residence_time_sem_ns=np.asarray([[0.5, 0.2, 0.1]], dtype=np.float64),
            residence_time_event_counts=np.asarray([[3, 2, 1]], dtype=np.int64),
        )
    sidecar = store.write_npz_sidecar(
        "sidecars/contact_profiles.npz",
        metadata={
            "kind": "contact_profiles",
            "layout": "condition_profile_table",
            "compute_residence_times": residence_enabled,
        },
        **profile_arrays,
    )
    coverage_values = [1.0 for _replicate in replicates]
    contact_values = [float(np.mean(row)) for row in contact_by_replicate]
    artifact = ConditionArtifact(
        analysis_name="contacts",
        condition_label=label,
        replicates=list(replicates),
        payload={
            "metrics": {
                "coverage": {
                    "name": "coverage",
                    "values": coverage_values,
                    "mean": float(np.mean(coverage_values)),
                    "sem": 0.0,
                    "std": 0.0,
                    "n": len(coverage_values),
                },
                "mean_contact_fraction": {
                    "name": "mean_contact_fraction",
                    "values": contact_values,
                    "mean": float(np.mean(contact_values)),
                    "sem": 0.0,
                    "std": 0.0,
                    "n": len(contact_values),
                },
            },
            "replicate_metrics": {
                str(replicate): {
                    "coverage": coverage_values[index],
                    "mean_contact_fraction": contact_values[index],
                }
                for index, replicate in enumerate(replicates)
            },
            "n_residues": 3,
            "polymer_types": ["PEG"],
            "profile_sidecar": sidecar.path,
        },
        sidecars=[sidecar],
        provenance={"profile_sidecar": sidecar.path},
        metadata={
            "contacts_detection_fingerprint": contacts_detection_fingerprint(settings),
            "compute_residence_times": residence_enabled,
            "equilibration": "10ns",
        },
    )
    # Simulate framework-owned condition artifact persistence for artifact-only plots
    store.write_condition_result(artifact)
    return aggregated_dir / "result.json"


def _contacts_plot_context(tmp_path, settings, plot_settings=None):
    """Build a contacts PlotContext for artifact-only plot tests."""

    from polyzymd.analyses.base import Condition, PlotContext
    from polyzymd.config.comparison import PlotSettings

    condition = Condition(
        label="A",
        config_path=Path("/tmp/a/config.yaml"),
        replicates=(1, 2),
        sim_config=object(),
    )
    return PlotContext(
        conditions=[condition],
        analysis_dirs={"A": tmp_path / "A" / "contacts"},
        results_dir=tmp_path / "results",
        output_dir=tmp_path / "plots",
        settings=settings,
        plot_settings=plot_settings or PlotSettings(),
        equilibration="10ns",
    )


class TestPlot:
    """Contacts plot tests using only condition artifacts and profile sidecars."""

    def test_all_plot_families_generate_files_from_artifacts(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        settings = ContactsSettings(
            protein_groups={"active_site": [1, 2], "surface": [3]},
            protein_partitions={"regions": ["active_site", "surface"]},
        )
        ctx = _contacts_plot_context(tmp_path, settings)
        _write_contacts_plot_artifact(ctx.analysis_dirs["A"], settings)

        paths = ContactsAnalysis().plot(ctx)

        assert len(paths) == 6
        assert {path.stem for path in paths} == {
            "contact_fraction_profile",
            "residence_time_profile",
            "cf_by_aa_class_bars",
            "cf_by_partition_regions_bars",
            "rt_by_aa_class_bars",
            "rt_by_partition_regions_bars",
        }
        assert all(path.exists() for path in paths)

    def test_plot_ignores_noncanonical_comparison_result_json(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        settings = ContactsSettings(
            protein_groups={"active_site": [1, 2], "surface": [3]},
            protein_partitions={"regions": ["active_site", "surface"]},
        )
        ctx = _contacts_plot_context(tmp_path, settings)
        _write_contacts_plot_artifact(ctx.analysis_dirs["A"], settings)
        ctx.results_dir.mkdir(parents=True)
        (ctx.results_dir / "result.json").write_text("not a comparison artifact", encoding="utf-8")

        with patch(
            "polyzymd.analyses.contacts._plotters.ArtifactStore.read_comparison_result",
            side_effect=AssertionError("comparison result should not be loaded"),
        ):
            paths = ContactsAnalysis().plot(ctx)

        assert paths
        assert all(path.exists() for path in paths)

    def test_plot_skips_residence_time_families_when_disabled(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        settings = ContactsSettings(
            compute_residence_times=False,
            protein_groups={"active_site": [1, 2], "surface": [3]},
            protein_partitions={"regions": ["active_site", "surface"]},
        )
        ctx = _contacts_plot_context(tmp_path, settings)
        _write_contacts_plot_artifact(ctx.analysis_dirs["A"], settings)

        paths = ContactsAnalysis().plot(ctx)

        assert {path.stem for path in paths} == {
            "contact_fraction_profile",
            "cf_by_aa_class_bars",
            "cf_by_partition_regions_bars",
        }

    @pytest.mark.parametrize(
        ("flag_name", "disabled_stem"),
        [
            ("generate_contact_fraction_profile", "contact_fraction_profile"),
            ("generate_residence_time_profile", "residence_time_profile"),
            ("generate_cf_by_aa_class_bars", "cf_by_aa_class_bars"),
            ("generate_cf_by_partition_bars", "cf_by_partition_regions_bars"),
            ("generate_rt_by_aa_class_bars", "rt_by_aa_class_bars"),
            ("generate_rt_by_partition_bars", "rt_by_partition_regions_bars"),
        ],
    )
    def test_plot_respects_disabled_plot_family_flags(self, tmp_path, flag_name, disabled_stem):
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.config.comparison import PlotSettings

        settings = ContactsSettings(
            protein_groups={"active_site": [1, 2], "surface": [3]},
            protein_partitions={"regions": ["active_site", "surface"]},
        )
        plot_settings = PlotSettings(contacts={flag_name: False})
        ctx = _contacts_plot_context(tmp_path, settings, plot_settings=plot_settings)
        _write_contacts_plot_artifact(ctx.analysis_dirs["A"], settings)

        paths = ContactsAnalysis().plot(ctx)

        stems = {path.stem for path in paths}
        assert disabled_stem not in stems
        assert stems == {
            "contact_fraction_profile",
            "residence_time_profile",
            "cf_by_aa_class_bars",
            "cf_by_partition_regions_bars",
            "rt_by_aa_class_bars",
            "rt_by_partition_regions_bars",
        } - {disabled_stem}

    def test_plot_skips_artifact_loading_when_all_plot_families_disabled(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.config.comparison import PlotSettings

        settings = ContactsSettings()
        plot_settings = PlotSettings(
            contacts={
                "generate_contact_fraction_profile": False,
                "generate_residence_time_profile": False,
                "generate_cf_by_aa_class_bars": False,
                "generate_cf_by_partition_bars": False,
                "generate_rt_by_aa_class_bars": False,
                "generate_rt_by_partition_bars": False,
            }
        )
        ctx = _contacts_plot_context(tmp_path, settings, plot_settings=plot_settings)

        with patch(
            "polyzymd.analyses.contacts._plotters.load_contacts_plot_data",
            side_effect=AssertionError("artifact loading should be skipped"),
        ):
            paths = ContactsAnalysis().plot(ctx)

        assert paths == []
        assert not ctx.output_dir.exists()

    def test_plot_skips_artifact_loading_when_only_rt_families_enabled_but_rt_disabled(
        self,
        tmp_path,
    ):
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.config.comparison import PlotSettings

        settings = ContactsSettings(compute_residence_times=False)
        plot_settings = PlotSettings(
            contacts={
                "generate_contact_fraction_profile": False,
                "generate_residence_time_profile": True,
                "generate_cf_by_aa_class_bars": False,
                "generate_cf_by_partition_bars": False,
                "generate_rt_by_aa_class_bars": True,
                "generate_rt_by_partition_bars": True,
            }
        )
        ctx = _contacts_plot_context(tmp_path, settings, plot_settings=plot_settings)

        with patch(
            "polyzymd.analyses.contacts._plotters.load_contacts_plot_data",
            side_effect=AssertionError("artifact loading should be skipped"),
        ):
            paths = ContactsAnalysis().plot(ctx)

        assert paths == []
        assert not ctx.output_dir.exists()

    def test_plot_rejects_tampered_profile_sidecar_without_fallback(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        settings = ContactsSettings()
        ctx = _contacts_plot_context(tmp_path, settings)
        _write_contacts_plot_artifact(ctx.analysis_dirs["A"], settings)
        profile_path = ctx.analysis_dirs["A"] / "aggregated" / "sidecars" / "contact_profiles.npz"
        profile_path.write_bytes(b"stale")

        with pytest.raises(ValueError, match="invalid profile sidecar"):
            ContactsAnalysis().plot(ctx)

    def test_plot_rejects_missing_profile_sidecar_without_fallback(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        settings = ContactsSettings()
        ctx = _contacts_plot_context(tmp_path, settings)
        _write_contacts_plot_artifact(ctx.analysis_dirs["A"], settings)
        profile_path = ctx.analysis_dirs["A"] / "aggregated" / "sidecars" / "contact_profiles.npz"
        profile_path.unlink()

        with pytest.raises(ValueError, match="invalid profile sidecar"):
            ContactsAnalysis().plot(ctx)

    def test_noncanonical_only_json_does_not_plot(self, tmp_path, caplog):
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        settings = ContactsSettings()
        ctx = _contacts_plot_context(tmp_path, settings)
        aggregate_path = ctx.analysis_dirs["A"] / "aggregated" / "result.json"
        aggregate_path.parent.mkdir(parents=True)
        aggregate_path.write_text('{"analysis_type": "contacts_aggregated", "residue_stats": []}')

        with caplog.at_level("INFO", logger="polyzymd.analyses.contacts"):
            paths = ContactsAnalysis().plot(ctx)

        assert paths == []
        assert "not a canonical ConditionArtifact" in caplog.text

    def test_plot_does_not_attempt_trajectory_loading(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.shared.loader import TrajectoryLoader

        settings = ContactsSettings()
        ctx = _contacts_plot_context(tmp_path, settings)
        _write_contacts_plot_artifact(ctx.analysis_dirs["A"], settings)

        with patch.object(TrajectoryLoader, "load_universe", side_effect=AssertionError):
            paths = ContactsAnalysis().plot(ctx)

        assert paths


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

        groups, partitions = _load_partition_definitions(settings, all_resids={1, 2, 3, 4, 5, 6})

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

        groups, partitions = _load_partition_definitions(settings, all_resids={1, 2, 3, 4, 5, 6})

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
        groups_1, partitions_1 = _load_partition_definitions(
            settings, all_resids={1, 2, 3, 4, 5, 6}
        )
        groups_2, partitions_2 = _load_partition_definitions(
            settings, all_resids={1, 2, 3, 4, 5, 6}
        )

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
        assert callable(analysis.build_mda_jobs)
        assert callable(analysis.build_mda_collector)
        assert callable(analysis.aggregate)
        assert callable(analysis.compare)
        assert callable(analysis.plot)
        assert callable(analysis.filter_conditions)


class TestContactsCanonicalArtifacts:
    """Contacts uses canonical artifact store paths only."""

    def test_noncanonical_cache_wrappers_are_not_public_on_plugin(self):
        """Non-canonical sidecar/cache path wrappers should not remain on the facade."""
        from polyzymd.analyses.contacts import ContactsAnalysis

        facade_names = set(ContactsAnalysis.__dict__)

        assert not {name for name in facade_names if "sidecar" in name and "path" in name}
        assert not {name for name in facade_names if "cache" in name and "candidate" in name}
        assert not {name for name in facade_names if "cache" in name and "context" in name}
        assert not {name for name in facade_names if name.startswith("_compute_contacts")}
        assert not {name for name in facade_names if name.startswith("_apply_fdr")}

    def test_condition_artifact_loads_from_canonical_result_json(self, tmp_path):
        """Comparison should read ``aggregated/result.json`` through ArtifactStore."""
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.analyses.mda import ArtifactStore

        settings = ContactsSettings()
        condition = Condition(
            label="A",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        analysis_dir = tmp_path / "A" / "contacts"
        aggregated_dir = analysis_dir / "aggregated"
        aggregated_dir.mkdir(parents=True)
        # Simulate framework-owned condition artifact persistence for compare
        ArtifactStore(aggregated_dir).write_condition_result(
            _make_condition_artifact("A", settings, replicates=(1, 2))
        )
        ctx = ComparisonContext(
            name="test",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"A": analysis_dir},
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=settings,
        )

        result = ContactsAnalysis().compare(ctx)

        assert result is not None
        assert [summary.label for summary in result.conditions] == ["A"]

    def test_noncanonical_sidecar_json_is_ignored_when_canonical_artifact_is_missing(
        self, tmp_path
    ):
        """Comparison should not discover non-canonical contacts JSON sidecars."""
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        settings = ContactsSettings()
        condition = Condition(
            label="A",
            config_path=tmp_path / "config.yaml",
            replicates=(1,),
            sim_config=MagicMock(),
        )
        analysis_dir = tmp_path / "A" / "contacts"
        analysis_dir.mkdir(parents=True)
        (analysis_dir / "contacts_eq10ns_cut4.5_rep1.json").write_text("{}")
        ctx = ComparisonContext(
            name="test",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"A": analysis_dir},
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=settings,
        )

        assert ContactsAnalysis().compare(ctx) is None
