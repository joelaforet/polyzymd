"""Tests for the secondary structure analysis plugin.

Tests the SecondaryStructureAnalysis class: discovery, aliases, settings,
compute_replicate, aggregate, extract_metrics, AggregatedResultClass,
_make_aggregated_filename, plot delegation, and the full lifecycle via
the orchestrator.

Heavy dependencies (mdtraj, trajectories) are mocked.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from pydantic import BaseModel

from polyzymd.analyses.base import (
    AggregateContext,
    ComparisonContext,
    Condition,
    MetricValue,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.secondary_structure import (
    SecondaryStructureAnalysis,
    SecondaryStructureSettings,
)

# ============================================================================
# Fixtures
# ============================================================================


@pytest.fixture
def ss_analysis():
    """Return a fresh SecondaryStructureAnalysis instance."""
    return SecondaryStructureAnalysis()


@pytest.fixture
def default_settings():
    """Return default SecondaryStructureSettings."""
    return SecondaryStructureSettings()


@pytest.fixture
def custom_settings():
    """Return SecondaryStructureSettings with a custom chain_id."""
    return SecondaryStructureSettings(chain_id="B")


@pytest.fixture
def condition():
    """Return a mock Condition."""
    return Condition(
        label="No Polymer",
        config_path=Path("/fake/config.yaml"),
        replicates=(1, 2, 3),
        sim_config=MagicMock(),
    )


def _make_mock_ss_result(
    replicate: int,
    overall_helix: float = 0.22,
    overall_strand: float = 0.18,
    overall_coil: float = 0.60,
    n_frames: int = 4800,
    n_residues: int = 5,
) -> MagicMock:
    """Create a mock SecondaryStructureResult for one replicate."""
    result = MagicMock()
    result.replicate = replicate
    result.n_frames = n_frames
    result.n_residues = n_residues
    result.residue_ids = list(range(1, n_residues + 1))
    result.residue_names = ["ALA"] * n_residues
    result.config_hash = "abc123"
    result.equilibration_time = 200.0
    result.equilibration_unit = "ns"
    result.selection_string = "chainid 0 (chain A)"

    # Per-residue persistence: slightly different per replicate
    rng = np.random.default_rng(seed=replicate)
    result.persistence_helix = (rng.uniform(0.15, 0.30, n_residues)).tolist()
    result.persistence_strand = (rng.uniform(0.10, 0.25, n_residues)).tolist()
    # Coil = 1 - helix - strand
    result.persistence_coil = [
        1.0 - h - s for h, s in zip(result.persistence_helix, result.persistence_strand)
    ]

    # Overall fractions (with slight variation per replicate)
    result.overall_helix_fraction = overall_helix + 0.005 * (replicate - 2)
    result.overall_strand_fraction = overall_strand + 0.003 * (replicate - 2)
    result.overall_coil_fraction = (
        1.0 - result.overall_helix_fraction - result.overall_strand_fraction
    )

    return result


# ============================================================================
# Discovery & class-level tests
# ============================================================================


class TestDiscovery:
    """Test that the plugin is discoverable by the framework."""

    def test_discovery_by_name(self):
        """The framework should find SecondaryStructureAnalysis by name."""
        from polyzymd.analyses.discovery import get_analysis

        cls = get_analysis("secondary_structure")
        assert cls is SecondaryStructureAnalysis

    def test_discovery_by_alias(self):
        """The framework should find SecondaryStructureAnalysis by alias 'ss'."""
        from polyzymd.analyses.discovery import get_analysis

        cls = get_analysis("ss")
        assert cls is SecondaryStructureAnalysis

    def test_listed(self):
        """SecondaryStructureAnalysis should appear in list_analyses()."""
        from polyzymd.analyses.discovery import list_analyses

        all_analyses = list_analyses()
        assert "secondary_structure" in all_analyses

    def test_all_names_include_alias(self):
        """list_all_names() should include both 'secondary_structure' and 'ss'."""
        from polyzymd.analyses.discovery import list_all_names

        names = list_all_names()
        assert "secondary_structure" in names
        assert "ss" in names


class TestClassAttributes:
    """Test class-level attributes are set correctly."""

    def test_name(self):
        assert SecondaryStructureAnalysis.name == "secondary_structure"

    def test_settings_type(self):
        assert SecondaryStructureAnalysis.Settings is SecondaryStructureSettings

    def test_aliases(self):
        assert SecondaryStructureAnalysis.aliases == ("ss",)

    def test_dependencies(self):
        assert SecondaryStructureAnalysis.dependencies == ()

    def test_min_replicates(self):
        assert SecondaryStructureAnalysis.min_replicates == 2


# ============================================================================
# Settings tests
# ============================================================================


class TestSettings:
    """Test SecondaryStructureSettings model."""

    def test_defaults(self, default_settings):
        assert default_settings.chain_id == "A"

    def test_custom_chain_id(self, custom_settings):
        assert custom_settings.chain_id == "B"

    def test_from_dict(self):
        settings = SecondaryStructureSettings(**{"chain_id": "C"})
        assert settings.chain_id == "C"

    def test_serialization_roundtrip(self, default_settings):
        d = default_settings.model_dump()
        restored = SecondaryStructureSettings(**d)
        assert restored.chain_id == default_settings.chain_id


# ============================================================================
# compute_replicate tests
# ============================================================================


class TestComputeReplicate:
    """Test compute_replicate performs inline DSSP computation."""

    def _make_mock_mdtraj_traj(self, n_frames: int = 200, n_residues: int = 5, n_atoms: int = 50):
        """Create a mock mdtraj Trajectory object."""
        mock_traj = MagicMock()
        mock_traj.n_frames = n_frames
        mock_traj.n_atoms = n_atoms
        mock_traj.n_residues = n_residues

        # Mock topology
        mock_top = MagicMock()
        mock_top.select.return_value = np.arange(n_atoms)
        mock_residues = []
        for i in range(n_residues):
            res = MagicMock()
            res.resSeq = i + 1
            res.name = "ALA"
            mock_residues.append(res)
        mock_top.residues = mock_residues
        mock_traj.topology = mock_top

        # atom_slice returns another mock trajectory with same topology
        sliced = MagicMock()
        sliced.n_frames = n_frames
        sliced.n_atoms = n_atoms
        sliced.n_residues = n_residues
        sliced.topology = mock_top
        # Slicing (for equilibration skip) returns another mock
        sliced.__getitem__ = MagicMock(return_value=sliced)
        mock_traj.atom_slice.return_value = sliced

        return mock_traj, sliced

    @patch("polyzymd.analyses.shared.config_hash.validate_config_hash")
    @patch("polyzymd.analyses.secondary_structure.compute_config_hash", return_value="abc123")
    @patch("polyzymd.analyses.secondary_structure.TrajectoryLoader")
    def test_computes_ss_inline(
        self,
        MockLoader,
        mock_hash,
        mock_validate_hash,
        ss_analysis,
        default_settings,
        condition,
        tmp_path,
    ):
        """compute_replicate should perform DSSP computation inline."""
        mock_loader_inst = MagicMock()
        MockLoader.return_value = mock_loader_inst

        traj_info = MagicMock()
        traj_info.trajectory_files = [Path("/fake/traj.dcd")]
        traj_info.topology_file = Path("/fake/topology.pdb")
        mock_loader_inst.get_trajectory_info.return_value = traj_info
        mock_loader_inst.get_timestep.return_value = 10.0  # 10 ps

        mock_traj, mock_protein_traj = self._make_mock_mdtraj_traj(n_frames=200, n_residues=5)

        # Create DSSP output: (n_frames, n_residues) of characters
        dssp_raw = np.array([["H", "E", "C", "H", "C"]] * 200)
        # Encode: H=1, E=2, C=0
        ss_matrix = np.zeros((200, 5), dtype=np.int8)
        ss_matrix[:, 0] = 1  # H
        ss_matrix[:, 1] = 2  # E
        ss_matrix[:, 2] = 0  # C
        ss_matrix[:, 3] = 1  # H
        ss_matrix[:, 4] = 0  # C

        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="200ns",
            recompute=False,
            settings=default_settings,
        )

        with (
            patch(
                "mdtraj.load",
                return_value=mock_traj,
            ),
            patch(
                "mdtraj.compute_dssp",
                return_value=dssp_raw,
            ),
            patch(
                "polyzymd.analyses.secondary_structure._encode_dssp_matrix",
                return_value=ss_matrix,
            ),
            patch("polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.2.1"),
        ):
            result = ss_analysis.compute_replicate(ctx, replicate=1)

        # Verify result has SS data
        assert result.replicate == 1
        assert result.n_residues == 5
        assert result.overall_helix_fraction == pytest.approx(0.4)  # 2 of 5 residues
        assert result.overall_strand_fraction == pytest.approx(0.2)  # 1 of 5
        assert result.overall_coil_fraction == pytest.approx(0.4)  # 2 of 5
        # Verify TrajectoryLoader was used
        MockLoader.assert_called_once_with(condition.sim_config)

    @patch("polyzymd.analyses.shared.config_hash.validate_config_hash")
    @patch("polyzymd.analyses.secondary_structure.compute_config_hash", return_value="abc123")
    @patch("polyzymd.analyses.secondary_structure.TrajectoryLoader")
    def test_custom_chain_id(
        self,
        MockLoader,
        mock_hash,
        mock_validate_hash,
        ss_analysis,
        custom_settings,
        condition,
        tmp_path,
    ):
        """compute_replicate should use the custom chain_id for chain selection."""
        mock_loader_inst = MagicMock()
        MockLoader.return_value = mock_loader_inst

        traj_info = MagicMock()
        traj_info.trajectory_files = [Path("/fake/traj.dcd")]
        traj_info.topology_file = Path("/fake/topology.pdb")
        mock_loader_inst.get_trajectory_info.return_value = traj_info
        mock_loader_inst.get_timestep.return_value = 10.0

        mock_traj, mock_protein_traj = self._make_mock_mdtraj_traj(n_frames=200, n_residues=5)

        dssp_raw = np.array([["C", "C", "C", "C", "C"]] * 200)
        ss_matrix = np.zeros((200, 5), dtype=np.int8)

        ctx = ReplicateContext(
            condition=condition,
            replicate=2,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_2",
            equilibration="100ns",
            recompute=True,
            settings=custom_settings,
        )

        with (
            patch(
                "mdtraj.load",
                return_value=mock_traj,
            ),
            patch(
                "mdtraj.compute_dssp",
                return_value=dssp_raw,
            ),
            patch(
                "polyzymd.analyses.secondary_structure._encode_dssp_matrix",
                return_value=ss_matrix,
            ),
            patch("polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.2.1"),
        ):
            result = ss_analysis.compute_replicate(ctx, replicate=2)

        # Verify chain B (index 1) was used for selection
        mock_traj.topology.select.assert_any_call("chainid 1")
        assert result.selection_string == "chainid 1 (chain B)"


# ============================================================================
# aggregate tests
# ============================================================================


class TestAggregate:
    """Test aggregation of per-replicate SS results."""

    def test_aggregate_produces_correct_result(
        self, ss_analysis, default_settings, condition, tmp_path
    ):
        """aggregate() should compute mean/SEM across replicates."""
        results = [_make_mock_ss_result(r) for r in (1, 2, 3)]

        agg_ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2, 3),
            output_dir=tmp_path / "aggregated",
            equilibration="200ns",
            settings=default_settings,
        )

        with patch(
            "polyzymd.analyses._results_base.get_polyzymd_version",
            return_value="0.0.0-test",
        ):
            agg = ss_analysis.aggregate(agg_ctx, results)

        # Check basic fields
        assert agg.n_replicates == 3
        assert agg.replicates == [1, 2, 3]
        assert agg.config_hash == "abc123"
        assert agg.selection_string == results[0].selection_string

        # Check mean overall values are the mean of per-replicate values
        expected_helix_mean = np.mean([r.overall_helix_fraction for r in results])
        assert abs(agg.mean_overall_helix - expected_helix_mean) < 1e-10

        expected_strand_mean = np.mean([r.overall_strand_fraction for r in results])
        assert abs(agg.mean_overall_strand - expected_strand_mean) < 1e-10

        # Check SEM values
        expected_helix_sem = float(
            np.std([r.overall_helix_fraction for r in results], ddof=1) / np.sqrt(3)
        )
        assert abs(agg.sem_overall_helix - expected_helix_sem) < 1e-10

        # Per-residue persistence lengths should match
        assert len(agg.mean_persistence_helix) == results[0].n_residues
        assert len(agg.sem_persistence_helix) == results[0].n_residues

        # Per-replicate values should be stored
        assert len(agg.per_replicate_helix) == 3
        assert len(agg.per_replicate_strand) == 3
        assert len(agg.per_replicate_coil) == 3

    def test_aggregate_saves_result(self, ss_analysis, default_settings, condition, tmp_path):
        """aggregate() should save the result JSON to disk."""
        results = [_make_mock_ss_result(r) for r in (1, 2, 3)]

        agg_ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2, 3),
            output_dir=tmp_path / "aggregated",
            equilibration="200ns",
            settings=default_settings,
        )

        with patch(
            "polyzymd.analyses._results_base.get_polyzymd_version",
            return_value="0.0.0-test",
        ):
            ss_analysis.aggregate(agg_ctx, results)

        # Check that a JSON file was saved
        json_files = list((tmp_path / "aggregated").glob("*.json"))
        assert len(json_files) == 1
        assert "secondary_structure" in json_files[0].name

    def test_aggregate_two_replicates(self, ss_analysis, default_settings, condition, tmp_path):
        """aggregate() should work with the minimum of 2 replicates."""
        results = [_make_mock_ss_result(r) for r in (1, 2)]

        agg_ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2),
            output_dir=tmp_path / "aggregated",
            equilibration="200ns",
            settings=default_settings,
        )

        with patch(
            "polyzymd.analyses._results_base.get_polyzymd_version",
            return_value="0.0.0-test",
        ):
            agg = ss_analysis.aggregate(agg_ctx, results)

        assert agg.n_replicates == 2


# ============================================================================
# extract_metrics tests
# ============================================================================


class TestExtractMetrics:
    """Test metric extraction for the default scalar comparison."""

    def test_extracts_helix_fraction(self, ss_analysis):
        """extract_metrics() should return helix_fraction MetricValue."""
        summary = MagicMock()
        summary.mean_overall_helix = 0.223
        summary.sem_overall_helix = 0.008
        summary.per_replicate_helix = [0.215, 0.225, 0.229]

        metrics = ss_analysis.extract_metrics(summary)

        assert "helix_fraction" in metrics
        m = metrics["helix_fraction"]
        assert isinstance(m, MetricValue)
        assert m.name == "helix_fraction"
        assert m.mean == 0.223
        assert m.sem == 0.008
        assert m.replicate_values == [0.215, 0.225, 0.229]

    def test_higher_is_better(self, ss_analysis):
        """Helix fraction should be higher_is_better=True."""
        summary = MagicMock()
        summary.mean_overall_helix = 0.2
        summary.sem_overall_helix = 0.01
        summary.per_replicate_helix = [0.19, 0.21]

        metrics = ss_analysis.extract_metrics(summary)
        assert metrics["helix_fraction"].higher_is_better is True

    def test_direction_labels(self, ss_analysis):
        """Direction labels should be destabilizing/unchanged/stabilizing."""
        summary = MagicMock()
        summary.mean_overall_helix = 0.2
        summary.sem_overall_helix = 0.01
        summary.per_replicate_helix = [0.19, 0.21]

        metrics = ss_analysis.extract_metrics(summary)
        assert metrics["helix_fraction"].direction_labels == (
            "destabilizing",
            "unchanged",
            "stabilizing",
        )

    def test_single_metric_returned(self, ss_analysis):
        """SS uses the default compare() so should return exactly 1 metric."""
        summary = MagicMock()
        summary.mean_overall_helix = 0.2
        summary.sem_overall_helix = 0.01
        summary.per_replicate_helix = [0.19, 0.21]

        metrics = ss_analysis.extract_metrics(summary)
        assert len(metrics) == 1


# ============================================================================
# AggregatedResultClass and _deserialize_result tests
# ============================================================================


class TestDeserializeResult:
    """Test loading aggregated results via AggregatedResultClass."""

    def test_aggregated_result_class_set(self, ss_analysis):
        """AggregatedResultClass should be SecondaryStructureAggregatedResult."""
        from polyzymd.analyses.secondary_structure._results import (
            SecondaryStructureAggregatedResult,
        )

        assert ss_analysis.AggregatedResultClass is SecondaryStructureAggregatedResult

    def test_loads_via_result_class(self, ss_analysis, tmp_path):
        """_deserialize_result should use AggregatedResultClass.load()."""
        fake_path = tmp_path / "agg.json"

        mock_loaded = MagicMock()
        with patch.object(
            ss_analysis.AggregatedResultClass, "load", return_value=mock_loaded
        ) as mock_load:
            result = ss_analysis._deserialize_result(fake_path)

            mock_load.assert_called_once_with(fake_path)
            assert result is mock_loaded


# ============================================================================
# _make_aggregated_filename tests
# ============================================================================


class TestMakeAggregatedFilename:
    """Test filename generation for aggregated results."""

    def test_consecutive_replicates(self, ss_analysis):
        """Consecutive replicates should use range format."""
        first = MagicMock()
        first.equilibration_time = 200.0
        first.equilibration_unit = "ns"

        name = ss_analysis._make_aggregated_filename((1, 2, 3, 4, 5), first)
        assert name == "secondary_structure_reps1-5_eq200.00ns.json"

    def test_nonconsecutive_replicates(self, ss_analysis):
        """Non-consecutive replicates should list individually."""
        first = MagicMock()
        first.equilibration_time = 100.0
        first.equilibration_unit = "ns"

        name = ss_analysis._make_aggregated_filename((1, 3, 5), first)
        assert name == "secondary_structure_reps1_3_5_eq100.00ns.json"

    def test_ps_equilibration(self, ss_analysis):
        """Equilibration in ps should be formatted correctly."""
        first = MagicMock()
        first.equilibration_time = 5000.0
        first.equilibration_unit = "ps"

        name = ss_analysis._make_aggregated_filename((1, 2), first)
        assert name == "secondary_structure_reps1-2_eq5000.00ps.json"


# ============================================================================
# plot tests
# ============================================================================


class TestPlot:
    """Test plot delegation to inline plotting helpers."""

    def _make_plot_ctx(self, tmp_path, conditions=None):
        """Helper to build a PlotContext."""
        from polyzymd.config.comparison import PlotSettings

        if conditions is None:
            conditions = [
                Condition(
                    label="No Polymer",
                    config_path=Path("/fake/config1.yaml"),
                    replicates=(1, 2, 3),
                    sim_config=MagicMock(),
                ),
                Condition(
                    label="100% SBMA",
                    config_path=Path("/fake/config2.yaml"),
                    replicates=(1, 2, 3),
                    sim_config=MagicMock(),
                ),
            ]

        analysis_dirs = {}
        for cond in conditions:
            d = tmp_path / "analysis" / cond.label / "secondary_structure"
            d.mkdir(parents=True)
            analysis_dirs[cond.label] = d

        return PlotContext(
            conditions=conditions,
            analysis_dirs=analysis_dirs,
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=SecondaryStructureSettings(),
            plot_settings=PlotSettings(),
        )

    def test_delegates_to_all_four_plotters(self, ss_analysis, tmp_path):
        """plot() should attempt all 4 SS plotters."""
        ctx = self._make_plot_ctx(tmp_path)

        with (
            patch(
                "polyzymd.analyses.secondary_structure._plot_ss_timeline_heatmap"
            ) as MockTimeline,
            patch("polyzymd.analyses.secondary_structure._plot_ss_content_bars") as MockContent,
            patch(
                "polyzymd.analyses.secondary_structure._plot_ss_individual_bars"
            ) as MockIndividual,
            patch(
                "polyzymd.analyses.secondary_structure._plot_ss_persistence_diff_heatmap"
            ) as MockPersistence,
        ):
            # All plotters return a single path
            for mock_fn in [MockTimeline, MockContent, MockIndividual, MockPersistence]:
                mock_fn.return_value = [tmp_path / "plots" / "fake.png"]

            plots = ss_analysis.plot(ctx)

            # All 4 plotters should have been called
            assert MockTimeline.called
            assert MockContent.called
            assert MockIndividual.called
            assert MockPersistence.called

            # Total paths: 4 (one from each plotter)
            assert len(plots) == 4

    def test_plot_handles_plotter_failure(self, ss_analysis, tmp_path):
        """plot() should catch exceptions from individual plotters."""
        ctx = self._make_plot_ctx(tmp_path)

        with (
            patch(
                "polyzymd.analyses.secondary_structure._plot_ss_timeline_heatmap"
            ) as MockTimeline,
            patch("polyzymd.analyses.secondary_structure._plot_ss_content_bars") as MockContent,
            patch(
                "polyzymd.analyses.secondary_structure._plot_ss_individual_bars"
            ) as MockIndividual,
            patch(
                "polyzymd.analyses.secondary_structure._plot_ss_persistence_diff_heatmap"
            ) as MockPersistence,
        ):
            # Timeline raises, others succeed
            MockTimeline.side_effect = RuntimeError("boom")
            MockContent.return_value = [tmp_path / "plots" / "bars.png"]
            MockIndividual.return_value = [tmp_path / "plots" / "ind.png"]
            MockPersistence.return_value = [tmp_path / "plots" / "diff.png"]

            plots = ss_analysis.plot(ctx)

            # Should still get 3 paths from the working plotters
            assert len(plots) == 3

    def test_plot_empty_conditions(self, ss_analysis, tmp_path):
        """plot() with no conditions should return empty list."""
        ctx = self._make_plot_ctx(tmp_path, conditions=[])
        plots = ss_analysis.plot(ctx)
        assert plots == []

    def test_plot_creates_default_plot_settings(self, ss_analysis, tmp_path):
        """Orchestrator guarantees non-None plot_settings; plugin passes it through."""
        ctx = self._make_plot_ctx(tmp_path)

        with (
            patch(
                "polyzymd.analyses.secondary_structure._plot_ss_timeline_heatmap"
            ) as MockTimeline,
            patch("polyzymd.analyses.secondary_structure._plot_ss_content_bars") as MockContent,
            patch(
                "polyzymd.analyses.secondary_structure._plot_ss_individual_bars"
            ) as MockIndividual,
            patch(
                "polyzymd.analyses.secondary_structure._plot_ss_persistence_diff_heatmap"
            ) as MockPersistence,
        ):
            for mock_fn in [MockTimeline, MockContent, MockIndividual, MockPersistence]:
                mock_fn.return_value = []

            ss_analysis.plot(ctx)

            # Check that PlotSettings was created and passed as 4th arg
            from polyzymd.config.comparison import PlotSettings

            for mock_fn in [MockTimeline, MockContent, MockIndividual, MockPersistence]:
                assert mock_fn.called
                call_args = mock_fn.call_args
                assert call_args is not None
                assert isinstance(call_args[0][3], PlotSettings)

    def test_plot_builds_data_dict_correctly(self, ss_analysis, tmp_path):
        """plot() should build the data dict with analysis_dir and replicates."""
        conditions = [
            Condition(
                label="Test",
                config_path=Path("/fake/config.yaml"),
                replicates=(1, 2),
                sim_config=MagicMock(),
            ),
        ]
        ctx = self._make_plot_ctx(tmp_path, conditions=conditions)

        captured_data = {}

        with (
            patch(
                "polyzymd.analyses.secondary_structure._plot_ss_timeline_heatmap"
            ) as MockTimeline,
            patch("polyzymd.analyses.secondary_structure._plot_ss_content_bars"),
            patch("polyzymd.analyses.secondary_structure._plot_ss_individual_bars"),
            patch("polyzymd.analyses.secondary_structure._plot_ss_persistence_diff_heatmap"),
        ):

            def capture_plot(data, labels, output_dir, plot_settings):
                captured_data.update(data)
                return []

            MockTimeline.side_effect = capture_plot

            ss_analysis.plot(ctx)

        assert "Test" in captured_data
        assert "analysis_dir" in captured_data["Test"]
        assert "aggregated_dir" in captured_data["Test"]
        assert "replicates" in captured_data["Test"]
        assert captured_data["Test"]["replicates"] == [1, 2]
        assert "__meta__" in captured_data


# ============================================================================
# compare (default) tests
# ============================================================================


class TestCompare:
    """Test the default scalar comparison via extract_metrics."""

    def test_compare_calls_default_scalar_comparison(self, ss_analysis, tmp_path):
        """compare() should delegate to default_scalar_comparison."""
        cond1 = Condition(
            label="Control",
            config_path=Path("/fake/c1.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        cond2 = Condition(
            label="Treatment",
            config_path=Path("/fake/c2.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )

        # Create aggregated dirs with fake results
        for label in ("Control", "Treatment"):
            d = tmp_path / label / "aggregated"
            d.mkdir(parents=True)

        ctx = ComparisonContext(
            name="test_compare",
            conditions=[cond1, cond2],
            excluded_conditions=[],
            control_label="Control",
            analysis_dirs={
                "Control": tmp_path / "Control",
                "Treatment": tmp_path / "Treatment",
            },
            results_dir=tmp_path / "results",
            equilibration="200ns",
            settings=SecondaryStructureSettings(),
            recompute=False,
        )

        mock_agg_control = MagicMock()
        mock_agg_control.mean_overall_helix = 0.22
        mock_agg_control.sem_overall_helix = 0.008
        mock_agg_control.per_replicate_helix = [0.215, 0.225, 0.220]

        mock_agg_treatment = MagicMock()
        mock_agg_treatment.mean_overall_helix = 0.25
        mock_agg_treatment.sem_overall_helix = 0.005
        mock_agg_treatment.per_replicate_helix = [0.245, 0.255, 0.250]

        agg_map = {
            str(tmp_path / "Control" / "aggregated"): mock_agg_control,
            str(tmp_path / "Treatment" / "aggregated"): mock_agg_treatment,
        }

        def mock_load_agg(self_inner, agg_dir):
            return agg_map.get(str(agg_dir))

        with (
            patch.object(
                SecondaryStructureAnalysis,
                "_load_aggregated_result",
                mock_load_agg,
            ),
            patch("polyzymd.analyses.stats.default_scalar_comparison") as mock_dsc,
        ):
            mock_dsc.return_value = MagicMock()
            result = ss_analysis.compare(ctx)

            mock_dsc.assert_called_once()
            call_kwargs = mock_dsc.call_args.kwargs
            assert call_kwargs["analysis_name"] == "secondary_structure"
            assert call_kwargs["project_name"] == "test_compare"
            assert call_kwargs["control_label"] == "Control"
            assert "Control" in call_kwargs["metrics_by_condition"]
            assert "Treatment" in call_kwargs["metrics_by_condition"]

    def test_compare_returns_result_with_single_condition(self, ss_analysis, tmp_path):
        """compare() should return a result with one condition having metrics."""
        cond = Condition(
            label="Only",
            config_path=Path("/fake/c.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )

        (tmp_path / "Only" / "aggregated").mkdir(parents=True)

        ctx = ComparisonContext(
            name="test",
            conditions=[cond],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"Only": tmp_path / "Only"},
            results_dir=tmp_path / "results",
            equilibration="200ns",
            settings=SecondaryStructureSettings(),
            recompute=False,
        )

        mock_agg = MagicMock()
        mock_agg.mean_overall_helix = 0.22
        mock_agg.sem_overall_helix = 0.008
        mock_agg.per_replicate_helix = [0.215, 0.225, 0.220]

        with patch.object(
            SecondaryStructureAnalysis,
            "_load_aggregated_result",
            return_value=mock_agg,
        ):
            result = ss_analysis.compare(ctx)

        assert result is not None
        assert len(result.conditions) == 1
        assert result.pairwise_comparisons == []


# ============================================================================
# Full lifecycle test
# ============================================================================


class TestLifecycle:
    """Integration test: settings -> compute -> aggregate -> extract_metrics."""

    def test_settings_to_compute_to_aggregate_to_metrics(
        self, ss_analysis, default_settings, condition, tmp_path
    ):
        """Full pipeline: settings, compute, aggregate, extract_metrics work together."""
        # Step 1: Verify settings
        assert default_settings.chain_id == "A"

        # Step 2: Mock compute for 3 replicates
        mock_results = [_make_mock_ss_result(r) for r in (1, 2, 3)]

        # Step 3: Aggregate
        agg_ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2, 3),
            output_dir=tmp_path / "aggregated",
            equilibration="200ns",
            settings=default_settings,
        )

        with patch(
            "polyzymd.analyses._results_base.get_polyzymd_version",
            return_value="0.0.0-test",
        ):
            agg = ss_analysis.aggregate(agg_ctx, mock_results)

        # Step 4: Extract metrics
        metrics = ss_analysis.extract_metrics(agg)

        assert "helix_fraction" in metrics
        m = metrics["helix_fraction"]
        assert m.mean == agg.mean_overall_helix
        assert m.sem == agg.sem_overall_helix
        assert m.replicate_values == agg.per_replicate_helix
        assert m.higher_is_better is True

    def test_aggregated_result_can_be_loaded_back(
        self, ss_analysis, default_settings, condition, tmp_path
    ):
        """Saved aggregated result should be loadable via AggregatedResultClass."""
        results = [_make_mock_ss_result(r) for r in (1, 2, 3)]

        agg_ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2, 3),
            output_dir=tmp_path / "aggregated",
            equilibration="200ns",
            settings=default_settings,
        )

        with patch(
            "polyzymd.analyses._results_base.get_polyzymd_version",
            return_value="0.0.0-test",
        ):
            agg = ss_analysis.aggregate(agg_ctx, results)

        # Find the saved file
        json_files = list((tmp_path / "aggregated").glob("*.json"))
        assert len(json_files) == 1

        # Load it back
        loaded = ss_analysis._deserialize_result(json_files[0])

        assert loaded.n_replicates == 3
        assert abs(loaded.mean_overall_helix - agg.mean_overall_helix) < 1e-10
        assert abs(loaded.sem_overall_helix - agg.sem_overall_helix) < 1e-10
        assert loaded.per_replicate_helix == agg.per_replicate_helix


# ============================================================================
# filter_conditions tests
# ============================================================================


class TestFilterConditions:
    """Test that filter_conditions keeps all conditions (default behavior)."""

    def test_keeps_all_conditions(self, ss_analysis):
        conditions = [
            Condition(
                label="A",
                config_path=Path("/a.yaml"),
                replicates=(1,),
                sim_config=MagicMock(),
            ),
            Condition(
                label="B",
                config_path=Path("/b.yaml"),
                replicates=(1,),
                sim_config=MagicMock(),
            ),
        ]
        result = ss_analysis.filter_conditions(conditions)
        assert len(result) == 2
        assert result[0].label == "A"
        assert result[1].label == "B"
