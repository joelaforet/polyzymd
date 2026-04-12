"""Tests for the RMSF analysis plugin.

Tests the RMSFAnalysis class: discovery, settings, compute_replicate,
aggregate, extract_metrics, AggregatedResultClass, inlined plotting,
and the full lifecycle via the orchestrator.

Heavy dependencies (MDAnalysis, trajectories) are mocked.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from polyzymd.analyses.base import (
    AggregateContext,
    Condition,
    MetricValue,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.rmsf import RMSFAnalysis, RMSFSettings

# ============================================================================
# Fixtures
# ============================================================================


@pytest.fixture
def rmsf_analysis():
    """Return a fresh RMSFAnalysis instance."""
    return RMSFAnalysis()


@pytest.fixture
def default_settings():
    """Return default RMSFSettings."""
    return RMSFSettings()


@pytest.fixture
def condition():
    """Return a Condition test object."""
    return Condition(
        label="No Polymer",
        config_path=Path("/fake/config.yaml"),
        replicates=(1, 2, 3),
        sim_config=MagicMock(),
    )


def _make_mock_rmsf_result(replicate: int, mean_rmsf: float = 1.5) -> MagicMock:
    """Create a mock RMSFResult with realistic fields."""
    result = MagicMock()
    result.replicate = replicate
    result.rmsf_values = [1.0, 1.5, 2.0, 1.2, 1.8]
    result.residue_ids = [1, 2, 3, 4, 5]
    result.residue_names = ["ALA", "GLY", "VAL", "LEU", "ILE"]
    result.mean_rmsf = mean_rmsf
    result.std_rmsf = 0.35
    result.min_rmsf = 1.0
    result.max_rmsf = 2.0
    result.config_hash = "abc123"
    result.equilibration_time = 10.0
    result.equilibration_unit = "ns"
    result.selection_string = "protein and name CA"
    result.n_independent_frames = 50
    return result


# ============================================================================
# Test: Discovery
# ============================================================================


class TestRMSFDiscovery:
    """Test that RMSFAnalysis is auto-discovered by the plugin system."""

    def test_discovery_finds_rmsf(self):
        from polyzymd.analyses.discovery import clear_cache, list_analyses

        clear_cache()
        analyses = list_analyses()
        assert "rmsf" in analyses
        assert analyses["rmsf"] is RMSFAnalysis

    def test_get_analysis_returns_rmsf(self):
        from polyzymd.analyses.discovery import clear_cache, get_analysis

        clear_cache()
        cls = get_analysis("rmsf")
        assert cls is RMSFAnalysis


# ============================================================================
# Test: Class variables and Settings
# ============================================================================


class TestRMSFClassVars:
    """Test class variables and settings."""

    def test_name(self, rmsf_analysis):
        assert rmsf_analysis.name == "rmsf"

    def test_aliases_empty(self, rmsf_analysis):
        assert rmsf_analysis.aliases == ()

    def test_dependencies_empty(self, rmsf_analysis):
        assert rmsf_analysis.dependencies == ()

    def test_min_replicates(self, rmsf_analysis):
        assert rmsf_analysis.min_replicates == 2

    def test_repr(self, rmsf_analysis):
        assert repr(rmsf_analysis) == "<RMSFAnalysis(name='rmsf')>"


class TestRMSFSettings:
    """Test RMSFSettings validation and defaults."""

    def test_defaults(self):
        s = RMSFSettings()
        assert s.selection == "protein and name CA"
        assert s.reference_mode == "centroid"
        assert s.reference_frame is None
        assert s.reference_file is None

    def test_custom_selection(self):
        s = RMSFSettings(selection="protein and name N")
        assert s.selection == "protein and name N"

    def test_frame_mode(self):
        s = RMSFSettings(reference_mode="frame", reference_frame=500)
        assert s.reference_mode == "frame"
        assert s.reference_frame == 500

    def test_external_mode(self):
        s = RMSFSettings(reference_mode="external", reference_file="/path/to/ref.pdb")
        assert s.reference_mode == "external"

    def test_invalid_reference_mode(self):
        with pytest.raises(ValueError, match="reference_mode must be one of"):
            RMSFSettings(reference_mode="invalid")

    def test_average_mode(self):
        s = RMSFSettings(reference_mode="average")
        assert s.reference_mode == "average"


# ============================================================================
# Test: compute_replicate
# ============================================================================


class TestComputeReplicate:
    """Test RMSFAnalysis.compute_replicate performs inline RMSF computation."""

    def _make_mock_universe(self, n_frames: int = 200, n_atoms: int = 5):
        """Create a mock MDAnalysis Universe for RMSF tests."""
        import numpy as np

        mock_u = MagicMock()

        # Mock trajectory
        mock_traj = MagicMock()
        mock_traj.__len__ = MagicMock(return_value=n_frames)
        # Make trajectory iterable (for slicing and iteration)
        mock_frames = []
        for i in range(n_frames):
            ts = MagicMock()
            ts.frame = i
            mock_frames.append(ts)
        mock_traj.__getitem__ = MagicMock(
            side_effect=lambda x: mock_frames[x] if isinstance(x, int) else mock_frames
        )
        mock_traj.__iter__ = MagicMock(return_value=iter(mock_frames))
        mock_u.trajectory = mock_traj

        # Mock atom selection
        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=n_atoms)
        # Each call to .positions returns random positions
        mock_atoms.positions = np.random.rand(n_atoms, 3).astype(np.float32)
        mock_atoms.indices = np.arange(n_atoms)

        # Mock residues
        mock_residues = []
        for i in range(n_atoms):
            res = MagicMock()
            res.resid = i + 1
            res.resname = "ALA"
            mock_residues.append(res)
        mock_atoms.residues = mock_residues

        mock_u.select_atoms = MagicMock(return_value=mock_atoms)

        return mock_u

    @patch("polyzymd.analyses.rmsf.align_trajectory", return_value=0)
    @patch("polyzymd.analyses.rmsf.validate_equilibration_time", return_value=(True, ""))
    @patch("polyzymd.analyses.shared.config_hash.validate_config_hash")
    @patch("polyzymd.analyses.rmsf.compute_config_hash", return_value="abc123")
    @patch("polyzymd.analyses.rmsf.TrajectoryLoader")
    def test_computes_rmsf_inline(
        self,
        MockLoader,
        mock_hash,
        mock_validate_hash,
        mock_eq_validate,
        mock_align,
        rmsf_analysis,
        condition,
        tmp_path,
    ):
        """Test that compute_replicate performs RMSF calculation inline."""
        import numpy as np

        # Setup mock loader
        mock_loader_inst = MagicMock()
        MockLoader.return_value = mock_loader_inst

        mock_u = self._make_mock_universe(n_frames=200, n_atoms=5)
        mock_loader_inst.load_universe.return_value = mock_u

        traj_info = MagicMock()
        traj_info.trajectory_files = [Path("/fake/traj.dcd")]
        mock_loader_inst.get_trajectory_info.return_value = traj_info
        mock_loader_inst.get_timestep.return_value = 10.0  # 10 ps

        settings = RMSFSettings()
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=False,
            settings=settings,
        )

        # Patch _compute_rmsf and _compute_rmsd_timeseries to avoid real iteration
        with (
            patch(
                "polyzymd.analyses.rmsf._compute_rmsf",
                return_value=np.array([1.0, 1.5, 2.0, 1.2, 1.8]),
            ),
            patch(
                "polyzymd.analyses.rmsf._compute_rmsd_timeseries", return_value=np.random.rand(100)
            ),
            patch("polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.2.1"),
        ):
            result = rmsf_analysis.compute_replicate(ctx, 1)

        # Verify result has RMSF data
        assert result.replicate == 1
        assert len(result.rmsf_values) == 5
        assert result.selection_string == "protein and name CA"
        # Verify the loader was used (not the old calculator)
        MockLoader.assert_called_once_with(condition.sim_config)
        mock_loader_inst.load_universe.assert_called_once_with(1)

    @patch("polyzymd.analyses.rmsf.align_trajectory", return_value=0)
    @patch("polyzymd.analyses.rmsf.validate_equilibration_time", return_value=(True, ""))
    @patch("polyzymd.analyses.shared.config_hash.validate_config_hash")
    @patch("polyzymd.analyses.rmsf.compute_config_hash", return_value="abc123")
    @patch("polyzymd.analyses.rmsf.TrajectoryLoader")
    def test_passes_custom_settings(
        self,
        MockLoader,
        mock_hash,
        mock_validate_hash,
        mock_eq_validate,
        mock_align,
        rmsf_analysis,
        condition,
        tmp_path,
    ):
        """Custom settings (selection, reference_mode) are used by the inline computation."""
        import numpy as np

        mock_loader_inst = MagicMock()
        MockLoader.return_value = mock_loader_inst
        mock_u = self._make_mock_universe(n_frames=200, n_atoms=5)
        mock_loader_inst.load_universe.return_value = mock_u

        traj_info = MagicMock()
        traj_info.trajectory_files = [Path("/fake/traj.dcd")]
        mock_loader_inst.get_trajectory_info.return_value = traj_info
        mock_loader_inst.get_timestep.return_value = 10.0

        settings = RMSFSettings(
            selection="backbone",
            reference_mode="average",
        )
        ctx = ReplicateContext(
            condition=condition,
            replicate=2,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_2",
            equilibration="5ns",
            recompute=True,
            settings=settings,
        )

        with (
            patch(
                "polyzymd.analyses.rmsf._compute_rmsf",
                return_value=np.array([1.0, 1.5, 2.0, 1.2, 1.8]),
            ),
            patch(
                "polyzymd.analyses.rmsf._compute_rmsd_timeseries", return_value=np.random.rand(100)
            ),
            patch("polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.2.1"),
        ):
            result = rmsf_analysis.compute_replicate(ctx, 2)

        # Verify custom selection was used
        mock_u.select_atoms.assert_called_with("backbone")
        # Verify alignment was called with average mode
        mock_align.assert_called_once()
        call_args = mock_align.call_args
        alignment_config = call_args[0][1]
        assert alignment_config.reference_mode == "average"

    @patch("polyzymd.analyses.rmsf.align_trajectory", return_value=0)
    @patch("polyzymd.analyses.rmsf.validate_equilibration_time", return_value=(True, ""))
    @patch("polyzymd.analyses.shared.config_hash.validate_config_hash")
    @patch("polyzymd.analyses.rmsf.compute_config_hash", return_value="abc123")
    @patch("polyzymd.analyses.rmsf.TrajectoryLoader")
    def test_handles_legacy_settings(
        self,
        MockLoader,
        mock_hash,
        mock_validate_hash,
        mock_eq_validate,
        mock_align,
        rmsf_analysis,
        condition,
        tmp_path,
    ):
        """Legacy settings objects (e.g. MagicMock with attrs) should work via getattr."""
        import numpy as np

        mock_loader_inst = MagicMock()
        MockLoader.return_value = mock_loader_inst
        mock_u = self._make_mock_universe(n_frames=200, n_atoms=5)
        mock_loader_inst.load_universe.return_value = mock_u

        traj_info = MagicMock()
        traj_info.trajectory_files = [Path("/fake/traj.dcd")]
        mock_loader_inst.get_trajectory_info.return_value = traj_info
        mock_loader_inst.get_timestep.return_value = 10.0

        # Simulate legacy settings with same attributes
        legacy_settings = MagicMock()
        legacy_settings.selection = "protein and name CA"
        legacy_settings.reference_mode = "average"
        legacy_settings.reference_frame = None
        legacy_settings.reference_file = None
        legacy_settings.alignment_selection = "protein and name CA"
        legacy_settings.centroid_selection = "protein"

        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=False,
            settings=legacy_settings,
        )

        with (
            patch(
                "polyzymd.analyses.rmsf._compute_rmsf",
                return_value=np.array([1.0, 1.5, 2.0, 1.2, 1.8]),
            ),
            patch(
                "polyzymd.analyses.rmsf._compute_rmsd_timeseries", return_value=np.random.rand(100)
            ),
            patch("polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.2.1"),
        ):
            result = rmsf_analysis.compute_replicate(ctx, 1)

        # Verify it completed successfully with legacy settings
        assert result.replicate == 1
        assert result.selection_string == "protein and name CA"


class TestComputeReplicateNegativePaths:
    """Test RMSFAnalysis.compute_replicate error and failure paths."""

    @staticmethod
    def _make_nonempty_mock_universe(n_frames: int = 100, n_atoms: int = 5) -> MagicMock:
        """Create a mock Universe with a non-empty selection."""
        mock_u = MagicMock()
        mock_u.trajectory.__len__ = MagicMock(return_value=n_frames)

        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=n_atoms)
        mock_u.select_atoms = MagicMock(return_value=mock_atoms)
        return mock_u

    @patch("polyzymd.analyses.rmsf.TrajectoryLoader")
    def test_raises_when_selection_is_empty(
        self,
        MockLoader,
        rmsf_analysis,
        condition,
        tmp_path,
    ):
        """An empty atom selection should raise ValueError with diagnostics."""
        mock_loader_inst = MagicMock()
        MockLoader.return_value = mock_loader_inst

        mock_u = MagicMock()
        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=0)
        mock_u.select_atoms = MagicMock(return_value=mock_atoms)
        mock_loader_inst.load_universe.return_value = mock_u

        traj_info = MagicMock()
        traj_info.trajectory_files = [Path("/fake/traj.dcd")]
        mock_loader_inst.get_trajectory_info.return_value = traj_info

        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=True,
            settings=RMSFSettings(selection="protein and name ZZ"),
        )

        with patch(
            "polyzymd.analyses.rmsf.get_selection_diagnostics",
            return_value="Selection diagnostics",
        ):
            with pytest.raises(ValueError, match="matched no atoms"):
                rmsf_analysis.compute_replicate(ctx, 1)

    @patch("polyzymd.analyses.rmsf.validate_equilibration_time")
    @patch("polyzymd.analyses.rmsf.TrajectoryLoader")
    def test_raises_when_equilibration_invalid(
        self,
        MockLoader,
        mock_validate_equilibration,
        rmsf_analysis,
        condition,
        tmp_path,
    ):
        """Invalid equilibration time should raise the validation message."""
        mock_validate_equilibration.return_value = (
            False,
            "Equilibration time exceeds trajectory length",
        )

        mock_loader_inst = MagicMock()
        MockLoader.return_value = mock_loader_inst
        mock_loader_inst.load_universe.return_value = self._make_nonempty_mock_universe(n_frames=20)

        traj_info = MagicMock()
        traj_info.trajectory_files = [Path("/fake/traj.dcd")]
        mock_loader_inst.get_trajectory_info.return_value = traj_info
        mock_loader_inst.get_timestep.return_value = 10.0

        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=True,
            settings=RMSFSettings(),
        )

        with pytest.raises(ValueError, match="Equilibration time exceeds trajectory length"):
            rmsf_analysis.compute_replicate(ctx, 1)

    @patch("polyzymd.analyses.rmsf.TrajectoryLoader")
    def test_propagates_file_not_found_from_loader(
        self,
        MockLoader,
        rmsf_analysis,
        condition,
        tmp_path,
    ):
        """FileNotFoundError from trajectory loading should propagate."""
        mock_loader_inst = MagicMock()
        MockLoader.return_value = mock_loader_inst
        mock_loader_inst.load_universe.side_effect = FileNotFoundError("Missing trajectory file")

        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=True,
            settings=RMSFSettings(),
        )

        with pytest.raises(FileNotFoundError, match="Missing trajectory file"):
            rmsf_analysis.compute_replicate(ctx, 1)

    def test_raises_when_frame_mode_missing_reference_frame(
        self,
        rmsf_analysis,
        condition,
        tmp_path,
    ):
        """Frame reference mode requires an explicit reference_frame."""
        settings = RMSFSettings(reference_mode="frame", reference_frame=None)
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=True,
            settings=settings,
        )

        with pytest.raises(ValueError, match="reference_frame is required"):
            rmsf_analysis.compute_replicate(ctx, 1)

    def test_raises_when_external_reference_file_missing(
        self,
        rmsf_analysis,
        condition,
        tmp_path,
    ):
        """External reference mode should reject nonexistent reference files."""
        settings = RMSFSettings(
            reference_mode="external",
            reference_file=str(tmp_path / "does_not_exist.pdb"),
        )
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=True,
            settings=settings,
        )

        with pytest.raises(ValueError, match="reference_file does not exist"):
            rmsf_analysis.compute_replicate(ctx, 1)


# ============================================================================
# Test: aggregate
# ============================================================================


class TestAggregate:
    """Test RMSFAnalysis.aggregate."""

    def test_aggregates_results(self, rmsf_analysis, condition, tmp_path):
        """Test aggregate produces an RMSFAggregatedResult with correct stats."""
        results = [_make_mock_rmsf_result(i, mean_rmsf=1.0 + 0.1 * i) for i in range(1, 4)]
        agg_dir = tmp_path / "aggregated"

        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2, 3),
            output_dir=agg_dir,
            equilibration="10ns",
            settings=RMSFSettings(),
        )

        # Patch only get_polyzymd_version (non-essential) and let real functions run
        with patch("polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.2.1"):
            result = rmsf_analysis.aggregate(ctx, results)

        # Verify the result has correct attributes
        assert result.n_replicates == 3
        assert result.replicates == [1, 2, 3]
        assert len(result.mean_rmsf_per_residue) == 5
        assert len(result.sem_rmsf_per_residue) == 5
        assert len(result.per_replicate_mean_rmsf) == 3
        assert result.overall_mean_rmsf == pytest.approx(1.2, abs=0.01)  # mean of 1.1, 1.2, 1.3

    def test_aggregate_saves_result_file(self, rmsf_analysis, condition, tmp_path):
        """Verify aggregated result is saved to the output directory."""
        results = [_make_mock_rmsf_result(i, mean_rmsf=1.5) for i in range(1, 3)]
        agg_dir = tmp_path / "aggregated"

        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2),
            output_dir=agg_dir,
            equilibration="10ns",
            settings=RMSFSettings(),
        )

        with patch("polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.2.1"):
            rmsf_analysis.aggregate(ctx, results)

        # Check that a file was saved
        json_files = list(agg_dir.glob("*.json"))
        assert len(json_files) == 1
        assert "rmsf_reps1-2_eq10ns.json" in json_files[0].name

    def test_aggregate_uses_correct_replicate_means(self, rmsf_analysis, condition, tmp_path):
        """Verify that per_replicate_mean_rmsf contains the correct values."""
        results = [_make_mock_rmsf_result(i, mean_rmsf=1.0 + 0.1 * i) for i in range(1, 4)]

        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2, 3),
            output_dir=tmp_path / "aggregated",
            equilibration="10ns",
            settings=RMSFSettings(),
        )

        with patch("polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.2.1"):
            result = rmsf_analysis.aggregate(ctx, results)

        assert result.per_replicate_mean_rmsf == [1.1, 1.2, 1.3]


# ============================================================================
# Test: extract_metrics
# ============================================================================


class TestExtractMetrics:
    """Test RMSFAnalysis.extract_metrics for default comparison pipeline."""

    def test_returns_mean_rmsf_metric(self, rmsf_analysis):
        summary = MagicMock()
        summary.overall_mean_rmsf = 1.45
        summary.overall_sem_rmsf = 0.08
        summary.per_replicate_mean_rmsf = [1.4, 1.5, 1.45]

        metrics = rmsf_analysis.extract_metrics(summary)

        assert "mean_rmsf" in metrics
        m = metrics["mean_rmsf"]
        assert isinstance(m, MetricValue)
        assert m.name == "mean_rmsf"
        assert m.mean == 1.45
        assert m.sem == 0.08
        assert m.replicate_values == [1.4, 1.5, 1.45]
        assert m.higher_is_better is False
        assert m.direction_labels == ("stabilizing", "unchanged", "destabilizing")

    def test_returns_single_metric(self, rmsf_analysis):
        """RMSF should return exactly one metric."""
        summary = MagicMock()
        summary.overall_mean_rmsf = 1.0
        summary.overall_sem_rmsf = 0.05
        summary.per_replicate_mean_rmsf = [1.0, 1.0]

        metrics = rmsf_analysis.extract_metrics(summary)
        assert len(metrics) == 1


# ============================================================================
# Test: AggregatedResultClass and _deserialize_result
# ============================================================================


class TestDeserializeResult:
    """Test RMSFAnalysis.AggregatedResultClass and _deserialize_result."""

    def test_aggregated_result_class_set(self, rmsf_analysis):
        """AggregatedResultClass should be RMSFAggregatedResult."""
        from polyzymd.analyses.rmsf._results import RMSFAggregatedResult

        assert rmsf_analysis.AggregatedResultClass is RMSFAggregatedResult

    def test_loads_aggregated_result(self, rmsf_analysis, tmp_path):
        """Test that _deserialize_result loads via RMSFAggregatedResult.load()."""
        mock_loaded = MagicMock()
        with patch.object(
            rmsf_analysis.AggregatedResultClass, "load", return_value=mock_loaded
        ) as mock_load:
            result = rmsf_analysis._deserialize_result(tmp_path / "test.json")

            mock_load.assert_called_once_with(tmp_path / "test.json")
            assert result is mock_loaded


# ============================================================================
# Test: filter_conditions (default behavior)
# ============================================================================


class TestFilterConditions:
    """RMSF applies to all conditions — filter should keep all."""

    def test_keeps_all_conditions(self, rmsf_analysis):
        conditions = [
            Condition(
                label=f"Cond {i}",
                config_path=Path(f"/fake/config_{i}.yaml"),
                replicates=(1, 2),
                sim_config=MagicMock(),
            )
            for i in range(4)
        ]

        filtered = rmsf_analysis.filter_conditions(conditions)
        assert len(filtered) == 4


# ============================================================================
# Test: _make_aggregated_filename
# ============================================================================


class TestMakeAggregatedFilename:
    """Test backward-compatible filename generation."""

    def test_contiguous_replicates(self):
        result = MagicMock()
        result.equilibration_time = 10.0
        result.equilibration_unit = "ns"

        filename = RMSFAnalysis._make_aggregated_filename((1, 2, 3, 4, 5), result)
        assert filename == "rmsf_reps1-5_eq10ns.json"

    def test_non_contiguous_replicates(self):
        result = MagicMock()
        result.equilibration_time = 100.0
        result.equilibration_unit = "ns"

        filename = RMSFAnalysis._make_aggregated_filename((1, 3, 5), result)
        assert filename == "rmsf_reps1_3_5_eq100ns.json"

    def test_single_pair(self):
        result = MagicMock()
        result.equilibration_time = 0.0
        result.equilibration_unit = "ns"

        filename = RMSFAnalysis._make_aggregated_filename((1, 2), result)
        assert filename == "rmsf_reps1-2_eq0ns.json"


# ============================================================================
# Test: plot
# ============================================================================


class TestPlot:
    """Test RMSFAnalysis.plot calls inlined private plotting functions."""

    def test_plot_returns_empty_on_no_conditions(self, rmsf_analysis, tmp_path):
        from polyzymd.config.comparison import PlotSettings

        ctx = PlotContext(
            conditions=[],
            analysis_dirs={},
            results_dir=tmp_path / "comparison",
            output_dir=tmp_path / "figures",
            settings=RMSFSettings(),
            plot_settings=PlotSettings(),
        )

        plots = rmsf_analysis.plot(ctx)
        assert plots == []

    @patch("polyzymd.analyses.rmsf._plot_rmsf_profile")
    @patch("polyzymd.analyses.rmsf._plot_rmsf_comparison")
    def test_plot_delegates_to_private_functions(
        self,
        mock_comp_plot,
        mock_prof_plot,
        rmsf_analysis,
        condition,
        tmp_path,
    ):
        from polyzymd.config.comparison import PlotSettings

        # Setup mock returns
        mock_comp_plot.return_value = [tmp_path / "figures" / "rmsf_comparison.png"]
        mock_prof_plot.return_value = [tmp_path / "figures" / "rmsf_profile.png"]

        analysis_dir = tmp_path / "analysis" / "no_polymer" / "rmsf"
        analysis_dir.mkdir(parents=True)

        ctx = PlotContext(
            conditions=[condition],
            analysis_dirs={"No Polymer": analysis_dir},
            results_dir=tmp_path / "comparison",
            output_dir=tmp_path / "figures",
            settings=RMSFSettings(),
            plot_settings=PlotSettings(),
        )

        with patch("polyzymd.config.comparison.PlotSettings") as MockPlotSettings:
            MockPlotSettings.return_value = MagicMock()
            plots = rmsf_analysis.plot(ctx)

        assert len(plots) == 2
        mock_comp_plot.assert_called_once()
        mock_prof_plot.assert_called_once()

    @patch("polyzymd.analyses.rmsf._plot_rmsf_profile", side_effect=Exception("plot error"))
    @patch("polyzymd.analyses.rmsf._plot_rmsf_comparison", side_effect=Exception("plot error"))
    def test_plot_propagates_exceptions(
        self,
        mock_comp_plot,
        mock_prof_plot,
        rmsf_analysis,
        condition,
        tmp_path,
    ):
        """Plotting failures should propagate to orchestrator."""
        from polyzymd.config.comparison import PlotSettings

        analysis_dir = tmp_path / "analysis" / "no_polymer" / "rmsf"
        analysis_dir.mkdir(parents=True)

        ctx = PlotContext(
            conditions=[condition],
            analysis_dirs={"No Polymer": analysis_dir},
            results_dir=tmp_path / "comparison",
            output_dir=tmp_path / "figures",
            settings=RMSFSettings(),
            plot_settings=PlotSettings(),
        )

        with patch("polyzymd.config.comparison.PlotSettings") as MockPlotSettings:
            MockPlotSettings.return_value = MagicMock()
            with pytest.raises(Exception, match="plot error"):
                rmsf_analysis.plot(ctx)


# ============================================================================
# Test: Full lifecycle via orchestrator (mocked compute)
# ============================================================================


class TestRMSFLifecycle:
    """Test the full compute -> aggregate -> compare lifecycle.

    Uses the orchestrator's run_analysis() and verifies the data flow
    through the RMSF plugin. Heavy deps are mocked.
    """

    @patch("polyzymd.analyses.rmsf.align_trajectory", return_value=0)
    @patch("polyzymd.analyses.rmsf.validate_equilibration_time", return_value=(True, ""))
    @patch("polyzymd.analyses.shared.config_hash.validate_config_hash")
    @patch("polyzymd.analyses.rmsf.compute_config_hash", return_value="abc123")
    @patch("polyzymd.analyses.rmsf.TrajectoryLoader")
    def test_run_analysis_lifecycle(
        self,
        MockLoader,
        mock_hash,
        mock_validate_hash,
        mock_eq_validate,
        mock_align,
        rmsf_analysis,
        condition,
        tmp_path,
    ):
        """Test compute_replicate -> aggregate via run_analysis()."""
        import numpy as np

        from polyzymd.analyses.orchestrator import run_analysis

        # Create a fresh mock loader instance for each call
        mock_loader_inst = MagicMock()
        MockLoader.return_value = mock_loader_inst

        # Track which replicate we're on so we can vary mean_rmsf
        call_count = {"n": 0}
        rmsf_by_rep = {
            1: [1.0, 1.2, 1.4, 1.0, 0.9],
            2: [1.1, 1.3, 1.5, 1.1, 1.0],
            3: [1.2, 1.4, 1.6, 1.2, 1.1],
        }

        def make_mock_universe(*args, **kwargs):
            """Return a mock Universe with 200 frames and 5 atoms."""
            mock_u = MagicMock()
            mock_traj = MagicMock()
            mock_traj.__len__ = MagicMock(return_value=200)
            mock_frames = [MagicMock(frame=i) for i in range(200)]
            mock_traj.__getitem__ = MagicMock(
                side_effect=lambda x: mock_frames[x] if isinstance(x, int) else mock_frames
            )
            mock_traj.__iter__ = MagicMock(return_value=iter(mock_frames))
            mock_u.trajectory = mock_traj

            mock_atoms = MagicMock()
            mock_atoms.__len__ = MagicMock(return_value=5)
            mock_atoms.positions = np.random.rand(5, 3).astype(np.float32)
            mock_atoms.indices = np.arange(5)
            mock_residues = []
            for i in range(5):
                res = MagicMock()
                res.resid = i + 1
                res.resname = "ALA"
                mock_residues.append(res)
            mock_atoms.residues = mock_residues
            mock_u.select_atoms = MagicMock(return_value=mock_atoms)
            return mock_u

        mock_loader_inst.load_universe = MagicMock(side_effect=make_mock_universe)

        traj_info = MagicMock()
        traj_info.trajectory_files = [Path("/fake/traj.dcd")]
        mock_loader_inst.get_trajectory_info.return_value = traj_info
        mock_loader_inst.get_timestep.return_value = 10.0

        # _compute_rmsf returns different values per call to vary per-replicate means
        rmsf_values_sequence = [
            np.array(rmsf_by_rep[1]),
            np.array(rmsf_by_rep[2]),
            np.array(rmsf_by_rep[3]),
        ]
        rmsf_call_idx = {"n": 0}

        def mock_compute_rmsf(*args, **kwargs):
            idx = rmsf_call_idx["n"]
            rmsf_call_idx["n"] += 1
            return rmsf_values_sequence[idx]

        with (
            patch("polyzymd.analyses.rmsf._compute_rmsf", side_effect=mock_compute_rmsf),
            patch(
                "polyzymd.analyses.rmsf._compute_rmsd_timeseries",
                return_value=np.random.rand(100),
            ),
            patch("polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.2.1"),
        ):
            output_dir = tmp_path / "analysis" / "no_polymer" / "rmsf"
            result = run_analysis(
                rmsf_analysis,
                condition,
                settings=RMSFSettings(),
                equilibration="10ns",
                output_dir=output_dir,
            )

        # Verify we got an aggregated result
        assert result.n_replicates == 3
        # Per-replicate means should match np.mean of each rmsf array
        expected_means = [float(np.mean(rmsf_by_rep[r])) for r in (1, 2, 3)]
        for actual, expected in zip(result.per_replicate_mean_rmsf, expected_means):
            assert actual == pytest.approx(expected, abs=0.01)
        # Loader should have been called 3 times (once per replicate)
        assert MockLoader.call_count == 3

    def test_extract_metrics_feeds_default_compare(self, rmsf_analysis):
        """Verify extract_metrics output is compatible with default_scalar_comparison."""
        from polyzymd.analyses.stats import default_scalar_comparison

        # Simulate two conditions with extracted metrics
        metrics_by_condition = {}
        for label, mean, sem, vals in [
            ("No Polymer", 1.8, 0.12, [1.7, 1.9, 1.8]),
            ("100% SBMA", 1.3, 0.08, [1.25, 1.35, 1.3]),
        ]:
            summary = MagicMock()
            summary.overall_mean_rmsf = mean
            summary.overall_sem_rmsf = sem
            summary.per_replicate_mean_rmsf = vals
            metrics_by_condition[label] = rmsf_analysis.extract_metrics(summary)

        # Run default comparison — should not crash
        result = default_scalar_comparison(
            analysis_name="rmsf",
            project_name="test_project",
            metrics_by_condition=metrics_by_condition,
            control_label="No Polymer",
            equilibration="10ns",
        )

        assert result.analysis_type == "rmsf"
        assert len(result.conditions) == 2
        assert len(result.pairwise_comparisons) >= 1
        # Ranking: lower RMSF first (100% SBMA should rank before No Polymer)
        assert result.ranking[0] == "100% SBMA"


class TestAggregatePerResidue:
    """Tests for _aggregate_per_residue cross-chain correctness."""

    def test_duplicate_resid_across_chains_not_merged(self):
        """Residues with same resid in different chains must stay separate."""
        from polyzymd.analyses.rmsf import _aggregate_per_residue

        mock_atoms = MagicMock()
        mock_atoms.indices = np.array([0, 1, 2, 3])

        res_a = MagicMock()
        res_a.resid = 1
        res_a.atoms.indices = np.array([0, 1])

        res_b = MagicMock()
        res_b.resid = 1
        res_b.atoms.indices = np.array([2, 3])

        mock_atoms.residues = [res_a, res_b]

        atom_rmsf = np.array([1.0, 2.0, 10.0, 20.0])

        result = _aggregate_per_residue(mock_atoms, atom_rmsf)
        assert len(result) == 2
        np.testing.assert_allclose(result[0], 1.5)
        np.testing.assert_allclose(result[1], 15.0)
