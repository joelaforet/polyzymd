"""Tests for the distances analysis plugin.

Covers discovery, settings, compute_replicate, aggregate, compare (full override),
plot delegation, AggregatedResultClass, _make_aggregated_filename, and lifecycle.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


class TestDiscovery:
    """The plugin is found by the automatic discovery system."""

    def test_discovery_by_name(self):
        from polyzymd.analyses.discovery import get_analysis

        cls = get_analysis("distances")
        assert cls is not None
        assert cls.name == "distances"

    def test_listed(self):
        from polyzymd.analyses.discovery import list_analyses

        analyses = list_analyses()
        names = list(analyses.keys())
        assert "distances" in names

    def test_all_names(self):
        from polyzymd.analyses.discovery import list_all_names

        names = list_all_names()
        assert "distances" in names


# ---------------------------------------------------------------------------
# Class Attributes
# ---------------------------------------------------------------------------


class TestClassAttributes:
    """Verify ClassVar declarations on DistancesAnalysis."""

    def test_name(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        assert DistancesAnalysis.name == "distances"

    def test_settings_type(self):
        from polyzymd.analyses.distances import DistancesAnalysis, DistancesSettings

        assert DistancesAnalysis.Settings is DistancesSettings

    def test_aliases(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        assert DistancesAnalysis.aliases == ()

    def test_dependencies(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        assert DistancesAnalysis.dependencies == ()

    def test_min_replicates(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        assert DistancesAnalysis.min_replicates == 2


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class TestSettings:
    """Validate DistancesSettings and DistancePairSettings."""

    def test_minimal_valid(self):
        from polyzymd.analyses.distances import DistancePairSettings, DistancesSettings

        s = DistancesSettings(
            pairs=[
                DistancePairSettings(
                    label="Pair1",
                    selection_a="resid 10 and name CA",
                    selection_b="resid 20 and name CA",
                )
            ]
        )
        assert s.threshold == 3.5
        assert s.use_pbc is True
        assert s.align_trajectory is True
        assert len(s.pairs) == 1

    def test_empty_pairs_rejected(self):
        from polyzymd.analyses.distances import DistancesSettings

        with pytest.raises(ValueError):
            DistancesSettings(pairs=[])

    def test_get_pair_selections(self):
        from polyzymd.analyses.distances import DistancePairSettings, DistancesSettings

        s = DistancesSettings(
            pairs=[
                DistancePairSettings(label="A", selection_a="sel1", selection_b="sel2"),
                DistancePairSettings(label="B", selection_a="sel3", selection_b="sel4"),
            ]
        )
        assert s.get_pair_selections() == [("sel1", "sel2"), ("sel3", "sel4")]

    def test_get_pair_labels(self):
        from polyzymd.analyses.distances import DistancePairSettings, DistancesSettings

        s = DistancesSettings(
            pairs=[
                DistancePairSettings(label="Alpha", selection_a="a", selection_b="b"),
                DistancePairSettings(label="Beta", selection_a="c", selection_b="d"),
            ]
        )
        assert s.get_pair_labels() == ["Alpha", "Beta"]

    def test_get_pair_thresholds_fallback(self):
        from polyzymd.analyses.distances import DistancePairSettings, DistancesSettings

        s = DistancesSettings(
            threshold=5.0,
            pairs=[
                DistancePairSettings(label="A", selection_a="a", selection_b="b", threshold=3.0),
                DistancePairSettings(
                    label="B", selection_a="c", selection_b="d"
                ),  # No per-pair threshold
            ],
        )
        thresholds = s.get_pair_thresholds()
        assert thresholds == [3.0, 5.0]  # per-pair, then fallback

    def test_alignment_mode_validation(self):
        from polyzymd.analyses.distances import DistancePairSettings, DistancesSettings

        with pytest.raises(ValueError):
            DistancesSettings(
                pairs=[DistancePairSettings(label="A", selection_a="a", selection_b="b")],
                alignment_mode="invalid",
            )

    def test_alignment_frame_required(self):
        from polyzymd.analyses.distances import DistancePairSettings, DistancesSettings

        with pytest.raises(ValueError):
            DistancesSettings(
                pairs=[DistancePairSettings(label="A", selection_a="a", selection_b="b")],
                alignment_mode="frame",
                # alignment_frame not provided
            )

    def test_alignment_config(self):
        from polyzymd.analyses.distances import DistancePairSettings, DistancesSettings

        s = DistancesSettings(
            pairs=[DistancePairSettings(label="A", selection_a="a", selection_b="b")],
            align_trajectory=True,
            alignment_selection="backbone",
            alignment_mode="centroid",
        )
        cfg = s.get_alignment_config()
        assert cfg.enabled is True
        assert cfg.selection == "backbone"

    def test_serialization_roundtrip(self):
        from polyzymd.analyses.distances import DistancePairSettings, DistancesSettings

        s = DistancesSettings(
            threshold=4.0,
            pairs=[
                DistancePairSettings(
                    label="X",
                    selection_a="resid 1",
                    selection_b="resid 2",
                    threshold=3.0,
                    below_label="Bound",
                    above_label="Unbound",
                )
            ],
            use_pbc=False,
        )
        d = s.model_dump()
        s2 = DistancesSettings.model_validate(d)
        assert s2.threshold == 4.0
        assert s2.pairs[0].below_label == "Bound"
        assert s2.use_pbc is False


# ---------------------------------------------------------------------------
# compute_replicate
# ---------------------------------------------------------------------------


def _make_mock_distance_result(replicate: int = 1, n_pairs: int = 2):
    """Create a mock DistanceResult."""
    mock = MagicMock()
    mock.replicate = replicate
    mock.config_hash = "abc123"
    mock.equilibration_time = 100.0
    mock.equilibration_unit = "ns"

    pair_results = []
    for i in range(n_pairs):
        pr = MagicMock()
        pr.pair_label = f"pair{i}"
        pr.selection1 = f"sel_a_{i}"
        pr.selection2 = f"sel_b_{i}"
        pr.mean_distance = 3.5 + i * 0.5
        pr.std_distance = 0.5
        pr.median_distance = 3.4 + i * 0.5
        pr.min_distance = 2.0 + i * 0.3
        pr.max_distance = 6.0 + i * 0.3
        pr.sem_distance = 0.1
        pr.threshold = 3.5
        pr.fraction_below_threshold = 0.6 - i * 0.1
        pr.kde_peak = 3.3 + i * 0.5
        pr.per_replicate_means = [3.5 + i * 0.5]
        pair_results.append(pr)

    mock.pair_results = pair_results
    return mock


class TestComputeReplicate:
    """compute_replicate delegates to DistanceCalculator."""

    def test_delegates_to_calculator(self):
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )

        analysis = DistancesAnalysis()

        settings = DistancesSettings(
            pairs=[DistancePairSettings(label="P1", selection_a="sel_a", selection_b="sel_b")]
        )

        mock_sim_config = MagicMock()
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2, 3),
            sim_config=mock_sim_config,
        )
        ctx = ReplicateContext(
            condition=cond,
            replicate=1,
            sim_config=mock_sim_config,
            output_dir=Path("/tmp/out/run_1"),
            equilibration="10ns",
            recompute=False,
            settings=settings,
        )

        mock_result = _make_mock_distance_result(replicate=1, n_pairs=1)

        with patch("polyzymd.analyses.distances.DistanceCalculator") as MockCalc:
            MockCalc.return_value.compute.return_value = mock_result
            result = analysis.compute_replicate(ctx, 1)

        assert result is mock_result
        MockCalc.assert_called_once()
        MockCalc.return_value.compute.assert_called_once_with(
            replicate=1,
            save=True,
            output_dir=Path("/tmp/out/run_1"),
            recompute=False,
        )

    def test_passes_settings_to_calculator(self):
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )

        analysis = DistancesAnalysis()

        settings = DistancesSettings(
            threshold=5.0,
            pairs=[
                DistancePairSettings(
                    label="P1",
                    selection_a="resid 10",
                    selection_b="resid 20",
                    threshold=4.0,
                )
            ],
            use_pbc=False,
        )

        mock_sim_config = MagicMock()
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
            output_dir=Path("/tmp/out/run_1"),
            equilibration="10ns",
            recompute=False,
            settings=settings,
        )

        mock_result = _make_mock_distance_result(replicate=1, n_pairs=1)

        with patch("polyzymd.analyses.distances.DistanceCalculator") as MockCalc:
            MockCalc.return_value.compute.return_value = mock_result
            analysis.compute_replicate(ctx, 1)

        # Verify settings were passed correctly
        call_kwargs = MockCalc.call_args[1]
        assert call_kwargs["pairs"] == [("resid 10", "resid 20")]
        assert call_kwargs["thresholds"] == [4.0]  # per-pair threshold, not global
        assert call_kwargs["use_pbc"] is False


# ---------------------------------------------------------------------------
# aggregate
# ---------------------------------------------------------------------------


class TestAggregate:
    """aggregate produces a DistanceAggregatedResult."""

    def _make_mock_results(self, n_reps: int = 3, n_pairs: int = 2):
        """Create mock per-replicate results for aggregation."""
        results = []
        for rep in range(1, n_reps + 1):
            mock = MagicMock()
            mock.replicate = rep
            mock.config_hash = "hash123"
            mock.equilibration_time = 100.0
            mock.equilibration_unit = "ns"

            pair_results = []
            for pair_idx in range(n_pairs):
                pr = MagicMock()
                pr.pair_label = f"pair{pair_idx}"
                pr.selection1 = f"sel_a_{pair_idx}"
                pr.selection2 = f"sel_b_{pair_idx}"
                pr.mean_distance = 3.5 + pair_idx * 0.5 + rep * 0.01
                pr.std_distance = 0.5
                pr.median_distance = 3.4 + pair_idx * 0.5 + rep * 0.01
                pr.threshold = 3.5
                pr.fraction_below_threshold = 0.6 - pair_idx * 0.1
                pr.kde_peak = 3.3 + pair_idx * 0.5 + rep * 0.01
                pair_results.append(pr)

            mock.pair_results = pair_results
            results.append(mock)
        return results

    def test_aggregate_produces_result(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )

        analysis = DistancesAnalysis()
        settings = DistancesSettings(
            pairs=[
                DistancePairSettings(label="P0", selection_a="sel_a_0", selection_b="sel_b_0"),
                DistancePairSettings(label="P1", selection_a="sel_a_1", selection_b="sel_b_1"),
            ]
        )

        output_dir = tmp_path / "aggregated"
        mock_sim_config = MagicMock()
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
            equilibration="100ns",
            settings=settings,
        )

        mock_results = self._make_mock_results(n_reps=3, n_pairs=2)

        with patch(
            "polyzymd.analyses.shared.aggregation.aggregate_distance_pair_stats"
        ) as mock_agg:
            # Create mock return for aggregate_distance_pair_stats
            mock_stats = MagicMock()
            mock_stats.mean_stats.mean = 3.5
            mock_stats.mean_stats.sem = 0.05
            mock_stats.median_stats.mean = 3.4
            mock_stats.per_rep_means = [3.49, 3.50, 3.51]
            mock_stats.per_rep_stds = [0.5, 0.5, 0.5]
            mock_stats.per_rep_medians = [3.39, 3.40, 3.41]
            mock_stats.fraction_stats = MagicMock()
            mock_stats.fraction_stats.mean = 0.6
            mock_stats.fraction_stats.sem = 0.02
            mock_stats.per_rep_fractions = [0.58, 0.60, 0.62]
            mock_stats.kde_peak_stats = MagicMock()
            mock_stats.kde_peak_stats.mean = 3.3
            mock_stats.kde_peak_stats.sem = 0.03
            mock_stats.per_rep_kde_peaks = [3.29, 3.30, 3.31]
            mock_agg.return_value = mock_stats

            with patch(
                "polyzymd.analyses._results_base.get_polyzymd_version",
                return_value="1.0.0-test",
            ):
                result = analysis.aggregate(ctx, mock_results)

        assert result is not None
        assert result.n_replicates == 3
        assert len(result.pair_results) == 2

    def test_aggregate_saves_file(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )

        analysis = DistancesAnalysis()
        settings = DistancesSettings(
            pairs=[
                DistancePairSettings(label="P0", selection_a="a", selection_b="b"),
            ]
        )

        output_dir = tmp_path / "aggregated"
        mock_sim_config = MagicMock()
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=mock_sim_config,
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1, 2),
            output_dir=output_dir,
            equilibration="100ns",
            settings=settings,
        )

        mock_results = self._make_mock_results(n_reps=2, n_pairs=1)

        with patch(
            "polyzymd.analyses.shared.aggregation.aggregate_distance_pair_stats"
        ) as mock_agg:
            mock_stats = MagicMock()
            mock_stats.mean_stats.mean = 3.5
            mock_stats.mean_stats.sem = 0.05
            mock_stats.median_stats.mean = 3.4
            mock_stats.per_rep_means = [3.49, 3.51]
            mock_stats.per_rep_stds = [0.5, 0.5]
            mock_stats.per_rep_medians = [3.39, 3.41]
            mock_stats.fraction_stats = None
            mock_stats.per_rep_fractions = None
            mock_stats.kde_peak_stats = None
            mock_stats.per_rep_kde_peaks = None
            mock_agg.return_value = mock_stats

            with patch(
                "polyzymd.analyses._results_base.get_polyzymd_version",
                return_value="1.0.0-test",
            ):
                analysis.aggregate(ctx, mock_results)

        # Check that output file was created
        json_files = list(output_dir.glob("*.json"))
        assert len(json_files) == 1


# ---------------------------------------------------------------------------
# compare (full override)
# ---------------------------------------------------------------------------


def _make_mock_agg_result(n_pairs: int = 2, n_reps: int = 3, offset: float = 0.0):
    """Create a mock DistanceAggregatedResult for comparison tests."""
    mock = MagicMock()
    mock.n_replicates = n_reps

    pair_results = []
    for i in range(n_pairs):
        pr = MagicMock()
        pr.pair_label = f"pair{i}"
        pr.selection1 = f"sel_a_{i}"
        pr.selection2 = f"sel_b_{i}"
        pr.overall_mean = 3.5 + i * 0.5 + offset
        pr.overall_sem = 0.05
        pr.threshold = 3.5
        pr.overall_fraction_below = 0.6 - i * 0.1 - offset * 0.1
        pr.sem_fraction_below = 0.02
        pr.per_replicate_means = [3.5 + i * 0.5 + offset + j * 0.01 for j in range(n_reps)]
        pr.per_replicate_fractions_below = [
            0.6 - i * 0.1 - offset * 0.1 + j * 0.005 for j in range(n_reps)
        ]
        pair_results.append(pr)

    mock.pair_results = pair_results
    return mock


class TestCompare:
    """compare() produces a DistanceComparisonResult with per-pair stats."""

    def _make_analysis_and_ctx(self, tmp_path, n_conditions: int = 3, control: str = "Control"):
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )

        analysis = DistancesAnalysis()

        settings = DistancesSettings(
            pairs=[
                DistancePairSettings(
                    label="Alpha",
                    selection_a="sel_a_0",
                    selection_b="sel_b_0",
                ),
                DistancePairSettings(
                    label="Beta",
                    selection_a="sel_a_1",
                    selection_b="sel_b_1",
                ),
            ]
        )

        conditions = []
        analysis_dirs = {}
        labels = ["Control", "Treatment_A", "Treatment_B"][:n_conditions]

        for label in labels:
            cond = Condition(
                label=label,
                config_path=Path(f"/tmp/{label}/config.yaml"),
                replicates=(1, 2, 3),
                sim_config=MagicMock(),
            )
            conditions.append(cond)
            agg_dir = tmp_path / label / "distances"
            agg_dir.mkdir(parents=True, exist_ok=True)
            analysis_dirs[label] = agg_dir

        ctx = ComparisonContext(
            name="test_comparison",
            conditions=conditions,
            excluded_conditions=[],
            control_label=control,
            analysis_dirs=analysis_dirs,
            results_dir=tmp_path / "results",
            equilibration="100ns",
            settings=settings,
            recompute=False,
        )

        return analysis, ctx, labels

    def test_compare_returns_result(self, tmp_path):
        analysis, ctx, labels = self._make_analysis_and_ctx(tmp_path)

        # Mock _load_aggregated_result to return different values per condition
        # Path structure: tmp_path / label / "distances" / "aggregated"
        # So agg_dir.parent.parent.name gives the condition label
        offsets = {"Control": 0.0, "Treatment_A": 0.5, "Treatment_B": -0.3}

        def mock_load(agg_dir):
            label = agg_dir.parent.parent.name
            return _make_mock_agg_result(n_pairs=2, offset=offsets.get(label, 0.0))

        with patch.object(analysis, "_load_aggregated_result", side_effect=mock_load):
            result = analysis.compare(ctx)

        assert result is not None
        assert result.metric == "mean_distance"
        assert result.name == "test_comparison"
        assert result.n_pairs == 2
        assert len(result.pair_labels) == 2
        assert result.pair_labels == ["Alpha", "Beta"]
        assert len(result.conditions) == 3
        assert result.control_label == "Control"

    def test_compare_has_per_pair_rankings(self, tmp_path):
        analysis, ctx, labels = self._make_analysis_and_ctx(tmp_path)

        offsets = {"Control": 0.0, "Treatment_A": 0.5, "Treatment_B": -0.3}

        def mock_load(agg_dir):
            label = agg_dir.parent.parent.name
            return _make_mock_agg_result(n_pairs=2, offset=offsets.get(label, 0.0))

        with patch.object(analysis, "_load_aggregated_result", side_effect=mock_load):
            result = analysis.compare(ctx)

        assert "Alpha" in result.ranking_by_pair
        assert "Beta" in result.ranking_by_pair
        # Treatment_B has lowest offset (-0.3) so should rank first for each pair
        assert result.ranking_by_pair["Alpha"][0] == "Treatment_B"

    def test_compare_pairwise_comparisons(self, tmp_path):
        analysis, ctx, labels = self._make_analysis_and_ctx(tmp_path)

        offsets = {"Control": 0.0, "Treatment_A": 0.5, "Treatment_B": -0.3}

        def mock_load(agg_dir):
            label = agg_dir.parent.parent.name
            return _make_mock_agg_result(n_pairs=2, offset=offsets.get(label, 0.0))

        with patch.object(analysis, "_load_aggregated_result", side_effect=mock_load):
            result = analysis.compare(ctx)

        # With control, should have (n_conditions - 1) * n_pairs comparisons
        assert len(result.pairwise_comparisons) == 4  # 2 treatments * 2 pairs

        # Each comparison should have distance and fraction fields
        for comp in result.pairwise_comparisons:
            assert comp.condition_a == "Control"
            assert hasattr(comp, "distance_t_statistic")
            assert hasattr(comp, "distance_direction")

    def test_compare_anova_with_3_conditions(self, tmp_path):
        analysis, ctx, labels = self._make_analysis_and_ctx(tmp_path, n_conditions=3)

        offsets = {"Control": 0.0, "Treatment_A": 0.5, "Treatment_B": -0.3}

        def mock_load(agg_dir):
            label = agg_dir.parent.parent.name
            return _make_mock_agg_result(n_pairs=2, offset=offsets.get(label, 0.0))

        with patch.object(analysis, "_load_aggregated_result", side_effect=mock_load):
            result = analysis.compare(ctx)

        # With 3+ conditions, ANOVA should be computed
        assert result.anova_by_pair is not None
        assert len(result.anova_by_pair) == 2  # one per pair

    def test_compare_no_anova_with_2_conditions(self, tmp_path):
        analysis, ctx, labels = self._make_analysis_and_ctx(tmp_path, n_conditions=2)

        offsets = {"Control": 0.0, "Treatment_A": 0.5}

        def mock_load(agg_dir):
            label = agg_dir.parent.parent.name
            return _make_mock_agg_result(n_pairs=2, offset=offsets.get(label, 0.0))

        with patch.object(analysis, "_load_aggregated_result", side_effect=mock_load):
            result = analysis.compare(ctx)

        assert result.anova_by_pair is None

    def test_compare_returns_result_with_single_condition(self, tmp_path):
        analysis, ctx, labels = self._make_analysis_and_ctx(tmp_path, n_conditions=1)

        def mock_load(agg_dir):
            return _make_mock_agg_result(n_pairs=2)

        with patch.object(analysis, "_load_aggregated_result", side_effect=mock_load):
            result = analysis.compare(ctx)

        assert result is not None
        assert len(result.conditions) == 1
        assert result.pairwise_comparisons == []

    def test_compare_no_control(self, tmp_path):
        analysis, ctx, labels = self._make_analysis_and_ctx(tmp_path, n_conditions=3, control=None)
        # Override control_label to None
        from polyzymd.analyses.base import ComparisonContext

        ctx = ComparisonContext(
            name=ctx.name,
            conditions=ctx.conditions,
            excluded_conditions=ctx.excluded_conditions,
            control_label=None,
            analysis_dirs=ctx.analysis_dirs,
            results_dir=ctx.results_dir,
            equilibration=ctx.equilibration,
            settings=ctx.settings,
            recompute=ctx.recompute,
        )

        offsets = {"Control": 0.0, "Treatment_A": 0.5, "Treatment_B": -0.3}

        def mock_load(agg_dir):
            label = agg_dir.parent.parent.name
            return _make_mock_agg_result(n_pairs=2, offset=offsets.get(label, 0.0))

        with patch.object(analysis, "_load_aggregated_result", side_effect=mock_load):
            result = analysis.compare(ctx)

        assert result is not None
        # Without control, should have C(3,2) * n_pairs = 3 * 2 = 6 comparisons
        assert len(result.pairwise_comparisons) == 6

    def test_compare_fraction_ranking(self, tmp_path):
        analysis, ctx, labels = self._make_analysis_and_ctx(tmp_path, n_conditions=2)

        offsets = {"Control": 0.0, "Treatment_A": 0.5}

        def mock_load(agg_dir):
            label = agg_dir.parent.parent.name
            return _make_mock_agg_result(n_pairs=2, offset=offsets.get(label, 0.0))

        with patch.object(analysis, "_load_aggregated_result", side_effect=mock_load):
            result = analysis.compare(ctx)

        # fraction_ranking_by_pair should be present since mock has fractions
        assert result.fraction_ranking_by_pair is not None


# ---------------------------------------------------------------------------
# _compare_pair
# ---------------------------------------------------------------------------


class TestComparePair:
    """_compare_pair performs correct statistical comparison for one pair."""

    def test_distance_direction_closer(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        pair_a = MagicMock()
        pair_a.per_replicate_means = [4.0, 4.1, 3.9]
        pair_a.mean_distance = 4.0
        pair_a.per_replicate_fractions = None

        pair_b = MagicMock()
        pair_b.per_replicate_means = [3.0, 3.1, 2.9]
        pair_b.mean_distance = 3.0
        pair_b.per_replicate_fractions = None

        comp = DistancesAnalysis._compare_pair("test_pair", "A", "B", pair_a, pair_b)
        assert comp.distance_direction == "closer"
        assert comp.distance_percent_change < 0

    def test_distance_direction_farther(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        pair_a = MagicMock()
        pair_a.per_replicate_means = [3.0, 3.1, 2.9]
        pair_a.mean_distance = 3.0
        pair_a.per_replicate_fractions = None

        pair_b = MagicMock()
        pair_b.per_replicate_means = [5.0, 5.1, 4.9]
        pair_b.mean_distance = 5.0
        pair_b.per_replicate_fractions = None

        comp = DistancesAnalysis._compare_pair("test_pair", "A", "B", pair_a, pair_b)
        assert comp.distance_direction == "farther"

    def test_fraction_comparison_when_present(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        pair_a = MagicMock()
        pair_a.per_replicate_means = [4.0, 4.1, 3.9]
        pair_a.mean_distance = 4.0
        pair_a.per_replicate_fractions = [0.5, 0.55, 0.45]
        pair_a.fraction_below_threshold = 0.5

        pair_b = MagicMock()
        pair_b.per_replicate_means = [3.0, 3.1, 2.9]
        pair_b.mean_distance = 3.0
        pair_b.per_replicate_fractions = [0.7, 0.75, 0.65]
        pair_b.fraction_below_threshold = 0.7

        comp = DistancesAnalysis._compare_pair("test_pair", "A", "B", pair_a, pair_b)
        assert comp.fraction_t_statistic is not None
        assert comp.fraction_direction == "more_contact"

    def test_no_fraction_when_absent(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        pair_a = MagicMock()
        pair_a.per_replicate_means = [4.0, 4.1]
        pair_a.mean_distance = 4.0
        pair_a.per_replicate_fractions = None

        pair_b = MagicMock()
        pair_b.per_replicate_means = [3.0, 3.1]
        pair_b.mean_distance = 3.0
        pair_b.per_replicate_fractions = None

        comp = DistancesAnalysis._compare_pair("test_pair", "A", "B", pair_a, pair_b)
        assert comp.fraction_t_statistic is None


# ---------------------------------------------------------------------------
# AggregatedResultClass and _deserialize_result
# ---------------------------------------------------------------------------


class TestDeserializeResult:
    """AggregatedResultClass loads via DistanceAggregatedResult.load()."""

    def test_aggregated_result_class_set(self):
        from polyzymd.analyses.distances import DistancesAnalysis
        from polyzymd.analyses.distances._results import DistanceAggregatedResult

        assert DistancesAnalysis.AggregatedResultClass is DistanceAggregatedResult

    def test_loads_via_result_class(self, tmp_path):
        from polyzymd.analyses.distances import DistancesAnalysis

        analysis = DistancesAnalysis()

        mock_result = MagicMock()
        with patch.object(
            analysis.AggregatedResultClass, "load", return_value=mock_result
        ) as mock_load:
            result = analysis._deserialize_result(tmp_path / "test.json")

        mock_load.assert_called_once_with(tmp_path / "test.json")
        assert result is mock_result


# ---------------------------------------------------------------------------
# _make_aggregated_filename
# ---------------------------------------------------------------------------


class TestMakeAggregatedFilename:
    """Filename generation matches existing convention."""

    def test_consecutive_replicates(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        mock_result = MagicMock()
        mock_result.equilibration_time = 100.0
        mock_result.equilibration_unit = "ns"

        filename = DistancesAnalysis._make_aggregated_filename((1, 2, 3), mock_result, "a1b2c3d4")
        assert filename == "distances_reps1-3_eq100.00ns_sa1b2c3d4.json"

    def test_nonconsecutive_replicates(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        mock_result = MagicMock()
        mock_result.equilibration_time = 50.0
        mock_result.equilibration_unit = "ns"

        filename = DistancesAnalysis._make_aggregated_filename((1, 3, 5), mock_result, "a1b2c3d4")
        assert filename == "distances_reps1_3_5_eq50.00ns_sa1b2c3d4.json"

    def test_ps_equilibration(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        mock_result = MagicMock()
        mock_result.equilibration_time = 5000.0
        mock_result.equilibration_unit = "ps"

        filename = DistancesAnalysis._make_aggregated_filename((1, 2), mock_result, "a1b2c3d4")
        assert filename == "distances_reps1-2_eq5000.00ps_sa1b2c3d4.json"


class TestSettingsCacheTag:
    """Settings fingerprinting should differentiate cache identities."""

    def test_cache_tag_changes_with_use_pbc(self):
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )

        settings_true = DistancesSettings(
            use_pbc=True,
            pairs=[DistancePairSettings(label="A", selection_a="resid 1", selection_b="resid 2")],
        )
        settings_false = DistancesSettings(
            use_pbc=False,
            pairs=[DistancePairSettings(label="A", selection_a="resid 1", selection_b="resid 2")],
        )

        tag_true = DistancesAnalysis._make_settings_cache_tag(settings_true)
        tag_false = DistancesAnalysis._make_settings_cache_tag(settings_false)

        assert tag_true != tag_false

    def test_cache_tag_changes_with_pair_definitions(self):
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )

        settings_a = DistancesSettings(
            pairs=[DistancePairSettings(label="A", selection_a="resid 1", selection_b="resid 2")],
        )
        settings_b = DistancesSettings(
            pairs=[DistancePairSettings(label="B", selection_a="resid 10", selection_b="resid 20")],
        )

        tag_a = DistancesAnalysis._make_settings_cache_tag(settings_a)
        tag_b = DistancesAnalysis._make_settings_cache_tag(settings_b)

        assert tag_a != tag_b


# ---------------------------------------------------------------------------
# plot
# ---------------------------------------------------------------------------


class TestPlot:
    """plot() delegates to all three inlined plotting functions."""

    def test_delegates_to_three_plotters(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.distances import DistancesAnalysis
        from polyzymd.config.comparison import PlotSettings

        analysis = DistancesAnalysis()

        cond = Condition(
            label="Cond1",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )

        ctx = PlotContext(
            conditions=[cond],
            analysis_dirs={"Cond1": tmp_path / "analysis"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=MagicMock(),
            plot_settings=PlotSettings(),
        )

        with (
            patch("polyzymd.analyses.distances._plot_distance_kde") as MockKDE,
            patch("polyzymd.analyses.distances._plot_distance_threshold_bars") as MockThreshold,
            patch("polyzymd.analyses.distances._plot_distance_state_bars") as MockState,
        ):
            MockKDE.return_value = [Path("kde.png")]
            MockThreshold.return_value = [Path("bars.png")]
            MockState.return_value = [Path("state.png")]

            plots = analysis.plot(ctx)

        assert len(plots) == 3
        MockKDE.assert_called_once()
        MockThreshold.assert_called_once()
        MockState.assert_called_once()

    def test_plot_handles_failure(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.distances import DistancesAnalysis
        from polyzymd.config.comparison import PlotSettings

        analysis = DistancesAnalysis()

        cond = Condition(
            label="Cond1",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=MagicMock(),
        )

        ctx = PlotContext(
            conditions=[cond],
            analysis_dirs={"Cond1": tmp_path / "analysis"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=MagicMock(),
            plot_settings=PlotSettings(),
        )

        with (
            patch(
                "polyzymd.analyses.distances._plot_distance_kde",
            ) as MockKDE,
            patch("polyzymd.analyses.distances._plot_distance_threshold_bars") as MockThreshold,
            patch("polyzymd.analyses.distances._plot_distance_state_bars") as MockState,
        ):
            MockKDE.side_effect = RuntimeError("KDE failed")
            MockThreshold.return_value = [Path("bars.png")]
            MockState.return_value = [Path("state.png")]

            plots = analysis.plot(ctx)

        # KDE failed, but threshold and state should still produce
        assert len(plots) == 2

    def test_plot_empty_conditions(self, tmp_path):
        from polyzymd.analyses.base import PlotContext
        from polyzymd.analyses.distances import DistancesAnalysis
        from polyzymd.config.comparison import PlotSettings

        analysis = DistancesAnalysis()

        ctx = PlotContext(
            conditions=[],
            analysis_dirs={},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=MagicMock(),
            plot_settings=PlotSettings(),
        )

        plots = analysis.plot(ctx)
        assert plots == []

    def test_plot_creates_default_plot_settings(self, tmp_path):
        """Orchestrator guarantees non-None plot_settings; plugin passes it through."""
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.distances import DistancesAnalysis
        from polyzymd.config.comparison import PlotSettings

        analysis = DistancesAnalysis()

        cond = Condition(
            label="Cond1",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=MagicMock(),
        )

        plot_settings = PlotSettings()
        ctx = PlotContext(
            conditions=[cond],
            analysis_dirs={"Cond1": tmp_path / "analysis"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=MagicMock(),
            plot_settings=plot_settings,
        )

        with (
            patch("polyzymd.analyses.distances._plot_distance_kde") as MockKDE,
            patch("polyzymd.analyses.distances._plot_distance_threshold_bars") as MockThreshold,
            patch("polyzymd.analyses.distances._plot_distance_state_bars") as MockState,
        ):
            MockKDE.return_value = []
            MockThreshold.return_value = []
            MockState.return_value = []

            analysis.plot(ctx)

        # Verify PlotSettings was passed through
        MockKDE.assert_called_once()
        call_args = MockKDE.call_args
        assert isinstance(call_args[0][3], PlotSettings)


# ---------------------------------------------------------------------------
# Plot loader cache filtering
# ---------------------------------------------------------------------------


class TestPlotLoadersCacheFiltering:
    """Plot data loaders should only read distances cache files."""

    def test_pooled_loader_ignores_unrelated_json_files(self, tmp_path):
        import json

        from polyzymd.analyses.distances._plotters import _load_pooled_distances

        analysis_dir = tmp_path / "condition" / "distances"
        run_dir = analysis_dir / "run_1"
        run_dir.mkdir(parents=True)

        distances_payload = {
            "pair_results": [
                {
                    "pair_label": "Distance Pair",
                    "distances": [3.0, 3.2, 3.4],
                    "threshold": 3.5,
                }
            ]
        }
        contacts_payload = {
            "pair_results": [
                {
                    "pair_label": "Contacts Pair",
                    "distances": [10.0, 10.5],
                    "threshold": 11.0,
                }
            ]
        }

        (run_dir / "distances_result.json").write_text(
            json.dumps(distances_payload), encoding="utf-8"
        )
        (run_dir / "contacts_result.json").write_text(
            json.dumps(contacts_payload), encoding="utf-8"
        )
        (run_dir / "notes.json").write_text(json.dumps({"note": "ignore me"}), encoding="utf-8")

        pooled = _load_pooled_distances(analysis_dir, [1])

        assert set(pooled.keys()) == {"Distance Pair"}
        assert pooled["Distance Pair"]["threshold"] == pytest.approx(3.5)

    def test_aggregated_loader_ignores_unrelated_json_files(self, tmp_path):
        import json

        from polyzymd.analyses.distances._plotters import _load_distance_aggregated_results

        aggregated_dir = tmp_path / "condition" / "distances" / "aggregated"
        aggregated_dir.mkdir(parents=True)

        distances_payload = {"pair_results": [{"pair_label": "Distance Pair", "overall_mean": 3.3}]}
        contacts_payload = {"pair_results": [{"pair_label": "Contacts Pair", "overall_mean": 9.9}]}

        (aggregated_dir / "distances_result.json").write_text(
            json.dumps(distances_payload),
            encoding="utf-8",
        )
        (aggregated_dir / "contacts_result.json").write_text(
            json.dumps(contacts_payload), encoding="utf-8"
        )
        (aggregated_dir / "notes.json").write_text(
            json.dumps({"note": "ignore me"}), encoding="utf-8"
        )

        data = {"Cond1": {"aggregated_dir": aggregated_dir}}
        loaded = _load_distance_aggregated_results(data, ["Cond1"])

        assert "Cond1" in loaded
        pair_results = loaded["Cond1"]["pair_results"]
        assert len(pair_results) == 1
        assert pair_results[0]["pair_label"] == "Distance Pair"


# ---------------------------------------------------------------------------
# filter_conditions (default keeps all)
# ---------------------------------------------------------------------------


class TestFilterConditions:
    """Default filter_conditions keeps all conditions."""

    def test_keeps_all_conditions(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.distances import DistancesAnalysis

        analysis = DistancesAnalysis()

        conditions = [
            Condition(
                label="A",
                config_path=Path("/tmp/a.yaml"),
                replicates=(1,),
                sim_config=MagicMock(),
            ),
            Condition(
                label="B",
                config_path=Path("/tmp/b.yaml"),
                replicates=(1,),
                sim_config=MagicMock(),
            ),
        ]

        result = analysis.filter_conditions(conditions)
        assert len(result) == 2


# ---------------------------------------------------------------------------
# Lifecycle test
# ---------------------------------------------------------------------------


class TestLifecycle:
    """End-to-end lifecycle: settings -> compute -> aggregate -> compare."""

    def test_settings_to_compute_to_aggregate_to_compare(self, tmp_path):
        """Verify the full lifecycle flow works end-to-end with mocks."""
        from polyzymd.analyses.base import (
            AggregateContext,
            Condition,
            ReplicateContext,
        )
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )

        analysis = DistancesAnalysis()
        settings = DistancesSettings(
            pairs=[DistancePairSettings(label="Pair1", selection_a="sel_a", selection_b="sel_b")]
        )

        mock_sim_config = MagicMock()
        cond = Condition(
            label="Test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=mock_sim_config,
        )

        # 1. Compute replicates
        rep_results = []
        for rep in [1, 2]:
            ctx = ReplicateContext(
                condition=cond,
                replicate=rep,
                sim_config=mock_sim_config,
                output_dir=tmp_path / f"run_{rep}",
                equilibration="10ns",
                recompute=False,
                settings=settings,
            )
            mock_result = _make_mock_distance_result(replicate=rep, n_pairs=1)
            with patch("polyzymd.analyses.distances.DistanceCalculator") as MockCalc:
                MockCalc.return_value.compute.return_value = mock_result
                result = analysis.compute_replicate(ctx, rep)
                rep_results.append(result)

        assert len(rep_results) == 2

        # 2. Aggregate
        agg_ctx = AggregateContext(
            condition=cond,
            replicates=(1, 2),
            output_dir=tmp_path / "aggregated",
            equilibration="10ns",
            settings=settings,
        )

        with patch(
            "polyzymd.analyses.shared.aggregation.aggregate_distance_pair_stats"
        ) as mock_agg:
            mock_stats = MagicMock()
            mock_stats.mean_stats.mean = 3.5
            mock_stats.mean_stats.sem = 0.05
            mock_stats.median_stats.mean = 3.4
            mock_stats.per_rep_means = [3.49, 3.51]
            mock_stats.per_rep_stds = [0.5, 0.5]
            mock_stats.per_rep_medians = [3.39, 3.41]
            mock_stats.fraction_stats = None
            mock_stats.per_rep_fractions = None
            mock_stats.kde_peak_stats = None
            mock_stats.per_rep_kde_peaks = None
            mock_agg.return_value = mock_stats

            with patch(
                "polyzymd.analyses._results_base.get_polyzymd_version",
                return_value="1.0.0-test",
            ):
                agg_result = analysis.aggregate(agg_ctx, rep_results)

        assert agg_result is not None
        assert agg_result.n_replicates == 2

    def test_extract_metrics_not_needed(self):
        """Since compare() is overridden, extract_metrics returns empty."""
        from polyzymd.analyses.distances import DistancesAnalysis

        analysis = DistancesAnalysis()
        # extract_metrics uses the default (returns {})
        assert analysis.extract_metrics(MagicMock()) == {}
