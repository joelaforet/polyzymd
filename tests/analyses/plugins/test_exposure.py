"""Tests for ExposureAnalysis plugin — Phase C-5.

Coverage:
- Discovery (registered by name)
- Class attributes (name, Settings, aliases, dependencies, min_replicates)
- Settings validation
- compute_replicate / aggregate (no-op)
- filter_conditions (polymer detection)
- compare (full override with dual rankings)
- pairwise comparisons
- ANOVA
- plot delegation
- extract_metrics (empty)
- lifecycle
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from pydantic import ValidationError

# ===========================================================================
# Discovery
# ===========================================================================


class TestDiscovery:
    """Exposure plugin is discovered by pkgutil scanning."""

    def test_discovery_by_name(self):
        from polyzymd.analyses import get_analysis

        cls = get_analysis("exposure")
        assert cls.__name__ == "ExposureAnalysis"

    def test_listed(self):
        from polyzymd.analyses import list_analyses

        analyses = list_analyses()
        assert "exposure" in analyses

    def test_all_names(self):
        from polyzymd.analyses import list_all_names

        names = list_all_names()
        assert "exposure" in names


# ===========================================================================
# Class attributes
# ===========================================================================


class TestClassAttributes:
    def test_name(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        assert ExposureAnalysis.name == "exposure"

    def test_settings_type(self):
        from polyzymd.analyses.exposure import ExposureAnalysis, ExposureSettings

        assert ExposureAnalysis.Settings is ExposureSettings

    def test_aliases(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        assert ExposureAnalysis.aliases == ()

    def test_dependencies(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        assert ExposureAnalysis.dependencies == ("contacts",)

    def test_min_replicates(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        assert ExposureAnalysis.min_replicates == 2


# ===========================================================================
# Settings
# ===========================================================================


class TestSettings:
    def test_defaults(self):
        from polyzymd.analyses.exposure import ExposureSettings

        s = ExposureSettings()
        assert s.protein_selection == "protein"
        assert s.polymer_selection == "chainID C"
        assert s.exposure_threshold == pytest.approx(0.2)
        assert s.transient_lower == pytest.approx(0.2)
        assert s.transient_upper == pytest.approx(0.8)
        assert s.min_event_length == 1
        assert s.probe_radius_nm == pytest.approx(0.14)
        assert s.n_sphere_points == 960
        assert s.protein_chain == "A"
        assert s.polymer_resnames is None

    def test_custom_values(self):
        from polyzymd.analyses.exposure import ExposureSettings

        s = ExposureSettings(
            exposure_threshold=0.3,
            transient_lower=0.1,
            transient_upper=0.9,
            min_event_length=5,
            protein_chain="B",
            polymer_resnames=["SBMA", "OEGMA"],
        )
        assert s.exposure_threshold == pytest.approx(0.3)
        assert s.transient_lower == pytest.approx(0.1)
        assert s.transient_upper == pytest.approx(0.9)
        assert s.min_event_length == 5
        assert s.protein_chain == "B"
        assert s.polymer_resnames == ["SBMA", "OEGMA"]

    def test_invalid_threshold_zero(self):
        from polyzymd.analyses.exposure import ExposureSettings

        with pytest.raises(ValidationError):
            ExposureSettings(exposure_threshold=0.0)

    def test_invalid_threshold_one(self):
        from polyzymd.analyses.exposure import ExposureSettings

        with pytest.raises(ValidationError):
            ExposureSettings(exposure_threshold=1.0)

    def test_invalid_transient_bounds(self):
        from polyzymd.analyses.exposure import ExposureSettings

        with pytest.raises(ValidationError):
            ExposureSettings(transient_lower=0.8, transient_upper=0.2)

    def test_serialization_roundtrip(self):
        from polyzymd.analyses.exposure import ExposureSettings

        s = ExposureSettings(exposure_threshold=0.3, min_event_length=3)
        d = s.model_dump()
        s2 = ExposureSettings(**d)
        assert s2.exposure_threshold == pytest.approx(0.3)
        assert s2.min_event_length == 3


# ===========================================================================
# compute_replicate / aggregate (no-op)
# ===========================================================================


class TestNoOp:
    def test_compute_replicate_returns_none(self):
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.exposure import ExposureAnalysis, ExposureSettings

        analysis = ExposureAnalysis()
        mock_sim = MagicMock()
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim,
        )
        ctx = ReplicateContext(
            condition=cond,
            replicate=1,
            sim_config=mock_sim,
            output_dir=Path("/tmp/out"),
            equilibration="10ns",
            recompute=False,
            settings=ExposureSettings(),
        )
        assert analysis.compute_replicate(ctx, 1) is None

    def test_aggregate_returns_none(self):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.exposure import ExposureAnalysis, ExposureSettings

        analysis = ExposureAnalysis()
        mock_sim = MagicMock()
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim,
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1,),
            output_dir=Path("/tmp/agg"),
            equilibration="10ns",
            settings=ExposureSettings(),
        )
        assert analysis.aggregate(ctx, []) is None


# ===========================================================================
# filter_conditions
# ===========================================================================


class TestFilterConditions:
    def test_keeps_polymer_conditions(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()

        mock_sim = MagicMock()
        mock_sim.polymers = MagicMock()
        mock_sim.polymers.enabled = True

        cond = Condition(
            label="With Polymer",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim,
        )

        result = analysis.filter_conditions([cond])
        assert len(result) == 1
        assert result[0].label == "With Polymer"

    def test_excludes_no_polymer_none(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()

        mock_sim = MagicMock()
        mock_sim.polymers = None

        cond = Condition(
            label="No Polymer",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim,
        )

        result = analysis.filter_conditions([cond])
        assert len(result) == 0

    def test_excludes_disabled_polymer(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()

        mock_sim = MagicMock()
        mock_sim.polymers = MagicMock()
        mock_sim.polymers.enabled = False

        cond = Condition(
            label="Disabled",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim,
        )

        result = analysis.filter_conditions([cond])
        assert len(result) == 0

    def test_error_includes_condition(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()

        mock_sim = MagicMock()
        # Accessing .polymers raises
        type(mock_sim).polymers = property(lambda s: (_ for _ in ()).throw(RuntimeError("boom")))

        cond = Condition(
            label="ErrorCond",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim,
        )

        result = analysis.filter_conditions([cond])
        assert len(result) == 1  # fail-open


# ===========================================================================
# compare
# ===========================================================================


def _make_mock_dynamics(n_residues=10, n_transient=3, chaperone_fraction=0.4):
    """Build a mock ExposureDynamicsResult."""
    mock = MagicMock()
    mock.n_residues = n_residues
    mock.n_transient.return_value = n_transient
    mock.total_chaperone_events.return_value = 5
    mock.total_unassisted_events.return_value = 3

    # Build transient residues with chaperone_fraction
    transient = []
    for _ in range(n_transient):
        r = MagicMock()
        r.chaperone_fraction = chaperone_fraction
        transient.append(r)
    mock.transient_residues.return_value = transient

    return mock


def _make_mock_enrichment(polymer_types=("SBMA",), aa_groups=("charged", "polar")):
    """Build a mock ChaperoneEnrichmentResult."""
    mock = MagicMock()
    entries = []
    for pt in polymer_types:
        for ag in aa_groups:
            entry = MagicMock()
            entry.polymer_type = pt
            entry.aa_group = ag
            entry.enrichment_residue = 0.3
            entries.append(entry)
    mock.entries = entries

    def _get(ptype, ag):
        for e in entries:
            if e.polymer_type == ptype and e.aa_group == ag:
                return e
        return None

    mock.get = _get
    return mock


class TestCompare:
    def _make_context(self, n_conditions=2, n_reps=3, control="A"):
        """Build a ComparisonContext with mock conditions."""
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.exposure import ExposureSettings

        conditions = []
        analysis_dirs = {}
        for i in range(n_conditions):
            label = chr(65 + i)  # A, B, C, ...
            mock_sim = MagicMock()
            mock_sim.polymers = MagicMock()
            mock_sim.polymers.enabled = True
            cond = Condition(
                label=label,
                config_path=Path(f"/tmp/{label}/config.yaml"),
                replicates=tuple(range(1, n_reps + 1)),
                sim_config=mock_sim,
            )
            conditions.append(cond)
            analysis_dirs[label] = Path(f"/tmp/{label}/exposure")

        return ComparisonContext(
            name="test_comparison",
            conditions=conditions,
            excluded_conditions=[],
            control_label=control,
            analysis_dirs=analysis_dirs,
            results_dir=Path("/tmp/results"),
            equilibration="10ns",
            settings=ExposureSettings(),
            recompute=False,
        )

    def test_compare_returns_result(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        ctx = self._make_context(n_conditions=2)

        dynamics = _make_mock_dynamics()
        enrichment = _make_mock_enrichment()

        with patch.object(
            analysis,
            "_load_or_compute_replicate",
            return_value=(dynamics, enrichment),
        ):
            result = analysis.compare(ctx)

        assert result is not None
        assert result.metric == "chaperone_fraction"
        assert len(result.conditions) == 2
        assert len(result.ranking) == 2

    def test_compare_pairwise_with_control(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        ctx = self._make_context(n_conditions=3, control="A")

        dynamics = _make_mock_dynamics()
        enrichment = _make_mock_enrichment()

        with patch.object(
            analysis,
            "_load_or_compute_replicate",
            return_value=(dynamics, enrichment),
        ):
            result = analysis.compare(ctx)

        # With control A and 3 conditions, expect 2 pairwise (A vs B, A vs C)
        assert len(result.pairwise_comparisons) == 2
        for comp in result.pairwise_comparisons:
            assert comp.condition_a == "A"

    def test_compare_no_control(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        ctx = self._make_context(n_conditions=3, control=None)

        dynamics = _make_mock_dynamics()
        enrichment = _make_mock_enrichment()

        with patch.object(
            analysis,
            "_load_or_compute_replicate",
            return_value=(dynamics, enrichment),
        ):
            result = analysis.compare(ctx)

        # All pairs: AB, AC, BC = 3
        assert len(result.pairwise_comparisons) == 3

    def test_compare_anova_with_3_conditions(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        ctx = self._make_context(n_conditions=3)

        dynamics = _make_mock_dynamics()
        enrichment = _make_mock_enrichment()

        with patch.object(
            analysis,
            "_load_or_compute_replicate",
            return_value=(dynamics, enrichment),
        ):
            result = analysis.compare(ctx)

        assert result.anova is not None
        assert result.anova.metric == "chaperone_fraction"

    def test_compare_no_anova_with_2_conditions(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        ctx = self._make_context(n_conditions=2)

        dynamics = _make_mock_dynamics()
        enrichment = _make_mock_enrichment()

        with patch.object(
            analysis,
            "_load_or_compute_replicate",
            return_value=(dynamics, enrichment),
        ):
            result = analysis.compare(ctx)

        assert result.anova is None

    def test_compare_returns_none_single_condition(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        ctx = self._make_context(n_conditions=1, control=None)

        result = analysis.compare(ctx)
        assert result is None

    def test_compare_dual_rankings(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        ctx = self._make_context(n_conditions=2)

        dynamics = _make_mock_dynamics()
        enrichment = _make_mock_enrichment()

        with patch.object(
            analysis,
            "_load_or_compute_replicate",
            return_value=(dynamics, enrichment),
        ):
            result = analysis.compare(ctx)

        assert len(result.ranking) == 2
        assert len(result.ranking_by_transient_fraction) == 2

    def test_compare_enrichment_aggregated(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        ctx = self._make_context(n_conditions=2)

        dynamics = _make_mock_dynamics()
        enrichment = _make_mock_enrichment(polymer_types=("SBMA",), aa_groups=("charged", "polar"))

        with patch.object(
            analysis,
            "_load_or_compute_replicate",
            return_value=(dynamics, enrichment),
        ):
            result = analysis.compare(ctx)

        for cond_summary in result.conditions:
            assert "SBMA" in cond_summary.polymer_types
            assert "charged" in cond_summary.aa_groups
            assert "SBMA" in cond_summary.enrichment_by_polymer_type

    def test_compare_respects_context_fdr_alpha_for_pairwise_and_anova(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        ctx = self._make_context(n_conditions=3, control="A")
        ctx = type(ctx)(
            name=ctx.name,
            conditions=ctx.conditions,
            excluded_conditions=ctx.excluded_conditions,
            control_label=ctx.control_label,
            analysis_dirs=ctx.analysis_dirs,
            results_dir=ctx.results_dir,
            equilibration=ctx.equilibration,
            settings=ctx.settings,
            recompute=ctx.recompute,
            fdr_alpha=0.2,
            ttest_method="welch",
            posthoc_method="ttest_bh",
        )

        dynamics = _make_mock_dynamics()
        enrichment = _make_mock_enrichment()

        class _FakeTTest:
            t_statistic = 1.23
            p_value = 0.2

            @property
            def significant(self):
                return False

        class _FakeAnova:
            f_statistic = 4.2
            p_value = 0.2

            @property
            def significant(self):
                return False

        with (
            patch.object(
                analysis,
                "_load_or_compute_replicate",
                return_value=(dynamics, enrichment),
            ),
            patch(
                "polyzymd.analyses.shared.inferential_statistics.independent_ttest",
                return_value=_FakeTTest(),
            ) as mock_ttest,
            patch(
                "polyzymd.analyses.shared.inferential_statistics.one_way_anova",
                return_value=_FakeAnova(),
            ),
        ):
            result = analysis.compare(ctx)

        assert result is not None
        assert result.fdr_alpha == pytest.approx(0.2)
        assert result.ttest_method == "welch"
        assert result.posthoc_method == "ttest_bh"
        assert result.anova is not None
        assert result.anova.significant is True
        assert all(comp.significant is True for comp in result.pairwise_comparisons)
        assert mock_ttest.call_count == 2
        assert all(call.kwargs["method"] == "welch" for call in mock_ttest.call_args_list)

    def test_compare_populates_adjusted_p_values(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        ctx = self._make_context(n_conditions=3, control=None)

        dynamics = _make_mock_dynamics()
        enrichment = _make_mock_enrichment()

        class _FakeTTest:
            def __init__(self, p_value: float):
                self.t_statistic = 1.0
                self.p_value = p_value

        with (
            patch.object(
                analysis,
                "_load_or_compute_replicate",
                return_value=(dynamics, enrichment),
            ),
            patch(
                "polyzymd.analyses.shared.inferential_statistics.independent_ttest",
                side_effect=[_FakeTTest(0.001), _FakeTTest(0.049), _FakeTTest(0.9)],
            ),
        ):
            result = analysis.compare(ctx)

        assert result is not None
        assert len(result.pairwise_comparisons) == 3
        assert all(comp.p_value_adjusted is not None for comp in result.pairwise_comparisons)

    def test_compare_significance_uses_bh_corrected_p_values(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        ctx = self._make_context(n_conditions=3, control=None)

        dynamics = _make_mock_dynamics()
        enrichment = _make_mock_enrichment()

        class _FakeTTest:
            def __init__(self, p_value: float):
                self.t_statistic = 1.0
                self.p_value = p_value

        with (
            patch.object(
                analysis,
                "_load_or_compute_replicate",
                return_value=(dynamics, enrichment),
            ),
            patch(
                "polyzymd.analyses.shared.inferential_statistics.independent_ttest",
                side_effect=[_FakeTTest(0.001), _FakeTTest(0.049), _FakeTTest(0.9)],
            ),
        ):
            result = analysis.compare(ctx)

        assert result is not None
        raw_significant = [comp.p_value <= 0.05 for comp in result.pairwise_comparisons]
        corrected_significant = [comp.significant for comp in result.pairwise_comparisons]

        assert any(raw_significant)
        assert any(
            raw and (not corrected)
            for raw, corrected in zip(raw_significant, corrected_significant)
        )


# ===========================================================================
# Plot
# ===========================================================================


class TestPlot:
    @patch("polyzymd.analyses.exposure._plot_enrichment_heatmap")
    @patch("polyzymd.analyses.exposure._plot_chaperone_fraction")
    def test_plot_creates_output_dir(self, mock_chaperone_fn, mock_heatmap_fn, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.exposure import ExposureAnalysis, ExposureSettings
        from polyzymd.config.comparison import PlotSettings

        analysis = ExposureAnalysis()
        mock_chaperone_fn.return_value = []
        mock_heatmap_fn.return_value = []

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
            analysis_dirs={"A": tmp_path / "A" / "exposure"},
            results_dir=tmp_path / "results",
            output_dir=output_dir,
            settings=ExposureSettings(),
            plot_settings=PlotSettings(),
        )

        analysis.plot(ctx)

        assert output_dir.exists()

    @patch(
        "polyzymd.analyses.exposure._plot_enrichment_heatmap",
        side_effect=Exception("plot error"),
    )
    @patch(
        "polyzymd.analyses.exposure._plot_chaperone_fraction",
        side_effect=Exception("plot error"),
    )
    def test_plot_propagates_exceptions(self, mock_chaperone_fn, mock_heatmap_fn, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.exposure import ExposureAnalysis, ExposureSettings
        from polyzymd.config.comparison import PlotSettings

        analysis = ExposureAnalysis()

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
            analysis_dirs={"A": tmp_path / "A" / "exposure"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=ExposureSettings(),
            plot_settings=PlotSettings(),
        )

        with pytest.raises(Exception, match="plot error"):
            analysis.plot(ctx)


# ===========================================================================
# extract_metrics
# ===========================================================================


class TestExtractMetrics:
    def test_returns_empty(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        assert analysis.extract_metrics(MagicMock()) == {}


# ===========================================================================
# Lifecycle
# ===========================================================================


class TestLifecycle:
    def test_instantiation(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        assert analysis.name == "exposure"

    def test_is_subclass_of_analysis(self):
        from polyzymd.analyses.base import Analysis
        from polyzymd.analyses.exposure import ExposureAnalysis

        assert issubclass(ExposureAnalysis, Analysis)

    def test_has_all_required_methods(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        assert callable(analysis.compute_replicate)
        assert callable(analysis.aggregate)
        assert callable(analysis.compare)
        assert callable(analysis.filter_conditions)
        assert callable(analysis.plot)
        assert callable(analysis.extract_metrics)

    def test_repr(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        assert "exposure" in repr(analysis)


# ===========================================================================
# _find_contact_result
# ===========================================================================


class TestFindContactResult:
    def test_finds_in_contacts_dir(self, tmp_path):
        from polyzymd.analyses.exposure import ExposureAnalysis

        contacts_dir = tmp_path / "contacts"
        contacts_dir.mkdir()
        (contacts_dir / "contacts_rep1.json").write_text("{}")

        result = ExposureAnalysis._find_contact_result(contacts_dir, 1)
        assert result is not None
        assert result.name == "contacts_rep1.json"

    def test_returns_none_when_missing(self, tmp_path):
        from polyzymd.analyses.exposure import ExposureAnalysis

        contacts_dir = tmp_path / "contacts"
        contacts_dir.mkdir()

        result = ExposureAnalysis._find_contact_result(contacts_dir, 1)
        assert result is None

    def test_returns_none_when_dir_is_none(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        result = ExposureAnalysis._find_contact_result(None, 1)
        assert result is None

    def test_raises_on_ambiguous_parameterized_matches(self, tmp_path):
        from polyzymd.analyses.exposure import ExposureAnalysis

        contacts_dir = tmp_path / "contacts"
        contacts_dir.mkdir()
        (contacts_dir / "contacts_eq10ns_cut4.0_rep1.json").write_text("{}")
        (contacts_dir / "contacts_eq20ns_cut4.5_rep1.json").write_text("{}")

        with pytest.raises(ValueError, match="Ambiguous contacts cache"):
            ExposureAnalysis._find_contact_result(contacts_dir, 1)


# ===========================================================================
# _condition_has_polymer
# ===========================================================================


class TestConditionHasPolymer:
    def test_true_when_polymer_configured(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        mock_sim = MagicMock()
        mock_sim.polymers = MagicMock()
        mock_sim.polymers.enabled = True

        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim,
        )
        assert analysis._condition_has_polymer(cond) is True

    def test_false_when_polymers_none(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        mock_sim = MagicMock()
        mock_sim.polymers = None

        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim,
        )
        assert analysis._condition_has_polymer(cond) is False

    def test_false_when_polymers_disabled(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        mock_sim = MagicMock()
        mock_sim.polymers = MagicMock()
        mock_sim.polymers.enabled = False

        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim,
        )
        assert analysis._condition_has_polymer(cond) is False


# ===========================================================================
# Formatter enrichment markers
# ===========================================================================


class TestExposureFormatterEnrichmentMarkers:
    def _make_result_with_enrichment(self, enrichment_value: float):
        from datetime import datetime

        from polyzymd.analyses.exposure._comparison_results import (
            ExposureComparisonResult,
            ExposureConditionSummary,
        )

        condition = ExposureConditionSummary(
            label="A",
            config_path="/tmp/a/config.yaml",
            n_replicates=3,
            replicate_values=[0.4, 0.5, 0.6],
            mean_transient_fraction=0.3,
            sem_transient_fraction=0.01,
            mean_chaperone_fraction=0.5,
            sem_chaperone_fraction=0.02,
            mean_n_transient=10.0,
            mean_total_chaperone_events=8.0,
            mean_total_unassisted_events=2.0,
            enrichment_by_polymer_type={"SBMA": {"charged": enrichment_value}},
            polymer_types=["SBMA"],
            aa_groups=["charged"],
        )
        return ExposureComparisonResult(
            name="test",
            metric="chaperone_fraction",
            control_label=None,
            fdr_alpha=0.05,
            ttest_method="student",
            posthoc_method="ttest_bh",
            conditions=[condition],
            pairwise_comparisons=[],
            anova=None,
            ranking=["A"],
            ranking_by_transient_fraction=["A"],
            excluded_conditions=[],
            equilibration_time="10ns",
            created_at=datetime.now(),
            polyzymd_version="test",
        )

    @pytest.mark.parametrize(
        ("enrichment_value", "expected_marker"),
        [(0.5, "+"), (-0.3, "-"), (0.0, "")],
    )
    def test_markdown_formatter_uses_zero_centered_markers(
        self,
        enrichment_value: float,
        expected_marker: str,
    ):
        from polyzymd.analyses.exposure._formatters import format_exposure_markdown

        result = self._make_result_with_enrichment(enrichment_value)
        text = format_exposure_markdown(result, show_pairwise=False, show_anova=False)

        assert "> **Key:** + = enriched (>0), - = depleted (<0)" in text
        if expected_marker:
            assert f"{enrichment_value:.2f} {expected_marker}" in text
        else:
            assert f"{enrichment_value:.2f} +" not in text
            assert f"{enrichment_value:.2f} -" not in text

    @pytest.mark.parametrize(
        ("enrichment_value", "expected_marker"),
        [(0.5, "+"), (-0.3, "-"), (0.0, "")],
    )
    def test_console_formatter_uses_zero_centered_markers(
        self,
        enrichment_value: float,
        expected_marker: str,
    ):
        from polyzymd.analyses.exposure._formatters import format_exposure_console_table

        result = self._make_result_with_enrichment(enrichment_value)
        text = format_exposure_console_table(result, show_pairwise=False, show_anova=False)

        assert "+ = enriched (>0), - = depleted (<0)" in text
        if expected_marker:
            assert f"{enrichment_value:>10.2f}{expected_marker} " in text
        else:
            assert f"{enrichment_value:>10.2f}+ " not in text
            assert f"{enrichment_value:>10.2f}- " not in text


# ===========================================================================
# Cache fingerprinting and chain validation
# ===========================================================================


class TestExposureCacheFingerprinting:
    def test_sasa_cache_path_changes_with_probe_radius_and_sphere_points(self, tmp_path):
        from polyzymd.analyses.exposure._sasa_config import SASAConfig
        from polyzymd.analyses.exposure._sasa_trajectory import SASATrajectoryResult
        from polyzymd.analyses.shared.config_hash import settings_fingerprint

        cfg_a = SASAConfig(probe_radius_nm=0.14, n_sphere_points=960)
        cfg_b = SASAConfig(probe_radius_nm=0.16, n_sphere_points=960)
        cfg_c = SASAConfig(probe_radius_nm=0.14, n_sphere_points=1000)

        path_a = SASATrajectoryResult.cache_path(tmp_path, settings_fp=settings_fingerprint(cfg_a))
        path_b = SASATrajectoryResult.cache_path(tmp_path, settings_fp=settings_fingerprint(cfg_b))
        path_c = SASATrajectoryResult.cache_path(tmp_path, settings_fp=settings_fingerprint(cfg_c))

        assert path_a != path_b
        assert path_a != path_c

    def test_dynamics_cache_path_changes_with_event_settings(self, tmp_path):
        from polyzymd.analyses.exposure._config import ExposureConfig
        from polyzymd.analyses.exposure._dynamics import ExposureDynamicsResult
        from polyzymd.analyses.shared.config_hash import settings_fingerprint

        cfg_a = ExposureConfig(transient_lower=0.2, transient_upper=0.8, min_event_length=1)
        cfg_b = ExposureConfig(transient_lower=0.25, transient_upper=0.8, min_event_length=1)
        cfg_c = ExposureConfig(transient_lower=0.2, transient_upper=0.8, min_event_length=3)

        path_a = ExposureDynamicsResult.cache_path(
            tmp_path, settings_fp=settings_fingerprint(cfg_a)
        )
        path_b = ExposureDynamicsResult.cache_path(
            tmp_path, settings_fp=settings_fingerprint(cfg_b)
        )
        path_c = ExposureDynamicsResult.cache_path(
            tmp_path, settings_fp=settings_fingerprint(cfg_c)
        )

        assert path_a != path_b
        assert path_a != path_c


class TestSASAChainSelection:
    @patch("mdtraj.load")
    def test_missing_chain_raises_value_error(self, mock_load, tmp_path):
        from polyzymd.analyses.exposure._sasa_config import SASAConfig
        from polyzymd.analyses.exposure._sasa_trajectory import compute_trajectory_sasa

        mock_topology = MagicMock()
        mock_topology.select.return_value = []
        mock_chain = MagicMock()
        mock_chain.index = 0
        mock_topology.chains = [mock_chain]

        mock_traj = MagicMock()
        mock_traj.topology = mock_topology
        mock_load.return_value = mock_traj

        with pytest.raises(ValueError, match="Chain 'Z' not found in topology"):
            compute_trajectory_sasa(
                topology_path=tmp_path / "top.pdb",
                trajectory_path=tmp_path / "traj.xtc",
                config=SASAConfig(chain_id="Z"),
                analysis_dir=tmp_path / "analysis",
                recompute=True,
            )
