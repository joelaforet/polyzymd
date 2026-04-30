"""Tests for ExposureAnalysis plugin — Phase C-5.

Coverage:
- Discovery (registered by name)
- Class attributes (name, Settings, aliases, dependencies, min_replicates)
- Settings validation
- run_replicate / aggregate (no-op)
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
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
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
# run_replicate / aggregate (no-op)
# ===========================================================================


class TestNoOp:
    def test_run_replicate_returns_none(self):
        """run_replicate is a no-op for compare-only exposure analysis.

        Exposure sets ``has_compute_stage = False``, so the base Analysis
        implementation returns ``None`` instead of requiring per-replicate
        compute behavior.
        """
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
        assert analysis.run_replicate(ctx, 1) is None

    def test_aggregate_returns_none(self):
        """aggregate is a no-op because exposure is comparator-only.

        Exposure sets ``has_aggregate_stage = False``, so aggregation is
        intentionally skipped and handled inside ``compare()``.
        """
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

    def test_compare_method_overrides_default_analysis_compare(self):
        """Exposure must override compare for custom multi-metric flow."""
        from polyzymd.analyses.base import Analysis
        from polyzymd.analyses.exposure import ExposureAnalysis

        assert ExposureAnalysis.compare is not Analysis.compare


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
        type(mock_sim).polymers = property(lambda s: (_ for _ in ()).throw(OSError("boom")))

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
        assert [cond.label for cond in result.conditions] == ["A", "B"]

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
        assert {comp.condition_b for comp in result.pairwise_comparisons} == {"B", "C"}
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

    def test_compare_ranking_prefers_higher_chaperone_fraction(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        ctx = self._make_context(n_conditions=2, n_reps=3, control=None)

        def _by_condition(*, cond, **_kwargs):
            if cond.label == "A":
                return _make_mock_dynamics(chaperone_fraction=0.2), _make_mock_enrichment()
            return _make_mock_dynamics(chaperone_fraction=0.8), _make_mock_enrichment()

        with patch.object(analysis, "_load_or_compute_replicate", side_effect=_by_condition):
            result = analysis.compare(ctx)

        assert result is not None
        assert result.ranking == ["B", "A"]
        means = {cond.label: cond.mean_chaperone_fraction for cond in result.conditions}
        assert means["A"] == pytest.approx(0.2)
        assert means["B"] == pytest.approx(0.8)

    def test_compare_pairwise_percent_change_uses_condition_means(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        ctx = self._make_context(n_conditions=2, n_reps=3, control="A")

        def _by_condition(*, cond, **_kwargs):
            if cond.label == "A":
                return _make_mock_dynamics(chaperone_fraction=0.2), _make_mock_enrichment()
            return _make_mock_dynamics(chaperone_fraction=0.8), _make_mock_enrichment()

        with patch.object(analysis, "_load_or_compute_replicate", side_effect=_by_condition):
            result = analysis.compare(ctx)

        assert result is not None
        assert len(result.pairwise_comparisons) == 1
        pair = result.pairwise_comparisons[0]
        assert pair.condition_a == "A"
        assert pair.condition_b == "B"
        assert pair.percent_change == pytest.approx(300.0)

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
        side_effect=RuntimeError("plot error"),
    )
    @patch(
        "polyzymd.analyses.exposure._plot_chaperone_fraction",
        side_effect=RuntimeError("plot error"),
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

        with pytest.raises(RuntimeError, match="plot error"):
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

    def test_class_contract_flags(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        assert ExposureAnalysis.name == "exposure"
        assert ExposureAnalysis.aliases == ()
        assert ExposureAnalysis.dependencies == ("contacts",)
        assert ExposureAnalysis.min_replicates == 2
        assert ExposureAnalysis.has_compute_stage is False
        assert ExposureAnalysis.has_aggregate_stage is False
        assert ExposureAnalysis.execution_cost_hint == "high"

    def test_repr(self):
        from polyzymd.analyses.exposure import ExposureAnalysis

        analysis = ExposureAnalysis()
        assert "exposure" in repr(analysis)


# ===========================================================================
# _find_contact_result
# ===========================================================================


class TestFindContactResult:
    def test_prefers_canonical_run_result(self, tmp_path):
        from polyzymd.analyses.exposure import ExposureAnalysis

        contacts_dir = tmp_path / "contacts"
        run_dir = contacts_dir / "run_1"
        run_dir.mkdir(parents=True)
        canonical = run_dir / "result.json"
        canonical.write_text("{}")
        (contacts_dir / "contacts_rep1.json").write_text("{}")

        result = ExposureAnalysis._find_contact_result(contacts_dir, 1)

        assert result == canonical

    def test_finds_fingerprinted_sidecar_when_canonical_absent(self, tmp_path):
        from polyzymd.analyses.exposure import ExposureAnalysis

        contacts_dir = tmp_path / "contacts"
        contacts_dir.mkdir()
        sidecar = contacts_dir / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json"
        sidecar.write_text('{"contacts_settings_fingerprint": "deadbeef"}')

        result = ExposureAnalysis._find_contact_result(
            contacts_dir,
            1,
            settings_fp="deadbeef",
        )

        assert result == sidecar

    def test_find_contact_result_filters_by_equilibration(self, tmp_path):
        from polyzymd.analyses.exposure import ExposureAnalysis

        contacts_dir = tmp_path / "contacts"
        contacts_dir.mkdir()
        wrong = contacts_dir / "contacts_eq0ns_cut4.5_sdeadbeef_rep1.json"
        expected = contacts_dir / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json"
        wrong.write_text("{}")
        expected.write_text('{"contacts_settings_fingerprint": "deadbeef"}')

        result = ExposureAnalysis._find_contact_result(
            contacts_dir,
            1,
            settings_fp="deadbeef",
            equilibration="10ns",
        )

        assert result == expected

    def test_matching_fingerprinted_sidecar_precedes_stale_canonical(self, tmp_path):
        from polyzymd.analyses.exposure import ExposureAnalysis

        contacts_dir = tmp_path / "contacts"
        run_dir = contacts_dir / "run_1"
        run_dir.mkdir(parents=True)
        canonical = run_dir / "result.json"
        sidecar = contacts_dir / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json"
        canonical.write_text("{}")
        sidecar.write_text('{"contacts_settings_fingerprint": "deadbeef"}')

        result = ExposureAnalysis._find_contact_result(
            contacts_dir,
            1,
            settings_fp="deadbeef",
        )

        assert result == sidecar

    def test_load_replicate_accepts_actual_non_default_contacts_fingerprint(self, tmp_path):
        """Exposure should trust the upstream contacts artifact identity."""

        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint
        from polyzymd.analyses.exposure import ExposureAnalysis
        from polyzymd.analyses.shared.config_hash import settings_fingerprint

        analysis = ExposureAnalysis()
        settings = ExposureAnalysis.Settings()
        contacts_settings = ContactsSettings(
            protein_selection=settings.protein_selection,
            polymer_selection=settings.polymer_selection,
            cutoff=6.0,
            grouping="none",
        )
        contacts_fp = contacts_settings_fingerprint(contacts_settings)
        default_contacts_fp = contacts_settings_fingerprint(
            ContactsSettings(
                protein_selection=settings.protein_selection,
                polymer_selection=settings.polymer_selection,
            )
        )
        assert contacts_fp != default_contacts_fp

        contacts_dir = tmp_path / "contacts"
        contacts_dir.mkdir()
        contact_file = contacts_dir / f"contacts_eq10ns_cut6.0_s{contacts_fp}_rep1.json"
        contact_file.write_text(f'{{"contacts_settings_fingerprint": "{contacts_fp}"}}')

        mock_sim = MagicMock()
        mock_sim.polymers = MagicMock(enabled=True)
        cond = Condition(
            label="test",
            config_path=tmp_path / "config.yaml",
            replicates=(1,),
            sim_config=mock_sim,
        )

        with (
            patch("polyzymd.analyses.contacts._results.ContactResult") as mock_contact_result,
            patch("polyzymd.analyses.shared.loader.TrajectoryLoader") as mock_loader_cls,
            patch(
                "polyzymd.analyses.exposure._sasa_trajectory.compute_trajectory_sasa"
            ) as mock_sasa,
            patch(
                "polyzymd.analyses.exposure._dynamics.analyze_exposure_dynamics"
            ) as mock_dynamics,
            patch(
                "polyzymd.analyses.exposure._enrichment.compute_chaperone_enrichment"
            ) as mock_enrichment,
        ):
            mock_contact_result.load.return_value = SimpleNamespace(
                metadata={"settings_fingerprint": contacts_fp},
                equilibration_time=10.0,
                equilibration_unit="ns",
                replicate=1,
                n_frames=5,
            )
            mock_loader = MagicMock()
            traj_info = MagicMock()
            traj_info.topology_file = tmp_path / "top.pdb"
            traj_info.trajectory_files = [tmp_path / "traj.xtc"]
            mock_loader.get_trajectory_info.return_value = traj_info
            mock_loader_cls.return_value = mock_loader
            mock_sasa.return_value = MagicMock(n_frames=5)
            mock_dynamics.return_value = MagicMock()
            mock_enrichment.return_value = MagicMock()

            result = analysis._load_or_compute_replicate(
                cond=cond,
                replicate=1,
                settings=settings,
                exposure_analysis_dir=tmp_path / "exposure",
                contacts_analysis_dir=contacts_dir,
                contacts_settings_fp=analysis._infer_contacts_settings_fingerprint(
                    contacts_dir,
                    cond.replicates,
                ),
                recompute=True,
                equilibration="10ns",
            )

        assert result is not None
        mock_contact_result.load.assert_called_once_with(contact_file)
        expected_identity = analysis._contacts_artifact_identity(
            mock_contact_result.load.return_value,
            contact_file,
        )
        assert "contacts_artifact_identity" not in mock_sasa.call_args.kwargs
        assert mock_dynamics.call_args.kwargs["contacts_artifact_identity"] == expected_identity

    def test_contact_result_settings_mismatch_is_rejected(self, tmp_path):
        from polyzymd.analyses.exposure import ExposureAnalysis

        contact_result = SimpleNamespace(metadata={"settings_fingerprint": "badcafe0"})

        assert (
            ExposureAnalysis._contact_result_matches_settings(
                contact_result,
                "deadbeef",
                tmp_path / "result.json",
            )
            is False
        )

    def test_contact_result_missing_required_fingerprint_is_rejected(self, tmp_path):
        from polyzymd.analyses.exposure import ExposureAnalysis

        contact_result = SimpleNamespace(metadata={})

        assert (
            ExposureAnalysis._contact_result_matches_settings(
                contact_result,
                "deadbeef",
                tmp_path / "result.json",
            )
            is False
        )

    def test_unproven_legacy_contact_artifact_is_not_resolved_when_fingerprint_required(
        self,
        tmp_path,
    ):
        from polyzymd.analyses.exposure import ExposureAnalysis

        contacts_dir = tmp_path / "contacts"
        contacts_dir.mkdir()
        (contacts_dir / "contacts_rep1.json").write_text("{}")

        result = ExposureAnalysis._find_contact_result(
            contacts_dir,
            1,
            settings_fp="deadbeef",
        )

        assert result is None

    def test_contact_artifact_path_rejects_filename_only_settings_token(self, tmp_path):
        """Exposure strict path matching should not trust ``_s`` filename tokens."""
        from polyzymd.analyses.exposure import ExposureAnalysis

        sidecar = tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json"
        sidecar.write_text("{}")

        assert not ExposureAnalysis._contact_artifact_path_matches_settings(sidecar, "deadbeef")

    def test_find_contact_result_rejects_filename_only_sidecar_without_requested_fp(
        self,
        tmp_path,
    ):
        """Strict downstream discovery must not fall back to filename-only sidecars."""
        from polyzymd.analyses.exposure import ExposureAnalysis

        contacts_dir = tmp_path / "contacts"
        contacts_dir.mkdir()
        (contacts_dir / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json").write_text("{}")

        result = ExposureAnalysis._find_contact_result(
            contacts_dir,
            1,
            equilibration="10ns",
        )

        assert result is None

    def test_contact_artifact_path_rejects_filename_only_sidecar_without_requested_fp(
        self,
        tmp_path,
    ):
        """Exposure strict matching rejects unproven fingerprinted sources."""
        from polyzymd.analyses.exposure import ExposureAnalysis

        sidecar = tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json"
        sidecar.write_text("{}")

        assert not ExposureAnalysis._contact_artifact_path_matches_settings(sidecar, None)

    def test_contact_result_rejects_filename_only_source_token(self, tmp_path):
        """Exposure strict result matching should require loaded metadata."""
        from polyzymd.analyses.exposure import ExposureAnalysis

        contact_result = SimpleNamespace(metadata={})
        source = tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json"

        assert not ExposureAnalysis._contact_result_matches_settings(
            contact_result,
            "deadbeef",
            source,
        )

    def test_contact_result_rejects_fingerprinted_source_without_requested_fp(self, tmp_path):
        """Exposure rejects fingerprinted sources when no embedded identity exists."""
        from polyzymd.analyses.exposure import ExposureAnalysis

        contact_result = SimpleNamespace(metadata={})
        source = tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json"
        source.write_text("{}")

        assert not ExposureAnalysis._contact_result_matches_settings(
            contact_result,
            None,
            source,
        )

    def test_mismatched_fingerprinted_sidecar_is_not_resolved(self, tmp_path):
        from polyzymd.analyses.exposure import ExposureAnalysis

        contacts_dir = tmp_path / "contacts"
        contacts_dir.mkdir()
        (contacts_dir / "contacts_eq10ns_cut4.5_sbadcafe_rep1.json").write_text("{}")

        result = ExposureAnalysis._find_contact_result(
            contacts_dir,
            1,
            settings_fp="deadbeef",
        )

        assert result is None

    def test_contact_result_matching_fingerprint_is_accepted(self, tmp_path):
        from polyzymd.analyses.exposure import ExposureAnalysis

        contact_result = SimpleNamespace(metadata={"settings_fingerprint": "deadbeef"})

        assert (
            ExposureAnalysis._contact_result_matches_settings(
                contact_result,
                "deadbeef",
                tmp_path / "result.json",
            )
            is True
        )

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

    def test_dynamics_cache_path_changes_with_contacts_identity(self, tmp_path):
        from polyzymd.analyses.exposure._config import ExposureConfig
        from polyzymd.analyses.exposure._dynamics import ExposureDynamicsResult
        from polyzymd.analyses.shared.config_hash import settings_fingerprint

        settings_fp = settings_fingerprint(ExposureConfig())

        path_a = ExposureDynamicsResult.cache_path(
            tmp_path,
            settings_fp=settings_fp,
            contacts_artifact_identity="aaaabbbbcccc",
        )
        path_b = ExposureDynamicsResult.cache_path(
            tmp_path,
            settings_fp=settings_fp,
            contacts_artifact_identity="dddd11112222",
        )

        assert path_a != path_b
        assert "_caaaabbbbcccc" in path_a.name

    def test_stale_dynamics_contacts_identity_fails_loud(self, tmp_path):
        from polyzymd.analyses.exposure import ExposureAnalysis
        from polyzymd.analyses.exposure._dynamics import ExposureDynamicsResult

        dynamics = ExposureDynamicsResult(contacts_artifact_identity="old")

        with pytest.raises(RuntimeError, match="contacts identity mismatch"):
            ExposureAnalysis._validate_dynamics_contacts_identity(
                dynamics,
                "new",
                tmp_path / "exposure_dynamics.json",
            )

    def test_analyze_dynamics_rejects_stale_contacts_identity_cache(self, tmp_path):
        from polyzymd.analyses.exposure._config import ExposureConfig
        from polyzymd.analyses.exposure._dynamics import (
            ExposureDynamicsResult,
            analyze_exposure_dynamics,
        )
        from polyzymd.analyses.shared.config_hash import settings_fingerprint

        config = ExposureConfig()
        cache_path = ExposureDynamicsResult.cache_path(
            tmp_path,
            settings_fp=settings_fingerprint(config),
            equilibration="0ns",
            contacts_artifact_identity="current",
        )
        ExposureDynamicsResult(contacts_artifact_identity="stale").save(cache_path)

        with pytest.raises(RuntimeError, match="contacts identity mismatch"):
            analyze_exposure_dynamics(
                sasa_result=MagicMock(),
                contact_result=MagicMock(),
                config=config,
                analysis_dir=tmp_path,
                recompute=False,
                equilibration="0ns",
                contacts_artifact_identity="current",
            )


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
        mock_traj.n_frames = 10
        mock_traj.n_atoms = 100
        mock_traj.topology = mock_topology
        mock_traj.__getitem__.return_value = mock_traj
        mock_load.return_value = mock_traj

        with pytest.raises(ValueError, match="Chain 'Z' not found in topology"):
            compute_trajectory_sasa(
                topology_path=tmp_path / "top.pdb",
                trajectory_path=tmp_path / "traj.xtc",
                config=SASAConfig(chain_id="Z"),
                analysis_dir=tmp_path / "analysis",
                recompute=True,
            )


# ===========================================================================
# Sibling SASA artifact reuse (Phase 2)
# ===========================================================================


class TestSiblingSASAReuse:
    """Tests for sibling SASA artifact lookup in compute_trajectory_sasa."""

    def test_sibling_hit_returns_adapted_result(self, tmp_path):
        """When a compatible sibling artifact exists, it should be loaded and adapted."""
        import json

        from polyzymd.analyses.exposure._sasa_config import SASAConfig
        from polyzymd.analyses.exposure._sasa_trajectory import compute_trajectory_sasa
        from polyzymd.analyses.shared.sasa import (
            SASA_ARTIFACT_COMPATIBILITY_VERSION,
            SASA_ARTIFACT_SCHEMA_NAME,
            SASA_ARTIFACT_SCHEMA_VERSION,
            compute_sasa_artifact_compatibility_hash,
        )

        exposure_run = tmp_path / "analysis" / "cond" / "exposure" / "run_1"
        exposure_run.mkdir(parents=True)
        sasa_run = tmp_path / "analysis" / "cond" / "sasa" / "run_1"
        sasa_run.mkdir(parents=True)

        n_frames = 10
        n_atoms = 50
        n_residues = 5
        protein_indices = np.arange(n_atoms, dtype=np.int64)

        compat_hash = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="protein and chainid 0",
            context_selection="protein and chainid 0",
            equilibration="0ns",
        )

        metadata = {
            "artifact_schema": SASA_ARTIFACT_SCHEMA_NAME,
            "artifact_schema_version": SASA_ARTIFACT_SCHEMA_VERSION,
            "artifact_compatibility_version": SASA_ARTIFACT_COMPATIBILITY_VERSION,
            "probe_radius_nm": 0.14,
            "n_sphere_points": 960,
            "equilibration": "0ns",
            "units": "A^2",
            "sasa_mode": "atom",
            "compatibility_hash": compat_hash,
            "target_selection": "protein and chainid 0",
            "context_selection": "protein and chainid 0",
        }
        (sasa_run / "sasa_protein_isolated.json").write_text(json.dumps(metadata))

        np.savez_compressed(
            sasa_run / "sasa_protein_isolated.npz",
            atom_sasa_a2=np.random.rand(n_frames, n_atoms).astype(np.float64) * 100,
            residue_sasa_a2=np.random.rand(n_frames, n_residues).astype(np.float64) * 100,
            total_sasa_a2=np.random.rand(n_frames).astype(np.float64) * 500,
            frames=np.arange(n_frames),
            time_ns=np.arange(n_frames, dtype=np.float64),
            target_atom_indices=protein_indices,
            context_atom_indices=protein_indices,
            residue_keys=np.array(["A:1:ALA", "A:2:GLY", "A:3:VAL", "A:4:LEU", "A:5:ILE"]),
            residue_chainids=np.array(["A"] * n_residues),
            residue_resids=np.arange(1, n_residues + 1),
            residue_resnames=np.array(["ALA", "GLY", "VAL", "LEU", "ILE"]),
        )

        config = SASAConfig(probe_radius_nm=0.14, n_sphere_points=960, chain_id="A")

        with patch("mdtraj.load_topology") as mock_load_topology:
            mock_topology = MagicMock()
            mock_topology.select.return_value = protein_indices
            mock_load_topology.return_value = mock_topology

            with patch("mdtraj.load") as mock_load_trajectory:
                result = compute_trajectory_sasa(
                    topology_path=tmp_path / "top.pdb",
                    trajectory_path=tmp_path / "traj.xtc",
                    config=config,
                    analysis_dir=exposure_run,
                    recompute=False,
                    equilibration="0ns",
                )

        mock_load_trajectory.assert_not_called()
        assert result is not None
        assert result.n_frames == n_frames
        assert result.n_residues == n_residues

    def test_sibling_index_mismatch_falls_through(self, tmp_path):
        """When sibling atom indices don't match, should return None."""
        import json

        from polyzymd.analyses.exposure._sasa_config import SASAConfig
        from polyzymd.analyses.exposure._sasa_trajectory import _try_load_sibling_sasa
        from polyzymd.analyses.shared.sasa import (
            SASA_ARTIFACT_COMPATIBILITY_VERSION,
            SASA_ARTIFACT_SCHEMA_NAME,
            SASA_ARTIFACT_SCHEMA_VERSION,
            compute_sasa_artifact_compatibility_hash,
        )

        exposure_run = tmp_path / "analysis" / "cond" / "exposure" / "run_1"
        exposure_run.mkdir(parents=True)
        sasa_run = tmp_path / "analysis" / "cond" / "sasa" / "run_1"
        sasa_run.mkdir(parents=True)

        n_frames, n_residues = 10, 5
        sibling_target = np.arange(200, dtype=np.int64)
        sibling_context = np.arange(200, dtype=np.int64)

        compat_hash = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="protein and chainid 0",
            context_selection="protein and chainid 0",
            equilibration="0ns",
        )

        metadata = {
            "artifact_schema": SASA_ARTIFACT_SCHEMA_NAME,
            "artifact_schema_version": SASA_ARTIFACT_SCHEMA_VERSION,
            "artifact_compatibility_version": SASA_ARTIFACT_COMPATIBILITY_VERSION,
            "probe_radius_nm": 0.14,
            "n_sphere_points": 960,
            "equilibration": "0ns",
            "units": "A^2",
            "sasa_mode": "atom",
            "compatibility_hash": compat_hash,
            "target_selection": "protein and chainid 0",
            "context_selection": "protein and chainid 0",
        }
        (sasa_run / "sasa_protein_isolated.json").write_text(json.dumps(metadata))

        np.savez_compressed(
            sasa_run / "sasa_protein_isolated.npz",
            atom_sasa_a2=np.random.rand(n_frames, 200).astype(np.float64) * 100,
            residue_sasa_a2=np.random.rand(n_frames, n_residues).astype(np.float64) * 100,
            total_sasa_a2=np.random.rand(n_frames).astype(np.float64) * 500,
            frames=np.arange(n_frames),
            time_ns=np.arange(n_frames, dtype=np.float64),
            target_atom_indices=sibling_target,
            context_atom_indices=sibling_context,
            residue_keys=np.array(["A:1:ALA", "A:2:GLY", "A:3:VAL", "A:4:LEU", "A:5:ILE"]),
            residue_chainids=np.array(["A"] * n_residues),
            residue_resids=np.arange(1, n_residues + 1),
            residue_resnames=np.array(["ALA", "GLY", "VAL", "LEU", "ILE"]),
        )

        config = SASAConfig(probe_radius_nm=0.14, n_sphere_points=960, chain_id="A")
        expected_indices = np.arange(50, dtype=np.int64)

        with patch(
            "polyzymd.analyses.exposure._sasa_trajectory._resolve_protein_indices_from_topology",
            return_value=expected_indices,
        ):
            result = _try_load_sibling_sasa(
                topology_path=tmp_path / "top.pdb",
                config=config,
                analysis_dir=exposure_run,
                equilibration="0ns",
            )

        assert result is None

    def test_legacy_sibling_without_atom_indices_returns_none(self, tmp_path):
        """Legacy sibling artifacts without index arrays should be skipped."""
        import json

        from polyzymd.analyses.exposure._sasa_config import SASAConfig
        from polyzymd.analyses.exposure._sasa_trajectory import _try_load_sibling_sasa
        from polyzymd.analyses.shared.sasa import (
            SASA_ARTIFACT_COMPATIBILITY_VERSION,
            SASA_ARTIFACT_SCHEMA_NAME,
            SASA_ARTIFACT_SCHEMA_VERSION,
            compute_sasa_artifact_compatibility_hash,
        )

        exposure_run = tmp_path / "analysis" / "cond" / "exposure" / "run_1"
        exposure_run.mkdir(parents=True)
        sasa_run = tmp_path / "analysis" / "cond" / "sasa" / "run_1"
        sasa_run.mkdir(parents=True)

        n_frames = 3
        n_residues = 5

        compat_hash = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="protein and chainid 0",
            context_selection="protein and chainid 0",
            equilibration="0ns",
        )

        metadata = {
            "artifact_schema": SASA_ARTIFACT_SCHEMA_NAME,
            "artifact_schema_version": SASA_ARTIFACT_SCHEMA_VERSION,
            "artifact_compatibility_version": SASA_ARTIFACT_COMPATIBILITY_VERSION,
            "probe_radius_nm": 0.14,
            "n_sphere_points": 960,
            "equilibration": "0ns",
            "units": "A^2",
            "sasa_mode": "atom",
            "compatibility_hash": compat_hash,
            "target_selection": "protein and chainid 0",
            "context_selection": "protein and chainid 0",
        }
        (sasa_run / "sasa_legacy.json").write_text(json.dumps(metadata))

        np.savez_compressed(
            sasa_run / "sasa_legacy.npz",
            atom_sasa_a2=np.random.rand(n_frames, 50).astype(np.float64) * 100,
            residue_sasa_a2=np.random.rand(n_frames, n_residues).astype(np.float64) * 100,
            total_sasa_a2=np.random.rand(n_frames).astype(np.float64) * 500,
            frames=np.arange(n_frames),
            time_ns=np.arange(n_frames, dtype=np.float64),
            residue_keys=np.array(["A:1:ALA", "A:2:GLY", "A:3:VAL", "A:4:LEU", "A:5:ILE"]),
            residue_chainids=np.array(["A"] * n_residues),
            residue_resids=np.arange(1, n_residues + 1),
            residue_resnames=np.array(["ALA", "GLY", "VAL", "LEU", "ILE"]),
        )

        config = SASAConfig(probe_radius_nm=0.14, n_sphere_points=960, chain_id="A")

        with patch(
            "polyzymd.analyses.exposure._sasa_trajectory._resolve_protein_indices_from_topology",
            return_value=np.arange(50, dtype=np.int64),
        ):
            result = _try_load_sibling_sasa(
                topology_path=tmp_path / "top.pdb",
                config=config,
                analysis_dir=exposure_run,
                equilibration="0ns",
            )

        assert result is None

    def test_sibling_context_with_polymer_rejected_even_when_target_matches(self, tmp_path):
        """Sibling artifacts with polymer context should be rejected for exposure."""
        import json

        from polyzymd.analyses.exposure._sasa_config import SASAConfig
        from polyzymd.analyses.exposure._sasa_trajectory import _try_load_sibling_sasa
        from polyzymd.analyses.shared.sasa import (
            SASA_ARTIFACT_COMPATIBILITY_VERSION,
            SASA_ARTIFACT_SCHEMA_NAME,
            SASA_ARTIFACT_SCHEMA_VERSION,
            compute_sasa_artifact_compatibility_hash,
        )

        exposure_run = tmp_path / "analysis" / "cond" / "exposure" / "run_1"
        exposure_run.mkdir(parents=True)
        sasa_run = tmp_path / "analysis" / "cond" / "sasa" / "run_1"
        sasa_run.mkdir(parents=True)

        n_frames = 3
        n_residues = 5
        protein_indices = np.arange(50, dtype=np.int64)
        protein_plus_polymer_indices = np.arange(100, dtype=np.int64)

        compat_hash = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="protein and chainid 0",
            context_selection="protein and chainid 0",
            equilibration="0ns",
        )

        metadata = {
            "artifact_schema": SASA_ARTIFACT_SCHEMA_NAME,
            "artifact_schema_version": SASA_ARTIFACT_SCHEMA_VERSION,
            "artifact_compatibility_version": SASA_ARTIFACT_COMPATIBILITY_VERSION,
            "probe_radius_nm": 0.14,
            "n_sphere_points": 960,
            "equilibration": "0ns",
            "units": "A^2",
            "sasa_mode": "atom",
            "compatibility_hash": compat_hash,
            "target_selection": "protein and chainid 0",
            "context_selection": "protein and chainid 0",
        }
        (sasa_run / "sasa_protein_with_polymer_context.json").write_text(json.dumps(metadata))

        np.savez_compressed(
            sasa_run / "sasa_protein_with_polymer_context.npz",
            atom_sasa_a2=np.random.rand(n_frames, 50).astype(np.float64) * 100,
            residue_sasa_a2=np.random.rand(n_frames, n_residues).astype(np.float64) * 100,
            total_sasa_a2=np.random.rand(n_frames).astype(np.float64) * 500,
            frames=np.arange(n_frames),
            time_ns=np.arange(n_frames, dtype=np.float64),
            target_atom_indices=protein_indices,
            context_atom_indices=protein_plus_polymer_indices,
            residue_keys=np.array(["A:1:ALA", "A:2:GLY", "A:3:VAL", "A:4:LEU", "A:5:ILE"]),
            residue_chainids=np.array(["A"] * n_residues),
            residue_resids=np.arange(1, n_residues + 1),
            residue_resnames=np.array(["ALA", "GLY", "VAL", "LEU", "ILE"]),
        )

        config = SASAConfig(probe_radius_nm=0.14, n_sphere_points=960, chain_id="A")

        with patch(
            "polyzymd.analyses.exposure._sasa_trajectory._resolve_protein_indices_from_topology",
            return_value=protein_indices,
        ):
            result = _try_load_sibling_sasa(
                topology_path=tmp_path / "top.pdb",
                config=config,
                analysis_dir=exposure_run,
                equilibration="0ns",
            )

        assert result is None

    def test_sibling_hit_converts_residue_sasa_a2_to_nm2(self, tmp_path):
        """Sibling adaptation should convert residue SASA from Å² to nm²."""
        import json

        from polyzymd.analyses.exposure._sasa_config import SASAConfig
        from polyzymd.analyses.exposure._sasa_trajectory import _try_load_sibling_sasa
        from polyzymd.analyses.shared.sasa import (
            A2_TO_NM2,
            SASA_ARTIFACT_COMPATIBILITY_VERSION,
            SASA_ARTIFACT_SCHEMA_NAME,
            SASA_ARTIFACT_SCHEMA_VERSION,
            compute_sasa_artifact_compatibility_hash,
        )

        exposure_run = tmp_path / "analysis" / "cond" / "exposure" / "run_1"
        exposure_run.mkdir(parents=True)
        sasa_run = tmp_path / "analysis" / "cond" / "sasa" / "run_1"
        sasa_run.mkdir(parents=True)

        protein_indices = np.arange(50, dtype=np.int64)
        residue_sasa_a2 = np.asarray([[100.0, 200.0, 300.0]], dtype=np.float64)

        compat_hash = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="protein and chainid 0",
            context_selection="protein and chainid 0",
            equilibration="0ns",
        )

        metadata = {
            "artifact_schema": SASA_ARTIFACT_SCHEMA_NAME,
            "artifact_schema_version": SASA_ARTIFACT_SCHEMA_VERSION,
            "artifact_compatibility_version": SASA_ARTIFACT_COMPATIBILITY_VERSION,
            "probe_radius_nm": 0.14,
            "n_sphere_points": 960,
            "equilibration": "0ns",
            "units": "A^2",
            "sasa_mode": "atom",
            "compatibility_hash": compat_hash,
            "target_selection": "protein and chainid 0",
            "context_selection": "protein and chainid 0",
        }
        (sasa_run / "sasa_unit_conversion.json").write_text(json.dumps(metadata))

        np.savez_compressed(
            sasa_run / "sasa_unit_conversion.npz",
            atom_sasa_a2=np.random.rand(1, 50).astype(np.float64) * 100,
            residue_sasa_a2=residue_sasa_a2,
            total_sasa_a2=np.asarray([600.0], dtype=np.float64),
            frames=np.asarray([0], dtype=np.int64),
            time_ns=np.asarray([0.0], dtype=np.float64),
            target_atom_indices=protein_indices,
            context_atom_indices=protein_indices,
            residue_keys=np.array(["A:1:ALA", "A:2:GLY", "A:3:VAL"]),
            residue_chainids=np.array(["A", "A", "A"]),
            residue_resids=np.asarray([1, 2, 3], dtype=np.int64),
            residue_resnames=np.array(["ALA", "GLY", "VAL"]),
        )

        config = SASAConfig(probe_radius_nm=0.14, n_sphere_points=960, chain_id="A")

        with patch(
            "polyzymd.analyses.exposure._sasa_trajectory._resolve_protein_indices_from_topology",
            return_value=protein_indices,
        ):
            result = _try_load_sibling_sasa(
                topology_path=tmp_path / "top.pdb",
                config=config,
                analysis_dir=exposure_run,
                equilibration="0ns",
            )

        assert result is not None
        expected_nm2 = residue_sasa_a2.astype(np.float32) * A2_TO_NM2
        np.testing.assert_allclose(result.sasa_per_frame, expected_nm2)

    def test_sibling_missing_probe_radius_metadata_returns_none(self, tmp_path):
        """Sibling metadata missing required fields should be rejected."""
        import json

        from polyzymd.analyses.exposure._sasa_config import SASAConfig
        from polyzymd.analyses.exposure._sasa_trajectory import _try_load_sibling_sasa
        from polyzymd.analyses.shared.sasa import (
            SASA_ARTIFACT_COMPATIBILITY_VERSION,
            SASA_ARTIFACT_SCHEMA_NAME,
            SASA_ARTIFACT_SCHEMA_VERSION,
            compute_sasa_artifact_compatibility_hash,
        )

        exposure_run = tmp_path / "analysis" / "cond" / "exposure" / "run_1"
        exposure_run.mkdir(parents=True)
        sasa_run = tmp_path / "analysis" / "cond" / "sasa" / "run_1"
        sasa_run.mkdir(parents=True)

        compat_hash = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="protein and chainid 0",
            context_selection="protein and chainid 0",
            equilibration="0ns",
        )

        metadata = {
            "artifact_schema": SASA_ARTIFACT_SCHEMA_NAME,
            "artifact_schema_version": SASA_ARTIFACT_SCHEMA_VERSION,
            "artifact_compatibility_version": SASA_ARTIFACT_COMPATIBILITY_VERSION,
            "n_sphere_points": 960,
            "equilibration": "0ns",
            "units": "A^2",
            "sasa_mode": "atom",
            "compatibility_hash": compat_hash,
            "target_selection": "protein and chainid 0",
            "context_selection": "protein and chainid 0",
        }
        (sasa_run / "sasa_missing_probe.json").write_text(json.dumps(metadata))

        protein_indices = np.arange(50, dtype=np.int64)
        np.savez_compressed(
            sasa_run / "sasa_missing_probe.npz",
            atom_sasa_a2=np.random.rand(1, 50).astype(np.float64) * 100,
            residue_sasa_a2=np.random.rand(1, 3).astype(np.float64) * 100,
            total_sasa_a2=np.asarray([250.0], dtype=np.float64),
            frames=np.asarray([0], dtype=np.int64),
            time_ns=np.asarray([0.0], dtype=np.float64),
            target_atom_indices=protein_indices,
            context_atom_indices=protein_indices,
            residue_keys=np.array(["A:1:ALA", "A:2:GLY", "A:3:VAL"]),
            residue_chainids=np.array(["A", "A", "A"]),
            residue_resids=np.asarray([1, 2, 3], dtype=np.int64),
            residue_resnames=np.array(["ALA", "GLY", "VAL"]),
        )

        config = SASAConfig(probe_radius_nm=0.14, n_sphere_points=960, chain_id="A")
        result = _try_load_sibling_sasa(
            topology_path=tmp_path / "top.pdb",
            config=config,
            analysis_dir=exposure_run,
            equilibration="0ns",
        )

        assert result is None

    def test_corrupted_sibling_metadata_json_skipped(self, tmp_path):
        """Sibling metadata JSON parse failures should be skipped safely."""
        from polyzymd.analyses.exposure._sasa_config import SASAConfig
        from polyzymd.analyses.exposure._sasa_trajectory import _try_load_sibling_sasa

        exposure_run = tmp_path / "analysis" / "cond" / "exposure" / "run_1"
        exposure_run.mkdir(parents=True)
        sasa_run = tmp_path / "analysis" / "cond" / "sasa" / "run_1"
        sasa_run.mkdir(parents=True)

        np.savez_compressed(
            sasa_run / "sasa_bad_metadata.npz",
            atom_sasa_a2=np.random.rand(1, 50).astype(np.float64) * 100,
            residue_sasa_a2=np.random.rand(1, 3).astype(np.float64) * 100,
            total_sasa_a2=np.asarray([250.0], dtype=np.float64),
            frames=np.asarray([0], dtype=np.int64),
            time_ns=np.asarray([0.0], dtype=np.float64),
            target_atom_indices=np.arange(50, dtype=np.int64),
            context_atom_indices=np.arange(50, dtype=np.int64),
            residue_keys=np.array(["A:1:ALA", "A:2:GLY", "A:3:VAL"]),
            residue_chainids=np.array(["A", "A", "A"]),
            residue_resids=np.asarray([1, 2, 3], dtype=np.int64),
            residue_resnames=np.array(["ALA", "GLY", "VAL"]),
        )
        (sasa_run / "sasa_bad_metadata.json").write_bytes(b"not valid json {{{")

        config = SASAConfig(probe_radius_nm=0.14, n_sphere_points=960, chain_id="A")
        result = _try_load_sibling_sasa(
            topology_path=tmp_path / "top.pdb",
            config=config,
            analysis_dir=exposure_run,
            equilibration="0ns",
        )

        assert result is None

    def test_corrupted_sibling_skipped(self, tmp_path):
        """Corrupted sibling NPZ should be caught and skipped."""
        import json

        from polyzymd.analyses.exposure._sasa_config import SASAConfig
        from polyzymd.analyses.exposure._sasa_trajectory import _try_load_sibling_sasa
        from polyzymd.analyses.shared.sasa import (
            SASA_ARTIFACT_COMPATIBILITY_VERSION,
            SASA_ARTIFACT_SCHEMA_NAME,
            SASA_ARTIFACT_SCHEMA_VERSION,
        )

        exposure_run = tmp_path / "analysis" / "cond" / "exposure" / "run_1"
        exposure_run.mkdir(parents=True)
        sasa_run = tmp_path / "analysis" / "cond" / "sasa" / "run_1"
        sasa_run.mkdir(parents=True)

        metadata = {
            "artifact_schema": SASA_ARTIFACT_SCHEMA_NAME,
            "artifact_schema_version": SASA_ARTIFACT_SCHEMA_VERSION,
            "artifact_compatibility_version": SASA_ARTIFACT_COMPATIBILITY_VERSION,
            "probe_radius_nm": 0.14,
            "n_sphere_points": 960,
            "equilibration": "0ns",
            "units": "A^2",
            "sasa_mode": "atom",
        }
        (sasa_run / "sasa_corrupted.json").write_text(json.dumps(metadata))
        (sasa_run / "sasa_corrupted.npz").write_bytes(b"not a valid npz file")

        config = SASAConfig(probe_radius_nm=0.14, n_sphere_points=960, chain_id="A")

        with patch(
            "polyzymd.analyses.exposure._sasa_trajectory._resolve_protein_indices_from_topology",
            return_value=np.arange(50, dtype=np.int64),
        ):
            result = _try_load_sibling_sasa(
                topology_path=tmp_path / "top.pdb",
                config=config,
                analysis_dir=exposure_run,
                equilibration="0ns",
            )

        assert result is None

    def test_no_siblings_returns_none(self, tmp_path):
        """When no sibling sasa directory exists, should return None."""
        from polyzymd.analyses.exposure._sasa_config import SASAConfig
        from polyzymd.analyses.exposure._sasa_trajectory import _try_load_sibling_sasa

        exposure_run = tmp_path / "analysis" / "cond" / "exposure" / "run_1"
        exposure_run.mkdir(parents=True)

        config = SASAConfig(probe_radius_nm=0.14, n_sphere_points=960, chain_id="A")

        result = _try_load_sibling_sasa(
            topology_path=tmp_path / "top.pdb",
            config=config,
            analysis_dir=exposure_run,
            equilibration="0ns",
        )

        assert result is None

    def test_recompute_bypasses_sibling_lookup(self, tmp_path):
        """When recompute=True, sibling lookup should be skipped entirely."""
        from polyzymd.analyses.exposure._sasa_config import SASAConfig
        from polyzymd.analyses.exposure._sasa_trajectory import compute_trajectory_sasa

        exposure_run = tmp_path / "analysis" / "cond" / "exposure" / "run_1"
        exposure_run.mkdir(parents=True)

        config = SASAConfig(probe_radius_nm=0.14, n_sphere_points=960, chain_id="A")

        with (
            patch(
                "polyzymd.analyses.exposure._sasa_trajectory._try_load_sibling_sasa"
            ) as mock_sibling,
            patch("mdtraj.load") as mock_md_load,
            patch("mdtraj.shrake_rupley") as mock_shrake_rupley,
        ):
            mock_topology = MagicMock()
            protein_indices = np.arange(50, dtype=np.int64)
            mock_topology.select.return_value = protein_indices

            residues = []
            residue_names = ["ALA", "GLY", "VAL", "LEU", "ILE"]
            for i, name in enumerate(residue_names, start=1):
                residue = MagicMock()
                residue.resSeq = i
                residue.name = name
                residues.append(residue)

            mock_protein_topology = MagicMock()
            mock_protein_topology.residues = residues

            mock_protein_traj = MagicMock()
            mock_protein_traj.n_atoms = 50
            mock_protein_traj.n_residues = 5
            mock_protein_traj.topology = mock_protein_topology

            mock_traj = MagicMock()
            mock_traj.n_frames = 10
            mock_traj.n_atoms = 100
            mock_traj.topology = mock_topology
            mock_traj.__getitem__.return_value = mock_traj
            mock_traj.atom_slice.return_value = mock_protein_traj
            mock_md_load.return_value = mock_traj

            mock_shrake_rupley.return_value = np.random.rand(10, 5).astype(np.float32)

            compute_trajectory_sasa(
                topology_path=tmp_path / "top.pdb",
                trajectory_path=tmp_path / "traj.xtc",
                config=config,
                analysis_dir=exposure_run,
                recompute=True,
                equilibration="0ns",
            )

        mock_sibling.assert_not_called()

    def test_equilibration_propagated_to_call_sites(self, tmp_path):
        """_load_or_compute_replicate should pass equilibration to compute_trajectory_sasa."""
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

        exposure_dir = tmp_path / "exposure"
        contacts_dir = tmp_path / "contacts"
        contacts_dir.mkdir(parents=True)
        contact_file = contacts_dir / "contacts_rep1.json"
        contact_file.write_text("{}")

        with (
            patch.object(analysis, "_find_contact_result", return_value=contact_file),
            patch("polyzymd.analyses.contacts._results.ContactResult") as mock_contact_result,
            patch("polyzymd.analyses.shared.loader.TrajectoryLoader") as mock_loader_cls,
            patch(
                "polyzymd.analyses.exposure._sasa_trajectory.compute_trajectory_sasa"
            ) as mock_sasa,
            patch(
                "polyzymd.analyses.exposure._dynamics.analyze_exposure_dynamics"
            ) as mock_dynamics,
            patch(
                "polyzymd.analyses.exposure._enrichment.compute_chaperone_enrichment"
            ) as mock_enrichment,
        ):
            loaded_contacts = MagicMock()
            loaded_contacts.equilibration_time = 10.0
            loaded_contacts.equilibration_unit = "ns"
            loaded_contacts.replicate = 1
            loaded_contacts.n_frames = 5
            mock_contact_result.load.return_value = loaded_contacts
            mock_loader = MagicMock()
            traj_info = MagicMock()
            traj_info.topology_file = tmp_path / "top.pdb"
            traj_info.trajectory_files = [tmp_path / "traj.xtc"]
            mock_loader.get_trajectory_info.return_value = traj_info
            mock_loader_cls.return_value = mock_loader

            mock_sasa_result = MagicMock()
            mock_sasa_result.n_frames = 5
            mock_sasa.return_value = mock_sasa_result
            mock_dynamics.return_value = MagicMock()
            mock_enrichment.return_value = MagicMock()

            analysis._load_or_compute_replicate(
                cond=cond,
                replicate=1,
                settings=ExposureAnalysis.Settings(),
                exposure_analysis_dir=exposure_dir,
                contacts_analysis_dir=contacts_dir,
                contacts_settings_fp=None,
                recompute=True,
                equilibration="10ns",
            )

        assert mock_sasa.called
        assert all(call.kwargs.get("equilibration") == "10ns" for call in mock_sasa.call_args_list)

    def test_adaptation_failure_skips_sibling(self, tmp_path):
        """If adapt_canonical_sasa_to_exposure() raises, the sibling is skipped."""
        import json

        import numpy as np

        from polyzymd.analyses.exposure._sasa_config import SASAConfig
        from polyzymd.analyses.exposure._sasa_trajectory import _try_load_sibling_sasa
        from polyzymd.analyses.shared.sasa import (
            SASA_ARTIFACT_COMPATIBILITY_VERSION,
            SASA_ARTIFACT_SCHEMA_NAME,
            SASA_ARTIFACT_SCHEMA_VERSION,
        )

        exposure_run = tmp_path / "analysis" / "cond" / "exposure" / "run_1"
        exposure_run.mkdir(parents=True)
        sasa_run = tmp_path / "analysis" / "cond" / "sasa" / "run_1"
        sasa_run.mkdir(parents=True)

        n_atoms = 50
        protein_indices = np.arange(n_atoms, dtype=np.int64)

        metadata = {
            "artifact_schema": SASA_ARTIFACT_SCHEMA_NAME,
            "artifact_schema_version": SASA_ARTIFACT_SCHEMA_VERSION,
            "artifact_compatibility_version": SASA_ARTIFACT_COMPATIBILITY_VERSION,
            "probe_radius_nm": 0.14,
            "n_sphere_points": 960,
            "equilibration": "0ns",
            "units": "A^2",
            "sasa_mode": "atom",
        }
        (sasa_run / "sasa_protein_isolated.json").write_text(json.dumps(metadata))

        # NPZ is loadable but has mismatched residue arrays that will break adaptation
        np.savez_compressed(
            sasa_run / "sasa_protein_isolated.npz",
            sasa_a2=np.random.rand(10, n_atoms).astype(np.float32),
            residue_sasa_a2=np.random.rand(10, 5).astype(np.float32),
            frames=np.arange(10),
            time_ns=np.arange(10, dtype=np.float64),
            target_atom_indices=protein_indices,
            context_atom_indices=protein_indices,
            residue_keys=np.array(["A:1:ALA", "A:2:GLY", "A:3:VAL", "A:4:LEU", "A:5:ILE"]),
            residue_chainids=np.array(["A"] * 5),
            residue_resids=np.arange(1, 6),
            residue_resnames=np.array(["ALA", "GLY", "VAL", "LEU", "ILE"]),
        )

        config = SASAConfig(probe_radius_nm=0.14, n_sphere_points=960, chain_id="A")

        with (
            patch(
                "polyzymd.analyses.exposure._sasa_trajectory._resolve_protein_indices_from_topology",
                return_value=protein_indices,
            ),
            patch(
                "polyzymd.analyses.shared.sasa.adapt_canonical_sasa_to_exposure",
                side_effect=ValueError("shape mismatch during adaptation"),
            ),
        ):
            result = _try_load_sibling_sasa(
                topology_path=tmp_path / "top.pdb",
                config=config,
                analysis_dir=exposure_run,
                equilibration="0ns",
            )

        # Should return None because adaptation failed
        assert result is None


# ===========================================================================
# Equilibration-aware cache identity (Phase 3)
# ===========================================================================


class TestEquilibrationCacheIdentity:
    """Cache paths must differ when equilibration differs."""

    def test_different_equilibration_gives_different_cache_path(self, tmp_path):
        from polyzymd.analyses.exposure._sasa_config import SASAConfig
        from polyzymd.analyses.exposure._sasa_trajectory import SASATrajectoryResult
        from polyzymd.analyses.shared.config_hash import settings_fingerprint

        cfg = SASAConfig(probe_radius_nm=0.14, n_sphere_points=960)
        fp = settings_fingerprint(cfg)

        path_0ns = SASATrajectoryResult.cache_path(tmp_path, settings_fp=fp, equilibration="0ns")
        path_10ns = SASATrajectoryResult.cache_path(tmp_path, settings_fp=fp, equilibration="10ns")
        path_50ns = SASATrajectoryResult.cache_path(tmp_path, settings_fp=fp, equilibration="50ns")

        assert path_0ns != path_10ns
        assert path_0ns != path_50ns
        assert path_10ns != path_50ns

    def test_equilibration_none_omits_eq_segment(self, tmp_path):
        """When equilibration is None (legacy callers), no eq_ segment appears."""
        from polyzymd.analyses.exposure._sasa_trajectory import SASATrajectoryResult

        path_none = SASATrajectoryResult.cache_path(tmp_path, settings_fp="abc", equilibration=None)
        path_0ns = SASATrajectoryResult.cache_path(tmp_path, settings_fp="abc", equilibration="0ns")

        assert "eq_" not in str(path_none)
        assert "eq_0ns" in str(path_0ns)
        assert path_none != path_0ns

    def test_whitespace_normalized(self, tmp_path):
        """Whitespace in equilibration labels should be normalized."""
        from polyzymd.analyses.exposure._sasa_trajectory import SASATrajectoryResult

        path_clean = SASATrajectoryResult.cache_path(
            tmp_path, settings_fp="abc", equilibration="10ns"
        )
        path_padded = SASATrajectoryResult.cache_path(
            tmp_path, settings_fp="abc", equilibration="  10ns  "
        )

        assert path_clean == path_padded

    def test_cache_path_includes_eq_and_fp_segments(self, tmp_path):
        from polyzymd.analyses.exposure._sasa_trajectory import SASATrajectoryResult

        path = SASATrajectoryResult.cache_path(
            tmp_path, settings_fp="deadbeef", equilibration="10ns"
        )
        parts = path.parts
        # Should contain sasa/eq_10ns/fp_deadbeef
        assert "sasa" in parts
        assert "eq_10ns" in parts
        assert "fp_deadbeef" in parts

    def test_different_equilibration_not_cached_together(self, tmp_path):
        """A cached SASA for eq=0ns must not be returned for eq=10ns."""
        from polyzymd.analyses.exposure._sasa_config import SASAConfig
        from polyzymd.analyses.exposure._sasa_trajectory import (
            SASATrajectoryResult,
            compute_trajectory_sasa,
        )

        config = SASAConfig(probe_radius_nm=0.14, n_sphere_points=960, chain_id="A")
        analysis_dir = tmp_path / "analysis" / "run_1"
        analysis_dir.mkdir(parents=True)

        # Pre-populate a cache for equilibration="0ns"
        from polyzymd.analyses.shared.config_hash import settings_fingerprint

        fp = settings_fingerprint(config)
        cache_0ns = SASATrajectoryResult.cache_path(
            analysis_dir, settings_fp=fp, equilibration="0ns"
        )
        cache_0ns.mkdir(parents=True)

        import numpy as np

        np.savez_compressed(
            cache_0ns / "sasa_trajectory.npz",
            sasa_per_frame=np.ones((5, 3), dtype=np.float32),
            relative_sasa_per_frame=np.ones((5, 3), dtype=np.float32) * 0.5,
            resids=np.array([1, 2, 3], dtype=np.int32),
            max_sasa_nm2=np.array([2.0, 2.0, 2.0], dtype=np.float32),
        )
        import json

        (cache_0ns / "sasa_metadata.json").write_text(
            json.dumps(
                {
                    "resnames": ["ALA", "GLY", "VAL"],
                    "aa_classes": ["nonpolar", "nonpolar", "nonpolar"],
                    "n_frames": 5,
                    "n_residues": 3,
                    "exposure_threshold": 0.2,
                    "trajectory_path": "",
                    "topology_path": "",
                    "equilibration": "0ns",
                }
            )
        )

        # Loading with eq="0ns" should find the cache (no mdtraj needed)
        result_0ns = compute_trajectory_sasa(
            topology_path=tmp_path / "top.pdb",
            trajectory_path=tmp_path / "traj.xtc",
            config=config,
            analysis_dir=analysis_dir,
            recompute=False,
            equilibration="0ns",
        )
        assert result_0ns.n_frames == 5

        # Loading with eq="10ns" should NOT find the cache (different path)
        # and should attempt to compute (hitting mdtraj.load)
        mock_traj = MagicMock()
        mock_traj.n_frames = 20
        mock_traj.n_atoms = 100
        mock_traj.timestep = 1000.0
        mock_topology = MagicMock()
        mock_topology.select.return_value = np.arange(50)
        mock_window_traj = MagicMock()
        mock_window_traj.n_frames = 10
        mock_window_traj.n_atoms = 100
        mock_window_traj.topology = mock_topology
        mock_sub_traj = MagicMock()
        mock_sub_traj.n_atoms = 50
        mock_sub_traj.n_residues = 3
        mock_residues = []
        for i in range(1, 4):
            r = MagicMock()
            r.resSeq = i
            r.name = "ALA"
            mock_residues.append(r)
        mock_sub_traj.topology.residues = mock_residues
        mock_traj.topology = mock_topology
        mock_traj.__getitem__.return_value = mock_window_traj
        mock_window_traj.atom_slice.return_value = mock_sub_traj

        with (
            patch("mdtraj.load", return_value=mock_traj) as mock_load,
            patch(
                "mdtraj.shrake_rupley",
                return_value=np.random.rand(10, 3).astype(np.float32),
            ),
        ):
            result_10ns = compute_trajectory_sasa(
                topology_path=tmp_path / "top.pdb",
                trajectory_path=tmp_path / "traj.xtc",
                config=config,
                analysis_dir=analysis_dir,
                recompute=False,
                equilibration="10ns",
            )
            # mdtraj.load SHOULD have been called (no cache for 10ns)
            mock_load.assert_called_once()
            mock_traj.__getitem__.assert_called_once()
            assert result_10ns.n_frames == 10
