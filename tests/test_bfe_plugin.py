"""Tests for BindingFreeEnergyAnalysis plugin — Phase C-6.

Coverage:
- Discovery (registered by name and alias)
- Class attributes (name, Settings, aliases, dependencies, min_replicates)
- Settings validation (units, threshold, k_b)
- compute_replicate / aggregate (no-op)
- filter_conditions (keeps all)
- compare (full override with temperature-aware pairwise)
- pairwise comparisons (cross-temperature suppression)
- plot delegation
- extract_metrics (empty)
- lifecycle
- ΔG_sel computation helpers
"""

from __future__ import annotations

import math
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# ===========================================================================
# Discovery
# ===========================================================================


class TestDiscovery:
    """BFE plugin is discovered by pkgutil scanning."""

    def test_discovery_by_name(self):
        from polyzymd.analyses import get_analysis

        cls = get_analysis("binding_free_energy")
        assert cls.__name__ == "BindingFreeEnergyAnalysis"

    def test_discovery_by_alias(self):
        from polyzymd.analyses import get_analysis

        cls = get_analysis("bfe")
        assert cls.__name__ == "BindingFreeEnergyAnalysis"

    def test_listed(self):
        from polyzymd.analyses import list_analyses

        analyses = list_analyses()
        assert "binding_free_energy" in analyses

    def test_all_names_includes_alias(self):
        from polyzymd.analyses import list_all_names

        names = list_all_names()
        assert "binding_free_energy" in names
        assert "bfe" in names


# ===========================================================================
# Class attributes
# ===========================================================================


class TestClassAttributes:
    def test_name(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        assert BindingFreeEnergyAnalysis.name == "binding_free_energy"

    def test_settings_type(self):
        from polyzymd.analyses.binding_free_energy import (
            BFESettings,
            BindingFreeEnergyAnalysis,
        )

        assert BindingFreeEnergyAnalysis.Settings is BFESettings

    def test_aliases(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        assert BindingFreeEnergyAnalysis.aliases == ("bfe",)

    def test_dependencies(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        assert BindingFreeEnergyAnalysis.dependencies == ("contacts",)

    def test_min_replicates(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        assert BindingFreeEnergyAnalysis.min_replicates == 1


# ===========================================================================
# Settings
# ===========================================================================


class TestSettings:
    def test_defaults(self):
        from polyzymd.analyses.binding_free_energy import BFESettings

        s = BFESettings()
        assert s.units == "kT"
        assert s.compute_binding_preference is True
        assert s.surface_exposure_threshold == pytest.approx(0.2)
        assert s.enzyme_pdb_for_sasa is None
        assert s.include_default_aa_groups is True
        assert s.protein_groups is None
        assert s.protein_partitions is None
        assert s.polymer_type_selections is None
        assert s.fdr_alpha == pytest.approx(0.05)

    def test_units_kcal(self):
        from polyzymd.analyses.binding_free_energy import BFESettings

        s = BFESettings(units="kcal/mol")
        assert s.units == "kcal/mol"

    def test_units_kj(self):
        from polyzymd.analyses.binding_free_energy import BFESettings

        s = BFESettings(units="kJ/mol")
        assert s.units == "kJ/mol"

    def test_invalid_units(self):
        from polyzymd.analyses.binding_free_energy import BFESettings

        with pytest.raises(Exception):
            BFESettings(units="eV")

    def test_k_b_kT(self):
        from polyzymd.analyses.binding_free_energy import BFESettings

        s = BFESettings(units="kT")
        assert s.k_b() == 0.0

    def test_k_b_kcal(self):
        from polyzymd.analyses.binding_free_energy import BFESettings

        s = BFESettings(units="kcal/mol")
        assert s.k_b() == pytest.approx(0.0019872041)

    def test_k_b_kj(self):
        from polyzymd.analyses.binding_free_energy import BFESettings

        s = BFESettings(units="kJ/mol")
        assert s.k_b() == pytest.approx(0.0083144626)

    def test_serialization_roundtrip(self):
        from polyzymd.analyses.binding_free_energy import BFESettings

        s = BFESettings(units="kcal/mol", fdr_alpha=0.01)
        d = s.model_dump()
        s2 = BFESettings(**d)
        assert s2.units == "kcal/mol"
        assert s2.fdr_alpha == pytest.approx(0.01)


# ===========================================================================
# compute_replicate / aggregate (no-op)
# ===========================================================================


class TestNoOp:
    def test_compute_replicate_returns_none(self):
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.binding_free_energy import (
            BFESettings,
            BindingFreeEnergyAnalysis,
        )

        analysis = BindingFreeEnergyAnalysis()
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
            settings=BFESettings(),
        )
        assert analysis.compute_replicate(ctx, 1) is None

    def test_aggregate_returns_none(self):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.binding_free_energy import (
            BFESettings,
            BindingFreeEnergyAnalysis,
        )

        analysis = BindingFreeEnergyAnalysis()
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
            settings=BFESettings(),
        )
        assert analysis.aggregate(ctx, []) is None


# ===========================================================================
# filter_conditions
# ===========================================================================


class TestFilterConditions:
    def test_keeps_all_conditions(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        mock_sim = MagicMock()

        conditions = [
            Condition(
                label="A",
                config_path=Path("/tmp/a.yaml"),
                replicates=(1,),
                sim_config=mock_sim,
            ),
            Condition(
                label="B",
                config_path=Path("/tmp/b.yaml"),
                replicates=(1,),
                sim_config=mock_sim,
            ),
        ]
        result = analysis.filter_conditions(conditions)
        assert len(result) == 2

    def test_keeps_no_polymer_conditions(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        mock_sim = MagicMock()
        mock_sim.polymers = None  # no polymer

        cond = Condition(
            label="No Polymer",
            config_path=Path("/tmp/np.yaml"),
            replicates=(1,),
            sim_config=mock_sim,
        )
        result = analysis.filter_conditions([cond])
        assert len(result) == 1  # BFE keeps all conditions


# ===========================================================================
# compare
# ===========================================================================


def _make_mock_bp_entry(
    partition_element="aromatic",
    mean_contact_share=0.3,
    expected_share=0.15,
    sem_contact_share=0.02,
    per_replicate_enrichments=None,
):
    """Build a mock AggregatedPartitionBindingEntry."""
    if per_replicate_enrichments is None:
        per_replicate_enrichments = [1.0, 0.8, 1.2]  # enrichment values

    entry = MagicMock()
    entry.partition_element = partition_element
    entry.mean_contact_share = mean_contact_share
    entry.expected_share = expected_share
    entry.sem_contact_share = sem_contact_share
    entry.per_replicate_enrichments = per_replicate_enrichments
    entry.n_exposed_in_group = 10
    return entry


def _make_mock_bp_result(entries=None, user_partitions=None):
    """Build a mock AggregatedBindingPreferenceResult with minimal structure."""
    if entries is None:
        entries = [_make_mock_bp_entry()]

    bp = MagicMock()
    partition_result = MagicMock()
    partition_result.entries = entries
    bp.aa_class_binding = {"SBM": partition_result}
    bp.user_defined_partitions = user_partitions or {}

    result = MagicMock()
    result.binding_preference = bp

    # Make isinstance check work for AggregatedBindingPreferenceResult
    return result


def _make_context(n_conditions=2, n_reps=3, control="A", temperature=300.0):
    """Build a ComparisonContext with mock conditions."""
    from polyzymd.analyses.base import ComparisonContext, Condition
    from polyzymd.analyses.binding_free_energy import BFESettings

    conditions = []
    analysis_dirs = {}
    for i in range(n_conditions):
        label = chr(65 + i)  # A, B, C, ...
        mock_sim = MagicMock()
        mock_sim.thermodynamics.temperature = temperature
        mock_sim.polymers = MagicMock()
        mock_sim.polymers.enabled = True
        cond = Condition(
            label=label,
            config_path=Path(f"/tmp/{label}/config.yaml"),
            replicates=tuple(range(1, n_reps + 1)),
            sim_config=mock_sim,
        )
        conditions.append(cond)
        analysis_dirs[label] = Path(f"/tmp/{label}/binding_free_energy")

    return ComparisonContext(
        name="test_comparison",
        conditions=conditions,
        excluded_conditions=[],
        control_label=control,
        analysis_dirs=analysis_dirs,
        results_dir=Path("/tmp/results"),
        equilibration="10ns",
        settings=BFESettings(),
        recompute=False,
    )


class TestCompare:
    def test_compare_returns_result(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis
        from polyzymd.analyses.binding_free_energy._comparison_results import FreeEnergyEntry

        analysis = BindingFreeEnergyAnalysis()
        ctx = _make_context(n_conditions=2)
        bp = _make_mock_bp_result()

        real_entry = FreeEnergyEntry(
            polymer_type="SBM",
            protein_group="aromatic",
            partition_name="aa_class",
            contact_share=0.2,
            expected_share=0.1,
            enrichment_ratio=2.0,
            delta_G=-0.5,
            delta_G_uncertainty=0.05,
            delta_G_per_replicate=[-0.4, -0.5, -0.6],
            units="kT",
            temperature_K=300.0,
            n_replicates=3,
            n_exposed_in_group=10,
        )

        with patch.object(analysis, "_load_binding_preference", return_value=bp):
            with patch.object(analysis, "_compute_dg_entries", return_value=[real_entry]):
                result = analysis.compare(ctx)

        assert result is not None
        assert result.units == "kT"
        assert len(result.conditions) == 2

    def test_compare_returns_none_no_data(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        ctx = _make_context(n_conditions=2)

        with patch.object(
            analysis,
            "_build_condition_summary",
            side_effect=ValueError("no data"),
        ):
            result = analysis.compare(ctx)

        assert result is None

    def test_compare_formula_kT(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        ctx = _make_context()

        with patch.object(analysis, "_load_binding_preference", return_value=None):
            result = analysis.compare(ctx)

        # All conditions got empty entries, so result is still created
        assert result is not None
        assert "k_bT" in result.formula

    def test_compare_mixed_temperatures(self):
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.binding_free_energy import (
            BFESettings,
            BindingFreeEnergyAnalysis,
        )

        analysis = BindingFreeEnergyAnalysis()

        conditions = []
        analysis_dirs = {}
        for i, (label, temp) in enumerate([("A", 300.0), ("B", 350.0)]):
            mock_sim = MagicMock()
            mock_sim.thermodynamics.temperature = temp
            cond = Condition(
                label=label,
                config_path=Path(f"/tmp/{label}/config.yaml"),
                replicates=(1, 2, 3),
                sim_config=mock_sim,
            )
            conditions.append(cond)
            analysis_dirs[label] = Path(f"/tmp/{label}/bfe")

        ctx = ComparisonContext(
            name="test",
            conditions=conditions,
            excluded_conditions=[],
            control_label="A",
            analysis_dirs=analysis_dirs,
            results_dir=Path("/tmp/results"),
            equilibration="10ns",
            settings=BFESettings(),
            recompute=False,
        )

        with patch.object(analysis, "_load_binding_preference", return_value=None):
            result = analysis.compare(ctx)

        assert result is not None
        assert result.mixed_temperatures is True
        assert len(result.temperature_groups) == 2


# ===========================================================================
# ΔG_sel computation
# ===========================================================================


class TestDGComputation:
    def test_entry_from_agg_bp_entry_basic(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()

        entry = _make_mock_bp_entry(
            mean_contact_share=0.2,
            expected_share=0.1,
            sem_contact_share=0.01,
            per_replicate_enrichments=[1.0, 0.8, 1.2],
        )

        fe = analysis._entry_from_agg_bp_entry(entry, "SBM", "aa_class", 1.0, "kT", 300.0)

        assert fe is not None
        # enrichment_ratio = 0.2 / 0.1 = 2.0
        assert fe.enrichment_ratio == pytest.approx(2.0)
        # ΔG = -kT * ln(2.0) ≈ -0.693
        assert fe.delta_G == pytest.approx(-math.log(2.0))
        assert fe.n_replicates == 3
        assert len(fe.delta_G_per_replicate) == 3

    def test_entry_zero_contact_share_returns_none(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        entry = _make_mock_bp_entry(mean_contact_share=0.0, expected_share=0.1)
        fe = analysis._entry_from_agg_bp_entry(entry, "SBM", "aa_class", 1.0, "kT", 300.0)
        assert fe is None

    def test_entry_zero_expected_share_returns_none(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        entry = _make_mock_bp_entry(mean_contact_share=0.2, expected_share=0.0)
        fe = analysis._entry_from_agg_bp_entry(entry, "SBM", "aa_class", 1.0, "kT", 300.0)
        assert fe is None

    def test_entry_negative_enrichment_gives_nan(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        # enrichment = -2.0 → ratio = -1.0 → log undefined
        entry = _make_mock_bp_entry(
            mean_contact_share=0.2,
            expected_share=0.1,
            per_replicate_enrichments=[-2.0],
        )
        fe = analysis._entry_from_agg_bp_entry(entry, "SBM", "aa_class", 1.0, "kT", 300.0)
        assert fe is not None
        assert math.isnan(fe.delta_G_per_replicate[0])

    def test_entry_kcal_units(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        kT = 0.0019872041 * 300.0  # kcal/(mol·K) * 300K

        entry = _make_mock_bp_entry(
            mean_contact_share=0.2,
            expected_share=0.1,
            per_replicate_enrichments=[1.0],
        )
        fe = analysis._entry_from_agg_bp_entry(entry, "SBM", "aa_class", kT, "kcal/mol", 300.0)
        assert fe is not None
        assert fe.units == "kcal/mol"
        assert fe.delta_G == pytest.approx(-kT * math.log(2.0))

    def test_uncertainty_propagation(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        entry = _make_mock_bp_entry(
            mean_contact_share=0.2,
            expected_share=0.1,
            sem_contact_share=0.02,
            per_replicate_enrichments=[1.0],
        )
        fe = analysis._entry_from_agg_bp_entry(entry, "SBM", "aa_class", 1.0, "kT", 300.0)
        assert fe is not None
        assert fe.delta_G_uncertainty is not None
        # σ = kT * σ_cs / cs = 1.0 * 0.02 / 0.2 = 0.1
        assert fe.delta_G_uncertainty == pytest.approx(0.1)


# ===========================================================================
# Pairwise comparison
# ===========================================================================


class TestPairwise:
    def _make_summaries(self, labels, temperature=300.0, dg_values=None):
        """Build mock FreeEnergyConditionSummary objects."""
        from polyzymd.analyses.binding_free_energy._comparison_results import (
            FreeEnergyConditionSummary,
            FreeEnergyEntry,
        )

        if dg_values is None:
            dg_values = {label: -0.5 - 0.1 * i for i, label in enumerate(labels)}

        summaries = []
        for label in labels:
            dg = dg_values[label]
            entry = FreeEnergyEntry(
                polymer_type="SBM",
                protein_group="aromatic",
                partition_name="aa_class",
                contact_share=0.2,
                expected_share=0.1,
                enrichment_ratio=2.0,
                delta_G=dg,
                delta_G_uncertainty=0.05,
                delta_G_per_replicate=[dg - 0.1, dg, dg + 0.1],
                units="kT",
                temperature_K=temperature,
                n_replicates=3,
                n_exposed_in_group=10,
            )
            summaries.append(
                FreeEnergyConditionSummary(
                    label=label,
                    config_path=f"/tmp/{label}/config.yaml",
                    temperature_K=temperature,
                    n_replicates=3,
                    units="kT",
                    entries=[entry],
                    polymer_types=["SBM"],
                    protein_groups=["aromatic"],
                )
            )
        return summaries

    def test_pairwise_with_control(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        summaries = self._make_summaries(["A", "B", "C"])

        pairs = analysis._compute_pairwise(summaries, effective_control="A")

        # Control A vs B, A vs C = 2 pairs × 1 (polymer, group) = 2
        assert len(pairs) == 2
        for p in pairs:
            assert p.condition_a == "A"

    def test_pairwise_no_control(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        summaries = self._make_summaries(["A", "B", "C"])

        pairs = analysis._compute_pairwise(summaries, effective_control=None)

        # All pairs: AB, AC, BC = 3 × 1 = 3
        assert len(pairs) == 3

    def test_pairwise_control_no_data_falls_back(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis
        from polyzymd.analyses.binding_free_energy._comparison_results import (
            FreeEnergyConditionSummary,
        )

        analysis = BindingFreeEnergyAnalysis()

        # Control has no entries (no-polymer condition)
        control = FreeEnergyConditionSummary(
            label="Control",
            config_path="/tmp/ctrl.yaml",
            temperature_K=300.0,
            n_replicates=0,
            units="kT",
            entries=[],
            polymer_types=[],
            protein_groups=[],
        )
        summaries = self._make_summaries(["B", "C"])
        all_summaries = [control] + summaries

        pairs = analysis._compute_pairwise(all_summaries, effective_control="Control")

        # Falls back to all-pairs among B, C = 1
        assert len(pairs) == 1
        labels = {(p.condition_a, p.condition_b) for p in pairs}
        assert ("B", "C") in labels

    def test_cross_temperature_suppresses_stats(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()

        summaries_300 = self._make_summaries(["A"], temperature=300.0)
        summaries_350 = self._make_summaries(["B"], temperature=350.0)

        pairs = analysis._compare_condition_pair(summaries_300[0], summaries_350[0])

        assert len(pairs) == 1
        assert pairs[0].cross_temperature is True
        assert pairs[0].t_statistic is None
        assert pairs[0].p_value is None

    def test_same_temperature_has_stats(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        summaries = self._make_summaries(["A", "B"], temperature=300.0)

        pairs = analysis._compare_condition_pair(summaries[0], summaries[1])

        assert len(pairs) == 1
        assert pairs[0].cross_temperature is False
        assert pairs[0].t_statistic is not None
        assert pairs[0].p_value is not None

    def test_delta_delta_g_sign(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        # A has dg=-0.5, B has dg=-0.6
        summaries = self._make_summaries(["A", "B"], dg_values={"A": -0.5, "B": -0.6})

        pairs = analysis._compare_condition_pair(summaries[0], summaries[1])

        assert len(pairs) == 1
        # ΔΔG = dg_B - dg_A = -0.6 - (-0.5) = -0.1
        assert pairs[0].delta_delta_G == pytest.approx(-0.1)


# ===========================================================================
# Plot
# ===========================================================================


class TestPlot:
    def test_plot_creates_output_dir(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.binding_free_energy import (
            BFESettings,
            BindingFreeEnergyAnalysis,
        )
        from polyzymd.config.comparison import PlotSettings

        analysis = BindingFreeEnergyAnalysis()
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
            analysis_dirs={"A": tmp_path / "A" / "bfe"},
            results_dir=tmp_path / "results",
            output_dir=output_dir,
            settings=BFESettings(),
            plot_settings=PlotSettings(),
        )

        with (
            patch(
                "polyzymd.analyses.binding_free_energy._plot_bfe_heatmap",
                return_value=[],
            ) as mock_heatmap,
            patch(
                "polyzymd.analyses.binding_free_energy._plot_bfe_bars",
                return_value=[],
            ) as mock_bars,
        ):
            analysis.plot(ctx)

        assert output_dir.exists()
        assert mock_heatmap.called
        assert mock_bars.called

    def test_plot_catches_exceptions(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.binding_free_energy import (
            BFESettings,
            BindingFreeEnergyAnalysis,
        )
        from polyzymd.config.comparison import PlotSettings

        analysis = BindingFreeEnergyAnalysis()
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
            analysis_dirs={"A": tmp_path / "A" / "bfe"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=BFESettings(),
            plot_settings=PlotSettings(),
        )

        with (
            patch(
                "polyzymd.analyses.binding_free_energy._plot_bfe_heatmap",
                side_effect=RuntimeError("boom"),
            ),
            patch(
                "polyzymd.analyses.binding_free_energy._plot_bfe_bars",
                side_effect=RuntimeError("boom"),
            ),
        ):
            result = analysis.plot(ctx)

        assert result == []

    def test_plot_collects_paths(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.binding_free_energy import (
            BFESettings,
            BindingFreeEnergyAnalysis,
        )
        from polyzymd.config.comparison import PlotSettings

        analysis = BindingFreeEnergyAnalysis()
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
            analysis_dirs={"A": tmp_path / "A" / "bfe"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=BFESettings(),
            plot_settings=PlotSettings(),
        )

        heatmap_path = tmp_path / "plots" / "bfe_heatmap.png"
        bars_path = tmp_path / "plots" / "bfe_bars.png"

        with (
            patch(
                "polyzymd.analyses.binding_free_energy._plot_bfe_heatmap",
                return_value=[heatmap_path],
            ),
            patch(
                "polyzymd.analyses.binding_free_energy._plot_bfe_bars",
                return_value=[bars_path],
            ),
        ):
            result = analysis.plot(ctx)

        assert heatmap_path in result
        assert bars_path in result

    def test_plot_passes_plot_settings(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.binding_free_energy import (
            BFESettings,
            BindingFreeEnergyAnalysis,
        )
        from polyzymd.config.comparison import PlotSettings

        analysis = BindingFreeEnergyAnalysis()
        mock_sim = MagicMock()
        ps = PlotSettings()
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
            analysis_dirs={"A": tmp_path / "A" / "bfe"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=BFESettings(),
            plot_settings=ps,
        )

        with (
            patch(
                "polyzymd.analyses.binding_free_energy._plot_bfe_heatmap",
                return_value=[],
            ) as mock_heatmap,
            patch(
                "polyzymd.analyses.binding_free_energy._plot_bfe_bars",
                return_value=[],
            ) as mock_bars,
        ):
            analysis.plot(ctx)

        # 4th positional arg is plot_settings
        assert mock_heatmap.call_args[0][3] is ps
        assert mock_bars.call_args[0][3] is ps


# ===========================================================================
# extract_metrics
# ===========================================================================


class TestExtractMetrics:
    def test_returns_empty(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        assert analysis.extract_metrics(MagicMock()) == {}


# ===========================================================================
# Lifecycle
# ===========================================================================


class TestLifecycle:
    def test_instantiation(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        assert analysis.name == "binding_free_energy"

    def test_is_subclass_of_analysis(self):
        from polyzymd.analyses.base import Analysis
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        assert issubclass(BindingFreeEnergyAnalysis, Analysis)

    def test_has_all_required_methods(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        assert callable(analysis.compute_replicate)
        assert callable(analysis.aggregate)
        assert callable(analysis.compare)
        assert callable(analysis.filter_conditions)
        assert callable(analysis.plot)
        assert callable(analysis.extract_metrics)

    def test_repr(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        assert "binding_free_energy" in repr(analysis)


# ===========================================================================
# CondProxy (shared helper)
# ===========================================================================


class TestCondProxy:
    def test_cond_proxy_attributes(self):
        from polyzymd.analyses.shared.binding_preference_helpers import CondProxy

        proxy = CondProxy(label="A", config="/tmp/a.yaml")
        assert proxy.label == "A"
        assert proxy.config == "/tmp/a.yaml"


# ===========================================================================
# Temperature handling
# ===========================================================================


class TestTemperature:
    def test_get_temperature(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        mock_sim = MagicMock()
        mock_sim.thermodynamics.temperature = 350.0

        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim,
        )
        assert analysis._get_temperature(cond) == pytest.approx(350.0)
