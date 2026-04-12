"""Tests for PolymerAffinityAnalysis plugin — Phase C-8.

Coverage:
- Discovery (registered by name and alias)
- Class attributes (name, Settings, aliases, dependencies, min_replicates)
- Settings validation (threshold range, defaults)
- compute_replicate / aggregate (no-op)
- filter_conditions (excludes no-polymer)
- compare (full override with temperature-aware pairwise)
- affinity score computation helpers
- pairwise comparisons (cross-temperature suppression)
- per-replicate score computation
- aggregation helpers
- plot delegation
- extract_metrics (empty)
- lifecycle
"""

from __future__ import annotations

import math
import sys
import types
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

# ===========================================================================
# Discovery
# ===========================================================================


class TestDiscovery:
    """Polymer affinity plugin is discovered by pkgutil scanning."""

    def test_discovery_by_name(self):
        from polyzymd.analyses import get_analysis

        cls = get_analysis("polymer_affinity")
        assert cls.__name__ == "PolymerAffinityAnalysis"

    def test_discovery_by_alias(self):
        from polyzymd.analyses import get_analysis

        cls = get_analysis("pa")
        assert cls.__name__ == "PolymerAffinityAnalysis"

    def test_listed(self):
        from polyzymd.analyses import list_analyses

        analyses = list_analyses()
        assert "polymer_affinity" in analyses

    def test_all_names_includes_alias(self):
        from polyzymd.analyses import list_all_names

        names = list_all_names()
        assert "polymer_affinity" in names
        assert "pa" in names


# ===========================================================================
# Class attributes
# ===========================================================================


class TestClassAttributes:
    def test_name(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        assert PolymerAffinityAnalysis.name == "polymer_affinity"

    def test_settings_type(self):
        from polyzymd.analyses.polymer_affinity import (
            PolymerAffinityAnalysis,
            PolymerAffinitySettings,
        )

        assert PolymerAffinityAnalysis.Settings is PolymerAffinitySettings

    def test_aliases(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        assert PolymerAffinityAnalysis.aliases == ("pa",)

    def test_dependencies(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        assert PolymerAffinityAnalysis.dependencies == ("contacts",)

    def test_min_replicates(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        assert PolymerAffinityAnalysis.min_replicates == 1


# ===========================================================================
# Settings
# ===========================================================================


class TestSettings:
    def test_defaults(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinitySettings

        s = PolymerAffinitySettings()
        assert s.compute_binding_preference is True
        assert s.surface_exposure_threshold == pytest.approx(0.2)
        assert s.enzyme_pdb_for_sasa is None
        assert s.include_default_aa_groups is True
        assert s.protein_groups is None
        assert s.protein_partitions is None
        assert s.polymer_type_selections is None
        assert s.fdr_alpha == pytest.approx(0.05)

    def test_custom_threshold(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinitySettings

        s = PolymerAffinitySettings(surface_exposure_threshold=0.3)
        assert s.surface_exposure_threshold == pytest.approx(0.3)

    def test_threshold_boundary_zero(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinitySettings

        s = PolymerAffinitySettings(surface_exposure_threshold=0.0)
        assert s.surface_exposure_threshold == pytest.approx(0.0)

    def test_threshold_boundary_one(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinitySettings

        s = PolymerAffinitySettings(surface_exposure_threshold=1.0)
        assert s.surface_exposure_threshold == pytest.approx(1.0)

    def test_threshold_below_zero_raises(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinitySettings

        with pytest.raises(ValueError):
            PolymerAffinitySettings(surface_exposure_threshold=-0.1)

    def test_threshold_above_one_raises(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinitySettings

        with pytest.raises(ValueError):
            PolymerAffinitySettings(surface_exposure_threshold=1.1)

    def test_serialization_roundtrip(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinitySettings

        s = PolymerAffinitySettings(
            fdr_alpha=0.01,
            surface_exposure_threshold=0.3,
            compute_binding_preference=False,
        )
        d = s.model_dump()
        s2 = PolymerAffinitySettings(**d)
        assert s2.fdr_alpha == pytest.approx(0.01)
        assert s2.surface_exposure_threshold == pytest.approx(0.3)
        assert s2.compute_binding_preference is False

    def test_custom_protein_groups(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinitySettings

        groups = {"active_site": [1, 2, 3], "binding_loop": [50, 51, 52]}
        s = PolymerAffinitySettings(protein_groups=groups)
        assert s.protein_groups == groups


# ===========================================================================
# compute_replicate / aggregate (no-op)
# ===========================================================================


class TestNoOp:
    def test_compute_replicate_returns_none(self):
        """Compare-only plugins must no-op compute_replicate by contract."""
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.polymer_affinity import (
            PolymerAffinityAnalysis,
            PolymerAffinitySettings,
        )

        analysis = PolymerAffinityAnalysis()
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
            settings=PolymerAffinitySettings(),
        )
        assert analysis.compute_replicate(ctx, 1) is None

    def test_aggregate_returns_none(self):
        """Compare-only plugins must no-op aggregate by contract."""
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.polymer_affinity import (
            PolymerAffinityAnalysis,
            PolymerAffinitySettings,
        )

        analysis = PolymerAffinityAnalysis()
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
            settings=PolymerAffinitySettings(),
        )
        assert analysis.aggregate(ctx, []) is None

    def test_compare_only_stage_flags(self):
        """Polymer affinity disables compute and aggregate stages explicitly."""
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        assert PolymerAffinityAnalysis.has_compute_stage is False
        assert PolymerAffinityAnalysis.has_aggregate_stage is False

    def test_compare_is_overridden(self):
        """Plugin must override Analysis.compare for custom comparison logic."""
        from polyzymd.analyses.base import Analysis
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        assert PolymerAffinityAnalysis.compare is not Analysis.compare


# ===========================================================================
# filter_conditions
# ===========================================================================


class TestFilterConditions:
    def test_excludes_no_polymer_conditions(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()
        mock_sim_poly = MagicMock()
        mock_sim_poly.polymers = MagicMock(enabled=True)  # has polymer

        mock_sim_no_poly = MagicMock(spec=[])  # no polymers, no topology, no attrs

        conditions = [
            Condition(
                label="With Polymer",
                config_path=Path("/tmp/a.yaml"),
                replicates=(1,),
                sim_config=mock_sim_poly,
            ),
            Condition(
                label="No Polymer",
                config_path=Path("/tmp/b.yaml"),
                replicates=(1,),
                sim_config=mock_sim_no_poly,
            ),
        ]
        result = analysis.filter_conditions(conditions)
        assert len(result) == 1
        assert result[0].label == "With Polymer"

    def test_keeps_all_when_all_have_polymer(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()
        mock_sim = MagicMock()
        mock_sim.polymers = MagicMock(enabled=True)

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

    def test_respects_polymer_chain_from_settings(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.polymer_affinity import (
            PolymerAffinityAnalysis,
            PolymerAffinitySettings,
        )

        analysis = PolymerAffinityAnalysis()
        mock_sim = MagicMock()
        mock_sim.polymers = None

        chain_b = MagicMock()
        chain_b.chain_id = "B"
        mock_sim.topology.chains = [chain_b]

        condition = Condition(
            label="With Polymer B",
            config_path=Path("/tmp/b.yaml"),
            replicates=(1,),
            sim_config=mock_sim,
        )

        settings = PolymerAffinitySettings(polymer_chain="B")
        result = analysis.filter_conditions([condition], settings=settings)

        assert len(result) == 1
        assert result[0].label == "With Polymer B"


# ===========================================================================
# compare
# ===========================================================================


def _make_mock_agg_bp_entry(
    polymer_type="SBM",
    protein_group="aromatic",
    mean_contact_fraction=0.3,
    n_exposed_in_group=10,
    mean_contact_share=0.2,
    expected_share=0.1,
    sem_contact_fraction=0.02,
    sem_contact_share=0.01,
    per_replicate_enrichments=None,
    n_replicates=3,
):
    """Build a mock AggregatedBindingPreferenceEntry for affinity scoring."""
    if per_replicate_enrichments is None:
        per_replicate_enrichments = [1.0, 0.8, 1.2]

    return SimpleNamespace(
        polymer_type=polymer_type,
        protein_group=protein_group,
        mean_contact_fraction=mean_contact_fraction,
        n_exposed_in_group=n_exposed_in_group,
        mean_contact_share=mean_contact_share,
        expected_share=expected_share,
        sem_contact_fraction=sem_contact_fraction,
        sem_contact_share=sem_contact_share,
        per_replicate_enrichments=per_replicate_enrichments,
        n_replicates=n_replicates,
    )


def _make_mock_agg_bp_result(entries=None):
    """Build a mock AggregatedBindingPreferenceResult."""
    if entries is None:
        entries = [_make_mock_agg_bp_entry()]

    return SimpleNamespace(entries=entries)


def _make_mock_single_bp_entry(
    polymer_type="SBM",
    protein_group="aromatic",
    mean_contact_fraction=0.3,
    n_exposed_in_group=10,
    contact_share=0.2,
    expected_share=0.1,
):
    """Build a mock single-replicate BindingPreferenceEntry."""
    return SimpleNamespace(
        polymer_type=polymer_type,
        protein_group=protein_group,
        mean_contact_fraction=mean_contact_fraction,
        n_exposed_in_group=n_exposed_in_group,
        contact_share=contact_share,
        expected_share=expected_share,
    )


def _make_pa_context(n_conditions=2, n_reps=3, control="A", temperature=300.0):
    """Build a ComparisonContext with mock conditions for polymer affinity."""
    from polyzymd.analyses.base import ComparisonContext, Condition
    from polyzymd.analyses.polymer_affinity import PolymerAffinitySettings

    conditions = []
    analysis_dirs = {}
    for i in range(n_conditions):
        label = chr(65 + i)  # A, B, C, ...
        mock_sim = MagicMock()
        mock_sim.thermodynamics.temperature = temperature
        mock_sim.polymer = MagicMock()  # has polymer
        cond = Condition(
            label=label,
            config_path=Path(f"/tmp/{label}/config.yaml"),
            replicates=tuple(range(1, n_reps + 1)),
            sim_config=mock_sim,
        )
        conditions.append(cond)
        analysis_dirs[label] = Path(f"/tmp/{label}/polymer_affinity")

    return ComparisonContext(
        name="test_comparison",
        conditions=conditions,
        excluded_conditions=[],
        control_label=control,
        analysis_dirs=analysis_dirs,
        results_dir=Path("/tmp/results"),
        equilibration="10ns",
        settings=PolymerAffinitySettings(),
        recompute=False,
    )


def _make_condition_summary(label, temperature=300.0, entries=None):
    """Build a real AffinityScoreConditionSummary for compare-level tests."""
    from polyzymd.analyses.polymer_affinity._comparison_results import (
        AffinityScoreConditionSummary,
        AffinityScoreEntry,
    )

    if entries is None:
        entries = [
            AffinityScoreEntry(
                polymer_type="SBM",
                protein_group="aromatic",
                n_contacts=3.0,
                affinity_score=-2.0,
                affinity_score_per_replicate=[-1.8, -2.0, -2.2],
                mean_contact_fraction=0.3,
                n_exposed_in_group=10,
                contact_share=0.2,
                expected_share=0.1,
                temperature_K=temperature,
                n_replicates=3,
            ),
        ]

    total_score = sum(e.affinity_score for e in entries if e.affinity_score is not None)
    polymer_types = sorted({e.polymer_type for e in entries})
    protein_groups = sorted({e.protein_group for e in entries})

    return AffinityScoreConditionSummary(
        label=label,
        config_path=f"/tmp/{label}/config.yaml",
        temperature_K=temperature,
        n_replicates=3,
        total_score=total_score,
        total_score_per_replicate=[-1.8, -2.0, -2.2],
        entries=entries,
        polymer_types=polymer_types,
        protein_groups=protein_groups,
    )


class TestCompare:
    def test_compare_returns_result(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis
        from polyzymd.analyses.polymer_affinity._comparison_results import AffinityScoreEntry

        analysis = PolymerAffinityAnalysis()
        ctx = _make_pa_context(n_conditions=2)

        summary_a = _make_condition_summary("A")
        summary_b = _make_condition_summary(
            "B",
            entries=[
                AffinityScoreEntry(
                    polymer_type="SBM",
                    protein_group="aromatic",
                    n_contacts=2.0,
                    affinity_score=-1.0,
                    affinity_score_per_replicate=[-0.9, -1.0, -1.1],
                    mean_contact_fraction=0.2,
                    n_exposed_in_group=10,
                    contact_share=0.2,
                    expected_share=0.1,
                    temperature_K=300.0,
                    n_replicates=3,
                )
            ],
        )

        summaries = [summary_a, summary_b]
        side_effects = iter(summaries)

        with patch.object(
            analysis, "_build_condition_summary", side_effect=lambda *a, **kw: next(side_effects)
        ):
            result = analysis.compare(ctx)

        assert result is not None
        assert len(result.conditions) == 2
        assert result.name == "test_comparison"
        assert [condition.label for condition in result.conditions] == ["A", "B"]
        assert result.mixed_temperatures is False
        assert result.temperature_groups == {"300.0": ["A", "B"]}
        assert len(result.pairwise_comparisons) == 1
        pair = result.pairwise_comparisons[0]
        assert pair.condition_a == "A"
        assert pair.condition_b == "B"
        assert pair.delta_score == pytest.approx(1.0)

        # More negative total score ranks first
        assert [condition.label for condition in result.get_ranking()] == ["A", "B"]

    def test_compare_returns_none_no_data(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()
        ctx = _make_pa_context(n_conditions=2)

        with patch.object(
            analysis,
            "_build_condition_summary",
            side_effect=ValueError("no data"),
        ):
            result = analysis.compare(ctx)

        assert result is None

    def test_compare_returns_none_all_summaries_empty(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis
        from polyzymd.analyses.polymer_affinity._comparison_results import (
            AffinityScoreConditionSummary,
        )

        analysis = PolymerAffinityAnalysis()
        ctx = _make_pa_context(n_conditions=2)

        summary_a = AffinityScoreConditionSummary(
            label="A",
            config_path="/tmp/A/config.yaml",
            temperature_K=300.0,
            n_replicates=0,
            entries=[],
            polymer_types=[],
            protein_groups=[],
        )
        summary_b = summary_a.model_copy(update={"label": "B", "config_path": "/tmp/B/config.yaml"})

        with patch.object(
            analysis,
            "_build_condition_summary",
            side_effect=[summary_a, summary_b],
        ):
            result = analysis.compare(ctx)

        assert result is None

    def test_compare_propagates_unexpected_runtime_error(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()
        ctx = _make_pa_context(n_conditions=1)

        with patch.object(
            analysis,
            "_build_condition_summary",
            side_effect=RuntimeError("unexpected compare failure"),
        ):
            with pytest.raises(RuntimeError, match="unexpected compare failure"):
                analysis.compare(ctx)

    def test_compare_mixed_temperatures(self):
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.polymer_affinity import (
            PolymerAffinityAnalysis,
            PolymerAffinitySettings,
        )

        analysis = PolymerAffinityAnalysis()

        conditions = []
        analysis_dirs = {}
        for label, temp in [("A", 300.0), ("B", 350.0)]:
            mock_sim = MagicMock()
            mock_sim.thermodynamics.temperature = temp
            mock_sim.polymer = MagicMock()
            cond = Condition(
                label=label,
                config_path=Path(f"/tmp/{label}/config.yaml"),
                replicates=(1, 2, 3),
                sim_config=mock_sim,
            )
            conditions.append(cond)
            analysis_dirs[label] = Path(f"/tmp/{label}/polymer_affinity")

        ctx = ComparisonContext(
            name="test",
            conditions=conditions,
            excluded_conditions=[],
            control_label="A",
            analysis_dirs=analysis_dirs,
            results_dir=Path("/tmp/results"),
            equilibration="10ns",
            settings=PolymerAffinitySettings(),
            recompute=False,
        )

        summaries = [
            _make_condition_summary("A", temperature=300.0),
            _make_condition_summary("B", temperature=350.0),
        ]
        side_effects = iter(summaries)

        with patch.object(
            analysis, "_build_condition_summary", side_effect=lambda *a, **kw: next(side_effects)
        ):
            result = analysis.compare(ctx)

        assert result is not None
        assert result.mixed_temperatures is True
        assert len(result.temperature_groups) == 2
        assert result.temperature_groups == {"300.0": ["A"], "350.0": ["B"]}
        assert len(result.pairwise_comparisons) == 1
        pair = result.pairwise_comparisons[0]
        assert pair.condition_a == "A"
        assert pair.condition_b == "B"
        assert pair.cross_temperature is True
        assert pair.t_statistic is None
        assert pair.p_value is None

    def test_compare_collects_metadata(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis
        from polyzymd.analyses.polymer_affinity._comparison_results import AffinityScoreEntry

        analysis = PolymerAffinityAnalysis()
        ctx = _make_pa_context(n_conditions=2)

        entries_a = [
            AffinityScoreEntry(
                polymer_type="SBM",
                protein_group="aromatic",
                n_contacts=3.0,
                affinity_score=-2.0,
                affinity_score_per_replicate=[-1.8, -2.0, -2.2],
                mean_contact_fraction=0.3,
                n_exposed_in_group=10,
                contact_share=0.2,
                expected_share=0.1,
                temperature_K=300.0,
                n_replicates=3,
            ),
            AffinityScoreEntry(
                polymer_type="SBM",
                protein_group="charged",
                n_contacts=2.0,
                affinity_score=-1.0,
                affinity_score_per_replicate=[-0.8, -1.0, -1.2],
                mean_contact_fraction=0.2,
                n_exposed_in_group=10,
                contact_share=0.15,
                expected_share=0.1,
                temperature_K=300.0,
                n_replicates=3,
            ),
            AffinityScoreEntry(
                polymer_type="EGM",
                protein_group="aromatic",
                n_contacts=1.5,
                affinity_score=-0.5,
                affinity_score_per_replicate=[-0.4, -0.5, -0.6],
                mean_contact_fraction=0.15,
                n_exposed_in_group=10,
                contact_share=0.1,
                expected_share=0.08,
                temperature_K=300.0,
                n_replicates=3,
            ),
        ]

        summaries = [
            _make_condition_summary("A", entries=entries_a),
            _make_condition_summary("B", entries=entries_a),
        ]
        side_effects = iter(summaries)

        with patch.object(
            analysis, "_build_condition_summary", side_effect=lambda *a, **kw: next(side_effects)
        ):
            result = analysis.compare(ctx)

        assert result is not None
        assert result.polymer_types == ["EGM", "SBM"]
        assert result.protein_groups == ["aromatic", "charged"]
        assert [condition.label for condition in result.conditions] == ["A", "B"]
        assert result.conditions[0].total_score == pytest.approx(-3.5)
        assert result.conditions[1].total_score == pytest.approx(-3.5)


# ===========================================================================
# Affinity score computation
# ===========================================================================


class TestAffinityScoreComputation:
    def test_entry_from_agg_bp_entry_basic(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()

        agg_entry = _make_mock_agg_bp_entry(
            polymer_type="SBM",
            protein_group="aromatic",
            mean_contact_fraction=0.3,
            n_exposed_in_group=10,
            mean_contact_share=0.2,
            expected_share=0.1,
            per_replicate_enrichments=[1.0, 0.8, 1.2],
        )

        entry = analysis._entry_from_agg_bp_entry(agg_entry, 300.0, None)

        assert entry is not None
        # N = 0.3 * 10 = 3.0
        assert entry.n_contacts == pytest.approx(3.0)
        # ΔG = -ln(0.2/0.1) = -ln(2) ≈ -0.693
        assert entry.delta_G_per_contact == pytest.approx(-math.log(2.0))
        # S = N * ΔG = 3.0 * (-0.693) ≈ -2.079
        assert entry.affinity_score == pytest.approx(3.0 * -math.log(2.0))
        assert entry.polymer_type == "SBM"
        assert entry.protein_group == "aromatic"

    def test_entry_zero_contact_share(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()

        agg_entry = _make_mock_agg_bp_entry(mean_contact_share=0.0, expected_share=0.1)
        entry = analysis._entry_from_agg_bp_entry(agg_entry, 300.0, None)

        assert entry is not None
        assert entry.delta_G_per_contact is None
        assert entry.affinity_score is None

    def test_entry_zero_expected_share(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()

        agg_entry = _make_mock_agg_bp_entry(mean_contact_share=0.2, expected_share=0.0)
        entry = analysis._entry_from_agg_bp_entry(agg_entry, 300.0, None)

        assert entry is not None
        assert entry.delta_G_per_contact is None
        assert entry.affinity_score is None

    def test_per_replicate_from_enrichments(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()

        # enrichment values: ratio = enrichment + 1
        # ΔG_rep = -ln(enrichment + 1)
        # S_rep = N * ΔG_rep
        enrichments = [1.0, 0.5]  # ratios: 2.0, 1.5
        n_contacts = 3.0  # mcf=0.3, n_exposed=10

        agg_entry = _make_mock_agg_bp_entry(
            mean_contact_fraction=0.3,
            n_exposed_in_group=10,
            mean_contact_share=0.2,
            expected_share=0.1,
            per_replicate_enrichments=enrichments,
        )

        entry = analysis._entry_from_agg_bp_entry(agg_entry, 300.0, None)

        assert entry is not None
        assert len(entry.affinity_score_per_replicate) == 2
        expected_0 = n_contacts * (-math.log(2.0))  # enrichment=1.0 → ratio=2.0
        expected_1 = n_contacts * (-math.log(1.5))  # enrichment=0.5 → ratio=1.5
        assert entry.affinity_score_per_replicate[0] == pytest.approx(expected_0)
        assert entry.affinity_score_per_replicate[1] == pytest.approx(expected_1)

    def test_negative_enrichment_gives_nan_filtered(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()

        # enrichment = -2.0 → ratio = -1.0 → log undefined → NaN → filtered out
        agg_entry = _make_mock_agg_bp_entry(
            mean_contact_share=0.2,
            expected_share=0.1,
            per_replicate_enrichments=[-2.0],
        )

        entry = analysis._entry_from_agg_bp_entry(agg_entry, 300.0, None)

        assert entry is not None
        # NaN replicates are filtered out
        assert len(entry.affinity_score_per_replicate) == 0

    def test_uncertainty_from_replicates(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()

        enrichments = [1.0, 0.5, 1.5]  # 3 replicates → SEM
        agg_entry = _make_mock_agg_bp_entry(
            mean_contact_fraction=0.3,
            n_exposed_in_group=10,
            mean_contact_share=0.2,
            expected_share=0.1,
            per_replicate_enrichments=enrichments,
        )

        entry = analysis._entry_from_agg_bp_entry(agg_entry, 300.0, None)

        assert entry is not None
        assert entry.affinity_score_uncertainty is not None
        # SEM = std(reps, ddof=1) / sqrt(3)
        reps = entry.affinity_score_per_replicate
        expected_sem = float(np.std(reps, ddof=1) / np.sqrt(len(reps)))
        assert entry.affinity_score_uncertainty == pytest.approx(expected_sem)

    def test_uncertainty_analytical_fallback(self):
        """When only 1 valid replicate, use analytical error propagation."""
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()

        agg_entry = _make_mock_agg_bp_entry(
            mean_contact_fraction=0.3,
            n_exposed_in_group=10,
            mean_contact_share=0.2,
            expected_share=0.1,
            sem_contact_fraction=0.02,
            sem_contact_share=0.01,
            per_replicate_enrichments=[1.0],  # only 1 replicate
        )

        entry = analysis._entry_from_agg_bp_entry(agg_entry, 300.0, None)

        assert entry is not None
        # Only 1 replicate, so analytical fallback should be used
        assert entry.affinity_score_uncertainty is not None
        # σ(S) = √[(N·σ_ΔG)² + (ΔG·σ_N)²]
        n_contacts = 0.3 * 10  # = 3.0
        sigma_n = 0.02 * 10  # = 0.2
        sigma_dg = 0.01 / 0.2  # = 0.05
        delta_g = -math.log(0.2 / 0.1)  # = -ln(2)
        expected = math.sqrt((n_contacts * sigma_dg) ** 2 + (delta_g * sigma_n) ** 2)
        assert entry.affinity_score_uncertainty == pytest.approx(expected)

    def test_entries_from_single(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()

        bp_entry = _make_mock_single_bp_entry(
            polymer_type="SBM",
            protein_group="aromatic",
            mean_contact_fraction=0.3,
            n_exposed_in_group=10,
            contact_share=0.2,
            expected_share=0.1,
        )

        # Single-replicate result
        result = MagicMock()
        result.entries = [bp_entry]

        entries = analysis._entries_from_single(result, 300.0)

        assert len(entries) == 1
        assert entries[0].n_contacts == pytest.approx(3.0)
        assert entries[0].n_replicates == 1
        assert len(entries[0].affinity_score_per_replicate) == 1


# ===========================================================================
# Per-replicate from files
# ===========================================================================


class TestPerRepFromFiles:
    def test_per_rep_scores_matched(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()

        entry1 = MagicMock()
        entry1.polymer_type = "SBM"
        entry1.protein_group = "aromatic"
        entry1.mean_contact_fraction = 0.3
        entry1.n_exposed_in_group = 10
        entry1.contact_share = 0.2
        entry1.expected_share = 0.1

        entry2 = MagicMock()
        entry2.polymer_type = "SBM"
        entry2.protein_group = "aromatic"
        entry2.mean_contact_fraction = 0.25
        entry2.n_exposed_in_group = 10
        entry2.contact_share = 0.18
        entry2.expected_share = 0.1

        per_rep_data = {1: [entry1], 2: [entry2]}

        scores = analysis._per_rep_scores_from_files(per_rep_data, "SBM", "aromatic")

        assert len(scores) == 2
        # rep 1: N=3.0, ΔG=-ln(2.0), S=3.0*-ln(2.0)
        assert scores[0] == pytest.approx(3.0 * (-math.log(2.0)))
        # rep 2: N=2.5, ΔG=-ln(1.8), S=2.5*-ln(1.8)
        assert scores[1] == pytest.approx(2.5 * (-math.log(1.8)))

    def test_per_rep_scores_unmatched_gives_zero(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()

        entry = MagicMock()
        entry.polymer_type = "EGM"  # different polymer type
        entry.protein_group = "aromatic"

        per_rep_data = {1: [entry]}

        scores = analysis._per_rep_scores_from_files(per_rep_data, "SBM", "aromatic")

        assert len(scores) == 1
        assert scores[0] == 0.0

    def test_per_rep_scores_zero_shares_gives_zero(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()

        entry = MagicMock()
        entry.polymer_type = "SBM"
        entry.protein_group = "aromatic"
        entry.mean_contact_fraction = 0.3
        entry.n_exposed_in_group = 10
        entry.contact_share = 0.0  # zero → can't compute ΔG
        entry.expected_share = 0.1

        per_rep_data = {1: [entry]}

        scores = analysis._per_rep_scores_from_files(per_rep_data, "SBM", "aromatic")

        assert len(scores) == 1
        assert scores[0] == 0.0


# ===========================================================================
# Aggregation helpers
# ===========================================================================


class TestAggregation:
    def test_aggregate_polymer_type_scores(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis
        from polyzymd.analyses.polymer_affinity._comparison_results import AffinityScoreEntry

        analysis = PolymerAffinityAnalysis()

        entries = [
            AffinityScoreEntry(
                polymer_type="SBM",
                protein_group="aromatic",
                n_contacts=3.0,
                affinity_score=-2.0,
                affinity_score_per_replicate=[-1.8, -2.0, -2.2],
                mean_contact_fraction=0.3,
                n_exposed_in_group=10,
                contact_share=0.2,
                expected_share=0.1,
                temperature_K=300.0,
                n_replicates=3,
            ),
            AffinityScoreEntry(
                polymer_type="SBM",
                protein_group="charged",
                n_contacts=2.0,
                affinity_score=-1.0,
                affinity_score_per_replicate=[-0.8, -1.0, -1.2],
                mean_contact_fraction=0.2,
                n_exposed_in_group=10,
                contact_share=0.15,
                expected_share=0.1,
                temperature_K=300.0,
                n_replicates=3,
            ),
        ]

        scores = analysis._aggregate_polymer_type_scores(entries)

        assert len(scores) == 1  # one polymer type: SBM
        pts = scores[0]
        assert pts.polymer_type == "SBM"
        assert pts.total_score == pytest.approx(-3.0)
        assert pts.total_n_contacts == pytest.approx(5.0)
        # Per-replicate: [-1.8+-0.8, -2.0+-1.0, -2.2+-1.2] = [-2.6, -3.0, -3.4]
        assert len(pts.total_score_per_replicate) == 3
        assert pts.total_score_per_replicate[0] == pytest.approx(-2.6)
        assert pts.total_score_per_replicate[1] == pytest.approx(-3.0)
        assert pts.total_score_per_replicate[2] == pytest.approx(-3.4)

    def test_aggregate_multiple_polymer_types_sorted(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis
        from polyzymd.analyses.polymer_affinity._comparison_results import AffinityScoreEntry

        analysis = PolymerAffinityAnalysis()

        entries = [
            AffinityScoreEntry(
                polymer_type="SBM",
                protein_group="aromatic",
                n_contacts=3.0,
                affinity_score=-2.0,
                affinity_score_per_replicate=[-2.0],
                mean_contact_fraction=0.3,
                n_exposed_in_group=10,
                contact_share=0.2,
                expected_share=0.1,
                temperature_K=300.0,
                n_replicates=1,
            ),
            AffinityScoreEntry(
                polymer_type="EGM",
                protein_group="aromatic",
                n_contacts=2.0,
                affinity_score=-5.0,  # more negative → sorted first
                affinity_score_per_replicate=[-5.0],
                mean_contact_fraction=0.2,
                n_exposed_in_group=10,
                contact_share=0.2,
                expected_share=0.1,
                temperature_K=300.0,
                n_replicates=1,
            ),
        ]

        scores = analysis._aggregate_polymer_type_scores(entries)

        assert len(scores) == 2
        # Sorted by total_score ascending (most negative first)
        assert scores[0].polymer_type == "EGM"
        assert scores[1].polymer_type == "SBM"

    def test_sum_per_replicate_empty(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()
        assert analysis._sum_per_replicate_across_groups([]) == []

    def test_compute_total_per_replicate_scores_empty(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()
        assert analysis._compute_total_per_replicate_scores([]) == []


# ===========================================================================
# Pairwise comparisons
# ===========================================================================


class TestPairwise:
    def _make_summaries(self, labels, temperature=300.0, total_scores=None):
        """Build AffinityScoreConditionSummary objects."""
        from polyzymd.analyses.polymer_affinity._comparison_results import (
            AffinityScoreConditionSummary,
            AffinityScoreEntry,
        )

        if total_scores is None:
            total_scores = {label: -2.0 - 0.5 * i for i, label in enumerate(labels)}

        summaries = []
        for label in labels:
            score = total_scores[label]
            entry = AffinityScoreEntry(
                polymer_type="SBM",
                protein_group="aromatic",
                n_contacts=3.0,
                affinity_score=score,
                affinity_score_per_replicate=[score - 0.1, score, score + 0.1],
                mean_contact_fraction=0.3,
                n_exposed_in_group=10,
                contact_share=0.2,
                expected_share=0.1,
                temperature_K=temperature,
                n_replicates=3,
            )
            summaries.append(
                AffinityScoreConditionSummary(
                    label=label,
                    config_path=f"/tmp/{label}/config.yaml",
                    temperature_K=temperature,
                    n_replicates=3,
                    total_score=score,
                    total_score_per_replicate=[score - 0.1, score, score + 0.1],
                    entries=[entry],
                    polymer_types=["SBM"],
                    protein_groups=["aromatic"],
                )
            )
        return summaries

    def test_pairwise_with_control(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()
        summaries = self._make_summaries(["A", "B", "C"])

        pairs = analysis._compute_pairwise(summaries, effective_control="A")

        # Control A vs B, A vs C = 2 pairs
        assert len(pairs) == 2
        for p in pairs:
            assert p.condition_a == "A"

    def test_pairwise_no_control(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()
        summaries = self._make_summaries(["A", "B", "C"])

        pairs = analysis._compute_pairwise(summaries, effective_control=None)

        # All pairs: AB, AC, BC = 3
        assert len(pairs) == 3

    def test_pairwise_control_no_data_falls_back(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis
        from polyzymd.analyses.polymer_affinity._comparison_results import (
            AffinityScoreConditionSummary,
        )

        analysis = PolymerAffinityAnalysis()

        # Control has no entries
        control = AffinityScoreConditionSummary(
            label="Control",
            config_path="/tmp/ctrl.yaml",
            temperature_K=300.0,
            n_replicates=0,
            entries=[],
            polymer_types=[],
            protein_groups=[],
        )
        summaries = self._make_summaries(["B", "C"])
        all_summaries = [control] + summaries

        pairs = analysis._compute_pairwise(all_summaries, effective_control="Control")

        # Falls back to all-pairs among B, C = 1
        assert len(pairs) == 1
        pair_labels = {(p.condition_a, p.condition_b) for p in pairs}
        assert ("B", "C") in pair_labels

    def test_cross_temperature_suppresses_stats(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()

        summaries_300 = self._make_summaries(["A"], temperature=300.0)
        summaries_350 = self._make_summaries(["B"], temperature=350.0)

        pair = analysis._compare_total_scores(summaries_300[0], summaries_350[0])

        assert pair.cross_temperature is True
        assert pair.t_statistic is None
        assert pair.p_value is None

    def test_same_temperature_has_stats(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()
        summaries = self._make_summaries(["A", "B"], temperature=300.0)

        pair = analysis._compare_total_scores(summaries[0], summaries[1])

        assert pair.cross_temperature is False
        assert pair.t_statistic is not None
        assert pair.p_value is not None

    def test_delta_score_sign(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()
        # A has score=-2.0, B has score=-3.0
        summaries = self._make_summaries(["A", "B"], total_scores={"A": -2.0, "B": -3.0})

        pair = analysis._compare_total_scores(summaries[0], summaries[1])

        # delta = score_B - score_A = -3.0 - (-2.0) = -1.0
        assert pair.delta_score == pytest.approx(-1.0)


# ===========================================================================
# Plot
# ===========================================================================


class TestPlot:
    def test_plot_creates_output_dir(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.polymer_affinity import (
            PolymerAffinityAnalysis,
            PolymerAffinitySettings,
        )
        from polyzymd.config.comparison import PlotSettings

        analysis = PolymerAffinityAnalysis()
        mock_sim = MagicMock()
        mock_sim.polymer = MagicMock()
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
            analysis_dirs={"A": tmp_path / "A" / "polymer_affinity"},
            results_dir=tmp_path / "results",
            output_dir=output_dir,
            settings=PolymerAffinitySettings(),
            plot_settings=PlotSettings(),
        )

        with (
            patch(
                "polyzymd.analyses.polymer_affinity._plot_affinity_stacked_bars",
                return_value=[],
            ),
            patch(
                "polyzymd.analyses.polymer_affinity._plot_affinity_group_bars",
                return_value=[],
            ),
        ):
            analysis.plot(ctx)

        assert output_dir.exists()

    def test_plot_propagates_exceptions(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.polymer_affinity import (
            PolymerAffinityAnalysis,
            PolymerAffinitySettings,
        )
        from polyzymd.config.comparison import PlotSettings

        analysis = PolymerAffinityAnalysis()
        mock_sim = MagicMock()
        mock_sim.polymer = MagicMock()
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
            analysis_dirs={"A": tmp_path / "A" / "polymer_affinity"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=PolymerAffinitySettings(),
            plot_settings=PlotSettings(),
        )

        with (
            patch(
                "polyzymd.analyses.polymer_affinity._plot_affinity_stacked_bars",
                side_effect=RuntimeError("boom"),
            ),
            patch(
                "polyzymd.analyses.polymer_affinity._plot_affinity_group_bars",
                side_effect=RuntimeError("boom"),
            ),
        ):
            with pytest.raises(RuntimeError, match="boom"):
                analysis.plot(ctx)

    def test_plot_empty_labels(self, tmp_path):
        from polyzymd.analyses.base import PlotContext
        from polyzymd.analyses.polymer_affinity import (
            PolymerAffinityAnalysis,
            PolymerAffinitySettings,
        )
        from polyzymd.config.comparison import PlotSettings

        analysis = PolymerAffinityAnalysis()

        ctx = PlotContext(
            conditions=[],
            analysis_dirs={},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=PolymerAffinitySettings(),
            plot_settings=PlotSettings(),
        )

        result = analysis.plot(ctx)
        assert result == []

    def test_plot_collects_paths(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.polymer_affinity import (
            PolymerAffinityAnalysis,
            PolymerAffinitySettings,
        )
        from polyzymd.config.comparison import PlotSettings

        analysis = PolymerAffinityAnalysis()
        mock_sim = MagicMock()
        mock_sim.polymer = MagicMock()
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
            analysis_dirs={"A": tmp_path / "A" / "polymer_affinity"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=PolymerAffinitySettings(),
            plot_settings=PlotSettings(),
        )

        stacked_path = tmp_path / "stacked.png"
        group_path = tmp_path / "group.png"

        with (
            patch(
                "polyzymd.analyses.polymer_affinity._plot_affinity_stacked_bars",
                return_value=[stacked_path],
            ),
            patch(
                "polyzymd.analyses.polymer_affinity._plot_affinity_group_bars",
                return_value=[group_path],
            ),
        ):
            result = analysis.plot(ctx)

        assert stacked_path in result
        assert group_path in result

    def test_plot_passes_plot_settings(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.polymer_affinity import (
            PolymerAffinityAnalysis,
            PolymerAffinitySettings,
        )
        from polyzymd.config.comparison import PlotSettings

        analysis = PolymerAffinityAnalysis()
        mock_sim = MagicMock()
        mock_sim.polymer = MagicMock()
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
            analysis_dirs={"A": tmp_path / "A" / "polymer_affinity"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=PolymerAffinitySettings(),
            plot_settings=ps,
        )

        with (
            patch(
                "polyzymd.analyses.polymer_affinity._plot_affinity_stacked_bars",
                return_value=[],
            ) as mock_stacked,
            patch(
                "polyzymd.analyses.polymer_affinity._plot_affinity_group_bars",
                return_value=[],
            ) as mock_group,
        ):
            analysis.plot(ctx)

        # 4th positional arg (index 3) is plot_settings
        assert mock_stacked.call_args[0][3] is ps
        assert mock_group.call_args[0][3] is ps


# ===========================================================================
# extract_metrics
# ===========================================================================


class TestExtractMetrics:
    def test_returns_empty(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()
        assert analysis.extract_metrics(MagicMock()) == {}


# ===========================================================================
# Lifecycle
# ===========================================================================


class TestLifecycle:
    def test_instantiation(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()
        assert analysis.name == "polymer_affinity"

    def test_is_subclass_of_analysis(self):
        from polyzymd.analyses.base import Analysis
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        assert issubclass(PolymerAffinityAnalysis, Analysis)

    def test_plugin_contract_surface(self):
        """Lifecycle-facing contract uses explicit class variables and overrides."""
        from polyzymd.analyses.base import Analysis
        from polyzymd.analyses.polymer_affinity import (
            PolymerAffinityAnalysis,
            PolymerAffinitySettings,
        )

        analysis = PolymerAffinityAnalysis()
        assert PolymerAffinityAnalysis.name == "polymer_affinity"
        assert PolymerAffinityAnalysis.Settings is PolymerAffinitySettings
        assert PolymerAffinityAnalysis.has_compute_stage is False
        assert PolymerAffinityAnalysis.has_aggregate_stage is False
        assert PolymerAffinityAnalysis.compare is not Analysis.compare
        assert analysis.extract_metrics({"dummy": 1}) == {}
        assert "polymer_affinity" in str(analysis)

    def test_repr(self):
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()
        assert "polymer_affinity" in repr(analysis)


# ===========================================================================
# Binding preference loading
# ===========================================================================


class TestLoadBindingPreference:
    def test_non_cached_path_uses_real_condition(self, tmp_path):
        """Non-cached compute path should pass the full Condition object."""
        from polyzymd.analyses.polymer_affinity import (
            PolymerAffinityAnalysis,
            PolymerAffinitySettings,
        )

        analysis = PolymerAffinityAnalysis()
        ctx = _make_pa_context(n_conditions=1)
        cond = ctx.conditions[0]
        ctx = replace(ctx, recompute=True)

        settings = PolymerAffinitySettings(compute_binding_preference=True)
        enzyme_pdb = tmp_path / "enzyme.pdb"
        enzyme_pdb.write_text("ATOM\n")
        expected = object()

        with (
            patch(
                "polyzymd.analyses.shared.binding_preference_helpers.resolve_enzyme_pdb",
                return_value=enzyme_pdb,
            ),
            patch(
                "polyzymd.analyses.contacts._paths.find_contact_results_for_replicates",
                return_value={1: tmp_path / "contacts_rep1.json"},
            ),
            patch(
                "polyzymd.analyses.shared.binding_preference.compute_condition_binding_preference",
                return_value=expected,
            ) as mock_compute,
        ):
            loaded = analysis._load_binding_preference(cond, ctx, settings)

        assert loaded is expected
        assert mock_compute.call_args.kwargs["cond"] is cond


# ===========================================================================
# Temperature handling
# ===========================================================================


class TestTemperature:
    def test_get_temperature(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis

        analysis = PolymerAffinityAnalysis()
        mock_sim = MagicMock()
        mock_sim.thermodynamics.temperature = 350.0

        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim,
        )
        assert analysis._get_temperature(cond) == pytest.approx(350.0)


# ===========================================================================
# _condition_has_polymer
# ===========================================================================


class TestConditionHasPolymer:
    def test_has_polymer_attribute(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.polymer_affinity import _condition_has_polymer

        mock_sim = MagicMock()
        mock_sim.polymers = MagicMock(enabled=True)  # has polymer

        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim,
        )
        assert _condition_has_polymer(cond) is True

    def test_no_polymer_attribute(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.polymer_affinity import _condition_has_polymer

        mock_sim = MagicMock(spec=[])  # no attributes at all

        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim,
        )
        assert _condition_has_polymer(cond) is False

    def test_polymer_is_none(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.polymer_affinity import _condition_has_polymer

        mock_sim = MagicMock()
        mock_sim.polymers = None
        # Also check topology path
        mock_sim.topology = MagicMock()
        mock_sim.topology.chains = []

        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim,
        )
        assert _condition_has_polymer(cond) is False

    def test_polymer_via_chain_id(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.polymer_affinity import _condition_has_polymer

        mock_sim = MagicMock()
        mock_sim.polymers = None

        chain_c = MagicMock()
        chain_c.chain_id = "C"
        mock_sim.topology.chains = [chain_c]

        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim,
        )
        assert _condition_has_polymer(cond) is True

    def test_invalid_polymer_selection_raises(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.polymer_affinity import _condition_has_polymer

        mock_sim = MagicMock(spec=[])
        mock_sim.get_working_directory = MagicMock(return_value=Path("/tmp/run_1"))

        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim,
        )

        mock_universe = MagicMock()
        mock_universe.select_atoms.side_effect = ValueError("Invalid selection token")

        mda_module = types.ModuleType("MDAnalysis")
        mda_module.Universe = MagicMock(return_value=mock_universe)

        with (
            patch.dict(sys.modules, {"MDAnalysis": mda_module}),
            patch("polyzymd.analyses.shared.loader.TrajectoryLoader") as mock_loader_cls,
        ):
            mock_loader = MagicMock()
            mock_loader_cls.return_value = mock_loader
            mock_loader.find_topology.return_value = Path("/tmp/system.pdb")

            with pytest.raises(ValueError, match="Invalid polymer selection"):
                _condition_has_polymer(
                    cond,
                    polymer_type_selections={"bad": "not a valid selection"},
                )


# ===========================================================================
# Deserialization
# ===========================================================================


class TestDeserialization:
    def test_deserialize_result_fallback_to_json(self, tmp_path):
        """Without AggregatedResultClass, _deserialize_result returns a dict via json.loads."""
        from polyzymd.analyses.polymer_affinity import PolymerAffinityAnalysis
        from polyzymd.analyses.polymer_affinity._comparison_results import (
            PolymerAffinityScoreResult,
        )

        analysis = PolymerAffinityAnalysis()

        # AggregatedResultClass should not be set for compare-only plugins
        assert analysis.AggregatedResultClass is None

        result = PolymerAffinityScoreResult(
            name="test",
            conditions=[],
            pairwise_comparisons=[],
        )
        path = tmp_path / "result.json"
        result.save(path)

        loaded = analysis._deserialize_result(path)
        assert isinstance(loaded, dict)
        assert loaded["name"] == "test"
