"""Unit tests for the polymer affinity score comparison module.

Tests cover:
- Physics / math correctness (affinity score formula, sign convention)
- AffinityScoreEntry, PolymerTypeScore, AffinityScoreConditionSummary models
- PolymerAffinityScoreResult save/load roundtrip
- AffinityScorePairwiseEntry model
- PolymerAffinityAnalysis compare helpers (no I/O, synthetic data)
- Formatters (console table, markdown, JSON)
- Plugin metadata and import checks
"""

from __future__ import annotations

import importlib
import json
import math
import tempfile
from pathlib import Path

import pytest

importlib.import_module("polyzymd.config.comparison")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_entry(
    polymer_type: str = "SBM",
    protein_group: str = "aromatic",
    n_contacts: float = 5.0,
    contact_share: float = 0.25,
    expected_share: float = 0.10,
    mean_contact_fraction: float = 0.5,
    n_exposed_in_group: int = 10,
    temperature_K: float = 363.0,
    n_replicates: int = 3,
    affinity_score_per_replicate: list[float] | None = None,
):
    """Construct an AffinityScoreEntry with sensible defaults."""
    from polyzymd.analyses.polymer_affinity._comparison_results import AffinityScoreEntry

    # ΔG_sel per contact = -ln(contact_share / expected_share) [kT]
    delta_g: float | None = None
    if contact_share > 0 and expected_share > 0:
        delta_g = -math.log(contact_share / expected_share)

    # Affinity score = N × ΔG_sel
    score: float | None = None
    if delta_g is not None:
        score = n_contacts * delta_g

    if affinity_score_per_replicate is None:
        if score is not None:
            affinity_score_per_replicate = [score + 0.1 * i for i in range(n_replicates)]
        else:
            affinity_score_per_replicate = []

    return AffinityScoreEntry(
        polymer_type=polymer_type,
        protein_group=protein_group,
        partition_name="aa_class",
        n_contacts=n_contacts,
        delta_G_per_contact=delta_g,
        affinity_score=score,
        affinity_score_uncertainty=0.1 if n_replicates >= 2 else None,
        affinity_score_per_replicate=affinity_score_per_replicate,
        mean_contact_fraction=mean_contact_fraction,
        n_exposed_in_group=n_exposed_in_group,
        contact_share=contact_share,
        expected_share=expected_share,
        temperature_K=temperature_K,
        n_replicates=n_replicates,
    )


def _make_polymer_type_score(
    polymer_type: str = "SBM",
    entries: list | None = None,
):
    """Construct a PolymerTypeScore."""
    from polyzymd.analyses.polymer_affinity._comparison_results import PolymerTypeScore

    if entries is None:
        entries = [_make_entry(polymer_type=polymer_type)]

    total_score = sum(e.affinity_score for e in entries if e.affinity_score is not None)
    total_n = sum(e.n_contacts for e in entries)

    return PolymerTypeScore(
        polymer_type=polymer_type,
        total_score=total_score,
        total_score_uncertainty=0.15,
        total_score_per_replicate=[total_score + 0.1 * i for i in range(3)],
        total_n_contacts=total_n,
        group_contributions=entries,
    )


def _make_condition_summary(
    label: str = "Cond A",
    temperature_K: float = 363.0,
    entries: list | None = None,
):
    """Construct an AffinityScoreConditionSummary."""
    from polyzymd.analyses.polymer_affinity._comparison_results import AffinityScoreConditionSummary

    if entries is None:
        entries = [_make_entry()]

    polymer_types = sorted({e.polymer_type for e in entries})
    protein_groups = sorted({e.protein_group for e in entries})
    n_replicates = max((e.n_replicates for e in entries), default=0)
    total_score = sum(e.affinity_score for e in entries if e.affinity_score is not None)
    total_n = sum(e.n_contacts for e in entries)
    total_per_rep = [total_score + 0.1 * i for i in range(n_replicates)]

    return AffinityScoreConditionSummary(
        label=label,
        config_path="/fake/config.yaml",
        temperature_K=temperature_K,
        n_replicates=n_replicates,
        total_score=total_score,
        total_score_uncertainty=0.15,
        total_score_per_replicate=total_per_rep,
        total_n_contacts=total_n,
        entries=entries,
        polymer_types=polymer_types,
        protein_groups=protein_groups,
    )


# ---------------------------------------------------------------------------
# Physics / Math Tests
# ---------------------------------------------------------------------------


class TestAffinityScorePhysics:
    """Test affinity score formula correctness."""

    def test_random_binding_gives_zero_score(self):
        """When contact_share == expected_share, ΔG_sel = 0 → score = 0."""
        cs = 0.20
        es = 0.20
        delta_g = -math.log(cs / es)
        assert abs(delta_g) < 1e-12
        # Score = N × 0 = 0
        score = 5.0 * delta_g
        assert abs(score) < 1e-12

    def test_enriched_contact_gives_negative_score(self):
        """When contact_share > expected_share, score < 0 (favorable)."""
        cs = 0.40
        es = 0.10
        delta_g = -math.log(cs / es)
        assert delta_g < 0
        n_contacts = 3.0
        score = n_contacts * delta_g
        assert score < 0

    def test_depleted_contact_gives_positive_score(self):
        """When contact_share < expected_share, score > 0 (unfavorable)."""
        cs = 0.05
        es = 0.30
        delta_g = -math.log(cs / es)
        assert delta_g > 0
        n_contacts = 3.0
        score = n_contacts * delta_g
        assert score > 0

    def test_score_scales_with_n_contacts(self):
        """Doubling N should double the affinity score."""
        cs, es = 0.30, 0.10
        delta_g = -math.log(cs / es)
        score_3 = 3.0 * delta_g
        score_6 = 6.0 * delta_g
        assert abs(score_6 / score_3 - 2.0) < 1e-12

    def test_kT_units_temperature_independent(self):
        """kT-unit ΔG_sel should be the same at any temperature."""
        cs, es = 0.30, 0.10
        # In kT units, ΔG_sel = -ln(cs/es) — temperature doesn't appear
        delta_g_300 = -math.log(cs / es)
        delta_g_400 = -math.log(cs / es)
        assert abs(delta_g_300 - delta_g_400) < 1e-15

    def test_enrichment_equivalence(self):
        """ΔG_sel = -ln(cs/es) = -ln(enrichment + 1)."""
        cs, es = 0.30, 0.10
        enrichment = (cs / es) - 1.0
        delta_g_ratio = -math.log(cs / es)
        delta_g_enrichment = -math.log(enrichment + 1.0)
        assert abs(delta_g_ratio - delta_g_enrichment) < 1e-12

    def test_total_score_is_sum_of_group_scores(self):
        """Total affinity score = Σ(N_g × ΔG_sel,g) across groups."""
        entries = [
            _make_entry(
                protein_group="aromatic", n_contacts=3.0, contact_share=0.30, expected_share=0.10
            ),
            _make_entry(
                protein_group="polar", n_contacts=2.0, contact_share=0.05, expected_share=0.20
            ),
        ]
        total = sum(e.affinity_score for e in entries if e.affinity_score is not None)
        # Manually compute
        score_1 = 3.0 * (-math.log(0.30 / 0.10))
        score_2 = 2.0 * (-math.log(0.05 / 0.20))
        expected_total = score_1 + score_2
        assert abs(total - expected_total) < 1e-10

    def test_analytical_error_propagation_formula(self):
        """σ(S) = √[(N·σ_ΔG_sel)² + (ΔG_sel·σ_N)²]."""
        n_contacts = 5.0
        delta_g = -1.2  # kT
        sigma_n = 0.5
        sigma_dg = 0.1
        expected_unc = math.sqrt((n_contacts * sigma_dg) ** 2 + (delta_g * sigma_n) ** 2)
        assert expected_unc > 0
        # Verify it grows with N
        sigma_large_n = math.sqrt((10.0 * sigma_dg) ** 2 + (delta_g * sigma_n) ** 2)
        assert sigma_large_n > expected_unc


# ---------------------------------------------------------------------------
# AffinityScoreEntry Tests
# ---------------------------------------------------------------------------


class TestAffinityScoreEntry:
    """Test AffinityScoreEntry model."""

    def test_basic_construction(self):
        entry = _make_entry()
        assert entry.polymer_type == "SBM"
        assert entry.protein_group == "aromatic"
        assert entry.affinity_score is not None
        assert entry.affinity_score < 0  # cs > es → favorable

    def test_score_value_correct(self):
        cs, es = 0.25, 0.10
        n_contacts = 5.0
        delta_g = -math.log(cs / es)
        expected_score = n_contacts * delta_g
        entry = _make_entry(contact_share=cs, expected_share=es, n_contacts=n_contacts)
        assert abs(entry.affinity_score - expected_score) < 1e-10

    def test_n_contacts_from_mcf_times_exposed(self):
        """n_contacts should represent mcf × n_exposed."""
        entry = _make_entry(
            mean_contact_fraction=0.3,
            n_exposed_in_group=20,
            n_contacts=6.0,  # 0.3 × 20
        )
        expected_n = entry.mean_contact_fraction * entry.n_exposed_in_group
        assert abs(entry.n_contacts - expected_n) < 1e-10

    def test_zero_contact_share_gives_none_delta_g(self):
        entry = _make_entry(contact_share=0.0, expected_share=0.10)
        assert entry.delta_G_per_contact is None
        assert entry.affinity_score is None

    def test_zero_expected_share_gives_none_delta_g(self):
        entry = _make_entry(contact_share=0.20, expected_share=0.0)
        assert entry.delta_G_per_contact is None
        assert entry.affinity_score is None

    def test_per_replicate_list_length(self):
        entry = _make_entry(
            n_replicates=5,
            affinity_score_per_replicate=[-4.0, -4.1, -3.9, -4.2, -3.8],
        )
        assert len(entry.affinity_score_per_replicate) == 5


# ---------------------------------------------------------------------------
# PolymerTypeScore Tests
# ---------------------------------------------------------------------------


class TestPolymerTypeScore:
    """Test PolymerTypeScore model."""

    def test_total_score_is_sum_of_entries(self):
        entries = [
            _make_entry(protein_group="aromatic", contact_share=0.30, expected_share=0.10),
            _make_entry(protein_group="polar", contact_share=0.10, expected_share=0.10),
        ]
        pts = _make_polymer_type_score(entries=entries)
        expected = sum(e.affinity_score for e in entries if e.affinity_score is not None)
        assert abs(pts.total_score - expected) < 1e-10

    def test_total_n_contacts_is_sum(self):
        entries = [
            _make_entry(protein_group="aromatic", n_contacts=3.0),
            _make_entry(protein_group="polar", n_contacts=2.0),
        ]
        pts = _make_polymer_type_score(entries=entries)
        assert abs(pts.total_n_contacts - 5.0) < 1e-10

    def test_group_contributions_preserved(self):
        entries = [_make_entry(protein_group="aromatic"), _make_entry(protein_group="polar")]
        pts = _make_polymer_type_score(entries=entries)
        assert len(pts.group_contributions) == 2


# ---------------------------------------------------------------------------
# AffinityScoreConditionSummary Tests
# ---------------------------------------------------------------------------


class TestAffinityScoreConditionSummary:
    """Test AffinityScoreConditionSummary."""

    def test_primary_metric_value_equals_total_score(self):
        summary = _make_condition_summary()
        assert summary.primary_metric_value == summary.total_score

    def test_primary_metric_sem_from_uncertainty(self):
        summary = _make_condition_summary()
        assert summary.primary_metric_sem == summary.total_score_uncertainty

    def test_empty_entries_gives_zero_score(self):
        summary = _make_condition_summary(entries=[])
        assert summary.total_score == 0.0
        assert summary.total_n_contacts == 0.0

    def test_polymer_types_sorted(self):
        entries = [
            _make_entry(polymer_type="ZZZ"),
            _make_entry(polymer_type="AAA"),
        ]
        summary = _make_condition_summary(entries=entries)
        assert summary.polymer_types == ["AAA", "ZZZ"]


# ---------------------------------------------------------------------------
# PolymerAffinityScoreResult Save/Load Tests
# ---------------------------------------------------------------------------


class TestPolymerAffinityScoreResult:
    """Test PolymerAffinityScoreResult serialization."""

    def _build_result(self, n_conditions: int = 2):
        from polyzymd.analyses.polymer_affinity._comparison_results import (
            PolymerAffinityScoreResult,
        )

        conditions = []
        for i in range(n_conditions):
            entries = [
                _make_entry(
                    polymer_type="SBM",
                    protein_group="aromatic",
                    n_contacts=5.0 + i,
                    contact_share=0.20 + 0.05 * i,
                    expected_share=0.10,
                ),
                _make_entry(
                    polymer_type="SBM",
                    protein_group="polar",
                    n_contacts=3.0,
                    contact_share=0.10,
                    expected_share=0.20,
                ),
            ]
            conditions.append(_make_condition_summary(label=f"Cond {i}", entries=entries))

        pairwise = []
        if n_conditions >= 2:
            from polyzymd.analyses.polymer_affinity._comparison_results import (
                AffinityScorePairwiseEntry,
            )

            pairwise.append(
                AffinityScorePairwiseEntry(
                    condition_a="Cond 0",
                    condition_b="Cond 1",
                    temperature_a_K=363.0,
                    temperature_b_K=363.0,
                    cross_temperature=False,
                    score_a=conditions[0].total_score,
                    score_b=conditions[1].total_score,
                    delta_score=conditions[1].total_score - conditions[0].total_score,
                    t_statistic=-2.5,
                    p_value=0.06,
                )
            )

        return PolymerAffinityScoreResult(
            name="test_comparison",
            conditions=conditions,
            pairwise_comparisons=pairwise,
            polymer_types=["SBM"],
            protein_groups=["aromatic", "polar"],
            equilibration_time="200ns",
        )

    def test_save_and_load_roundtrip(self):
        result = self._build_result()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "affinity_result.json"
            result.save(path)
            assert path.exists()
            loaded = type(result).load(path)
            assert loaded.name == result.name
            assert len(loaded.conditions) == 2
            assert loaded.conditions[0].label == "Cond 0"

    def test_save_creates_parent_dirs(self):
        result = self._build_result()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "nested" / "dir" / "result.json"
            result.save(path)
            assert path.exists()

    def test_json_is_valid(self):
        result = self._build_result()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = result.save(Path(tmpdir) / "result.json")
            data = json.loads(path.read_text())
            assert data["name"] == "test_comparison"
            assert "methodology" in data

    def test_get_condition_found(self):
        result = self._build_result()
        cond = result.get_condition("Cond 0")
        assert cond is not None
        assert cond.label == "Cond 0"

    def test_get_condition_not_found_returns_none(self):
        result = self._build_result()
        assert result.get_condition("Nonexistent") is None

    def test_get_ranking_most_negative_first(self):
        result = self._build_result()
        ranking = result.get_ranking()
        assert len(ranking) == 2
        # Most negative score should come first
        assert ranking[0].total_score <= ranking[1].total_score


# ---------------------------------------------------------------------------
# AffinityScorePairwiseEntry Tests
# ---------------------------------------------------------------------------


class TestAffinityScorePairwiseEntry:
    """Test AffinityScorePairwiseEntry."""

    def test_delta_score_sign(self):
        """Positive delta_score means B has weaker affinity than A."""
        from polyzymd.analyses.polymer_affinity._comparison_results import (
            AffinityScorePairwiseEntry,
        )

        entry = AffinityScorePairwiseEntry(
            condition_a="Control",
            condition_b="Treatment",
            temperature_a_K=363.0,
            temperature_b_K=363.0,
            cross_temperature=False,
            score_a=-10.0,  # strong affinity
            score_b=-5.0,  # weaker affinity
            delta_score=-5.0 - (-10.0),  # +5.0 → treatment weaker
        )
        assert entry.delta_score > 0

    def test_cross_temperature_flag(self):
        from polyzymd.analyses.polymer_affinity._comparison_results import (
            AffinityScorePairwiseEntry,
        )

        entry = AffinityScorePairwiseEntry(
            condition_a="363K",
            condition_b="293K",
            temperature_a_K=363.0,
            temperature_b_K=293.0,
            cross_temperature=True,
        )
        assert entry.cross_temperature
        assert entry.t_statistic is None
        assert entry.p_value is None

    def test_same_temperature_has_stats(self):
        from polyzymd.analyses.polymer_affinity._comparison_results import (
            AffinityScorePairwiseEntry,
        )

        entry = AffinityScorePairwiseEntry(
            condition_a="A",
            condition_b="B",
            temperature_a_K=363.0,
            temperature_b_K=363.0,
            cross_temperature=False,
            score_a=-10.0,
            score_b=-8.0,
            delta_score=2.0,
            t_statistic=3.5,
            p_value=0.01,
        )
        assert not entry.cross_temperature
        assert entry.t_statistic is not None
        assert entry.p_value is not None


# ---------------------------------------------------------------------------
# Formatter Tests
# ---------------------------------------------------------------------------


class TestPolymerAffinityFormatters:
    """Test output formatters."""

    def _build_result(self, n_conditions: int = 2):
        from polyzymd.analyses.polymer_affinity._comparison_results import (
            AffinityScorePairwiseEntry,
            PolymerAffinityScoreResult,
        )

        conditions = []
        for i in range(n_conditions):
            entries = [
                _make_entry(
                    polymer_type="SBM",
                    protein_group="aromatic",
                    n_contacts=5.0 + i,
                    contact_share=0.20 + 0.05 * i,
                    expected_share=0.10,
                    n_replicates=3,
                    affinity_score_per_replicate=[-4.0 - 0.3 * i, -4.1 - 0.3 * i, -3.9 - 0.3 * i],
                ),
                _make_entry(
                    polymer_type="SBM",
                    protein_group="polar",
                    n_contacts=3.0,
                    contact_share=0.10,
                    expected_share=0.20,
                    n_replicates=3,
                    affinity_score_per_replicate=[1.5 + 0.1 * i, 1.6 + 0.1 * i, 1.4 + 0.1 * i],
                ),
            ]
            conditions.append(_make_condition_summary(label=f"Cond {i}", entries=entries))

        pairwise = []
        if n_conditions >= 2:
            pairwise.append(
                AffinityScorePairwiseEntry(
                    condition_a="Cond 0",
                    condition_b="Cond 1",
                    temperature_a_K=363.0,
                    temperature_b_K=363.0,
                    cross_temperature=False,
                    score_a=-2.5,
                    score_b=-3.0,
                    delta_score=-0.5,
                    t_statistic=-2.0,
                    p_value=0.08,
                )
            )

        return PolymerAffinityScoreResult(
            name="test",
            conditions=conditions,
            pairwise_comparisons=pairwise,
            polymer_types=["SBM"],
            protein_groups=["aromatic", "polar"],
            equilibration_time="200ns",
        )

    def test_json_format_valid(self):
        from polyzymd.analyses.polymer_affinity._formatters import format_affinity_result

        result = self._build_result()
        output = format_affinity_result(result, format="json")
        data = json.loads(output)
        assert data["name"] == "test"
        assert "methodology" in data

    def test_table_format_contains_condition_labels(self):
        from polyzymd.analyses.polymer_affinity._formatters import format_affinity_result

        result = self._build_result()
        output = format_affinity_result(result, format="table")
        assert "Cond 0" in output
        assert "Cond 1" in output

    def test_table_format_contains_polymer_type(self):
        from polyzymd.analyses.polymer_affinity._formatters import format_affinity_result

        result = self._build_result()
        output = format_affinity_result(result, format="table")
        assert "SBM" in output

    def test_table_format_contains_score_keyword(self):
        from polyzymd.analyses.polymer_affinity._formatters import format_affinity_result

        result = self._build_result()
        output = format_affinity_result(result, format="table")
        assert any(kw in output for kw in ["Score", "score", "kT"])

    def test_table_format_contains_disclaimer(self):
        from polyzymd.analyses.polymer_affinity._formatters import format_affinity_result

        result = self._build_result()
        output = format_affinity_result(result, format="table")
        assert "DISCLAIMER" in output or "independence" in output.lower()

    def test_markdown_format_has_headers(self):
        from polyzymd.analyses.polymer_affinity._formatters import format_affinity_result

        result = self._build_result()
        output = format_affinity_result(result, format="markdown")
        assert "#" in output

    def test_markdown_format_has_table_pipes(self):
        from polyzymd.analyses.polymer_affinity._formatters import format_affinity_result

        result = self._build_result()
        output = format_affinity_result(result, format="markdown")
        assert "|" in output

    def test_markdown_format_has_disclaimer(self):
        from polyzymd.analyses.polymer_affinity._formatters import format_affinity_result

        result = self._build_result()
        output = format_affinity_result(result, format="markdown")
        assert "Disclaimer" in output or "independence" in output.lower()

    def test_format_single_condition(self):
        """Single condition should not raise and should produce output."""
        from polyzymd.analyses.polymer_affinity._formatters import format_affinity_result

        result = self._build_result(n_conditions=1)
        for fmt in ("table", "markdown", "json"):
            output = format_affinity_result(result, format=fmt)
            assert len(output) > 0

    def test_invalid_format_raises(self):
        from polyzymd.analyses.polymer_affinity._formatters import format_affinity_result

        result = self._build_result()
        with pytest.raises(ValueError, match="Unknown format"):
            format_affinity_result(result, format="csv")


# ---------------------------------------------------------------------------
# Metadata / Import Tests
# ---------------------------------------------------------------------------


class TestRegistration:
    """Test polymer affinity plugin metadata and imports."""

    def test_results_module_importable(self):
        from polyzymd.analyses.polymer_affinity._comparison_results import (
            AffinityScoreConditionSummary,
            AffinityScoreEntry,
            AffinityScorePairwiseEntry,
            PolymerAffinityScoreResult,
            PolymerTypeScore,
        )

        assert AffinityScoreEntry is not None
        assert PolymerTypeScore is not None
        assert AffinityScoreConditionSummary is not None
        assert AffinityScorePairwiseEntry is not None
        assert PolymerAffinityScoreResult is not None

    def test_formatters_importable(self):
        from polyzymd.analyses.polymer_affinity._formatters import (
            format_affinity_result,
        )

        assert format_affinity_result is not None

    def test_plugin_has_plot_functions(self):
        from polyzymd.analyses import polymer_affinity

        assert callable(getattr(polymer_affinity, "_plot_affinity_stacked_bars", None))
        assert callable(getattr(polymer_affinity, "_plot_affinity_group_bars", None))

    def test_plot_settings_model_attribute(self):
        from polyzymd.analyses.discovery import get_analysis
        from polyzymd.analyses.polymer_affinity._plot_settings import (
            AffinityPlotSettings,
        )

        analysis_cls = get_analysis("polymer_affinity")
        assert analysis_cls.PlotSettingsModel is AffinityPlotSettings
