"""Unit tests for the polymer affinity score comparison module.

Tests cover:
- Physics / math correctness (affinity score formula, sign convention)
- AffinityScoreEntry, PolymerTypeScore, AffinityScoreConditionSummary models
- PolymerAffinityScoreResult save/load roundtrip
- AffinityScorePairwiseEntry model
- PolymerAffinityScoreComparator helpers (no I/O, synthetic data)
- Formatters (console table, markdown, JSON)
- Registry and import checks
"""

from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path

import pytest

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
    from polyzymd.compare.results.polymer_affinity import AffinityScoreEntry

    # ΔΔG per contact = -ln(contact_share / expected_share) [kT]
    delta_g: float | None = None
    if contact_share > 0 and expected_share > 0:
        delta_g = -math.log(contact_share / expected_share)

    # Affinity score = N × ΔΔG
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
    from polyzymd.compare.results.polymer_affinity import PolymerTypeScore

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
    from polyzymd.compare.results.polymer_affinity import AffinityScoreConditionSummary

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
        """When contact_share == expected_share, ΔΔG = 0 → score = 0."""
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
        """kT-unit ΔΔG should be the same at any temperature."""
        cs, es = 0.30, 0.10
        # In kT units, ΔΔG = -ln(cs/es) — temperature doesn't appear
        delta_g_300 = -math.log(cs / es)
        delta_g_400 = -math.log(cs / es)
        assert abs(delta_g_300 - delta_g_400) < 1e-15

    def test_enrichment_equivalence(self):
        """ΔΔG = -ln(cs/es) = -ln(enrichment + 1)."""
        cs, es = 0.30, 0.10
        enrichment = (cs / es) - 1.0
        delta_g_ratio = -math.log(cs / es)
        delta_g_enrichment = -math.log(enrichment + 1.0)
        assert abs(delta_g_ratio - delta_g_enrichment) < 1e-12

    def test_total_score_is_sum_of_group_scores(self):
        """Total affinity score = Σ(N_g × ΔΔG_g) across groups."""
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
        """σ(S) = √[(N·σ_ΔΔG)² + (ΔΔG·σ_N)²]."""
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
        from polyzymd.compare.results.polymer_affinity import PolymerAffinityScoreResult

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
            from polyzymd.compare.results.polymer_affinity import AffinityScorePairwiseEntry

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
        from polyzymd.compare.results.polymer_affinity import AffinityScorePairwiseEntry

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
        from polyzymd.compare.results.polymer_affinity import AffinityScorePairwiseEntry

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
        from polyzymd.compare.results.polymer_affinity import AffinityScorePairwiseEntry

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
# Comparator Helper Tests (no I/O)
# ---------------------------------------------------------------------------


class TestPolymerAffinityComparatorHelpers:
    """Test comparator helper methods using synthetic data."""

    def _make_comparator(self):
        """Build a comparator without requiring real config files."""
        from unittest.mock import MagicMock

        from polyzymd.compare.comparators.polymer_affinity import PolymerAffinityScoreComparator
        from polyzymd.compare.settings import PolymerAffinityScoreSettings

        config = MagicMock()
        config.name = "test"
        config.conditions = []
        config.control = None

        analysis_settings = PolymerAffinityScoreSettings()
        comparator = PolymerAffinityScoreComparator(config, analysis_settings)
        return comparator

    def test_metric_type_is_mean_based(self):
        from polyzymd.analysis.core.metric_type import MetricType

        comp = self._make_comparator()
        assert comp.metric_type == MetricType.MEAN_BASED

    def test_rank_summaries_most_negative_first(self):
        comp = self._make_comparator()
        entries_a = [_make_entry(contact_share=0.40, expected_share=0.10)]  # strongly negative
        entries_b = [_make_entry(contact_share=0.10, expected_share=0.10)]  # ~0
        s_a = _make_condition_summary(label="A (preferred)", entries=entries_a)
        s_b = _make_condition_summary(label="B (neutral)", entries=entries_b)
        ranked = comp._rank_summaries([s_b, s_a])
        # Most negative total score first
        assert ranked[0].label == "A (preferred)"

    def test_compare_total_scores_same_temperature(self):
        """Same-temperature pair should have t-stats when n_reps >= 2."""
        comp = self._make_comparator()
        entries_a = [
            _make_entry(
                affinity_score_per_replicate=[-8.0, -9.0, -7.5],
                n_replicates=3,
            )
        ]
        entries_b = [
            _make_entry(
                affinity_score_per_replicate=[-4.0, -5.0, -3.5],
                n_replicates=3,
            )
        ]
        sa = _make_condition_summary(label="A", temperature_K=363.0, entries=entries_a)
        sb = _make_condition_summary(label="B", temperature_K=363.0, entries=entries_b)

        pw = comp._compare_total_scores(sa, sb)
        assert not pw.cross_temperature
        assert pw.t_statistic is not None
        assert pw.p_value is not None

    def test_compare_total_scores_cross_temperature_suppresses_stats(self):
        """Cross-temperature pair should have no t-stats."""
        comp = self._make_comparator()
        entries_a = [
            _make_entry(
                affinity_score_per_replicate=[-8.0, -9.0, -7.5],
                n_replicates=3,
            )
        ]
        entries_b = [
            _make_entry(
                affinity_score_per_replicate=[-4.0, -5.0, -3.5],
                n_replicates=3,
            )
        ]
        sa = _make_condition_summary(label="A", temperature_K=363.0, entries=entries_a)
        sb = _make_condition_summary(label="B", temperature_K=293.0, entries=entries_b)

        pw = comp._compare_total_scores(sa, sb)
        assert pw.cross_temperature
        assert pw.t_statistic is None
        assert pw.p_value is None

    def test_aggregate_polymer_type_scores(self):
        """Should group entries by polymer type and sum scores."""
        comp = self._make_comparator()
        entries = [
            _make_entry(
                polymer_type="SBM",
                protein_group="aromatic",
                n_contacts=3.0,
                contact_share=0.30,
                expected_share=0.10,
            ),
            _make_entry(
                polymer_type="SBM",
                protein_group="polar",
                n_contacts=2.0,
                contact_share=0.15,
                expected_share=0.20,
            ),
            _make_entry(
                polymer_type="EGM",
                protein_group="aromatic",
                n_contacts=4.0,
                contact_share=0.25,
                expected_share=0.10,
            ),
        ]
        pts_list = comp._aggregate_polymer_type_scores(entries)
        assert len(pts_list) == 2

        # Find SBM and EGM
        sbm_pts = next(p for p in pts_list if p.polymer_type == "SBM")
        egm_pts = next(p for p in pts_list if p.polymer_type == "EGM")

        assert len(sbm_pts.group_contributions) == 2
        assert len(egm_pts.group_contributions) == 1

        # SBM total should be sum of its entries
        expected_sbm = sum(e.affinity_score for e in entries[:2] if e.affinity_score is not None)
        assert abs(sbm_pts.total_score - expected_sbm) < 1e-10

    def test_direction_labels(self):
        """Check direction label property exists."""
        comp = self._make_comparator()
        labels = comp._direction_labels
        assert len(labels) == 3
        assert "affinity" in labels[0].lower()

    def test_get_mean_value(self):
        """_get_mean_value should return total score."""
        comp = self._make_comparator()
        summary = _make_condition_summary()
        assert comp._get_mean_value(summary) == summary.total_score

    def test_get_replicate_values_with_data(self):
        """_get_replicate_values should return per-replicate total scores."""
        comp = self._make_comparator()
        summary = _make_condition_summary()
        reps = comp._get_replicate_values(summary)
        assert len(reps) == len(summary.total_score_per_replicate)


# ---------------------------------------------------------------------------
# Formatter Tests
# ---------------------------------------------------------------------------


class TestPolymerAffinityFormatters:
    """Test output formatters."""

    def _build_result(self, n_conditions: int = 2):
        from polyzymd.compare.results.polymer_affinity import (
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
        from polyzymd.compare.polymer_affinity_formatters import format_affinity_result

        result = self._build_result()
        output = format_affinity_result(result, format="json")
        data = json.loads(output)
        assert data["name"] == "test"
        assert "methodology" in data

    def test_table_format_contains_condition_labels(self):
        from polyzymd.compare.polymer_affinity_formatters import format_affinity_result

        result = self._build_result()
        output = format_affinity_result(result, format="table")
        assert "Cond 0" in output
        assert "Cond 1" in output

    def test_table_format_contains_polymer_type(self):
        from polyzymd.compare.polymer_affinity_formatters import format_affinity_result

        result = self._build_result()
        output = format_affinity_result(result, format="table")
        assert "SBM" in output

    def test_table_format_contains_score_keyword(self):
        from polyzymd.compare.polymer_affinity_formatters import format_affinity_result

        result = self._build_result()
        output = format_affinity_result(result, format="table")
        assert any(kw in output for kw in ["Score", "score", "kT"])

    def test_table_format_contains_disclaimer(self):
        from polyzymd.compare.polymer_affinity_formatters import format_affinity_result

        result = self._build_result()
        output = format_affinity_result(result, format="table")
        assert "DISCLAIMER" in output or "independence" in output.lower()

    def test_markdown_format_has_headers(self):
        from polyzymd.compare.polymer_affinity_formatters import format_affinity_result

        result = self._build_result()
        output = format_affinity_result(result, format="markdown")
        assert "#" in output

    def test_markdown_format_has_table_pipes(self):
        from polyzymd.compare.polymer_affinity_formatters import format_affinity_result

        result = self._build_result()
        output = format_affinity_result(result, format="markdown")
        assert "|" in output

    def test_markdown_format_has_disclaimer(self):
        from polyzymd.compare.polymer_affinity_formatters import format_affinity_result

        result = self._build_result()
        output = format_affinity_result(result, format="markdown")
        assert "Disclaimer" in output or "independence" in output.lower()

    def test_format_single_condition(self):
        """Single condition should not raise and should produce output."""
        from polyzymd.compare.polymer_affinity_formatters import format_affinity_result

        result = self._build_result(n_conditions=1)
        for fmt in ("table", "markdown", "json"):
            output = format_affinity_result(result, format=fmt)
            assert len(output) > 0

    def test_invalid_format_raises(self):
        from polyzymd.compare.polymer_affinity_formatters import format_affinity_result

        result = self._build_result()
        with pytest.raises(ValueError, match="Unknown format"):
            format_affinity_result(result, format="csv")


# ---------------------------------------------------------------------------
# Registry / Import Tests
# ---------------------------------------------------------------------------


class TestRegistration:
    """Test that polymer affinity types are registered properly."""

    def test_comparator_registered(self):
        from polyzymd.compare.core.registry import ComparatorRegistry

        assert ComparatorRegistry.is_registered("polymer_affinity")

    def test_analysis_settings_registered(self):
        from polyzymd.analysis.core.registry import AnalysisSettingsRegistry

        assert AnalysisSettingsRegistry.is_registered("polymer_affinity")

    def test_comparison_settings_registered(self):
        from polyzymd.analysis.core.registry import ComparisonSettingsRegistry

        assert ComparisonSettingsRegistry.is_registered("polymer_affinity")

    def test_plotter_stacked_bars_registered(self):
        from polyzymd.compare.plotter import PlotterRegistry

        assert PlotterRegistry.is_registered("affinity_stacked_bars")

    def test_plotter_group_bars_registered(self):
        from polyzymd.compare.plotter import PlotterRegistry

        assert PlotterRegistry.is_registered("affinity_group_bars")

    def test_plot_settings_registered(self):
        from polyzymd.analysis.core.registry import PlotSettingsRegistry

        assert PlotSettingsRegistry.is_registered("polymer_affinity")

    def test_comparators_init_exports(self):
        from polyzymd.compare.comparators import PolymerAffinityScoreComparator

        assert PolymerAffinityScoreComparator is not None

    def test_results_module_importable(self):
        from polyzymd.compare.results.polymer_affinity import (
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

    def test_results_init_exports(self):
        from polyzymd.compare.results import (
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
        from polyzymd.compare.polymer_affinity_formatters import (
            format_affinity_result,
        )

        assert format_affinity_result is not None

    def test_compare_init_exports(self):
        from polyzymd.compare import (
            PolymerAffinityScoreComparator,
            PolymerAffinityScoreComparisonSettings,
            PolymerAffinityScoreResult,
            PolymerAffinityScoreSettings,
            format_affinity_result,
        )

        assert PolymerAffinityScoreComparator is not None
        assert PolymerAffinityScoreResult is not None
        assert PolymerAffinityScoreSettings is not None
        assert PolymerAffinityScoreComparisonSettings is not None
        assert format_affinity_result is not None

    def test_plotters_init_imports_module(self):
        from polyzymd.compare.plotters import polymer_affinity

        assert polymer_affinity is not None


# ---------------------------------------------------------------------------
# Settings Tests
# ---------------------------------------------------------------------------


class TestPolymerAffinitySettings:
    """Test PolymerAffinityScoreSettings and ComparisonSettings."""

    def test_default_settings_construction(self):
        from polyzymd.compare.settings import PolymerAffinityScoreSettings

        s = PolymerAffinityScoreSettings()
        assert s.surface_exposure_threshold > 0

    def test_comparison_settings_defaults(self):
        from polyzymd.compare.settings import PolymerAffinityScoreComparisonSettings

        cs = PolymerAffinityScoreComparisonSettings()
        assert 0 < cs.fdr_alpha <= 1.0

    def test_analysis_type_string(self):
        from polyzymd.compare.settings import PolymerAffinityScoreSettings

        assert PolymerAffinityScoreSettings.analysis_type() == "polymer_affinity"

    def test_comparison_analysis_type_string(self):
        from polyzymd.compare.settings import PolymerAffinityScoreComparisonSettings

        assert PolymerAffinityScoreComparisonSettings.analysis_type() == "polymer_affinity"
