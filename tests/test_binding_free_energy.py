"""Unit tests for the binding free energy analysis module.

Tests cover:
- Physics / math correctness (ΔΔG formula, delta-method uncertainty)
- FreeEnergyEntry, FreeEnergyConditionSummary, BindingFreeEnergyResult models
- BindingFreeEnergyAnalysisSettings and BindingFreeEnergyComparisonSettings
- BindingFreeEnergyComparator helpers (no I/O, uses synthetic data)
- Formatters (console table, markdown, JSON)
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

KB_KCAL = 0.0019872041  # kcal mol⁻¹ K⁻¹
KB_KJ = 0.0083144626  # kJ mol⁻¹ K⁻¹
KT_KT = 1.0  # dimensionless (kT units)
T_REF = 300.0  # K
KT_KCAL = KB_KCAL * T_REF  # 0.5962 kcal/mol


def _make_entry(
    polymer_type: str = "SBM",
    protein_group: str = "aromatic",
    contact_share: float = 0.25,
    expected_share: float = 0.10,
    units: str = "kT",
    temperature_K: float = T_REF,
    n_replicates: int = 3,
    delta_G_per_replicate: list[float] | None = None,
    delta_G_uncertainty: float | None = None,
):
    """Construct a FreeEnergyEntry with sensible defaults."""
    from polyzymd.compare.results.binding_free_energy import FreeEnergyEntry

    enrichment_ratio = contact_share / expected_share
    if units == "kT":
        kT = KT_KT
    elif units == "kcal/mol":
        kT = KB_KCAL * temperature_K
    else:
        kT = KB_KJ * temperature_K
    delta_G = -kT * math.log(enrichment_ratio)

    if delta_G_per_replicate is None:
        delta_G_per_replicate = [delta_G + 0.01 * i for i in range(n_replicates)]

    return FreeEnergyEntry(
        polymer_type=polymer_type,
        protein_group=protein_group,
        partition_name="aa_class",
        contact_share=contact_share,
        expected_share=expected_share,
        enrichment_ratio=enrichment_ratio,
        delta_G=delta_G,
        delta_G_uncertainty=delta_G_uncertainty,
        delta_G_per_replicate=delta_G_per_replicate,
        units=units,
        temperature_K=temperature_K,
        n_replicates=n_replicates,
    )


def _make_condition_summary(
    label: str = "Cond A",
    temperature_K: float = T_REF,
    entries: list | None = None,
    units: str = "kT",
):
    """Construct a FreeEnergyConditionSummary."""
    from polyzymd.compare.results.binding_free_energy import FreeEnergyConditionSummary

    if entries is None:
        entries = [_make_entry()]

    polymer_types = sorted({e.polymer_type for e in entries})
    protein_groups = sorted({e.protein_group for e in entries})
    n_replicates = max((e.n_replicates for e in entries), default=0)

    return FreeEnergyConditionSummary(
        label=label,
        config_path="/fake/config.yaml",
        temperature_K=temperature_K,
        n_replicates=n_replicates,
        units=units,
        entries=entries,
        polymer_types=polymer_types,
        protein_groups=protein_groups,
    )


# ---------------------------------------------------------------------------
# Physics / Math Tests
# ---------------------------------------------------------------------------


class TestBoltzmannInversion:
    """Test ΔΔG formula correctness."""

    def test_random_binding_gives_zero_ddg(self):
        """When contact_share == expected_share, ΔΔG = 0."""
        cs = 0.20
        es = 0.20
        kT = KB_KCAL * T_REF
        ddg = -kT * math.log(cs / es)
        assert abs(ddg) < 1e-12

    def test_enriched_contact_gives_negative_ddg(self):
        """When contact_share > expected_share, ΔΔG < 0 (favorable)."""
        cs = 0.40  # polymer prefers this group
        es = 0.10  # only 10% surface area
        kT = KB_KCAL * T_REF
        ddg = -kT * math.log(cs / es)
        assert ddg < 0

    def test_depleted_contact_gives_positive_ddg(self):
        """When contact_share < expected_share, ΔΔG > 0 (unfavorable)."""
        cs = 0.05
        es = 0.30
        kT = KB_KCAL * T_REF
        ddg = -kT * math.log(cs / es)
        assert ddg > 0

    def test_equivalence_with_enrichment(self):
        """ΔΔG = -kT·ln(enrichment + 1) (exact, not approximation)."""
        cs = 0.30
        es = 0.10
        enrichment = (cs / es) - 1.0  # = 2.0
        kT = KB_KCAL * T_REF
        ddg_ratio = -kT * math.log(cs / es)
        ddg_enrichment = -kT * math.log(enrichment + 1.0)
        assert abs(ddg_ratio - ddg_enrichment) < 1e-12

    def test_delta_method_uncertainty(self):
        """σ(ΔΔG) = kT·√[(σ_cs/cs)² + (σ_es/es)²]."""
        cs, es = 0.30, 0.10
        sem_cs, sem_es = 0.02, 0.005
        kT = KB_KCAL * T_REF
        expected_unc = kT * math.sqrt((sem_cs / cs) ** 2 + (sem_es / es) ** 2)
        # Sanity check: uncertainty is positive
        assert expected_unc > 0
        # With zero σ_es (single SASA computation)
        unc_zero_es = kT * math.sqrt((sem_cs / cs) ** 2)
        assert abs(unc_zero_es - kT * sem_cs / cs) < 1e-12

    def test_kcal_vs_kj_scale(self):
        """kJ/mol values should be 4.184× kcal/mol values."""
        cs, es = 0.25, 0.10
        ddg_kcal = -(KB_KCAL * T_REF) * math.log(cs / es)
        ddg_kj = -(KB_KJ * T_REF) * math.log(cs / es)
        assert abs(ddg_kj / ddg_kcal - 4.184) < 0.01

    def test_kT_units_ddg_equals_negative_log(self):
        """In kT units, ΔΔG = -ln(contact_share / expected_share) exactly."""
        cs, es = 0.30, 0.10
        ddg_kT = -KT_KT * math.log(cs / es)
        expected = -math.log(cs / es)
        assert abs(ddg_kT - expected) < 1e-12

    def test_kT_units_temperature_independent(self):
        """kT-unit ΔΔG should be the same at any temperature."""
        cs, es = 0.30, 0.10
        ddg_300 = -KT_KT * math.log(cs / es)  # T=300K
        ddg_400 = -KT_KT * math.log(cs / es)  # T=400K — same formula
        assert abs(ddg_300 - ddg_400) < 1e-15


# ---------------------------------------------------------------------------
# Settings Tests
# ---------------------------------------------------------------------------


class TestBindingFreeEnergySettings:
    """Test BindingFreeEnergyAnalysisSettings and ComparisonSettings."""

    def test_default_units_kT(self):
        from polyzymd.compare.settings import BindingFreeEnergyAnalysisSettings

        s = BindingFreeEnergyAnalysisSettings()
        assert s.units == "kT"

    def test_k_b_kT_returns_zero(self):
        from polyzymd.compare.settings import BindingFreeEnergyAnalysisSettings

        s = BindingFreeEnergyAnalysisSettings(units="kT")
        assert s.k_b() == 0.0

    def test_k_b_kcal(self):
        from polyzymd.compare.settings import BindingFreeEnergyAnalysisSettings

        s = BindingFreeEnergyAnalysisSettings(units="kcal/mol")
        assert abs(s.k_b() - KB_KCAL) < 1e-12

    def test_k_b_kj(self):
        from polyzymd.compare.settings import BindingFreeEnergyAnalysisSettings

        s = BindingFreeEnergyAnalysisSettings(units="kJ/mol")
        assert abs(s.k_b() - KB_KJ) < 1e-12

    def test_invalid_units_raises(self):
        from pydantic import ValidationError

        from polyzymd.compare.settings import BindingFreeEnergyAnalysisSettings

        with pytest.raises(ValidationError):
            BindingFreeEnergyAnalysisSettings(units="eV")

    def test_comparison_settings_defaults(self):
        from polyzymd.compare.settings import BindingFreeEnergyComparisonSettings

        cs = BindingFreeEnergyComparisonSettings()
        assert 0 < cs.fdr_alpha <= 1.0

    def test_surface_exposure_threshold_positive(self):
        from pydantic import ValidationError

        from polyzymd.compare.settings import BindingFreeEnergyAnalysisSettings

        with pytest.raises(ValidationError):
            BindingFreeEnergyAnalysisSettings(surface_exposure_threshold=-0.1)


# ---------------------------------------------------------------------------
# FreeEnergyEntry Tests
# ---------------------------------------------------------------------------


class TestFreeEnergyEntry:
    """Test FreeEnergyEntry model."""

    def test_basic_construction(self):
        entry = _make_entry()
        assert entry.polymer_type == "SBM"
        assert entry.protein_group == "aromatic"
        assert entry.delta_G is not None
        assert entry.delta_G < 0  # cs > es → favorable

    def test_delta_G_value_correct(self):
        cs, es = 0.25, 0.10
        kT = KT_KT  # default units are now "kT" → kT = 1.0
        expected = -kT * math.log(cs / es)
        entry = _make_entry(contact_share=cs, expected_share=es)
        assert abs(entry.delta_G - expected) < 1e-10

    def test_enrichment_ratio_stored(self):
        cs, es = 0.30, 0.12
        entry = _make_entry(contact_share=cs, expected_share=es)
        assert abs(entry.enrichment_ratio - cs / es) < 1e-12

    def test_n_replicates_matches_per_rep_list(self):
        entry = _make_entry(n_replicates=4, delta_G_per_replicate=[-1.0, -1.1, -0.9, -1.05])
        assert entry.n_replicates == 4
        assert len(entry.delta_G_per_replicate) == 4

    def test_kj_units(self):
        cs, es = 0.25, 0.10
        kT = KB_KJ * T_REF
        expected = -kT * math.log(cs / es)
        entry = _make_entry(contact_share=cs, expected_share=es, units="kJ/mol")
        assert abs(entry.delta_G - expected) < 1e-10


# ---------------------------------------------------------------------------
# FreeEnergyConditionSummary Tests
# ---------------------------------------------------------------------------


class TestFreeEnergyConditionSummary:
    """Test FreeEnergyConditionSummary."""

    def test_primary_metric_value_mean_of_entries(self):
        entries = [
            _make_entry(contact_share=0.25, expected_share=0.10),  # negative ΔΔG
            _make_entry(
                protein_group="polar", contact_share=0.10, expected_share=0.20
            ),  # positive ΔΔG
        ]
        summary = _make_condition_summary(entries=entries)
        dg_vals = [e.delta_G for e in entries]
        expected_mean = sum(dg_vals) / len(dg_vals)
        assert abs(summary.primary_metric_value - expected_mean) < 1e-10

    def test_primary_metric_value_empty_entries(self):
        summary = _make_condition_summary(entries=[])
        assert summary.primary_metric_value == 0.0

    def test_get_entry_found(self):
        entry = _make_entry(polymer_type="SBM", protein_group="aromatic")
        summary = _make_condition_summary(entries=[entry])
        found = summary.get_entry("SBM", "aromatic")
        assert found is not None
        assert found.polymer_type == "SBM"
        assert found.protein_group == "aromatic"

    def test_get_entry_not_found(self):
        summary = _make_condition_summary(entries=[_make_entry()])
        assert summary.get_entry("SBM", "nonexistent") is None

    def test_polymer_types_sorted(self):
        entries = [
            _make_entry(polymer_type="ZZZ"),
            _make_entry(polymer_type="AAA"),
        ]
        summary = _make_condition_summary(entries=entries)
        assert summary.polymer_types == ["AAA", "ZZZ"]


# ---------------------------------------------------------------------------
# BindingFreeEnergyResult Save/Load Tests
# ---------------------------------------------------------------------------


class TestBindingFreeEnergyResult:
    """Test BindingFreeEnergyResult serialization."""

    def _build_result(self, units: str = "kT"):
        from polyzymd.compare.results.binding_free_energy import BindingFreeEnergyResult

        cond = _make_condition_summary(label="Cond A")
        return BindingFreeEnergyResult(
            name="test_comparison",
            units=units,
            conditions=[cond],
            polymer_types=["SBM"],
            protein_groups=["aromatic"],
            equilibration_time="10ns",
        )

    def test_save_and_load_roundtrip(self):
        result = self._build_result()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "bfe_result.json"
            result.save(path)
            assert path.exists()
            loaded = type(result).load(path)
            assert loaded.name == result.name
            assert loaded.units == result.units
            assert len(loaded.conditions) == 1
            assert loaded.conditions[0].label == "Cond A"

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
            assert data["units"] == "kT"

    def test_get_condition_found(self):
        result = self._build_result()
        cond = result.get_condition("Cond A")
        assert cond.label == "Cond A"

    def test_get_condition_not_found_raises(self):
        result = self._build_result()
        with pytest.raises(KeyError):
            result.get_condition("Nonexistent")

    def test_kj_units_preserved_in_roundtrip(self):
        result = self._build_result(units="kJ/mol")
        with tempfile.TemporaryDirectory() as tmpdir:
            path = result.save(Path(tmpdir) / "result_kj.json")
            loaded = type(result).load(path)
            assert loaded.units == "kJ/mol"


# ---------------------------------------------------------------------------
# FreeEnergyPairwiseEntry Tests
# ---------------------------------------------------------------------------


class TestFreeEnergyPairwiseEntry:
    """Test FreeEnergyPairwiseEntry."""

    def test_delta_delta_g_sign(self):
        """ΔΔG_B - ΔΔG_A: positive means B less favorable than A."""
        from polyzymd.compare.results.binding_free_energy import FreeEnergyPairwiseEntry

        entry = FreeEnergyPairwiseEntry(
            polymer_type="SBM",
            protein_group="aromatic",
            condition_a="Control",
            condition_b="Treatment",
            temperature_a_K=T_REF,
            temperature_b_K=T_REF,
            cross_temperature=False,
            delta_G_a=-1.5,
            delta_G_b=-0.5,
            delta_delta_G=-0.5 - (-1.5),  # +1.0 → treatment less favorable
        )
        assert entry.delta_delta_G > 0

    def test_cross_temperature_flag(self):
        from polyzymd.compare.results.binding_free_energy import FreeEnergyPairwiseEntry

        entry = FreeEnergyPairwiseEntry(
            polymer_type="SBM",
            protein_group="aromatic",
            condition_a="300K",
            condition_b="320K",
            temperature_a_K=300.0,
            temperature_b_K=320.0,
            cross_temperature=True,
        )
        assert entry.cross_temperature
        assert entry.t_statistic is None
        assert entry.p_value is None


# ---------------------------------------------------------------------------
# BindingFreeEnergyComparator Helper Tests (no I/O)
# ---------------------------------------------------------------------------


class TestBindingFreeEnergyComparatorHelpers:
    """Test comparator helper methods using synthetic data."""

    def _make_comparator(self, units: str = "kT"):
        """Build a comparator without requiring real config files."""
        from unittest.mock import MagicMock

        from polyzymd.compare.comparators.binding_free_energy import BindingFreeEnergyComparator
        from polyzymd.compare.settings import BindingFreeEnergyAnalysisSettings

        config = MagicMock()
        config.name = "test"
        config.conditions = []
        config.control = None

        analysis_settings = BindingFreeEnergyAnalysisSettings(units=units)
        comparator = BindingFreeEnergyComparator(config, analysis_settings)
        return comparator

    def test_metric_type_is_mean_based(self):
        from polyzymd.analysis.core.metric_type import MetricType

        comp = self._make_comparator()
        assert comp.metric_type == MetricType.MEAN_BASED

    def test_rank_summaries_most_negative_first(self):
        comp = self._make_comparator()
        entries_a = [_make_entry(contact_share=0.30, expected_share=0.10)]  # strongly negative
        entries_b = [_make_entry(contact_share=0.10, expected_share=0.10)]  # ~0
        s_a = _make_condition_summary(label="A (preferred)", entries=entries_a)
        s_b = _make_condition_summary(label="B (neutral)", entries=entries_b)
        ranked = comp._rank_summaries([s_b, s_a])
        assert ranked[0].label == "A (preferred)"

    def test_interpret_direction_negative(self):
        comp = self._make_comparator()
        assert "favorable" in comp._interpret_direction(-0.5)

    def test_interpret_direction_positive(self):
        comp = self._make_comparator()
        assert "favorable" in comp._interpret_direction(0.5)

    def test_interpret_direction_zero(self):
        comp = self._make_comparator()
        assert "unchanged" in comp._interpret_direction(0.0)

    def test_compare_condition_pair_same_temperature(self):
        """Same-temperature pair should have t-stats when n_reps >= 2."""
        comp = self._make_comparator()

        entries_a = [
            _make_entry(
                delta_G_per_replicate=[-1.2, -1.4, -1.1],
                n_replicates=3,
            )
        ]
        entries_b = [
            _make_entry(
                delta_G_per_replicate=[-0.5, -0.6, -0.7],
                n_replicates=3,
            )
        ]
        sa = _make_condition_summary(label="A", temperature_K=300.0, entries=entries_a)
        sb = _make_condition_summary(label="B", temperature_K=300.0, entries=entries_b)

        pairwise = comp._compare_condition_pair(sa, sb)
        assert len(pairwise) == 1
        pw = pairwise[0]
        assert not pw.cross_temperature
        assert pw.t_statistic is not None
        assert pw.p_value is not None

    def test_compare_condition_pair_cross_temperature_suppresses_stats(self):
        """Cross-temperature pair should have no t-stats."""
        comp = self._make_comparator()

        entries_a = [_make_entry(delta_G_per_replicate=[-1.2, -1.4, -1.1], n_replicates=3)]
        entries_b = [_make_entry(delta_G_per_replicate=[-0.5, -0.6, -0.7], n_replicates=3)]
        sa = _make_condition_summary(label="A", temperature_K=300.0, entries=entries_a)
        sb = _make_condition_summary(label="B", temperature_K=320.0, entries=entries_b)

        pairwise = comp._compare_condition_pair(sa, sb)
        assert len(pairwise) == 1
        pw = pairwise[0]
        assert pw.cross_temperature
        assert pw.t_statistic is None
        assert pw.p_value is None

    def test_entry_from_agg_bp_entry_zero_contact_share_returns_none(self):
        """Entry with contact_share=0 should return None."""
        from unittest.mock import MagicMock

        comp = self._make_comparator()
        entry = MagicMock()
        entry.mean_contact_share = 0.0
        entry.expected_share = 0.10
        entry.sem_contact_share = 0.01
        entry.per_replicate_enrichments = []
        entry.partition_element = "aromatic"

        result = comp._entry_from_agg_bp_entry(
            entry, "SBM", "aa_class", KB_KCAL * T_REF, "kcal/mol", T_REF
        )
        assert result is None

    def test_entry_from_agg_bp_entry_zero_expected_share_returns_none(self):
        """Entry with expected_share=0 should return None."""
        from unittest.mock import MagicMock

        comp = self._make_comparator()
        entry = MagicMock()
        entry.mean_contact_share = 0.20
        entry.expected_share = 0.0
        entry.sem_contact_share = 0.01
        entry.per_replicate_enrichments = []
        entry.partition_element = "aromatic"

        result = comp._entry_from_agg_bp_entry(
            entry, "SBM", "aa_class", KB_KCAL * T_REF, "kcal/mol", T_REF
        )
        assert result is None

    def test_entry_from_agg_bp_entry_per_replicate_nan_for_negative_ratio(self):
        """Enrichment_rep + 1 <= 0 should produce NaN in delta_G_per_replicate."""
        from unittest.mock import MagicMock

        comp = self._make_comparator()
        entry = MagicMock()
        entry.mean_contact_share = 0.20
        entry.expected_share = 0.10
        entry.sem_contact_share = 0.0
        entry.per_replicate_enrichments = [-1.5]  # ratio = -0.5 → NaN
        entry.partition_element = "aromatic"
        entry.n_exposed_in_group = 0

        result = comp._entry_from_agg_bp_entry(
            entry, "SBM", "aa_class", KB_KCAL * T_REF, "kcal/mol", T_REF
        )
        assert result is not None
        assert math.isnan(result.delta_G_per_replicate[0])


# ---------------------------------------------------------------------------
# Formatter Tests
# ---------------------------------------------------------------------------


class TestBindingFreeEnergyFormatters:
    """Test output formatters."""

    def _build_result(self, n_conditions: int = 2):
        from polyzymd.compare.results.binding_free_energy import BindingFreeEnergyResult

        conditions = []
        for i in range(n_conditions):
            entries = [
                _make_entry(
                    polymer_type="SBM",
                    protein_group="aromatic",
                    contact_share=0.20 + 0.05 * i,
                    expected_share=0.10,
                    n_replicates=3,
                    delta_G_per_replicate=[-1.0 - 0.1 * i, -0.9 - 0.1 * i, -1.1 - 0.1 * i],
                ),
                _make_entry(
                    polymer_type="SBM",
                    protein_group="polar",
                    contact_share=0.10,
                    expected_share=0.20,
                    n_replicates=3,
                    delta_G_per_replicate=[0.3 + 0.05 * i, 0.35 + 0.05 * i, 0.28 + 0.05 * i],
                ),
            ]
            conditions.append(_make_condition_summary(label=f"Cond {i}", entries=entries))

        from polyzymd.compare.results.binding_free_energy import FreeEnergyPairwiseEntry

        pairwise = []
        if n_conditions >= 2:
            pairwise.append(
                FreeEnergyPairwiseEntry(
                    polymer_type="SBM",
                    protein_group="aromatic",
                    condition_a="Cond 0",
                    condition_b="Cond 1",
                    temperature_a_K=T_REF,
                    temperature_b_K=T_REF,
                    cross_temperature=False,
                    delta_G_a=-1.0,
                    delta_G_b=-1.1,
                    delta_delta_G=-0.1,
                    t_statistic=-2.5,
                    p_value=0.06,
                )
            )

        return BindingFreeEnergyResult(
            name="test",
            units="kT",
            conditions=conditions,
            pairwise_comparisons=pairwise,
            polymer_types=["SBM"],
            protein_groups=["aromatic", "polar"],
            equilibration_time="10ns",
        )

    def test_json_format_valid(self):
        from polyzymd.compare.binding_free_energy_formatters import format_bfe_result

        result = self._build_result()
        output = format_bfe_result(result, format="json")
        data = json.loads(output)
        assert data["name"] == "test"
        assert data["units"] == "kT"

    def test_table_format_contains_condition_labels(self):
        from polyzymd.compare.binding_free_energy_formatters import format_bfe_result

        result = self._build_result()
        output = format_bfe_result(result, format="table")
        assert "Cond 0" in output
        assert "Cond 1" in output

    def test_table_format_contains_polymer_type(self):
        from polyzymd.compare.binding_free_energy_formatters import format_bfe_result

        result = self._build_result()
        output = format_bfe_result(result, format="table")
        assert "SBM" in output

    def test_table_format_contains_ddg_symbol(self):
        from polyzymd.compare.binding_free_energy_formatters import format_bfe_result

        result = self._build_result()
        output = format_bfe_result(result, format="table")
        # The output should mention ΔΔG or dG somewhere
        assert any(sym in output for sym in ["ΔΔG", "dG", "delta_G", "DDG", "ΔG"])

    def test_markdown_format_has_headers(self):
        from polyzymd.compare.binding_free_energy_formatters import format_bfe_result

        result = self._build_result()
        output = format_bfe_result(result, format="markdown")
        assert "#" in output  # Markdown headers

    def test_markdown_format_has_table_pipes(self):
        from polyzymd.compare.binding_free_energy_formatters import format_bfe_result

        result = self._build_result()
        output = format_bfe_result(result, format="markdown")
        assert "|" in output  # Markdown table delimiters

    def test_format_single_condition(self):
        """Single condition should not raise and should produce output."""
        from polyzymd.compare.binding_free_energy_formatters import format_bfe_result

        result = self._build_result(n_conditions=1)
        for fmt in ("table", "markdown", "json"):
            output = format_bfe_result(result, format=fmt)
            assert len(output) > 0


# ---------------------------------------------------------------------------
# Registry / Import Tests
# ---------------------------------------------------------------------------


class TestRegistration:
    """Test that BFE types are registered properly."""

    def test_comparator_registered(self):
        from polyzymd.compare.core.registry import ComparatorRegistry

        assert ComparatorRegistry.is_registered("binding_free_energy")

    def test_analysis_settings_registered(self):
        from polyzymd.analysis.core.registry import AnalysisSettingsRegistry

        assert AnalysisSettingsRegistry.is_registered("binding_free_energy")

    def test_comparison_settings_registered(self):
        from polyzymd.analysis.core.registry import ComparisonSettingsRegistry

        assert ComparisonSettingsRegistry.is_registered("binding_free_energy")

    def test_comparators_init_exports(self):
        from polyzymd.compare.comparators import BindingFreeEnergyComparator

        assert BindingFreeEnergyComparator is not None

    def test_results_module_importable(self):
        from polyzymd.compare.results.binding_free_energy import (
            BindingFreeEnergyResult,
            FreeEnergyConditionSummary,
            FreeEnergyEntry,
            FreeEnergyPairwiseEntry,
        )

        assert BindingFreeEnergyResult is not None
        assert FreeEnergyConditionSummary is not None
        assert FreeEnergyEntry is not None
        assert FreeEnergyPairwiseEntry is not None

    def test_formatters_importable(self):
        from polyzymd.compare.binding_free_energy_formatters import (
            format_bfe_result,
        )

        assert format_bfe_result is not None
