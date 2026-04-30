"""Tests for BindingFreeEnergyAnalysis plugin — Phase C-6.

Coverage:
- Discovery (registered by name and alias)
- Class attributes (name, Settings, aliases, dependencies, min_replicates)
- Settings validation (units, threshold, k_b)
- run_replicate / aggregate (no-op)
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
from dataclasses import dataclass, replace
from pathlib import Path
from types import SimpleNamespace
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

        with pytest.raises(ValueError):
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

    def test_binding_preference_fingerprint_ignores_bfe_only_settings(self):
        from polyzymd.analyses.binding_free_energy import BFESettings
        from polyzymd.analyses.shared.binding_preference_helpers import (
            binding_preference_settings_fingerprint,
        )

        default_units = BFESettings(units="kT", fdr_alpha=0.05)
        changed_bfe_only = BFESettings(units="kcal/mol", fdr_alpha=0.01)

        assert binding_preference_settings_fingerprint(default_units) == (
            binding_preference_settings_fingerprint(changed_bfe_only)
        )

    def test_binding_preference_fingerprint_changes_with_bp_settings(self):
        from polyzymd.analyses.binding_free_energy import BFESettings
        from polyzymd.analyses.shared.binding_preference_helpers import (
            binding_preference_settings_fingerprint,
        )

        default = BFESettings(surface_exposure_threshold=0.2)
        changed = BFESettings(surface_exposure_threshold=0.35)

        assert binding_preference_settings_fingerprint(default) != (
            binding_preference_settings_fingerprint(changed)
        )


class TestBindingPreferenceCacheIdentity:
    """Regression tests for BP-specific cache identity."""

    def test_bp_cache_rejects_mismatched_bp_settings_fingerprint(self, tmp_path):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.binding_free_energy import BFESettings
        from polyzymd.analyses.shared.binding_preference import BindingPreferenceResult
        from polyzymd.analyses.shared.binding_preference_helpers import (
            binding_preference_settings_fingerprint,
            try_load_cached_binding_preference,
        )

        matching_fp = binding_preference_settings_fingerprint(
            BFESettings(surface_exposure_threshold=0.2)
        )
        mismatched_fp = binding_preference_settings_fingerprint(
            BFESettings(surface_exposure_threshold=0.35)
        )
        contact_fp = "deadbeef"
        result = BindingPreferenceResult(
            n_frames=10,
            metadata={
                "binding_preference_settings_fingerprint": matching_fp,
                "contacts_settings_fingerprint": contact_fp,
                "replicate": 1,
                "equilibration": "10ns",
            },
        )
        result.save(tmp_path / "binding_preference_rep1.json")
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=MagicMock(),
        )

        loaded = try_load_cached_binding_preference(
            cond,
            tmp_path,
            settings_fp=matching_fp,
            contact_settings_fp=contact_fp,
            equilibration="10ns",
            successful_replicates=(1,),
        )
        rejected = try_load_cached_binding_preference(
            cond,
            tmp_path,
            settings_fp=mismatched_fp,
            contact_settings_fp=contact_fp,
            equilibration="10ns",
            successful_replicates=(1,),
        )

        assert loaded is not None
        assert rejected is None

    def test_contact_result_replicate_mismatch_is_rejected(self, tmp_path):
        from polyzymd.analyses.shared.binding_preference._orchestration import (
            BindingPreferenceContactIdentityError,
            _validate_contact_result_contract,
        )

        contact_result = SimpleNamespace(
            replicate=2,
            metadata={"settings_fingerprint": "deadbeef"},
        )

        with pytest.raises(BindingPreferenceContactIdentityError, match="replicate mismatch"):
            _validate_contact_result_contract(
                contact_result,
                tmp_path / "contacts_rep1.json",
                expected_replicate=1,
                contact_settings_fp="deadbeef",
            )

    def test_contact_result_settings_mismatch_is_rejected(self, tmp_path):
        from polyzymd.analyses.shared.binding_preference._orchestration import (
            BindingPreferenceContactIdentityError,
            _validate_contact_result_contract,
        )

        contact_result = SimpleNamespace(
            replicate=1,
            metadata={"settings_fingerprint": "badcafe0"},
        )

        with pytest.raises(
            BindingPreferenceContactIdentityError,
            match="settings fingerprint mismatch",
        ):
            _validate_contact_result_contract(
                contact_result,
                tmp_path / "contacts_rep1.json",
                expected_replicate=1,
                contact_settings_fp="deadbeef",
            )

    def test_contact_result_filename_only_fingerprint_is_rejected(self, tmp_path):
        from polyzymd.analyses.shared.binding_preference._orchestration import (
            BindingPreferenceContactIdentityError,
            _validate_contact_result_contract,
        )

        contact_path = tmp_path / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json"
        contact_path.write_text("{}")
        contact_result = SimpleNamespace(replicate=1, metadata={})

        with pytest.raises(
            BindingPreferenceContactIdentityError,
            match="settings fingerprint mismatch",
        ):
            _validate_contact_result_contract(
                contact_result,
                contact_path,
                expected_replicate=1,
                contact_settings_fp="deadbeef",
            )


# ===========================================================================
# run_replicate / aggregate (no-op)
# ===========================================================================


class TestNoOp:
    def test_run_replicate_returns_none(self):
        """Compare-only plugins must no-op compute stage by returning ``None``.

        BindingFreeEnergyAnalysis sets ``has_compute_stage = False``, so the
        inherited Analysis.run_replicate contract is to return ``None``.
        """
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
        assert analysis.run_replicate(ctx, 1) is None

    def test_aggregate_returns_none(self):
        """Compare-only plugins must no-op aggregate stage by returning ``None``.

        BindingFreeEnergyAnalysis sets ``has_aggregate_stage = False``, so the
        inherited Analysis.aggregate contract is to return ``None``.
        """
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

    def test_compare_is_overridden(self):
        """BindingFreeEnergyAnalysis must override Analysis.compare.

        BFE performs custom compare-only logic, so it cannot use the base
        scalar comparison implementation from Analysis.
        """
        from polyzymd.analyses.base import Analysis
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        assert BindingFreeEnergyAnalysis.compare is not Analysis.compare
        assert "compare" in BindingFreeEnergyAnalysis.__dict__


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


@dataclass
class _MockBPEntry:
    partition_element: str
    mean_contact_share: float
    expected_share: float
    sem_contact_share: float
    per_replicate_enrichments: list[float]
    n_exposed_in_group: int


@dataclass
class _MockPartitionResult:
    entries: list[_MockBPEntry]


@dataclass
class _MockBindingPreference:
    aa_class_binding: dict[str, _MockPartitionResult]
    user_defined_partitions: dict[str, dict[str, _MockPartitionResult]]


@dataclass
class _MockBPResult:
    binding_preference: _MockBindingPreference


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

    return _MockBPEntry(
        partition_element=partition_element,
        mean_contact_share=mean_contact_share,
        expected_share=expected_share,
        sem_contact_share=sem_contact_share,
        per_replicate_enrichments=per_replicate_enrichments,
        n_exposed_in_group=10,
    )


def _make_mock_bp_result(entries=None, user_partitions=None):
    """Build a mock AggregatedBindingPreferenceResult with minimal structure."""
    if entries is None:
        entries = [_make_mock_bp_entry()]

    partition_result = _MockPartitionResult(entries=entries)
    bp = _MockBindingPreference(
        aa_class_binding={"SBM": partition_result},
        user_defined_partitions=user_partitions or {},
    )
    return _MockBPResult(binding_preference=bp)


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

        # All conditions got empty entries, so compare should return None
        assert result is None

    def test_compare_returns_none_all_summaries_empty(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis
        from polyzymd.analyses.binding_free_energy._comparison_results import (
            FreeEnergyConditionSummary,
        )

        analysis = BindingFreeEnergyAnalysis()
        ctx = _make_context(n_conditions=2)
        empty_summary = FreeEnergyConditionSummary(
            label="A",
            config_path="/tmp/A/config.yaml",
            temperature_K=300.0,
            n_replicates=0,
            units="kT",
            entries=[],
            polymer_types=[],
            protein_groups=[],
        )
        empty_summary_2 = empty_summary.model_copy(update={"label": "B"})

        with patch.object(
            analysis,
            "_build_condition_summary",
            side_effect=[empty_summary, empty_summary_2],
        ):
            result = analysis.compare(ctx)

        assert result is None

    def test_compare_propagates_unexpected_summary_runtime_error(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        ctx = _make_context(n_conditions=1)

        with patch.object(
            analysis,
            "_build_condition_summary",
            side_effect=RuntimeError("unexpected summary failure"),
        ):
            with pytest.raises(RuntimeError, match="unexpected summary failure"):
                analysis.compare(ctx)

    def test_compare_mixed_temperatures(self):
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.binding_free_energy import (
            BFESettings,
            BindingFreeEnergyAnalysis,
        )
        from polyzymd.analyses.binding_free_energy._comparison_results import (
            FreeEnergyConditionSummary,
            FreeEnergyEntry,
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

        summary_a = FreeEnergyConditionSummary(
            label="A",
            config_path="/tmp/A/config.yaml",
            temperature_K=300.0,
            n_replicates=3,
            units="kT",
            entries=[
                FreeEnergyEntry(
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
            ],
            polymer_types=["SBM"],
            protein_groups=["aromatic"],
        )
        summary_b = summary_a.model_copy(
            update={
                "label": "B",
                "config_path": "/tmp/B/config.yaml",
                "temperature_K": 350.0,
                "entries": [
                    summary_a.entries[0].model_copy(
                        update={
                            "delta_G": -0.6,
                            "delta_G_per_replicate": [-0.5, -0.6, -0.7],
                            "temperature_K": 350.0,
                        }
                    )
                ],
            }
        )

        with patch.object(
            analysis,
            "_build_condition_summary",
            side_effect=[summary_a, summary_b],
        ):
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

    def test_pairwise_keeps_distinct_partition_identity(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis
        from polyzymd.analyses.binding_free_energy._comparison_results import (
            FreeEnergyConditionSummary,
            FreeEnergyEntry,
        )

        analysis = BindingFreeEnergyAnalysis()

        entries_a = [
            FreeEnergyEntry(
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
            ),
            FreeEnergyEntry(
                polymer_type="SBM",
                protein_group="aromatic",
                partition_name="custom_partition",
                contact_share=0.25,
                expected_share=0.1,
                enrichment_ratio=2.5,
                delta_G=-0.7,
                delta_G_uncertainty=0.05,
                delta_G_per_replicate=[-0.6, -0.7, -0.8],
                units="kT",
                temperature_K=300.0,
                n_replicates=3,
                n_exposed_in_group=10,
            ),
        ]
        entries_b = [
            entries_a[0].model_copy(
                update={"delta_G": -0.6, "delta_G_per_replicate": [-0.5, -0.6, -0.7]}
            ),
            entries_a[1].model_copy(
                update={"delta_G": -0.9, "delta_G_per_replicate": [-0.8, -0.9, -1.0]}
            ),
        ]

        summary_a = FreeEnergyConditionSummary(
            label="A",
            config_path="/tmp/A/config.yaml",
            temperature_K=300.0,
            n_replicates=3,
            units="kT",
            entries=entries_a,
            polymer_types=["SBM"],
            protein_groups=["aromatic"],
        )
        summary_b = FreeEnergyConditionSummary(
            label="B",
            config_path="/tmp/B/config.yaml",
            temperature_K=300.0,
            n_replicates=3,
            units="kT",
            entries=entries_b,
            polymer_types=["SBM"],
            protein_groups=["aromatic"],
        )

        pairs = analysis._compare_condition_pair(summary_a, summary_b)

        assert len(pairs) == 2
        partition_names = {pair.partition_name for pair in pairs}
        assert partition_names == {"aa_class", "custom_partition"}


class TestPairwiseSerialization:
    def test_free_energy_pairwise_entry_partition_roundtrip(self):
        from polyzymd.analyses.binding_free_energy._comparison_results import (
            FreeEnergyPairwiseEntry,
        )

        entry = FreeEnergyPairwiseEntry(
            polymer_type="SBM",
            protein_group="aromatic",
            partition_name="custom_partition",
            condition_a="A",
            condition_b="B",
            temperature_a_K=300.0,
            temperature_b_K=300.0,
            cross_temperature=False,
            delta_G_a=-0.5,
            delta_G_b=-0.7,
            delta_delta_G=-0.2,
            t_statistic=1.2,
            p_value=0.05,
        )

        loaded = FreeEnergyPairwiseEntry.model_validate_json(entry.model_dump_json())
        assert loaded.partition_name == "custom_partition"


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

    def test_plot_propagates_exceptions(self, tmp_path):
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
            with pytest.raises(RuntimeError, match="boom"):
                analysis.plot(ctx)

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
    def test_class_variables_match_compare_only_contract(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        assert BindingFreeEnergyAnalysis.name == "binding_free_energy"
        assert BindingFreeEnergyAnalysis.aliases == ("bfe",)
        assert BindingFreeEnergyAnalysis.dependencies == ("contacts",)
        assert BindingFreeEnergyAnalysis.min_replicates == 1
        assert BindingFreeEnergyAnalysis.has_compute_stage is False
        assert BindingFreeEnergyAnalysis.has_aggregate_stage is False

    def test_is_subclass_of_analysis(self):
        from polyzymd.analyses.base import Analysis
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        assert issubclass(BindingFreeEnergyAnalysis, Analysis)

    def test_extract_metrics_returns_empty_for_custom_compare(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        assert analysis.extract_metrics(MagicMock()) == {}

    def test_repr(self):
        from polyzymd.analyses.binding_free_energy import BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        assert "binding_free_energy" in repr(analysis)


# ===========================================================================
# Binding preference loading
# ===========================================================================


class TestLoadBindingPreference:
    def test_non_cached_path_uses_real_condition(self, tmp_path):
        """Non-cached compute path should pass the full Condition object."""
        from polyzymd.analyses.binding_free_energy import (
            BFESettings,
            BindingFreeEnergyAnalysis,
        )
        from polyzymd.analyses.shared.binding_preference_helpers import (
            binding_preference_settings_fingerprint,
        )

        analysis = BindingFreeEnergyAnalysis()
        ctx = _make_context(n_conditions=1)
        cond = ctx.conditions[0]
        ctx = replace(ctx, recompute=True)

        settings = BFESettings(compute_binding_preference=True)
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
            ) as mock_find_contacts,
            patch(
                "polyzymd.analyses.shared.binding_preference.compute_condition_binding_preference",
                return_value=expected,
            ) as mock_compute,
        ):
            loaded = analysis._load_binding_preference(cond, ctx, settings)

        assert loaded is expected
        assert mock_compute.call_args.kwargs["cond"] is cond
        assert "settings_fp" not in mock_find_contacts.call_args.kwargs
        assert mock_compute.call_args.kwargs[
            "settings_fp"
        ] == binding_preference_settings_fingerprint(settings)

    def test_non_cached_path_discovers_canonical_contacts_results(self, tmp_path):
        """Binding free energy should consume canonical contacts run paths."""
        from polyzymd.analyses.binding_free_energy import BFESettings, BindingFreeEnergyAnalysis
        from polyzymd.analyses.contacts._results import ContactResult

        analysis = BindingFreeEnergyAnalysis()
        ctx = _make_context(n_conditions=1, n_reps=2)
        cond = ctx.conditions[0]
        contacts_dir = tmp_path / cond.label / "contacts"
        run_1 = contacts_dir / "run_1"
        run_2 = contacts_dir / "run_2"
        run_1.mkdir(parents=True)
        run_2.mkdir(parents=True)
        ContactResult(
            replicate=1,
            residue_contacts=[],
            n_frames=10,
            timestep_ps=10.0,
            criteria_label="any_atom_4.5A",
            criteria_cutoff=4.5,
            equilibration_time=10.0,
            equilibration_unit="ns",
        ).save(run_1 / "result.json")
        ContactResult(
            replicate=2,
            residue_contacts=[],
            n_frames=10,
            timestep_ps=10.0,
            criteria_label="any_atom_4.5A",
            criteria_cutoff=4.5,
            equilibration_time=10.0,
            equilibration_unit="ns",
        ).save(run_2 / "result.json")
        ctx = replace(
            ctx,
            recompute=True,
            analysis_dirs={cond.label: tmp_path / cond.label / "binding_free_energy"},
        )

        settings = BFESettings(compute_binding_preference=True)
        enzyme_pdb = tmp_path / "enzyme.pdb"
        enzyme_pdb.write_text("ATOM\n")

        with (
            patch(
                "polyzymd.analyses.shared.binding_preference_helpers.resolve_enzyme_pdb",
                return_value=enzyme_pdb,
            ),
            patch(
                "polyzymd.analyses.shared.binding_preference.compute_condition_binding_preference",
                return_value=object(),
            ) as mock_compute,
        ):
            analysis._load_binding_preference(cond, ctx, settings)

        assert mock_compute.call_args.kwargs["contact_results_by_replicate"] == {
            1: run_1 / "result.json",
            2: run_2 / "result.json",
        }

    def test_non_cached_path_discovers_fingerprinted_contacts_without_canonical(
        self,
        tmp_path,
    ):
        """Binding free energy should consume matching fingerprinted contacts sidecars."""

        from polyzymd.analyses.binding_free_energy import BFESettings, BindingFreeEnergyAnalysis
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint
        from polyzymd.analyses.shared.binding_preference_helpers import (
            binding_preference_settings_fingerprint,
        )
        from polyzymd.analyses.shared.config_hash import settings_fingerprint

        analysis = BindingFreeEnergyAnalysis()
        ctx = _make_context(n_conditions=1, n_reps=1)
        cond = ctx.conditions[0]
        contacts_dir = tmp_path / cond.label / "contacts"
        contacts_dir.mkdir(parents=True)
        ctx = replace(
            ctx,
            recompute=True,
            analysis_dirs={cond.label: tmp_path / cond.label / "binding_free_energy"},
        )

        settings = BFESettings(compute_binding_preference=True, units="kcal/mol")
        contacts_settings = ContactsSettings(
            compute_binding_preference=True,
            cutoff=6.0,
            grouping="none",
        )
        fp = contacts_settings_fingerprint(contacts_settings)
        assert fp != settings_fingerprint(settings)
        sidecar = contacts_dir / f"contacts_eq10ns_cut6.0_s{fp}_rep1.json"
        sidecar.write_text(f'{{"contacts_settings_fingerprint": "{fp}"}}')
        enzyme_pdb = tmp_path / "enzyme.pdb"
        enzyme_pdb.write_text("ATOM\n")

        with (
            patch(
                "polyzymd.analyses.shared.binding_preference_helpers.resolve_enzyme_pdb",
                return_value=enzyme_pdb,
            ),
            patch(
                "polyzymd.analyses.shared.binding_preference.compute_condition_binding_preference",
                return_value=object(),
            ) as mock_compute,
        ):
            analysis._load_binding_preference(cond, ctx, settings)

        assert mock_compute.call_args.kwargs["contact_results_by_replicate"] == {1: sidecar}
        assert mock_compute.call_args.kwargs["contact_settings_fp"] == fp
        assert mock_compute.call_args.kwargs[
            "settings_fp"
        ] == binding_preference_settings_fingerprint(settings)
        assert mock_compute.call_args.kwargs["settings_fp"] != fp

    def test_non_cached_path_uses_current_sidecar_over_stale_canonical(self, tmp_path):
        """Fingerprint inference should let current sidecars override stale canonical output."""

        from polyzymd.analyses.binding_free_energy import BFESettings, BindingFreeEnergyAnalysis
        from polyzymd.analyses.contacts._results import ContactResult

        analysis = BindingFreeEnergyAnalysis()
        ctx = _make_context(n_conditions=1, n_reps=1)
        cond = ctx.conditions[0]
        contacts_dir = tmp_path / cond.label / "contacts"
        run_dir = contacts_dir / "run_1"
        run_dir.mkdir(parents=True)
        current_fp = "deadbeef"
        stale_fp = "badcafe0"
        canonical = run_dir / "result.json"
        ContactResult(
            replicate=1,
            residue_contacts=[],
            n_frames=10,
            timestep_ps=10.0,
            criteria_label="any_atom_4.5A",
            criteria_cutoff=4.5,
            equilibration_time=10.0,
            equilibration_unit="ns",
            metadata={"contacts_settings_fingerprint": stale_fp},
        ).save(canonical)
        ctx = replace(
            ctx,
            recompute=True,
            analysis_dirs={cond.label: tmp_path / cond.label / "binding_free_energy"},
        )

        settings = BFESettings(compute_binding_preference=True, units="kcal/mol")
        sidecar = contacts_dir / f"contacts_eq10ns_cut6.0_s{current_fp}_rep1.json"
        sidecar.write_text(f'{{"contacts_settings_fingerprint": "{current_fp}"}}')
        enzyme_pdb = tmp_path / "enzyme.pdb"
        enzyme_pdb.write_text("ATOM\n")

        with (
            patch(
                "polyzymd.analyses.shared.binding_preference_helpers.resolve_enzyme_pdb",
                return_value=enzyme_pdb,
            ),
            patch(
                "polyzymd.analyses.shared.binding_preference.compute_condition_binding_preference",
                return_value=object(),
            ) as mock_compute,
        ):
            analysis._load_binding_preference(cond, ctx, settings)

        assert mock_compute.call_args.kwargs["contact_results_by_replicate"] == {1: sidecar}
        assert mock_compute.call_args.kwargs["contact_settings_fp"] == current_fp
        assert mock_compute.call_args.kwargs["settings_fp"] != current_fp

    def test_non_cached_path_preserves_metadata_proven_legacy_contacts_sidecar(
        self,
        tmp_path,
    ):
        """Identity inference should re-resolve proven legacy contacts sidecars."""

        from polyzymd.analyses.binding_free_energy import BFESettings, BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        ctx = _make_context(n_conditions=1, n_reps=1)
        cond = ctx.conditions[0]
        contacts_dir = tmp_path / cond.label / "contacts"
        contacts_dir.mkdir(parents=True)
        ctx = replace(
            ctx,
            recompute=True,
            analysis_dirs={cond.label: tmp_path / cond.label / "binding_free_energy"},
        )

        contacts_fp = "deadbeef"
        legacy = contacts_dir / "contacts_eq10ns_cut6.0_rep1.json"
        legacy.write_text(
            f'{{"metadata": {{"contacts_settings_fingerprint": "{contacts_fp}", '
            '"equilibration": "10ns"}}'
        )
        enzyme_pdb = tmp_path / "enzyme.pdb"
        enzyme_pdb.write_text("ATOM\n")

        with (
            patch(
                "polyzymd.analyses.shared.binding_preference_helpers.resolve_enzyme_pdb",
                return_value=enzyme_pdb,
            ),
            patch(
                "polyzymd.analyses.shared.binding_preference.compute_condition_binding_preference",
                return_value=object(),
            ) as mock_compute,
        ):
            analysis._load_binding_preference(
                cond,
                ctx,
                BFESettings(compute_binding_preference=True),
            )

        assert mock_compute.call_args.kwargs["contact_results_by_replicate"] == {1: legacy}
        assert mock_compute.call_args.kwargs["contact_settings_fp"] == contacts_fp

    def test_cached_lookup_uses_actual_non_default_contacts_fingerprint(self, tmp_path):
        """Cached BP lookup should use the upstream contacts artifact fingerprint."""

        from polyzymd.analyses.binding_free_energy import BFESettings, BindingFreeEnergyAnalysis
        from polyzymd.analyses.contacts import ContactsSettings
        from polyzymd.analyses.contacts._identity import contacts_settings_fingerprint

        analysis = BindingFreeEnergyAnalysis()
        ctx = _make_context(n_conditions=1, n_reps=1)
        cond = ctx.conditions[0]
        contacts_dir = tmp_path / cond.label / "contacts"
        contacts_dir.mkdir(parents=True)
        ctx = replace(
            ctx,
            recompute=False,
            analysis_dirs={cond.label: tmp_path / cond.label / "binding_free_energy"},
        )
        contacts_settings = ContactsSettings(
            compute_binding_preference=True,
            cutoff=6.0,
            grouping="none",
        )
        fp = contacts_settings_fingerprint(contacts_settings)
        (contacts_dir / f"contacts_eq10ns_cut6.0_s{fp}_rep1.json").write_text(
            f'{{"contacts_settings_fingerprint": "{fp}"}}'
        )
        expected = object()

        with patch(
            "polyzymd.analyses.shared.binding_preference_helpers."
            "try_load_cached_binding_preference",
            return_value=expected,
        ) as mock_cached:
            loaded = analysis._load_binding_preference(
                cond,
                ctx,
                BFESettings(compute_binding_preference=True),
            )

        assert loaded is expected
        assert mock_cached.call_args.kwargs["contact_settings_fp"] == fp
        assert mock_cached.call_args.kwargs["settings_fp"] != fp
        assert mock_cached.call_args.kwargs["equilibration"] == "10ns"
        assert mock_cached.call_args.kwargs["successful_replicates"] == (1,)

    def test_cached_lookup_ignores_stale_bp_fingerprint_when_contacts_identify_current(
        self,
        tmp_path,
    ):
        """Stale BP sidecars must not poison contacts fingerprint inference."""

        from polyzymd.analyses.binding_free_energy import BFESettings, BindingFreeEnergyAnalysis
        from polyzymd.analyses.shared.binding_preference import AggregatedBindingPreferenceResult

        analysis = BindingFreeEnergyAnalysis()
        ctx = _make_context(n_conditions=1, n_reps=1)
        cond = ctx.conditions[0]
        contacts_dir = tmp_path / cond.label / "contacts"
        contacts_dir.mkdir(parents=True)
        ctx = replace(
            ctx,
            recompute=False,
            analysis_dirs={cond.label: tmp_path / cond.label / "binding_free_energy"},
        )

        current_fp = "deadbeef"
        stale_fp = "badcafe0"
        run_dir = contacts_dir / "run_1"
        run_dir.mkdir()
        (run_dir / "result.json").write_text(
            f'{{"metadata": {{"contacts_settings_fingerprint": "{stale_fp}", '
            '"equilibration": "10ns"}}}}'
        )
        (contacts_dir / f"contacts_eq10ns_cut6.0_s{current_fp}_rep1.json").write_text(
            f'{{"contacts_settings_fingerprint": "{current_fp}"}}'
        )
        AggregatedBindingPreferenceResult(
            n_replicates=1,
            metadata={"settings_fingerprint": stale_fp, "equilibration": "10ns"},
        ).save(contacts_dir / f"binding_preference_aggregated_s{stale_fp}_reps1-1.json")
        expected = object()

        with patch(
            "polyzymd.analyses.shared.binding_preference_helpers."
            "try_load_cached_binding_preference",
            return_value=expected,
        ) as mock_cached:
            loaded = analysis._load_binding_preference(
                cond,
                ctx,
                BFESettings(compute_binding_preference=True),
            )

        assert loaded is expected
        assert mock_cached.call_args.kwargs["contact_settings_fp"] == current_fp
        assert mock_cached.call_args.kwargs["settings_fp"] != current_fp

    def test_filename_only_contacts_sidecar_blocks_cached_bp_load(self, tmp_path):
        """Unproven contacts sidecars should prevent downstream BP cache fallback."""
        from polyzymd.analyses.binding_free_energy import BFESettings, BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        ctx = _make_context(n_conditions=1, n_reps=1)
        cond = ctx.conditions[0]
        contacts_dir = tmp_path / cond.label / "contacts"
        contacts_dir.mkdir(parents=True)
        (contacts_dir / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json").write_text("{}")
        ctx = replace(
            ctx,
            recompute=False,
            analysis_dirs={cond.label: tmp_path / cond.label / "binding_free_energy"},
        )

        with patch(
            "polyzymd.analyses.shared.binding_preference_helpers."
            "try_load_cached_binding_preference",
            return_value=object(),
        ) as mock_cached:
            loaded = analysis._load_binding_preference(
                cond,
                ctx,
                BFESettings(compute_binding_preference=True),
            )

        assert loaded is None
        mock_cached.assert_not_called()

    def test_filename_only_contacts_sidecar_blocks_bp_recompute(self, tmp_path):
        """Unproven contacts sidecars should prevent recompute from contacts."""
        from polyzymd.analyses.binding_free_energy import BFESettings, BindingFreeEnergyAnalysis

        analysis = BindingFreeEnergyAnalysis()
        ctx = _make_context(n_conditions=1, n_reps=1)
        cond = ctx.conditions[0]
        contacts_dir = tmp_path / cond.label / "contacts"
        contacts_dir.mkdir(parents=True)
        (contacts_dir / "contacts_eq10ns_cut4.5_sdeadbeef_rep1.json").write_text("{}")
        ctx = replace(
            ctx,
            recompute=True,
            analysis_dirs={cond.label: tmp_path / cond.label / "binding_free_energy"},
        )

        with patch(
            "polyzymd.analyses.shared.binding_preference.compute_condition_binding_preference",
            return_value=object(),
        ) as mock_compute:
            loaded = analysis._load_binding_preference(
                cond,
                ctx,
                BFESettings(compute_binding_preference=True),
            )

        assert loaded is None
        mock_compute.assert_not_called()


class TestBindingPreferenceCacheHelper:
    """Regression tests for settings-aware binding-preference cache loading."""

    def test_rejects_legacy_canonical_cache_when_settings_fp_known(self, tmp_path):
        from polyzymd.analyses.shared.binding_preference import AggregatedBindingPreferenceResult
        from polyzymd.analyses.shared.binding_preference_helpers import (
            try_load_cached_binding_preference,
        )

        cond = SimpleNamespace(label="Cond", replicates=(1, 2))
        AggregatedBindingPreferenceResult(n_replicates=2).save(
            tmp_path / "binding_preference_aggregated.json"
        )

        result = try_load_cached_binding_preference(cond, tmp_path, settings_fp="deadbeef")

        assert result is None

    def test_accepts_fingerprinted_cache_when_settings_fp_known(self, tmp_path):
        from polyzymd.analyses.shared.binding_preference import AggregatedBindingPreferenceResult
        from polyzymd.analyses.shared.binding_preference_helpers import (
            try_load_cached_binding_preference,
        )

        cond = SimpleNamespace(label="Cond", replicates=(1, 2))
        AggregatedBindingPreferenceResult(
            n_replicates=2,
            metadata={"replicates": [1, 2]},
        ).save(tmp_path / "binding_preference_aggregated_sdeadbeef.json")

        result = try_load_cached_binding_preference(cond, tmp_path, settings_fp="deadbeef")

        assert result is not None
        assert result.n_replicates == 2

    def test_rejects_fingerprinted_aggregate_with_wrong_successful_subset(self, tmp_path):
        """Aggregate BP caches must prove the requested successful replicate set."""

        from polyzymd.analyses.shared.binding_preference import AggregatedBindingPreferenceResult
        from polyzymd.analyses.shared.binding_preference_helpers import (
            try_load_cached_binding_preference,
        )

        cond = SimpleNamespace(label="Cond", replicates=(1, 2, 3))
        AggregatedBindingPreferenceResult(
            n_replicates=2,
            metadata={"replicates": [1, 2], "settings_fingerprint": "deadbeef"},
        ).save(tmp_path / "binding_preference_aggregated_sdeadbeef_reps1-2.json")

        result = try_load_cached_binding_preference(
            cond,
            tmp_path,
            settings_fp="deadbeef",
            successful_replicates=(2, 3),
        )

        assert result is None

    def test_rejects_legacy_aggregate_without_replicate_identity(self, tmp_path):
        """Legacy aggregate BP caches cannot prove which replicates they contain."""

        from polyzymd.analyses.shared.binding_preference import AggregatedBindingPreferenceResult
        from polyzymd.analyses.shared.binding_preference_helpers import (
            try_load_cached_binding_preference,
        )

        cond = SimpleNamespace(label="Cond", replicates=(1, 2))
        AggregatedBindingPreferenceResult(n_replicates=2).save(
            tmp_path / "binding_preference_aggregated_reps1-2.json"
        )

        result = try_load_cached_binding_preference(cond, tmp_path)

        assert result is None

    def test_rejects_cached_binding_preference_from_wrong_window(self, tmp_path):
        from polyzymd.analyses.shared.binding_preference import AggregatedBindingPreferenceResult
        from polyzymd.analyses.shared.binding_preference_helpers import (
            try_load_cached_binding_preference,
        )

        cond = SimpleNamespace(label="Cond", replicates=(1, 2))
        AggregatedBindingPreferenceResult(
            n_replicates=2,
            metadata={"settings_fingerprint": "deadbeef", "equilibration": "0ns"},
        ).save(tmp_path / "binding_preference_aggregated_sdeadbeef.json")

        result = try_load_cached_binding_preference(
            cond,
            tmp_path,
            settings_fp="deadbeef",
            equilibration="10ns",
        )

        assert result is None


class TestBindingPreferenceAggregationArtifacts:
    """Regression tests for binding-preference aggregate cache identity."""

    def test_aggregate_metadata_and_filename_use_successful_replicates(self, tmp_path):
        """Skipped contact replicates should be excluded from aggregate identity."""

        from polyzymd.analyses.shared.binding_preference import (
            AggregatedBindingPreferenceResult,
            BindingPreferenceResult,
            compute_condition_binding_preference,
        )

        cond = SimpleNamespace(label="Cond", replicates=(1, 2, 3))
        sim_config = MagicMock()
        sim_config.get_working_directory.return_value = tmp_path / "run_1"
        enzyme_pdb = tmp_path / "enzyme.pdb"
        enzyme_pdb.write_text("ATOM\n")

        def _compute_bp(**_kwargs):
            return BindingPreferenceResult(n_frames=10, metadata={})

        def _aggregate_bp(results):
            return AggregatedBindingPreferenceResult(n_replicates=len(results))

        with (
            patch.dict("sys.modules", {"MDAnalysis": MagicMock()}),
            patch("polyzymd.analyses.shared.loader.TrajectoryLoader") as mock_loader_cls,
            patch(
                "polyzymd.analyses.shared.surface_exposure.SurfaceExposureFilter"
            ) as mock_filter_cls,
            patch(
                "polyzymd.analyses.shared.binding_preference._orchestration."
                "resolve_protein_groups_from_surface_exposure",
                return_value={"aromatic": [1]},
            ),
            patch(
                "polyzymd.analyses.shared.binding_preference._orchestration."
                "compute_binding_preference",
                side_effect=_compute_bp,
            ),
            patch(
                "polyzymd.analyses.shared.binding_preference._aggregate."
                "aggregate_binding_preference",
                side_effect=_aggregate_bp,
            ),
        ):
            mock_loader_cls.return_value.find_topology.side_effect = FileNotFoundError
            mock_filter_cls.return_value.calculate.return_value = SimpleNamespace(
                exposed_count=1,
                total_count=1,
            )

            result = compute_condition_binding_preference(
                cond,
                sim_config,
                tmp_path,
                enzyme_pdb=enzyme_pdb,
                contact_results_by_replicate={
                    1: tmp_path / "contacts_rep1.json",
                    3: tmp_path / "contacts_rep3.json",
                },
                load_contact_result=lambda _path: SimpleNamespace(
                    metadata={"equilibration": "10ns"}
                ),
                threshold=0.2,
                include_default_aa_groups=True,
                custom_protein_groups=None,
                protein_partitions=None,
                polymer_type_selections=None,
                polymer_chain="C",
                settings_fp="deadbeef",
                equilibration="10ns",
            )

        assert result is not None
        assert result.n_replicates == 2
        assert result.metadata["replicates"] == [1, 3]
        assert (tmp_path / "binding_preference_aggregated_sdeadbeef_reps1_3.json").exists()
        assert not (tmp_path / "binding_preference_aggregated_sdeadbeef_reps1-3.json").exists()


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
