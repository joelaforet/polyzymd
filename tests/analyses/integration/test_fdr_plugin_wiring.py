"""Focused tests for FDR/effect-size/top-residue wiring in analysis plugins."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from polyzymd.analyses import get_analysis
from polyzymd.analyses.contacts._aggregator import AggregatedContactResult, AggregatedResidueStats
from polyzymd.analyses.contacts._comparison_results import (
    AggregateComparisonResult,
    ContactsANOVASummary,
    ContactsComparisonResult,
    ContactsConditionSummary,
    ContactsPairwiseComparison,
)
from polyzymd.analyses.contacts._formatters import format_contacts_console_table


def _make_contacts_result_for_formatting() -> ContactsComparisonResult:
    """Build a minimal contacts comparison result with adjusted p-values."""
    conditions = [
        ContactsConditionSummary(
            label="Control",
            config_path="/tmp/control.yaml",
            n_replicates=3,
            n_residues=100,
            coverage_mean=0.20,
            coverage_sem=0.01,
            mean_contact_fraction=0.05,
            mean_contact_fraction_sem=0.005,
        ),
        ContactsConditionSummary(
            label="Treatment",
            config_path="/tmp/treatment.yaml",
            n_replicates=3,
            n_residues=100,
            coverage_mean=0.25,
            coverage_sem=0.01,
            mean_contact_fraction=0.07,
            mean_contact_fraction_sem=0.004,
        ),
    ]
    pairwise = [
        ContactsPairwiseComparison(
            condition_a="Control",
            condition_b="Treatment",
            aggregate_comparisons=[
                AggregateComparisonResult(
                    metric="coverage",
                    condition_a="Control",
                    condition_b="Treatment",
                    condition_a_mean=0.20,
                    condition_a_sem=0.01,
                    condition_b_mean=0.25,
                    condition_b_sem=0.01,
                    t_statistic=2.5,
                    p_value=0.03,
                    p_value_adjusted=0.04,
                    cohens_d=0.8,
                    effect_size_interpretation="large",
                    significant=True,
                    meets_effect_size_threshold=True,
                    percent_change=25.0,
                    direction="increased",
                )
            ],
        )
    ]
    return ContactsComparisonResult(
        name="test_contacts",
        contacts_name="contacts",
        polymer_selection="chainID C",
        protein_selection="protein",
        cutoff=4.5,
        contact_criteria="distance_cutoff",
        conditions=conditions,
        pairwise_comparisons=pairwise,
        ranking_by_coverage=["Treatment", "Control"],
        ranking_by_contact_fraction=["Treatment", "Control"],
        equilibration_time="10ns",
        created_at=datetime.now(),
        control_label="Control",
        fdr_alpha=0.05,
        min_effect_size=0.5,
        top_residues=3,
    )


def test_contacts_fdr_correction_applied() -> None:
    """Contacts BH correction should produce exact adjusted p-values."""
    cls = get_analysis("contacts")
    analysis = cls()

    raw_p_values = [0.001, 0.01, 0.03, 0.031, 0.08, 0.2, 0.21, 0.9]
    comparisons: list[ContactsPairwiseComparison] = []

    for pair_idx in range(4):
        p_cov = raw_p_values[2 * pair_idx]
        p_cf = raw_p_values[2 * pair_idx + 1]

        comparisons.append(
            ContactsPairwiseComparison(
                condition_a=f"A{pair_idx}",
                condition_b=f"B{pair_idx}",
                aggregate_comparisons=[
                    AggregateComparisonResult(
                        metric="coverage",
                        condition_a=f"A{pair_idx}",
                        condition_b=f"B{pair_idx}",
                        condition_a_mean=0.2,
                        condition_a_sem=0.01,
                        condition_b_mean=0.3,
                        condition_b_sem=0.01,
                        t_statistic=2.0,
                        p_value=p_cov,
                        cohens_d=0.6,
                        effect_size_interpretation="medium",
                        significant=True,
                        percent_change=50.0,
                        direction="increased",
                    ),
                    AggregateComparisonResult(
                        metric="mean_contact_fraction",
                        condition_a=f"A{pair_idx}",
                        condition_b=f"B{pair_idx}",
                        condition_a_mean=0.1,
                        condition_a_sem=0.01,
                        condition_b_mean=0.11,
                        condition_b_sem=0.01,
                        t_statistic=0.8,
                        p_value=p_cf,
                        cohens_d=0.2,
                        effect_size_interpretation="small",
                        significant=True,
                        percent_change=10.0,
                        direction="increased",
                    ),
                ],
            )
        )

    analysis._apply_fdr_correction(comparisons, [], fdr_alpha=0.05)

    expected_adjusted = [0.008, 0.04, 0.062, 0.062, 0.128, 0.24, 0.24, 0.9]

    flattened = [agg for comp in comparisons for agg in comp.aggregate_comparisons]
    for idx, agg in enumerate(flattened):
        assert agg.p_value_adjusted == pytest.approx(expected_adjusted[idx])

    assert [agg.significant for agg in flattened] == [
        True,
        True,
        False,
        False,
        False,
        False,
        False,
        False,
    ]


def test_contacts_effect_size_threshold() -> None:
    """Contacts effect-size threshold should tag entries by |d| cutoff."""
    cls = get_analysis("contacts")
    analysis = cls()

    comparisons = [
        ContactsPairwiseComparison(
            condition_a="A",
            condition_b="B",
            aggregate_comparisons=[
                AggregateComparisonResult(
                    metric="coverage",
                    condition_a="A",
                    condition_b="B",
                    condition_a_mean=0.2,
                    condition_a_sem=0.01,
                    condition_b_mean=0.3,
                    condition_b_sem=0.01,
                    t_statistic=2.0,
                    p_value=0.04,
                    cohens_d=0.6,
                    effect_size_interpretation="medium",
                    significant=True,
                    percent_change=50.0,
                    direction="increased",
                ),
                AggregateComparisonResult(
                    metric="mean_contact_fraction",
                    condition_a="A",
                    condition_b="B",
                    condition_a_mean=0.1,
                    condition_a_sem=0.01,
                    condition_b_mean=0.11,
                    condition_b_sem=0.01,
                    t_statistic=0.8,
                    p_value=0.50,
                    cohens_d=0.4,
                    effect_size_interpretation="small",
                    significant=False,
                    percent_change=10.0,
                    direction="increased",
                ),
            ],
        )
    ]

    analysis._apply_effect_size_threshold(comparisons, min_effect_size=0.5)

    first, second = comparisons[0].aggregate_comparisons
    assert first.meets_effect_size_threshold is True
    assert second.meets_effect_size_threshold is False


def test_contacts_top_residues_selection() -> None:
    """Contacts top-residue helper should work with real aggregated models."""
    cls = get_analysis("contacts")
    analysis = cls()

    agg_a = AggregatedContactResult(
        residue_stats=[
            AggregatedResidueStats(
                protein_resid=10,
                protein_resname="ALA",
                protein_group="nonpolar",
                contact_fraction_mean=0.40,
                contact_fraction_per_replicate=[0.40],
            ),
            AggregatedResidueStats(
                protein_resid=11,
                protein_resname="GLY",
                protein_group="nonpolar",
                contact_fraction_mean=0.80,
                contact_fraction_per_replicate=[0.80],
            ),
            AggregatedResidueStats(
                protein_resid=12,
                protein_resname="SER",
                protein_group="polar",
                contact_fraction_mean=0.20,
                contact_fraction_per_replicate=[0.20],
            ),
            AggregatedResidueStats(
                protein_resid=13,
                protein_resname="TYR",
                protein_group="aromatic",
                contact_fraction_mean=0.60,
                contact_fraction_per_replicate=[0.60],
            ),
        ],
        n_replicates=1,
        total_frames_per_replicate=[100],
        timestep_ps=1.0,
        criteria_label="any_atom_4.5A",
        criteria_cutoff=4.5,
        coverage_mean=0.0,
        coverage_sem=0.0,
        mean_contact_fraction=0.0,
        mean_contact_fraction_sem=0.0,
    )
    agg_b = AggregatedContactResult(
        residue_stats=[
            AggregatedResidueStats(
                protein_resid=20,
                protein_resname="LEU",
                protein_group="nonpolar",
                contact_fraction_mean=0.10,
                contact_fraction_per_replicate=[0.10],
            ),
            AggregatedResidueStats(
                protein_resid=21,
                protein_resname="VAL",
                protein_group="nonpolar",
                contact_fraction_mean=0.30,
                contact_fraction_per_replicate=[0.30],
            ),
            AggregatedResidueStats(
                protein_resid=22,
                protein_resname="ASP",
                protein_group="charged_negative",
                contact_fraction_mean=0.50,
                contact_fraction_per_replicate=[0.50],
            ),
            AggregatedResidueStats(
                protein_resid=23,
                protein_resname="LYS",
                protein_group="charged_positive",
                contact_fraction_mean=0.90,
                contact_fraction_per_replicate=[0.90],
            ),
        ],
        n_replicates=1,
        total_frames_per_replicate=[100],
        timestep_ps=1.0,
        criteria_label="any_atom_4.5A",
        criteria_cutoff=4.5,
        coverage_mean=0.0,
        coverage_sem=0.0,
        mean_contact_fraction=0.0,
        mean_contact_fraction_sem=0.0,
    )

    condition_data = [
        (SimpleNamespace(label="A"), {"agg_result": agg_a}),
        (SimpleNamespace(label="B"), {"agg_result": agg_b}),
    ]

    top = analysis._compute_top_contacted_residues(condition_data, top_n=3)

    assert top is not None
    assert top["A"] == [(11, "GLY", 0.8), (13, "TYR", 0.6), (10, "ALA", 0.4)]
    assert top["B"] == [(23, "LYS", 0.9), (22, "ASP", 0.5), (21, "VAL", 0.3)]


def test_contacts_settings_wiring_to_helper_arguments() -> None:
    """Contacts settings values should flow to helper arguments via ctx.settings."""
    cls = get_analysis("contacts")
    analysis = cls()

    ctx = SimpleNamespace(
        settings=cls.Settings(
            fdr_alpha=0.10,
            min_effect_size=0.3,
            top_residues=2,
        )
    )

    comparisons = [
        ContactsPairwiseComparison(
            condition_a="Control",
            condition_b="Treatment",
            aggregate_comparisons=[
                AggregateComparisonResult(
                    metric="coverage",
                    condition_a="Control",
                    condition_b="Treatment",
                    condition_a_mean=0.2,
                    condition_a_sem=0.01,
                    condition_b_mean=0.3,
                    condition_b_sem=0.01,
                    t_statistic=2.0,
                    p_value=0.08,
                    cohens_d=0.25,
                    effect_size_interpretation="small",
                    significant=False,
                    percent_change=50.0,
                    direction="increased",
                ),
                AggregateComparisonResult(
                    metric="mean_contact_fraction",
                    condition_a="Control",
                    condition_b="Treatment",
                    condition_a_mean=0.1,
                    condition_a_sem=0.01,
                    condition_b_mean=0.11,
                    condition_b_sem=0.01,
                    t_statistic=1.1,
                    p_value=0.08,
                    cohens_d=0.35,
                    effect_size_interpretation="small",
                    significant=False,
                    percent_change=10.0,
                    direction="increased",
                ),
            ],
        )
    ]

    analysis._apply_fdr_correction(comparisons, [], ctx.settings.fdr_alpha)
    analysis._apply_effect_size_threshold(comparisons, ctx.settings.min_effect_size)

    coverage_comp, fraction_comp = comparisons[0].aggregate_comparisons
    assert coverage_comp.p_value_adjusted == pytest.approx(0.08)
    assert fraction_comp.p_value_adjusted == pytest.approx(0.08)
    assert coverage_comp.significant is True
    assert fraction_comp.significant is True
    assert coverage_comp.meets_effect_size_threshold is False
    assert fraction_comp.meets_effect_size_threshold is True

    condition_data = [
        (
            SimpleNamespace(label="Treatment"),
            {
                "agg_result": AggregatedContactResult(
                    residue_stats=[
                        AggregatedResidueStats(
                            protein_resid=1,
                            protein_resname="ALA",
                            protein_group="nonpolar",
                            contact_fraction_mean=0.9,
                            contact_fraction_per_replicate=[0.9],
                        ),
                        AggregatedResidueStats(
                            protein_resid=2,
                            protein_resname="GLY",
                            protein_group="nonpolar",
                            contact_fraction_mean=0.4,
                            contact_fraction_per_replicate=[0.4],
                        ),
                        AggregatedResidueStats(
                            protein_resid=3,
                            protein_resname="SER",
                            protein_group="polar",
                            contact_fraction_mean=0.2,
                            contact_fraction_per_replicate=[0.2],
                        ),
                    ],
                    n_replicates=1,
                    total_frames_per_replicate=[100],
                    timestep_ps=1.0,
                    criteria_label="any_atom_4.5A",
                    criteria_cutoff=4.5,
                    coverage_mean=0.0,
                    coverage_sem=0.0,
                    mean_contact_fraction=0.0,
                    mean_contact_fraction_sem=0.0,
                )
            },
        )
    ]
    top = analysis._compute_top_contacted_residues(condition_data, ctx.settings.top_residues)

    assert top is not None
    assert top["Treatment"] == [(1, "ALA", 0.9), (2, "GLY", 0.4)]


class TestContactsCompareWiring:
    """Test that fdr_alpha, min_effect_size, and top_residues flow through compare()."""

    @staticmethod
    def _make_real_agg_result(condition_idx: int) -> AggregatedContactResult:
        """Construct a real AggregatedContactResult with deterministic replicate data.

        Parameters
        ----------
        condition_idx : int
            Condition index: 0=Control, 1=TreatmentStrong, 2=TreatmentMild.

        Returns
        -------
        AggregatedContactResult
            Aggregated result instance suitable for ContactsAnalysis.compare().
        """
        if condition_idx == 0:
            residue_replicates = [
                [0.20, 0.20, 0.20],
                [0.18, 0.20, 0.22],
                [0.16, 0.18, 0.20],
                [0.14, 0.16, 0.18],
                [0.12, 0.14, 0.16],
                [0.00, 0.00, 0.20],
                [0.00, 0.00, 0.00],
                [0.00, 0.00, 0.00],
                [0.00, 0.00, 0.00],
                [0.00, 0.00, 0.00],
            ]
        elif condition_idx == 1:
            residue_replicates = [
                [0.90, 0.85, 0.90],
                [0.80, 0.75, 0.80],
                [0.70, 0.65, 0.70],
                [0.60, 0.55, 0.60],
                [0.50, 0.45, 0.50],
                [0.40, 0.35, 0.40],
                [0.30, 0.25, 0.30],
                [0.20, 0.00, 0.20],
                [0.10, 0.00, 0.10],
                [0.00, 0.00, 0.00],
            ]
        else:
            residue_replicates = [
                [0.22, 0.20, 0.22],
                [0.20, 0.18, 0.20],
                [0.18, 0.16, 0.18],
                [0.16, 0.14, 0.16],
                [0.14, 0.12, 0.14],
                [0.10, 0.00, 0.20],
                [0.00, 0.00, 0.00],
                [0.00, 0.00, 0.00],
                [0.00, 0.00, 0.00],
                [0.00, 0.00, 0.00],
            ]

        residue_stats = []
        for resid, per_rep in enumerate(residue_replicates, start=1):
            residue_stats.append(
                AggregatedResidueStats(
                    protein_resid=resid,
                    protein_resname="ALA",
                    protein_group="nonpolar",
                    contact_fraction_mean=sum(per_rep) / len(per_rep),
                    contact_fraction_per_replicate=per_rep,
                )
            )

        coverage_per_rep = []
        mean_cf_per_rep = []
        for rep_idx in range(3):
            rep_vals = [r[rep_idx] for r in residue_replicates]
            contacted = sum(1 for v in rep_vals if v > 0.0)
            coverage_per_rep.append(contacted / len(residue_replicates))
            mean_cf_per_rep.append(sum(rep_vals) / len(rep_vals))

        return AggregatedContactResult(
            residue_stats=residue_stats,
            n_replicates=3,
            total_frames_per_replicate=[1000, 1000, 1000],
            timestep_ps=1.0,
            criteria_label="any_atom_4.5A",
            criteria_cutoff=4.5,
            coverage_mean=sum(coverage_per_rep) / len(coverage_per_rep),
            coverage_sem=0.0,
            mean_contact_fraction=sum(mean_cf_per_rep) / len(mean_cf_per_rep),
            mean_contact_fraction_sem=0.0,
            residence_time_by_polymer_type={"SBM": (10.0, 1.0)},
        )

    @staticmethod
    def _make_compare_ctx(tmp_path: Path, settings: BaseModel):
        """Build a 3-condition ComparisonContext for contacts compare tests.

        Parameters
        ----------
        tmp_path : Path
            Temporary pytest directory.
        settings : BaseModel
            Contacts settings instance used by the context.

        Returns
        -------
        ComparisonContext
            Context ready for ContactsAnalysis.compare().
        """
        from polyzymd.analyses.base import ComparisonContext, Condition

        labels = ["Control", "TreatmentStrong", "TreatmentMild"]
        conditions = []
        analysis_dirs = {}

        for label in labels:
            cond = Condition(
                label=label,
                config_path=tmp_path / label / "config.yaml",
                replicates=(1, 2, 3),
                sim_config=MagicMock(spec=["temperature"]),
            )
            conditions.append(cond)

            analysis_dir = tmp_path / label / "contacts"
            (analysis_dir / "aggregated").mkdir(parents=True, exist_ok=True)
            analysis_dirs[label] = analysis_dir

        return ComparisonContext(
            name="test_contacts_compare",
            conditions=conditions,
            excluded_conditions=[],
            control_label="Control",
            analysis_dirs=analysis_dirs,
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=settings,
            recompute=False,
        )

    def test_compare_with_nondefault_fdr_alpha(self, tmp_path: Path) -> None:
        """Non-default fdr_alpha should be used by compare() output and BH correction."""
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings(fdr_alpha=0.10, min_effect_size=0.0, top_residues=0)
        ctx = self._make_compare_ctx(tmp_path, settings)

        labels = [c.label for c in ctx.conditions]
        mock_results = {label: self._make_real_agg_result(i) for i, label in enumerate(labels)}

        def load_side_effect(agg_dir: Path) -> AggregatedContactResult | None:
            label = agg_dir.parent.parent.name
            return mock_results.get(label)

        with patch.object(analysis, "_load_aggregated_result", side_effect=load_side_effect):
            result = analysis.compare(ctx)

        assert result is not None
        assert result.fdr_alpha == pytest.approx(0.10)

        all_pairwise = [
            agg for pair in result.pairwise_comparisons for agg in pair.aggregate_comparisons
        ]
        assert all(entry.p_value_adjusted is not None for entry in all_pairwise)
        assert all(entry.p_value_adjusted >= 0.0 for entry in all_pairwise)

        assert len(result.anova) == 2
        assert all(a.p_value_adjusted is not None for a in result.anova)

    def test_compare_with_effect_size_threshold(self, tmp_path: Path) -> None:
        """Effect-size threshold should tag both passing and failing comparisons."""
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings(fdr_alpha=0.05, min_effect_size=0.5, top_residues=0)
        ctx = self._make_compare_ctx(tmp_path, settings)

        labels = [c.label for c in ctx.conditions]
        mock_results = {label: self._make_real_agg_result(i) for i, label in enumerate(labels)}

        def load_side_effect(agg_dir: Path) -> AggregatedContactResult | None:
            label = agg_dir.parent.parent.name
            return mock_results.get(label)

        with patch.object(analysis, "_load_aggregated_result", side_effect=load_side_effect):
            result = analysis.compare(ctx)

        assert result is not None
        assert result.min_effect_size == pytest.approx(0.5)

        all_pairwise = [
            agg for pair in result.pairwise_comparisons for agg in pair.aggregate_comparisons
        ]
        flags = [entry.meets_effect_size_threshold for entry in all_pairwise]
        assert any(flags)
        assert not all(flags)

    def test_compare_with_top_residues(self, tmp_path: Path) -> None:
        """top_residues should control emitted top-contacted residue lists."""
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings(fdr_alpha=0.05, min_effect_size=0.0, top_residues=3)
        ctx = self._make_compare_ctx(tmp_path, settings)

        labels = [c.label for c in ctx.conditions]
        mock_results = {label: self._make_real_agg_result(i) for i, label in enumerate(labels)}

        def load_side_effect(agg_dir: Path) -> AggregatedContactResult | None:
            label = agg_dir.parent.parent.name
            return mock_results.get(label)

        with patch.object(analysis, "_load_aggregated_result", side_effect=load_side_effect):
            result = analysis.compare(ctx)

        assert result is not None
        assert result.top_residues == 3
        assert result.top_contacted_residues is not None

        top_dict = result.top_contacted_residues
        assert set(top_dict.keys()) == {"Control", "TreatmentStrong", "TreatmentMild"}

        for entries in top_dict.values():
            assert len(entries) <= 3
            if len(entries) > 1:
                fracs = [value[2] for value in entries]
                assert fracs == sorted(fracs, reverse=True)

    def test_compare_fdr_changes_significance(self, tmp_path: Path) -> None:
        """Permissive alpha should retain at least as many significant pairwise tests."""
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()

        strict_settings = ContactsSettings(fdr_alpha=0.001, min_effect_size=0.0, top_residues=0)
        strict_ctx = self._make_compare_ctx(tmp_path / "strict", strict_settings)

        permissive_settings = ContactsSettings(fdr_alpha=0.50, min_effect_size=0.0, top_residues=0)
        permissive_ctx = self._make_compare_ctx(tmp_path / "permissive", permissive_settings)

        strict_labels = [c.label for c in strict_ctx.conditions]
        strict_results = {
            label: self._make_real_agg_result(i) for i, label in enumerate(strict_labels)
        }
        permissive_labels = [c.label for c in permissive_ctx.conditions]
        permissive_results = {
            label: self._make_real_agg_result(i) for i, label in enumerate(permissive_labels)
        }

        def strict_side_effect(agg_dir: Path) -> AggregatedContactResult | None:
            return strict_results.get(agg_dir.parent.parent.name)

        def permissive_side_effect(agg_dir: Path) -> AggregatedContactResult | None:
            return permissive_results.get(agg_dir.parent.parent.name)

        with patch.object(analysis, "_load_aggregated_result", side_effect=strict_side_effect):
            strict_result = analysis.compare(strict_ctx)

        with patch.object(analysis, "_load_aggregated_result", side_effect=permissive_side_effect):
            permissive_result = analysis.compare(permissive_ctx)

        assert strict_result is not None
        assert permissive_result is not None

        strict_sig = sum(
            1
            for pair in strict_result.pairwise_comparisons
            for entry in pair.aggregate_comparisons
            if entry.significant
        )
        permissive_sig = sum(
            1
            for pair in permissive_result.pairwise_comparisons
            for entry in pair.aggregate_comparisons
            if entry.significant
        )

        assert strict_sig <= permissive_sig


class TestYAMLToSettingsWiring:
    """Test that YAML and parser inputs reach plugin settings correctly."""

    def test_contacts_settings_from_yaml_dict(self) -> None:
        """Non-default FDR/effect-size/top-residue values parse into ContactsSettings."""
        from polyzymd.analyses.contacts import ContactsSettings

        settings = ContactsSettings(
            fdr_alpha=0.10,
            min_effect_size=0.3,
            top_residues=5,
        )
        assert settings.fdr_alpha == pytest.approx(0.10)
        assert settings.min_effect_size == pytest.approx(0.3)
        assert settings.top_residues == 5

    def test_contacts_settings_from_dict_construction(self) -> None:
        """ContactsSettings should accept plain dict values as from YAML parsing."""
        from polyzymd.analyses.contacts import ContactsSettings

        yaml_dict = {
            "fdr_alpha": 0.10,
            "min_effect_size": 0.3,
            "top_residues": 5,
        }
        settings = ContactsSettings(**yaml_dict)

        assert settings.fdr_alpha == pytest.approx(0.10)
        assert settings.min_effect_size == pytest.approx(0.3)
        assert settings.top_residues == 5

    def test_contacts_settings_from_comparison_yaml(self, tmp_path: Path) -> None:
        """ComparisonConfig.from_yaml should preserve non-default contacts settings."""
        from polyzymd.config.comparison import ComparisonConfig

        cfg_path = tmp_path / "comparison.yaml"
        cfg_path.write_text(
            "\n".join(
                [
                    'name: "wiring-test"',
                    "control: Control",
                    "conditions:",
                    '  - label: "Control"',
                    '    config: "control/config.yaml"',
                    "    replicates: [1, 2, 3]",
                    '  - label: "Treatment"',
                    '    config: "treatment/config.yaml"',
                    "    replicates: [1, 2, 3]",
                    "plugins:",
                    "  contacts:",
                    "    fdr_alpha: 0.10",
                    "    min_effect_size: 0.3",
                    "    top_residues: 5",
                ]
            )
        )

        cfg = ComparisonConfig.from_yaml(cfg_path)
        contacts_settings = cfg.plugins.get("contacts")

        assert contacts_settings is not None
        assert contacts_settings.fdr_alpha == pytest.approx(0.10)
        assert contacts_settings.min_effect_size == pytest.approx(0.3)
        assert contacts_settings.top_residues == 5


def test_contacts_anova_fdr_correction() -> None:
    """Contacts ANOVA entries should receive BH-adjusted p-values."""
    cls = get_analysis("contacts")
    analysis = cls()

    anova = [
        ContactsANOVASummary(metric="coverage", f_statistic=6.0, p_value=0.01, significant=True),
        ContactsANOVASummary(
            metric="mean_contact_fraction",
            f_statistic=0.5,
            p_value=0.20,
            significant=False,
        ),
    ]

    analysis._apply_fdr_correction([], anova, fdr_alpha=0.05)

    assert anova[0].p_value_adjusted == pytest.approx(0.02)
    assert anova[1].p_value_adjusted == pytest.approx(0.20)
    assert anova[0].significant is True
    assert anova[1].significant is False


def test_contacts_result_backward_compat() -> None:
    """Contacts result should deserialize old JSON missing new fields."""
    payload = {
        "name": "legacy_contacts",
        "contacts_name": "contacts",
        "contacts_description": None,
        "polymer_selection": "chainID C",
        "protein_selection": "protein",
        "cutoff": 4.5,
        "contact_criteria": "distance_cutoff",
        "control_label": "Control",
        "conditions": [
            {
                "label": "Control",
                "config_path": "/tmp/control.yaml",
                "n_replicates": 3,
                "n_residues": 100,
                "coverage_mean": 0.2,
                "coverage_sem": 0.01,
                "mean_contact_fraction": 0.05,
                "mean_contact_fraction_sem": 0.005,
            }
        ],
        "pairwise_comparisons": [
            {
                "condition_a": "Control",
                "condition_b": "Treatment",
                "aggregate_comparisons": [
                    {
                        "metric": "coverage",
                        "condition_a": "Control",
                        "condition_b": "Treatment",
                        "condition_a_mean": 0.2,
                        "condition_a_sem": 0.01,
                        "condition_b_mean": 0.25,
                        "condition_b_sem": 0.01,
                        "t_statistic": 2.2,
                        "p_value": 0.03,
                        "cohens_d": 0.8,
                        "effect_size_interpretation": "large",
                        "significant": True,
                        "percent_change": 25.0,
                        "direction": "increased",
                    }
                ],
            }
        ],
        "anova": [],
        "ranking_by_coverage": ["Control"],
        "ranking_by_contact_fraction": ["Control"],
        "equilibration_time": "10ns",
        "created_at": "2026-01-01T00:00:00",
        "polyzymd_version": "1.3.0",
    }

    result = ContactsComparisonResult.model_validate_json(json.dumps(payload))

    assert result.fdr_alpha == pytest.approx(0.05)
    assert result.min_effect_size == pytest.approx(0.0)
    assert result.top_residues == 0


def test_contacts_formatter_shows_adjusted_pvalues() -> None:
    """Contacts formatter should include both raw and adjusted p-values."""
    result = _make_contacts_result_for_formatting()
    formatted = format_contacts_console_table(result)

    assert "0.0300" in formatted
    assert "0.0400" in formatted
