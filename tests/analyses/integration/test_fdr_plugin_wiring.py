"""Focused tests for FDR/effect-size/top-residue wiring in analysis plugins."""

from __future__ import annotations

import math
from collections.abc import Sequence
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
import pytest
from pydantic import BaseModel

from polyzymd.analyses import get_analysis
from polyzymd.analyses.base import ComparisonContext, Condition
from polyzymd.analyses.contacts._comparison import (
    apply_effect_size_threshold,
    apply_fdr_correction,
    compute_top_contacted_residues,
)
from polyzymd.analyses.contacts._comparison_results import (
    AggregateComparisonResult,
    ContactsANOVASummary,
    ContactsComparisonResult,
    ContactsConditionSummary,
    ContactsPairwiseComparison,
)
from polyzymd.analyses.contacts._formatters import format_contacts_console_table
from polyzymd.analyses.contacts._identity import contacts_detection_fingerprint
from polyzymd.analyses.mda import ConditionArtifact


def _metric_summary(name: str, values: Sequence[float]) -> dict[str, object]:
    """Build a canonical aggregate metric summary from replicate values."""

    metric_values = [float(value) for value in values]
    std = float(np.std(metric_values, ddof=1)) if len(metric_values) > 1 else 0.0
    sem = float(std / math.sqrt(len(metric_values))) if len(metric_values) > 1 else 0.0
    mean = float(np.mean(metric_values)) if metric_values else 0.0
    return {
        "name": name,
        "values": metric_values,
        "mean": mean,
        "sem": sem,
        "std": std,
        "n": len(values),
    }


def _make_contacts_artifact(
    label: str,
    residue_rows: Sequence[dict[str, object]],
    settings: BaseModel,
    *,
    replicates: Sequence[int] = (1,),
) -> ConditionArtifact:
    """Build a canonical contacts condition artifact for wiring tests.

    Parameters
    ----------
    label : str
        Condition label stored in the artifact.
    residue_rows : sequence of dict
        Residue summary rows containing ``contact_fraction_per_replicate`` values.
    settings : BaseModel
        Contacts settings used to produce the detection fingerprint.
    replicates : sequence of int, optional
        Replicate IDs represented by each row's value vector, by default ``(1,)``.

    Returns
    -------
    ConditionArtifact
        Canonical contacts aggregate artifact accepted by ``ContactsAnalysis.compare()``.
    """

    replicate_ids = [int(replicate) for replicate in replicates]
    normalized_rows = []
    for row in residue_rows:
        values = [float(value) for value in row["contact_fraction_per_replicate"]]
        normalized_rows.append(
            {
                "protein_resid": int(row["protein_resid"]),
                "protein_resname": str(row["protein_resname"]),
                "protein_chain_id": str(row.get("protein_chain_id", "A")),
                "protein_group": str(row.get("protein_group", "nonpolar")),
                "contact_fraction_mean": float(np.mean(values)) if values else 0.0,
                "contact_fraction_per_replicate": values,
            }
        )

    coverage_values = []
    contact_values = []
    for rep_idx, _replicate in enumerate(replicate_ids):
        rep_values = [row["contact_fraction_per_replicate"][rep_idx] for row in normalized_rows]
        coverage_values.append(
            sum(1 for value in rep_values if float(value) > 0.0) / len(normalized_rows)
            if normalized_rows
            else 0.0
        )
        contact_values.append(float(np.mean(rep_values)) if rep_values else 0.0)

    replicate_metrics = {
        str(replicate): {
            "coverage": coverage_values[idx],
            "mean_contact_fraction": contact_values[idx],
        }
        for idx, replicate in enumerate(replicate_ids)
    }
    return ConditionArtifact(
        analysis_name="contacts",
        condition_label=label,
        replicates=replicate_ids,
        payload={
            "metrics": {
                "coverage": _metric_summary("coverage", coverage_values),
                "mean_contact_fraction": _metric_summary("mean_contact_fraction", contact_values),
            },
            "replicate_metrics": replicate_metrics,
            "n_replicates": len(replicate_ids),
            "n_residues": len(normalized_rows),
            "total_frames_per_replicate": [1000 for _ in replicate_ids],
            "criteria_cutoff": float(getattr(settings, "cutoff", 4.5)),
            "residue_stats": normalized_rows,
            "residence_time_by_polymer_type": {"SBM": {"mean_ns": 10.0, "sem_ns": 1.0}},
        },
        provenance={
            "source": "contacts_fdr_wiring_test",
            "frame_selection": {"equilibration": "10ns"},
        },
        metadata={
            "contacts_detection_fingerprint": contacts_detection_fingerprint(settings),
            "compute_residence_times": bool(getattr(settings, "compute_residence_times", True)),
            "equilibration": "10ns",
        },
    )


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
        polymer_selection="chainid C",
        protein_selection="chainid A",
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

    apply_fdr_correction(comparisons, [], fdr_alpha=0.05)

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

    apply_effect_size_threshold(comparisons, min_effect_size=0.5)

    first, second = comparisons[0].aggregate_comparisons
    assert first.meets_effect_size_threshold is True
    assert second.meets_effect_size_threshold is False


def test_contacts_top_residues_selection() -> None:
    """Contacts top-residue helper should work with canonical artifacts."""
    cls = get_analysis("contacts")
    settings = cls.Settings()

    agg_a = _make_contacts_artifact(
        "A",
        [
            {
                "protein_resid": 10,
                "protein_resname": "ALA",
                "protein_group": "nonpolar",
                "contact_fraction_per_replicate": [0.40],
            },
            {
                "protein_resid": 11,
                "protein_resname": "GLY",
                "protein_group": "nonpolar",
                "contact_fraction_per_replicate": [0.80],
            },
            {
                "protein_resid": 12,
                "protein_resname": "SER",
                "protein_group": "polar",
                "contact_fraction_per_replicate": [0.20],
            },
            {
                "protein_resid": 13,
                "protein_resname": "TYR",
                "protein_group": "aromatic",
                "contact_fraction_per_replicate": [0.60],
            },
        ],
        settings,
    )
    agg_b = _make_contacts_artifact(
        "B",
        [
            {
                "protein_resid": 20,
                "protein_resname": "LEU",
                "protein_group": "nonpolar",
                "contact_fraction_per_replicate": [0.10],
            },
            {
                "protein_resid": 21,
                "protein_resname": "VAL",
                "protein_group": "nonpolar",
                "contact_fraction_per_replicate": [0.30],
            },
            {
                "protein_resid": 22,
                "protein_resname": "ASP",
                "protein_group": "charged_negative",
                "contact_fraction_per_replicate": [0.50],
            },
            {
                "protein_resid": 23,
                "protein_resname": "LYS",
                "protein_group": "charged_positive",
                "contact_fraction_per_replicate": [0.90],
            },
        ],
        settings,
    )

    condition_data = [
        (SimpleNamespace(label="A"), {"agg_result": agg_a}),
        (SimpleNamespace(label="B"), {"agg_result": agg_b}),
    ]

    top = compute_top_contacted_residues(condition_data, top_n=3)

    assert top is not None
    assert top["A"] == [(11, "GLY", 0.8), (13, "TYR", 0.6), (10, "ALA", 0.4)]
    assert top["B"] == [(23, "LYS", 0.9), (22, "ASP", 0.5), (21, "VAL", 0.3)]


def test_contacts_settings_wiring_to_helper_arguments() -> None:
    """Contacts settings values should flow to helper arguments via ctx.settings."""
    ctx = SimpleNamespace(
        settings=get_analysis("contacts").Settings(
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

    apply_fdr_correction(comparisons, [], ctx.settings.fdr_alpha)
    apply_effect_size_threshold(comparisons, ctx.settings.min_effect_size)

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
                "agg_result": _make_contacts_artifact(
                    "Treatment",
                    [
                        {
                            "protein_resid": 1,
                            "protein_resname": "ALA",
                            "protein_group": "nonpolar",
                            "contact_fraction_per_replicate": [0.9],
                        },
                        {
                            "protein_resid": 2,
                            "protein_resname": "GLY",
                            "protein_group": "nonpolar",
                            "contact_fraction_per_replicate": [0.4],
                        },
                        {
                            "protein_resid": 3,
                            "protein_resname": "SER",
                            "protein_group": "polar",
                            "contact_fraction_per_replicate": [0.2],
                        },
                    ],
                    ctx.settings,
                )
            },
        )
    ]
    top = compute_top_contacted_residues(condition_data, ctx.settings.top_residues)

    assert top is not None
    assert top["Treatment"] == [(1, "ALA", 0.9), (2, "GLY", 0.4)]


class TestContactsCompareWiring:
    """Test that fdr_alpha, min_effect_size, and top_residues flow through compare()."""

    @staticmethod
    def _make_real_agg_result(condition_idx: int, settings: BaseModel) -> ConditionArtifact:
        """Construct a canonical contacts artifact with deterministic replicate data.

        Parameters
        ----------
        condition_idx : int
            Condition index: 0=Control, 1=TreatmentStrong, 2=TreatmentMild.
        settings : BaseModel
            Contacts settings used to fingerprint the artifact.

        Returns
        -------
        ConditionArtifact
            Canonical aggregate artifact suitable for ContactsAnalysis.compare().
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

        residue_rows = []
        for resid, per_rep in enumerate(residue_replicates, start=1):
            residue_rows.append(
                {
                    "protein_resid": resid,
                    "protein_resname": "ALA",
                    "protein_group": "nonpolar",
                    "protein_chain_id": "A",
                    "contact_fraction_per_replicate": per_rep,
                }
            )

        return _make_contacts_artifact(
            ["Control", "TreatmentStrong", "TreatmentMild"][condition_idx],
            residue_rows,
            settings,
            replicates=(1, 2, 3),
        )

    @staticmethod
    def _make_compare_ctx(tmp_path: Path, settings: BaseModel) -> ComparisonContext:
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
            aggregated_results={
                label: TestContactsCompareWiring._make_real_agg_result(idx, settings)
                for idx, label in enumerate(labels)
            },
        )

    def test_compare_with_nondefault_fdr_alpha(self, tmp_path: Path) -> None:
        """Non-default fdr_alpha should be used by compare() output and BH correction."""
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings(fdr_alpha=0.10, min_effect_size=0.0, top_residues=0)
        ctx = self._make_compare_ctx(tmp_path, settings)

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

        strict_result = analysis.compare(strict_ctx)
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
    anova = [
        ContactsANOVASummary(metric="coverage", f_statistic=6.0, p_value=0.01, significant=True),
        ContactsANOVASummary(
            metric="mean_contact_fraction",
            f_statistic=0.5,
            p_value=0.20,
            significant=False,
        ),
    ]

    apply_fdr_correction([], anova, fdr_alpha=0.05)

    assert anova[0].p_value_adjusted == pytest.approx(0.02)
    assert anova[1].p_value_adjusted == pytest.approx(0.20)
    assert anova[0].significant is True
    assert anova[1].significant is False


def test_contacts_formatter_shows_adjusted_pvalues() -> None:
    """Contacts formatter should include both raw and adjusted p-values."""
    result = _make_contacts_result_for_formatting()
    formatted = format_contacts_console_table(result)

    assert "0.0300" in formatted
    assert "0.0400" in formatted
