"""Tests for the contacts analysis plugin.

Covers discovery, settings, compute_replicate, aggregate, compare (full override),
filter_conditions, binding preference sub-pipeline, plot delegation,
AggregatedResultClass, and per-replicate metric helpers.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from pydantic import BaseModel


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


class TestDiscovery:
    """The plugin is found by the automatic discovery system."""

    def test_discovery_by_name(self):
        from polyzymd.analyses.discovery import get_analysis

        cls = get_analysis("contacts")
        assert cls is not None
        assert cls.name == "contacts"

    def test_listed(self):
        from polyzymd.analyses.discovery import list_analyses

        analyses = list_analyses()
        names = list(analyses.keys())
        assert "contacts" in names

    def test_all_names(self):
        from polyzymd.analyses.discovery import list_all_names

        names = list_all_names()
        assert "contacts" in names


# ---------------------------------------------------------------------------
# Class Attributes
# ---------------------------------------------------------------------------


class TestClassAttributes:
    """Verify ClassVar declarations on ContactsAnalysis."""

    def test_name(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        assert ContactsAnalysis.name == "contacts"

    def test_settings_type(self):
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        assert ContactsAnalysis.Settings is ContactsSettings

    def test_aliases(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        assert ContactsAnalysis.aliases == ()

    def test_dependencies(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        assert ContactsAnalysis.dependencies == ()

    def test_min_replicates(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        assert ContactsAnalysis.min_replicates == 2


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class TestSettings:
    """Validate ContactsSettings."""

    def test_defaults(self):
        from polyzymd.analyses.contacts import ContactsSettings

        s = ContactsSettings()
        assert s.polymer_selection == "chainID C"
        assert s.protein_selection == "protein"
        assert s.cutoff == 4.5
        assert s.grouping == "aa_class"
        assert s.compute_residence_times is True
        assert s.compute_binding_preference is False
        assert s.fdr_alpha == 0.05

    def test_custom_values(self):
        from polyzymd.analyses.contacts import ContactsSettings

        s = ContactsSettings(
            polymer_selection="chainID D",
            protein_selection="protein and name CA",
            cutoff=5.0,
            grouping="none",
            compute_binding_preference=True,
            surface_exposure_threshold=0.3,
            fdr_alpha=0.01,
        )
        assert s.polymer_selection == "chainID D"
        assert s.cutoff == 5.0
        assert s.grouping == "none"
        assert s.compute_binding_preference is True
        assert s.surface_exposure_threshold == 0.3
        assert s.fdr_alpha == 0.01

    def test_invalid_grouping(self):
        from polyzymd.analyses.contacts import ContactsSettings

        with pytest.raises(Exception):
            ContactsSettings(grouping="invalid")

    def test_invalid_fdr_alpha(self):
        from polyzymd.analyses.contacts import ContactsSettings

        with pytest.raises(Exception):
            ContactsSettings(fdr_alpha=1.5)

    def test_protein_partitions_validation(self):
        from polyzymd.analyses.contacts import ContactsSettings

        # Requires protein_groups to be defined
        with pytest.raises(Exception):
            ContactsSettings(
                protein_partitions={"part1": ["group1"]},
                # protein_groups not defined
            )

    def test_protein_partitions_undefined_group(self):
        from polyzymd.analyses.contacts import ContactsSettings

        with pytest.raises(Exception):
            ContactsSettings(
                protein_groups={"existing": [1, 2]},
                protein_partitions={"part1": ["nonexistent"]},
            )

    def test_protein_partitions_overlapping_groups(self):
        from polyzymd.analyses.contacts import ContactsSettings

        with pytest.raises(Exception):
            ContactsSettings(
                protein_groups={"g1": [1, 2, 3], "g2": [3, 4, 5]},
                protein_partitions={"part1": ["g1", "g2"]},  # resid 3 overlaps
            )

    def test_valid_protein_partitions(self):
        from polyzymd.analyses.contacts import ContactsSettings

        s = ContactsSettings(
            protein_groups={"g1": [1, 2, 3], "g2": [4, 5, 6]},
            protein_partitions={"part1": ["g1", "g2"]},
        )
        assert s.protein_partitions is not None

    def test_serialization_roundtrip(self):
        from polyzymd.analyses.contacts import ContactsSettings

        s = ContactsSettings(
            cutoff=5.0,
            polymer_types=["SBM", "EGM"],
            compute_binding_preference=True,
        )
        d = s.model_dump()
        s2 = ContactsSettings.model_validate(d)
        assert s2.cutoff == 5.0
        assert s2.polymer_types == ["SBM", "EGM"]
        assert s2.compute_binding_preference is True


# ---------------------------------------------------------------------------
# compute_replicate
# ---------------------------------------------------------------------------


def _make_mock_contact_result(replicate: int = 1):
    """Create a mock ContactResult."""
    mock = MagicMock()
    mock.replicate = replicate
    mock.save = MagicMock()
    mock.has_per_residue_statistics = MagicMock(return_value=True)
    return mock


class TestComputeReplicate:
    """compute_replicate delegates to ParallelContactAnalyzer."""

    def test_delegates_to_analyzer(self, tmp_path):
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings()

        mock_sim_config = MagicMock()
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2, 3),
            sim_config=mock_sim_config,
        )
        output_dir = tmp_path / "run_1"
        ctx = ReplicateContext(
            condition=cond,
            replicate=1,
            sim_config=mock_sim_config,
            output_dir=output_dir,
            equilibration="10ns",
            recompute=False,
            settings=settings,
        )

        mock_result = _make_mock_contact_result(replicate=1)

        with (
            patch("polyzymd.analyses.contacts.ParallelContactAnalyzer") as MockAnalyzer,
            patch("polyzymd.analyses.shared.loader.TrajectoryLoader") as MockLoader,
            patch("MDAnalysis.Universe") as MockUniverse,
        ):
            # Set up loader
            mock_traj_info = MagicMock()
            mock_traj_info.trajectory_files = [Path("/tmp/traj.xtc")]
            mock_traj_info.topology_file = Path("/tmp/top.pdb")
            MockLoader.return_value.get_trajectory_info.return_value = mock_traj_info

            # Set up universe
            mock_universe = MagicMock()
            mock_universe.trajectory.dt = 10.0
            MockUniverse.return_value = mock_universe

            # Set up analyzer
            MockAnalyzer.return_value.run.return_value = mock_result

            result = analysis.compute_replicate(ctx, 1)

        assert result is mock_result
        MockAnalyzer.assert_called_once()
        call_kwargs = MockAnalyzer.call_args[1]
        assert call_kwargs["cutoff"] == 4.5

    def test_returns_none_on_missing_trajectory(self, tmp_path):
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings()

        mock_sim_config = MagicMock()
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim_config,
        )
        ctx = ReplicateContext(
            condition=cond,
            replicate=1,
            sim_config=mock_sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=False,
            settings=settings,
        )

        with patch("polyzymd.analyses.shared.loader.TrajectoryLoader") as MockLoader:
            MockLoader.return_value.get_trajectory_info.side_effect = FileNotFoundError(
                "No trajectory"
            )
            result = analysis.compute_replicate(ctx, 1)

        assert result is None


# ---------------------------------------------------------------------------
# aggregate
# ---------------------------------------------------------------------------


def _make_mock_agg_result(n_replicates: int = 3, n_residues: int = 5):
    """Create a mock AggregatedContactResult."""
    mock = MagicMock()
    mock.n_replicates = n_replicates
    mock.n_residues = n_residues
    mock.coverage_mean = 0.8
    mock.coverage_sem = 0.05
    mock.mean_contact_fraction = 0.3
    mock.mean_contact_fraction_sem = 0.02
    mock.residence_time_by_polymer_type = {"SBM": (10.0, 1.0)}
    mock.save = MagicMock()

    residue_stats = []
    for i in range(n_residues):
        rs = MagicMock()
        rs.protein_resid = i + 1
        # Per-replicate contact fractions
        fracs = [0.2 + i * 0.05 + rep * 0.01 for rep in range(n_replicates)]
        rs.contact_fraction_per_replicate = fracs
        residue_stats.append(rs)
    mock.residue_stats = residue_stats

    return mock


class TestAggregate:
    """aggregate delegates to aggregate_contact_results."""

    def test_aggregate_delegates(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        output_dir = tmp_path / "aggregated"

        mock_sim_config = MagicMock()
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2, 3),
            sim_config=mock_sim_config,
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1, 2, 3),
            output_dir=output_dir,
            equilibration="10ns",
            settings=settings,
        )

        mock_results = [_make_mock_contact_result(i) for i in range(1, 4)]
        mock_agg = _make_mock_agg_result(n_replicates=3)

        with patch(
            "polyzymd.analyses.contacts._aggregator.aggregate_contact_results",
            return_value=mock_agg,
        ) as mock_fn:
            result = analysis.aggregate(ctx, mock_results)

        mock_fn.assert_called_once()
        assert result is mock_agg
        mock_agg.save.assert_called_once()

    def test_aggregate_saves_file(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings()
        output_dir = tmp_path / "aggregated"

        mock_sim_config = MagicMock()
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 3),
            sim_config=mock_sim_config,
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1, 3),
            output_dir=output_dir,
            equilibration="10ns",
            settings=settings,
        )

        mock_results = [_make_mock_contact_result(1), _make_mock_contact_result(3)]
        mock_agg = _make_mock_agg_result(n_replicates=2)

        with patch(
            "polyzymd.analyses.contacts._aggregator.aggregate_contact_results",
            return_value=mock_agg,
        ):
            analysis.aggregate(ctx, mock_results)

        # Check save was called with correct filename
        save_call = mock_agg.save.call_args[0][0]
        assert "reps1-3" in str(save_call)


# ---------------------------------------------------------------------------
# Per-replicate metric computation
# ---------------------------------------------------------------------------


class TestPerReplicateMetrics:
    """Test _compute_coverage_per_replicate and _compute_contact_fraction_per_replicate."""

    def test_coverage_per_replicate(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        agg = _make_mock_agg_result(n_replicates=3, n_residues=5)
        # By default all residues have >0 fractions, so coverage should be 1.0
        coverages = ContactsAnalysis._compute_coverage_per_replicate(agg)
        assert len(coverages) == 3
        for c in coverages:
            assert c == 1.0  # all residues have >0 contact fraction

    def test_coverage_with_zero_fractions(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        agg = MagicMock()
        agg.n_replicates = 2
        agg.n_residues = 3

        # Residue 1: contacted in both reps
        # Residue 2: contacted only in rep 0
        # Residue 3: never contacted
        rs1 = MagicMock()
        rs1.contact_fraction_per_replicate = [0.5, 0.3]
        rs2 = MagicMock()
        rs2.contact_fraction_per_replicate = [0.1, 0.0]
        rs3 = MagicMock()
        rs3.contact_fraction_per_replicate = [0.0, 0.0]
        agg.residue_stats = [rs1, rs2, rs3]

        coverages = ContactsAnalysis._compute_coverage_per_replicate(agg)
        assert len(coverages) == 2
        assert coverages[0] == pytest.approx(2 / 3)  # rs1 + rs2 contacted in rep 0
        assert coverages[1] == pytest.approx(1 / 3)  # only rs1 contacted in rep 1

    def test_contact_fraction_per_replicate(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        agg = MagicMock()
        agg.n_replicates = 2

        rs1 = MagicMock()
        rs1.contact_fraction_per_replicate = [0.4, 0.6]
        rs2 = MagicMock()
        rs2.contact_fraction_per_replicate = [0.2, 0.8]
        agg.residue_stats = [rs1, rs2]

        fractions = ContactsAnalysis._compute_contact_fraction_per_replicate(agg)
        assert len(fractions) == 2
        assert fractions[0] == pytest.approx(0.3)  # (0.4 + 0.2) / 2
        assert fractions[1] == pytest.approx(0.7)  # (0.6 + 0.8) / 2


# ---------------------------------------------------------------------------
# Residue set validation
# ---------------------------------------------------------------------------


class TestResidueSetValidation:
    """Test _validate_residue_sets."""

    def test_matching_residue_sets(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        cond_a = MagicMock()
        cond_a.label = "A"
        data_a = {"agg_result": MagicMock()}
        rs_a = [MagicMock(protein_resid=i) for i in [1, 2, 3]]
        data_a["agg_result"].residue_stats = rs_a

        cond_b = MagicMock()
        cond_b.label = "B"
        data_b = {"agg_result": MagicMock()}
        rs_b = [MagicMock(protein_resid=i) for i in [1, 2, 3]]
        data_b["agg_result"].residue_stats = rs_b

        # Should not raise
        ContactsAnalysis._validate_residue_sets([(cond_a, data_a), (cond_b, data_b)])

    def test_mismatched_residue_sets(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        cond_a = MagicMock()
        cond_a.label = "A"
        data_a = {"agg_result": MagicMock()}
        data_a["agg_result"].residue_stats = [MagicMock(protein_resid=i) for i in [1, 2, 3]]

        cond_b = MagicMock()
        cond_b.label = "B"
        data_b = {"agg_result": MagicMock()}
        data_b["agg_result"].residue_stats = [MagicMock(protein_resid=i) for i in [1, 2, 4]]

        with pytest.raises(ValueError, match="Residue set mismatch"):
            ContactsAnalysis._validate_residue_sets([(cond_a, data_a), (cond_b, data_b)])


# ---------------------------------------------------------------------------
# filter_conditions
# ---------------------------------------------------------------------------


class TestFilterConditions:
    """Test filter_conditions polymer detection logic."""

    def test_condition_with_cached_results_included(self, tmp_path):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()

        # Create a cached result file
        projects_dir = tmp_path / "projects"
        contacts_dir = projects_dir / "analysis" / "contacts"
        contacts_dir.mkdir(parents=True)
        (contacts_dir / "contacts_rep1.json").write_text("{}")

        mock_sim_config = MagicMock()
        mock_sim_config.output.projects_directory = projects_dir

        cond = Condition(
            label="With Polymer",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim_config,
        )

        result = analysis.filter_conditions([cond])
        assert len(result) == 1
        assert result[0].label == "With Polymer"

    def test_condition_without_polymer_excluded(self, tmp_path):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()

        # No cached results, topology exists but no polymer atoms
        projects_dir = tmp_path / "projects"
        projects_dir.mkdir()
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()

        mock_sim_config = MagicMock()
        mock_sim_config.output.projects_directory = projects_dir
        mock_sim_config.get_working_directory.return_value = run_dir

        # Create a dummy topology file
        (run_dir / "solvated_system.pdb").write_text("ATOM ...")

        cond = Condition(
            label="No Polymer",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim_config,
        )

        # Mock MDAnalysis to return empty selection
        with patch("MDAnalysis.Universe") as MockUniverse:
            mock_universe = MagicMock()
            mock_universe.select_atoms.return_value = MagicMock(__len__=lambda s: 0)
            MockUniverse.return_value = mock_universe

            result = analysis.filter_conditions([cond])

        assert len(result) == 0

    def test_error_during_check_includes_condition(self, tmp_path):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()

        mock_sim_config = MagicMock()
        mock_sim_config.output.projects_directory = tmp_path / "projects"
        mock_sim_config.get_working_directory.side_effect = Exception("boom")

        cond = Condition(
            label="ErrorCond",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=mock_sim_config,
        )

        result = analysis.filter_conditions([cond])
        assert len(result) == 1  # Included because of error (fail-open)


# ---------------------------------------------------------------------------
# compare — full override
# ---------------------------------------------------------------------------


def _make_condition_and_data(label: str, n_reps: int = 3, n_residues: int = 5):
    """Build a mock condition + data dict for compare() tests."""
    from polyzymd.analyses.base import Condition

    mock_sim_config = MagicMock()
    cond = Condition(
        label=label,
        config_path=Path(f"/tmp/{label}/config.yaml"),
        replicates=tuple(range(1, n_reps + 1)),
        sim_config=mock_sim_config,
    )

    agg = _make_mock_agg_result(n_replicates=n_reps, n_residues=n_residues)
    data = {
        "agg_result": agg,
        "coverage_per_replicate": [0.8 + i * 0.01 for i in range(n_reps)],
        "contact_fraction_per_replicate": [0.3 + i * 0.005 for i in range(n_reps)],
    }
    return cond, data


class TestCompare:
    """Test compare() full override."""

    def _make_ctx(
        self,
        tmp_path,
        n_conditions: int = 3,
        control: str | None = "Control",
    ):
        """Create a ComparisonContext with mock data."""
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsSettings

        labels = ["Control", "Treatment A", "Treatment B"][:n_conditions]
        conditions = []
        for label in labels:
            mock_sim = MagicMock()
            conditions.append(
                Condition(
                    label=label,
                    config_path=Path(f"/tmp/{label}/config.yaml"),
                    replicates=(1, 2, 3),
                    sim_config=mock_sim,
                )
            )

        analysis_dirs = {}
        for cond in conditions:
            agg_dir = tmp_path / cond.label / "contacts" / "aggregated"
            agg_dir.mkdir(parents=True, exist_ok=True)
            analysis_dirs[cond.label] = tmp_path / cond.label / "contacts"

        return ComparisonContext(
            name="test_comparison",
            conditions=conditions,
            excluded_conditions=[],
            control_label=control,
            analysis_dirs=analysis_dirs,
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=ContactsSettings(),
            recompute=False,
        )

    def test_compare_returns_result(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        ctx = self._make_ctx(tmp_path, n_conditions=3, control="Control")

        # Mock _load_aggregated_result to return proper aggregated results
        mock_agg_results = {}
        for cond in ctx.conditions:
            mock_agg = _make_mock_agg_result(n_replicates=3, n_residues=5)
            mock_agg_results[cond.label] = mock_agg

        def side_effect(agg_dir):
            # Extract label from path
            label = agg_dir.parent.parent.name
            return mock_agg_results.get(label)

        with (
            patch.object(analysis, "_load_aggregated_result", side_effect=side_effect),
            patch.object(analysis, "_load_or_compute_binding_preference", return_value=None),
        ):
            result = analysis.compare(ctx)

        assert result is not None
        assert result.name == "test_comparison"
        assert len(result.conditions) == 3
        assert result.control_label == "Control"
        assert len(result.ranking_by_coverage) == 3
        assert len(result.ranking_by_contact_fraction) == 3

    def test_compare_pairwise_with_control(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        ctx = self._make_ctx(tmp_path, n_conditions=3, control="Control")

        mock_agg_results = {}
        for cond in ctx.conditions:
            mock_agg_results[cond.label] = _make_mock_agg_result(3, 5)

        def side_effect(agg_dir):
            label = agg_dir.parent.parent.name
            return mock_agg_results.get(label)

        with (
            patch.object(analysis, "_load_aggregated_result", side_effect=side_effect),
            patch.object(analysis, "_load_or_compute_binding_preference", return_value=None),
        ):
            result = analysis.compare(ctx)

        # With 3 conditions and a control, expect 2 pairwise comparisons
        assert len(result.pairwise_comparisons) == 2
        for comp in result.pairwise_comparisons:
            assert comp.condition_a == "Control"
            assert len(comp.aggregate_comparisons) == 2  # coverage + contact_fraction

    def test_compare_pairwise_without_control(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        ctx = self._make_ctx(tmp_path, n_conditions=3, control=None)

        mock_agg_results = {}
        for cond in ctx.conditions:
            mock_agg_results[cond.label] = _make_mock_agg_result(3, 5)

        def side_effect(agg_dir):
            label = agg_dir.parent.parent.name
            return mock_agg_results.get(label)

        with (
            patch.object(analysis, "_load_aggregated_result", side_effect=side_effect),
            patch.object(analysis, "_load_or_compute_binding_preference", return_value=None),
        ):
            result = analysis.compare(ctx)

        # Without control: all pairs = C(3,2) = 3
        assert len(result.pairwise_comparisons) == 3

    def test_compare_anova_with_three_conditions(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        ctx = self._make_ctx(tmp_path, n_conditions=3)

        mock_agg_results = {}
        for cond in ctx.conditions:
            mock_agg_results[cond.label] = _make_mock_agg_result(3, 5)

        def side_effect(agg_dir):
            label = agg_dir.parent.parent.name
            return mock_agg_results.get(label)

        with (
            patch.object(analysis, "_load_aggregated_result", side_effect=side_effect),
            patch.object(analysis, "_load_or_compute_binding_preference", return_value=None),
        ):
            result = analysis.compare(ctx)

        # ANOVA with 3+ conditions: 2 summaries (coverage + contact_fraction)
        assert len(result.anova) == 2
        metrics = {a.metric for a in result.anova}
        assert metrics == {"coverage", "mean_contact_fraction"}

    def test_compare_no_anova_with_two_conditions(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        ctx = self._make_ctx(tmp_path, n_conditions=2, control="Control")

        mock_agg_results = {}
        for cond in ctx.conditions:
            mock_agg_results[cond.label] = _make_mock_agg_result(3, 5)

        def side_effect(agg_dir):
            label = agg_dir.parent.parent.name
            return mock_agg_results.get(label)

        with (
            patch.object(analysis, "_load_aggregated_result", side_effect=side_effect),
            patch.object(analysis, "_load_or_compute_binding_preference", return_value=None),
        ):
            result = analysis.compare(ctx)

        # No ANOVA with < 3 conditions
        assert len(result.anova) == 0

    def test_compare_returns_none_with_insufficient_conditions(self, tmp_path):
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()

        mock_sim = MagicMock()
        cond = Condition(
            label="Only",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=mock_sim,
        )
        ctx = ComparisonContext(
            name="test",
            conditions=[cond],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"Only": tmp_path / "Only" / "contacts"},
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=ContactsSettings(),
            recompute=False,
        )

        result = analysis.compare(ctx)
        assert result is None

    def test_compare_excluded_conditions_recorded(self, tmp_path):
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()

        mock_sim = MagicMock()
        conditions = [
            Condition(
                label=label,
                config_path=Path(f"/tmp/{label}/config.yaml"),
                replicates=(1, 2, 3),
                sim_config=mock_sim,
            )
            for label in ["A", "B"]
        ]
        excluded = [
            Condition(
                label="No Polymer",
                config_path=Path("/tmp/np/config.yaml"),
                replicates=(1,),
                sim_config=mock_sim,
            )
        ]

        analysis_dirs = {}
        for cond in conditions:
            agg_dir = tmp_path / cond.label / "contacts" / "aggregated"
            agg_dir.mkdir(parents=True, exist_ok=True)
            analysis_dirs[cond.label] = tmp_path / cond.label / "contacts"

        ctx = ComparisonContext(
            name="test",
            conditions=conditions,
            excluded_conditions=excluded,
            control_label=None,
            analysis_dirs=analysis_dirs,
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=ContactsSettings(),
            recompute=False,
        )

        mock_agg = _make_mock_agg_result(3, 5)

        with (
            patch.object(analysis, "_load_aggregated_result", return_value=mock_agg),
            patch.object(analysis, "_load_or_compute_binding_preference", return_value=None),
        ):
            result = analysis.compare(ctx)

        assert "No Polymer" in result.excluded_conditions


# ---------------------------------------------------------------------------
# Pairwise comparison details
# ---------------------------------------------------------------------------


class TestPairwiseComparison:
    """Test _compare_contacts_pair static method."""

    def test_pair_produces_two_metrics(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        summary_a = MagicMock()
        summary_a.coverage_mean = 0.8
        summary_a.coverage_sem = 0.05
        summary_a.mean_contact_fraction = 0.3
        summary_a.mean_contact_fraction_sem = 0.02

        summary_b = MagicMock()
        summary_b.coverage_mean = 0.9
        summary_b.coverage_sem = 0.04
        summary_b.mean_contact_fraction = 0.4
        summary_b.mean_contact_fraction_sem = 0.03

        data_a = {
            "coverage_per_replicate": [0.78, 0.80, 0.82],
            "contact_fraction_per_replicate": [0.28, 0.30, 0.32],
        }
        data_b = {
            "coverage_per_replicate": [0.88, 0.90, 0.92],
            "contact_fraction_per_replicate": [0.38, 0.40, 0.42],
        }

        comp = ContactsAnalysis._compare_contacts_pair(
            "A", summary_a, data_a, "B", summary_b, data_b
        )

        assert comp.condition_a == "A"
        assert comp.condition_b == "B"
        assert len(comp.aggregate_comparisons) == 2

        metrics = {ac.metric for ac in comp.aggregate_comparisons}
        assert metrics == {"coverage", "mean_contact_fraction"}

    def test_pair_direction_labels(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        summary_a = MagicMock(
            coverage_mean=0.5,
            coverage_sem=0.05,
            mean_contact_fraction=0.2,
            mean_contact_fraction_sem=0.02,
        )
        summary_b = MagicMock(
            coverage_mean=0.9,
            coverage_sem=0.04,
            mean_contact_fraction=0.5,
            mean_contact_fraction_sem=0.03,
        )

        data_a = {
            "coverage_per_replicate": [0.48, 0.50, 0.52],
            "contact_fraction_per_replicate": [0.18, 0.20, 0.22],
        }
        data_b = {
            "coverage_per_replicate": [0.88, 0.90, 0.92],
            "contact_fraction_per_replicate": [0.48, 0.50, 0.52],
        }

        comp = ContactsAnalysis._compare_contacts_pair(
            "Control", summary_a, data_a, "Treatment", summary_b, data_b
        )

        for ac in comp.aggregate_comparisons:
            assert ac.direction in ("increased", "decreased", "unchanged")


# ---------------------------------------------------------------------------
# ANOVA
# ---------------------------------------------------------------------------


class TestANOVA:
    """Test _compute_contacts_anova."""

    def test_anova_returns_two_summaries(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        condition_data = []
        for i in range(3):
            cond = MagicMock()
            cond.label = f"Cond{i}"
            data = {
                "coverage_per_replicate": [0.7 + i * 0.05 + j * 0.01 for j in range(3)],
                "contact_fraction_per_replicate": [0.2 + i * 0.03 + j * 0.005 for j in range(3)],
            }
            condition_data.append((cond, data))

        results = ContactsAnalysis._compute_contacts_anova(condition_data)
        assert len(results) == 2
        metrics = {r.metric for r in results}
        assert metrics == {"coverage", "mean_contact_fraction"}

        for r in results:
            assert hasattr(r, "f_statistic")
            assert hasattr(r, "p_value")
            assert hasattr(r, "significant")


# ---------------------------------------------------------------------------
# AggregatedResultClass and _deserialize_result
# ---------------------------------------------------------------------------


class TestDeserializeResult:
    """Test AggregatedResultClass and _deserialize_result hook."""

    def test_aggregated_result_class_set(self):
        from polyzymd.analyses.contacts._aggregator import AggregatedContactResult
        from polyzymd.analyses.contacts import ContactsAnalysis

        assert ContactsAnalysis.AggregatedResultClass is AggregatedContactResult

    def test_deserialize_delegates_to_load(self, tmp_path):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()

        mock_result = MagicMock()
        with patch.object(analysis.AggregatedResultClass, "load", return_value=mock_result):
            result = analysis._deserialize_result(tmp_path / "test.json")

        assert result is mock_result


# ---------------------------------------------------------------------------
# plot delegation
# ---------------------------------------------------------------------------


_PLOT_FUNCTIONS = [
    "_plot_contact_fraction_profile",
    "_plot_residence_time_profile",
    "_plot_cf_by_aa_class_bars",
    "_plot_cf_by_partition_bars",
    "_plot_rt_by_aa_class_bars",
    "_plot_rt_by_partition_bars",
    "_plot_user_partition_bars",
    "_plot_system_coverage_bars",
    "_plot_system_coverage_heatmap",
    "_plot_binding_preference_bars",
    "_plot_binding_preference_heatmap",
]


class TestPlot:
    """Test plot() delegates to private module-level plotting functions."""

    def test_plot_creates_output_dir(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.compare.config import PlotSettings

        analysis = ContactsAnalysis()

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
            analysis_dirs={"A": tmp_path / "A" / "contacts"},
            results_dir=tmp_path / "results",
            output_dir=output_dir,
            settings=ContactsSettings(),
            plot_settings=PlotSettings(),
        )

        # Mock all 11 private plot functions to return empty lists
        patches = {
            fn: patch(f"polyzymd.analyses.contacts.{fn}", return_value=[]) for fn in _PLOT_FUNCTIONS
        }
        mocks = {name: p.start() for name, p in patches.items()}
        try:
            analysis.plot(ctx)
        finally:
            for p in patches.values():
                p.stop()

        assert output_dir.exists()
        # Verify all 11 functions were called
        for name, mock_fn in mocks.items():
            assert mock_fn.called, f"{name} was not called"

    def test_plot_returns_combined_paths(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.compare.config import PlotSettings

        analysis = ContactsAnalysis()

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
            analysis_dirs={"A": tmp_path / "A" / "contacts"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=ContactsSettings(),
            plot_settings=PlotSettings(),
        )

        # Make first two functions return paths, rest return empty
        path_a = tmp_path / "plot_a.png"
        path_b = tmp_path / "plot_b.png"

        patches = {}
        for fn in _PLOT_FUNCTIONS:
            patches[fn] = patch(f"polyzymd.analyses.contacts.{fn}", return_value=[])

        patches[_PLOT_FUNCTIONS[0]] = patch(
            f"polyzymd.analyses.contacts.{_PLOT_FUNCTIONS[0]}",
            return_value=[path_a],
        )
        patches[_PLOT_FUNCTIONS[1]] = patch(
            f"polyzymd.analyses.contacts.{_PLOT_FUNCTIONS[1]}",
            return_value=[path_b],
        )

        for p in patches.values():
            p.start()
        try:
            result = analysis.plot(ctx)
        finally:
            for p in patches.values():
                p.stop()

        assert path_a in result
        assert path_b in result
        assert len(result) == 2

    def test_plot_catches_plotter_exceptions(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.compare.config import PlotSettings

        analysis = ContactsAnalysis()

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
            analysis_dirs={"A": tmp_path / "A" / "contacts"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=ContactsSettings(),
            plot_settings=PlotSettings(),
        )

        # Make all plot functions raise
        patches = {
            fn: patch(
                f"polyzymd.analyses.contacts.{fn}",
                side_effect=RuntimeError("plot failed"),
            )
            for fn in _PLOT_FUNCTIONS
        }
        for p in patches.values():
            p.start()
        try:
            result = analysis.plot(ctx)
        finally:
            for p in patches.values():
                p.stop()

        # Should return empty list, not raise
        assert result == []

    def test_plot_passes_plot_settings(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
        from polyzymd.compare.config import PlotSettings

        analysis = ContactsAnalysis()
        ps = PlotSettings()

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
            analysis_dirs={"A": tmp_path / "A" / "contacts"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=ContactsSettings(),
            plot_settings=ps,
        )

        patches = {
            fn: patch(f"polyzymd.analyses.contacts.{fn}", return_value=[]) for fn in _PLOT_FUNCTIONS
        }
        mocks = {name: p.start() for name, p in patches.items()}
        try:
            analysis.plot(ctx)
        finally:
            for p in patches.values():
                p.stop()

        # Verify PlotSettings was passed as 4th positional arg to each function
        for name, mock_fn in mocks.items():
            assert mock_fn.call_args[0][3] is ps, f"{name} did not receive PlotSettings as 4th arg"


# ---------------------------------------------------------------------------
# Binding preference sub-pipeline
# ---------------------------------------------------------------------------


class TestBindingPreference:
    """Test _load_or_compute_binding_preference and helpers."""

    def test_returns_none_when_disabled_and_no_cached(self, tmp_path):
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings

        analysis = ContactsAnalysis()
        settings = ContactsSettings(compute_binding_preference=False)

        mock_sim = MagicMock()
        conditions = [
            Condition(
                label=label,
                config_path=Path(f"/tmp/{label}/config.yaml"),
                replicates=(1, 2),
                sim_config=mock_sim,
            )
            for label in ["A", "B"]
        ]
        analysis_dirs = {
            "A": tmp_path / "A" / "contacts",
            "B": tmp_path / "B" / "contacts",
        }
        for d in analysis_dirs.values():
            d.mkdir(parents=True, exist_ok=True)

        ctx = ComparisonContext(
            name="test",
            conditions=conditions,
            excluded_conditions=[],
            control_label=None,
            analysis_dirs=analysis_dirs,
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=settings,
            recompute=False,
        )

        condition_data = [(cond, {"agg_result": MagicMock()}) for cond in conditions]

        result = analysis._load_or_compute_binding_preference(ctx, condition_data)
        assert result is None

    def test_try_load_cached_bp_not_found(self, tmp_path):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.shared.binding_preference_helpers import (
            try_load_cached_binding_preference,
        )

        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=MagicMock(),
        )

        # Empty directory — no cached results
        result = try_load_cached_binding_preference(cond, tmp_path)
        assert result is None


# ---------------------------------------------------------------------------
# Binding preference pairwise p-values
# ---------------------------------------------------------------------------


class TestBPPairwisePValues:
    """Test _compute_bp_pairwise_pvalues."""

    def test_returns_empty_with_one_condition(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        result = ContactsAnalysis._compute_bp_pairwise_pvalues({"A": [1.0, 2.0, 3.0]})
        assert result == {}

    def test_returns_pvalues_with_two_conditions(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        data = {
            "A": [1.0, 1.1, 1.2, 0.9, 1.05],
            "B": [2.0, 2.1, 2.2, 1.9, 2.05],
        }
        result = ContactsAnalysis._compute_bp_pairwise_pvalues(data)
        assert "A_vs_B" in result
        assert 0.0 <= result["A_vs_B"] <= 1.0

    def test_skips_insufficient_data(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        data = {
            "A": [1.0],  # Only 1 value — not enough for t-test
            "B": [2.0, 2.1, 2.2],
        }
        result = ContactsAnalysis._compute_bp_pairwise_pvalues(data)
        assert result == {}  # Skipped due to insufficient data


# ---------------------------------------------------------------------------
# Integration: full lifecycle
# ---------------------------------------------------------------------------


class TestLifecycle:
    """Verify the plugin can be instantiated and key methods are callable."""

    def test_instantiation(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        assert analysis.name == "contacts"
        assert repr(analysis) == "<ContactsAnalysis(name='contacts')>"

    def test_extract_metrics_returns_empty(self):
        """contacts overrides compare() entirely, so extract_metrics() returns {}."""
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        result = analysis.extract_metrics(MagicMock())
        assert result == {}

    def test_is_subclass_of_analysis(self):
        from polyzymd.analyses.base import Analysis
        from polyzymd.analyses.contacts import ContactsAnalysis

        assert issubclass(ContactsAnalysis, Analysis)

    def test_has_all_required_methods(self):
        from polyzymd.analyses.contacts import ContactsAnalysis

        analysis = ContactsAnalysis()
        assert callable(analysis.compute_replicate)
        assert callable(analysis.aggregate)
        assert callable(analysis.compare)
        assert callable(analysis.plot)
        assert callable(analysis.filter_conditions)
