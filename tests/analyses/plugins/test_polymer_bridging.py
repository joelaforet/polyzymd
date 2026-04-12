"""Tests for the polymer bridging analysis plugin."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from pydantic import BaseModel

from polyzymd.analyses.base import AggregateContext, Condition, PlotContext, ReplicateContext
from polyzymd.analyses.polymer_bridging import (
    PolymerBridgingAggregatedResult,
    PolymerBridgingAnalysis,
    PolymerBridgingObservation,
    PolymerBridgingReplicateResult,
    PolymerBridgingSettings,
    _compute_bridging_statistics_from_frames,
)


class TestDiscovery:
    def test_discovery_by_name(self):
        from polyzymd.analyses.discovery import get_analysis

        assert get_analysis("polymer_bridging") is PolymerBridgingAnalysis

    def test_discovery_by_alias(self):
        from polyzymd.analyses.discovery import get_analysis

        assert get_analysis("bridging") is PolymerBridgingAnalysis


class TestSettings:
    def test_defaults(self):
        settings = PolymerBridgingSettings()
        assert settings.polymer_selection == "chainID C"
        assert settings.min_ca_distance_angstrom == pytest.approx(0.0)

    def test_negative_distance_rejected(self):
        with pytest.raises(ValueError):
            PolymerBridgingSettings(min_ca_distance_angstrom=-1.0)


class TestCacheFingerprint:
    """Verify cache filenames include settings fingerprint."""

    def test_cache_tag_changes_with_settings(self):
        """Different settings must produce different cache tags."""
        analysis = PolymerBridgingAnalysis()
        tag_default = analysis._make_settings_cache_tag(PolymerBridgingSettings())
        tag_custom = analysis._make_settings_cache_tag(
            PolymerBridgingSettings(cutoff=5.0, min_ca_distance_angstrom=10.0)
        )
        assert tag_default != tag_custom
        assert len(tag_default) == 8
        assert len(tag_custom) == 8


class TestCoreComputation:
    def test_bridging_stats_without_distance_threshold(self):
        frame_contacts = [{10}, {10, 35}, {10, 35}, {60}]
        stats = _compute_bridging_statistics_from_frames(
            frame_contacts,
            min_ca_distance_angstrom=0.0,
            ca_distances={},
        )
        assert stats["contacting_observations"] == 4
        assert stats["multisite_observations"] == 2
        assert stats["high_valency_observations"] == 0
        assert stats["multisite_fraction"] == pytest.approx(0.5)
        assert stats["mean_contacts_per_contacting_oligomer"] == pytest.approx(1.5)

    def test_bridging_stats_with_distance_threshold(self):
        frame_contacts = [{10}, {10, 35}, {10, 35}, {60}]
        stats = _compute_bridging_statistics_from_frames(
            frame_contacts,
            min_ca_distance_angstrom=100.0,
            ca_distances={(10, 35): 95.0, (10, 60): 190.0, (35, 60): 95.0},
        )
        assert stats["multisite_fraction"] == pytest.approx(0.0)
        assert stats["mean_contacts_per_contacting_oligomer"] == pytest.approx(1.0)

    def test_bridging_stats_use_dynamic_observation_distances(self):
        frame_contacts = [
            ({10, 35}, {(10, 35): 7.0}),
            ({10, 35}, {(10, 35): 9.0}),
            ({60}, {}),
        ]
        stats = _compute_bridging_statistics_from_frames(
            frame_contacts,
            min_ca_distance_angstrom=8.0,
        )
        assert stats["contacting_observations"] == 3
        assert stats["multisite_observations"] == 1
        assert stats["multisite_fraction"] == pytest.approx(1.0 / 3.0)
        assert stats["mean_contacts_per_contacting_oligomer"] == pytest.approx(4.0 / 3.0)

    def test_bridging_stats_capture_anchor_and_signature_metadata(self):
        observation = PolymerBridgingObservation(
            protein_residues={10, 35, 60},
            protein_resnames={10: "PHE", 35: "SER", 60: "LEU"},
            protein_groups={10: "aromatic", 35: "polar", 60: "nonpolar"},
            contacting_polymer_resids={101, 103},
            polymer_resnames={101: "SBM", 102: "EGM", 103: "SBM", 104: "EGM", 105: "SBM"},
            fragment_signature=("SBM", "EGM", "SBM", "EGM", "SBM"),
            ca_distances={(10, 35): 9.0, (10, 60): 12.0, (35, 60): 7.0},
            pair_min_distances={(101, 10): 3.2, (103, 35): 3.8, (103, 60): 4.2},
        )
        stats = _compute_bridging_statistics_from_frames(
            [observation],
            min_ca_distance_angstrom=8.0,
        )

        assert stats["anchor_protein_group_probabilities"]["aromatic"] == pytest.approx(1.0)
        assert stats["polymer_anchor_type_probabilities"]["SBM"] == pytest.approx(1.0)
        assert stats["multivalent_protein_group_probabilities"]["nonpolar"] > 0.0
        assert "SBM-EGM-SBM-EGM-SBM" in stats["fragment_signature_probabilities"]


class TestCAValidation:
    """Verify loud failure when CA atoms are missing for distance filtering."""

    def test_no_ca_atoms_with_distance_filter_raises(self):
        """Raise ValueError when CA filtering is requested but no CA atoms exist."""
        condition = MagicMock()
        condition.sim_config = MagicMock()

        mock_loader = MagicMock()
        mock_universe = MagicMock()

        atom = MagicMock()
        atom.resid = 1
        atom.resname = "ALA"

        mock_protein = MagicMock()
        mock_protein.__len__.return_value = 1
        mock_protein.atoms = [atom]
        mock_residue = MagicMock()
        mock_residue.resid = 1
        mock_ca_selection = MagicMock()
        mock_ca_selection.__len__.return_value = 0
        mock_residue.atoms.select_atoms.return_value = mock_ca_selection
        mock_protein.residues = [mock_residue]

        mock_polymer = MagicMock()
        mock_polymer.__len__.return_value = 1

        mock_universe.select_atoms.side_effect = [mock_protein, mock_polymer]
        mock_loader.load_universe.return_value = mock_universe

        with patch(
            "polyzymd.analyses.shared.loader.TrajectoryLoader",
            return_value=mock_loader,
        ):
            with pytest.raises(ValueError, match="(?i)no CA atoms"):
                PolymerBridgingAnalysis._compute_frame_contacts(
                    condition,
                    1,
                    protein_selection="protein",
                    polymer_selection="chainID C",
                    cutoff=4.5,
                    equilibration="0ns",
                    min_ca_distance_angstrom=8.0,
                )

    def test_no_ca_atoms_without_distance_filter_succeeds(self):
        """Allow execution when CA filtering is disabled."""
        condition = MagicMock()
        condition.sim_config = MagicMock()

        mock_loader = MagicMock()
        mock_universe = MagicMock()

        atom = MagicMock()
        atom.resid = 1
        atom.resname = "ALA"

        mock_protein = MagicMock()
        mock_protein.__len__.return_value = 1
        mock_protein.atoms = [atom]
        mock_residue = MagicMock()
        mock_residue.resid = 1
        mock_ca_selection = MagicMock()
        mock_ca_selection.__len__.return_value = 0
        mock_residue.atoms.select_atoms.return_value = mock_ca_selection
        mock_protein.residues = [mock_residue]

        mock_polymer = MagicMock()
        mock_polymer.__len__.return_value = 1
        mock_polymer.fragments = []

        mock_universe.select_atoms.side_effect = [mock_protein, mock_polymer]
        mock_universe.trajectory = []
        mock_loader.load_universe.return_value = mock_universe
        mock_loader.get_timestep.return_value = 10.0

        with patch(
            "polyzymd.analyses.shared.loader.TrajectoryLoader",
            return_value=mock_loader,
        ):
            observations, n_frames, timestep_ps = PolymerBridgingAnalysis._compute_frame_contacts(
                condition,
                1,
                protein_selection="protein",
                polymer_selection="chainID C",
                cutoff=4.5,
                equilibration="0ns",
                min_ca_distance_angstrom=0.0,
            )

        assert observations == []
        assert n_frames == 0
        assert timestep_ps == pytest.approx(10.0)


class TestLifecycle:
    def test_compute_replicate_uses_trajectory_contacts(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        condition = Condition(
            label="Cond",
            config_path=tmp_path / "condition" / "config.yaml",
            replicates=(1,),
            sim_config=MagicMock(),
        )
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "out",
            equilibration="0ns",
            recompute=False,
            settings=PolymerBridgingSettings(),
        )
        analysis._compute_frame_contacts = MagicMock(
            return_value=(
                [
                    ({10}, {}),
                    ({10, 35}, {(10, 35): 20.0}),
                    ({10, 35}, {(10, 35): 20.0}),
                    ({60}, {}),
                ],
                4,
                10.0,
            )
        )

        result = analysis.compute_replicate(ctx, 1)

        assert result.replicate == 1
        assert result.multisite_fraction == pytest.approx(0.5)

    def test_compute_replicate_uses_cached_result_when_available(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        condition = Condition(
            label="Cond",
            config_path=tmp_path / "condition" / "config.yaml",
            replicates=(1,),
            sim_config=MagicMock(),
        )
        settings_tag = analysis._make_settings_cache_tag(settings)
        result_path = tmp_path / "out" / f"polymer_bridging_{settings_tag}.json"
        cached = PolymerBridgingReplicateResult(
            replicate=1,
            n_frames=4,
            timestep_ps=10.0,
            min_ca_distance_angstrom=0.0,
            contacting_observations=4,
            multisite_observations=2,
            high_valency_observations=0,
            mean_contacts_per_contacting_oligomer=1.5,
            multisite_fraction=0.5,
            high_valency_fraction=0.0,
            valency_probabilities={"1": 0.5, "2": 0.5, "3+": 0.0},
        )
        result_path.parent.mkdir(parents=True, exist_ok=True)
        cached.save(result_path)

        analysis._compute_frame_contacts = MagicMock(side_effect=AssertionError("should not run"))

        loaded = analysis.compute_replicate(
            ReplicateContext(
                condition=condition,
                replicate=1,
                sim_config=condition.sim_config,
                output_dir=tmp_path / "out",
                equilibration="0ns",
                recompute=False,
                settings=settings,
            ),
            1,
        )

        assert loaded == cached

    def test_aggregate_builds_typed_result(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        condition = Condition(
            label="Cond",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2),
            output_dir=tmp_path / "agg",
            equilibration="0ns",
            settings=PolymerBridgingSettings(),
            result_path=tmp_path / "agg" / "result.json",
        )
        analysis._compute_frame_contacts = MagicMock(
            return_value=(
                [
                    ({10}, {}),
                    ({10, 35}, {(10, 35): 20.0}),
                    ({10, 35}, {(10, 35): 20.0}),
                    ({60}, {}),
                ],
                4,
                10.0,
            )
        )
        results = [
            analysis.compute_replicate(
                ReplicateContext(
                    condition=condition,
                    replicate=1,
                    sim_config=condition.sim_config,
                    output_dir=tmp_path / "r1",
                    equilibration="0ns",
                    recompute=False,
                    settings=PolymerBridgingSettings(),
                ),
                1,
            ),
            analysis.compute_replicate(
                ReplicateContext(
                    condition=condition,
                    replicate=2,
                    sim_config=condition.sim_config,
                    output_dir=tmp_path / "r2",
                    equilibration="0ns",
                    recompute=False,
                    settings=PolymerBridgingSettings(),
                ),
                2,
            ),
        ]

        aggregated = analysis.aggregate(ctx, results)

        assert isinstance(aggregated, PolymerBridgingAggregatedResult)
        assert aggregated.n_replicates == 2

    def test_aggregate_uses_cached_result_when_available(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()
        condition = Condition(
            label="Cond",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        settings_tag = analysis._make_settings_cache_tag(settings)
        rep_str = analysis._format_replicate_range((1, 2))
        result_path = tmp_path / "agg" / f"polymer_bridging_{rep_str}_{settings_tag}.json"
        cached = PolymerBridgingAggregatedResult(
            n_replicates=2,
            replicates=[1, 2],
            min_ca_distance_angstrom=0.0,
            mean_contacts_per_contacting_oligomer=1.5,
            mean_contacts_sem=0.1,
            multisite_fraction=0.5,
            multisite_fraction_sem=0.05,
            high_valency_fraction=0.1,
            high_valency_fraction_sem=0.02,
            mean_contacts_per_contacting_oligomer_replicates=[1.4, 1.6],
            multisite_fraction_replicates=[0.45, 0.55],
            high_valency_fraction_replicates=[0.08, 0.12],
            valency_probabilities_mean={"1": 0.5, "2": 0.4, "3+": 0.1},
            valency_probabilities_sem={"1": 0.01, "2": 0.02, "3+": 0.01},
            valency_probabilities_per_replicate={
                "1": [0.52, 0.48],
                "2": [0.38, 0.42],
                "3+": [0.1, 0.1],
            },
        )
        result_path.parent.mkdir(parents=True, exist_ok=True)
        cached.save(result_path)

        loaded = analysis.aggregate(
            AggregateContext(
                condition=condition,
                replicates=(1, 2),
                output_dir=tmp_path / "agg",
                equilibration="0ns",
                settings=settings,
            ),
            [],
        )

        assert loaded == cached

    def test_plot_returns_paths(self, tmp_path):
        analysis = PolymerBridgingAnalysis()
        aggregated = PolymerBridgingAggregatedResult(
            n_replicates=2,
            replicates=[1, 2],
            min_ca_distance_angstrom=0.0,
            mean_contacts_per_contacting_oligomer=1.5,
            mean_contacts_sem=0.1,
            multisite_fraction=0.5,
            multisite_fraction_sem=0.05,
            high_valency_fraction=0.1,
            high_valency_fraction_sem=0.02,
            mean_contacts_per_contacting_oligomer_replicates=[1.4, 1.6],
            multisite_fraction_replicates=[0.45, 0.55],
            high_valency_fraction_replicates=[0.08, 0.12],
            valency_probabilities_mean={"1": 0.5, "2": 0.4, "3+": 0.1},
            valency_probabilities_sem={"1": 0.01, "2": 0.02, "3+": 0.01},
            valency_probabilities_per_replicate={
                "1": [0.52, 0.48],
                "2": [0.38, 0.42],
                "3+": [0.1, 0.1],
            },
        )
        comparison_json = tmp_path / "comparison" / "polymer_bridging" / "result.json"
        comparison_json.parent.mkdir(parents=True, exist_ok=True)
        comparison_json.write_text(
            '{"analysis_type":"polymer_bridging","name":"Test","control_label":null,'
            '"conditions":[{"label":"A","n_replicates":2,'
            '"multisite_fraction_mean":0.5,"multisite_fraction_sem":0.05,'
            '"multisite_fraction_replicate_values":[0.45,0.55],'
            '"mean_contacts_per_contacting_oligomer_mean":1.5,'
            '"mean_contacts_per_contacting_oligomer_sem":0.1,'
            '"mean_contacts_per_contacting_oligomer_replicate_values":[1.4,1.6],'
            '"high_valency_fraction_mean":0.1,"high_valency_fraction_sem":0.02,'
            '"high_valency_fraction_replicate_values":[0.08,0.12]},'
            '{"label":"B","n_replicates":2,'
            '"multisite_fraction_mean":0.7,"multisite_fraction_sem":0.04,'
            '"multisite_fraction_replicate_values":[0.66,0.74],'
            '"mean_contacts_per_contacting_oligomer_mean":1.8,'
            '"mean_contacts_per_contacting_oligomer_sem":0.08,'
            '"mean_contacts_per_contacting_oligomer_replicate_values":[1.72,1.88],'
            '"high_valency_fraction_mean":0.2,"high_valency_fraction_sem":0.03,'
            '"high_valency_fraction_replicate_values":[0.17,0.23]}],'
            '"pairwise_comparisons":[],"anova":null,"ranking":["B","A"],'
            '"rankings_by_metric":{"multisite_fraction":["B","A"],'
            '"mean_contacts_per_contacting_oligomer":["B","A"],'
            '"high_valency_fraction":["B","A"]},"equilibration_time":"0ns",'
            '"created_at":"now","polyzymd_version":"test"}'
        )
        cond_a = Condition(
            label="A", config_path=Path("/tmp/a.yaml"), replicates=(1, 2), sim_config=MagicMock()
        )
        cond_b = Condition(
            label="B", config_path=Path("/tmp/b.yaml"), replicates=(1, 2), sim_config=MagicMock()
        )
        analysis_dir_a = tmp_path / "analysis" / "A" / "polymer_bridging" / "aggregated"
        analysis_dir_b = tmp_path / "analysis" / "B" / "polymer_bridging" / "aggregated"
        analysis_dir_a.mkdir(parents=True, exist_ok=True)
        analysis_dir_b.mkdir(parents=True, exist_ok=True)
        aggregated.save(analysis_dir_a / "result.json")
        aggregated.model_copy(
            update={
                "multisite_fraction": 0.7,
                "mean_contacts_per_contacting_oligomer": 1.8,
                "high_valency_fraction": 0.2,
                "valency_probabilities_mean": {"1": 0.3, "2": 0.5, "3+": 0.2},
            }
        ).save(analysis_dir_b / "result.json")
        from polyzymd.config.comparison import PlotSettings

        ctx = PlotContext(
            conditions=[cond_a, cond_b],
            analysis_dirs={"A": analysis_dir_a.parent, "B": analysis_dir_b.parent},
            results_dir=comparison_json.parent,
            output_dir=tmp_path / "figures",
            settings=PolymerBridgingSettings(),
            plot_settings=PlotSettings(),
            comparison_path=comparison_json,
        )

        paths = analysis.plot(ctx)

        assert len(paths) >= 3
        assert all(path.exists() for path in paths)

    def test_format_emits_sections(self):
        analysis = PolymerBridgingAnalysis()
        from polyzymd.analyses.base import ComparisonResult, ConditionSummary

        result = ComparisonResult(
            analysis_type="polymer_bridging",
            name="Test",
            conditions=[
                ConditionSummary(
                    label="A",
                    n_replicates=2,
                    multisite_fraction_mean=0.5,
                    multisite_fraction_sem=0.05,
                    multisite_fraction_replicate_values=[0.45, 0.55],
                    mean_contacts_per_contacting_oligomer_mean=1.5,
                    mean_contacts_per_contacting_oligomer_sem=0.1,
                    mean_contacts_per_contacting_oligomer_replicate_values=[1.4, 1.6],
                    high_valency_fraction_mean=0.1,
                    high_valency_fraction_sem=0.02,
                    high_valency_fraction_replicate_values=[0.08, 0.12],
                ),
                ConditionSummary(
                    label="B",
                    n_replicates=2,
                    multisite_fraction_mean=0.7,
                    multisite_fraction_sem=0.04,
                    multisite_fraction_replicate_values=[0.66, 0.74],
                    mean_contacts_per_contacting_oligomer_mean=1.8,
                    mean_contacts_per_contacting_oligomer_sem=0.08,
                    mean_contacts_per_contacting_oligomer_replicate_values=[1.72, 1.88],
                    high_valency_fraction_mean=0.2,
                    high_valency_fraction_sem=0.03,
                    high_valency_fraction_replicate_values=[0.17, 0.23],
                ),
            ],
            pairwise_comparisons=[],
            anova=None,
            ranking=["B", "A"],
            rankings_by_metric={
                "multisite_fraction": ["B", "A"],
                "mean_contacts_per_contacting_oligomer": ["B", "A"],
                "high_valency_fraction": ["B", "A"],
            },
            equilibration_time="0ns",
            created_at="now",
            polyzymd_version="test",
        )

        text = analysis.format(result, output_format="text")

        assert "WARNING: Experimental analysis" in text
        assert "Polymer Bridging Comparison" in text
        assert "Average Oligomer Valency" in text
        assert "High-Valency Oligomers" in text

    def test_compare_sanitizes_nan_stats(self):
        analysis = PolymerBridgingAnalysis()
        condition_a = Condition(
            label="A",
            config_path=Path("/tmp/a.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        condition_b = Condition(
            label="B",
            config_path=Path("/tmp/b.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        aggregated_a = PolymerBridgingAggregatedResult(
            n_replicates=3,
            replicates=[1, 2, 3],
            min_ca_distance_angstrom=0.0,
            mean_contacts_per_contacting_oligomer=1.0,
            mean_contacts_sem=0.0,
            multisite_fraction=0.0,
            multisite_fraction_sem=0.0,
            high_valency_fraction=0.0,
            high_valency_fraction_sem=0.0,
            mean_contacts_per_contacting_oligomer_replicates=[1.0, 1.0, 1.0],
            multisite_fraction_replicates=[0.0, 0.0, 0.0],
            high_valency_fraction_replicates=[0.0, 0.0, 0.0],
            valency_probabilities_mean={"1": 1.0, "2": 0.0, "3+": 0.0},
            valency_probabilities_sem={"1": 0.0, "2": 0.0, "3+": 0.0},
            valency_probabilities_per_replicate={
                "1": [1.0, 1.0, 1.0],
                "2": [0.0, 0.0, 0.0],
                "3+": [0.0, 0.0, 0.0],
            },
        )
        aggregated_b = aggregated_a.model_copy(
            update={
                "mean_contacts_per_contacting_oligomer": 2.0,
                "mean_contacts_per_contacting_oligomer_replicates": [2.0, 2.0, 2.0],
                "multisite_fraction": 1.0,
                "multisite_fraction_replicates": [1.0, 1.0, 1.0],
                "high_valency_fraction": 1.0,
                "high_valency_fraction_replicates": [1.0, 1.0, 1.0],
                "valency_probabilities_mean": {"1": 0.0, "2": 0.0, "3+": 1.0},
                "valency_probabilities_per_replicate": {
                    "1": [0.0, 0.0, 0.0],
                    "2": [0.0, 0.0, 0.0],
                    "3+": [1.0, 1.0, 1.0],
                },
            }
        )
        from polyzymd.analyses.base import ComparisonContext

        ctx = ComparisonContext(
            name="Test",
            conditions=[condition_a, condition_b],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"A": Path("/tmp/a"), "B": Path("/tmp/b")},
            results_dir=Path("/tmp/results"),
            equilibration="0ns",
            settings=PolymerBridgingSettings(),
            recompute=False,
            aggregated_results={"A": aggregated_a, "B": aggregated_b},
        )

        result = analysis.compare(ctx)

        assert result is not None
        assert all(pair.p_value >= 0.0 for pair in result.pairwise_comparisons)

    def test_compare_passes_welch_method_to_pairwise_helper(self, monkeypatch):
        analysis = PolymerBridgingAnalysis()
        condition_a = Condition(
            label="A",
            config_path=Path("/tmp/a.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        condition_b = Condition(
            label="B",
            config_path=Path("/tmp/b.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        aggregated_a = PolymerBridgingAggregatedResult(
            n_replicates=2,
            replicates=[1, 2],
            min_ca_distance_angstrom=0.0,
            mean_contacts_per_contacting_oligomer=1.0,
            mean_contacts_sem=0.1,
            multisite_fraction=0.4,
            multisite_fraction_sem=0.05,
            high_valency_fraction=0.1,
            high_valency_fraction_sem=0.02,
            mean_contacts_per_contacting_oligomer_replicates=[0.9, 1.1],
            multisite_fraction_replicates=[0.35, 0.45],
            high_valency_fraction_replicates=[0.09, 0.11],
            valency_probabilities_mean={"1": 0.6, "2": 0.3, "3+": 0.1},
            valency_probabilities_sem={"1": 0.01, "2": 0.01, "3+": 0.01},
            valency_probabilities_per_replicate={
                "1": [0.59, 0.61],
                "2": [0.31, 0.29],
                "3+": [0.1, 0.1],
            },
        )
        aggregated_b = aggregated_a.model_copy(
            update={
                "mean_contacts_per_contacting_oligomer": 1.4,
                "mean_contacts_per_contacting_oligomer_replicates": [1.3, 1.5],
                "multisite_fraction": 0.6,
                "multisite_fraction_replicates": [0.55, 0.65],
                "high_valency_fraction": 0.2,
                "high_valency_fraction_replicates": [0.18, 0.22],
            }
        )

        captured: dict[str, object] = {}

        def _fake_pairwise(metrics_by_condition, control_label, **kwargs):
            captured["ttest_method"] = kwargs.get("ttest_method")
            captured["posthoc_method"] = kwargs.get("posthoc_method")
            captured["fdr_alpha"] = kwargs.get("fdr_alpha")
            return []

        monkeypatch.setattr(
            "polyzymd.analyses.polymer_bridging.pairwise_comparisons",
            _fake_pairwise,
        )

        from polyzymd.analyses.base import ComparisonContext

        ctx = ComparisonContext(
            name="Test",
            conditions=[condition_a, condition_b],
            excluded_conditions=[],
            control_label="A",
            analysis_dirs={"A": Path("/tmp/a"), "B": Path("/tmp/b")},
            results_dir=Path("/tmp/results"),
            equilibration="0ns",
            settings=PolymerBridgingSettings(),
            recompute=False,
            aggregated_results={"A": aggregated_a, "B": aggregated_b},
            ttest_method="welch",
            posthoc_method="ttest_bh",
            fdr_alpha=0.1,
        )

        result = analysis.compare(ctx)

        assert result is not None
        assert captured["ttest_method"] == "welch"
        assert captured["posthoc_method"] == "ttest_bh"
        assert captured["fdr_alpha"] == pytest.approx(0.1)
        assert result.ttest_method == "welch"

    def test_compare_honors_tukey_hsd_posthoc_method(self):
        analysis = PolymerBridgingAnalysis()
        condition_a = Condition(
            label="A",
            config_path=Path("/tmp/a.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        condition_b = Condition(
            label="B",
            config_path=Path("/tmp/b.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        condition_c = Condition(
            label="C",
            config_path=Path("/tmp/c.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        aggregated_a = PolymerBridgingAggregatedResult(
            n_replicates=3,
            replicates=[1, 2, 3],
            min_ca_distance_angstrom=0.0,
            mean_contacts_per_contacting_oligomer=1.0,
            mean_contacts_sem=0.1,
            multisite_fraction=0.2,
            multisite_fraction_sem=0.02,
            high_valency_fraction=0.1,
            high_valency_fraction_sem=0.01,
            mean_contacts_per_contacting_oligomer_replicates=[0.9, 1.0, 1.1],
            multisite_fraction_replicates=[0.18, 0.2, 0.22],
            high_valency_fraction_replicates=[0.08, 0.1, 0.12],
            valency_probabilities_mean={"1": 0.7, "2": 0.2, "3+": 0.1},
            valency_probabilities_sem={"1": 0.01, "2": 0.01, "3+": 0.01},
            valency_probabilities_per_replicate={
                "1": [0.7, 0.69, 0.71],
                "2": [0.2, 0.21, 0.19],
                "3+": [0.1, 0.1, 0.1],
            },
        )
        aggregated_b = aggregated_a.model_copy(
            update={
                "mean_contacts_per_contacting_oligomer": 1.4,
                "mean_contacts_per_contacting_oligomer_replicates": [1.3, 1.4, 1.5],
                "multisite_fraction": 0.5,
                "multisite_fraction_replicates": [0.48, 0.5, 0.52],
                "high_valency_fraction": 0.2,
                "high_valency_fraction_replicates": [0.18, 0.2, 0.22],
            }
        )
        aggregated_c = aggregated_a.model_copy(
            update={
                "mean_contacts_per_contacting_oligomer": 1.8,
                "mean_contacts_per_contacting_oligomer_replicates": [1.7, 1.8, 1.9],
                "multisite_fraction": 0.8,
                "multisite_fraction_replicates": [0.78, 0.8, 0.82],
                "high_valency_fraction": 0.4,
                "high_valency_fraction_replicates": [0.38, 0.4, 0.42],
            }
        )

        from polyzymd.analyses.base import ComparisonContext

        ctx = ComparisonContext(
            name="Test",
            conditions=[condition_a, condition_b, condition_c],
            excluded_conditions=[],
            control_label="A",
            analysis_dirs={"A": Path("/tmp/a"), "B": Path("/tmp/b"), "C": Path("/tmp/c")},
            results_dir=Path("/tmp/results"),
            equilibration="0ns",
            settings=PolymerBridgingSettings(),
            recompute=False,
            aggregated_results={"A": aggregated_a, "B": aggregated_b, "C": aggregated_c},
            ttest_method="welch",
            posthoc_method="tukey_hsd",
            fdr_alpha=0.05,
        )

        result = analysis.compare(ctx)

        assert result is not None
        assert result.posthoc_method == "tukey_hsd"
        assert result.ttest_method == "welch"
        assert result.pairwise_comparisons
        assert all(pair.posthoc_method == "tukey_hsd" for pair in result.pairwise_comparisons)


class TestConfigCompatibility:
    def test_legacy_analysis_settings_promoted_to_plugins(self, tmp_path):
        yaml_path = tmp_path / "comparison.yaml"
        cond_dir = tmp_path / "cond"
        cond_dir.mkdir()
        (cond_dir / "config.yaml").write_text("name: test\n")
        yaml_path.write_text(
            "name: Legacy\n"
            "defaults:\n  equilibration_time: 0ns\n"
            "conditions:\n"
            f"  - label: A\n    config: {cond_dir / 'config.yaml'}\n    replicates: [1, 2]\n"
            f"  - label: B\n    config: {cond_dir / 'config.yaml'}\n    replicates: [1, 2]\n"
            "analysis_settings:\n"
            "  polymer_bridging:\n"
            "    min_ca_distance_angstrom: 12.0\n"
        )

        from polyzymd.config.comparison import ComparisonConfig

        cfg = ComparisonConfig.from_yaml(yaml_path)

        assert "polymer_bridging" in cfg.plugins.get_enabled_plugins()
        settings = cfg.plugins.get("polymer_bridging")
        assert settings.min_ca_distance_angstrom == pytest.approx(12.0)


class TestSharedPlotSettingsReexport:
    def test_plotsettings_can_be_imported_from_shared(self):
        from polyzymd.analyses.shared import PlotSettings

        assert issubclass(PlotSettings, BaseModel)
