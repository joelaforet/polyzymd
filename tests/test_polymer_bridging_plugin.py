"""Tests for the polymer bridging analysis plugin."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest
from pydantic import BaseModel

from polyzymd.analyses.base import AggregateContext, Condition, PlotContext, ReplicateContext
from polyzymd.analyses.polymer_bridging import (
    PolymerBridgingObservation,
    PolymerBridgingAggregatedResult,
    PolymerBridgingAnalysis,
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
        condition = Condition(
            label="Cond",
            config_path=tmp_path / "condition" / "config.yaml",
            replicates=(1,),
            sim_config=MagicMock(),
        )
        result_path = tmp_path / "out" / "result.json"
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
                settings=PolymerBridgingSettings(),
                result_path=result_path,
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
        condition = Condition(
            label="Cond",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        result_path = tmp_path / "agg" / "result.json"
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
                settings=PolymerBridgingSettings(),
                result_path=result_path,
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
