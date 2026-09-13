"""Tests for the agent-facing analysis protocol."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any, ClassVar, Sequence

import pytest
from pydantic import BaseModel

from polyzymd.analyses.base import Condition
from polyzymd.analyses.contract import Observable, contract_analysis, reduce_replicate
from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.mda import ReplicateArtifact
from polyzymd.analyses.protocols import (
    VERDICT_LARGER,
    VERDICT_NO_DIFFERENCE,
    VERDICT_NOT_TESTABLE,
    ConditionReport,
    PairwiseReport,
    ProtocolReport,
    analyze,
    build_report,
)
from polyzymd.analyses.stats import interpret_direction

# Replicate values the toy plugin reports, keyed by condition label. Tests set
# this before calling analyze().
REPLICATE_VALUES: dict[str, list[float]] = {}


class ToyProtocolSettings(BaseModel):
    """Settings for the toy protocol plugin."""

    scale: float = 1.0


class ToyProtocolPlugin:
    """Contract plugin whose observables the test supplies directly."""

    name = "toy_protocol"
    Settings = ToyProtocolSettings
    references: tuple[str, ...] = ()

    #: Observables reported per replicate; the second is added by the multi
    #: metric variant below.
    observable_names: tuple[str, ...] = ("mean_value",)

    def compute(self, universe: Any, frames: Any, settings: Any) -> list[Observable]:
        """Never called; the test drives the compute stage directly."""
        raise AssertionError("the toy plugin computes through _run_compute_stage")


class ToyMultiMetricPlugin(ToyProtocolPlugin):
    """Toy plugin that reports two observables."""

    name = "toy_multi"
    observable_names: tuple[str, ...] = ("mean_value", "doubled_value")


class ToyProtocolAnalysis(contract_analysis(ToyProtocolPlugin)):  # type: ignore[misc]
    """Analysis that skips the universe and reports the test's values."""

    protocol_version: ClassVar[str] = "1"

    def _run_compute_stage(self, ctx: Any, replicate: int) -> ReplicateArtifact:
        """Write one replicate artifact from the values the test assigned.

        The contract lifecycle owns aggregation, testing and formatting from
        here, so the toy only has to produce the observables a real plugin's
        ``compute()`` would have measured.
        """
        value = float(REPLICATE_VALUES[ctx.condition.label][replicate - 1]) * ctx.settings.scale
        scales = {"mean_value": 1.0, "doubled_value": 2.0}
        observables = [
            Observable(
                name=name,
                kind="mean_of_timeseries",
                values=[value * scales[name]],
                unit="A",
            )
            for name in self.plugin.observable_names
        ]
        return ReplicateArtifact(
            analysis_name=self.name,
            condition_label=ctx.condition.label,
            replicate=replicate,
            payload={
                "observables": [
                    estimate.model_dump(mode="json") for estimate in reduce_replicate(observables)
                ]
            },
        )

    def plot(self, ctx: Any) -> list[Path]:
        """Write a placeholder figure so the pipeline does not need matplotlib."""
        out = ctx.output_dir / "toy_protocol.txt"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("figure")
        return [out]


class ToyMultiMetricAnalysis(ToyProtocolAnalysis):
    """Toy analysis that reports two observables."""

    name: ClassVar[str] = "toy_multi"
    protocol_version: ClassVar[str] = "2"
    plugin = ToyMultiMetricPlugin()


def _install_toy(monkeypatch: pytest.MonkeyPatch, analysis_cls: type = ToyProtocolAnalysis) -> None:
    """Make the toy plugin discoverable and skip simulation config loading.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Patching fixture.
    analysis_cls : type, optional
        Plugin class to register, by default :class:`ToyProtocolAnalysis`.
    """
    registry = {analysis_cls.name: analysis_cls}

    def _get_analysis(name: str) -> type:
        if name in registry:
            return registry[name]
        raise KeyError(f"Unknown analysis {name!r}")

    monkeypatch.setattr("polyzymd.analyses.discovery.get_analysis", _get_analysis)
    monkeypatch.setattr("polyzymd.analyses.discovery.list_all_names", lambda: sorted(registry))
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.Condition.from_condition_config",
        lambda cond: Condition(
            cond.label, Path(cond.config), tuple(cond.replicates), SimpleNamespace()
        ),
    )


def _write_configs(tmp_path: Path, labels: Sequence[str]) -> list[Path]:
    """Create one placeholder simulation config per label.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory.
    labels : sequence of str
        Condition directory names.

    Returns
    -------
    list of Path
        Config paths in the given order.
    """
    paths = []
    for label in labels:
        directory = tmp_path / label
        directory.mkdir(parents=True, exist_ok=True)
        config = directory / "config.yaml"
        config.write_text("placeholder: true\n")
        paths.append(config)
    return paths


def _run(
    tmp_path: Path,
    values: dict[str, list[float]],
    *,
    name: str = "toy_protocol",
    equilibration: str = "10ns",
) -> ProtocolReport:
    """Run the toy protocol over the given per-condition replicate values.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory used for configs and outputs.
    values : dict
        Replicate values keyed by condition label.
    name : str, optional
        Analysis name, by default ``"toy_protocol"``.
    equilibration : str, optional
        Equilibration window, by default ``"10ns"``.

    Returns
    -------
    ProtocolReport
        The report.
    """
    REPLICATE_VALUES.clear()
    REPLICATE_VALUES.update(values)
    labels = list(values)
    configs = _write_configs(tmp_path, labels)
    replicates = list(range(1, len(next(iter(values.values()))) + 1))
    return analyze(
        name,
        configs,
        replicates=replicates,
        equilibration=equilibration,
        output_dir=tmp_path / "out",
    )


class TestReportFields:
    """The report states what every number is."""

    def test_two_conditions_report_is_fully_typed(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """Every documented field is present and carries the declared type."""
        _install_toy(monkeypatch)
        report = _run(
            tmp_path,
            {"A": [10.0, 10.1, 10.2], "B": [12.0, 12.1, 12.2]},
        )

        assert isinstance(report, ProtocolReport)
        assert report.analysis == "toy_protocol"
        assert report.protocol_version == "1"
        assert report.metric == "mean_value"
        assert report.unit == "A"
        assert report.all_metrics == ["mean_value"]
        assert report.equilibration == "10ns"
        assert set(report.frames_per_replicate) == {"A", "B"}

        assert [condition.label for condition in report.conditions] == ["A", "B"]
        first = report.conditions[0]
        assert isinstance(first, ConditionReport)
        assert first.n_replicates == 3
        assert first.mean == pytest.approx(10.1)
        assert first.sem is not None and first.sem > 0.0
        assert first.ci95 is not None and first.ci95[0] < first.mean < first.ci95[1]
        assert first.ci_method == "student_t"
        assert first.replicate_values == [10.0, 10.1, 10.2]

        assert len(report.pairwise) == 1
        pair = report.pairwise[0]
        assert isinstance(pair, PairwiseReport)
        assert (pair.a, pair.b) == ("A", "B")
        assert pair.delta == pytest.approx(2.0)
        assert pair.delta_ci95 is not None
        assert pair.p is not None and pair.p_adjusted is not None
        assert pair.test == "student_t"
        assert pair.correction == "benjamini_hochberg"
        assert pair.cohens_d is not None
        # The report orients the effect size like delta, so both are positive
        # when the second condition is larger.
        assert pair.cohens_d > 0.0
        assert pair.hedges_g is not None and 0.0 < pair.hedges_g < pair.cohens_d
        assert pair.testable is True

        assert report.provenance.polyzymd_version
        assert set(report.provenance.config_hashes) == {"A", "B"}
        assert "comparison_result" in report.provenance.output_paths
        assert report.verdict

    def test_primary_metric_is_the_first_and_the_rest_are_listed(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """A multi-metric plugin reports its first metric and names the others."""
        _install_toy(monkeypatch, ToyMultiMetricAnalysis)
        report = _run(
            tmp_path,
            {"A": [10.0, 10.1, 10.2], "B": [12.0, 12.1, 12.2]},
            name="toy_multi",
        )

        assert report.metric == "mean_value"
        assert report.all_runs == ["mean_value", "doubled_value"]
        assert report.protocol_version == "2"
        assert report.conditions[0].mean == pytest.approx(10.1)
        assert len(report.pairwise) == 1

    def test_json_round_trips_through_the_model(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """The JSON form validates back into an equal report."""
        _install_toy(monkeypatch)
        report = _run(tmp_path, {"A": [10.0, 10.1, 10.2], "B": [12.0, 12.1, 12.2]})

        restored = ProtocolReport.model_validate_json(report.model_dump_json())

        assert restored == report


class TestVerdict:
    """The verdict answers the question in one sentence."""

    def test_significant_difference_names_direction_and_evidence(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """A clear difference is reported as larger, with delta, CI, p and n."""
        _install_toy(monkeypatch)
        report = _run(tmp_path, {"A": [10.0, 10.1, 10.2], "B": [12.0, 12.1, 12.2]})

        assert len(report.verdict) == 1
        sentence = report.verdict[0]
        assert sentence.startswith(f"B {VERDICT_LARGER} mean_value than A")
        assert "delta +2" in sentence
        assert "95% CI" in sentence
        assert "p_adj" in sentence
        assert "n 3 vs 3" in sentence
        assert report.pairwise[0].significant is True

    def test_overlapping_conditions_report_no_significant_difference(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """Conditions that overlap get the no-difference sentence."""
        _install_toy(monkeypatch)
        report = _run(tmp_path, {"A": [10.0, 11.0, 12.0], "B": [10.2, 11.1, 11.9]})

        sentence = report.verdict[0]
        assert sentence.startswith(f"{VERDICT_NO_DIFFERENCE} in mean_value between A and B")
        assert "n 3 vs 3" in sentence
        assert report.pairwise[0].significant is False

    def test_single_condition_has_no_comparison(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """One config gives an empty pairwise list and a summary verdict."""
        _install_toy(monkeypatch)
        report = _run(tmp_path, {"A": [10.0, 10.1, 10.2]})

        assert report.pairwise == []
        assert len(report.verdict) == 1
        assert report.verdict[0].startswith("A mean_value 10.1 A (95% CI")
        assert "n 3" in report.verdict[0]

    def test_single_replicate_is_not_testable_rather_than_not_different(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """One replicate per condition makes the test undefined, and says so."""
        _install_toy(monkeypatch)
        report = _run(tmp_path, {"A": [10.0], "B": [12.0]})

        assert report.pairwise[0].testable is False
        assert report.pairwise[0].significant is False
        assert report.pairwise[0].delta_ci95 is None
        assert report.verdict[0].startswith(VERDICT_NOT_TESTABLE)
        assert any("one replicate" in warning for warning in report.warnings)


class TestAgentText:
    """The agent rendering stays inside its line budget."""

    def test_two_condition_report_fits_in_25_lines(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """A two-condition comparison renders compactly and carries the verdict."""
        _install_toy(monkeypatch)
        report = _run(tmp_path, {"A": [10.0, 10.1, 10.2], "B": [12.0, 12.1, 12.2]})

        text = report.to_agent_text()
        lines = text.strip().split("\n")

        assert len(lines) <= 25
        assert lines[0].startswith("# polyzymd analyze toy_protocol")
        assert "metric mean_value" in lines[0]
        assert "unit A" in lines[0]
        assert "eq 10ns" in lines[0]
        assert any(line.startswith("A  n 3") for line in lines)
        assert any(line.startswith("A vs B") for line in lines)
        assert any(line.startswith("verdict:") for line in lines)
        assert "|" not in text
        assert "" not in [line.strip() for line in lines]

    def test_many_conditions_still_fit_the_budget(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """A wide comparison drops lines and says how many it dropped."""
        _install_toy(monkeypatch)
        values = {f"C{index}": [10.0 + index, 10.1 + index, 10.2 + index] for index in range(12)}
        report = _run(tmp_path, values)

        lines = report.to_agent_text().strip().split("\n")

        assert len(lines) <= 25
        assert any("omitted" in line for line in lines)


class TestErrors:
    """Setup failures raise typed errors that say how to fix them."""

    def test_unknown_analysis_lists_the_known_ones(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """An unknown name raises ProtocolError naming the available plugins."""
        _install_toy(monkeypatch)

        with pytest.raises(ProtocolError) as excinfo:
            analyze("not_an_analysis", _write_configs(tmp_path, ["A"]))

        assert "Unknown analysis" in str(excinfo.value)
        assert "toy_protocol" in (excinfo.value.hint or "")

    def test_missing_config_names_the_path(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """A config path that does not exist raises before any computation."""
        _install_toy(monkeypatch)

        with pytest.raises(ProtocolError) as excinfo:
            analyze("toy_protocol", [tmp_path / "nope" / "config.yaml"], replicates=[1])

        assert "not found" in str(excinfo.value)
        assert excinfo.value.hint

    def test_no_configs_is_rejected(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Calling with an empty config list raises a typed error."""
        _install_toy(monkeypatch)

        with pytest.raises(ProtocolError) as excinfo:
            analyze("toy_protocol", [])

        assert "No simulation configs" in str(excinfo.value)

    def test_label_count_must_match_config_count(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """One label per config is required when labels are given at all."""
        _install_toy(monkeypatch)
        configs = _write_configs(tmp_path, ["A", "B"])

        with pytest.raises(ProtocolError) as excinfo:
            analyze("toy_protocol", configs, labels=["only_one"], replicates=[1])

        assert "label(s) for" in str(excinfo.value)


class TestArtifactShape:
    """The MDA comparison artifact shape is normalized the same way."""

    def test_condition_artifact_payload_is_read(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """A ComparisonArtifact payload yields the same report fields."""
        from polyzymd.analyses.mda.artifacts import ComparisonArtifact
        from polyzymd.config.comparison import ComparisonConfig

        _install_toy(monkeypatch)
        configs = _write_configs(tmp_path, ["A", "B"])
        config = ComparisonConfig(
            name="artifact_case",
            control="A",
            conditions=[
                {"label": "A", "config": configs[0], "replicates": [1, 2, 3]},
                {"label": "B", "config": configs[1], "replicates": [1, 2, 3]},
            ],
        )
        artifact = ComparisonArtifact(
            analysis_name="toy_protocol",
            conditions=["A", "B"],
            payload={
                "condition_summaries": [
                    {
                        "label": "A",
                        "n_replicates": 3,
                        "mean_rg_mean": 18.4,
                        "mean_rg_sem": 0.05,
                        "mean_rg_replicate_values": [18.4, 18.5, 18.3],
                        "mean_rg_unit": "A",
                        "mean_rg_ci95_low": 18.2,
                        "mean_rg_ci95_high": 18.6,
                        "mean_rg_ci_method": "student_t",
                    },
                    {
                        "label": "B",
                        "n_replicates": 3,
                        "mean_rg_mean": 18.71,
                        "mean_rg_sem": 0.06,
                        "mean_rg_replicate_values": [18.7, 18.8, 18.63],
                        "mean_rg_unit": "A",
                        "mean_rg_ci95_low": 18.5,
                        "mean_rg_ci95_high": 18.9,
                        "mean_rg_ci_method": "student_t",
                    },
                ],
                "pairwise_comparisons": [
                    {
                        "condition_a": "A",
                        "condition_b": "B",
                        "metric": "mean_rg",
                        "t_statistic": 5.2,
                        "p_value": 0.006,
                        "p_value_adjusted": 0.006,
                        "cohens_d": 4.2,
                        "effect_size_interpretation": "large",
                        "direction": "increased",
                        "significant": True,
                        "percent_change": 1.7,
                        "testable": True,
                    }
                ],
                "statistical_parameters": {
                    "ttest_method": "welch",
                    "posthoc_method": "ttest_bh",
                },
            },
        )
        pipeline_result = {
            "comparison": artifact,
            "aggregated": {},
            "comparison_path": tmp_path / "comparison.json",
            "plots": [],
        }

        report = build_report(ToyProtocolAnalysis(), config, pipeline_result)

        assert report.metric == "mean_rg"
        assert report.unit == "A"
        assert report.pairwise[0].test == "welch_t"
        assert report.pairwise[0].correction == "BH"
        assert report.pairwise[0].delta == pytest.approx(0.31, abs=1e-9)
        assert report.pairwise[0].delta_ci95 is not None
        assert report.pairwise[0].cohens_d == pytest.approx(-4.2)
        assert report.verdict[0].startswith("B larger mean_rg than A")

    @staticmethod
    def _observable_payload(*, significant: bool, p_adjusted: float) -> dict[str, Any]:
        """Contract comparison payload for two conditions and one observable."""

        def aggregate(mean: float, values: list[float]) -> dict[str, Any]:
            return {
                "name": "rg_protein",
                "kind": "mean_of_timeseries",
                "unit": "A",
                "n_replicates": 3,
                "replicate_values": values,
                "mean": mean,
                "sem": 0.05,
                "ci95_low": mean - 0.13,
                "ci95_high": mean + 0.13,
                "ci_method": "student_t",
            }

        return {
            "conditions": {
                "A": [
                    aggregate(18.4, [18.4, 18.5, 18.3]),
                    {
                        "name": "rg_protein_fragments",
                        "kind": "profile",
                        "unit": "A",
                        "n_replicates": 3,
                        "profile_mean": [1.0, 2.0],
                        "index": [0.0, 1.0],
                    },
                ],
                "B": [aggregate(18.71, [18.7, 18.8, 18.63])],
            },
            "comparisons": [
                {
                    "name": "rg_protein",
                    "kind": "mean_of_timeseries",
                    "unit": "A",
                    "control": "A",
                    "condition": "B",
                    "n_control": 3,
                    "n_condition": 3,
                    "delta": 0.31,
                    "percent_change": 1.7,
                    "test": "welch_t",
                    "p_value": 0.006,
                    "p_adjusted": p_adjusted,
                    "correction": "benjamini_hochberg",
                    "cohens_d": 4.2,
                    "significant": significant,
                    "testable": True,
                }
            ],
        }

    def _report(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
        payload: dict[str, Any],
    ) -> Any:
        """Build a report from one comparison payload, whichever shape it has."""
        from polyzymd.analyses.mda.artifacts import ComparisonArtifact
        from polyzymd.config.comparison import ComparisonConfig

        _install_toy(monkeypatch)
        configs = _write_configs(tmp_path, ["A", "B"])
        config = ComparisonConfig(
            name="payload_case",
            control="A",
            conditions=[
                {"label": "A", "config": configs[0], "replicates": [1, 2, 3]},
                {"label": "B", "config": configs[1], "replicates": [1, 2, 3]},
            ],
        )
        artifact = ComparisonArtifact(
            analysis_name="toy_protocol", conditions=["A", "B"], payload=payload
        )
        return build_report(
            ToyProtocolAnalysis(),
            config,
            {
                "comparison": artifact,
                "aggregated": {},
                "comparison_path": tmp_path / "comparison.json",
                "plots": [],
            },
        )

    def test_observable_contract_payload_is_read(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """A contract plugin's comparison payload reports its observables."""
        report = self._report(
            monkeypatch,
            tmp_path,
            self._observable_payload(significant=True, p_adjusted=0.006),
        )

        assert report.run == "rg_protein"
        assert report.all_runs == ["rg_protein"]
        assert report.unit == "A"
        assert [condition.label for condition in report.conditions] == ["A", "B"]
        assert report.pairwise[0].test == "welch_t"
        assert report.pairwise[0].delta == pytest.approx(0.31)
        assert report.pairwise[0].direction == "increased"
        assert report.pairwise[0].significant is True

    def test_observable_contract_claims_no_direction_without_significance(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """The same difference with a p value above alpha names no direction."""
        from polyzymd.analyses.shared.inferential_statistics import NO_SIGNIFICANT_CHANGE

        report = self._report(
            monkeypatch,
            tmp_path,
            self._observable_payload(significant=False, p_adjusted=0.42),
        )

        assert report.pairwise[0].delta == pytest.approx(0.31)
        assert report.pairwise[0].significant is False
        assert report.pairwise[0].direction == NO_SIGNIFICANT_CHANGE

    def test_both_payload_shapes_use_one_direction_vocabulary(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """The same rise reads the same whether a plugin or the framework tested it.

        ``direction`` is a public field of the report, so a plugin that still
        owns its comparison and one on the observable contract must not offer
        an agent two words for the same finding.
        """
        legacy = self._report(monkeypatch, tmp_path / "legacy", _LEGACY_PAYLOAD)
        contract = self._report(
            monkeypatch,
            tmp_path / "contract",
            self._observable_payload(significant=True, p_adjusted=0.006),
        )

        assert legacy.pairwise[0].direction == contract.pairwise[0].direction


_LEGACY_PAYLOAD: dict[str, Any] = {
    "condition_summaries": [
        {
            "label": "A",
            "n_replicates": 3,
            "mean_rg_mean": 18.4,
            "mean_rg_sem": 0.05,
            "mean_rg_replicate_values": [18.4, 18.5, 18.3],
            "mean_rg_unit": "A",
            "mean_rg_ci95_low": 18.27,
            "mean_rg_ci95_high": 18.53,
            "mean_rg_ci_method": "student_t",
        },
        {
            "label": "B",
            "n_replicates": 3,
            "mean_rg_mean": 18.71,
            "mean_rg_sem": 0.05,
            "mean_rg_replicate_values": [18.7, 18.8, 18.63],
            "mean_rg_unit": "A",
            "mean_rg_ci95_low": 18.58,
            "mean_rg_ci95_high": 18.84,
            "mean_rg_ci_method": "student_t",
        },
    ],
    "pairwise_comparisons": [
        {
            "condition_a": "A",
            "condition_b": "B",
            "metric": "mean_rg",
            "p_value": 0.006,
            "p_value_adjusted": 0.006,
            "cohens_d": 4.2,
            "direction": interpret_direction(1.7),
            "significant": True,
            "percent_change": 1.7,
            "testable": True,
        }
    ],
    "statistical_parameters": {"ttest_method": "welch", "posthoc_method": "ttest_bh"},
}
