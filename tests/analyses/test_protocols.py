"""Tests for the agent-facing analysis protocol: report fields, verdicts, text and errors.

The reports come from :class:`~polyzymd.analyses.timeseries.ReplicateValues`
built from given replicate values, the path every analysis of
``polyzymd analyze`` reports through.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Sequence

import pytest

from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.protocols import (
    FUNCTION_ANALYSES,
    VERDICT_LARGER,
    VERDICT_NO_DIFFERENCE,
    VERDICT_NOT_TESTABLE,
    ConditionReport,
    PairwiseReport,
    ProtocolReport,
    _difference_ci,
    analyze,
)


class _Study:
    """The parts of a Study that ReplicateValues reads."""

    def __init__(self, labels: Sequence[str]) -> None:
        self.labels, self.control = list(labels), labels[0]

    def __getitem__(self, label: str) -> SimpleNamespace:
        return SimpleNamespace(equilibration="10ns", stride=1, config_hash=f"hash-{label}")


def _report(tmp_path: Path, values: dict[str, list[float]]) -> ProtocolReport:
    """Compare the given per-condition replicate values against the first condition.

    Parameters
    ----------
    tmp_path : Path
        Folder named as the results path of the report.
    values : dict
        Replicate values keyed by condition label, control first.

    Returns
    -------
    ProtocolReport
        The summary for one condition, the comparison for several.
    """
    from polyzymd.analyses.timeseries import ReplicateValues, Source

    rows = {
        label: [
            (index + 1, value, None, None, 100, None, None) for index, value in enumerate(items)
        ]
        for label, items in values.items()
    }
    source = Source("rg", "A", _Study(list(values)), {label: [] for label in values}, tmp_path)
    result = ReplicateValues(source, "mean_rg", "A", False, rows)
    return result.compare() if len(values) > 1 else result.summary()


class TestReportFields:
    """The report states what every number is."""

    def test_two_conditions_report_is_fully_typed(self, tmp_path: Path) -> None:
        """Every documented field is present and carries the declared type."""
        report = _report(tmp_path, {"A": [10.0, 10.1, 10.2], "B": [12.0, 12.1, 12.2]})

        assert isinstance(report, ProtocolReport)
        assert report.analysis == "rg"
        assert report.metric == "mean_rg"
        assert report.unit == "A"
        assert report.equilibration == "10ns"
        assert report.frames_per_replicate == {"A": [100, 100, 100], "B": [100, 100, 100]}

        assert [condition.label for condition in report.conditions] == ["A", "B"]
        first = report.conditions[0]
        assert isinstance(first, ConditionReport)
        assert first.n_replicates == 3
        assert first.mean == pytest.approx(10.1)
        assert first.sem is not None and first.sem > 0.0
        assert first.ci95 is not None and first.ci95[0] < first.mean < first.ci95[1]
        assert first.ci_method == "student_t"
        assert first.replicate_values == [10.0, 10.1, 10.2]
        assert first.replicates == [1, 2, 3]

        assert len(report.pairwise) == 1
        pair = report.pairwise[0]
        assert isinstance(pair, PairwiseReport)
        assert (pair.a, pair.b) == ("A", "B")
        assert pair.delta == pytest.approx(2.0)
        assert pair.delta_ci95 is not None
        assert pair.p is not None and pair.p_adjusted is not None
        assert pair.test == "welch_t"
        assert pair.correction == "BH"
        assert pair.family_size == 1
        # Effect sizes are oriented like delta, so both are positive when the
        # second condition is larger.
        assert pair.cohens_d is not None and pair.cohens_d > 0.0
        assert pair.hedges_g is not None and 0.0 < pair.hedges_g < pair.cohens_d
        assert pair.testable is True

        assert report.provenance.polyzymd_version
        assert report.provenance.config_hashes == {"A": "hash-A", "B": "hash-B"}
        assert report.provenance.output_paths == {"results": str(tmp_path)}
        assert report.verdict

    def test_json_round_trips_through_the_model(self, tmp_path: Path) -> None:
        """The JSON form validates back into an equal report."""
        report = _report(tmp_path, {"A": [10.0, 10.1, 10.2], "B": [12.0, 12.1, 12.2]})

        assert ProtocolReport.model_validate_json(report.model_dump_json()) == report


class TestVerdict:
    """The verdict answers the question in one sentence."""

    def test_significant_difference_names_direction_and_evidence(self, tmp_path: Path) -> None:
        """A clear difference is reported as larger, with delta, CI, p and n."""
        report = _report(tmp_path, {"A": [10.0, 10.1, 10.2], "B": [12.0, 12.1, 12.2]})

        assert len(report.verdict) == 1
        sentence = report.verdict[0]
        assert sentence.startswith(f"B {VERDICT_LARGER} mean_rg than A")
        assert "delta +2" in sentence
        assert "95% CI" in sentence
        assert "p_adj" in sentence
        assert "n 3 vs 3" in sentence
        assert report.pairwise[0].significant is True

    def test_overlapping_conditions_report_no_significant_difference(self, tmp_path: Path) -> None:
        """Conditions that overlap get the no-difference sentence."""
        report = _report(tmp_path, {"A": [10.0, 11.0, 12.0], "B": [10.2, 11.1, 11.9]})

        sentence = report.verdict[0]
        assert sentence.startswith(f"{VERDICT_NO_DIFFERENCE} in mean_rg between A and B")
        assert "n 3 vs 3" in sentence
        assert report.pairwise[0].significant is False

    def test_single_condition_has_no_comparison(self, tmp_path: Path) -> None:
        """One condition gives an empty pairwise list and a summary verdict."""
        report = _report(tmp_path, {"A": [10.0, 10.1, 10.2]})

        assert report.pairwise == []
        assert len(report.verdict) == 1
        assert report.verdict[0].startswith("A mean_rg 10.1 A (95% CI")
        assert "n 3" in report.verdict[0]

    def test_single_replicate_is_not_testable_rather_than_not_different(
        self, tmp_path: Path
    ) -> None:
        """One replicate per condition makes the test undefined, and says so."""
        report = _report(tmp_path, {"A": [10.0], "B": [12.0]})

        assert report.pairwise[0].testable is False
        assert report.pairwise[0].significant is False
        assert report.pairwise[0].delta_ci95 is None
        assert report.verdict[0].startswith(VERDICT_NOT_TESTABLE)
        assert any("one replicate" in warning for warning in report.warnings)


class TestAgentText:
    """The agent rendering prints every item on its own line."""

    def test_two_condition_report_fits_in_25_lines(self, tmp_path: Path) -> None:
        """A two-condition comparison renders compactly and carries the verdict."""
        report = _report(tmp_path, {"A": [10.0, 10.1, 10.2], "B": [12.0, 12.1, 12.2]})

        text = report.to_agent_text()
        lines = text.strip().split("\n")

        assert len(lines) <= 25
        assert lines[0].startswith("# polyzymd analyze rg")
        assert "metric mean_rg" in lines[0]
        assert "unit A" in lines[0]
        assert "eq 10ns" in lines[0]
        assert any(line.startswith("A  n 3") for line in lines)
        assert any(line.startswith("A vs B") for line in lines)
        assert any(line.startswith("verdict:") for line in lines)
        assert "|" not in text
        assert "" not in [line.strip() for line in lines]

    def test_many_conditions_print_every_line(self, tmp_path: Path) -> None:
        """A wide comparison prints every condition and comparison."""
        values = {f"C{index}": [10.0 + index, 10.1 + index, 10.2 + index] for index in range(12)}
        lines = _report(tmp_path, values).to_agent_text().strip().split("\n")

        conditions = [line for line in lines if line.startswith("C") and "  n 3  " in line]
        comparisons = [line for line in lines if line.startswith("C0 vs ")]
        assert len(conditions) == 12
        assert len(comparisons) == 11
        assert not any("omitted" in line for line in lines)

    def test_every_replicate_value_is_printed(self, tmp_path: Path) -> None:
        """A condition with many replicates prints all of their values."""
        values = {
            "A": [10.0 + 0.1 * index for index in range(8)],
            "B": [12.0 + 0.1 * index for index in range(8)],
        }
        text = _report(tmp_path, values).to_agent_text()

        line = next(line for line in text.split("\n") if line.startswith("A  n 8"))
        printed = line.split("  values ", 1)[1].split(", ")
        assert len(printed) == 8
        assert "more" not in line


class TestDifferenceInterval:
    """The interval on a difference matches the reported test."""

    A = [1.0, 2.0, 3.0]
    B = [2.0, 4.0, 6.0, 8.0]

    def test_student_uses_the_pooled_variance(self) -> None:
        """Pooled variance 4.4 with 5 degrees of freedom gives 3 +/- 4.118."""
        low, high = _difference_ci(self.A, self.B, "student_t")
        assert (low, high) == pytest.approx((-1.118283, 7.118283), abs=1e-6)

    def test_welch_uses_separate_variances(self) -> None:
        """Welch's interval is narrower here than the pooled one."""
        low, high = _difference_ci(self.A, self.B, "welch_t")
        assert (low, high) == pytest.approx((-0.897975, 6.897975), abs=1e-6)

    @pytest.mark.parametrize(
        "values_a, values_b, test",
        [
            (A, B, "tukey_hsd"),
            ([1.0], B, "student_t"),
            ([1.0, 1.0], [2.0, 2.0], "student_t"),
            ([1.0, 1.0], [2.0, 2.0], "welch_t"),
        ],
    )
    def test_no_interval_when_none_can_be_estimated(
        self, values_a: list[float], values_b: list[float], test: str
    ) -> None:
        """Tukey, one replicate and zero variance give no interval."""
        assert _difference_ci(values_a, values_b, test) is None


class TestErrors:
    """Setup failures raise typed errors that say how to fix them."""

    def test_unknown_analysis_lists_the_analyses_and_the_function_api(self, tmp_path: Path) -> None:
        """An unknown name names every analysis and the page on writing a function."""
        with pytest.raises(ProtocolError) as excinfo:
            analyze("not_an_analysis", [tmp_path / "A" / "config.yaml"])

        assert str(excinfo.value) == "No analysis named 'not_an_analysis'."
        hint = excinfo.value.hint or ""
        assert hint.startswith(f"Use one of {', '.join(FUNCTION_ANALYSES)}.")
        assert "Study.timeseries or Study.per_replicate" in hint
        assert "how_to/study_api.html" in hint

    @pytest.mark.parametrize("name", ["toy_protocol", "radius_of_gyration", "catalytic_triad_v1"])
    def test_names_outside_the_analyses_are_refused_before_any_config_is_read(
        self, name: str, tmp_path: Path
    ) -> None:
        """No name outside FUNCTION_ANALYSES reaches a config, even a missing one."""
        with pytest.raises(ProtocolError, match="No analysis named"):
            analyze(name, [tmp_path / "missing" / "config.yaml"], stride=5)

    def test_catalytic_triad_keeps_its_refusal(self, tmp_path: Path) -> None:
        """The retired triad analysis points at the routine on the study API."""
        with pytest.raises(ProtocolError, match="no longer a polyzymd analyze analysis") as info:
            analyze("catalytic_triad", [tmp_path / "A" / "config.yaml"])

        assert "how_to/analysis_triad_quickstart.html" in (info.value.hint or "")

    def test_missing_config_names_the_path(self, tmp_path: Path) -> None:
        """A config path that does not exist raises before any computation."""
        missing = tmp_path / "nope" / "config.yaml"
        with pytest.raises(ProtocolError) as excinfo:
            analyze("rg", [missing], replicates=[1])

        assert "not found" in str(excinfo.value)
        assert str(missing) in str(excinfo.value)
        assert excinfo.value.hint

    def test_no_configs_is_rejected(self) -> None:
        """Calling with an empty config list raises a typed error."""
        with pytest.raises(ProtocolError, match="none given"):
            analyze("rg", [])

    def test_label_count_must_match_config_count(self, tmp_path: Path) -> None:
        """One label per config is required when labels are given at all."""
        configs = []
        for label in ("A", "B"):
            (tmp_path / label).mkdir()
            configs.append(tmp_path / label / "config.yaml")
            configs[-1].write_text("placeholder: true\n")

        with pytest.raises(ProtocolError, match="label\\(s\\) for"):
            analyze("rg", configs, labels=["only_one"], replicates=[1])


def test_repeated_pair_labels_are_refused() -> None:
    """Two distance pairs with one label would overwrite each other's values."""
    from polyzymd.analyses.protocols import _analyze_pairs

    pair = {"label": "d", "selection_a": "name C1", "selection_b": "name C2"}
    with pytest.raises(ProtocolError, match="repeat"):
        _analyze_pairs(
            "distances",
            None,
            {"pairs": [pair, pair]},
            None,
            recompute=False,
            output_dir=None,
            eq_check=False,
            plots=False,
        )


@pytest.mark.parametrize(
    ("limits", "text"),
    [((1.9999995, 2.0000004), "1.999999 to 2"), ((0.7733, 1.227), "0.7733 to 1.227"), (None, "na")],
)
def test_a_narrow_interval_is_not_printed_as_one_number(limits, text) -> None:
    """A real, narrow interval prints with enough digits to show both ends, not as 'ci95 2 to 2'."""
    from polyzymd.analyses.protocols import _interval

    assert _interval(limits) == text


def test_the_untestable_reason_is_the_real_one() -> None:
    """Two replicates with no variance are said to have no variance, not to lack replicates."""
    from polyzymd.analyses.protocols import _verdict

    conditions = [
        ConditionReport(label=label, n_replicates=2, mean=0.0, replicate_values=[0.0, 0.0])
        for label in ("A", "B")
    ]
    pair = PairwiseReport(a="A", b="B", delta=0.0, testable=False)
    text = " ".join(_verdict("m", None, conditions, [pair]))
    assert "same value in every replicate" in text and "at least two" not in text


def test_a_selection_on_a_missing_topology_attribute_is_a_protocol_error() -> None:
    """'chainid A' on a topology without chain IDs names the selection, not a bare AttributeError."""
    mda = pytest.importorskip("MDAnalysis")
    from polyzymd.analyses.protocols import _empty_selections

    universe = mda.Universe.empty(3, trajectory=True)
    replicate = SimpleNamespace(index=1, universe=lambda: universe)
    study = [SimpleNamespace(label="Water", replicates=[replicate])]
    with pytest.raises(ProtocolError, match="'chainid A'.*chainIDs") as info:
        _empty_selections(study, {"protein_selection": "chainid A"})
    assert "resname" in (info.value.hint or "")


def test_a_test_with_two_replicates_and_a_constant_control_warns_of_low_power() -> None:
    """A significant verdict at n 2 vs 2 against a control with no variance says it has little power."""
    from polyzymd.analyses.protocols import _verdict

    conditions = [
        ConditionReport(label="Water", n_replicates=2, mean=0.0, replicate_values=[0.0, 0.0]),
        ConditionReport(label="SDS", n_replicates=2, mean=0.4, replicate_values=[0.39, 0.41]),
    ]
    pair = PairwiseReport(
        a="Water", b="SDS", delta=0.4, p=0.03, p_adjusted=0.03, significant=True, testable=True
    )
    (text,) = _verdict("coverage", None, conditions, [pair])
    assert text.startswith("SDS larger coverage than Water")
    assert "little power: Water and SDS have fewer than 3 replicates" in text
    assert "Water has the same value in every replicate" in text


def test_a_test_with_three_varying_replicates_has_no_power_warning() -> None:
    from polyzymd.analyses.protocols import _verdict

    conditions = [
        ConditionReport(label=label, n_replicates=3, mean=1.0, replicate_values=[0.9, 1.0, 1.1])
        for label in ("A", "B")
    ]
    pair = PairwiseReport(a="A", b="B", delta=0.0, p=0.9, p_adjusted=0.9, testable=True)
    assert "power" not in " ".join(_verdict("m", None, conditions, [pair]))


def test_a_single_replicate_is_not_said_to_have_the_same_value_in_every_replicate() -> None:
    """At n 1 the power note names only the replicate count."""
    from polyzymd.analyses.protocols import _verdict

    conditions = [
        ConditionReport(label="Water", n_replicates=1, mean=0.0, replicate_values=[0.0]),
        ConditionReport(label="SDS", n_replicates=3, mean=0.4, replicate_values=[0.3, 0.4, 0.5]),
    ]
    pair = PairwiseReport(a="Water", b="SDS", delta=0.4, p=0.03, p_adjusted=0.03, testable=True)
    (text,) = _verdict("coverage", None, conditions, [pair])
    assert "little power: Water has fewer than 3 replicates" in text
    assert "same value in every replicate" not in text
