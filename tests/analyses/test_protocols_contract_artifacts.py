"""Build a report from the comparison artifact every plugin now writes.

The nine stored campaign artifacts this file replaces covered three shapes a
plugin could write: the MDAnalysis comparison artifact, the framework scalar
result, and a plugin's own result grouped by run or pair label. There is one
shape left, so the cases here are generated rather than stored: each one runs
the real ``Analysis.compare()`` over condition aggregates built in memory and
checks that ``build_report`` states the metric, the interval, the test and the
verdict it should.

The numbers are chosen so the expected answer is arithmetic a reader can check:
the control replicates are 1.0, 1.1 and 1.2, the treated replicates are the same
values plus a fixed offset, so the difference in means is exactly that offset.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any, Sequence

import pytest

from polyzymd.analyses._framework.contexts import ComparisonContext, Condition
from polyzymd.analyses.contract import ObservableEstimate, aggregate_observables
from polyzymd.analyses.discovery import get_analysis, list_all_names
from polyzymd.analyses.mda.artifacts import ComparisonArtifact, ConditionArtifact
from polyzymd.analyses.protocols import ProtocolReport, build_report
from tests.analyses.conftest import make_simulation_config

CONTROL = "No Polymer (Control)"
TREATED = "SBMA-EGMA 0:100"
CONTROL_VALUES = (1.0, 1.1, 1.2)
OFFSET = 0.5


def _aggregate(
    analysis_name: str,
    label: str,
    values: Sequence[float],
    *,
    observable: str,
    kind: str = "mean_of_timeseries",
    unit: str | None = "A",
) -> ConditionArtifact:
    """Build one condition aggregate holding a single observable."""
    replicates = [
        [
            ObservableEstimate(
                name=observable, kind=kind, unit=unit, value=float(value), n_frames=100
            )
        ]
        for value in values
    ]
    aggregates = aggregate_observables(replicates)
    return ConditionArtifact(
        analysis_name=analysis_name,
        condition_label=label,
        replicates=list(range(1, len(values) + 1)),
        payload={"observables": [item.model_dump(mode="json") for item in aggregates]},
        metadata={"n_replicates": len(values)},
    )


def _comparison(
    analysis_name: str, tmp_path: Path, *, observable: str, kind: str, unit: str | None
) -> ComparisonArtifact:
    """Run the real comparison stage over two conditions of one plugin."""
    analysis = get_analysis(analysis_name)()
    treated = [value + OFFSET for value in CONTROL_VALUES]
    aggregates = {
        CONTROL: _aggregate(
            analysis_name, CONTROL, CONTROL_VALUES, observable=observable, kind=kind, unit=unit
        ),
        TREATED: _aggregate(
            analysis_name, TREATED, treated, observable=observable, kind=kind, unit=unit
        ),
    }
    conditions = [
        Condition(
            label=label,
            config_path=tmp_path / f"{label}.yaml",
            replicates=(1, 2, 3),
            sim_config=make_simulation_config(label),
        )
        for label in aggregates
    ]
    comparison = analysis.compare(
        ComparisonContext(
            name="campaign",
            conditions=conditions,
            excluded_conditions=[],
            control_label=CONTROL,
            analysis_dirs={
                label: tmp_path / "analysis" / label / analysis_name for label in aggregates
            },
            results_dir=tmp_path / "comparison" / analysis_name,
            equilibration="10ns",
            settings=None,
            aggregated_results=aggregates,
        )
    )
    assert isinstance(comparison, ComparisonArtifact)
    return comparison


def _report(analysis_name: str, comparison: ComparisonArtifact, tmp_path: Path) -> ProtocolReport:
    """Normalize one comparison artifact into a report."""
    config = SimpleNamespace(
        name="campaign",
        source_path=tmp_path / "comparison.yaml",
        defaults=SimpleNamespace(equilibration_time="10ns"),
        control=CONTROL,
        conditions=[
            SimpleNamespace(label=label, config=tmp_path / f"{label}.yaml", replicates=[1, 2, 3])
            for label in (CONTROL, TREATED)
        ],
    )
    return build_report(
        get_analysis(analysis_name)(),
        config,
        {"comparison": comparison, "aggregated": {}, "plots": []},
    )


@pytest.mark.parametrize("analysis_name", sorted(list_all_names()))
def test_every_plugin_artifact_normalizes_into_a_report(analysis_name: str, tmp_path: Path) -> None:
    """Each plugin's comparison artifact gives a fully typed report."""
    comparison = _comparison(
        analysis_name,
        tmp_path,
        observable=f"{analysis_name}_value",
        kind="mean_of_timeseries",
        unit="A",
    )

    report = _report(analysis_name, comparison, tmp_path)

    assert report.analysis == analysis_name
    assert report.metric == f"{analysis_name}_value"
    assert report.unit == "A"
    assert [condition.label for condition in report.conditions] == [CONTROL, TREATED]
    assert report.conditions[0].mean == pytest.approx(1.1)
    assert report.conditions[0].ci95 is not None
    assert report.conditions[0].ci_method == "student_t"
    assert len(report.pairwise) == 1
    pair = report.pairwise[0]
    assert (pair.a, pair.b) == (CONTROL, TREATED)
    assert pair.delta == pytest.approx(OFFSET)
    assert pair.test == "student_t"
    assert pair.correction == "benjamini_hochberg"
    assert pair.testable is True
    assert pair.hedges_g is not None


@pytest.mark.parametrize(
    ("kind", "unit"),
    [("mean_of_timeseries", "A"), ("fluctuation", "A"), ("fraction", "fraction")],
)
def test_every_scalar_kind_reaches_the_report(kind: str, unit: str, tmp_path: Path) -> None:
    """A report states the unit of whichever scalar kind the plugin reported."""
    comparison = _comparison("rg", tmp_path, observable="probe", kind=kind, unit=unit)

    report = _report("rg", comparison, tmp_path)

    assert report.metric == "probe"
    assert report.unit == unit
    assert report.verdict


def test_a_second_observable_is_listed_rather_than_dropped(tmp_path: Path) -> None:
    """The report names the group it renders and lists the others."""
    analysis = get_analysis("rg")()
    aggregates: dict[str, Any] = {}
    for label, offset in ((CONTROL, 0.0), (TREATED, OFFSET)):
        replicates = [
            [
                ObservableEstimate(
                    name=name,
                    kind="mean_of_timeseries",
                    unit="A",
                    value=float(value + offset) * scale,
                    n_frames=100,
                )
                for name, scale in (("first", 1.0), ("second", 2.0))
            ]
            for value in CONTROL_VALUES
        ]
        aggregates[label] = ConditionArtifact(
            analysis_name="rg",
            condition_label=label,
            replicates=[1, 2, 3],
            payload={
                "observables": [
                    item.model_dump(mode="json") for item in aggregate_observables(replicates)
                ]
            },
        )
    conditions = [
        Condition(
            label=label,
            config_path=tmp_path / f"{label}.yaml",
            replicates=(1, 2, 3),
            sim_config=make_simulation_config(label),
        )
        for label in aggregates
    ]
    comparison = analysis.compare(
        ComparisonContext(
            name="campaign",
            conditions=conditions,
            excluded_conditions=[],
            control_label=CONTROL,
            analysis_dirs={label: tmp_path / "analysis" / label / "rg" for label in aggregates},
            results_dir=tmp_path / "comparison" / "rg",
            equilibration="10ns",
            settings=None,
            aggregated_results=aggregates,
        )
    )

    report = _report("rg", comparison, tmp_path)

    assert report.run == "first"
    assert report.all_runs == ["first", "second"]
    assert len(report.pairwise) == 1


def test_an_artifact_from_before_the_contract_is_a_typed_error(tmp_path: Path) -> None:
    """A campaign tree still holding a pre-contract result says so."""
    from polyzymd.analyses.exceptions import ProtocolError

    stale = ComparisonArtifact(
        analysis_name="rg",
        conditions=[CONTROL, TREATED],
        control_label=CONTROL,
        payload={"condition_summaries": [{"label": CONTROL, "mean_rg_mean": 1.1}]},
    )

    with pytest.raises(ProtocolError, match="observable-contract artifact"):
        _report("rg", stale, tmp_path)


def test_the_report_json_round_trips(tmp_path: Path) -> None:
    """The JSON form validates back into an equal report."""
    comparison = _comparison(
        "rmsd", tmp_path, observable="mean_rmsd", kind="mean_of_timeseries", unit="A"
    )

    report = _report("rmsd", comparison, tmp_path)

    assert ProtocolReport.model_validate_json(report.model_dump_json()) == report
