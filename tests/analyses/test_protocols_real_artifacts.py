"""Build reports from the nine stored comparison artifacts of a real campaign.

The files under ``tests/data/comparison_artifacts/`` are the comparison results
of a LipA campaign, trimmed to the fields the report builder reads. They cover
all three shapes a plugin can store: the MDAnalysis comparison artifact
(catalytic_triad, hydrogen_bonds, rmsf, secondary_structure), a custom result
grouped by run or pair label (rg, rmsd, sasa, distances), and a custom result
that nests its statistics one level deeper (contacts).

The catalytic_triad and distances files were written by the implementations
those two plugins had before they moved to the observable contract, and they
are kept exactly as they were: the report builder still has to read artifacts a
user already has on disk. Two things in them no longer match a fresh run. The
triad file states its simultaneous contact fraction as a percentage, from 0.0
to 26.65, because the old code multiplied by 100; a new run stores the same
quantity as a fraction with unit "fraction". The distances file groups by pair
label under a custom result, where a new run stores one observable per pair.

Each case pins one hand-checked pair against the numbers in its own file, so a
regression in the normalizer shows up as a wrong verdict rather than as a
silently different number.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from polyzymd.analyses.mda.artifacts import ComparisonArtifact
from polyzymd.analyses.protocols import ProtocolReport, build_report

ARTIFACT_DIR = Path(__file__).resolve().parents[1] / "data" / "comparison_artifacts"

CONTROL = "No Polymer (Control)"
TREATED = "SBMA-EGMA 0:100"

# analysis, expected metric, expected run, then the hand-checked pair: the two
# condition labels, the adjusted p value in the file, whether the file calls it
# significant, the replicate count per arm, and the sign of mean(b) - mean(a).
CASES = [
    (
        "catalytic_triad",
        "simultaneous_contact_fraction",
        None,
        CONTROL,
        TREATED,
        0.3465935070873343,
        False,
        5,
        +1,
    ),
    ("contacts", "coverage", None, TREATED, "SBMA-EGMA 25:75", 0.0007176137549934729, True, 5, +1),
    (
        "distances",
        "mean_distance",
        "Ile12(N)-Substrate(carbonyl C)",
        CONTROL,
        TREATED,
        None,
        False,
        5,
        -1,
    ),
    (
        "hydrogen_bonds",
        "mean_hbonds_protein_polymer",
        None,
        CONTROL,
        TREATED,
        1.2333522411286452e-07,
        True,
        5,
        +1,
    ),
    ("rg", "mean_rg", "Protein", CONTROL, TREATED, 0.9345684211892554, False, 5, +1),
    ("rmsd", "mean_rmsd", "Whole Protein CA", CONTROL, TREATED, 0.06274019478822104, False, 5, -1),
    ("rmsf", "mean_rmsf", None, CONTROL, TREATED, 0.01060537998278552, True, 5, -1),
    ("sasa", "mean_sasa", "protein_isolated", CONTROL, TREATED, 0.5308012869205083, False, 5, +1),
    (
        "secondary_structure",
        "helix_fraction",
        None,
        CONTROL,
        TREATED,
        0.008399959791450017,
        True,
        5,
        +1,
    ),
]


class _StoredResult:
    """A plugin comparison result loaded from its saved JSON."""

    def __init__(self, data: dict[str, Any]) -> None:
        self._data = data
        self.warnings = list(data.get("warnings") or [])

    def model_dump(self) -> dict[str, Any]:
        return self._data


def _load(analysis: str) -> Any:
    """Load one stored comparison as the object the pipeline would have returned."""
    data = json.loads((ARTIFACT_DIR / f"{analysis}.json").read_text())
    if data.get("artifact_type") == "comparison":
        return ComparisonArtifact.model_validate(data)
    return _StoredResult(data)


def _labels(comparison: Any) -> list[str]:
    """List the condition labels a stored comparison covers."""
    if isinstance(comparison, ComparisonArtifact):
        return [
            str(item.get("label", "")) for item in comparison.payload.get("condition_summaries", [])
        ]
    return [str(item.get("label", "")) for item in comparison.model_dump().get("conditions", [])]


def _report(analysis: str, run: str | None = None) -> ProtocolReport:
    """Build a report from one stored comparison without touching a trajectory."""
    comparison = _load(analysis)
    config = SimpleNamespace(
        defaults=SimpleNamespace(equilibration_time="200ns"),
        conditions=[
            SimpleNamespace(label=label, config=Path("/nonexistent/config.yaml"))
            for label in _labels(comparison)
        ],
        plugins=SimpleNamespace(get=lambda name: None),
    )
    plugin = SimpleNamespace(
        name=analysis,
        protocol_version="1",
        Settings=type("Settings", (), {}),
        aggregate_settings_fingerprint=lambda settings: None,
    )
    return build_report(
        plugin,
        config,
        {"comparison": comparison, "aggregated": {}, "plots": []},
        run=run,
    )


@pytest.mark.parametrize(
    "analysis,metric,run,label_a,label_b,p_adjusted,significant,n,sign",
    CASES,
    ids=[case[0] for case in CASES],
)
def test_stored_comparison_reports_its_own_numbers(
    analysis: str,
    metric: str,
    run: str | None,
    label_a: str,
    label_b: str,
    p_adjusted: float | None,
    significant: bool,
    n: int,
    sign: int,
) -> None:
    """Every stored comparison builds a report whose hand-checked pair matches the file."""
    report = _report(analysis)

    assert report.analysis == analysis
    assert report.metric == metric
    assert report.run == run
    assert report.conditions, "no condition summaries were recovered"
    assert all(condition.n_replicates == n for condition in report.conditions)
    assert report.pairwise, "no comparisons were recovered"

    pair = next(p for p in report.pairwise if (p.a, p.b) == (label_a, label_b))
    if p_adjusted is None:
        assert pair.p_adjusted is None
    else:
        assert pair.p_adjusted == pytest.approx(p_adjusted)
    assert pair.significant is significant
    assert pair.delta * sign > 0.0

    verdict = next(text for text in report.verdict if label_a in text and label_b in text)
    if p_adjusted is None:
        assert verdict.startswith("no test recorded")
    elif significant:
        assert verdict.startswith(f"{label_b} {'larger' if sign > 0 else 'smaller'} {metric}")
    else:
        assert verdict.startswith("no significant difference")
    assert f"n {n} vs {n}" in verdict


@pytest.mark.parametrize("analysis", [case[0] for case in CASES])
def test_stored_comparison_renders_agent_text(analysis: str) -> None:
    """Agent text stays inside its budget and carries a verdict for every stored result."""
    text = _report(analysis).to_agent_text()
    lines = text.strip().split("\n")

    assert len(lines) <= 25
    assert lines[0].startswith(f"# polyzymd analyze {analysis}")
    assert any(line.startswith("verdict:") for line in lines)
    assert "" not in [line.strip() for line in lines]


@pytest.mark.parametrize("analysis", [case[0] for case in CASES])
def test_stored_comparison_round_trips_as_json(analysis: str) -> None:
    """Every report validates back from its own JSON."""
    report = _report(analysis)
    assert ProtocolReport.model_validate_json(report.model_dump_json()) == report


def test_contacts_nested_statistics_are_not_reported_as_a_null_result() -> None:
    """Contacts nests its statistics one level deeper; they must still be read.

    The stored file records p_adjusted 7.176e-4 and Cohen's d -3.863 for the
    coverage of SBMA-EGMA 0:100 against SBMA-EGMA 25:75. Reading only the outer
    rows loses every statistic and turns a real difference into a null result.
    """
    report = _report("contacts")
    pair = next(
        item
        for item in report.pairwise
        if (item.a, item.b) == ("SBMA-EGMA 0:100", "SBMA-EGMA 25:75")
    )

    assert pair.p_adjusted == pytest.approx(0.0007176137549934729)
    assert pair.significant is True
    # The report orients Cohen's d like delta, so the file's -3.863 is reported
    # as +3.863 next to a positive difference.
    assert pair.cohens_d == pytest.approx(3.86281027985395)
    assert pair.delta > 0.0
    assert not any("no significant difference" in text for text in report.verdict[:1])


def test_run_grouped_result_reports_one_run_and_lists_the_rest() -> None:
    """rg measures Rg on two selections; the report covers one and names the other."""
    default = _report("rg")

    assert default.run == "Protein"
    assert default.all_runs == ["Protein", "Polymer Oligomers"]
    assert len(default.pairwise) == 5
    assert any("also reported runs" in text for text in default.warnings)

    other = _report("rg", run="Polymer Oligomers")

    assert other.run == "Polymer Oligomers"
    assert other.all_runs == ["Polymer Oligomers", "Protein"]
    # The control has no polymer, so this run covers the five polymer
    # conditions and compares every pair of them rather than against a control.
    assert [item.label for item in other.conditions] == [
        item.label for item in default.conditions if item.label != CONTROL
    ]
    assert len(other.pairwise) == 10


def test_unknown_run_is_a_typed_error() -> None:
    """Asking for a run the plugin did not report names the ones it did."""
    from polyzymd.analyses.exceptions import ProtocolError

    with pytest.raises(ProtocolError) as excinfo:
        _report("rg", run="Nonexistent Selection")

    assert "no run or metric named" in str(excinfo.value)
    assert "Protein" in (excinfo.value.hint or "")


def test_result_without_replicate_values_still_reports_an_interval() -> None:
    """Contacts stores only means and standard errors; intervals are rebuilt from them."""
    report = _report("contacts")
    condition = report.conditions[0]

    assert condition.replicate_values == []
    assert condition.ci95 is not None
    assert condition.ci_method == "student_t_from_sem"
    assert condition.ci95[0] < condition.mean < condition.ci95[1]
    assert report.pairwise[0].delta_ci95 is None
    assert any("rebuilt from stored standard errors" in text for text in report.warnings)
