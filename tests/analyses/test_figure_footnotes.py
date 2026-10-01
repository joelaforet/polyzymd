"""Every study figure that draws an uncertainty says what it is and what it spans.

Grossfield et al. (2018, LiveCoMS 1:5067, doi:10.33011/livecoms.1.1.5067) ask
that every figure describe the meaning and basis of its uncertainties. Each
test draws one figure type of :mod:`polyzymd.analyses.figures` from values
built without a trajectory and checks that its footnote names the interval
the figure draws, what it is the interval of, and the replicates it is
computed across. The last test checks that the autouse audit in
``tests/analyses/conftest.py`` fails a figure saved without one.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

from polyzymd.analyses import figures as study_figures  # noqa: E402
from polyzymd.analyses.shared import plotting  # noqa: E402
from polyzymd.analyses.timeseries import ReplicateValues  # noqa: E402
from polyzymd.config.analysis_settings import PlotSettings  # noqa: E402
from tests._support.analysis_testkit import replicate_values  # noqa: E402

PER_CONDITION = {"A": [1.0, 1.2, 1.4], "B": [2.0, 2.3, 2.5]}


@pytest.fixture()
def saved(monkeypatch) -> dict[str, object]:
    """Keep every saved figure open, by file stem, so its footnote can be read."""
    kept: dict[str, object] = {}
    original = plotting.save_figure

    def keep(fig, output_path, plot_settings, **kwargs):
        kept[Path(output_path).stem] = fig
        return original(fig, output_path, plot_settings, close=False)

    monkeypatch.setattr(plotting, "save_figure", keep)
    yield kept
    plt.close("all")


def _note(fig) -> str:
    """Return the footnote of ``fig``, unwrapped into one line."""
    text = next(item for item in fig.texts if item.get_gid() == plotting.FIGURE_NOTE_GID)
    return " ".join(text.get_text().split())


def _profile(per_condition: dict[str, list[float]]) -> ReplicateValues:
    """Labelled values on residues 1 to 3: replicate ``v`` holds ``v, 2v, 3v``."""
    source = replicate_values(per_condition).source
    rows = {
        label: [
            (index, np.array([v, 2 * v, 3 * v]), None, None, 20, None, None)
            for index, v in enumerate(values)
        ]
        for label, values in per_condition.items()
    }
    return ReplicateValues(source, "rmsf", "A", False, rows, labels=[1, 2, 3])


def test_condition_bars_name_the_t_interval_of_the_mean(tmp_path, saved) -> None:
    replicate_values(PER_CONDITION).plot(tmp_path, name="bars")
    assert _note(saved["bars"]) == (
        "Error bars: 95% Student t confidence interval of the condition mean across n = 3 "
        "replicates; production window t >= 10ns. Bars: condition means; points: "
        "per-replicate values."
    )


def test_different_counts_are_given_per_condition(tmp_path, saved) -> None:
    replicate_values({"A": [1.0, 1.2, 1.4], "B": [2.0, 2.3]}).plot(tmp_path, name="mixed")
    note = _note(saved["mixed"])
    assert "across the replicates of each condition (n per condition under each bar)" in note
    assert "n = " not in note


def test_a_one_replicate_condition_is_said_to_have_no_interval(tmp_path, saved) -> None:
    replicate_values({"A": [1.0, 1.2, 1.4], "B": [2.0]}).plot(tmp_path, name="one")
    assert _note(saved["one"]).endswith("a condition with one replicate has no interval.")


def test_grouped_bars_name_the_t_interval_of_the_mean(tmp_path, saved) -> None:
    values = replicate_values(PER_CONDITION)
    study_figures.plot_values([values, values], ["x", "y"], tmp_path, "grouped")
    assert _note(saved["grouped"]).startswith(
        "Error bars: 95% Student t confidence interval of the condition mean across n = 3 "
        "replicates;"
    )


def test_timeseries_band_is_the_interval_of_the_mean_at_each_time(tmp_path, saved) -> None:
    replicate_values(PER_CONDITION).source.plot(tmp_path, "series")
    note = _note(saved["series"])
    assert note.startswith(
        "Band: 95% Student t confidence interval of the condition mean at each time across "
        "n = 3 replicates; production window t >= 10ns."
    )
    assert "grey: the equilibration window" in note


def test_profile_band_is_the_interval_of_the_mean_at_each_label(tmp_path, saved) -> None:
    study_figures.plot_profile(_profile(PER_CONDITION), tmp_path, "profile", xlabel="Residue")
    assert _note(saved["profile"]).startswith(
        "Band: 95% Student t confidence interval of the condition mean at each residue across "
        "n = 3 replicates;"
    )


def test_decomposition_bands_are_intervals_of_each_line(tmp_path, saved) -> None:
    profile = _profile(PER_CONDITION)
    parts = {"RMSF": profile, "offset": profile}
    study_figures.plot_decomposition(parts, tmp_path, "parts", xlabel="Residue")
    assert _note(saved["parts"]).startswith(
        "Bands: 95% Student t confidence interval of each line's condition mean at each "
        "residue across n = 3 replicates;"
    )


@pytest.mark.parametrize(
    ("test", "method"), [("welch", "Welch t"), ("student", "Student t (pooled variance)")]
)
def test_difference_bands_name_the_test_of_the_difference(tmp_path, saved, test, method) -> None:
    profile = _profile({"A": [1.0, 1.2, 1.4], "B": [2.0, 2.3]})
    report = profile.compare(test=test)
    study_figures.plot_differences(profile, report, tmp_path, "delta", xlabel="Residue")
    note = _note(saved["delta"])
    assert note.startswith(
        f"Bands: 95% {method} confidence interval of the difference of condition means at "
        "each residue across the replicates of each condition (n of both conditions in each "
        "panel title); production window t >= 10ns."
    )
    assert "Benjamini-Hochberg correction over a family of 3 tests" in note


def test_distributions_draw_no_interval_and_claim_none(tmp_path, saved) -> None:
    source = replicate_values({"A": [1.0, 1.2, 1.4], "B": [2.0, 2.3, 2.5]}).source
    for item in (*source.series["A"], *source.series["B"]):
        item.values = item.values + np.linspace(0.0, 0.1, item.values.size)
    source.plot_distribution(output_dir=tmp_path, name="kde")
    note = _note(saved["kde"])
    assert note.startswith("Thick lines: Gaussian KDE of all replicates' frames pooled")
    assert "95%" not in note


def test_the_audit_fails_a_band_saved_without_a_footnote(tmp_path) -> None:
    fig, ax = plt.subplots()
    ax.fill_between([0, 1], [0.9, 1.9], [1.1, 2.1])
    with pytest.raises(AssertionError, match="without a footnote"):
        plotting.save_figure(fig, tmp_path / "bare.png", PlotSettings())
    plt.close(fig)
