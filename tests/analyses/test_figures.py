"""Tests for the figures of stored study results and for ``polyzymd analyze`` figure output.

Each replicate is an OpenMM run directory with one DCD segment of four
unit-mass atoms on a cross scaled on frame ``k`` by ``base + 0.01 * k``, so its
radius of gyration is that scale and frame ``k`` is at ``0.1 * k`` ns. The
0.25 ns window leaves frames 3 to 9. Each test checks the data drawn on the
matplotlib Axes against the stored values.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml
from click.testing import CliRunner
from matplotlib.collections import PathCollection
from matplotlib.container import BarContainer, ErrorbarContainer
from scipy.stats import gaussian_kde

import polyzymd as pz
from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.figures import reflected_kde
from polyzymd.analyses.shared import plotting
from polyzymd.analyses.shared.statistics import mean_sem_ci
from polyzymd.analyses.timeseries import _safe
from polyzymd.cli.analyze import analyze_command
from tests._support.analysis_testkit import (
    replicate_values,
    write_openmm_replicate,
    write_simulation_config,
)

pytest.importorskip("MDAnalysis")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

EQUILIBRATION = "0.25ns"
PAIRS = [
    {"label": "C1-C2", "selection_a": "name C1", "selection_b": "name C2", "threshold": 2.31},
    {"label": "mid-C3", "selection_a": "midpoint(name C1 C2)", "selection_b": "name C3"},
]


def rg(atoms):
    """Mass-weighted radius of gyration of ``atoms``."""
    return atoms.radius_of_gyration()


CALLS: list[int] = []


def counted_rg(atoms):
    """Mass-weighted radius of gyration of ``atoms``, counting each call in ``CALLS``."""
    CALLS.append(1)
    return atoms.radius_of_gyration()


@pytest.fixture()
def configs(tmp_path: Path) -> dict[str, Path]:
    """Two conditions of three replicates, B larger than A by one Å."""
    paths = {}
    for label, offset in (("A", 1.0), ("B", 2.0)):
        config = write_simulation_config(tmp_path / label, scratch=tmp_path / label / "scratch")
        for replicate in (1, 2, 3):
            scales = [offset + 0.1 * replicate + 0.01 * k for k in range(10)]
            write_openmm_replicate(config, replicate, scales)
        paths[label] = config
    return paths


@pytest.fixture()
def series(configs, tmp_path):
    """The per-frame Rg of every replicate, stored under ``tmp_path``."""
    study = pz.Study.from_configs(configs, equilibration=EQUILIBRATION)
    return study.timeseries(rg, pz.select("all"), unit="A", name="rg", output_dir=tmp_path)


@pytest.fixture()
def figures(monkeypatch) -> dict[str, object]:
    """Keep every saved figure open, by file stem, so its Axes can be read."""
    saved: dict[str, object] = {}
    original = plotting.save_figure

    def keep(fig, output_path, plot_settings, **kwargs):
        saved[Path(output_path).stem] = fig
        return original(fig, output_path, plot_settings, close=False)

    monkeypatch.setattr(plotting, "save_figure", keep)
    yield saved
    import matplotlib.pyplot as plt

    plt.close("all")


def _footnote(fig) -> str:
    """Return the footnote text of ``fig``, the text drawn at the bottom left."""
    return next(text.get_text() for text in fig.texts if text.get_position() == (0.01, 0.01))


def test_timeseries_plot_draws_every_replicate_and_the_mean(series, tmp_path, figures) -> None:
    """Each replicate is a line of its stored values, and the window is shaded from 0 ns."""
    path = series.plot()
    assert path == tmp_path / "figures" / "rg_timeseries.png" and path.is_file()
    ax = figures["rg_timeseries"].axes[0]
    lines = [(line.get_xdata(), line.get_ydata()) for line in ax.get_lines()]
    for items in series.series.values():
        for item in items:
            assert any(
                np.array_equal(x, item.times) and np.array_equal(y, item.values) for x, y in lines
            )
        mean = np.mean([item.values for item in items], axis=0)
        assert any(np.allclose(y, mean) and len(y) == len(mean) for _, y in lines)
    assert len(lines) == 6 + 2
    span = next(patch for patch in ax.patches)
    assert span.get_x() == pytest.approx(0.0) and span.get_width() == pytest.approx(0.25)
    legend = [text.get_text() for text in ax.get_legend().get_texts()]
    assert legend == ["A (n = 3)", "B (n = 3)", "equilibration window"]
    assert ax.get_xlabel() == "Time (ns)" and ax.get_ylabel() == "rg (Å)"
    note = _footnote(figures["rg_timeseries"])
    assert "Band: 95% CI (Student t) across n = 3 replicates" in note
    assert "production window t >= 0.25ns" in note


def test_distribution_plot_draws_pooled_and_replicate_kdes(series, tmp_path, figures) -> None:
    """The thick line is the KDE of the pooled frames, and the threshold line sits at 1.5 Å."""
    path = series.plot_distribution(1.5, output_dir=tmp_path / "out")
    assert path == tmp_path / "out" / "rg_distribution.png" and path.is_file()
    ax = figures["rg_distribution"].axes[0]
    lines = ax.get_lines()
    assert len(lines) == 6 + 2 + 1
    for label, items in series.series.items():
        pooled = np.concatenate([item.values for item in items])
        thick = next(line for line in lines if line.get_label() == f"{label} (n = 3)")
        assert np.allclose(thick.get_ydata(), gaussian_kde(pooled)(thick.get_xdata()))
        for item in items:
            assert any(
                np.allclose(line.get_ydata(), gaussian_kde(item.values)(line.get_xdata()))
                for line in lines
                if line.get_linewidth() == 0.8
            )
    threshold = next(line for line in lines if line.get_label() == "threshold 1.5 Å")
    assert list(threshold.get_xdata()) == [1.5, 1.5]


def test_values_plot_draws_means_intervals_and_every_replicate(series, tmp_path, figures) -> None:
    """Bar heights, error bars and points equal the summary means, intervals and values."""
    values = series.reduce("mean")
    path = values.plot(tmp_path, name="rg_comparison")
    assert path == tmp_path / "rg_comparison.png" and path.is_file()
    fig = figures["rg_comparison"]
    ax = fig.axes[0]
    bars = next(item for item in ax.containers if isinstance(item, BarContainer))
    dots = [item for item in ax.collections if isinstance(item, PathCollection)]
    summary = {row.label: row for row in values.summary().conditions}
    for index, label in enumerate(["A", "B"]):
        stats = mean_sem_ci(values.values[label])
        assert bars[index].get_height() == pytest.approx(summary[label].mean)
        low, high = bars.errorbar.lines[2][0].get_segments()[index][:, 1]
        assert (low, high) == pytest.approx((stats.ci_low, stats.ci_high))
        assert (low, high) == pytest.approx(tuple(summary[label].ci95))
        points = dots[index].get_offsets()[:, 1]
        assert list(points) == pytest.approx(values.values[label])
    ticks = [text.get_text() for text in ax.get_xticklabels()]
    assert ticks == ["A\nn = 3", "B\nn = 3"]
    note = _footnote(fig)
    assert "Error bars: 95% CI (Student t) across n = 3 replicates" in note
    assert "production window t >= 0.25ns" in note


def test_single_replicates_get_no_interval_claim(tmp_path, figures) -> None:
    """With one replicate per condition the footnote says no interval is drawn."""
    replicate_values({"A": [1.0], "B": [2.0]}).plot(tmp_path, name="single")
    ax = figures["single"].axes[0]
    assert not any(isinstance(item, ErrorbarContainer) for item in ax.containers)
    assert _footnote(figures["single"]).startswith("No interval: every condition has one replicate")


def test_cli_rg_writes_the_legacy_figures_and_records_them(configs, tmp_path) -> None:
    """rg writes its time series, comparison and distribution; --no-plots writes none."""
    base = ["rg", "-c", str(configs["A"]), "-c", str(configs["B"]), "--eq", EQUILIBRATION]
    base += ["--set", "selection=all", "--format", "json"]
    result = CliRunner().invoke(analyze_command, [*base, "--output-dir", str(tmp_path / "plots")])
    assert result.exit_code == 0, result.output
    folder = tmp_path / "plots" / "figures" / "rg"
    assert json.loads(result.stdout)["provenance"]["output_paths"]["figures"] == str(folder)
    assert sorted(path.name for path in folder.iterdir()) == [
        "rg_comparison.png",
        "rg_distribution.png",
        "rg_timeseries.png",
    ]

    quiet = CliRunner().invoke(
        analyze_command, [*base, "--output-dir", str(tmp_path / "quiet"), "--no-plots"]
    )
    assert quiet.exit_code == 0, quiet.output
    assert not (tmp_path / "quiet" / "figures").exists()
    assert "figures" not in json.loads(quiet.stdout)["provenance"]["output_paths"]


def test_cli_triad_draws_each_pair_with_its_threshold(configs, tmp_path, figures) -> None:
    """The triad draws each pair's distribution and every fraction, with each pair's threshold."""
    pairs = tmp_path / "pairs.yaml"
    pairs.write_text(yaml.safe_dump(PAIRS))
    arguments = ["catalytic_triad", "-c", str(configs["A"]), "-c", str(configs["B"])]
    arguments += ["--eq", EQUILIBRATION, "--set", f"pairs={pairs}", "--set", "threshold=1.135"]
    arguments += ["--output-dir", str(tmp_path)]
    result = CliRunner().invoke(analyze_command, arguments)
    assert result.exit_code == 0, result.output
    folder = tmp_path / "figures" / "catalytic_triad"
    assert sorted(path.name for path in folder.iterdir()) == [
        "triad_fraction_C1-C2_below_2.31_A.png",
        "triad_fraction_mid-C3_below_1.135_A.png",
        "triad_fraction_simultaneous.png",
        "triad_kde_C1-C2.png",
        "triad_kde_mid-C3.png",
        "triad_kde_panel.png",
        "triad_threshold_bars.png",
    ]
    for stem, threshold in (("triad_kde_C1-C2", 2.31), ("triad_kde_mid-C3", 1.135)):
        lines = figures[stem].axes[0].get_lines()
        assert list(lines[-1].get_xdata()) == [threshold, threshold]
    stored = np.load(
        tmp_path / "polyzymd_results/catalytic_triad_simultaneous/A/replicate_1/series.npz"
    )
    ax = figures["triad_fraction_simultaneous"].axes[0]
    dots = [item for item in ax.collections if isinstance(item, PathCollection)]
    assert dots[0].get_offsets()[0, 1] == pytest.approx(stored["values"].mean())
    assert ax.get_ylim() == (0, 1.05)


def _stored_means(root: Path, name: str, label: str) -> list[float]:
    """Return the per-replicate means of the stored series ``name`` of condition ``label``."""
    folder = root / "polyzymd_results" / _safe(name) / label
    return [
        float(np.load(folder / f"replicate_{r}/series.npz")["values"].mean()) for r in (1, 2, 3)
    ]


def _bars(ax) -> list[tuple[list[float], list[tuple[float, float]], np.ndarray]]:
    """Return, per condition, the bar heights, the error bar ends and the point values."""
    bars = [item for item in ax.containers if isinstance(item, BarContainer)]
    dots = [item.get_offsets()[:, 1] for item in ax.collections if isinstance(item, PathCollection)]
    rows = []
    for index, container in enumerate(bars):
        ends = [tuple(seg[:, 1]) for seg in container.errorbar.lines[2][0].get_segments()]
        points = dots[index * len(container) : (index + 1) * len(container)]
        rows.append(([patch.get_height() for patch in container], ends, points))
    return rows


def test_cli_triad_grouped_figures_hold_the_stored_values(configs, tmp_path, figures) -> None:
    """The threshold bars hold every pair's fraction and All pairs; the panel holds each threshold."""
    pairs = tmp_path / "pairs.yaml"
    pairs.write_text(yaml.safe_dump(PAIRS))
    arguments = ["catalytic_triad", "-c", str(configs["A"]), "-c", str(configs["B"])]
    arguments += ["--eq", EQUILIBRATION, "--set", f"pairs={pairs}", "--set", "threshold=2.2"]
    arguments += ["--output-dir", str(tmp_path)]
    assert CliRunner().invoke(analyze_command, arguments).exit_code == 0
    names = [
        "catalytic_triad_C1-C2_below 2.31 A",
        "catalytic_triad_mid-C3_below 2.2 A",
        "catalytic_triad_simultaneous",
    ]
    ax = figures["triad_threshold_bars"].axes[0]
    ticks = [text.get_text() for text in ax.get_xticklabels()]
    assert ticks == ["C1-C2 below 2.31 A", "mid-C3 below 2.2 A", "All pairs"]
    rows = _bars(ax)
    assert len(rows) == 2
    for label, (heights, ends, points) in zip(["A", "B"], rows, strict=True):
        for group, name in enumerate(names):
            means = _stored_means(tmp_path, name, label)
            stats = mean_sem_ci(means)
            assert heights[group] == pytest.approx(stats.mean)
            if stats.sem:
                assert ends[group] == pytest.approx((stats.ci_low, stats.ci_high))
            assert sorted(points[group]) == pytest.approx(sorted(means))
    assert "across n = 3 replicates" in _footnote(figures["triad_threshold_bars"])
    panel = figures["triad_kde_panel"].axes
    assert [ax.get_title() for ax in panel] == ["C1-C2 distance", "mid-C3 distance"]
    for ax, threshold in zip(panel, (2.31, 2.2), strict=True):
        assert list(ax.get_lines()[-1].get_xdata()) == [threshold, threshold]
    quiet = CliRunner().invoke(
        analyze_command, [*arguments[:-1], str(tmp_path / "q"), "--no-plots"]
    )
    assert quiet.exit_code == 0 and not (tmp_path / "q" / "figures").exists()


def test_cli_distances_writes_the_grouped_figures(configs, tmp_path) -> None:
    """distances writes its threshold bars and KDE panel next to the per-pair figures."""
    pairs = tmp_path / "pairs.yaml"
    pairs.write_text(yaml.safe_dump(PAIRS))
    arguments = ["distances", "-c", str(configs["A"]), "--eq", EQUILIBRATION]
    arguments += ["--set", f"pairs={pairs}", "--output-dir", str(tmp_path)]
    assert CliRunner().invoke(analyze_command, arguments).exit_code == 0
    written = {path.name for path in (tmp_path / "figures" / "distances").iterdir()}
    assert {"distance_threshold_bars.png", "distance_kde_panel.png"} <= written


def test_grouped_plots_refuse_mixed_units(series, tmp_path, figures) -> None:
    """Results or series of different units are refused; one unit draws one panel per series."""
    doubled = series.transform(lambda values: 2 * values, name="double")
    below = series.transform(lambda values: values < 1.5, unit=None, name="below")
    with pytest.raises(ProtocolError, match="share one unit"):
        pz.plot_values([series.reduce("mean"), below.reduce("fraction")], output_dir=tmp_path)
    with pytest.raises(ProtocolError, match="share one unit"):
        pz.plot_distributions([series, below], output_dir=tmp_path)
    path = pz.plot_distributions([series, doubled], [1.5, 3.0], output_dir=tmp_path, name="two")
    assert path == tmp_path / "two.png" and path.is_file()
    panel = figures["two"].axes
    for ax, source, threshold in zip(panel, (series, doubled), (1.5, 3.0), strict=True):
        pooled = np.concatenate([item.values for item in source.series["A"]])
        thick = next(line for line in ax.get_lines() if line.get_label() == "A (n = 3)")
        assert np.allclose(thick.get_ydata(), gaussian_kde(pooled)(thick.get_xdata()))
        assert list(ax.get_lines()[-1].get_xdata()) == [threshold, threshold]
    assert panel[-1].get_xlabel() == "value (Å)"


def test_importing_polyzymd_does_not_import_matplotlib() -> None:
    """matplotlib loads only when a figure is drawn."""
    code = (
        "import sys, polyzymd, polyzymd.analyses, polyzymd.analyses.figures, "
        "polyzymd.analyses.timeseries, polyzymd.cli.analyze; "
        "sys.exit('matplotlib' in sys.modules)"
    )
    assert subprocess.run([sys.executable, "-c", code], check=False).returncode == 0


class TestReflectedKDE:
    """Boundary correction by reflection (Schuster 1985; Silverman 1986, section 2.10)."""

    @staticmethod
    def _bandwidth(values: np.ndarray) -> float:
        return gaussian_kde(values).factor * float(np.std(values, ddof=1))

    def test_lower_bound_keeps_the_density_inside_the_support(self) -> None:
        """A half-normal sample gives a curve from 0 that integrates to 1 and is plain far away."""
        values = np.abs(np.random.default_rng(1).normal(size=2000))
        grid, density = reflected_kde(values, (0.0, None), points=20001)
        assert grid[0] == 0.0 and (density >= 0).all()
        assert np.trapezoid(density, grid) == pytest.approx(1.0, abs=1e-3)
        far = grid > values.min() + 8 * self._bandwidth(values)
        assert np.allclose(density[far], gaussian_kde(values)(grid[far]), rtol=0, atol=1e-12)
        assert np.trapezoid(gaussian_kde(values)(grid), grid) < 0.99

    def test_upper_bound_and_far_or_missing_bounds(self) -> None:
        """An upper bound caps the grid; a far bound or no bound leaves the plain KDE."""
        values = 1.0 - np.abs(np.random.default_rng(2).normal(scale=0.1, size=2000))
        grid, density = reflected_kde(values, (None, 1.0), points=20001)
        assert grid[-1] == 1.0 and np.trapezoid(density, grid) == pytest.approx(1.0, abs=1e-3)
        far = np.random.default_rng(3).normal(loc=20.0, size=500)
        width = 3 * self._bandwidth(far)
        plain_grid = np.linspace(far.min() - width, far.max() + width, 200)
        for bounds in ((0.0, None), (None, None)):
            grid, density = reflected_kde(far, bounds)
            assert np.array_equal(grid, plain_grid)
            assert np.allclose(density, gaussian_kde(far)(plain_grid), rtol=0, atol=1e-12)

    def test_distribution_figure_uses_the_series_bounds(self, series, tmp_path, figures) -> None:
        """A series with bounds (0, None) is drawn from its reflected KDE with the footnote."""
        assert series.bounds == (None, None)
        record = json.loads((series.series["A"][0].path / "record.json").read_text())
        assert record["bounds"] == [None, None]
        bounded = series.transform(lambda values: values - 1.0, name="shifted", bounds=(0.0, None))
        assert bounded.transform(lambda values: values, name="same").bounds == (0.0, None)
        bounded.plot_distribution(output_dir=tmp_path)
        fig = figures["shifted_distribution"]
        pooled = np.concatenate([item.values for item in bounded.series["A"]])
        thick = next(line for line in fig.axes[0].get_lines() if line.get_label() == "A (n = 3)")
        grid, density = reflected_kde(pooled, (0.0, None))
        assert np.array_equal(thick.get_xdata(), grid) and np.allclose(thick.get_ydata(), density)
        assert "corrected by reflection (Schuster 1985; Silverman 1986)" in _footnote(fig)
        series.plot_distribution(output_dir=tmp_path)
        assert "reflection" not in _footnote(figures["rg_distribution"])


def test_new_bounds_reuse_the_stored_series_and_update_its_record(configs, tmp_path) -> None:
    """Bounds change only the drawing, so a stored series is reused and its record updated."""
    study = pz.Study.from_configs(configs, equilibration=EQUILIBRATION)
    arguments = (counted_rg, pz.select("all"))
    CALLS.clear()
    first = study.timeseries(*arguments, unit="A", name="rg", output_dir=tmp_path)
    measured = len(CALLS)
    assert measured == 6 * 7
    again = study.timeseries(
        *arguments, unit="A", name="rg", output_dir=tmp_path, bounds=(0.0, None)
    )
    assert len(CALLS) == measured and again.bounds == (0.0, None)
    for label in ("A", "B"):
        for old, new in zip(first.series[label], again.series[label], strict=True):
            assert np.array_equal(old.values, new.values)
            record = json.loads((new.path / "record.json").read_text())
            assert record["bounds"] == [0, None] and "versions" in record
