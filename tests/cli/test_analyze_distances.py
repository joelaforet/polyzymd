"""Tests for ``polyzymd analyze distances``, the pair distance, and the retired ``catalytic_triad``.

Each replicate is an OpenMM run directory with one DCD segment of four atoms
on a cross scaled on frame ``k`` by ``s = base + 0.01 * k``, so atoms C1 and C2
are ``2 s`` apart and the midpoint of C1 and C2, the origin, is ``s`` from C3.
The trajectory has no box, so no periodic image is used.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import yaml
from click.testing import CliRunner

from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.functions import pair_distance
from polyzymd.cli.analyze import EXIT_ANALYSIS_ERROR, analyze_command
from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config

mda = pytest.importorskip("MDAnalysis")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

EQUILIBRATION = "0.25ns"
PRODUCTION = np.arange(3, 10)
PAIRS = [
    {"label": "C1-C2", "selection_a": "name C1", "selection_b": "name C2", "threshold": 2.31},
    {
        "label": "mid-C3",
        "selection_a": "midpoint(name C1 C2)",
        "selection_b": "name C3",
        "below_label": "bound",
    },
]


def _scales(offset: float, replicate: int) -> np.ndarray:
    """Per-frame scale of the cross in one replicate's production frames."""
    return offset + 0.1 * replicate + 0.01 * PRODUCTION


@pytest.fixture()
def configs(tmp_path: Path) -> dict[str, Path]:
    """Two conditions of three replicates, B's cross larger than A's by one."""
    paths = {}
    for label, offset in (("A", 1.0), ("B", 2.0)):
        config = write_simulation_config(tmp_path / label, scratch=tmp_path / label / "scratch")
        for replicate in (1, 2, 3):
            write_openmm_replicate(
                config, replicate, [offset + 0.1 * replicate + 0.01 * k for k in range(10)]
            )
        paths[label] = config
    return paths


def test_pair_distance_points_and_minimum_image() -> None:
    """Single atoms, midpoint and center of mass, with and without the periodic box."""
    universe = mda.Universe.empty(3, trajectory=True)
    universe.add_TopologyAttr("masses", [1.0, 3.0, 1.0])
    universe.atoms.positions = [[1.0, 0.0, 0.0], [3.0, 0.0, 0.0], [9.0, 0.0, 0.0]]
    universe.dimensions = [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]
    first, pair, last = universe.atoms[[0]], universe.atoms[:2], universe.atoms[[2]]
    assert pair_distance(first, last) == pytest.approx(2.0)
    assert pair_distance(first, last, pbc=False) == pytest.approx(8.0)
    assert pair_distance(pair, last, mode_a="midpoint") == pytest.approx(3.0)
    assert pair_distance(pair, last, mode_a="com", pbc=False) == pytest.approx(6.5)
    universe.dimensions = None
    with pytest.warns(UserWarning, match="no valid box"):
        assert pair_distance(first, last) == pytest.approx(8.0)
    with pytest.raises(ValueError, match="exactly one atom"):
        pair_distance(pair, last)


def test_cli_distances_reports_each_pair(configs, tmp_path: Path) -> None:
    """--set pairs=FILE measures every pair, and --run picks the mean or the fraction."""
    pairs = tmp_path / "pairs.yaml"
    pairs.write_text(yaml.safe_dump(PAIRS))
    base = ["distances", "-c", str(configs["A"]), "-c", str(configs["B"]), "--eq", EQUILIBRATION]
    base += ["--set", f"pairs={pairs}", "--output-dir", str(tmp_path)]

    result = CliRunner().invoke(analyze_command, [*base, "--format", "json"])
    assert result.exit_code == 0, result.output
    report = json.loads(result.stdout)
    assert report["analysis"] == "distances" and report["run"] == "C1-C2"
    assert report["all_runs"] == ["C1-C2", "C1-C2 below 2.31 A", "mid-C3", "mid-C3 bound"]
    assert report["metric"] == "mean_distance" and report["unit"] == "A"
    values = report["conditions"][0]["replicate_values"]
    assert values == pytest.approx([2 * _scales(1.0, r).mean() for r in (1, 2, 3)], abs=1e-5)

    fraction = CliRunner().invoke(analyze_command, [*base, "--run", "C1-C2 below 2.31 A"])
    assert fraction.exit_code == 0, fraction.output
    lines = fraction.stdout.strip().split("\n")
    assert lines[0].startswith("# polyzymd analyze distances  metric fraction_below_threshold")
    assert "run C1-C2 below 2.31 A" in lines[0]
    # Replicate A1 has C1-C2 distances 2.26 to 2.38, three of them below 2.31.
    assert lines[1].startswith("A  n 3  mean 0.1429")
    assert lines[2].startswith("B  n 3  mean 0 ")

    unknown = CliRunner().invoke(analyze_command, [*base, "--run", "C1-C3"])
    assert unknown.exit_code == EXIT_ANALYSIS_ERROR
    assert "no result named 'C1-C3'" in unknown.stderr


def test_all_below_gives_the_fraction_with_every_pair_below_its_threshold(
    configs, tmp_path
) -> None:
    """Two stored pair distances combine, as in the triad routine, into one fraction."""
    import polyzymd as pz
    from polyzymd.analyses.functions import all_below

    study = pz.Study.from_configs({"A": configs["A"]}, equilibration=EQUILIBRATION)
    c1_c2 = study.timeseries(
        pair_distance,
        pz.select("name C1"),
        pz.select("name C2"),
        unit="A",
        name="c1_c2",
        output_dir=tmp_path,
    )
    mid_c3 = study.timeseries(
        pair_distance,
        pz.select("name C1 C2"),
        pz.select("name C3"),
        mode_a="midpoint",
        unit="A",
        name="mid_c3",
        output_dir=tmp_path,
    )
    both = c1_c2.transform(
        all_below, mid_c3, unit=None, bounds=(0.0, 1.0), name="both", thresholds=[2.31, 1.135]
    )
    # Replicate A1: C1-C2 below 2.31 in frames 3 to 5, mid-C3 below 1.135 in frame 3 only.
    assert both.reduce("fraction").values["A"] == pytest.approx([1 / 7, 0.0, 0.0])
    record = json.loads((tmp_path / "polyzymd_results/both/A/replicate_1/record.json").read_text())
    assert record["transform"]["qualname"] == "all_below"
    assert record["kwargs"] == {"thresholds": [2.31, 1.135]}
    assert len(record["inputs"]) == 2


def test_catalytic_triad_is_refused_with_the_routine_and_distances(configs, tmp_path) -> None:
    """catalytic_triad is a routine now: Python and the CLI name its page and distances."""
    from polyzymd.analyses import analyze

    url = "https://polyzymd.readthedocs.io/en/latest/how_to/analysis_triad_quickstart.html"
    with pytest.raises(ProtocolError, match="routine on the study API") as raised:
        analyze("catalytic_triad", [configs["A"]], equilibration=EQUILIBRATION)
    assert url in raised.value.hint
    assert "polyzymd analyze distances -c <config.yaml> --set pairs=<pairs.yaml>" in (
        raised.value.hint
    )
    assert ".claude/skills/polyzymd-analyze/SKILL.md" in raised.value.hint

    result = CliRunner().invoke(analyze_command, ["catalytic_triad", "-c", str(configs["A"])])
    assert result.exit_code == EXIT_ANALYSIS_ERROR
    assert "error: catalytic_triad is no longer a polyzymd analyze analysis" in result.stderr
    assert f"fix: Follow {url}" in result.stderr
    assert not (tmp_path / "polyzymd_results").exists()


@pytest.mark.parametrize(
    ("settings", "match"),
    [
        ({}, "needs pairs"),
        ({"pairs": [{"label": "x", "selection_a": "name C1"}]}, "needs pairs"),
        ({"pairs": PAIRS, "cutoff": 3.0}, "no setting other than pairs, threshold, use_pbc"),
        ({"pairs": "/nonexistent/pairs.yaml"}, "cannot read the pairs file"),
    ],
)
def test_bad_pair_settings_are_rejected(configs, tmp_path, settings, match) -> None:
    """Missing or malformed pairs and unknown settings raise ProtocolError."""
    from polyzymd.analyses import analyze

    with pytest.raises(ProtocolError, match=match):
        analyze(
            "distances",
            [configs["A"]],
            equilibration=EQUILIBRATION,
            settings=settings,
            output_dir=tmp_path,
        )
