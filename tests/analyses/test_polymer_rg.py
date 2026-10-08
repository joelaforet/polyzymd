"""Tests for ``polyzymd analyze polymer_rg`` on OpenMM run directories.

Each replicate has a three-atom protein in chain A and two polymer chains in
chain C, each a ring of eight carbon atoms whose radius changes from frame to
frame, so each chain's radius of gyration on each frame is its ring radius.
The second ring sits across the box face on x and is written wrapped, as a
simulation writes it. Bonds come from the PDB's CONECT records.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

from polyzymd.analyses import analyze
from polyzymd.analyses.exceptions import NoMatchingAtomsError, ProtocolError
from polyzymd.cli.analyze import analyze_command
from tests._support.analysis_testkit import write_openmm_frames, write_simulation_config

mda = pytest.importorskip("MDAnalysis")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

BOX = [30.0, 30.0, 30.0, 90.0, 90.0, 90.0]
RING = 8
N_FRAMES = 6
PROTEIN = np.array([[15.0, 15.0, 5.0], [16.5, 15.0, 5.0], [18.0, 15.0, 5.0]])
CENTERS = ([10.0, 15.0, 15.0], [29.0, 15.0, 22.0])


def _ring(radius: float, center) -> np.ndarray:
    angles = 2 * np.pi * np.arange(RING) / RING
    return np.asarray(center) + radius * np.stack(
        [np.cos(angles), np.sin(angles), np.zeros(RING)], axis=1
    )


def _radii(label: str, replicate: int) -> np.ndarray:
    """Ring radii, ``(frames, 2)``: B's chains are 1.5 A larger than A's."""
    base = 3.0 + (1.5 if label == "B" else 0.0) + 0.1 * replicate
    steps = 0.2 * np.arange(N_FRAMES)
    return np.stack([base + steps, base + 1.0 - steps / 2], axis=1)


def _write(folder: Path, label: str, *, polymer: bool = True) -> Path:
    config = write_simulation_config(folder / label, scratch=folder / label / "scratch")
    n_poly = 2 * RING if polymer else 0
    resindex = [0, 0, 0] + [1 + i // RING for i in range(n_poly)]
    bonds = [(0, 1), (1, 2)]
    for chain in range(n_poly // RING):
        start = 3 + chain * RING
        bonds += [(start + i, start + i + 1) for i in range(RING - 1)]
    for replicate in (1, 2, 3):
        radii = _radii(label, replicate)
        frames = []
        for row in radii:
            rings = [_ring(r, c) for r, c in zip(row, CENTERS)] if polymer else []
            frame = np.vstack([PROTEIN, *rings])
            frame[:, 0] %= BOX[0]
            frames.append(frame)
        write_openmm_frames(
            config,
            replicate,
            np.asarray(frames),
            resindex,
            resids=[1, 1, 2][: max(resindex) + 1],
            names=["C"] * len(resindex),
            resnames=["ALA", "SBM", "SBM"][: max(resindex) + 1],
            elements=["C"] * len(resindex),
            chain_ids=["A"] * 3 + ["C"] * n_poly,
            dimensions=BOX,
            bonds=bonds,
        )
    return config


@pytest.fixture()
def configs(tmp_path: Path) -> dict[str, Path]:
    return {label: _write(tmp_path / "runs", label) for label in ("A", "B")}


def _values(report) -> dict[str, list[float]]:
    return {row.label: row.replicate_values for row in report.conditions}


def test_mean_rg_is_each_replicate_mean_over_frames_of_the_chain_mean(
    configs, tmp_path: Path
) -> None:
    report = analyze(
        "polymer_rg",
        [configs["A"], configs["B"]],
        equilibration="0ns",
        output_dir=tmp_path / "out",
        plots=False,
    )
    expected = {label: [float(_radii(label, r).mean()) for r in (1, 2, 3)] for label in ("A", "B")}
    for label, values in _values(report).items():
        np.testing.assert_allclose(values, expected[label], atol=1e-4)
    assert report.run == "mean_rg" and report.all_runs == ["mean_rg", "chain_rg"]
    assert report.unit == "A"
    (pair,) = report.pairwise
    assert pair.delta == pytest.approx(1.5, abs=1e-4)
    assert report.provenance.settings["selection"] == "chainid C"
    assert set(report.provenance.settings["chains"].values()) == {2}


def test_chain_rg_reports_each_chain_mean_by_position(configs, tmp_path: Path) -> None:
    report = analyze(
        "polymer_rg",
        [configs["A"]],
        equilibration="0ns",
        output_dir=tmp_path / "out",
        plots=False,
        run="chain_rg",
    )
    rows = {(row.label, row.entry): row.replicate_values for row in report.conditions}
    for chain in (1, 2):
        expected = [float(_radii("A", r)[:, chain - 1].mean()) for r in (1, 2, 3)]
        np.testing.assert_allclose(rows[("A", str(chain))], expected, atol=1e-4)


def test_a_control_without_polymer_is_left_out_with_a_warning(tmp_path: Path) -> None:
    control = _write(tmp_path / "runs", "A", polymer=False)
    treated = _write(tmp_path / "runs", "B")
    report = analyze(
        "polymer_rg",
        [control, treated],
        equilibration="0ns",
        output_dir=tmp_path / "out",
        plots=False,
    )
    assert list(_values(report)) == ["B"]
    assert not report.pairwise
    assert any("matched no atoms" in warning for warning in report.warnings)


def test_a_study_without_polymer_raises(tmp_path: Path) -> None:
    control = _write(tmp_path / "runs", "A", polymer=False)
    with pytest.raises(NoMatchingAtomsError):
        analyze("polymer_rg", [control], equilibration="0ns", output_dir=tmp_path, plots=False)


def test_an_unknown_run_is_refused(configs, tmp_path: Path) -> None:
    with pytest.raises(ProtocolError, match="no result named"):
        analyze("polymer_rg", [configs["A"]], equilibration="0ns", output_dir=tmp_path, run="x")


def test_polyzymd_analyze_polymer_rg_draws_its_figures(configs, tmp_path: Path) -> None:
    arguments = ["polymer_rg", "-c", str(configs["A"]), "-c", str(configs["B"])]
    arguments += ["--eq", "0ns", "--output-dir", str(tmp_path / "out")]
    result = CliRunner().invoke(analyze_command, arguments)
    assert result.exit_code == 0, result.output
    figures = {path.stem for path in (tmp_path / "out" / "figures" / "polymer_rg").iterdir()}
    assert {"polymer_rg_timeseries", "polymer_rg_comparison"} <= figures
