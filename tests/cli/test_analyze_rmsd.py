"""Tests for ``polyzymd analyze rmsd`` and the rmsd function it runs.

Each replicate is an OpenMM run directory with one DCD segment of four
atoms on a cross, scaled on frame ``k`` by ``base + 0.01 * k``, so its RMSD
from the cross scaled by ``c`` is ``|base + 0.01 * k - c|``.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

from polyzymd.analyses.functions import rmsd
from polyzymd.analyses.protocols import ProtocolReport
from polyzymd.cli.analyze import EXIT_ANALYSIS_ERROR, analyze_command
from tests._support.analysis_testkit import (
    CROSS,
    write_openmm_replicate,
    write_simulation_config,
)

mda = pytest.importorskip("MDAnalysis")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

EQUILIBRATION = "0.25ns"
GMX = shutil.which("gmx") or shutil.which(str(Path.home() / ".pixi/bin/gmx"))


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
def reference_pdb(tmp_path: Path) -> Path:
    """The unit cross as a PDB file."""
    universe = mda.Universe.empty(4, n_residues=1, atom_resindex=[0] * 4, trajectory=True)
    universe.add_TopologyAttr("names", ["C1", "C2", "C3", "C4"])
    universe.add_TopologyAttr("resnames", ["MOL"])
    universe.atoms.positions = np.asarray(CROSS)
    universe.atoms.write(str(tmp_path / "cross.pdb"))
    return tmp_path / "cross.pdb"


def _moving_frames(n_frames: int = 12) -> np.ndarray:
    """Random coordinates of 20 atoms that drift, rotate and translate over the frames."""
    from scipy.spatial.transform import Rotation

    rng = np.random.default_rng(7)
    shape = rng.normal(scale=4.0, size=(20, 3))
    return np.array(
        [
            Rotation.random(random_state=rng).apply(shape + rng.normal(scale=0.4, size=shape.shape))
            + rng.normal(scale=6.0, size=3)
            for _ in range(n_frames)
        ],
        dtype=np.float32,
    )


def test_rmsd_equals_mdanalysis_rmsd_analysis() -> None:
    """rmsd gives what MDAnalysis.analysis.rms.RMSD, which legacy rmsd ran, gives per frame."""
    from MDAnalysis.analysis.rms import RMSD

    coordinates = _moving_frames()
    universe = mda.Universe.empty(20, trajectory=True)
    universe.add_TopologyAttr("masses", [1.0] * 20)
    universe.load_new(coordinates, format="MEMORY")
    reference = mda.Merge(universe.atoms)
    reference.load_new(coordinates[4][np.newaxis], format="MEMORY")
    expected = RMSD(universe.atoms, reference.atoms, select="all").run().results.rmsd[:, 2]
    measured = [rmsd(universe.atoms, reference.atoms) for _ in universe.trajectory]
    assert measured == pytest.approx(expected, abs=1e-6)
    assert measured[4] == pytest.approx(0.0, abs=1e-5)


def test_cli_rmsd_prints_the_report(configs, reference_pdb: Path, tmp_path: Path) -> None:
    """External mode measures each frame from the unit cross and compares with the first config."""
    arguments = ["rmsd", "-c", str(configs["A"]), "-c", str(configs["B"]), "--eq", EQUILIBRATION]
    arguments += ["--set", "selection=all", "--set", "reference_mode=external"]
    arguments += ["--set", f"reference_file={reference_pdb}", "--output-dir", str(tmp_path)]
    result = CliRunner().invoke(analyze_command, arguments)
    assert result.exit_code == 0, result.output
    lines = result.stdout.strip().split("\n")
    assert lines[0].startswith("# polyzymd analyze rmsd  metric mean_rmsd  unit A  eq 0.25ns")
    assert lines[1].startswith("A  n 3  mean 0.26")
    assert lines[2].startswith("B  n 3  mean 1.26")
    assert lines[-1].startswith("verdict: B larger mean_rmsd than A")

    default = CliRunner().invoke(
        analyze_command, ["rmsd", "-c", str(configs["A"]), "--eq", EQUILIBRATION]
    )
    assert default.exit_code == EXIT_ANALYSIS_ERROR
    assert "matched no atoms" in default.stderr


def test_python_analyze_rmsd_gives_the_cli_report(configs, tmp_path: Path) -> None:
    """analyze("rmsd", ...) with the average reference gives the CLI report."""
    from polyzymd.analyses import analyze

    settings = {"selection": "all", "alignment_selection": "all", "reference_mode": "average"}
    report = analyze(
        "rmsd",
        [configs["A"], configs["B"]],
        equilibration=EQUILIBRATION,
        replicates=[1, 2],
        settings=settings,
        output_dir=tmp_path,
    )
    # Frames 3 to 9 average to the cross at base + 0.06, 0.03, 0.02, 0.01, 0, ... Å away.
    assert report.conditions[0].replicate_values == pytest.approx([0.12 / 7] * 2, abs=1e-5)
    assert report.frames_per_replicate == {"A": [7, 7], "B": [7, 7]}

    arguments = ["rmsd", "-c", str(configs["A"]), "-c", str(configs["B"]), "--eq", EQUILIBRATION]
    arguments += ["--replicates", "1-2", "--output-dir", str(tmp_path), "--format", "json"]
    for key, value in settings.items():
        arguments += ["--set", f"{key}={value}"]
    result = CliRunner().invoke(analyze_command, arguments)
    assert ProtocolReport.model_validate_json(result.stdout) == report


def test_python_analyze_rmsd_refuses_other_settings(configs) -> None:
    from polyzymd.analyses import analyze
    from polyzymd.analyses.exceptions import ProtocolError

    with pytest.raises(ProtocolError, match="no setting other than selection, alignment_sel"):
        analyze("rmsd", [configs["A"]], equilibration=EQUILIBRATION, settings={"runs": []})
    with pytest.raises(ProtocolError) as excinfo:
        analyze("not_an_analysis", [configs["A"]])
    assert "rmsd" in excinfo.value.hint.split("Use one of: ")[1].split(", ")


@pytest.mark.skipif(GMX is None, reason="needs the GROMACS gmx binary")
def test_rmsd_matches_gmx_rms(tmp_path: Path) -> None:
    """rmsd equals gmx rms with a least-squares fit on the same atoms, in nm times 10."""
    coordinates = _moving_frames()
    universe = mda.Universe.empty(20, n_residues=1, atom_resindex=[0] * 20, trajectory=True)
    universe.add_TopologyAttr("names", [f"C{i}" for i in range(20)])
    universe.add_TopologyAttr("resnames", ["MOL"])
    universe.add_TopologyAttr("masses", [1.0] * 20)
    universe.load_new(coordinates, format="MEMORY", dimensions=[99.0, 99.0, 99.0, 90, 90, 90])
    universe.atoms.write(str(tmp_path / "ref.pdb"))
    with mda.Writer(str(tmp_path / "traj.xtc"), n_atoms=20) as writer:
        for _ in universe.trajectory:
            writer.write(universe.atoms)
    subprocess.run(
        [GMX, "rms", "-s", "ref.pdb", "-f", "traj.xtc", "-o", "rmsd.xvg", "-nomw", "-xvg", "none"],
        cwd=tmp_path,
        input="0\n0\n",
        text=True,
        capture_output=True,
        check=True,
    )
    gmx = np.loadtxt(tmp_path / "rmsd.xvg")[:, 1] * 10.0
    xtc = mda.Universe(str(tmp_path / "ref.pdb"), str(tmp_path / "traj.xtc"))
    reference = mda.Universe(str(tmp_path / "ref.pdb"))
    measured = [rmsd(xtc.atoms, reference.atoms) for _ in xtc.trajectory]
    assert measured == pytest.approx(gmx, abs=1e-5)
