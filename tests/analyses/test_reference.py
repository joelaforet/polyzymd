"""Known-answer tests for pz.reference, build_reference and the centroid frame.

The study replicates are the four-atom crosses of
:func:`tests._support.analysis_testkit.write_openmm_replicate`: frame ``k``
is the unit cross scaled by ``scales[k]``, so the RMSD between frames ``j``
and ``k`` after superposition is ``|scales[j] - scales[k]|``.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import polyzymd as pz
from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.reference import build_reference
from polyzymd.analyses.shared.centroid import _find_frame_closest_to_aligned_mean
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

# Frame k is at 0.1 * k ns, and a 0.25 ns window leaves frames 3 to 9.
EQUILIBRATION = "0.25ns"
SCALES = [1.0 + 0.01 * k for k in range(10)]


def rmsd_to(atoms, reference):
    """RMSD after superposition, written here so the test does not use the shipped function."""
    from MDAnalysis.analysis.rms import rmsd

    return rmsd(atoms.positions, reference.positions, center=True, superposition=True)


@pytest.fixture()
def study(tmp_path: Path) -> pz.Study:
    """One condition, one replicate of scaled crosses."""
    config = write_simulation_config(tmp_path / "A", scratch=tmp_path / "A" / "scratch")
    write_openmm_replicate(config, 1, SCALES)
    return pz.Study.from_configs({"A": config}, equilibration=EQUILIBRATION)


def _write_cross(path: Path, scale: float) -> Path:
    """Write the cross scaled by ``scale`` as a PDB file."""
    universe = mda.Universe.empty(4, n_residues=1, atom_resindex=[0] * 4, trajectory=True)
    universe.add_TopologyAttr("names", ["C1", "C2", "C3", "C4"])
    universe.add_TopologyAttr("resnames", ["MOL"])
    universe.atoms.positions = np.asarray(CROSS) * scale
    universe.atoms.write(str(path))
    return path


def _rotated(shape: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Rotate ``shape`` by a random rotation and translate it at random."""
    from scipy.spatial.transform import Rotation

    return Rotation.random(random_state=rng).apply(shape) + rng.normal(scale=5.0, size=3)


def test_reference_checks_its_arguments(tmp_path: Path) -> None:
    with pytest.raises(ProtocolError, match="Unknown reference mode"):
        pz.reference("first", "all")
    with pytest.raises(ProtocolError, match="production frame from 1"):
        pz.reference("frame", "all", frame=0)
    with pytest.raises(ProtocolError, match="existing structure file"):
        pz.reference("external", "all", file=tmp_path / "missing.pdb")
    assert pz.reference("average", "name C1 C2").alignment == "name C1 C2"
    assert pz.reference("centroid", "all", frame=3).frame is None


@pytest.mark.parametrize(
    ("mode", "options", "expected"),
    [
        ("frame", {"frame": 2}, [abs(s - SCALES[4]) for s in SCALES[3:]]),
        ("average", {}, [abs(s - 1.06) for s in SCALES[3:]]),
        ("centroid", {}, [abs(s - 1.06) for s in SCALES[3:]]),
    ],
)
def test_modes_on_scaled_crosses(study, tmp_path, mode, options, expected) -> None:
    """Frame 2 is trajectory frame 4, and the mean and centroid are the cross at 1.06."""
    series = study.timeseries(
        rmsd_to,
        pz.select("all"),
        pz.reference(mode, "all", **options),
        unit="A",
        name=mode,
        output_dir=tmp_path,
    ).series["A"][0]
    assert series.values == pytest.approx(expected, abs=1e-5)
    record = json.loads((series.path / "record.json").read_text())
    assert record["arguments"]["args"][1]["reference"]["mode"] == mode
    if mode != "average":
        assert record["chosen"]["1"] == {
            "frame": 4 if mode == "centroid" else 2,
            "trajectory_frame": 6 if mode == "centroid" else 4,
        }


def test_external_file_is_hashed_and_a_changed_file_recomputes(study, tmp_path) -> None:
    path = _write_cross(tmp_path / "ref.pdb", 1.0)

    def run():
        return study.timeseries(
            rmsd_to,
            pz.select("all"),
            pz.reference("external", "all", file=path),
            unit="A",
            output_dir=tmp_path,
        ).series["A"][0]

    first = run()
    record = json.loads((first.path / "record.json").read_text())
    stored = record["arguments"]["args"][1]["reference"]["file"]
    assert stored["path"] == str(path.resolve()) and len(stored["sha256"]) == 64
    assert first.values == pytest.approx([s - 1.0 for s in SCALES[3:]], abs=1e-5)
    _write_cross(path, 2.0)
    assert run().values == pytest.approx([2.0 - s for s in SCALES[3:]], abs=1e-5)


def test_plain_path_argument_is_hashed(study, tmp_path) -> None:
    path = tmp_path / "offset.txt"
    path.write_text("0.5")

    def offset_rg(atoms, file):
        return atoms.radius_of_gyration() + float(Path(file).read_text())

    def run():
        return study.timeseries(
            offset_rg, pz.select("all"), str(path), unit="A", output_dir=tmp_path
        ).series["A"][0]

    first = run()
    record = json.loads((first.path / "record.json").read_text())
    assert record["arguments"]["args"][1]["file"]["path"] == str(path)
    path.write_text("1.5")
    assert run().values == pytest.approx(first.values + 1.0, abs=1e-5)


def test_frame_past_the_window_and_atom_count_mismatch(study, tmp_path) -> None:
    with pytest.raises(ProtocolError, match="past the last of 7 production frames"):
        study.timeseries(
            rmsd_to,
            pz.select("all"),
            pz.reference("frame", "all", frame=8),
            unit="A",
            output_dir=tmp_path,
        )
    path = _write_cross(tmp_path / "ref.pdb", 1.0)
    path.write_text("".join(line for line in path.read_text().splitlines(True) if "C4" not in line))
    with pytest.raises(ProtocolError, match="has 3 atoms"):
        study.timeseries(
            rmsd_to,
            pz.select("all"),
            pz.reference("external", "all", file=path),
            unit="A",
            output_dir=tmp_path,
        )


def test_rotated_frames_leave_the_live_trajectory_unchanged() -> None:
    """Average and centroid undo random rigid motions, and the trajectory is not moved."""
    rng = np.random.default_rng(3)
    shape = rng.normal(scale=3.0, size=(8, 3))
    bend = rng.normal(size=(8, 3))
    amounts = [0.0, 0.1, 0.45, 0.9, 1.0]
    coordinates = np.array([_rotated(shape + a * bend, rng) for a in [9.0, *amounts]])
    universe = mda.Universe.empty(8, trajectory=True)
    universe.load_new(coordinates.astype(np.float32), format="MEMORY")
    before = universe.trajectory.timeseries(order="fac").copy()
    frames = np.arange(1, 6)

    average, _ = build_reference(pz.reference("average", "all"), universe, frames)
    target = shape + np.mean(amounts) * bend
    from MDAnalysis.analysis.rms import rmsd

    assert rmsd(average.positions, target, center=True, superposition=True) < 0.05
    centroid, chosen = build_reference(pz.reference("centroid", "all"), universe, frames)
    assert chosen == {"frame": 3, "trajectory_frame": 3}
    assert centroid.positions == pytest.approx(coordinates[3], abs=1e-5)
    assert centroid.universe is not universe
    assert np.array_equal(universe.trajectory.timeseries(order="fac"), before)


def test_centroid_rmsd_is_per_atom() -> None:
    """Crosses at 0.8, 1.0 and 1.3 average to 1.0333, 0.0333 Å from the frame at 1.0."""
    coordinates = np.array([np.asarray(CROSS) * scale for scale in (0.8, 1.0, 1.3)])
    index, value = _find_frame_closest_to_aligned_mean(coordinates)
    assert index == 1
    assert value == pytest.approx(0.1 / 3, abs=1e-5)
