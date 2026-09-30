"""Writers of small on-disk simulations and replicate values for the analysis tests.

:func:`write_simulation_config`, :func:`write_openmm_replicate` and
:func:`write_openmm_frames` write a minimal OpenMM simulation that
:class:`~polyzymd.analyses.study.Study` can load; :func:`replicate_values`
builds a :class:`~polyzymd.analyses.timeseries.ReplicateValues` from given
numbers.
"""

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from polyzymd.config.schema import SimulationConfig

# ---------------------------------------------------------------------------
# On-disk simulation outputs
# ---------------------------------------------------------------------------

#: Four unit-mass atoms on a unit cross. The radius of gyration is the scale.
CROSS = ((1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0))


def write_simulation_config(directory: Path, *, scratch: Path, name: str = "toy") -> Path:
    """Write a minimal OpenMM ``config.yaml`` whose run directories go in ``scratch``.

    Parameters
    ----------
    directory : Path
        Folder to write ``config.yaml`` into. It is created when missing.
    scratch : Path
        Scratch directory the config puts its run directories in.
    name : str, optional
        Simulation name.

    Returns
    -------
    Path
        Path of the written ``config.yaml``.
    """
    import yaml

    directory.mkdir(parents=True, exist_ok=True)
    data = {
        "name": name,
        "engine": "openmm",
        "enzyme": {"name": "TestEnzyme", "pdb_path": "test.pdb"},
        "thermodynamics": {"temperature": 300.0},
        "simulation_phases": {
            "equilibration_stages": [
                {"name": "eq1", "duration": 0.1, "temperature": 300.0, "ensemble": "NVT"}
            ],
            "production": {
                "ensemble": "NPT",
                "duration": 1.0,
                "samples": 10,
                "checkpoint_interval": 60.0,
            },
        },
        "output": {"projects_directory": str(directory), "scratch_directory": str(scratch)},
    }
    path = directory / "config.yaml"
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    return path


def write_openmm_replicate(
    config_path: Path, replicate: int, scales: Sequence[float], *, dt_ps: float = 100.0
) -> Path:
    """Write an OpenMM run directory with one DCD production segment.

    Frame ``k`` is the four-atom cross of :data:`CROSS` scaled by
    ``scales[k]``, so its mass-weighted radius of gyration is ``scales[k]``
    and its time is ``k * dt_ps`` ps.

    Parameters
    ----------
    config_path : Path
        Simulation ``config.yaml`` written by :func:`write_simulation_config`.
    replicate : int
        Replicate number.
    scales : sequence of float
        Radius of gyration of each frame.
    dt_ps : float, optional
        Time between frames in ps.

    Returns
    -------
    Path
        The run directory.
    """
    mda = importlib.import_module("MDAnalysis")

    run_dir = SimulationConfig.from_yaml(config_path).get_working_directory(replicate)
    segment = run_dir / "production_0"
    segment.mkdir(parents=True, exist_ok=True)
    universe = mda.Universe.empty(4, n_residues=1, atom_resindex=[0] * 4, trajectory=True)
    universe.add_TopologyAttr("names", ["C1", "C2", "C3", "C4"])
    universe.add_TopologyAttr("resnames", ["MOL"])
    universe.add_TopologyAttr("masses", [1.0] * 4)
    cross = np.asarray(CROSS, dtype=np.float32)
    universe.atoms.positions = cross
    universe.atoms.write(str(run_dir / "solvated_system.pdb"))
    path = segment / "production_0_trajectory.dcd"
    with mda.Writer(str(path), n_atoms=4, dt=dt_ps, istart=0, nsavc=1) as writer:
        for scale in scales:
            universe.atoms.positions = cross * float(scale)
            writer.write(universe.atoms)
    return run_dir


def write_openmm_frames(
    config_path: Path,
    replicate: int,
    coordinates: Any,
    atom_resindex: Sequence[int],
    *,
    resids: Sequence[int] | None = None,
    names: Sequence[str] | None = None,
    resnames: Sequence[str] | None = None,
    elements: Sequence[str] | None = None,
    chain_ids: Sequence[str] | None = None,
    dimensions: Sequence[float] | None = None,
    bonds: Sequence[tuple[int, int]] | None = None,
    dt_ps: float = 100.0,
) -> Path:
    """Write an OpenMM run directory whose DCD holds ``coordinates`` frame by frame.

    Atom ``i`` is named ``names[i]``, by default ``C<i>``, and belongs to
    residue ``atom_resindex[i]``, whose residue ID is
    ``resids[atom_resindex[i]]``, by default that index plus one, and whose
    name is ``resnames[atom_resindex[i]]``, by default ``ALA``. ``elements``
    and ``chain_ids``, one per atom, fill the PDB's element and chain ID
    columns when given. ``dimensions``, the box as ``[a, b, c, alpha, beta,
    gamma]`` in Å and degrees, is written on every frame and in the PDB when
    given; without it the DCD holds no box. ``bonds``, pairs of atom indices,
    are written as CONECT records, so the loaded topology has bonds and
    bonded fragments; without it the loaded topology has no bonds. The
    topology holds the first frame.

    Returns
    -------
    Path
        The run directory.
    """
    mda = importlib.import_module("MDAnalysis")

    coordinates = np.asarray(coordinates, dtype=np.float32)
    n_atoms, n_residues = coordinates.shape[1], max(atom_resindex) + 1
    run_dir = SimulationConfig.from_yaml(config_path).get_working_directory(replicate)
    segment = run_dir / "production_0"
    segment.mkdir(parents=True, exist_ok=True)
    universe = mda.Universe.empty(
        n_atoms, n_residues=n_residues, atom_resindex=list(atom_resindex), trajectory=True
    )
    universe.add_TopologyAttr("names", list(names or [f"C{i}" for i in range(n_atoms)]))
    universe.add_TopologyAttr("resnames", list(resnames or ["ALA"] * n_residues))
    universe.add_TopologyAttr("resids", list(resids or range(1, n_residues + 1)))
    universe.add_TopologyAttr("masses", [1.0] * n_atoms)
    if elements is not None:
        universe.add_TopologyAttr("elements", list(elements))
    if chain_ids is not None:
        universe.add_TopologyAttr("chainIDs", list(chain_ids))
    if dimensions is not None:
        universe.dimensions = np.asarray(dimensions, dtype=np.float32)
    if bonds is not None:
        universe.add_TopologyAttr("bonds", [tuple(pair) for pair in bonds])
    universe.atoms.positions = coordinates[0]
    universe.atoms.write(str(run_dir / "solvated_system.pdb"), bonds="all")
    path = segment / "production_0_trajectory.dcd"
    with mda.Writer(str(path), n_atoms=n_atoms, dt=dt_ps, istart=0, nsavc=1) as writer:
        for frame in coordinates:
            universe.atoms.positions = frame
            writer.write(universe.atoms)
    return run_dir


def replicate_values(
    per_condition: dict[str, list[float]],
    how: Any = "mean",
    series_bounds: Any = (None, None),
    **reduce_kwargs: Any,
) -> Any:
    """Build :class:`~polyzymd.analyses.timeseries.ReplicateValues` without a trajectory.

    Each replicate is a constant series of 20 frames, so ``"mean"`` and
    ``"fraction"`` reduce it to the given value, its statistical inefficiency
    is 1 and its effective sample size is 20.

    Parameters
    ----------
    per_condition : dict of str to list of float
        Replicate values of each condition, control first.
    how : str, optional
        Reduction passed to ``Timeseries.reduce``.

    Returns
    -------
    ReplicateValues
        Values named ``<how>_rg`` in unit ``A``, equilibration ``10ns``.
    """
    from types import SimpleNamespace

    from polyzymd.analyses.timeseries import ReplicateSeries, Timeseries

    class _Study:
        control = next(iter(per_condition))

        def __getitem__(self, label: str) -> Any:
            return SimpleNamespace(equilibration="10ns", config_hash="hash")

    frames = np.arange(20)
    series = {
        label: [
            ReplicateSeries(label, index, np.full(20, value), frames, frames * 0.1, Path())
            for index, value in enumerate(values, start=1)
        ]
        for label, values in per_condition.items()
    }
    return Timeseries("rg", "A", _Study(), series, Path(), tuple(series_bounds)).reduce(
        how, **reduce_kwargs
    )
