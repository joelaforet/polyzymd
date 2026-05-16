"""Shared helpers for GROMACS smoke tests.

This module centralizes lightweight fixtures used by analysis plugin smoke tests
that exercise real TrajectoryLoader/GROMACS path resolution while mocking only
MDAnalysis-heavy behavior.
"""

from __future__ import annotations

from pathlib import Path
from types import ModuleType
from typing import Sequence
from unittest.mock import MagicMock

import numpy as np

from polyzymd.analyses.base import Condition


def make_gromacs_config(tmp_path: Path) -> MagicMock:
    """Create a minimal simulation config mock that dispatches to GROMACS.

    Parameters
    ----------
    tmp_path : Path
        Temporary root directory used for replicate run folders.

    Returns
    -------
    MagicMock
        Config mock with attributes required by ``TrajectoryLoader`` and
        ``GromacsEngine`` layout resolution.
    """
    config = MagicMock()
    config.engine = "gromacs"

    config.gromacs = MagicMock()
    config.gromacs.gmx_binary = "gmx"

    config.generate_system_name = MagicMock(return_value="test_system")
    config.get_working_directory = MagicMock(
        side_effect=lambda replicate: tmp_path / f"run_{replicate}"
    )

    config.output = MagicMock()
    config.output.effective_scratch_directory = tmp_path
    return config


def create_gromacs_layout(run_dir: Path) -> None:
    """Create a minimal GROMACS run directory layout.

    Parameters
    ----------
    run_dir : Path
        Replicate directory such as ``run_1``.
    """
    gromacs_dir = run_dir / "gromacs"
    gromacs_dir.mkdir(parents=True, exist_ok=True)
    (gromacs_dir / "solvated_system.pdb").write_bytes(b"\n")
    (gromacs_dir / "prod.xtc").write_bytes(b"\n")


class MockTrajectory:
    """Simple trajectory object that emulates minimal MDAnalysis behavior.

    Parameters
    ----------
    n_frames : int, optional
        Number of frames available in the trajectory, by default 200.
    dt_ps : float, optional
        Timestep in picoseconds exposed via ``dt``, by default 10.0.
    """

    def __init__(self, n_frames: int = 200, dt_ps: float = 10.0) -> None:
        self._n_frames = n_frames
        self._dt_ps = dt_ps
        self.time = 0.0

    @property
    def dt(self) -> float:
        """Return the trajectory timestep in picoseconds."""
        return self._dt_ps

    def __len__(self) -> int:
        """Return the number of frames."""
        return self._n_frames

    def __getitem__(self, item: int | slice):
        """Return one frame or a list of frames for slices."""
        if isinstance(item, slice):
            start, stop, step = item.indices(self._n_frames)
            return [self[idx] for idx in range(start, stop, step)]

        frame = int(item)
        if frame < 0:
            frame += self._n_frames
        self.time = frame * self._dt_ps

        ts = MagicMock()
        ts.frame = frame
        ts.time = self.time
        return ts

    def __iter__(self):
        """Yield all trajectory frames in order."""
        for frame in range(self._n_frames):
            yield self[frame]


def make_mock_universe(
    n_frames: int = 200,
    n_atoms: int = 5,
    dt_ps: float = 10.0,
    n_residues: int | None = None,
    chain_ids: Sequence[str] | None = None,
    atom_names: Sequence[str] | None = None,
    residue_names: Sequence[str] | None = None,
    residue_ids: Sequence[int] | None = None,
    seed: int = 7,
) -> MagicMock:
    """Build a mock MDAnalysis Universe for smoke tests.

    Parameters
    ----------
    n_frames : int, optional
        Number of trajectory frames, by default 200.
    n_atoms : int, optional
        Number of atoms in the default atom group, by default 5.
    dt_ps : float, optional
        Trajectory timestep in picoseconds, by default 10.0.
    n_residues : int | None, optional
        Number of residues to create. If ``None``, uses ``n_atoms``.
    chain_ids : Sequence[str] | None, optional
        Chain IDs to assign to residues. Values are cycled if shorter than
        ``n_residues``. Defaults to ``("A",)``.
    atom_names : Sequence[str] | None, optional
        Atom names for the default atom group. Values are cycled if shorter
        than ``n_atoms``. Defaults to ``("CA",)``.
    residue_names : Sequence[str] | None, optional
        Residue names to assign. Values are cycled if shorter than
        ``n_residues``. Defaults to ``("ALA",)``.
    residue_ids : Sequence[int] | None, optional
        Residue IDs to assign. If omitted, IDs are ``1..n_residues``.
    seed : int, optional
        Random seed used for deterministic coordinate generation.

    Returns
    -------
    MagicMock
        Universe-like object with ``trajectory`` and ``select_atoms`` behavior
        suitable for smoke tests.
    """
    residue_count = n_atoms if n_residues is None else n_residues
    resolved_chain_ids = tuple(chain_ids) if chain_ids else ("A",)
    resolved_atom_names = tuple(atom_names) if atom_names else ("CA",)
    resolved_residue_names = tuple(residue_names) if residue_names else ("ALA",)
    resolved_residue_ids = tuple(residue_ids) if residue_ids else tuple(range(1, residue_count + 1))

    universe = MagicMock()
    universe.trajectory = MockTrajectory(n_frames=n_frames, dt_ps=dt_ps)

    atoms = MagicMock()
    atoms.__len__ = MagicMock(return_value=n_atoms)
    atoms.n_atoms = n_atoms

    rng = np.random.default_rng(seed)
    atoms.positions = rng.random((n_atoms, 3)).astype(np.float32)
    atoms.indices = np.arange(n_atoms)
    atoms.names = np.array(
        [resolved_atom_names[i % len(resolved_atom_names)] for i in range(n_atoms)]
    )

    residues: list[MagicMock] = []
    for idx in range(residue_count):
        residue = MagicMock()
        residue.ix = idx
        residue.resid = resolved_residue_ids[idx % len(resolved_residue_ids)]
        residue.resname = resolved_residue_names[idx % len(resolved_residue_names)]
        residue.segid = resolved_chain_ids[idx % len(resolved_chain_ids)]
        residue.chainID = resolved_chain_ids[idx % len(resolved_chain_ids)]
        residue.segment = MagicMock()
        residue.segment.segid = residue.segid
        residue.atoms = MagicMock()
        residue.atoms.indices = np.asarray([idx % n_atoms], dtype=np.int64)
        residues.append(residue)

    atoms.residues = residues
    universe.select_atoms = MagicMock(return_value=atoms)
    return universe


def install_fake_mdanalysis() -> ModuleType:
    """Create a fake MDAnalysis module for ``sys.modules`` patching.

    Returns
    -------
    ModuleType
        Fake module exposing a mock ``Universe`` constructor.
    """
    fake_mda = ModuleType("MDAnalysis")
    fake_mda.Universe = MagicMock()
    return fake_mda


def make_condition(
    tmp_path: Path,
    config: MagicMock,
    label: str = "GROMACS Smoke",
    replicates: tuple[int, ...] = (1, 2),
) -> Condition:
    """Build a test condition for GROMACS smoke tests.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory used for config path placement.
    config : MagicMock
        Simulation config mock associated with the condition.
    label : str, optional
        Human-readable condition label, by default ``"GROMACS Smoke"``.
    replicates : tuple[int, ...], optional
        Replicate IDs included in the condition, by default ``(1, 2)``.

    Returns
    -------
    Condition
        Constructed condition object for analysis contexts.
    """
    return Condition(
        label=label,
        config_path=tmp_path / "config.yaml",
        replicates=replicates,
        sim_config=config,
    )
