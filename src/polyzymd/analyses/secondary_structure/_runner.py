"""Runner-backed DSSP computation for secondary-structure analysis."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class SecondaryStructureRunnerResult:
    """Trajectory-native secondary-structure payload.

    The plugin converts this payload into the public
    :class:`SecondaryStructureResult` schema during ``summarize_replicate()``.
    """

    selection_string: str
    n_frames: int
    n_residues: int
    residue_ids: list[int]
    residue_names: list[str]
    ss_matrix: np.ndarray
    persistence_coil: list[float]
    persistence_helix: list[float]
    persistence_strand: list[float]
    overall_helix_fraction: float
    overall_strand_fraction: float
    overall_coil_fraction: float


class SecondaryStructureReplicateRunner:
    """Compute DSSP assignments for one replicate with mdtraj.

    Parameters
    ----------
    trajectory_files : Sequence[Path]
        Trajectory files for the replicate.
    topology_file : Path
        Topology file used to load the trajectory.
    chain_id : str
        Protein chain letter following the PolyzyMD chain convention.
    chain_index_func : Callable[[str], int]
        Function converting chain letters to mdtraj chain indices.
    encode_dssp_func : Callable[[np.ndarray], np.ndarray]
        Function converting DSSP character assignments to integer codes.
    """

    def __init__(
        self,
        *,
        trajectory_files: Sequence[Path],
        topology_file: Path,
        chain_id: str,
        chain_index_func: Callable[[str], int],
        encode_dssp_func: Callable[[np.ndarray], np.ndarray],
    ) -> None:
        self.trajectory_files = tuple(Path(path) for path in trajectory_files)
        self.topology_file = Path(topology_file)
        self.chain_id = chain_id
        self._chain_index_func = chain_index_func
        self._encode_dssp_func = encode_dssp_func
        self.results: SecondaryStructureRunnerResult | None = None

    def run(self, start: int, stop: int, step: int) -> SecondaryStructureReplicateRunner:
        """Execute DSSP computation over the requested frame window.

        Parameters
        ----------
        start : int
            Inclusive start frame.
        stop : int
            Exclusive stop frame.
        step : int
            Frame stride.

        Returns
        -------
        SecondaryStructureReplicateRunner
            The executed runner with ``results`` populated.
        """

        import mdtraj as md

        traj = md.load([str(path) for path in self.trajectory_files], top=str(self.topology_file))
        chain_idx = self._chain_index_func(self.chain_id)
        protein_indices = traj.topology.select(f"chainid {chain_idx}")
        if len(protein_indices) == 0:
            raise ValueError(
                f"Chain '{self.chain_id}' (index {chain_idx}) not found in topology. "
                f"Available chains: {[chain.index for chain in traj.topology.chains]}. "
                "Check your Settings.chain_id configuration."
            )

        protein_traj = traj.atom_slice(protein_indices)
        protein_traj = protein_traj[start:stop:step]
        dssp_raw = md.compute_dssp(protein_traj, simplified=True)
        ss_matrix = self._encode_dssp_func(dssp_raw)

        residue_ids = [residue.resSeq for residue in protein_traj.topology.residues]
        residue_names = [residue.name.upper() for residue in protein_traj.topology.residues]
        self.results = self._summarize_matrix(
            ss_matrix=ss_matrix,
            residue_ids=residue_ids,
            residue_names=residue_names,
            selection_string=f"chainid {chain_idx} (chain {self.chain_id})",
        )
        return self

    @staticmethod
    def _summarize_matrix(
        *,
        ss_matrix: np.ndarray,
        residue_ids: list[int],
        residue_names: list[str],
        selection_string: str,
    ) -> SecondaryStructureRunnerResult:
        """Compute persistence and content fractions from encoded DSSP data.

        Parameters
        ----------
        ss_matrix : np.ndarray
            Integer-encoded secondary-structure matrix.
        residue_ids : list[int]
            Residue IDs matching matrix columns.
        residue_names : list[str]
            Residue names matching matrix columns.
        selection_string : str
            Human-readable chain selection used for the analysis.

        Returns
        -------
        SecondaryStructureRunnerResult
            Summary payload for plugin-level serialization.
        """

        total_entries = ss_matrix.size
        n_frames, n_residues = ss_matrix.shape
        return SecondaryStructureRunnerResult(
            selection_string=selection_string,
            n_frames=n_frames,
            n_residues=n_residues,
            residue_ids=residue_ids,
            residue_names=residue_names,
            ss_matrix=ss_matrix,
            persistence_coil=(ss_matrix == 0).mean(axis=0).tolist(),
            persistence_helix=(ss_matrix == 1).mean(axis=0).tolist(),
            persistence_strand=(ss_matrix == 2).mean(axis=0).tolist(),
            overall_helix_fraction=float(np.sum(ss_matrix == 1)) / total_entries,
            overall_strand_fraction=float(np.sum(ss_matrix == 2)) / total_entries,
            overall_coil_fraction=float(np.sum(ss_matrix == 0)) / total_entries,
        )
