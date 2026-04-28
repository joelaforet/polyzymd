"""Runner-backed trajectory logic for polymer bridging.

This module owns the per-trajectory MDAnalysis iteration for oligomer bridging.
The public plugin keeps lifecycle orchestration, aggregation, comparison, and
plotting in ``polymer_bridging.__init__``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

from polyzymd.analyses.shared.groupings import ProteinAAClassification

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from polyzymd.analyses.base import Condition
    from polyzymd.analyses.polymer_bridging import (
        FrameContactObservation,
        ResiduePairDistances,
    )


LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class PolymerBridgingRunnerResult:
    """Trajectory-native polymer bridging observations for one replicate.

    Parameters
    ----------
    observations : list[FrameContactObservation]
        One observation per contacting fragment per selected frame.
    n_frames : int
        Number of selected trajectory frames inspected.
    timestep_ps : float
        Effective time spacing between selected frames in picoseconds.
    """

    observations: list[FrameContactObservation]
    n_frames: int
    timestep_ps: float


@dataclass
class PolymerBridgingReplicateRunner:
    """Execute polymer bridging observation assembly for one replicate."""

    universe: Any
    protein_selection: str
    polymer_selection: str
    cutoff: float
    min_ca_distance_angstrom: float = 0.0
    results: PolymerBridgingRunnerResult | None = field(default=None, init=False)

    def run(self, start: int, stop: int, step: int = 1) -> PolymerBridgingReplicateRunner:
        """Run observation assembly over the selected trajectory window.

        Parameters
        ----------
        start : int
            Inclusive start frame.
        stop : int
            Exclusive stop frame.
        step : int, optional
            Frame stride, by default 1.

        Returns
        -------
        PolymerBridgingReplicateRunner
            This runner with ``results`` populated.
        """

        observations = compute_observations_for_universe(
            self.universe,
            protein_selection=self.protein_selection,
            polymer_selection=self.polymer_selection,
            cutoff=self.cutoff,
            min_ca_distance_angstrom=self.min_ca_distance_angstrom,
            start=start,
            stop=stop,
            step=step,
        )
        self.results = PolymerBridgingRunnerResult(
            observations=observations,
            n_frames=len(range(start, stop, step)),
            timestep_ps=_trajectory_timestep_ps(self.universe) * step,
        )
        return self


def compute_frame_contacts(
    condition: "Condition",
    replicate: int,
    *,
    protein_selection: str,
    polymer_selection: str,
    cutoff: float,
    equilibration: str,
    min_ca_distance_angstrom: float = 0.0,
) -> tuple[list[FrameContactObservation], int, float]:
    """Compute per-fragment, per-frame contact observations from a trajectory.

    Parameters
    ----------
    condition : Condition
        Simulation condition providing the simulation configuration.
    replicate : int
        Replicate index.
    protein_selection : str
        MDAnalysis selection string for protein atoms.
    polymer_selection : str
        MDAnalysis selection string for polymer atoms.
    cutoff : float
        Atom-level contact cutoff in Angstroms.
    equilibration : str
        Equilibration time string used to resolve the start frame.
    min_ca_distance_angstrom : float, optional
        Minimum frame-wise CA distance for multisite eligibility, by default
        0.0.

    Returns
    -------
    tuple[list[FrameContactObservation], int, float]
        Observations, selected frame count, and effective timestep in
        picoseconds.
    """

    from polyzymd.analyses.shared.loader import TrajectoryLoader
    from polyzymd.analyses.shared.window import resolve_replicate_trajectory_window

    loader = TrajectoryLoader(condition.sim_config)
    universe = loader.load_universe(replicate, cache=False)

    n_frames_total = len(universe.trajectory)
    if n_frames_total == 0:
        return [], 0, 0.0

    protein = universe.select_atoms(protein_selection)
    polymer = universe.select_atoms(polymer_selection)
    if len(protein) == 0 or len(polymer) == 0:
        return [], 0, 0.0

    window = resolve_replicate_trajectory_window(
        loader=loader,
        replicate=replicate,
        equilibration=equilibration,
        n_frames_total=n_frames_total,
    )
    observations = compute_observations_for_groups(
        universe=universe,
        protein=protein,
        polymer=polymer,
        cutoff=cutoff,
        min_ca_distance_angstrom=min_ca_distance_angstrom,
        start=window.start,
        stop=window.stop,
        step=window.step,
    )
    return observations, window.n_frames_selected, window.timestep_ps * window.step


def compute_observations_for_universe(
    universe: Any,
    *,
    protein_selection: str,
    polymer_selection: str,
    cutoff: float,
    min_ca_distance_angstrom: float,
    start: int,
    stop: int,
    step: int = 1,
) -> list[FrameContactObservation]:
    """Assemble bridging observations for selected frames in one Universe.

    Parameters
    ----------
    universe : Any
        MDAnalysis Universe containing topology and trajectory.
    protein_selection : str
        Protein atom selection string.
    polymer_selection : str
        Polymer atom selection string.
    cutoff : float
        Atom-level contact cutoff in Angstroms.
    min_ca_distance_angstrom : float
        CA distance threshold used for validation.
    start : int
        Inclusive start frame.
    stop : int
        Exclusive stop frame.
    step : int, optional
        Frame stride, by default 1.

    Returns
    -------
    list[FrameContactObservation]
        One observation per contacting polymer fragment per frame.
    """

    protein = universe.select_atoms(protein_selection)
    polymer = universe.select_atoms(polymer_selection)
    if len(protein) == 0 or len(polymer) == 0:
        return []
    return compute_observations_for_groups(
        universe=universe,
        protein=protein,
        polymer=polymer,
        cutoff=cutoff,
        min_ca_distance_angstrom=min_ca_distance_angstrom,
        start=start,
        stop=stop,
        step=step,
    )


def compute_observations_for_groups(
    *,
    universe: Any,
    protein: Any,
    polymer: Any,
    cutoff: float,
    min_ca_distance_angstrom: float,
    start: int,
    stop: int,
    step: int = 1,
) -> list[FrameContactObservation]:
    """Assemble bridging observations for pre-selected atom groups.

    Parameters
    ----------
    universe : Any
        MDAnalysis Universe that owns the trajectory.
    protein : Any
        Selected protein atom group.
    polymer : Any
        Selected polymer atom group.
    cutoff : float
        Atom-level contact cutoff in Angstroms.
    min_ca_distance_angstrom : float
        CA distance threshold used for validation.
    start : int
        Inclusive start frame.
    stop : int
        Exclusive stop frame.
    step : int, optional
        Frame stride, by default 1.

    Returns
    -------
    list[FrameContactObservation]
        One observation per contacting polymer fragment per frame.
    """

    from MDAnalysis.lib.distances import capped_distance

    from polyzymd.analyses.polymer_bridging import PolymerBridgingObservation

    protein_atom_to_resid = np.array([int(atom.resid) for atom in protein.atoms], dtype=np.int64)
    protein_atom_to_resname = np.array([str(atom.resname) for atom in protein.atoms], dtype=object)
    protein_grouping = ProteinAAClassification()
    ca_atom_index_by_resid = _ca_atom_index_by_resid(protein)
    if min_ca_distance_angstrom > 0.0 and not ca_atom_index_by_resid:
        raise ValueError(
            "No CA atoms were found in the selected protein topology while "
            "min_ca_distance_angstrom > 0; fix topology files to include CA labels"
        )

    observations: list[FrameContactObservation] = []
    fragments = _fragments_or_single(polymer, context="polymer bridging fragment detection")
    for ts in universe.trajectory[start:stop:step]:
        box = ts.dimensions
        all_positions = ts.positions
        for fragment in fragments:
            pairs, atom_distances = capped_distance(
                fragment.positions,
                protein.positions,
                max_cutoff=cutoff,
                box=box,
                return_distances=True,
            )
            if len(pairs) == 0:
                continue

            protein_atom_indices = np.asarray(pairs[:, 1], dtype=np.int64)
            polymer_atom_indices = np.asarray(pairs[:, 0], dtype=np.int64)
            residue_ids = {int(protein_atom_to_resid[idx]) for idx in protein_atom_indices}
            protein_resnames = {
                int(resid): str(resname)
                for resid, resname in zip(
                    protein_atom_to_resid[protein_atom_indices],
                    protein_atom_to_resname[protein_atom_indices],
                    strict=False,
                )
            }
            protein_groups = {
                resid: protein_grouping.classify(resname)
                for resid, resname in protein_resnames.items()
            }
            contacting_polymer_resids = {
                int(fragment.atoms[int(local_idx)].resid) for local_idx in polymer_atom_indices
            }
            polymer_resnames = {int(res.resid): str(res.resname) for res in fragment.residues}
            pair_min_distances = _compute_pair_min_distances(
                fragment,
                protein,
                polymer_atom_indices,
                protein_atom_indices,
                protein_atom_to_resid,
                atom_distances,
            )
            observations.append(
                PolymerBridgingObservation(
                    protein_residues=residue_ids,
                    protein_resnames=protein_resnames,
                    protein_groups=protein_groups,
                    contacting_polymer_resids=contacting_polymer_resids,
                    polymer_resnames=polymer_resnames,
                    fragment_signature=tuple(str(res.resname) for res in fragment.residues),
                    ca_distances=_observation_ca_distances(
                        residue_ids,
                        all_positions,
                        ca_atom_index_by_resid,
                        box,
                    ),
                    pair_min_distances=pair_min_distances,
                )
            )
    return observations


def _fragments_or_single(atom_group: Any, *, context: str) -> list[Any]:
    """Return MDAnalysis fragments with a no-bond topology fallback.

    Parameters
    ----------
    atom_group : Any
        MDAnalysis atom group to split into connected fragments.
    context : str
        Diagnostic context for warning messages.

    Returns
    -------
    list[Any]
        Fragment atom groups, or the input atom group as a single fragment
        when bond topology is unavailable.
    """

    from MDAnalysis.exceptions import NoDataError

    try:
        fragments = atom_group.fragments
    except NoDataError:
        LOGGER.warning(
            "%s: topology has no bond information; treating the selected polymer as one fragment",
            context,
        )
        return [atom_group]

    return list(fragments) if fragments else [atom_group]


def _ca_atom_index_by_resid(protein: Any) -> dict[int, int]:
    """Map protein residue IDs to their global CA atom indices."""

    ca_atom_index_by_resid: dict[int, int] = {}
    for residue in protein.residues:
        ca_atoms = residue.atoms.select_atoms("name CA")
        if len(ca_atoms) == 0:
            continue
        ca_atom_index_by_resid[int(residue.resid)] = int(ca_atoms.indices[0])
    return ca_atom_index_by_resid


def _trajectory_timestep_ps(universe: Any) -> float:
    """Return the trajectory timestep in picoseconds with a safe fallback."""

    try:
        timestep_ps = float(universe.trajectory.dt)
    except (AttributeError, TypeError, ValueError):
        timestep_ps = 1.0
    return timestep_ps if timestep_ps > 0.0 else 1.0


def _observation_ca_distances(
    residues: set[int],
    all_positions: "NDArray[np.float64]",
    ca_atom_index_by_resid: dict[int, int],
    box: Any,
) -> "ResiduePairDistances":
    """Compute frame-wise CA distances for one contacted residue set.

    Parameters
    ----------
    residues : set[int]
        Contacted protein residue IDs.
    all_positions : NDArray[np.float64]
        Frame atom positions indexed by global atom index.
    ca_atom_index_by_resid : dict[int, int]
        Mapping from protein residue ID to global CA atom index.
    box : Any
        Periodic box dimensions for minimum-image distance calculations.

    Returns
    -------
    ResiduePairDistances
        Pairwise CA distances keyed by sorted residue ID pairs.
    """

    from MDAnalysis.lib.distances import distance_array

    distances: dict[tuple[int, int], float] = {}
    residues_sorted = sorted(residues)
    for i, resid_i in enumerate(residues_sorted):
        for resid_j in residues_sorted[i + 1 :]:
            key = (resid_i, resid_j)
            idx_i = ca_atom_index_by_resid.get(resid_i)
            idx_j = ca_atom_index_by_resid.get(resid_j)
            if idx_i is None or idx_j is None:
                continue
            pos_i = all_positions[idx_i]
            pos_j = all_positions[idx_j]
            distance = distance_array(
                np.asarray(pos_i, dtype=float).reshape(1, 3),
                np.asarray(pos_j, dtype=float).reshape(1, 3),
                box=box,
            )[0, 0]
            distances[key] = float(distance)
    return distances


def _compute_pair_min_distances(
    fragment: Any,
    protein: Any,
    polymer_atom_indices: "NDArray[np.int64]",
    protein_atom_indices: "NDArray[np.int64]",
    protein_atom_to_resid: "NDArray[np.int64]",
    atom_distances: "NDArray[np.float64]",
) -> dict[tuple[int, int], float]:
    """Compute minimum atom distance for each polymer/protein residue pair.

    Parameters
    ----------
    fragment : Any
        Contacting polymer fragment.
    protein : Any
        Protein atom group. It is retained for API symmetry with the contact
        pair indices.
    polymer_atom_indices : NDArray[np.int64]
        Local polymer atom indices from the PBC-aware contact search.
    protein_atom_indices : NDArray[np.int64]
        Local protein atom indices from the PBC-aware contact search.
    protein_atom_to_resid : NDArray[np.int64]
        Protein atom index to residue ID mapping.
    atom_distances : NDArray[np.float64]
        PBC-aware atom distances returned by ``capped_distance``.

    Returns
    -------
    dict[tuple[int, int], float]
        Minimum PBC-aware distance for each polymer/protein residue pair.
    """

    del protein
    distances: dict[tuple[int, int], float] = {}
    for local_poly_idx, protein_idx, atom_distance in zip(
        polymer_atom_indices, protein_atom_indices, atom_distances, strict=False
    ):
        polymer_atom = fragment.atoms[int(local_poly_idx)]
        polymer_resid = int(polymer_atom.resid)
        protein_resid = int(protein_atom_to_resid[int(protein_idx)])
        key = (polymer_resid, protein_resid)
        distance = float(atom_distance)
        if key not in distances or distance < distances[key]:
            distances[key] = distance
    return distances
