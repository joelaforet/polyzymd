"""Solvent-accessible surface area, written against the observable contract.

The plugin measures one target selection inside one or more contexts. A context
is the set of atoms allowed to block the surface, so ``protein`` inside
``protein`` is the isolated protein surface and ``protein`` inside
``protein or resname SBM EGM`` is the protein surface with the polymer present.
Each context reports two observables: the total area of the target per frame as
a ``mean_of_timeseries`` in square angstrom, and the per-residue mean relative
area as a ``profile``, the residue area divided by the maximum accessible area
of that residue type.

Areas come from the Shrake-Rupley algorithm as implemented by MDTraj, with a
probe radius of 0.14 nm and 960 sphere points by default. MDTraj works in
square nanometre, so every area is multiplied by 100 to reach square angstrom.
Solvent and ions are excluded because they are never part of a target or a
context selection.

References
----------
Shrake, A. & Rupley, J. A. (1973). Environment and exposure to solvent of
protein atoms. Lysozyme and insulin. *Journal of Molecular Biology*, 79(2),
351-371. doi:10.1016/0022-2836(73)90011-9

Tien, M. Z., Meyer, A. G., Sydykova, D. K., Spielman, S. J. & Wilke, C. O.
(2013). Maximum allowed solvent accessibilities of residues in proteins.
*PLoS ONE*, 8(11), e80635. doi:10.1371/journal.pone.0080635

McGibbon, R. T., Beauchamp, K. A., Harrigan, M. P., Klein, C., Swails, J. M.,
Hernandez, C. X., Schwantes, C. R., Wang, L.-P., Lane, T. J. & Pande, V. S.
(2015). MDTraj: a modern open library for the analysis of molecular dynamics
trajectories. *Biophysical Journal*, 109(8), 1528-1532.
doi:10.1016/j.bpj.2015.08.015
"""

from __future__ import annotations

import tempfile
import warnings
from typing import Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, Field, model_validator

from polyzymd.analyses.base import SlurmResourceHint
from polyzymd.analyses.contract import Observable, iter_frames
from polyzymd.analyses.contract_runner import contract_analysis
from polyzymd.analyses.exceptions import ReplicateError, SelectionError
from polyzymd.analyses.shared.aa_classification import get_max_asa

#: MDTraj reports square nanometre; every area is reported in square angstrom.
NM2_TO_A2 = 100.0


class SASARun(BaseModel):
    """One context in which the area of a target selection is measured."""

    label: str = Field(min_length=1, description="Context name, for example 'protein_isolated'")
    target_selection: str = Field(min_length=1, description="Atoms whose area is reported")
    context_selection: str | None = Field(
        default=None,
        description="Atoms allowed to block the surface; defaults to the target selection",
    )
    stride: int = Field(default=1, ge=1, description="Deprecated and ignored since v1.3")

    @model_validator(mode="after")
    def _reject_per_run_stride(self) -> SASARun:
        """Warn that a per-context stride no longer selects frames."""
        if self.stride != 1:
            warnings.warn(
                f"sasa run {self.label!r} sets stride={self.stride}; the framework now resolves "
                "one frame window for every observable, so the setting is ignored. Remove it and "
                "set the window with --eq-time instead. It will be rejected in v1.4.",
                DeprecationWarning,
                stacklevel=2,
            )
        return self


class SASASettings(BaseModel):
    """Settings for the solvent-accessible surface area analysis."""

    runs: list[SASARun] = Field(min_length=1, description="Contexts to measure")
    probe_radius_nm: float = Field(default=0.14, gt=0.0, description="Probe radius in nm")
    n_sphere_points: int = Field(default=960, ge=100, description="Sphere points per atom")
    chunk_size: int = Field(default=100, ge=1, description="Frames per MDTraj call")

    @model_validator(mode="after")
    def _require_unique_labels(self) -> SASASettings:
        """Reject two contexts that would produce the same observable name."""
        labels = [run.label for run in self.runs]
        if len(set(labels)) != len(labels):
            raise ValueError(f"sasa run labels must be unique, got {labels}")
        return self


class SASA:
    """Solvent-accessible surface area of a selection in one or more contexts."""

    name: ClassVar[str] = "sasa"
    Settings: ClassVar[type[BaseModel]] = SASASettings
    execution_cost_hint: ClassVar[str] = "high"
    slurm_resource_hint: ClassVar[SlurmResourceHint] = SlurmResourceHint(mem="8G", time="02:00:00")
    references: ClassVar[tuple[str, ...]] = (
        "Shrake & Rupley 1973, J Mol Biol 79:351, doi:10.1016/0022-2836(73)90011-9",
        "Tien et al. 2013, PLoS ONE 8:e80635, doi:10.1371/journal.pone.0080635",
        "McGibbon et al. 2015, Biophys J 109:1528, doi:10.1016/j.bpj.2015.08.015",
    )

    def compute(self, universe: Any, frames: Any, settings: SASASettings) -> Sequence[Observable]:
        """Measure the surface area of every configured context.

        Parameters
        ----------
        universe : MDAnalysis.Universe
            Universe loaded by the framework.
        frames : FrameSelection
            Production window resolved by the framework.
        settings : SASASettings
            Contexts to measure and the Shrake-Rupley settings.

        Returns
        -------
        Sequence[Observable]
            Two observables per context: ``sasa_<label>`` in square angstrom and
            ``relative_sasa_<label>`` over the residues of the target.

        Raises
        ------
        SelectionError
            If a target or context selection matches no atoms.
        ReplicateError
            If a target is not a subset of its context, or if a target residue
            has no maximum accessible area in the Tien et al. 2013 table.
        """
        observables: list[Observable] = []
        for run in settings.runs:
            target = universe.select_atoms(run.target_selection)
            context = universe.select_atoms(run.context_selection or run.target_selection)
            _check_selections(run, target, context)
            target_local, groups = _index_target(target, context)
            totals, residue_areas = _sasa_series(
                universe, frames, context, target_local, groups, settings
            )
            observables.append(
                Observable(
                    name=f"sasa_{run.label}",
                    kind="mean_of_timeseries",
                    unit="A^2",
                    values=totals,
                    higher_is_better=False,
                )
            )
            observables.append(
                Observable(
                    name=f"relative_sasa_{run.label}",
                    kind="profile",
                    unit="fraction",
                    index=[float(resid) for resid, _, _ in groups],
                    values=np.mean(residue_areas, axis=0) / _max_asa(run.label, groups),
                    higher_is_better=False,
                )
            )
        return observables


def _check_selections(run: SASARun, target: Any, context: Any) -> None:
    """Reject an empty selection or a target that its context does not contain."""
    empty = [name for name, group in (("target", target), ("context", context)) if len(group) == 0]
    if empty:
        raise SelectionError(
            f"sasa run {run.label!r}: {' and '.join(empty)} selection matched no atoms "
            f"(target_selection={run.target_selection!r}, "
            f"context_selection={run.context_selection or run.target_selection!r})"
        )
    missing = np.setdiff1d(target.indices, context.indices)
    if missing.size:
        raise ReplicateError(
            f"sasa run {run.label!r}: {missing.size} target atoms are outside the context "
            "selection, so their area would be computed without their own neighbours. "
            "Widen context_selection so it contains the target."
        )


def _index_target(target: Any, context: Any) -> tuple[np.ndarray, list[tuple[int, str, list[int]]]]:
    """Locate the target inside its context and group it into residues.

    Returns the context-local index of every target atom, and one
    ``(resid, resname, context-local indices)`` entry per target residue in
    trajectory order. Grouping keys on the topology residue index, so two
    residues that share a chain, a residue ID and a residue name stay separate.
    """
    local = {int(index): position for position, index in enumerate(context.indices.tolist())}
    target_local = [local[int(index)] for index in target.indices.tolist()]
    grouped: dict[int, list[int]] = {}
    for resindex, position in zip(target.resindices.tolist(), target_local):
        grouped.setdefault(int(resindex), []).append(position)
    residues = {int(residue.resindex): residue for residue in target.residues}
    return np.asarray(target_local), [
        (int(residues[resindex].resid), str(residues[resindex].resname), positions)
        for resindex, positions in grouped.items()
    ]


def _max_asa(label: str, groups: Sequence[tuple[int, str, list[int]]]) -> np.ndarray:
    """Maximum accessible area of every target residue, in square angstrom."""
    maxima = [get_max_asa(resname) for _, resname, _ in groups]
    unknown = sorted({resname for (_, resname, _), value in zip(groups, maxima) if value is None})
    if unknown:
        raise ReplicateError(
            f"sasa run {label!r}: residues {unknown} have no maximum accessible area in the "
            "Tien et al. 2013 table, so their relative area is undefined. Restrict "
            "target_selection to standard amino acids."
        )
    return np.asarray(maxima, dtype=np.float64)


def _sasa_series(
    universe: Any,
    frames: Any,
    context: Any,
    target_local: np.ndarray,
    groups: Sequence[tuple[int, str, list[int]]],
    settings: SASASettings,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-frame total target area and per-frame per-residue area, both in A^2.

    Coordinates are buffered ``chunk_size`` frames at a time so a long window
    does not hold the whole trajectory in memory. Chunking does not change the
    numbers, because Shrake-Rupley treats every frame independently.
    """
    import mdtraj as md

    topology = _mdtraj_topology(context)
    residue_indices = [positions for _, _, positions in groups]
    totals: list[np.ndarray] = []
    residues: list[np.ndarray] = []
    buffer: list[np.ndarray] = []

    def flush() -> None:
        if not buffer:
            return
        trajectory = md.Trajectory(
            xyz=np.asarray(buffer, dtype=np.float32) / 10.0, topology=topology
        )
        atom_nm2 = np.asarray(
            md.shrake_rupley(
                trajectory,
                mode="atom",
                probe_radius=settings.probe_radius_nm,
                n_sphere_points=settings.n_sphere_points,
            ),
            dtype=np.float64,
        )
        totals.append(np.sum(atom_nm2[:, target_local] * NM2_TO_A2, axis=1))
        residues.append(
            np.column_stack(
                [np.sum(atom_nm2[:, indices], axis=1) * NM2_TO_A2 for indices in residue_indices]
            )
        )
        buffer.clear()

    for _ in iter_frames(universe, frames):
        buffer.append(np.asarray(context.positions, dtype=np.float32).copy())
        if len(buffer) >= settings.chunk_size:
            flush()
    flush()
    if not totals:
        raise ReplicateError("sasa: the resolved frame window contains no frames")
    return np.concatenate(totals), np.concatenate(residues, axis=0)


def _mdtraj_topology(atoms: Any) -> Any:
    """MDTraj topology for an MDAnalysis atom group, through a temporary PDB file."""
    import mdtraj as md

    with tempfile.NamedTemporaryFile(suffix=".pdb") as handle:
        atoms.write(handle.name)
        return md.load(handle.name).topology


SASAAnalysis = contract_analysis(SASA)
