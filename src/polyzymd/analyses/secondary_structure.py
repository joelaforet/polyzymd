"""Secondary structure content from DSSP, written against the observable contract.

The plugin assigns a simplified DSSP class to every protein residue in every
production frame with ``mdtraj.compute_dssp(simplified=True)`` and reports four
``fraction`` observables (the per-frame fraction of residues in helix, in
strand, in coil, and unassigned) plus two ``profile`` observables (the fraction
of the window each residue spends in helix and in strand).

Only ``ss_helix`` and ``ss_strand`` are tested across conditions. The four
fractions sum to one in every frame, so testing all four would put two
dependent tests into the multiple-comparison family and weaken the adjusted
p-values of the two that carry independent information. ``ss_coil`` and
``ss_unassigned`` are still aggregated and reported with their uncertainty.

``unassigned`` counts residues mdtraj returns as ``NA`` because they carry no
usable backbone or an unrecognised residue name. The previous implementation
encoded those residues as coil, which inflated the coil fraction and hid a
broken selection. They now have their own fraction, so a non-zero
``ss_unassigned`` says the selection needs attention.

References
----------
Kabsch, W. & Sander, C. (1983). Dictionary of protein secondary structure:
pattern recognition of hydrogen-bonded and geometrical features. *Biopolymers*,
22(12), 2577-2637. doi:10.1002/bip.360221211

McGibbon, R. T., Beauchamp, K. A., Harrigan, M. P., Klein, C., Swails, J. M.,
Hernandez, C. X., Schwantes, C. R., Wang, L.-P., Lane, T. J. & Pande, V. S.
(2015). MDTraj: a modern open library for the analysis of molecular dynamics
trajectories. *Biophysical Journal*, 109(8), 1528-1532.
doi:10.1016/j.bpj.2015.08.015
"""

from __future__ import annotations

from typing import Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, Field, field_validator

from polyzymd.analyses.base import SlurmResourceHint
from polyzymd.analyses.contract import Observable, iter_frames
from polyzymd.analyses.contract_runner import contract_analysis
from polyzymd.analyses.exceptions import ReplicateError

#: Simplified DSSP characters mdtraj returns, mapped to observable names.
DSSP_CLASSES: dict[str, str] = {"H": "helix", "E": "strand", "C": "coil", "NA": "unassigned"}


class SecondaryStructureSettings(BaseModel):
    """Which protein residues to assign.

    Chain ``A`` is the protein chain under the PolyzyMD chain convention. Set
    ``selection`` when the topology does not preserve chain IDs.
    """

    chain_id: str = Field(default="A", description="Chain letter for the protein chain")
    selection: str | None = Field(
        default=None,
        description="Explicit MDAnalysis protein selection. Overrides chain_id when provided.",
    )

    @field_validator("chain_id", "selection")
    @classmethod
    def _non_blank(cls, value: str | None) -> str | None:
        """Reject a blank chain ID or selection."""
        if value is None:
            return None
        stripped = value.strip()
        if not stripped:
            raise ValueError("secondary_structure chain_id and selection must not be blank")
        return stripped


class SecondaryStructure:
    """Simplified DSSP content and per-residue occupancy of one protein chain."""

    name: ClassVar[str] = "secondary_structure"
    Settings: ClassVar[type[BaseModel]] = SecondaryStructureSettings
    references: ClassVar[tuple[str, ...]] = (
        "Kabsch and Sander 1983, Biopolymers 22:2577, doi:10.1002/bip.360221211",
        "McGibbon et al. 2015, Biophys J 109:1528, doi:10.1016/j.bpj.2015.08.015",
    )
    #: DSSP holds the whole window of protein coordinates in memory.
    slurm_resource_hint: ClassVar[SlurmResourceHint | None] = SlurmResourceHint(mem="16G")

    def compute(
        self, universe: Any, frames: Any, settings: SecondaryStructureSettings
    ) -> Sequence[Observable]:
        """Assign DSSP classes over the production window.

        Parameters
        ----------
        universe : MDAnalysis.Universe
            Universe loaded by the framework.
        frames : FrameSelection
            Production window resolved by the framework.
        settings : SecondaryStructureSettings
            Chain or selection to assign.

        Returns
        -------
        Sequence[Observable]
            Four ``fraction`` observables named ``ss_helix``, ``ss_strand``,
            ``ss_coil`` and ``ss_unassigned``, each a per-frame fraction of the
            selected residues, and two ``profile`` observables named
            ``helix_occupancy`` and ``strand_occupancy`` indexed by residue ID.
            ``ss_helix`` and ``ss_strand`` declare ``higher_is_better=True``
            and are the two that enter the cross-condition tests; ``ss_coil``
            and ``ss_unassigned`` are ``tested=False`` because they are
            determined by the other two.

        Raises
        ------
        ReplicateError
            If the selection matches no atoms or covers partial residues.
        """
        import mdtraj as md

        group = self._select(universe, settings)
        topology, residue_ids = _mdtraj_topology(group)
        positions = np.asarray(
            [group.positions.copy() for _ in iter_frames(universe, frames)], dtype=np.float32
        )
        if positions.size == 0:
            raise ReplicateError("secondary_structure: the production window selected no frames")
        classes = md.compute_dssp(
            md.Trajectory(xyz=positions / 10.0, topology=topology), simplified=True
        )
        fractions = [
            Observable(
                name=f"ss_{label}",
                kind="fraction",
                unit="fraction",
                values=np.asarray((classes == char).mean(axis=1), dtype=np.float64).tolist(),
                higher_is_better=True if label in ("helix", "strand") else None,
                # Coil and unassigned are what the other two are not, so testing
                # all four would add two dependent tests to the correction
                # family. They are still aggregated and reported.
                tested=label in ("helix", "strand"),
            )
            for char, label in DSSP_CLASSES.items()
        ]
        profiles = [
            Observable(
                name=f"{label}_occupancy",
                kind="profile",
                unit="fraction",
                values=np.asarray((classes == char).mean(axis=0), dtype=np.float64).tolist(),
                index=[float(residue_id) for residue_id in residue_ids],
            )
            for char, label in (("H", "helix"), ("E", "strand"))
        ]
        return fractions + profiles

    def _select(self, universe: Any, settings: SecondaryStructureSettings) -> Any:
        """Return the selected atoms, refusing a selection DSSP cannot use."""
        selection = settings.selection or f"protein and chainid {settings.chain_id}"
        try:
            group = universe.select_atoms(selection)
        except AttributeError as exc:
            raise ReplicateError(
                f"secondary_structure: selection {selection!r} failed",
                hint="Topologies without chain IDs need an explicit selection such as 'protein'.",
            ) from exc
        if len(group) == 0:
            raise ReplicateError(
                f"secondary_structure: selection {selection!r} matched no atoms",
                hint="Check plugins.secondary_structure.chain_id or set selection explicitly.",
            )
        if sum(len(residue.atoms) for residue in group.residues) != len(group):
            raise ReplicateError(
                f"secondary_structure: selection {selection!r} covers partial residues",
                hint="DSSP needs whole residues; use 'protein' or 'protein and resid A:B'.",
            )
        return group


def _mdtraj_topology(group: Any) -> tuple[Any, list[int]]:
    """Build an mdtraj topology from an MDAnalysis atom group.

    Returns the topology and the residue IDs in selection order, which index
    the occupancy profiles.
    """
    import mdtraj as md

    topology = md.Topology()
    chain = topology.add_chain()
    residue_ids: list[int] = []
    for residue in group.residues:
        residue_id = int(residue.resid)
        md_residue = topology.add_residue(str(residue.resname).upper(), chain, resSeq=residue_id)
        residue_ids.append(residue_id)
        for atom in residue.atoms:
            topology.add_atom(str(atom.name), _element(md, atom), md_residue)
    return topology, residue_ids


def _element(md: Any, atom: Any) -> Any:
    """Best-effort mdtraj element for one atom, falling back to carbon."""
    symbol = str(getattr(atom, "element", "") or "").strip()
    if not symbol:
        symbol = "".join(char for char in str(atom.name) if char.isalpha())[:1]
    try:
        return md.element.get_by_symbol(symbol.capitalize() or "C")
    except KeyError:
        return md.element.carbon


SecondaryStructureAnalysis = contract_analysis(SecondaryStructure)
