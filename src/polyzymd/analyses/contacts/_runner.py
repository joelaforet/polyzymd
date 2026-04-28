"""Runner-backed contact execution helpers.

This module keeps the trajectory-native contact work behind the runner seam so
the contacts plugin can preserve its legacy result models and side outputs while
PolyzyMD owns replicate, aggregation, and comparison orchestration.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from MDAnalysis.core.groups import AtomGroup, Residue
    from MDAnalysis.core.universe import Universe
    from numpy.typing import NDArray

    from polyzymd.analyses.contacts._results import ContactResult

LOGGER = logging.getLogger(__name__)


class _IdentityResidueGrouping:
    """Group residues by residue name without biochemical coarsening."""

    @property
    def available_groups(self) -> list[str]:
        """Return an empty group list for open-ended residue-name labels."""

        return []

    def classify(self, resname: str) -> str:
        """Return the normalized residue name as the group label.

        Parameters
        ----------
        resname : str
            Residue name to normalize.

        Returns
        -------
        str
            Uppercase residue name, or ``"unknown"`` when blank.
        """

        normalized = str(resname).upper().strip()
        return normalized or "unknown"


class _SecondaryStructureGrouping:
    """Group residues by attached secondary-structure annotations."""

    _DSSP_LABELS = {
        "H": "helix",
        "G": "helix",
        "I": "helix",
        "E": "sheet",
        "B": "sheet",
        "T": "turn",
        "S": "turn",
        "C": "coil",
        "-": "coil",
        "": "coil",
    }

    @property
    def available_groups(self) -> list[str]:
        """Return supported secondary-structure group labels."""

        return ["helix", "sheet", "turn", "coil", "unknown"]

    def classify(self, resname: str) -> str:
        """Classify residue names without structural annotations.

        Parameters
        ----------
        resname : str
            Residue name. It is unused because secondary structure requires
            residue-level annotations.

        Returns
        -------
        str
            ``"unknown"`` because no residue annotation was provided.
        """

        del resname
        return "unknown"

    def classify_residue(self, residue: Any) -> str:
        """Classify one residue from common secondary-structure attributes.

        Parameters
        ----------
        residue : Any
            MDAnalysis residue or residue-like test double.

        Returns
        -------
        str
            Coarse secondary-structure group label.
        """

        raw_label = self._extract_secondary_structure_label(residue)
        if raw_label is None:
            return "unknown"
        label = str(raw_label).strip()
        if not label:
            return "coil"
        return self._DSSP_LABELS.get(label.upper(), label.lower())

    @staticmethod
    def _extract_secondary_structure_label(residue: Any) -> Any | None:
        """Extract a secondary-structure label from common residue attributes."""

        for attr in ("secondary_structure", "secstruct", "ss", "dssp"):
            value = getattr(residue, attr, None)
            if value is not None:
                return value

        atoms = getattr(residue, "atoms", None)
        if atoms is None:
            return None
        for attr in ("secondary_structures", "secstructs", "ss", "dssp"):
            values = getattr(atoms, attr, None)
            if values is None:
                continue
            try:
                if len(values) == 0:
                    continue
                return values[0]
            except TypeError:
                return values
        return None


def build_contact_grouping(grouping: str) -> Any:
    """Build a residue grouping strategy from contacts settings.

    Parameters
    ----------
    grouping : str
        Contacts grouping mode: ``"aa_class"``, ``"secondary_structure"``,
        or ``"none"``.

    Returns
    -------
    Any
        Grouping object with a ``classify`` method.
    """

    from polyzymd.analyses.shared.groupings import ProteinAAClassification

    if grouping == "aa_class":
        return ProteinAAClassification()
    if grouping == "secondary_structure":
        return _SecondaryStructureGrouping()
    if grouping == "none":
        return _IdentityResidueGrouping()
    raise ValueError(f"Unsupported contacts grouping mode: {grouping}")


def identify_polymer_chains(
    query_residues: "AtomGroup",
) -> list[list[tuple[int, str, "Residue"]]]:
    """Identify polymer chains and their segments from an atom group.

    Groups residues by connected fragments and returns structured chain
    information.

    Parameters
    ----------
    query_residues : AtomGroup
        MDAnalysis atom group containing the polymer residues to classify.

    Returns
    -------
    list[list[tuple[int, str, Residue]]]
        List of chains, where each chain is a list of
        ``(chain_idx, resname, residue)`` tuples.
    """

    all_atoms = query_residues.atoms
    fragments = _fragments_or_single(all_atoms, context="contacts polymer chain detection")

    chains: list[list[tuple[int, str, "Residue"]]] = []
    for chain_idx, fragment in enumerate(fragments):
        chain_residues = [(chain_idx, residue.resname, residue) for residue in fragment.residues]
        if chain_residues:
            chains.append(chain_residues)

    return chains


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


def _get_contact_analysis_base_cls() -> type:
    """Return a cached ``AnalysisBase`` subclass for residue contacts.

    MDAnalysis is imported lazily because it is a heavy dependency.
    """

    cached = getattr(_get_contact_analysis_base_cls, "_cls", None)
    if cached is not None:
        return cached

    from MDAnalysis.analysis.base import AnalysisBase

    class _ContactAnalysisBase(AnalysisBase):  # type: ignore[misc]
        """MDAnalysis ``AnalysisBase`` for residue-level contact detection."""

        def __init__(
            self,
            target_atoms: "AtomGroup",
            query_atoms: "AtomGroup",
            target_residue_indices: "NDArray[np.int64]",
            query_residue_indices: "NDArray[np.int64]",
            cutoff: float,
            n_target_residues: int,
            n_query_residues: int,
            **kwargs: Any,
        ):
            super().__init__(target_atoms.universe.trajectory, **kwargs)

            self.target_atoms = target_atoms
            self.query_atoms = query_atoms
            self.target_residue_indices = target_residue_indices
            self.query_residue_indices = query_residue_indices
            self.cutoff = cutoff
            self.n_target_residues = n_target_residues
            self.n_query_residues = n_query_residues

        def _prepare(self) -> None:
            self._contact_matrix = np.zeros(
                (self.n_frames, self.n_target_residues, self.n_query_residues),
                dtype=np.uint8,
            )

        def _single_frame(self) -> None:
            from MDAnalysis.lib.distances import capped_distance

            pairs, _ = capped_distance(
                self.query_atoms.positions,
                self.target_atoms.positions,
                max_cutoff=self.cutoff,
                box=self._ts.dimensions,
                return_distances=True,
            )

            if len(pairs) == 0:
                return

            query_atom_indices = pairs[:, 0]
            target_atom_indices = pairs[:, 1]
            query_res_indices = self.query_residue_indices[query_atom_indices]
            target_res_indices = self.target_residue_indices[target_atom_indices]
            self._contact_matrix[self._frame_index, target_res_indices, query_res_indices] = 1

        def _conclude(self) -> None:
            self.results.contact_matrix = self._contact_matrix

    _get_contact_analysis_base_cls._cls = _ContactAnalysisBase  # type: ignore[attr-defined]
    return _ContactAnalysisBase


class ParallelContactAnalyzer:
    """Optimized contact analyzer using capped-distance neighbour searching.

    Parameters
    ----------
    target_selector : Any
        Selector for target residues, typically protein atoms.
    query_selector : Any
        Selector for query residues, typically polymer atoms.
    cutoff : float, optional
        Contact distance cutoff in Å, by default 4.0.
    grouping : Any | None, optional
        Residue grouping strategy, by default ``ProteinAAClassification()``.
    """

    def __init__(
        self,
        target_selector: Any,
        query_selector: Any,
        cutoff: float = 4.0,
        grouping: Any | None = None,
    ) -> None:
        self.target_selector = target_selector
        self.query_selector = query_selector
        self.cutoff = cutoff
        self.grouping = grouping or build_contact_grouping("aa_class")

    def run(
        self,
        universe: "Universe",
        start: int = 0,
        stop: int | None = None,
        step: int = 1,
        verbose: bool = False,
    ) -> "ContactResult":
        """Run contact analysis over one trajectory window.

        Parameters
        ----------
        universe : Universe
            MDAnalysis universe with loaded trajectory.
        start : int, optional
            Inclusive start frame, by default 0.
        stop : int | None, optional
            Exclusive stop frame, by default ``None``.
        step : int, optional
            Frame stride, by default 1.
        verbose : bool, optional
            Emit progress logging, by default ``False``.

        Returns
        -------
        ContactResult
            Per-replicate contact result in the legacy contacts schema.
        """

        from polyzymd.analyses.contacts._results import (
            ContactResult,
            PolymerSegmentContacts,
            ResidueContactData,
            compress_contact_array,
        )

        target_result = self.target_selector.select(universe)
        query_result = self.query_selector.select(universe)

        target_atoms = target_result.atoms
        query_atoms = query_result.atoms
        target_residues = target_result.residues
        query_residues = query_result.residues

        if verbose:
            LOGGER.info(
                "Analyzing contacts: %d query residues (%d atoms) -> %d target residues (%d atoms)",
                len(query_residues),
                len(query_atoms),
                len(target_residues),
                len(target_atoms),
            )
            LOGGER.info("Cutoff: %.1f A", self.cutoff)

        target_atom_to_res = self._build_atom_to_residue_map(target_atoms, target_residues)
        query_atom_to_res = self._build_atom_to_residue_map(query_atoms, query_residues)

        try:
            timestep_ps = universe.trajectory.dt
        except (AttributeError, TypeError, ValueError):
            timestep_ps = 1.0

        analysis = _get_contact_analysis_base_cls()(
            target_atoms=target_atoms,
            query_atoms=query_atoms,
            target_residue_indices=target_atom_to_res,
            query_residue_indices=query_atom_to_res,
            cutoff=self.cutoff,
            n_target_residues=len(target_residues),
            n_query_residues=len(query_residues),
        )
        analysis.run(start=start, stop=stop, step=step, verbose=verbose)

        contact_matrix = analysis.results.contact_matrix
        n_frames = contact_matrix.shape[0]

        if verbose:
            LOGGER.info("Processing %d frames of contact data", n_frames)

        polymer_chains = identify_polymer_chains(query_residues)
        residue_lookup = {
            self._unique_residue_key(residue): i for i, residue in enumerate(query_residues)
        }
        res_to_chain_seg: dict[int, tuple[int, int, str, int]] = {}
        for chain_idx, chain in enumerate(polymer_chains):
            for seg_idx, (_, _, residue) in enumerate(chain):
                query_idx = residue_lookup.get(self._unique_residue_key(residue))
                if query_idx is not None:
                    res_to_chain_seg[query_idx] = (
                        chain_idx,
                        seg_idx,
                        residue.resname,
                        residue.resid,
                    )

        residue_contacts: list[ResidueContactData] = []
        for target_idx, target_residue in enumerate(target_residues):
            group = self._classify_target_residue(target_residue)
            target_contacts = contact_matrix[:, target_idx, :]
            segment_contacts: list[PolymerSegmentContacts] = []

            for query_idx in range(target_contacts.shape[1]):
                contact_array = target_contacts[:, query_idx].astype(bool)
                if not np.any(contact_array):
                    continue

                if query_idx in res_to_chain_seg:
                    chain_idx, _, resname, resid = res_to_chain_seg[query_idx]
                else:
                    query_residue = query_residues[query_idx]
                    chain_idx = 0
                    resname = query_residue.resname
                    resid = query_residue.resid

                events = compress_contact_array(contact_array)
                if events:
                    segment_contacts.append(
                        PolymerSegmentContacts(
                            polymer_resname=resname,
                            polymer_resid=resid,
                            polymer_chain_idx=chain_idx,
                            events=events,
                        )
                    )

            residue_contacts.append(
                ResidueContactData(
                    protein_resid=target_residue.resid,
                    protein_resname=target_residue.resname,
                    protein_group=group,
                    segment_contacts=segment_contacts,
                )
            )

        result = ContactResult(
            residue_contacts=residue_contacts,
            n_frames=n_frames,
            timestep_ps=timestep_ps * step,
            criteria_label=f"any_atom_{self.cutoff:.1f}A",
            criteria_cutoff=self.cutoff,
            start_frame=start,
            metadata={
                "target_selector": self.target_selector.label,
                "query_selector": self.query_selector.label,
                "n_polymer_chains": len(polymer_chains),
                "n_polymer_segments": sum(len(chain) for chain in polymer_chains),
                "optimized": True,
                "algorithm": "capped_distance",
            },
        )
        result.compute_per_residue_statistics()

        if verbose:
            LOGGER.info(
                "Analysis complete: %d/%d residues contacted (%.1f%%)",
                result.n_contacted_residues,
                result.n_protein_residues,
                result.coverage_fraction() * 100.0,
            )

        return result

    def _classify_target_residue(self, residue: Any) -> str:
        """Classify one target residue with the configured grouping strategy.

        Parameters
        ----------
        residue : Any
            MDAnalysis residue or residue-like object.

        Returns
        -------
        str
            Group label for the residue.
        """

        classify_residue = getattr(self.grouping, "classify_residue", None)
        if callable(classify_residue):
            return str(classify_residue(residue))
        return str(self.grouping.classify(residue.resname))

    @staticmethod
    def _unique_residue_key(residue: Any) -> int | tuple[str, int | Any, str]:
        """Build a stable identity key for one residue."""

        residue_ix = getattr(residue, "ix", None)
        if residue_ix is not None:
            try:
                return int(residue_ix)
            except (TypeError, ValueError):
                pass

        chain_id = getattr(residue, "chainID", None) or getattr(residue, "segid", "")
        return (str(chain_id), getattr(residue, "resid", None), getattr(residue, "resname", ""))

    @staticmethod
    def _build_atom_to_residue_map(
        atoms: Any,
        residues: Any,
    ) -> "NDArray[np.int64]":
        """Map each atom to its parent-residue index."""

        residue_lookup = {
            ParallelContactAnalyzer._unique_residue_key(residue): index
            for index, residue in enumerate(residues)
        }
        atom_to_residue = np.zeros(len(atoms), dtype=np.int64)
        for atom_index, atom in enumerate(atoms):
            residue_key = ParallelContactAnalyzer._unique_residue_key(atom.residue)
            atom_to_residue[atom_index] = residue_lookup[residue_key]
        return atom_to_residue


@dataclass
class ContactsReplicateRunner:
    """Execute contacts analysis for one replicate through the runner seam."""

    universe: Any
    target_selector: Any
    query_selector: Any
    cutoff: float
    grouping: Any | None = None
    analyzer_cls: type[ParallelContactAnalyzer] = ParallelContactAnalyzer
    verbose: bool = False
    results: Any = field(default=None, init=False)

    def run(self, start: int, stop: int, step: int = 1) -> ContactsReplicateRunner:
        """Execute contacts over the provided trajectory window."""

        analyzer = self.analyzer_cls(
            target_selector=self.target_selector,
            query_selector=self.query_selector,
            cutoff=self.cutoff,
            grouping=self.grouping,
        )
        self.results = analyzer.run(
            self.universe,
            start=start,
            stop=stop,
            step=step,
            verbose=self.verbose,
        )
        return self
