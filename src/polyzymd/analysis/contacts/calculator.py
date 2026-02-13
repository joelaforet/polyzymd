"""Contact analyzer for computing polymer-protein contacts.

This module provides the main ContactAnalyzer class that computes
contacts between polymer residues and protein residues frame-by-frame.

Key features:
- Frame-by-frame computation for proper autocorrelation analysis
- Compressed event storage for space efficiency
- Support for multiple polymer chains
- Configurable contact criteria (Strategy pattern)
- Progress reporting for long trajectories

Examples
--------
>>> from MDAnalysis import Universe
>>> from polyzymd.analysis.contacts import ContactAnalyzer, AnyAtomWithinCutoff
>>> from polyzymd.analysis.common.selectors import ProteinResidues, PolymerChains
>>> from polyzymd.analysis.common.groupings import ProteinAAClassification
>>>
>>> # Load trajectory
>>> u = Universe("topology.pdb", "trajectory.dcd")
>>>
>>> # Create analyzer
>>> analyzer = ContactAnalyzer(
...     target_selector=ProteinResidues(),
...     query_selector=PolymerChains(),
...     criteria=AnyAtomWithinCutoff(cutoff=4.0),
...     grouping=ProteinAAClassification(),
... )
>>>
>>> # Run analysis
>>> result = analyzer.run(u, start=1000, verbose=True)
>>> print(f"Coverage: {result.coverage_fraction():.1%}")
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Callable

import numpy as np
from numpy.typing import NDArray

from polyzymd.analysis.common.groupings import ProteinAAClassification, ResidueGrouping
from polyzymd.analysis.common.selectors import PolymerChains, ProteinResidues
from polyzymd.analysis.common.selectors.base import MolecularSelector
from polyzymd.analysis.contacts.criteria import AnyAtomWithinCutoff, ContactCriteria
from polyzymd.analysis.contacts.results import (
    ContactResult,
    ResidueContactData,
    PolymerSegmentContacts,
    compress_contact_array,
)

if TYPE_CHECKING:
    from MDAnalysis.core.groups import AtomGroup, Residue
    from MDAnalysis.core.universe import Universe

logger = logging.getLogger(__name__)


class ContactAnalyzer:
    """Analyzer for computing molecular contacts frame-by-frame.

    This class computes contacts between two molecular groups (typically
    polymer and protein) across a trajectory. Results are stored in
    compressed event format for space efficiency while preserving
    frame-by-frame information for autocorrelation analysis.

    Parameters
    ----------
    target_selector : MolecularSelector
        Selector for target residues (typically protein).
        Contacts are computed TO these residues.
    query_selector : MolecularSelector
        Selector for query residues (typically polymer).
        Contacts are computed FROM these residues.
    criteria : ContactCriteria
        Definition of what constitutes a "contact".
        Default: AnyAtomWithinCutoff(cutoff=4.0)
    grouping : ResidueGrouping, optional
        Classification scheme for target residues.
        Default: ProteinAAClassification()
    per_residue : bool
        If True, compute contacts for each target residue separately.
        If False, treat all targets as one group. Default True.

    Examples
    --------
    >>> # Basic usage with defaults
    >>> analyzer = ContactAnalyzer(
    ...     target_selector=ProteinResidues(),
    ...     query_selector=PolymerChains(),
    ... )
    >>> result = analyzer.run(universe)
    >>>
    >>> # Custom criteria and grouping
    >>> analyzer = ContactAnalyzer(
    ...     target_selector=ProteinResiduesNearReference("resid 77 133 156", cutoff=10),
    ...     query_selector=PolymerResiduesByType(["SBM"]),
    ...     criteria=AnyAtomToCOM(cutoff=5.0),
    ...     grouping=ProteinAAClassification(),
    ... )
    """

    def __init__(
        self,
        target_selector: MolecularSelector,
        query_selector: MolecularSelector,
        criteria: ContactCriteria | None = None,
        grouping: ResidueGrouping | None = None,
        per_residue: bool = True,
    ):
        self.target_selector = target_selector
        self.query_selector = query_selector
        self.criteria = criteria or AnyAtomWithinCutoff(cutoff=4.0)
        self.grouping = grouping or ProteinAAClassification()
        self.per_residue = per_residue

    def run(
        self,
        universe: "Universe",
        start: int = 0,
        stop: int | None = None,
        step: int = 1,
        verbose: bool = False,
        progress_callback: Callable[[int, int], None] | None = None,
    ) -> ContactResult:
        """Run contact analysis on a trajectory.

        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe with loaded trajectory
        start : int
            First frame to analyze (0-indexed). Default 0.
        stop : int, optional
            Last frame to analyze (exclusive). Default None (all frames).
        step : int
            Frame stride. Default 1 (every frame).
        verbose : bool
            Print progress information. Default False.
        progress_callback : callable, optional
            Function called with (current_frame, total_frames) for progress updates.

        Returns
        -------
        ContactResult
            Complete contact analysis results
        """
        # Get selections
        target_result = self.target_selector.select(universe)
        query_result = self.query_selector.select(universe)

        target_residues = target_result.residues
        query_residues = query_result.residues

        if verbose:
            logger.info(
                f"Analyzing contacts: {len(query_residues)} query residues -> "
                f"{len(target_residues)} target residues"
            )
            logger.info(f"Criteria: {self.criteria}")

        # Get trajectory info
        trajectory = universe.trajectory
        n_total_frames = len(trajectory)
        stop = stop or n_total_frames

        # Validate frame range
        if start < 0 or start >= n_total_frames:
            raise ValueError(f"start={start} out of range [0, {n_total_frames})")
        if stop <= start:
            raise ValueError(f"stop={stop} must be > start={start}")

        # Get timestep (in ps)
        try:
            timestep_ps = trajectory.dt
        except (AttributeError, ValueError):
            timestep_ps = 1.0  # Default if not available
            logger.warning("Could not determine timestep, using 1.0 ps")

        # Identify polymer chains and their segments
        polymer_chains = self._identify_polymer_chains(query_residues)

        # Initialize contact arrays for each (target_residue, polymer_segment) pair
        # Shape: (n_target_residues, n_polymer_segments, n_frames)
        frame_indices = list(range(start, stop, step))
        n_frames = len(frame_indices)

        contact_data = self._initialize_contact_data(target_residues, polymer_chains, n_frames)

        # Frame-by-frame analysis
        if verbose:
            logger.info(f"Analyzing {n_frames} frames...")

        for frame_idx, ts_idx in enumerate(frame_indices):
            trajectory[ts_idx]
            box = trajectory.ts.dimensions

            # Compute contacts for this frame
            self._compute_frame_contacts(
                target_residues,
                polymer_chains,
                contact_data,
                frame_idx,
                box,
            )

            if progress_callback:
                progress_callback(frame_idx + 1, n_frames)

            if verbose and (frame_idx + 1) % max(1, n_frames // 10) == 0:
                logger.info(f"  Frame {frame_idx + 1}/{n_frames}")

        # Convert to compressed events
        result = self._build_result(
            target_residues,
            polymer_chains,
            contact_data,
            n_frames,
            timestep_ps * step,  # Adjusted for step
            start,
        )

        if verbose:
            logger.info(
                f"Analysis complete: {result.n_contacted_residues}/{result.n_protein_residues} "
                f"residues contacted ({result.coverage_fraction():.1%})"
            )

        return result

    def _identify_polymer_chains(
        self, query_residues: "AtomGroup"
    ) -> list[list[tuple[int, str, "Residue"]]]:
        """Identify polymer chains and their segments.

        Returns list of chains, where each chain is a list of
        (chain_idx, resname, residue) tuples.
        """
        # Group residues by fragment (connected component)
        all_atoms = query_residues.atoms
        fragments = all_atoms.fragments if all_atoms.fragments else [all_atoms]

        chains = []
        for chain_idx, frag in enumerate(fragments):
            chain_residues = []
            for res in frag.residues:
                chain_residues.append((chain_idx, res.resname, res))
            if chain_residues:
                chains.append(chain_residues)

        return chains

    def _initialize_contact_data(
        self,
        target_residues: "AtomGroup",
        polymer_chains: list[list[tuple[int, str, "Residue"]]],
        n_frames: int,
    ) -> dict:
        """Initialize data structure for storing contact arrays.

        Returns dict with structure:
        {
            target_resid: {
                (chain_idx, segment_idx): np.zeros(n_frames, dtype=bool),
                ...
            },
            ...
        }
        """
        data = {}
        for res in target_residues:
            data[res.resid] = {}
            for chain_idx, chain in enumerate(polymer_chains):
                for seg_idx, (_, _, _) in enumerate(chain):
                    data[res.resid][(chain_idx, seg_idx)] = np.zeros(n_frames, dtype=bool)
        return data

    def _compute_frame_contacts(
        self,
        target_residues: "AtomGroup",
        polymer_chains: list[list[tuple[int, str, "Residue"]]],
        contact_data: dict,
        frame_idx: int,
        box: NDArray[np.float64] | None,
    ) -> None:
        """Compute contacts for a single frame."""
        for target_res in target_residues:
            target_atoms = target_res.atoms

            for chain_idx, chain in enumerate(polymer_chains):
                for seg_idx, (_, _, polymer_res) in enumerate(chain):
                    polymer_atoms = polymer_res.atoms

                    result = self.criteria.check_contact(
                        query_atoms=polymer_atoms,
                        target_atoms=target_atoms,
                        box=box,
                    )

                    if result.is_contact:
                        contact_data[target_res.resid][(chain_idx, seg_idx)][frame_idx] = True

    def _build_result(
        self,
        target_residues: "AtomGroup",
        polymer_chains: list[list[tuple[int, str, "Residue"]]],
        contact_data: dict,
        n_frames: int,
        timestep_ps: float,
        start_frame: int,
    ) -> ContactResult:
        """Convert contact arrays to compressed ContactResult."""
        residue_contacts = []

        for target_res in target_residues:
            # Classify the target residue
            group = self.grouping.classify(target_res.resname)

            # Build segment contacts
            segment_contacts = []
            for chain_idx, chain in enumerate(polymer_chains):
                for seg_idx, (_, resname, polymer_res) in enumerate(chain):
                    contact_array = contact_data[target_res.resid][(chain_idx, seg_idx)]

                    # Compress to events
                    events = compress_contact_array(contact_array)

                    if events:  # Only include if there are contacts
                        segment_contacts.append(
                            PolymerSegmentContacts(
                                polymer_resname=resname,
                                polymer_resid=polymer_res.resid,
                                polymer_chain_idx=chain_idx,
                                events=events,
                            )
                        )

            residue_contacts.append(
                ResidueContactData(
                    protein_resid=target_res.resid,
                    protein_resname=target_res.resname,
                    protein_group=group,
                    segment_contacts=segment_contacts,
                )
            )

        result = ContactResult(
            residue_contacts=residue_contacts,
            n_frames=n_frames,
            timestep_ps=timestep_ps,
            criteria_label=self.criteria.label,
            criteria_cutoff=self.criteria.cutoff,
            start_frame=start_frame,
            metadata={
                "target_selector": self.target_selector.label,
                "query_selector": self.query_selector.label,
                "n_polymer_chains": len(polymer_chains),
                "n_polymer_segments": sum(len(c) for c in polymer_chains),
            },
        )

        # Compute per-residue statistical inefficiency for proper uncertainty quantification
        # This follows LiveCoMS best practices (Grossfield et al. 2018)
        logger.debug("Computing per-residue statistical inefficiency...")
        result.compute_per_residue_statistics()

        return result

    def validate(self, universe: "Universe") -> dict:
        """Validate selectors and criteria against a universe.

        Returns diagnostic information about whether the analysis
        can be performed successfully.

        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe to validate against

        Returns
        -------
        dict
            Validation results with keys:
            - valid: bool
            - target_info: dict
            - query_info: dict
            - errors: list[str]
            - warnings: list[str]
        """
        errors = []
        warnings = []

        # Validate target selector
        target_info = self.target_selector.validate(universe)
        if not target_info.get("valid", False):
            errors.append(f"Target selector invalid: {target_info.get('error', 'unknown')}")

        # Validate query selector
        query_info = self.query_selector.validate(universe)
        if not query_info.get("valid", False):
            errors.append(f"Query selector invalid: {query_info.get('error', 'unknown')}")

        # Check trajectory
        if len(universe.trajectory) < 10:
            warnings.append(
                f"Short trajectory ({len(universe.trajectory)} frames). "
                "Results may have limited statistical significance."
            )

        return {
            "valid": len(errors) == 0,
            "target_info": target_info,
            "query_info": query_info,
            "errors": errors,
            "warnings": warnings,
        }


def run_contact_analysis(
    topology: str,
    trajectory: str,
    cutoff: float = 4.0,
    start: int = 0,
    stop: int | None = None,
    step: int = 1,
    output_path: str | None = None,
    verbose: bool = True,
) -> ContactResult:
    """Convenience function for running contact analysis.

    Parameters
    ----------
    topology : str
        Path to topology file (PDB, GRO, etc.)
    trajectory : str
        Path to trajectory file (DCD, XTC, etc.)
    cutoff : float
        Contact distance cutoff in Angstroms. Default 4.0.
    start : int
        First frame to analyze. Default 0.
    stop : int, optional
        Last frame to analyze. Default None (all frames).
    step : int
        Frame stride. Default 1.
    output_path : str, optional
        Path to save JSON results. Default None (don't save).
    verbose : bool
        Print progress. Default True.

    Returns
    -------
    ContactResult
        Analysis results
    """
    import MDAnalysis as mda

    # Load universe
    u = mda.Universe(topology, trajectory)

    if verbose:
        logger.info(f"Loaded: {topology}")
        logger.info(f"Trajectory: {trajectory} ({len(u.trajectory)} frames)")

    # Create analyzer with defaults
    analyzer = ContactAnalyzer(
        target_selector=ProteinResidues(),
        query_selector=PolymerChains(),
        criteria=AnyAtomWithinCutoff(cutoff=cutoff),
        grouping=ProteinAAClassification(),
    )

    # Run analysis
    result = analyzer.run(u, start=start, stop=stop, step=step, verbose=verbose)

    # Save if requested
    if output_path:
        result.save(output_path)
        if verbose:
            logger.info(f"Results saved to: {output_path}")

    return result
