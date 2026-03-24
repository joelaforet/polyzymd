"""Polymer chain and residue selectors.

This module provides selectors for polymer chains and residues:

- PolymerChains: Select all polymer chains
- PolymerResiduesByType: Select polymer residues by residue name (monomer type)

For systems built with PolyzyMD, use chain_id="C" (the default) to select
polymers based on the PolyzyMD chain convention:
- Chain A: Protein/Enzyme
- Chain B: Substrate/Ligand
- Chain C: Polymers
- Chain D+: Solvent (water, ions, co-solvents)

Examples
--------
>>> # Select polymer chain C (PolyzyMD default)
>>> selector = PolymerChains()
>>> result = selector.select(universe)
>>>
>>> # Select by residue names (for non-PolyzyMD systems)
>>> selector = PolymerChains(chain_id=None, residue_names=["SBM", "EGP"])
>>>
>>> # Select specific polymer types within chain C
>>> selector = PolymerResiduesByType(residue_names=["SBM"])
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from polyzymd.analysis.common.selectors.base import MolecularSelector, SelectionResult

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe


# PolyzyMD chain convention
POLYZYMD_POLYMER_CHAIN_ID = "C"

# Common polymer residue names used in PolyzyMD simulations
# Users can extend this or provide their own lists
DEFAULT_POLYMER_RESNAMES = [
    # Sulfobetaine methacrylate variants
    "SBM",
    "SBMA",
    "SB",
    # Ethylene glycol methacrylate variants
    "EGM",
    "EGMA",
    "EGP",
    "EGPMA",
    "OEGMA",
    # Phosphorylcholine
    "MPC",
    "PC",
    # Generic polymer names
    "MON",  # Monomer
    "POL",  # Polymer
    "PLY",
]


class PolymerChains(MolecularSelector):
    """Select polymer chains from the system.

    For PolyzyMD-built systems, polymers are assigned to Chain C by convention.
    This selector uses chain ID selection by default, which is more reliable
    than residue name matching.

    Parameters
    ----------
    chain_id : str, optional
        Chain ID for polymer selection. Default "C" (PolyzyMD convention).
        Set to None to use residue_names instead.
    residue_names : list[str], optional
        Residue names that identify polymer residues.
        Only used when chain_id is None, or as a filter within the chain.
        Default uses common PolyzyMD polymer names.
    chain_indices : list[int], optional
        If provided, select only these polymer chain indices (0-indexed)
        from within the selected atoms. Useful when analyzing specific
        polymer chains in multi-chain systems.
    segids : list[str], optional
        If provided, select only polymers with these segment IDs.

    Notes
    -----
    The PolyzyMD chain convention is:
    - Chain A: Protein/Enzyme
    - Chain B: Substrate/Ligand
    - Chain C: Polymers
    - Chain D+: Solvent (water, ions, co-solvents)

    For systems not built with PolyzyMD, set chain_id=None and provide
    residue_names explicitly.

    Examples
    --------
    >>> # PolyzyMD system (recommended)
    >>> selector = PolymerChains()  # Uses chain C
    >>>
    >>> # Non-PolyzyMD system
    >>> selector = PolymerChains(chain_id=None, residue_names=["SBM", "EGM"])
    >>>
    >>> # PolyzyMD system with specific polymer types
    >>> selector = PolymerChains(residue_names=["SBM"])  # SBM in chain C only
    """

    def __init__(
        self,
        chain_id: str | None = POLYZYMD_POLYMER_CHAIN_ID,
        residue_names: list[str] | None = None,
        chain_indices: list[int] | None = None,
        segids: list[str] | None = None,
    ):
        self.chain_id = chain_id
        self.residue_names = residue_names or DEFAULT_POLYMER_RESNAMES
        self.chain_indices = chain_indices
        self.segids = segids

    def select(self, universe: "Universe") -> SelectionResult:
        """Select polymer atoms/residues."""
        selection_parts = []

        # Primary selection: by chain ID or residue names
        if self.chain_id is not None:
            selection_parts.append(f"chainID {self.chain_id}")
            # If residue_names also provided, use as additional filter
            if self.residue_names and self.residue_names != DEFAULT_POLYMER_RESNAMES:
                resname_str = " ".join(self.residue_names)
                selection_parts.append(f"resname {resname_str}")
        else:
            # Fall back to residue name selection
            resname_str = " ".join(self.residue_names)
            selection_parts.append(f"resname {resname_str}")

        if self.segids:
            segid_str = " ".join(self.segids)
            selection_parts.append(f"segid {segid_str}")

        selection = " and ".join(f"({part})" for part in selection_parts)
        atoms = universe.select_atoms(selection)

        if len(atoms) == 0:
            if self.chain_id is not None:
                raise ValueError(
                    f"No polymer atoms found in chain '{self.chain_id}'. "
                    "If this is not a PolyzyMD system, use chain_id=None and "
                    "provide residue_names explicitly."
                )
            else:
                raise ValueError(
                    f"No polymer atoms found with residue names: {self.residue_names}. "
                    "Check that polymer residue names match your topology."
                )

        # If chain_indices specified, filter to those chains
        if self.chain_indices is not None:
            # Group residues by fragment (connected component)
            fragments = atoms.fragments
            if not fragments:
                # Fallback: use residue groups
                fragments = [atoms]

            selected_atoms = None
            for idx in self.chain_indices:
                if idx >= len(fragments):
                    raise ValueError(
                        f"Chain index {idx} out of range. Found {len(fragments)} polymer chains."
                    )
                if selected_atoms is None:
                    selected_atoms = fragments[idx]
                else:
                    selected_atoms = selected_atoms | fragments[idx]

            atoms = selected_atoms

        return SelectionResult(
            atoms=atoms,
            residues=atoms.residues,
            label=self.label,
            metadata={
                "chain_id": self.chain_id,
                "residue_names": self.residue_names,
                "chain_indices": self.chain_indices,
                "segids": self.segids,
                "n_chains": len(atoms.fragments) if atoms.fragments else 1,
            },
        )

    @property
    def label(self) -> str:
        if self.chain_indices:
            return f"polymer_chains_{'-'.join(str(i) for i in self.chain_indices)}"
        if self.chain_id:
            return f"polymer_chain{self.chain_id}"
        return "polymer"


class PolymerResiduesByType(MolecularSelector):
    """Select polymer residues by monomer type (residue name).

    This selector groups polymer residues by their residue names,
    allowing analysis of specific monomer types within copolymers.

    Parameters
    ----------
    residue_names : list[str]
        Residue names to select (e.g., ["SBM", "EGP"] for SBMA-EGMA copolymer)
    exclude : bool, optional
        If True, select polymer residues NOT matching these names. Default False.

    Examples
    --------
    >>> # Select SBMA monomers only
    >>> selector = PolymerResiduesByType(residue_names=["SBM", "SBMA"])
    >>>
    >>> # Select non-SBMA monomers
    >>> selector = PolymerResiduesByType(residue_names=["SBM", "SBMA"], exclude=True)
    """

    def __init__(
        self,
        residue_names: list[str],
        exclude: bool = False,
    ):
        if not residue_names:
            raise ValueError("residue_names cannot be empty")

        self.residue_names = residue_names
        self.exclude = exclude

    def select(self, universe: "Universe") -> SelectionResult:
        """Select polymer residues by type."""
        resname_str = " ".join(self.residue_names)

        if self.exclude:
            # Select all polymers EXCEPT these types
            # First get all polymer residues
            all_polymer_str = " ".join(DEFAULT_POLYMER_RESNAMES)
            selection = f"resname {all_polymer_str} and not resname {resname_str}"
        else:
            selection = f"resname {resname_str}"

        atoms = universe.select_atoms(selection)

        if len(atoms) == 0:
            mode = "excluding" if self.exclude else "with"
            raise ValueError(f"No polymer residues found {mode} names: {self.residue_names}")

        return SelectionResult(
            atoms=atoms,
            residues=atoms.residues,
            label=self.label,
            metadata={
                "residue_names": self.residue_names,
                "exclude": self.exclude,
            },
        )

    @property
    def label(self) -> str:
        prefix = "not_" if self.exclude else ""
        return f"{prefix}{'_'.join(self.residue_names)}"


class PolymerSegments(MolecularSelector):
    """Select individual segments (residues) within polymer chains.

    This selector provides fine-grained access to polymer segments,
    useful for per-segment contact analysis.

    Parameters
    ----------
    residue_names : list[str], optional
        Residue names that identify polymer residues.
    chain_index : int, optional
        Specific chain to select segments from (0-indexed).
        If None, selects from all chains.
    segment_indices : list[int], optional
        Specific segment indices within chains to select.
        Uses 0-indexed positions within each chain.

    Notes
    -----
    A "segment" in this context refers to a single residue/monomer unit
    within a polymer chain, not MDAnalysis segments.
    """

    def __init__(
        self,
        residue_names: list[str] | None = None,
        chain_index: int | None = None,
        segment_indices: list[int] | None = None,
    ):
        self.residue_names = residue_names or DEFAULT_POLYMER_RESNAMES
        self.chain_index = chain_index
        self.segment_indices = segment_indices

    def select(self, universe: "Universe") -> SelectionResult:
        """Select polymer segments."""
        # First get all polymer atoms
        resname_str = " ".join(self.residue_names)
        all_polymer = universe.select_atoms(f"resname {resname_str}")

        if len(all_polymer) == 0:
            raise ValueError(f"No polymer atoms found with residue names: {self.residue_names}")

        # Get fragments (chains)
        fragments = all_polymer.fragments
        if not fragments:
            fragments = [all_polymer]

        # Select specific chain if requested
        if self.chain_index is not None:
            if self.chain_index >= len(fragments):
                raise ValueError(
                    f"Chain index {self.chain_index} out of range. Found {len(fragments)} chains."
                )
            fragments = [fragments[self.chain_index]]

        # Collect residues, optionally filtering by segment index
        selected_residues = []
        for frag in fragments:
            residues = frag.residues
            if self.segment_indices is not None:
                for idx in self.segment_indices:
                    if idx < len(residues):
                        selected_residues.append(residues[idx])
            else:
                selected_residues.extend(residues)

        if not selected_residues:
            raise ValueError("No polymer segments matched the selection criteria")

        # Combine into single AtomGroup
        atoms = selected_residues[0].atoms
        for res in selected_residues[1:]:
            atoms = atoms | res.atoms

        return SelectionResult(
            atoms=atoms,
            residues=atoms.residues,
            label=self.label,
            metadata={
                "residue_names": self.residue_names,
                "chain_index": self.chain_index,
                "segment_indices": self.segment_indices,
                "n_segments": len(selected_residues),
            },
        )

    @property
    def label(self) -> str:
        parts = ["polymer_segments"]
        if self.chain_index is not None:
            parts.append(f"chain{self.chain_index}")
        return "_".join(parts)
