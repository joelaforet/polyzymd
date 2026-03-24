"""Solvent, cosolvent, and substrate selectors.

This module provides selectors for non-protein, non-polymer molecules:

- SolventMolecules: Select water molecules
- CosolventMolecules: Select cosolvent (e.g., DMSO)
- SubstrateMolecule: Select substrate/ligand molecules

Examples
--------
>>> # Select water molecules
>>> selector = SolventMolecules()
>>>
>>> # Select DMSO cosolvent
>>> selector = CosolventMolecules(residue_names=["DMSO", "DMS"])
>>>
>>> # Select substrate by residue name
>>> selector = SubstrateMolecule(residue_name="RBU")  # Resorufin butyrate
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from polyzymd.analysis.common.selectors.base import MolecularSelector, SelectionResult

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe


# Common water residue names across force fields
DEFAULT_WATER_RESNAMES = [
    "HOH",  # PDB standard
    "WAT",  # Common
    "TIP3",  # TIP3P
    "TIP4",  # TIP4P
    "TIP5",  # TIP5P
    "SOL",  # GROMACS
    "SPC",  # SPC water
    "T3P",  # TIP3P variant
    "W",  # Coarse-grained
]

# Common cosolvent residue names
# NOTE: ACE is intentionally excluded - it's commonly used for N-terminal
# acetyl caps on proteins, not acetone cosolvent. Use "ACO" for acetone.
DEFAULT_COSOLVENT_RESNAMES = [
    "DMSO",  # Dimethyl sulfoxide
    "DMS",  # DMSO variant
    "ACN",  # Acetonitrile
    "ACO",  # Acetone (not ACE - that's often acetyl cap)
    "MeOH",  # Methanol
    "MEOH",
    "EtOH",  # Ethanol
    "ETOH",
    "THF",  # Tetrahydrofuran
    "GLYC",  # Glycerol
    "GOL",  # Glycerol variant
]


class SolventMolecules(MolecularSelector):
    """Select solvent (water) molecules.

    Parameters
    ----------
    residue_names : list[str], optional
        Residue names for water. Default uses common water names.
    exclude_near : str, optional
        Exclude waters within a cutoff of this selection.
        E.g., "protein" to exclude waters in first hydration shell.
    exclude_cutoff : float, optional
        Cutoff in Angstroms for exclude_near. Default 3.0.

    Examples
    --------
    >>> # Select all water
    >>> selector = SolventMolecules()
    >>>
    >>> # Select bulk water (exclude first shell around protein)
    >>> selector = SolventMolecules(
    ...     exclude_near="protein",
    ...     exclude_cutoff=5.0
    ... )
    """

    def __init__(
        self,
        residue_names: list[str] | None = None,
        exclude_near: str | None = None,
        exclude_cutoff: float = 3.0,
    ):
        self.residue_names = residue_names or DEFAULT_WATER_RESNAMES
        self.exclude_near = exclude_near
        self.exclude_cutoff = exclude_cutoff

    def select(self, universe: "Universe") -> SelectionResult:
        """Select water molecules."""
        resname_str = " ".join(self.residue_names)
        selection = f"resname {resname_str}"

        if self.exclude_near:
            # Exclude waters near the specified selection
            selection = (
                f"({selection}) and not (around {self.exclude_cutoff} ({self.exclude_near}))"
            )

        atoms = universe.select_atoms(selection)

        if len(atoms) == 0:
            raise ValueError(f"No solvent molecules found with residue names: {self.residue_names}")

        return SelectionResult(
            atoms=atoms,
            residues=atoms.residues,
            label=self.label,
            metadata={
                "residue_names": self.residue_names,
                "exclude_near": self.exclude_near,
                "exclude_cutoff": self.exclude_cutoff,
                "n_molecules": len(atoms.residues),
            },
        )

    @property
    def label(self) -> str:
        if self.exclude_near:
            return f"bulk_water_{self.exclude_cutoff:.1f}A"
        return "water"


class CosolventMolecules(MolecularSelector):
    """Select cosolvent molecules (e.g., DMSO, acetonitrile).

    Parameters
    ----------
    residue_names : list[str], optional
        Residue names for cosolvent. Default uses common names.
        You should typically specify this for your system.

    Examples
    --------
    >>> # Select DMSO molecules
    >>> selector = CosolventMolecules(residue_names=["DMSO", "DMS"])
    """

    def __init__(
        self,
        residue_names: list[str] | None = None,
    ):
        self.residue_names = residue_names or DEFAULT_COSOLVENT_RESNAMES

    def select(self, universe: "Universe") -> SelectionResult:
        """Select cosolvent molecules."""
        resname_str = " ".join(self.residue_names)
        selection = f"resname {resname_str}"

        atoms = universe.select_atoms(selection)

        if len(atoms) == 0:
            raise ValueError(
                f"No cosolvent molecules found with residue names: {self.residue_names}. "
                "Specify the correct residue names for your cosolvent."
            )

        return SelectionResult(
            atoms=atoms,
            residues=atoms.residues,
            label=self.label,
            metadata={
                "residue_names": self.residue_names,
                "n_molecules": len(atoms.residues),
            },
        )

    @property
    def label(self) -> str:
        return "cosolvent"


class SubstrateMolecule(MolecularSelector):
    """Select substrate or ligand molecules.

    Parameters
    ----------
    residue_name : str
        Residue name of the substrate.
    n_molecules : int, optional
        Expected number of substrate molecules. If provided, validates
        that exactly this many are found. Default None (no validation).

    Examples
    --------
    >>> # Select resorufin butyrate substrate
    >>> selector = SubstrateMolecule(residue_name="RBU")
    >>>
    >>> # Select single substrate, validate count
    >>> selector = SubstrateMolecule(residue_name="RBU", n_molecules=1)
    """

    def __init__(
        self,
        residue_name: str,
        n_molecules: int | None = None,
    ):
        self.residue_name = residue_name
        self.n_molecules = n_molecules

    def select(self, universe: "Universe") -> SelectionResult:
        """Select substrate molecules."""
        selection = f"resname {self.residue_name}"
        atoms = universe.select_atoms(selection)

        if len(atoms) == 0:
            raise ValueError(f"No substrate molecules found with residue name: {self.residue_name}")

        n_found = len(atoms.residues)
        if self.n_molecules is not None and n_found != self.n_molecules:
            raise ValueError(
                f"Expected {self.n_molecules} substrate molecule(s), but found {n_found}"
            )

        return SelectionResult(
            atoms=atoms,
            residues=atoms.residues,
            label=self.label,
            metadata={
                "residue_name": self.residue_name,
                "n_molecules": n_found,
            },
        )

    @property
    def label(self) -> str:
        return f"substrate_{self.residue_name}"


class IonSelector(MolecularSelector):
    """Select ion molecules (Na+, Cl-, etc.).

    Parameters
    ----------
    residue_names : list[str], optional
        Residue names for ions. Default includes common ions.
    ion_type : {"all", "cations", "anions"}, optional
        Filter to specific ion types. Default "all".

    Examples
    --------
    >>> # Select all ions
    >>> selector = IonSelector()
    >>>
    >>> # Select only sodium ions
    >>> selector = IonSelector(residue_names=["NA", "Na+", "SOD"])
    """

    # Common ion residue names
    DEFAULT_CATIONS = ["NA", "Na+", "SOD", "K", "K+", "POT", "MG", "Mg2+", "CA", "Ca2+"]
    DEFAULT_ANIONS = ["CL", "Cl-", "CLA", "BR", "Br-"]

    def __init__(
        self,
        residue_names: list[str] | None = None,
        ion_type: str = "all",
    ):
        if ion_type not in ("all", "cations", "anions"):
            raise ValueError(f"ion_type must be 'all', 'cations', or 'anions', got {ion_type}")

        self.ion_type = ion_type

        if residue_names is not None:
            self.residue_names = residue_names
        elif ion_type == "cations":
            self.residue_names = self.DEFAULT_CATIONS
        elif ion_type == "anions":
            self.residue_names = self.DEFAULT_ANIONS
        else:
            self.residue_names = self.DEFAULT_CATIONS + self.DEFAULT_ANIONS

    def select(self, universe: "Universe") -> SelectionResult:
        """Select ion molecules."""
        resname_str = " ".join(self.residue_names)
        selection = f"resname {resname_str}"

        atoms = universe.select_atoms(selection)

        if len(atoms) == 0:
            raise ValueError(f"No ions found with residue names: {self.residue_names}")

        return SelectionResult(
            atoms=atoms,
            residues=atoms.residues,
            label=self.label,
            metadata={
                "residue_names": self.residue_names,
                "ion_type": self.ion_type,
                "n_ions": len(atoms.residues),
            },
        )

    @property
    def label(self) -> str:
        return f"ions_{self.ion_type}"
