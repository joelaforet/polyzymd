"""Modifier-to-protein linker abstractions for POC construction workflows."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Sequence
from pathlib import Path

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation.contracts import (
    ExplicitLinkageContract,
    LinkageBond,
    PdbAtomSelector,
    ReactiveEndpoint,
    ResolvedAttachmentPlan,
    resolve_explicit_linkage_contract,
)
from polyzymd.builders.conjugation.pdb_assembly import NhsLysPdbAttachment, PdbAtomRecord
from polyzymd.builders.conjugation.polymer_fragment import GeneratedPolymerFragment

_ATOM_RECORD_PREFIXES = ("ATOM", "HETATM")


class ProteinLinkSite(BaseModel):
    """Resolved protein atom site used for modifier placement."""

    atom: PdbAtomRecord
    removable_hydrogens: tuple[PdbAtomRecord, ...] = Field(default_factory=tuple)
    product_residue_name: str
    linker_atom_name: str


class ModifierLinkageSpec(BaseModel):
    """Resolved modifier and protein linkage names for Pablo and PDB assembly."""

    protein_residue_name: str
    modifier_residue_name: str
    protein_atom_name: str
    modifier_atom_name: str
    protein_leaving_atom_names: tuple[str, ...] = Field(default_factory=tuple)
    modifier_leaving_atom_names: tuple[str, ...] = Field(default_factory=tuple)
    bond_order: int | float = 1


class ModifierLinker(ABC):
    """Base class for modifier-to-protein covalent linker strategies."""

    @abstractmethod
    def resolve_site(self, protein_pdb_path: Path | str) -> ProteinLinkSite:
        """Resolve the protein atom site from a PDB file.

        Parameters
        ----------
        protein_pdb_path : pathlib.Path or str
            Protein PDB containing the target residue.

        Returns
        -------
        ProteinLinkSite
            Resolved protein atom and leaving hydrogens.
        """

    @abstractmethod
    def linkage_spec(self, modifier: GeneratedPolymerFragment) -> ModifierLinkageSpec:
        """Return product residue and atom names for a generated modifier.

        Parameters
        ----------
        modifier : GeneratedPolymerFragment
            Generated modifier fragment with a resolved reactive atom.

        Returns
        -------
        ModifierLinkageSpec
            Crosslink specification aligned with Pablo settings.
        """

    def resolve_plan(
        self, protein_pdb_path: Path | str, modifier: GeneratedPolymerFragment
    ) -> ResolvedAttachmentPlan:
        """Resolve this linker to a generic attachment plan.

        Parameters
        ----------
        protein_pdb_path : pathlib.Path or str
            Protein PDB containing the target residue.
        modifier : GeneratedPolymerFragment
            Generated modifier fragment with a resolved reactive atom.

        Returns
        -------
        ResolvedAttachmentPlan
            Generic atom-level linkage plan.
        """
        raise NotImplementedError("This linker does not provide a generic attachment plan")


class NhsLysModifierLinker(ModifierLinker):
    """NHS ester to lysine linker realization for modifier placement."""

    def __init__(
        self,
        *,
        target_chain: str = "A",
        target_residue_number: int,
        target_residue_name: str = "LYS",
        target_insertion_code: str = "",
        target_atom_name: str = "NZ",
        lysine_target_resname: str = "LYX",
        modifier_target_resname: str = "NHX",
        max_nz_hydrogens_to_remove: int = 2,
    ) -> None:
        """Initialize the NHS-Lys linker.

        Parameters
        ----------
        target_chain : str, optional
            Protein chain containing the target lysine, by default ``"A"``.
        target_residue_number : int
            PDB residue number of the target lysine.
        target_residue_name : str, optional
            Input protein residue name for the target site, by default ``"LYS"``.
        target_insertion_code : str, optional
            Input protein insertion code for the target site, by default ``""``.
        target_atom_name : str, optional
            Lysine atom joined to the modifier, by default ``"NZ"``.
        lysine_target_resname : str, optional
            Product residue name for linked lysine, by default ``"LYX"``.
        modifier_target_resname : str, optional
            Product residue name for linked NHS modifier residue, by default ``"NHX"``.
        max_nz_hydrogens_to_remove : int, optional
            Maximum number of lysine NZ hydrogens to remove, by default 2.
        """
        if len(target_chain) != 1:
            raise ValueError("target_chain must be a single PDB chain ID")
        self.target_chain = target_chain.upper()
        self.target_residue_number = target_residue_number
        self.target_residue_name = target_residue_name.upper()
        self.target_insertion_code = target_insertion_code.strip().upper()
        self.target_atom_name = target_atom_name.upper()
        self.lysine_target_resname = lysine_target_resname.upper()
        self.modifier_target_resname = modifier_target_resname.upper()
        self.max_nz_hydrogens_to_remove = max_nz_hydrogens_to_remove

    def resolve_site(self, protein_pdb_path: Path | str) -> ProteinLinkSite:
        """Resolve the target lysine NZ atom and removable hydrogens."""
        atoms = _parse_pdb_atoms(Path(protein_pdb_path))
        target_atoms = [
            atom
            for atom in atoms
            if atom.chain_id == self.target_chain
            and atom.residue_name.upper() == self.target_residue_name
            and atom.residue_number == self.target_residue_number
            and (atom.insertion_code or "").upper() == self.target_insertion_code
            and atom.atom_name.upper() == self.target_atom_name
        ]
        if len(target_atoms) != 1:
            raise ValueError(
                "Expected exactly one protein linker atom "
                f"{self.target_chain}:{self.target_residue_number}:{self.target_atom_name}, "
                f"found {len(target_atoms)}"
            )

        leaving_names = tuple(
            f"HZ{index}" for index in range(2, 2 + self.max_nz_hydrogens_to_remove)
        )
        removable = tuple(
            _resolve_target_residue_atom(atoms, self, atom_name) for atom_name in leaving_names
        )
        return ProteinLinkSite(
            atom=target_atoms[0],
            removable_hydrogens=removable,
            product_residue_name=self.lysine_target_resname,
            linker_atom_name=self.target_atom_name,
        )

    def linkage_spec(self, modifier: GeneratedPolymerFragment) -> ModifierLinkageSpec:
        """Return LYX-NHX crosslink names through the generic contract path."""
        contract = self.generic_contract(modifier)
        return ModifierLinkageSpec(
            protein_residue_name=contract.protein_endpoint.product_residue_name,
            modifier_residue_name=contract.modifier_endpoint.product_residue_name,
            protein_atom_name=contract.protein_endpoint.selector.atom_name,
            modifier_atom_name=contract.modifier_endpoint.selector.atom_name,
            protein_leaving_atom_names=contract.protein_endpoint.leaving_atom_names,
            modifier_leaving_atom_names=contract.modifier_endpoint.leaving_atom_names,
            bond_order=contract.bond.bond_order,
        )

    def generic_contract(self, modifier: GeneratedPolymerFragment) -> ExplicitLinkageContract:
        """Return the NHS-Lys convenience contract in generic form.

        Parameters
        ----------
        modifier : GeneratedPolymerFragment
            Generated modifier fragment with explicit reactive and leaving
            atom selectors.

        Returns
        -------
        ExplicitLinkageContract
            Generic linkage contract equivalent to the legacy NHS-Lys helper.
        """
        reactive_atom = resolve_modifier_reactive_atom(modifier)
        leaving_atoms = resolve_modifier_leaving_atoms(modifier)
        protein_selector = PdbAtomSelector(
            chain_id=self.target_chain,
            residue_name=self.target_residue_name,
            residue_number=self.target_residue_number,
            insertion_code=self.target_insertion_code,
            atom_name=self.target_atom_name,
        )
        modifier_selector = PdbAtomSelector(
            chain_id=reactive_atom.chain_id,
            residue_name=reactive_atom.residue_name,
            residue_number=reactive_atom.residue_number,
            insertion_code=reactive_atom.insertion_code,
            atom_name=reactive_atom.atom_name,
            atom_serial=reactive_atom.serial,
            atom_index=reactive_atom.atom_index,
        )
        return ExplicitLinkageContract(
            protein_endpoint=ReactiveEndpoint(
                participant="protein",
                selector=protein_selector,
                product_residue_name=self.lysine_target_resname,
                leaving_atom_names=tuple(
                    f"HZ{index}" for index in range(2, 2 + self.max_nz_hydrogens_to_remove)
                ),
            ),
            modifier_endpoint=ReactiveEndpoint(
                participant="modifier",
                selector=modifier_selector,
                product_residue_name=self.modifier_target_resname,
                leaving_atom_names=tuple(atom.atom_name.upper() for atom in leaving_atoms),
                leaving_atom_selectors=tuple(
                    PdbAtomSelector(
                        chain_id=atom.chain_id,
                        residue_name=atom.residue_name,
                        residue_number=atom.residue_number,
                        insertion_code=atom.insertion_code,
                        atom_name=atom.atom_name,
                        atom_serial=atom.serial,
                        atom_index=atom.atom_index,
                    )
                    for atom in leaving_atoms
                ),
            ),
            bond=LinkageBond(
                protein_atom_name=self.target_atom_name,
                modifier_atom_name=reactive_atom.atom_name,
                bond_order=1,
            ),
            mechanism_name="nhs_lys",
        )

    def resolve_plan(
        self, protein_pdb_path: Path | str, modifier: GeneratedPolymerFragment
    ) -> ResolvedAttachmentPlan:
        """Resolve the NHS-Lys helper through the generic PDB contract."""
        return resolve_explicit_linkage_contract(
            protein_pdb_path,
            modifier,
            self.generic_contract(modifier),
        )

    def attachment(self, protein_pdb_path: Path | str) -> NhsLysPdbAttachment:
        """Build the existing NHS-Lys PDB assembly attachment model.

        Parameters
        ----------
        protein_pdb_path : pathlib.Path or str
            Protein PDB containing the target lysine.

        Returns
        -------
        NhsLysPdbAttachment
            Attachment data compatible with :func:`write_crosslinked_pdb`.
        """
        site = self.resolve_site(protein_pdb_path)
        return NhsLysPdbAttachment(
            target_chain=self.target_chain,
            target_residue_name=self.target_residue_name,
            target_residue_number=self.target_residue_number,
            target_insertion_code=self.target_insertion_code,
            target_atom_name=self.target_atom_name,
            nz_hydrogen_atom_names_to_remove=tuple(
                atom.atom_name for atom in site.removable_hydrogens
            ),
            nz_hydrogen_atom_indices_to_remove=tuple(
                atom.atom_index for atom in site.removable_hydrogens if atom.atom_index is not None
            ),
            nz_hydrogen_atom_serials_to_remove=tuple(
                atom.serial for atom in site.removable_hydrogens if atom.serial is not None
            ),
            max_nz_hydrogens_to_remove=self.max_nz_hydrogens_to_remove,
            lysine_target_resname=self.lysine_target_resname,
            polymer_target_resname=self.modifier_target_resname,
        )


def resolve_modifier_reactive_atom(modifier: GeneratedPolymerFragment) -> PdbAtomRecord:
    """Resolve the single reactive atom in a generated modifier.

    Parameters
    ----------
    modifier : GeneratedPolymerFragment
        Generated modifier fragment with reactive selectors.

    Returns
    -------
    PdbAtomRecord
        Reactive atom converted to the PDB assembly atom model.
    """
    matches = [
        atom.to_pdb_atom()
        for atom in modifier.atoms
        if (
            modifier.reactive_atom_serial is not None
            and atom.serial == modifier.reactive_atom_serial
        )
        or (
            modifier.reactive_atom_index is not None
            and atom.atom_index == modifier.reactive_atom_index
        )
        or (
            modifier.reactive_atom_name is not None
            and atom.atom_name.upper() == modifier.reactive_atom_name.upper()
        )
    ]
    if len(matches) != 1:
        raise ValueError(f"Modifier reactive atom selector resolved {len(matches)} atoms")
    return matches[0]


def resolve_modifier_leaving_atoms(modifier: GeneratedPolymerFragment) -> tuple[PdbAtomRecord, ...]:
    """Resolve modifier leaving atoms from serial, index, or name selectors."""
    serials = set(modifier.leaving_atom_serials)
    indices = set(modifier.leaving_atom_indices)
    names = {name.upper() for name in modifier.leaving_atom_names}
    matches = [
        atom.to_pdb_atom()
        for atom in modifier.atoms
        if (atom.serial is not None and atom.serial in serials)
        or atom.atom_index in indices
        or atom.atom_name.upper() in names
    ]
    return tuple(
        sorted(matches, key=lambda atom: atom.atom_index if atom.atom_index is not None else -1)
    )


def _parse_pdb_atoms(path: Path) -> tuple[PdbAtomRecord, ...]:
    """Parse atom records from a PDB file."""
    atoms: list[PdbAtomRecord] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(_ATOM_RECORD_PREFIXES):
                atoms.append(PdbAtomRecord.from_pdb_line(line, atom_index=len(atoms)))
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records found in {path}")
    return tuple(atoms)


def _is_hydrogen_atom(atom: PdbAtomRecord) -> bool:
    """Return whether an atom record represents hydrogen."""
    return (atom.element or "").upper() == "H" or atom.atom_name.upper().startswith("H")


def _resolve_target_residue_atom(
    atoms: Sequence[PdbAtomRecord], linker: NhsLysModifierLinker, atom_name: str
) -> PdbAtomRecord:
    """Resolve one explicitly named atom on the NHS-Lys target residue."""
    matches = [
        atom
        for atom in atoms
        if atom.chain_id == linker.target_chain
        and atom.residue_name.upper() == linker.target_residue_name
        and atom.residue_number == linker.target_residue_number
        and (atom.insertion_code or "").upper() == linker.target_insertion_code
        and atom.atom_name.upper() == atom_name.upper()
    ]
    if len(matches) != 1:
        raise ValueError(
            "Expected exactly one protein leaving atom "
            f"{linker.target_chain}:{linker.target_residue_name}:"
            f"{linker.target_residue_number}{linker.target_insertion_code}:{atom_name}, "
            f"found {len(matches)}"
        )
    return matches[0]


def atom_names(atoms: Sequence[PdbAtomRecord]) -> tuple[str, ...]:
    """Return uppercase atom names for diagnostics and tests."""
    return tuple(atom.atom_name.upper() for atom in atoms)
