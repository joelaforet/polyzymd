"""Modifier-to-protein linker abstractions for POC construction workflows."""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from collections.abc import Sequence
from math import dist
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation.contracts import (
    ExplicitLinkageContract,
    LinkageBond,
    PdbAtomSelector,
    ReactiveEndpoint,
    ResolvedAttachmentPlan,
    resolve_explicit_linkage_contract,
)
from polyzymd.builders.conjugation.nhs_lys import detect_nhs_reactive_group
from polyzymd.builders.conjugation.pdb_assembly import NhsLysPdbAttachment, PdbAtomRecord
from polyzymd.builders.conjugation.polymer_fragment import GeneratedPolymerFragment
from polyzymd.builders.conjugation.rdkit_inputs import AtomIdentity, load_pdb_as_rdkit_input

LOGGER = logging.getLogger(__name__)

_ATOM_RECORD_PREFIXES = ("ATOM", "HETATM")


class ProteinLinkSite(BaseModel):
    """Resolved protein atom site used for modifier placement."""

    atom: PdbAtomRecord
    removable_hydrogens: tuple[PdbAtomRecord, ...] = Field(default_factory=tuple)
    product_residue_name: str
    linker_atom_name: str
    warnings: tuple[str, ...] = Field(default_factory=tuple)


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
                "Expected exactly one protein link atom "
                f"{self.target_chain}:{self.target_residue_number}:{self.target_atom_name}, "
                f"found {len(target_atoms)}"
            )

        removable, warnings = _resolve_lysine_leaving_hydrogens(
            protein_pdb_path,
            atoms,
            self,
            target_atoms[0],
        )
        return ProteinLinkSite(
            atom=target_atoms[0],
            removable_hydrogens=removable,
            product_residue_name=self.lysine_target_resname,
            linker_atom_name=self.target_atom_name,
            warnings=warnings,
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

    def generic_contract(
        self,
        modifier: GeneratedPolymerFragment,
        *,
        protein_pdb_path: Path | str | None = None,
    ) -> ExplicitLinkageContract:
        """Return the NHS-Lys convenience contract in generic form.

        Parameters
        ----------
        modifier : GeneratedPolymerFragment
            Generated modifier fragment with explicit reactive and leaving
            atom selectors.
        protein_pdb_path : pathlib.Path, str, or None, optional
            Optional source PDB used to resolve lysine leaving hydrogens by
            RDKit graph chemistry, by default ``None``. When omitted, no
            protein-side leaving atoms are emitted because atom names are only
            metadata hints.

        Returns
        -------
        ExplicitLinkageContract
            Generic linkage contract equivalent to the legacy NHS-Lys helper.
        """
        reactive_atom, leaving_atoms = resolve_modifier_nhs_atoms(modifier)
        protein_leaving_selectors: tuple[PdbAtomSelector, ...] = ()
        if protein_pdb_path is not None:
            site = self.resolve_site(protein_pdb_path)
            protein_leaving_selectors = tuple(
                _selector_from_pdb_atom(atom) for atom in site.removable_hydrogens
            )
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
                leaving_atom_selectors=protein_leaving_selectors,
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
            self.generic_contract(modifier, protein_pdb_path=protein_pdb_path),
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


def resolve_modifier_nhs_atoms(
    modifier: GeneratedPolymerFragment,
) -> tuple[PdbAtomRecord, tuple[PdbAtomRecord, ...]]:
    """Resolve modifier reactive and leaving atoms, preferring RDKit chemistry.

    Parameters
    ----------
    modifier : GeneratedPolymerFragment
        Generated fragment. Current generated POC fragments usually carry
        explicit selectors but not a complete RDKit topology.

    Returns
    -------
    tuple
        Reactive atom and leaving atoms as PDB records.
    """
    rdkit_failure: ValueError | None = None
    mol = getattr(modifier, "rdkit_mol", None) or getattr(modifier, "mol", None)
    if mol is not None:
        try:
            group = detect_nhs_reactive_group(mol)
        except ValueError as exc:
            rdkit_failure = exc
        else:
            atoms_by_index = {atom.atom_index: atom.to_pdb_atom() for atom in modifier.atoms}
            reactive_atom = atoms_by_index.get(group.reactive_carbon_index)
            leaving_atoms = tuple(
                atoms_by_index[index]
                for index in group.leaving_group_atom_indices
                if index in atoms_by_index
            )
            if reactive_atom is not None and len(leaving_atoms) == len(
                group.leaving_group_atom_indices
            ):
                return reactive_atom, leaving_atoms
            rdkit_failure = ValueError(
                "RDKit NHS detection returned atom indices that could not be mapped to the "
                "generated fragment: reactive carbon "
                f"{group.reactive_carbon_index}, leaving atoms {group.leaving_group_atom_indices}"
            )

    # Generated-fragment fallback: explicit selectors are provenance from the
    # polymer construction layer until complete modifier topology is retained
    try:
        fallback_reactive_atom = resolve_modifier_reactive_atom(modifier)
        fallback_leaving_atoms = resolve_modifier_leaving_atoms(modifier)
    except ValueError as fallback_exc:
        if rdkit_failure is not None:
            raise ValueError(
                "RDKit NHS detection failed and explicit generated-fragment fallback also "
                "failed. RDKit failure: "
                f"{rdkit_failure}. Fallback failure: {fallback_exc}"
            ) from fallback_exc
        raise
    if rdkit_failure is not None:
        LOGGER.warning(
            "RDKit NHS detection failed for modifier %s; using explicit generated-fragment "
            "fallback reactive atom %s and leaving atoms %s. RDKit failure: %s",
            modifier.name,
            fallback_reactive_atom.atom_name,
            tuple(atom.atom_name for atom in fallback_leaving_atoms),
            rdkit_failure,
        )
    return fallback_reactive_atom, fallback_leaving_atoms


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


def _resolve_lysine_leaving_hydrogens(
    protein_pdb_path: Path | str,
    atoms: Sequence[PdbAtomRecord],
    linker: NhsLysModifierLinker,
    target_atom: PdbAtomRecord,
) -> tuple[tuple[PdbAtomRecord, ...], tuple[str, ...]]:
    """Resolve NZ-bound hydrogens with RDKit graph chemistry and local fallback."""
    warnings: list[str] = []
    bundle = load_pdb_as_rdkit_input(protein_pdb_path, sanitize=False, proximity_bonding=True)
    target_identity = _identity_for_atom(bundle.atom_identities, target_atom)
    hydrogens = _rdkit_bound_hydrogens(bundle, target_identity, atoms)
    if not hydrogens:
        hydrogens = _geometry_bound_hydrogens(atoms, target_atom)
        if hydrogens:
            warnings.append(
                "RDKit did not infer NZ-hydrogen bonds; used same-residue distance fallback"
            )

    removable = _select_leaving_hydrogens(hydrogens, linker)
    if len(hydrogens) == 3:
        kept = _sort_hydrogen_records(hydrogens)[0]
        warnings.append(
            "Protonated lysine NZ had three bound hydrogens; kept "
            f"{kept.atom_name} and removed the other two by deterministic atom order"
        )
    return removable, tuple(warnings)


def _identity_for_atom(identities: Sequence[AtomIdentity], atom: PdbAtomRecord) -> AtomIdentity:
    """Return the RDKit identity matching a parsed PDB atom record."""
    matches = [
        identity
        for identity in identities
        if identity.atom_index == atom.atom_index
        and (atom.serial is None or identity.serial == atom.serial)
        and identity.atom_name.upper() == atom.atom_name.upper()
        and identity.same_residue_as(atom)
    ]
    if len(matches) != 1:
        raise ValueError(
            "Expected exactly one RDKit identity for lysine target atom "
            f"{atom.chain_id}:{atom.residue_name}:{atom.residue_number}:{atom.atom_name}, "
            f"found {len(matches)}"
        )
    return matches[0]


def _rdkit_bound_hydrogens(
    bundle: Any,
    target_identity: AtomIdentity,
    atoms: Sequence[PdbAtomRecord],
) -> tuple[PdbAtomRecord, ...]:
    """Return hydrogens bonded to NZ according to the RDKit molecule graph."""
    atom_by_index = {atom.atom_index: atom for atom in atoms if atom.atom_index is not None}
    target_atom = bundle.mol.GetAtomWithIdx(target_identity.rdkit_index)
    hydrogens: list[PdbAtomRecord] = []
    for neighbor in target_atom.GetNeighbors():
        if neighbor.GetAtomicNum() != 1:
            continue
        identity = bundle.identity_for_rdkit_index(neighbor.GetIdx())
        if not _same_identity_residue(identity, target_identity):
            continue
        pdb_atom = atom_by_index.get(identity.atom_index)
        if pdb_atom is not None:
            hydrogens.append(pdb_atom)
    return tuple(_sort_hydrogen_records(hydrogens))


def _geometry_bound_hydrogens(
    atoms: Sequence[PdbAtomRecord], target_atom: PdbAtomRecord, *, cutoff: float = 1.35
) -> tuple[PdbAtomRecord, ...]:
    """Return same-residue hydrogens within a local N-H distance cutoff."""
    target_xyz = (target_atom.x, target_atom.y, target_atom.z)
    hydrogens = [
        atom
        for atom in atoms
        if _same_pdb_residue(atom, target_atom)
        and _is_hydrogen_atom(atom)
        and dist(target_xyz, (atom.x, atom.y, atom.z)) <= cutoff
    ]
    return tuple(_sort_hydrogen_records(hydrogens))


def _select_leaving_hydrogens(
    hydrogens: Sequence[PdbAtomRecord], linker: NhsLysModifierLinker
) -> tuple[PdbAtomRecord, ...]:
    """Apply the initial NHS-Lys product proton policy to NZ hydrogens."""
    sorted_hydrogens = _sort_hydrogen_records(hydrogens)
    count = len(sorted_hydrogens)
    if count == 3:
        remove_count = 2
    elif count == 2:
        remove_count = 1
    else:
        raise ValueError(
            "NHS-Lys convenience resolution requires explicit lysine NZ hydrogens. "
            f"Expected 2 neutral or 3 protonated N-bound hydrogens, found {count}. "
            "Automatic hydrogen addition or protonation normalization is not implemented."
        )
    if remove_count > linker.max_nz_hydrogens_to_remove:
        raise ValueError(
            "NHS-Lys product policy requires removing "
            f"{remove_count} NZ hydrogens, but max_nz_hydrogens_to_remove is "
            f"{linker.max_nz_hydrogens_to_remove}"
        )
    return tuple(sorted_hydrogens[1 : 1 + remove_count])


def _sort_hydrogen_records(hydrogens: Sequence[PdbAtomRecord]) -> list[PdbAtomRecord]:
    """Return hydrogens in deterministic source order."""
    return sorted(
        hydrogens,
        key=lambda atom: (
            atom.atom_index if atom.atom_index is not None else 10**9,
            atom.serial if atom.serial is not None else 10**9,
            atom.atom_name.upper(),
        ),
    )


def _selector_from_pdb_atom(atom: PdbAtomRecord) -> PdbAtomSelector:
    """Build an exact selector for a resolved PDB atom."""
    return PdbAtomSelector(
        chain_id=atom.chain_id,
        residue_name=atom.residue_name,
        residue_number=atom.residue_number,
        insertion_code=atom.insertion_code,
        atom_name=atom.atom_name,
        atom_serial=atom.serial,
        atom_index=atom.atom_index,
    )


def _same_identity_residue(left: AtomIdentity, right: AtomIdentity) -> bool:
    """Return whether two atom identities share residue identity."""
    return (
        left.chain_id == right.chain_id
        and left.residue_name == right.residue_name
        and left.residue_number == right.residue_number
        and left.insertion_code == right.insertion_code
    )


def _same_pdb_residue(left: PdbAtomRecord, right: PdbAtomRecord) -> bool:
    """Return whether two PDB atoms belong to the same residue."""
    return (
        left.chain_id.upper() == right.chain_id.upper()
        and left.residue_name.upper() == right.residue_name.upper()
        and left.residue_number == right.residue_number
        and left.insertion_code.upper() == right.insertion_code.upper()
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
