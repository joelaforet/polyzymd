"""Internal linkage contracts, crosslink validation, and linker strategies."""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from collections.abc import Sequence
from math import dist
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field, computed_field, field_validator, model_validator

from polyzymd.builders.conjugation.polymer.fragment import (
    GeneratedPolymerFragment,
    PreparedFragment,
)
from polyzymd.builders.conjugation.structure.parsing import parse_pdb_atom_records as _parse_atoms
from polyzymd.builders.conjugation.structure.pdb import (
    NhsLysPdbAttachment,
    PdbAtomRecord,
    PdbLinkageAttachment,
    PlacedPolymerFragment,
)


class PdbAtomSelector(BaseModel):
    """Fully qualified PDB atom selector for linkage contracts."""

    chain_id: str = Field(..., min_length=1, max_length=1)
    residue_name: str = Field(..., min_length=1)
    residue_number: int
    atom_name: str = Field(..., min_length=1)
    insertion_code: str = Field("", max_length=1)
    atom_serial: int | None = Field(None, ge=1)
    atom_index: int | None = Field(None, ge=0)

    @field_validator("chain_id", "residue_name", "atom_name", "insertion_code")
    @classmethod
    def normalize_text_fields(cls, value: str) -> str:
        """Normalize PDB selector text to uppercase values."""
        return value.strip().upper()

    @model_validator(mode="after")
    def validate_required_text_fields(self) -> PdbAtomSelector:
        """Validate that normalized selector fields remain non-empty."""
        if not self.chain_id or not self.residue_name or not self.atom_name:
            raise ValueError("PDB atom selectors require non-empty chain, residue, and atom names")
        return self

    def matches(self, atom: PdbAtomRecord) -> bool:
        """Return whether an atom matches all supplied selector fields.

        Parameters
        ----------
        atom : PdbAtomRecord
            Parsed PDB atom record to test.

        Returns
        -------
        bool
            Whether the atom is selected by this fully qualified selector.
        """
        return (
            atom.chain_id.upper() == self.chain_id
            and atom.residue_name.upper() == self.residue_name
            and atom.residue_number == self.residue_number
            and (atom.insertion_code or "").upper() == self.insertion_code
            and atom.atom_name.upper() == self.atom_name
            and (self.atom_serial is None or atom.serial == self.atom_serial)
            and (self.atom_index is None or atom.atom_index == self.atom_index)
        )

    def same_residue(self, atom: PdbAtomRecord) -> bool:
        """Return whether an atom belongs to this selector's residue.

        Parameters
        ----------
        atom : PdbAtomRecord
            Parsed PDB atom record to test.

        Returns
        -------
        bool
            Whether the atom has the same chain, residue name, number, and
            insertion code as this selector.
        """
        return (
            atom.chain_id.upper() == self.chain_id
            and atom.residue_name.upper() == self.residue_name
            and atom.residue_number == self.residue_number
            and (atom.insertion_code or "").upper() == self.insertion_code
        )


class ReactiveEndpoint(BaseModel):
    """One participant endpoint in an explicit PDB linkage contract."""

    participant: Literal["protein", "modifier"]
    selector: PdbAtomSelector
    product_residue_name: str = Field(..., min_length=1, max_length=4)
    leaving_atom_selectors: tuple[PdbAtomSelector, ...] = Field(default_factory=tuple)
    leaving_atom_names: tuple[str, ...] = Field(default_factory=tuple)
    allow_external_leaving_residue: bool = False

    @field_validator("product_residue_name")
    @classmethod
    def normalize_product_residue_name(cls, value: str) -> str:
        """Normalize product residue names for downstream Pablo matching."""
        return value.strip().upper()

    @field_validator("leaving_atom_names", mode="before")
    @classmethod
    def normalize_leaving_atom_names(cls, value: object) -> tuple[str, ...]:
        """Normalize leaving atom name selectors to uppercase tuples."""
        if value is None:
            return ()
        return tuple(str(item).strip().upper() for item in value)

    @model_validator(mode="after")
    def validate_selector_scope(self) -> ReactiveEndpoint:
        """Validate that explicit leaving selectors remain in endpoint scope."""
        for selector in self.leaving_atom_selectors:
            if not _same_residue_selector(selector, self.selector):
                if self.allow_external_leaving_residue and selector.residue_name == "ROH":
                    continue
                raise ValueError(
                    "Leaving atom selectors must use the same chain, residue name, "
                    "residue number, and insertion code as the endpoint selector"
                )
        return self


class LinkageBond(BaseModel):
    """Bond instructions connecting the protein and modifier endpoints."""

    protein_atom_selector: PdbAtomSelector | None = None
    modifier_atom_selector: PdbAtomSelector | None = None
    protein_atom_name: str | None = None
    modifier_atom_name: str | None = None
    bond_order: int | float = 1
    target_bond_length_angstrom: float = Field(1.33, gt=0)

    @field_validator("protein_atom_name", "modifier_atom_name")
    @classmethod
    def normalize_optional_atom_name(cls, value: str | None) -> str | None:
        """Normalize optional atom-name aliases."""
        return value.strip().upper() if value is not None else None


class PabloCrosslinkRequirement(BaseModel):
    """Explicit Pablo crosslink requirement derived from a linkage contract."""

    residues: tuple[str, str]
    linking_atoms: tuple[str, str]
    leaving_atoms: tuple[tuple[str, ...], tuple[str, ...]] = Field(default_factory=tuple)
    bond_order: int | float = 1

    @field_validator("residues", "linking_atoms", mode="before")
    @classmethod
    def normalize_pair(cls, value: object) -> tuple[str, str]:
        """Normalize a two-name Pablo field."""
        normalized = tuple(str(item).strip().upper() for item in value)  # type: ignore[arg-type]
        if len(normalized) != 2:
            raise ValueError("Pablo crosslink fields must contain exactly two values")
        return normalized

    @field_validator("leaving_atoms", mode="before")
    @classmethod
    def normalize_leaving_atom_groups(
        cls, value: object
    ) -> tuple[tuple[str, ...], tuple[str, ...]]:
        """Normalize two leaving-atom groups for Pablo matching."""
        if value is None:
            return ((), ())
        groups = tuple(
            tuple(str(atom).strip().upper() for atom in group)
            for group in value  # type: ignore[union-attr]
        )
        if len(groups) != 2:
            raise ValueError("Pablo leaving_atoms must contain exactly two groups")
        return groups


class ExplicitLinkageContract(BaseModel):
    """Complete explicit PDB-only linkage contract."""

    protein_endpoint: ReactiveEndpoint
    modifier_endpoint: ReactiveEndpoint
    bond: LinkageBond = Field(default_factory=LinkageBond)
    mechanism_name: str | None = None

    @model_validator(mode="after")
    def validate_participants(self) -> ExplicitLinkageContract:
        """Validate endpoint participant roles."""
        if self.protein_endpoint.participant != "protein":
            raise ValueError("protein_endpoint participant must be 'protein'")
        if self.modifier_endpoint.participant != "modifier":
            raise ValueError("modifier_endpoint participant must be 'modifier'")
        return self


class ReactionProduct(BaseModel):
    """Canonical resolved product state shared by construction stages."""

    model_config = {"arbitrary_types_allowed": True}

    contract: ExplicitLinkageContract
    protein_link_atom: PdbAtomRecord
    modifier_link_atom: PdbAtomRecord
    protein_leaving_atoms: tuple[PdbAtomRecord, ...] = Field(default_factory=tuple)
    modifier_leaving_atoms: tuple[PdbAtomRecord, ...] = Field(default_factory=tuple)
    fragment: PreparedFragment
    attachment_id: str = ""
    attachment_index: int = Field(1, ge=1)
    attachment_config: Any = Field(default=None, exclude=True)
    reaction_name: str = ""
    product_residue_mappings: dict[str, dict[str, int | str]] = Field(
        default_factory=dict,
        exclude=True,
    )
    attachment_force_field_owner: str = "modifier"
    attachment_force_field_domain: str = ""
    endpoint_provenance: dict[str, Any] = Field(default_factory=dict)
    scoped_residue_aliases: dict[str, str] = Field(default_factory=dict)
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)

    @computed_field
    @property
    def protein_product_residue_name(self) -> str:
        """Return the contract-owned protein product residue name."""
        return self.contract.protein_endpoint.product_residue_name

    @computed_field
    @property
    def modifier_product_residue_name(self) -> str:
        """Return the contract-owned modifier product residue name."""
        return self.contract.modifier_endpoint.product_residue_name

    @computed_field
    @property
    def target_bond_length_angstrom(self) -> float:
        """Return the contract-owned target linkage distance."""
        return self.contract.bond.target_bond_length_angstrom

    @computed_field
    @property
    def pablo_crosslink_requirement(self) -> PabloCrosslinkRequirement:
        """Derive the Pablo requirement from the resolved atoms and contract."""
        return PabloCrosslinkRequirement(
            residues=(
                self.protein_product_residue_name,
                self.modifier_product_residue_name,
            ),
            linking_atoms=(self.protein_link_atom.atom_name, self.modifier_link_atom.atom_name),
            leaving_atoms=(
                tuple(atom.atom_name for atom in self.protein_leaving_atoms),
                tuple(atom.atom_name for atom in self.modifier_leaving_atoms),
            ),
            bond_order=self.contract.bond.bond_order,
        )

    def to_nhs_lys_pdb_attachment(self) -> NhsLysPdbAttachment:
        """Convert this resolved plan to the legacy NHS-Lys assembly adapter.

        Returns
        -------
        NhsLysPdbAttachment
            Legacy assembly attachment model populated from generic endpoint
            selectors for POC-compatible NHS-Lys construction paths.
        """
        selector = self.contract.protein_endpoint.selector
        return NhsLysPdbAttachment(
            target_chain=selector.chain_id,
            target_residue_name=selector.residue_name,
            target_residue_number=selector.residue_number,
            target_insertion_code=selector.insertion_code,
            target_atom_name=self.protein_link_atom.atom_name.upper(),
            nz_hydrogen_atom_serials_to_remove=tuple(
                atom.serial for atom in self.protein_leaving_atoms if atom.serial is not None
            ),
            nz_hydrogen_atom_indices_to_remove=tuple(
                atom.atom_index
                for atom in self.protein_leaving_atoms
                if atom.atom_index is not None
            ),
            nz_hydrogen_atom_names_to_remove=tuple(
                atom.atom_name for atom in self.protein_leaving_atoms
            ),
            max_nz_hydrogens_to_remove=len(self.protein_leaving_atoms),
            lysine_target_resname=self.protein_product_residue_name,
            polymer_target_resname=self.modifier_product_residue_name,
        )

    def to_pdb_linkage_attachment(self) -> PdbLinkageAttachment:
        """Convert this resolved plan to a generic PDB assembly attachment.

        Returns
        -------
        PdbLinkageAttachment
            Generic assembly attachment retaining exact resolved protein-side
            leaving atom identities and configured product residue names.
        """
        selector = self.contract.protein_endpoint.selector
        return PdbLinkageAttachment(
            target_chain=selector.chain_id,
            target_residue_name=selector.residue_name,
            target_residue_number=selector.residue_number,
            target_insertion_code=selector.insertion_code,
            target_atom_name=self.protein_link_atom.atom_name.upper(),
            protein_leaving_atoms_to_remove=self.protein_leaving_atoms,
            modifier_link_atom=self.modifier_link_atom,
            modifier_leaving_atoms_to_remove=self.modifier_leaving_atoms,
            protein_target_resname=self.protein_product_residue_name,
            modifier_target_resname=self.modifier_product_residue_name,
        )


def explicit_linkage_contract_from_config(attachment: Any) -> ExplicitLinkageContract:
    """Build an explicit PDB linkage contract from an attachment config.

    Parameters
    ----------
    attachment : Any
        Duck-typed ``ConjugationAttachmentConfig`` with validated explicit
        linkage fields.

    Returns
    -------
    ExplicitLinkageContract
        Generic explicit linkage contract used by construction and placement.
    """
    mechanism = attachment.mechanism
    site = attachment.site
    moiety_site = attachment.moiety.link_site
    if moiety_site is None:
        raise ValueError("explicit_linkage requires moiety.link_site")
    product_residues = mechanism.product_residues
    if product_residues.site is None or product_residues.moiety is None:
        raise ValueError("explicit_linkage requires both product residue names")

    protein_selector = PdbAtomSelector(
        chain_id=site.chain_id,
        residue_name=_required_config_value(site.residue_name, "site.residue_name"),
        residue_number=_required_config_value(site.residue_number, "site.residue_number"),
        atom_name=_required_config_value(site.atom_name, "site.atom_name"),
        insertion_code=getattr(site, "insertion_code", ""),
        atom_serial=getattr(site, "atom_serial", None),
        atom_index=getattr(site, "atom_index", None),
    )
    modifier_selector = PdbAtomSelector(
        chain_id=moiety_site.chain_id,
        residue_name=moiety_site.residue_name,
        residue_number=moiety_site.residue_number,
        atom_name=moiety_site.atom_name,
        insertion_code=moiety_site.insertion_code,
        atom_serial=moiety_site.atom_serial,
        atom_index=moiety_site.atom_index,
    )

    bond = mechanism.bond
    site_atom = bond.site_atom or mechanism.site_atom or protein_selector.atom_name
    moiety_atom = bond.moiety_atom or mechanism.moiety_atom or modifier_selector.atom_name
    return ExplicitLinkageContract(
        protein_endpoint=ReactiveEndpoint(
            participant="protein",
            selector=protein_selector,
            product_residue_name=product_residues.site,
            leaving_atom_names=tuple(mechanism.leaving_atoms.site),
        ),
        modifier_endpoint=ReactiveEndpoint(
            participant="modifier",
            selector=modifier_selector,
            product_residue_name=product_residues.moiety,
            leaving_atom_names=tuple(mechanism.leaving_atoms.moiety),
        ),
        bond=LinkageBond(
            protein_atom_name=site_atom,
            modifier_atom_name=moiety_atom,
            bond_order=bond.order,
            target_bond_length_angstrom=bond.target_bond_length_angstrom,
        ),
        mechanism_name=mechanism.name,
    )


def resolve_explicit_linkage_contract(
    protein_pdb_path: Path | str,
    modifier: (
        Path | str | GeneratedPolymerFragment | PlacedPolymerFragment | Sequence[PdbAtomRecord]
    ),
    contract: ExplicitLinkageContract,
    *,
    fragment: PreparedFragment,
    attachment_id: str = "",
    attachment_index: int = 1,
    attachment_config: Any = None,
    reaction_name: str = "",
    attachment_force_field_domain: str = "",
    diagnostics: tuple[str, ...] = (),
) -> ReactionProduct:
    """Resolve an explicit PDB linkage contract against protein and modifier atoms.

    Parameters
    ----------
    protein_pdb_path : pathlib.Path or str
        Protein PDB containing the selected protein endpoint atom.
    modifier : path-like, GeneratedPolymerFragment, PlacedPolymerFragment, or sequence
        Modifier PDB source or already parsed/generated modifier atom records.
    contract : ExplicitLinkageContract
        User-supplied explicit linkage contract.

    Returns
    -------
    ReactionProduct
        Atom-level plan with selected atoms, scoped leaving atoms, and the Pablo
        crosslink requirement.
    """
    protein_atoms = parse_pdb_atom_records(protein_pdb_path)
    modifier_atoms = _modifier_atom_records(modifier)

    protein_selector = contract.bond.protein_atom_selector or contract.protein_endpoint.selector
    modifier_selector = contract.bond.modifier_atom_selector or contract.modifier_endpoint.selector
    protein_link_atom = _resolve_single_atom(protein_atoms, protein_selector, "protein link atom")
    modifier_link_atom = _resolve_single_atom(
        modifier_atoms, modifier_selector, "modifier link atom"
    )

    _validate_bond_atom_name(contract.bond.protein_atom_name, protein_link_atom, "protein")
    _validate_bond_atom_name(contract.bond.modifier_atom_name, modifier_link_atom, "modifier")

    protein_leaving_atoms = _resolve_leaving_atoms(
        protein_atoms,
        contract.protein_endpoint,
        label="protein leaving atom",
    )
    modifier_leaving_atoms = _resolve_leaving_atoms(
        modifier_atoms,
        contract.modifier_endpoint,
        label="modifier leaving atom",
    )

    return ReactionProduct(
        contract=contract,
        protein_link_atom=protein_link_atom,
        modifier_link_atom=modifier_link_atom,
        protein_leaving_atoms=protein_leaving_atoms,
        modifier_leaving_atoms=modifier_leaving_atoms,
        fragment=fragment,
        attachment_id=attachment_id,
        attachment_index=attachment_index,
        attachment_config=attachment_config,
        reaction_name=reaction_name or contract.mechanism_name or "",
        attachment_force_field_domain=attachment_force_field_domain,
        diagnostics=diagnostics,
    )


def parse_pdb_atom_records(path: Path | str) -> tuple[PdbAtomRecord, ...]:
    """Parse all PDB ATOM/HETATM records from a file.

    Parameters
    ----------
    path : pathlib.Path or str
        PDB file path.

    Returns
    -------
    tuple of PdbAtomRecord
        Parsed atom records with zero-based ``atom_index`` values.
    """
    return _parse_atoms(path, require_atoms=True)


def _modifier_atom_records(
    modifier: (
        Path | str | GeneratedPolymerFragment | PlacedPolymerFragment | Sequence[PdbAtomRecord]
    ),
) -> tuple[PdbAtomRecord, ...]:
    """Normalize supported modifier atom sources to PDB atom records."""
    if isinstance(modifier, (str, Path)):
        return parse_pdb_atom_records(modifier)
    if isinstance(modifier, GeneratedPolymerFragment):
        return tuple(atom.to_pdb_atom() for atom in modifier.atoms)
    if isinstance(modifier, PlacedPolymerFragment):
        return tuple(modifier.atoms)
    return tuple(modifier)


def _resolve_single_atom(
    atoms: Sequence[PdbAtomRecord], selector: PdbAtomSelector, label: str
) -> PdbAtomRecord:
    """Resolve exactly one atom from a fully qualified selector."""
    matches = [atom for atom in atoms if selector.matches(atom)]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one {label} for selector "
            f"{selector.chain_id}:{selector.residue_name}:{selector.residue_number}"
            f"{selector.insertion_code}:{selector.atom_name}, found {len(matches)}"
        )
    return matches[0]


def _resolve_leaving_atoms(
    atoms: Sequence[PdbAtomRecord], endpoint: ReactiveEndpoint, *, label: str
) -> tuple[PdbAtomRecord, ...]:
    """Resolve endpoint leaving atoms scoped to the selected residue."""
    resolved: list[PdbAtomRecord] = []
    for selector in endpoint.leaving_atom_selectors:
        resolved.append(_resolve_single_atom(atoms, selector, label))

    scoped_atoms = [atom for atom in atoms if endpoint.selector.same_residue(atom)]
    for atom_name in endpoint.leaving_atom_names:
        matches = [atom for atom in scoped_atoms if atom.atom_name.upper() == atom_name]
        if len(matches) != 1:
            raise ValueError(
                f"Expected exactly one {label} named {atom_name} in selected residue "
                f"{endpoint.selector.chain_id}:{endpoint.selector.residue_name}:"
                f"{endpoint.selector.residue_number}, found {len(matches)}"
            )
        resolved.append(matches[0])

    unique: list[PdbAtomRecord] = []
    seen: set[tuple[int | None, int | None, str, int, str, str]] = set()
    for atom in resolved:
        key = (
            atom.atom_index,
            atom.serial,
            atom.atom_name.upper(),
            atom.residue_number,
            atom.insertion_code,
            atom.chain_id.upper(),
        )
        if key not in seen:
            seen.add(key)
            unique.append(atom)
    return tuple(unique)


def _validate_bond_atom_name(
    expected_name: str | None, atom: PdbAtomRecord, participant: str
) -> None:
    """Validate an optional bond atom-name alias against the resolved atom."""
    if expected_name is None:
        return
    if atom.atom_name.upper() != expected_name:
        raise ValueError(
            f"{participant.capitalize()} bond atom name {expected_name} does not match "
            f"resolved selector atom {atom.atom_name.upper()}"
        )


def _required_config_value(value: Any, field_name: str) -> Any:
    """Return a required config value or raise a clear adapter error."""
    if value is None:
        raise ValueError(f"explicit_linkage requires {field_name}")
    return value


def _same_residue_selector(left: PdbAtomSelector, right: PdbAtomSelector) -> bool:
    """Return whether two selectors address the same PDB residue."""
    return (
        left.chain_id == right.chain_id
        and left.residue_name == right.residue_name
        and left.residue_number == right.residue_number
        and left.insertion_code == right.insertion_code
    )


class CrosslinkValidationResult(BaseModel):
    """Resolved explicit Pablo crosslink validation details."""

    matched_index: int
    residues: tuple[str, str]
    linking_atoms: tuple[str, str]
    leaving_atoms: tuple[tuple[str, ...], tuple[str, ...]]
    bond_order: int | float = 1


class MissingPabloCrosslinkError(ValueError):
    """Raised when required Pablo crosslink configuration is absent."""

    def __init__(self, requirement: ModifierLinkageSpec | PabloCrosslinkRequirement) -> None:
        """Build an actionable missing-crosslink error."""
        super().__init__(_missing_crosslink_message(_as_requirement(requirement)))


def require_explicit_pablo_crosslink(
    policy: Any,
    linker: ModifierLinker,
    modifier: GeneratedPolymerFragment,
) -> CrosslinkValidationResult:
    """Require an explicit Pablo crosslink matching the linker realization.

    Parameters
    ----------
    policy : Any
        Duck-typed ``ccd_pablo`` policy containing a ``crosslinks`` field.
    linker : ModifierLinker
        Linker strategy that defines product residue and atom names.
    modifier : GeneratedPolymerFragment
        Generated modifier used to resolve the reactive atom name.

    Returns
    -------
    CrosslinkValidationResult
        Matched crosslink details.

    Raises
    ------
    MissingPabloCrosslinkError
        If no configured crosslink matches the product residue and atom names.
    """
    requirement = _as_requirement(linker.linkage_spec(modifier))
    return require_pablo_crosslink_requirement(policy, requirement)


def require_pablo_crosslink_requirement(
    policy: Any,
    requirement: PabloCrosslinkRequirement,
) -> CrosslinkValidationResult:
    """Require an explicit Pablo crosslink matching a resolved requirement.

    Parameters
    ----------
    policy : Any
        Duck-typed ``ccd_pablo`` policy containing a ``crosslinks`` field.
    requirement : PabloCrosslinkRequirement
        Generic resolved crosslink requirement.

    Returns
    -------
    CrosslinkValidationResult
        Matched crosslink details.

    Raises
    ------
    MissingPabloCrosslinkError
        If no configured crosslink matches the required product residue and
        atom names.
    """
    crosslinks = list(getattr(policy, "crosslinks", ()) if policy is not None else ())
    for index, crosslink in enumerate(crosslinks):
        normalized = _normalize_crosslink(crosslink)
        if _matches_requirement(normalized, requirement):
            return CrosslinkValidationResult(matched_index=index, **normalized)
    raise MissingPabloCrosslinkError(requirement)


def _normalize_crosslink(crosslink: Any) -> dict[str, Any]:
    """Normalize a duck-typed crosslink config."""
    return {
        "residues": _two_upper(_crosslink_attr(crosslink, "residues", ())),
        "linking_atoms": _two_upper(_crosslink_attr(crosslink, "linking_atoms", ())),
        "leaving_atoms": tuple(
            tuple(str(atom).strip().upper() for atom in group)
            for group in _crosslink_attr(crosslink, "leaving_atoms", ((), ()))
        ),
        "bond_order": _crosslink_attr(crosslink, "bond_order", 1),
    }


def _crosslink_attr(crosslink: Any, name: str, default: Any) -> Any:
    """Return a field from a mapping or object crosslink config."""
    if isinstance(crosslink, dict):
        return crosslink.get(name, default)
    return getattr(crosslink, name, default)


def _matches_requirement(crosslink: dict[str, Any], requirement: PabloCrosslinkRequirement) -> bool:
    """Return whether a normalized crosslink matches either residue order."""
    expected_forward = (
        requirement.residues,
        requirement.linking_atoms,
        requirement.leaving_atoms,
    )
    expected_reverse = (
        tuple(reversed(expected_forward[0])),
        tuple(reversed(expected_forward[1])),
        tuple(reversed(expected_forward[2])),
    )
    observed = (crosslink["residues"], crosslink["linking_atoms"], crosslink["leaving_atoms"])
    return (
        observed in {expected_forward, expected_reverse}
        and crosslink["bond_order"] == requirement.bond_order
    )


def _as_requirement(
    requirement: ModifierLinkageSpec | PabloCrosslinkRequirement,
) -> PabloCrosslinkRequirement:
    """Normalize legacy linkage specs to generic Pablo requirements."""
    if isinstance(requirement, PabloCrosslinkRequirement):
        return requirement
    return PabloCrosslinkRequirement(
        residues=(requirement.protein_residue_name, requirement.modifier_residue_name),
        linking_atoms=(requirement.protein_atom_name, requirement.modifier_atom_name),
        leaving_atoms=(
            requirement.protein_leaving_atom_names,
            requirement.modifier_leaving_atom_names,
        ),
        bond_order=requirement.bond_order,
    )


def _two_upper(values: Any) -> tuple[str, str]:
    """Normalize a two-name config field."""
    normalized = tuple(str(item).strip().upper() for item in values)
    if len(normalized) != 2:
        return ("", "")
    return normalized


def _missing_crosslink_message(requirement: PabloCrosslinkRequirement) -> str:
    """Return an actionable missing-crosslink diagnostic."""
    modifier_leaving = ", ".join(f'"{name}"' for name in requirement.leaving_atoms[1])
    if not modifier_leaving:
        modifier_leaving = '"<modifier leaving atom>"'
    return (
        "Required explicit ccd_pablo.crosslinks entry was not found for modifier linking to "
        "protein. Add a Pablo crosslink before running construction, for example:\n\n"
        "conjugation:\n"
        "  ccd_pablo:\n"
        "    crosslinks:\n"
        "      - residues: "
        f'["{requirement.residues[0]}", "{requirement.residues[1]}"]\n'
        "        linking_atoms: "
        f'["{requirement.linking_atoms[0]}", "{requirement.linking_atoms[1]}"]\n'
        "        leaving_atoms:\n"
        f"          - {list(requirement.leaving_atoms[0])}\n          - [{modifier_leaving}]\n"
        f"        bond_order: {requirement.bond_order}\n"
    )


class ExplicitCrosslinkRequirement(BaseModel):
    """Small serializable record documenting the required crosslink."""

    residues: tuple[str, str]
    linking_atoms: tuple[str, str]
    leaving_atoms: tuple[tuple[str, ...], tuple[str, ...]] = Field(default_factory=tuple)
    bond_order: int | float = 1


LOGGER = logging.getLogger(__name__)


def detect_nhs_reactive_group(mol: Any) -> Any:
    """Lazily resolve NHS reactive-group detection for tests and callers."""
    from polyzymd.builders.conjugation.reactions.nhs_lys import (
        detect_nhs_reactive_group as detector,
    )

    return detector(mol)


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
    ) -> ReactionProduct:
        """Resolve this linker to a generic attachment plan.

        Parameters
        ----------
        protein_pdb_path : pathlib.Path or str
            Protein PDB containing the target residue.
        modifier : GeneratedPolymerFragment
            Generated modifier fragment with a resolved reactive atom.

        Returns
        -------
        ReactionProduct
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
        """Return legacy LYX-NHX linkage names without protein leaving atoms.

        This helper has no protein PDB context and therefore cannot know the
        actual lysine hydrogens present in a source structure. Pablo/product
        workflows must use :meth:`resolve_plan` and its resolved crosslink
        requirement instead.
        """
        contract = self.generic_contract(modifier)
        return ModifierLinkageSpec(
            protein_residue_name=contract.protein_endpoint.product_residue_name,
            modifier_residue_name=contract.modifier_endpoint.product_residue_name,
            protein_atom_name=contract.protein_endpoint.selector.atom_name,
            modifier_atom_name=contract.modifier_endpoint.selector.atom_name,
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
    ) -> ReactionProduct:
        """Resolve the NHS-Lys helper through the generic PDB contract."""
        return resolve_explicit_linkage_contract(
            protein_pdb_path,
            modifier,
            self.generic_contract(modifier, protein_pdb_path=protein_pdb_path),
            fragment=PreparedFragment.from_generated_fragment(
                modifier,
                source_identity=modifier.name,
                source_kind="polymer",
            ),
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
    from polyzymd.builders.conjugation.structure.pdb import load_pdb_as_rdkit_input

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


def _identity_for_atom(identities: Sequence[Any], atom: PdbAtomRecord) -> Any:
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
    target_identity: Any,
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


def _same_identity_residue(left: Any, right: Any) -> bool:
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
    return _parse_atoms(path, require_atoms=True)


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
