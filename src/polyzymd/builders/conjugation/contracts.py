"""Generic PDB linkage contract models and resolvers."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.builders.conjugation.pdb_assembly import (
    NhsLysPdbAttachment,
    PdbAtomRecord,
    PdbLinkageAttachment,
    PlacedPolymerFragment,
)
from polyzymd.builders.conjugation.polymer_fragment import GeneratedPolymerFragment

_ATOM_RECORD_PREFIXES = ("ATOM", "HETATM")


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


class ResolvedAttachmentPlan(BaseModel):
    """Resolved atom-level attachment plan for PDB assembly and Pablo."""

    contract: ExplicitLinkageContract
    protein_link_atom: PdbAtomRecord
    modifier_link_atom: PdbAtomRecord
    protein_leaving_atoms: tuple[PdbAtomRecord, ...] = Field(default_factory=tuple)
    modifier_leaving_atoms: tuple[PdbAtomRecord, ...] = Field(default_factory=tuple)
    protein_product_residue_name: str
    modifier_product_residue_name: str
    pablo_crosslink_requirement: PabloCrosslinkRequirement
    target_bond_length_angstrom: float = Field(1.33, gt=0)

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
    modifier: Path
    | str
    | GeneratedPolymerFragment
    | PlacedPolymerFragment
    | Sequence[PdbAtomRecord],
    contract: ExplicitLinkageContract,
) -> ResolvedAttachmentPlan:
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
    ResolvedAttachmentPlan
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

    requirement = PabloCrosslinkRequirement(
        residues=(
            contract.protein_endpoint.product_residue_name,
            contract.modifier_endpoint.product_residue_name,
        ),
        linking_atoms=(protein_link_atom.atom_name, modifier_link_atom.atom_name),
        leaving_atoms=(
            tuple(atom.atom_name for atom in protein_leaving_atoms),
            tuple(atom.atom_name for atom in modifier_leaving_atoms),
        ),
        bond_order=contract.bond.bond_order,
    )
    return ResolvedAttachmentPlan(
        contract=contract,
        protein_link_atom=protein_link_atom,
        modifier_link_atom=modifier_link_atom,
        protein_leaving_atoms=protein_leaving_atoms,
        modifier_leaving_atoms=modifier_leaving_atoms,
        protein_product_residue_name=contract.protein_endpoint.product_residue_name,
        modifier_product_residue_name=contract.modifier_endpoint.product_residue_name,
        pablo_crosslink_requirement=requirement,
        target_bond_length_angstrom=contract.bond.target_bond_length_angstrom,
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
    pdb_path = Path(path)
    atoms: list[PdbAtomRecord] = []
    with pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(_ATOM_RECORD_PREFIXES):
                atoms.append(PdbAtomRecord.from_pdb_line(line, atom_index=len(atoms)))
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records found in {pdb_path}")
    return tuple(atoms)


def placed_fragment_from_resolved_plan(
    fragment: PlacedPolymerFragment,
    plan: ResolvedAttachmentPlan,
) -> PlacedPolymerFragment:
    """Return a placed fragment with resolved generic linkage selectors.

    Parameters
    ----------
    fragment : PlacedPolymerFragment
        Existing placed fragment whose atom identities match the resolved plan.
    plan : ResolvedAttachmentPlan
        Resolved generic linkage plan.

    Returns
    -------
    PlacedPolymerFragment
        Fragment whose reactive and leaving selectors are scoped by resolved
        atom serials/indices instead of global names.
    """
    return fragment.model_copy(
        update={
            "reactive_atom_serial": plan.modifier_link_atom.serial,
            "reactive_atom_index": plan.modifier_link_atom.atom_index,
            "reactive_atom_name": None,
            "leaving_atom_serials": tuple(
                atom.serial for atom in plan.modifier_leaving_atoms if atom.serial is not None
            ),
            "leaving_atom_indices": tuple(
                atom.atom_index
                for atom in plan.modifier_leaving_atoms
                if atom.atom_index is not None
            ),
            "leaving_atom_names": (),
        }
    )


def _modifier_atom_records(
    modifier: Path
    | str
    | GeneratedPolymerFragment
    | PlacedPolymerFragment
    | Sequence[PdbAtomRecord],
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
