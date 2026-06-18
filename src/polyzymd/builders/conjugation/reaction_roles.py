"""Generic atom-mapped reaction role helpers for covalent conjugation POCs."""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from typing import Any, Literal

from pydantic import BaseModel, Field, field_validator, model_validator

ParticipantRole = Literal["site", "moiety"]
AtomRole = Literal["linking", "leaving", "retained", "geometry_anchor"]
BondChangeKind = Literal["added", "removed", "order_changed"]
STRUCTURE_MATCHING_BLOCKER_MESSAGE = (
    "Generic atom-mapped coordinate surgery is not implemented: this preflight only derives "
    "atom-map roles and bond changes from SMARTS plus optional caller-supplied identity maps. "
    "TODO: implement structure matching from selected protein and modifier atoms before using "
    "generic mechanisms for coordinate edits."
)


class ReactionParticipant(BaseModel):
    """Reactant participant in an atom-mapped covalent mechanism."""

    name: str = Field(..., min_length=1, description="Human-readable participant name")
    role: ParticipantRole = Field(..., description="Participant role in the reaction")
    reactant_index: int | None = Field(
        None,
        ge=0,
        description="Optional index into the Pablo/RDKit reactant SMARTS list",
    )


class PdbAtomIdentity(BaseModel):
    """Concrete PDB atom identity resolved for an atom-map number."""

    chain_id: str = Field(..., min_length=1, description="PDB chain identifier")
    residue_name: str = Field(..., min_length=1, description="PDB residue name")
    residue_number: int = Field(..., description="PDB residue sequence number")
    atom_name: str = Field(..., min_length=1, description="PDB atom name")
    serial: int | None = Field(None, ge=1, description="Optional PDB atom serial")
    insertion_code: str | None = Field(None, description="Optional PDB insertion code")

    @field_validator("chain_id", "residue_name", "atom_name", "insertion_code")
    @classmethod
    def normalize_labels(cls, value: str | None) -> str | None:
        """Normalize fixed-width PDB labels without assuming residue chemistry."""
        if value is None:
            return None
        normalized = value.strip()
        if not normalized:
            raise ValueError("PDB atom identity labels must be non-empty")
        return normalized.upper()


class AtomRoleSpec(BaseModel):
    """Role assigned to one atom-map number by reaction metadata."""

    map_number: int = Field(..., ge=1, description="Atom-map number in reaction SMARTS")
    participant: ParticipantRole = Field(..., description="Participant owning the mapped atom")
    role: AtomRole = Field(..., description="Role of this atom in the reaction")
    label: str | None = Field(None, description="Optional mechanism-local label")
    required: bool = Field(True, description="Whether resolution must find a concrete atom")

    @field_validator("label")
    @classmethod
    def normalize_label(cls, value: str | None) -> str | None:
        """Normalize optional role labels while keeping them chemistry-agnostic."""
        if value is None:
            return None
        normalized = value.strip()
        if not normalized:
            raise ValueError("Role labels must be non-empty when provided")
        return normalized


class MappedBond(BaseModel):
    """Bond between two atom-map numbers."""

    atom_maps: tuple[int, int] = Field(..., description="Sorted atom-map numbers")
    order: float = Field(..., ge=0, description="SMARTS-derived bond order; 0 means query/any")

    @field_validator("atom_maps")
    @classmethod
    def normalize_atom_maps(cls, value: tuple[int, int]) -> tuple[int, int]:
        """Store map pairs deterministically."""
        if len(value) != 2:
            raise ValueError("Mapped bonds require exactly two atom-map numbers")
        a, b = value
        if a < 1 or b < 1:
            raise ValueError("Atom-map numbers must be positive")
        if a == b:
            raise ValueError("Mapped bonds cannot connect an atom to itself")
        return tuple(sorted((a, b)))

    @property
    def key(self) -> tuple[int, int]:
        """Return the atom-map pair as a dictionary key."""
        return self.atom_maps


class BondChange(BaseModel):
    """Mapped bond change between reactant and product SMARTS."""

    kind: BondChangeKind = Field(..., description="Type of mapped bond change")
    atom_maps: tuple[int, int] = Field(..., description="Sorted changed atom-map pair")
    reactant_order: float | None = Field(None, description="Reactant bond order, when present")
    product_order: float | None = Field(None, description="Product bond order, when present")

    @field_validator("atom_maps")
    @classmethod
    def normalize_atom_maps(cls, value: tuple[int, int]) -> tuple[int, int]:
        """Store changed map pairs deterministically."""
        return MappedBond(atom_maps=value, order=1).atom_maps


class MappedSmartsValidation(BaseModel):
    """Validation summary for reactant and product atom-mapped SMARTS."""

    reactant_maps: tuple[int, ...]
    product_maps: tuple[int, ...]
    missing_product_maps: tuple[int, ...] = Field(
        default_factory=tuple,
        description="Reactant maps absent from products; usually leaving fragments",
    )
    new_product_maps: tuple[int, ...] = Field(
        default_factory=tuple,
        description="Product maps absent from reactants",
    )

    @property
    def is_balanced(self) -> bool:
        """Return whether product SMARTS retain every mapped reactant atom."""
        return not self.missing_product_maps and not self.new_product_maps


class BondChangeSet(BaseModel):
    """Mapped bond changes derived from atom-mapped reactant and product SMARTS."""

    added_bonds: tuple[BondChange, ...] = Field(default_factory=tuple)
    removed_bonds: tuple[BondChange, ...] = Field(default_factory=tuple)
    order_changes: tuple[BondChange, ...] = Field(default_factory=tuple)


class AtomMappedReaction(BaseModel):
    """POC-level atom-mapped covalent mechanism declaration."""

    name: str = Field(..., min_length=1, description="Mechanism name")
    reactant_smarts: tuple[str, ...] = Field(..., min_length=1)
    product_smarts: tuple[str, ...] = Field(..., min_length=1)
    participants: tuple[ReactionParticipant, ...] = Field(default_factory=tuple)
    atom_roles: tuple[AtomRoleSpec, ...] = Field(default_factory=tuple)
    description: str | None = None

    @field_validator("reactant_smarts", "product_smarts")
    @classmethod
    def validate_smarts_entries(cls, values: tuple[str, ...]) -> tuple[str, ...]:
        """Reject empty SMARTS fragments early."""
        normalized = tuple(value.strip() for value in values)
        if any(not value for value in normalized):
            raise ValueError("SMARTS entries must be non-empty")
        return normalized

    @model_validator(mode="after")
    def validate_roles_reference_maps(self) -> "AtomMappedReaction":
        """Ensure role metadata references mapped atoms in the declaration."""
        validation = validate_mapped_smarts(self.reactant_smarts, self.product_smarts)
        known_maps = set(validation.reactant_maps) | set(validation.product_maps)
        unknown = sorted(
            {role.map_number for role in self.atom_roles if role.map_number not in known_maps}
        )
        if unknown:
            raise ValueError(f"Atom roles reference unknown atom-map numbers: {unknown}")
        return self

    @classmethod
    def from_reaction_smarts(
        cls,
        name: str,
        reaction_smarts: str,
        **kwargs: Any,
    ) -> "AtomMappedReaction":
        """Create a declaration from ``reactants>>products`` SMARTS."""
        if ">>" not in reaction_smarts:
            raise ValueError("Reaction SMARTS must contain '>>'")
        reactants, products = reaction_smarts.split(">>", 1)
        return cls(
            name=name,
            reactant_smarts=tuple(_split_smarts_fragments(reactants)),
            product_smarts=tuple(_split_smarts_fragments(products)),
            **kwargs,
        )


class ResolvedAtomRole(BaseModel):
    """Resolved atom role with optional concrete PDB identity."""

    spec: AtomRoleSpec
    pdb_identity: PdbAtomIdentity | None = None


class ResolvedBondChange(BaseModel):
    """Bond change resolved to concrete identities when available."""

    change: BondChange
    atom_identities: tuple[PdbAtomIdentity | None, PdbAtomIdentity | None]


class GenericResolvedReactionPlan(BaseModel):
    """POC resolved reaction plan from mapped SMARTS plus atom identities."""

    reaction: AtomMappedReaction
    validation: MappedSmartsValidation
    bond_changes: BondChangeSet
    resolved_roles: tuple[ResolvedAtomRole, ...]
    resolved_added_bonds: tuple[ResolvedBondChange, ...]
    resolved_removed_bonds: tuple[ResolvedBondChange, ...]
    resolved_order_changes: tuple[ResolvedBondChange, ...]
    inferred_leaving_maps: tuple[int, ...]


def validate_mapped_smarts(
    reactant_smarts: Sequence[str],
    product_smarts: Sequence[str],
    *,
    require_all_reactant_maps_in_products: bool = False,
) -> MappedSmartsValidation:
    """Validate and summarize atom maps in reactant/product SMARTS.

    Parameters
    ----------
    reactant_smarts : Sequence[str]
        Atom-mapped reactant SMARTS fragments.
    product_smarts : Sequence[str]
        Atom-mapped product SMARTS fragments.
    require_all_reactant_maps_in_products : bool, optional
        Pablo requires every reactant map to appear in product SMARTS. PolyzyMD's
        generic POC path permits missing product maps so leaving groups can be
        represented explicitly, by default ``False``.
    """
    reactant_maps = _mapped_atom_numbers(reactant_smarts)
    product_maps = _mapped_atom_numbers(product_smarts)
    if not reactant_maps:
        raise ValueError("Reactant SMARTS must contain atom-map numbers")
    if not product_maps:
        raise ValueError("Product SMARTS must contain atom-map numbers")

    missing = tuple(sorted(set(reactant_maps) - set(product_maps)))
    new = tuple(sorted(set(product_maps) - set(reactant_maps)))
    if require_all_reactant_maps_in_products and missing:
        raise ValueError(f"Product SMARTS are missing reactant atom-map numbers: {missing}")
    return MappedSmartsValidation(
        reactant_maps=tuple(sorted(reactant_maps)),
        product_maps=tuple(sorted(product_maps)),
        missing_product_maps=missing,
        new_product_maps=new,
    )


def derive_bond_changes(
    reactant_smarts: Sequence[str],
    product_smarts: Sequence[str],
) -> BondChangeSet:
    """Derive added, removed, and order-changed bonds between mapped atoms."""
    reactant_bonds = _mapped_bond_orders(reactant_smarts)
    product_bonds = _mapped_bond_orders(product_smarts)

    reactant_keys = set(reactant_bonds)
    product_keys = set(product_bonds)
    added = tuple(
        BondChange(kind="added", atom_maps=key, product_order=product_bonds[key])
        for key in sorted(product_keys - reactant_keys)
    )
    removed = tuple(
        BondChange(kind="removed", atom_maps=key, reactant_order=reactant_bonds[key])
        for key in sorted(reactant_keys - product_keys)
    )
    order_changes = tuple(
        BondChange(
            kind="order_changed",
            atom_maps=key,
            reactant_order=reactant_bonds[key],
            product_order=product_bonds[key],
        )
        for key in sorted(reactant_keys & product_keys)
        if reactant_bonds[key] != product_bonds[key]
    )
    return BondChangeSet(added_bonds=added, removed_bonds=removed, order_changes=order_changes)


def atom_mapped_reaction_from_mechanism_config(mechanism: Any) -> AtomMappedReaction:
    """Build an atom-mapped reaction declaration from a mechanism config-like object.

    This is a schema-to-reaction adapter only. It validates SMARTS, role metadata,
    and bond-change derivation inputs, but it does not match those mapped atoms to
    a protein PDB or modifier structure.
    """
    reaction_smarts = getattr(mechanism, "reaction_smarts", None)
    if reaction_smarts is None:
        raise ValueError("mechanism.reaction_smarts is required for atom-mapped reactions")
    atom_roles = tuple(_coerce_atom_role(role) for role in getattr(mechanism, "atom_roles", ()))
    participants = _participants_from_atom_roles(atom_roles)
    return AtomMappedReaction.from_reaction_smarts(
        name=str(getattr(mechanism, "name", "atom_mapped_reaction")),
        reaction_smarts=str(reaction_smarts),
        participants=participants,
        atom_roles=atom_roles,
    )


def resolve_reaction_roles_from_identity_map(
    reaction: AtomMappedReaction,
    map_number_to_pdb_identity: Mapping[int, PdbAtomIdentity | Mapping[str, Any]],
    *,
    require_required_identities: bool = True,
) -> GenericResolvedReactionPlan:
    """Resolve roles using a caller-supplied atom-map to PDB-identity map.

    This function is intentionally not a structure matcher. It does not inspect
    protein coordinates, modifier coordinates, or RDKit substructure matches; the
    caller must supply any concrete atom identities. Use
    ``STRUCTURE_MATCHING_BLOCKER_MESSAGE`` when a workflow reaches this boundary
    without a mechanism-specific coordinate backend.
    """
    identities = {
        int(map_number): _coerce_identity(identity)
        for map_number, identity in map_number_to_pdb_identity.items()
    }
    validation = validate_mapped_smarts(reaction.reactant_smarts, reaction.product_smarts)
    bond_changes = derive_bond_changes(reaction.reactant_smarts, reaction.product_smarts)
    inferred_leaving = tuple(
        sorted(
            set(validation.missing_product_maps)
            | {role.map_number for role in reaction.atom_roles if role.role == "leaving"}
        )
    )

    missing_required = sorted(
        role.map_number
        for role in reaction.atom_roles
        if role.required and role.map_number not in identities
    )
    if require_required_identities and missing_required:
        raise ValueError(f"No PDB identity was provided for required atom maps: {missing_required}")

    resolved_roles = tuple(
        ResolvedAtomRole(spec=role, pdb_identity=identities.get(role.map_number))
        for role in reaction.atom_roles
    )
    return GenericResolvedReactionPlan(
        reaction=reaction,
        validation=validation,
        bond_changes=bond_changes,
        resolved_roles=resolved_roles,
        resolved_added_bonds=_resolve_bond_changes(bond_changes.added_bonds, identities),
        resolved_removed_bonds=_resolve_bond_changes(bond_changes.removed_bonds, identities),
        resolved_order_changes=_resolve_bond_changes(bond_changes.order_changes, identities),
        inferred_leaving_maps=inferred_leaving,
    )


def _coerce_atom_role(role: AtomRoleSpec | Mapping[str, Any] | Any) -> AtomRoleSpec:
    """Coerce config-style atom role metadata to the reaction role model."""
    if isinstance(role, AtomRoleSpec):
        return role
    if hasattr(role, "model_dump"):
        return AtomRoleSpec.model_validate(role.model_dump())
    return AtomRoleSpec.model_validate(role)


def _participants_from_atom_roles(atom_roles: Sequence[AtomRoleSpec]) -> tuple[ReactionParticipant, ...]:
    """Infer participant declarations from role metadata without chemistry assumptions."""
    participants: list[ReactionParticipant] = []
    seen: set[ParticipantRole] = set()
    for role in atom_roles:
        if role.participant in seen:
            continue
        seen.add(role.participant)
        participants.append(ReactionParticipant(name=role.participant, role=role.participant))
    return tuple(participants)


def _coerce_identity(identity: PdbAtomIdentity | Mapping[str, Any]) -> PdbAtomIdentity:
    """Coerce mapping-style identities to the typed model."""
    if isinstance(identity, PdbAtomIdentity):
        return identity
    return PdbAtomIdentity.model_validate(identity)


def _resolve_bond_changes(
    changes: Sequence[BondChange],
    identities: Mapping[int, PdbAtomIdentity],
) -> tuple[ResolvedBondChange, ...]:
    """Resolve changed mapped bonds to PDB identities where possible."""
    return tuple(
        ResolvedBondChange(
            change=change,
            atom_identities=(
                identities.get(change.atom_maps[0]),
                identities.get(change.atom_maps[1]),
            ),
        )
        for change in changes
    )


def _mapped_atom_numbers(smarts_entries: Sequence[str]) -> tuple[int, ...]:
    """Return all unique atom-map numbers in SMARTS entries."""
    maps: set[int] = set()
    for smarts in smarts_entries:
        maps.update(int(match) for match in re.findall(r":(\d+)(?=[^\]]*\])", smarts))
    return tuple(sorted(maps))


def _mapped_bond_orders(smarts_entries: Sequence[str]) -> dict[tuple[int, int], float]:
    """Return mapped bond orders, preferring RDKit and falling back to a light parser."""
    try:
        return _rdkit_mapped_bond_orders(smarts_entries)
    except Exception:
        return _fallback_mapped_bond_orders(smarts_entries)


def _rdkit_mapped_bond_orders(smarts_entries: Sequence[str]) -> dict[tuple[int, int], float]:
    """Use RDKit to extract mapped bond orders from SMARTS entries."""
    from rdkit import Chem

    bonds: dict[tuple[int, int], float] = {}
    for smarts in smarts_entries:
        mol = Chem.MolFromSmarts(smarts)
        if mol is None:
            raise ValueError(f"RDKit could not parse SMARTS: {smarts}")
        for bond in mol.GetBonds():
            begin = bond.GetBeginAtom().GetAtomMapNum()
            end = bond.GetEndAtom().GetAtomMapNum()
            if begin and end:
                key = tuple(sorted((int(begin), int(end))))
                bonds[key] = _rdkit_bond_order(bond)
    return bonds


def _rdkit_bond_order(bond: Any) -> float:
    """Convert an RDKit bond into a compact numeric order."""
    from rdkit import Chem

    bond_type = bond.GetBondType()
    if bond_type == Chem.BondType.UNSPECIFIED:
        return 0.0
    if bond_type == Chem.BondType.AROMATIC:
        return 1.5
    return float(bond.GetBondTypeAsDouble())


def _fallback_mapped_bond_orders(smarts_entries: Sequence[str]) -> dict[tuple[int, int], float]:
    """Small SMARTS parser for bracketed mapped atoms used by lightweight tests."""
    bonds: dict[tuple[int, int], float] = {}
    for smarts in smarts_entries:
        bonds.update(_fallback_fragment_bonds(smarts))
    return bonds


def _fallback_fragment_bonds(smarts: str) -> dict[tuple[int, int], float]:
    """Parse bracketed mapped atoms, branches, dots, and simple ring closures."""
    bonds: dict[tuple[int, int], float] = {}
    branch_stack: list[int | None] = []
    ring_open: dict[str, tuple[int, float]] = {}
    previous_map: int | None = None
    pending_order = 1.0
    index = 0
    while index < len(smarts):
        char = smarts[index]
        if char in "-=#:~":
            pending_order = _bond_symbol_order(char)
            index += 1
            continue
        if char == ".":
            previous_map = None
            pending_order = 1.0
            index += 1
            continue
        if char == "(":
            branch_stack.append(previous_map)
            index += 1
            continue
        if char == ")":
            previous_map = branch_stack.pop() if branch_stack else None
            pending_order = 1.0
            index += 1
            continue
        if char != "[":
            index += 1
            continue

        end = smarts.find("]", index)
        if end == -1:
            raise ValueError(f"Unclosed bracket atom in SMARTS: {smarts}")
        atom_text = smarts[index + 1 : end]
        match = re.search(r":(\d+)", atom_text)
        current_map = int(match.group(1)) if match else None
        if previous_map is not None and current_map is not None:
            bonds[tuple(sorted((previous_map, current_map)))] = pending_order
        previous_map = current_map
        pending_order = 1.0
        index = end + 1

        while index < len(smarts) and smarts[index].isdigit():
            ring_digit = smarts[index]
            if current_map is not None:
                if ring_digit in ring_open:
                    other_map, order = ring_open.pop(ring_digit)
                    bonds[tuple(sorted((current_map, other_map)))] = pending_order or order
                    pending_order = 1.0
                else:
                    ring_open[ring_digit] = (current_map, pending_order)
                    pending_order = 1.0
            index += 1
    return bonds


def _bond_symbol_order(symbol: str) -> float:
    """Map simple SMARTS bond symbols to numeric order."""
    return {"-": 1.0, "=": 2.0, "#": 3.0, ":": 1.5, "~": 0.0}[symbol]


def _split_smarts_fragments(smarts: str) -> list[str]:
    """Split dot-separated SMARTS fragments for simple reaction declarations."""
    return [fragment for fragment in (part.strip() for part in smarts.split(".")) if fragment]


__all__ = [
    "AtomMappedReaction",
    "AtomRoleSpec",
    "BondChange",
    "BondChangeSet",
    "GenericResolvedReactionPlan",
    "MappedBond",
    "MappedSmartsValidation",
    "PdbAtomIdentity",
    "ReactionParticipant",
    "ResolvedAtomRole",
    "ResolvedBondChange",
    "STRUCTURE_MATCHING_BLOCKER_MESSAGE",
    "atom_mapped_reaction_from_mechanism_config",
    "derive_bond_changes",
    "resolve_reaction_roles_from_identity_map",
    "validate_mapped_smarts",
]
