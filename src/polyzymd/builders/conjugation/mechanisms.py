"""Declarative reaction mechanism models for covalent modification planning."""

from __future__ import annotations

from pydantic import BaseModel, Field, field_validator, model_validator


def _normalize_atom_label(value: str) -> str:
    """Normalize an atom label used by declarative graph edits."""
    normalized = value.strip().upper()
    if not normalized:
        raise ValueError("Atom labels must be non-empty strings")
    return normalized


def _normalize_identifier(value: str) -> str:
    """Normalize a mechanism identifier."""
    normalized = value.strip().lower()
    if not normalized:
        raise ValueError("Mechanism identifiers must be non-empty strings")
    return normalized


def _reject_duplicate_values(values: list[str], field_name: str) -> None:
    """Reject duplicate labels while preserving readable error messages."""
    duplicates = sorted({value for value in values if values.count(value) > 1})
    if duplicates:
        duplicate_text = ", ".join(duplicates)
        raise ValueError(f"Duplicate {field_name} labels are not allowed: {duplicate_text}")


class SiteAtomRule(BaseModel):
    """Allowed biomolecular attachment site for a declarative mechanism."""

    residue_name: str = Field(..., min_length=1, description="Allowed residue name")
    atom_name: str = Field(..., min_length=1, description="Allowed reactive atom name")
    role: str = Field("nucleophile", description="Chemical role of the site atom")
    rationale: str | None = Field(None, description="Human-readable rationale for this rule")

    @field_validator("residue_name", "atom_name")
    @classmethod
    def normalize_labels(cls, value: str) -> str:
        """Normalize residue and atom names for deterministic matching."""
        return _normalize_atom_label(value)


class MoietyReactiveGroup(BaseModel):
    """Reactive moiety labels needed for a future graph edit."""

    group: str = Field(..., min_length=1, description="Reactive group identifier")
    anchor_atom: str = Field(..., min_length=1, description="Moiety atom bonded to the site")
    leaving_atoms: list[str] = Field(
        default_factory=list,
        description="Moiety atom labels removed during future graph surgery",
    )
    rationale: str | None = Field(None, description="Human-readable rationale for this group")

    @field_validator("group")
    @classmethod
    def normalize_group(cls, value: str) -> str:
        """Normalize a reactive group identifier."""
        normalized = value.strip().lower()
        if not normalized:
            raise ValueError("Reactive group identifiers must be non-empty strings")
        return normalized

    @field_validator("anchor_atom")
    @classmethod
    def normalize_anchor_atom(cls, value: str) -> str:
        """Normalize the moiety anchor atom label."""
        return _normalize_atom_label(value)

    @field_validator("leaving_atoms")
    @classmethod
    def normalize_leaving_atoms(cls, values: list[str]) -> list[str]:
        """Normalize leaving atom labels and reject duplicates."""
        normalized = [_normalize_atom_label(value) for value in values]
        _reject_duplicate_values(normalized, "leaving atom")
        return normalized

    @model_validator(mode="after")
    def validate_anchor_not_leaving(self) -> "MoietyReactiveGroup":
        """Prevent the future bond anchor from also being removed."""
        if self.anchor_atom in self.leaving_atoms:
            raise ValueError("Moiety anchor atom cannot also be listed as a leaving atom")
        return self


class BondSpec(BaseModel):
    """Declarative covalent bond to add during a future graph edit."""

    site_atom: str = Field(..., min_length=1, description="Biomolecular atom label")
    moiety_atom: str = Field(..., min_length=1, description="Moiety atom label")
    order: int = Field(1, ge=1, description="Placeholder bond order")

    @field_validator("site_atom", "moiety_atom")
    @classmethod
    def normalize_atom_names(cls, value: str) -> str:
        """Normalize atom labels for duplicate bond detection."""
        return _normalize_atom_label(value)

    def key(self) -> tuple[str, str, int]:
        """Return a stable key for duplicate detection."""
        return (self.site_atom, self.moiety_atom, self.order)


class GraphEditPlan(BaseModel):
    """Non-executable graph edit placeholders for a reaction mechanism."""

    delete_site_atoms: list[str] = Field(
        default_factory=list,
        description="Biomolecular atom labels removed in future graph surgery",
    )
    delete_moiety_atoms: list[str] = Field(
        default_factory=list,
        description="Moiety atom labels removed in future graph surgery",
    )
    remove_site_hydrogens: list[str] = Field(
        default_factory=list,
        description="Site hydrogen labels removed in future proton handling",
    )
    remove_moiety_protons: list[str] = Field(
        default_factory=list,
        description="Moiety proton labels removed in future proton handling",
    )
    add_bonds: list[BondSpec] = Field(
        default_factory=list,
        description="Covalent bonds added in future graph surgery",
    )

    @field_validator(
        "delete_site_atoms",
        "delete_moiety_atoms",
        "remove_site_hydrogens",
        "remove_moiety_protons",
    )
    @classmethod
    def normalize_atom_lists(cls, values: list[str]) -> list[str]:
        """Normalize graph edit atom labels and reject duplicates."""
        normalized = [_normalize_atom_label(value) for value in values]
        _reject_duplicate_values(normalized, "graph edit atom")
        return normalized

    @model_validator(mode="after")
    def validate_unique_bonds(self) -> "GraphEditPlan":
        """Reject duplicate future bond specifications."""
        bond_keys = [bond.key() for bond in self.add_bonds]
        duplicates = sorted({key for key in bond_keys if bond_keys.count(key) > 1})
        if duplicates:
            duplicate_text = ", ".join(
                f"{site}-{moiety}:{order}" for site, moiety, order in duplicates
            )
            raise ValueError(f"Duplicate bond specifications are not allowed: {duplicate_text}")
        return self


class ChargePatchHint(BaseModel):
    """Placeholder for future local charge patching around a covalent junction."""

    strategy: str = Field("defer", min_length=1, description="Charge patch strategy hint")
    rationale: str | None = Field(None, description="Human-readable charge patch rationale")
    patch_radius_bonds: int | None = Field(
        None,
        ge=0,
        description="Optional local patch radius for future charge workflows",
    )

    @field_validator("strategy")
    @classmethod
    def normalize_strategy(cls, value: str) -> str:
        """Normalize a charge patch strategy label."""
        normalized = value.strip().lower()
        if not normalized:
            raise ValueError("Charge patch strategy must be a non-empty string")
        return normalized


class ReactionMechanism(BaseModel):
    """Serializable declaration of a named covalent modification mechanism."""

    identifier: str = Field(..., min_length=1, description="Stable mechanism identifier")
    display_name: str = Field(..., min_length=1, description="Human-readable mechanism name")
    allowed_sites: list[SiteAtomRule] = Field(
        ..., min_length=1, description="Allowed biomolecular site residue and atom rules"
    )
    moiety_reactive_group: MoietyReactiveGroup = Field(
        ..., description="Expected reactive group labels on the moiety"
    )
    graph_edits: GraphEditPlan = Field(
        default_factory=GraphEditPlan,
        description="Non-executable graph edit placeholders",
    )
    charge_patch_hint: ChargePatchHint = Field(
        default_factory=ChargePatchHint,
        description="Future charge patching hint",
    )
    rationale: str = Field(..., min_length=1, description="Mechanism rationale and intended use")
    notes: list[str] = Field(default_factory=list, description="Additional implementation notes")

    @field_validator("identifier")
    @classmethod
    def normalize_identifier(cls, value: str) -> str:
        """Normalize a mechanism identifier."""
        return _normalize_identifier(value)

    @model_validator(mode="after")
    def validate_unique_site_rules(self) -> "ReactionMechanism":
        """Reject duplicate allowed site residue and atom rules."""
        keys = [(rule.residue_name, rule.atom_name) for rule in self.allowed_sites]
        duplicates = sorted({key for key in keys if keys.count(key) > 1})
        if duplicates:
            duplicate_text = ", ".join(f"{residue}:{atom}" for residue, atom in duplicates)
            raise ValueError(f"Duplicate allowed site rules are not allowed: {duplicate_text}")
        return self

    @property
    def name(self) -> str:
        """Return the stable mechanism identifier for config compatibility."""
        return self.identifier
