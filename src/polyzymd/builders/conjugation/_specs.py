"""Internal resolved attachment build specifications.

These records bridge reaction/linkage resolution and the current construction
functions without changing the public conjugation API.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation._linkage import ResolvedAttachmentPlan
from polyzymd.builders.conjugation.polymer.fragment import (
    GeneratedPolymerFragment,
    PolymerFragmentAtom,
    PolymerFragmentResidue,
)
from polyzymd.builders.conjugation.polymer.moiety import GeneratedMoietyFragment


class ConjugationFragment(BaseModel):
    """Generic internal fragment shape for resolved conjugation attachments."""

    model_config = {"arbitrary_types_allowed": True}

    atoms: tuple[PolymerFragmentAtom, ...] = Field(..., min_length=1)
    residues: tuple[PolymerFragmentResidue, ...] = Field(default_factory=tuple)
    bonds: tuple[tuple[int | str, int | str], ...] = Field(default_factory=tuple)
    bond_orders: tuple[tuple[int | str, int | str, float], ...] = Field(default_factory=tuple)
    sequence: str | None = None
    name: str = "fragment"
    source_kind: Literal["moiety", "polymer"] = "polymer"
    reactive_atom_serial: int | None = None
    reactive_atom_index: int | None = None
    reactive_atom_name: str | None = None
    leaving_atom_serials: tuple[int, ...] = Field(default_factory=tuple)
    leaving_atom_indices: tuple[int, ...] = Field(default_factory=tuple)
    leaving_atom_names: tuple[str, ...] = Field(default_factory=tuple)
    sidecars: dict[str, Path] = Field(default_factory=dict)
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)

    @classmethod
    def from_generated_polymer_fragment(
        cls,
        fragment: GeneratedPolymerFragment,
        *,
        source_kind: Literal["moiety", "polymer"] = "polymer",
        sidecars: dict[str, Path] | None = None,
        diagnostics: tuple[str, ...] = (),
    ) -> "ConjugationFragment":
        """Create a generic fragment from the existing polymer fragment model."""
        return cls(
            atoms=fragment.atoms,
            residues=fragment.residues,
            bonds=fragment.bonds,
            bond_orders=fragment.bond_orders,
            sequence=fragment.sequence,
            name=fragment.name,
            source_kind=source_kind,
            reactive_atom_serial=fragment.reactive_atom_serial,
            reactive_atom_index=fragment.reactive_atom_index,
            reactive_atom_name=fragment.reactive_atom_name,
            leaving_atom_serials=fragment.leaving_atom_serials,
            leaving_atom_indices=fragment.leaving_atom_indices,
            leaving_atom_names=fragment.leaving_atom_names,
            sidecars=sidecars or {},
            diagnostics=diagnostics,
        )

    def to_generated_polymer_fragment(self) -> GeneratedPolymerFragment:
        """Convert to the existing construction fragment adapter."""
        return GeneratedPolymerFragment(
            atoms=self.atoms,
            bonds=self.bonds,
            bond_orders=self.bond_orders,
            residues=self.residues,
            sequence=self.sequence,
            reactive_atom_serial=self.reactive_atom_serial,
            reactive_atom_index=self.reactive_atom_index,
            reactive_atom_name=self.reactive_atom_name,
            leaving_atom_serials=self.leaving_atom_serials,
            leaving_atom_indices=self.leaving_atom_indices,
            leaving_atom_names=self.leaving_atom_names,
            name=self.name,
        )


class AttachmentBuildSpec(BaseModel):
    """Resolved build input for one enabled conjugation attachment."""

    model_config = {"arbitrary_types_allowed": True}

    attachment_id: str
    attachment_index: int = Field(..., ge=1)
    attachment_config: Any = Field(exclude=True)
    reaction_name: str
    fragment: ConjugationFragment
    source_fragment: Any | None = Field(default=None, exclude=True)
    generated_fragment: GeneratedPolymerFragment = Field(exclude=True)
    resolved_plan: ResolvedAttachmentPlan
    source_sidecars: dict[str, Path] = Field(default_factory=dict)
    product_residue_mappings: dict[str, dict[str, int | str]] = Field(
        default_factory=dict,
        exclude=True,
    )
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)


def attachment_spec_from_moiety_plan(
    fragment: GeneratedMoietyFragment,
    plan: ResolvedAttachmentPlan,
    *,
    attachment_config: Any,
    attachment_index: int,
    reaction_name: str,
) -> AttachmentBuildSpec:
    """Adapt a one-residue moiety and resolved plan into a build spec."""
    generated_fragment = _generated_fragment_from_moiety_plan(fragment, plan)
    sidecars = _moiety_sidecars(fragment)
    generic_fragment = ConjugationFragment.from_generated_polymer_fragment(
        generated_fragment,
        source_kind="moiety",
        sidecars=sidecars,
        diagnostics=("Adapted GeneratedMoietyFragment to ConjugationFragment",),
    )
    return AttachmentBuildSpec(
        attachment_id=_attachment_id(attachment_config, attachment_index),
        attachment_index=attachment_index,
        attachment_config=attachment_config,
        reaction_name=reaction_name,
        fragment=generic_fragment,
        source_fragment=fragment,
        generated_fragment=generated_fragment,
        resolved_plan=plan,
        source_sidecars=sidecars,
        diagnostics=("Resolved moiety attachment build spec",),
    )


def attachment_spec_from_generated_polymer_plan(
    fragment: GeneratedPolymerFragment,
    sdf_path: Path | str | None,
    plan: ResolvedAttachmentPlan,
    *,
    attachment_config: Any,
    attachment_index: int,
    reaction_name: str,
    charged_sdf_path: Path | str | None = None,
    source_kind: Literal["moiety", "polymer"] = "polymer",
    source_fragment: Any | None = None,
    sidecars: dict[str, Path] | None = None,
) -> AttachmentBuildSpec:
    """Adapt a generated polymer fragment and resolved plan into a build spec."""
    resolved_sidecars = dict(sidecars or {})
    if sdf_path is not None:
        resolved_sidecars["sdf"] = Path(sdf_path)
        resolved_sidecars["bond_sdf"] = Path(sdf_path)
    if charged_sdf_path is not None:
        resolved_sidecars["charged_sdf"] = Path(charged_sdf_path)
    generic_fragment = ConjugationFragment.from_generated_polymer_fragment(
        fragment,
        source_kind=source_kind,
        sidecars=resolved_sidecars,
        diagnostics=("Adapted GeneratedPolymerFragment to ConjugationFragment",),
    )
    return AttachmentBuildSpec(
        attachment_id=_attachment_id(attachment_config, attachment_index),
        attachment_index=attachment_index,
        attachment_config=attachment_config,
        reaction_name=reaction_name,
        fragment=generic_fragment,
        source_fragment=source_fragment or fragment,
        generated_fragment=fragment,
        resolved_plan=plan,
        source_sidecars=resolved_sidecars,
        diagnostics=(f"Resolved {source_kind} attachment build spec",),
    )


def _generated_fragment_from_moiety_plan(
    fragment: GeneratedMoietyFragment,
    plan: ResolvedAttachmentPlan,
) -> GeneratedPolymerFragment:
    """Adapt a one-residue moiety fragment to the existing placement fragment model."""
    return GeneratedPolymerFragment(
        atoms=fragment.atoms,
        bonds=fragment.bonds,
        bond_orders=fragment.bond_orders,
        residues=fragment.residues,
        sequence=None,
        reactive_atom_serial=plan.modifier_link_atom.serial,
        reactive_atom_index=plan.modifier_link_atom.atom_index,
        reactive_atom_name=None,
        leaving_atom_serials=tuple(
            atom.serial for atom in plan.modifier_leaving_atoms if atom.serial is not None
        ),
        leaving_atom_indices=tuple(
            atom.atom_index for atom in plan.modifier_leaving_atoms if atom.atom_index is not None
        ),
        leaving_atom_names=(),
        name=fragment.name,
    )


def _moiety_sidecars(fragment: GeneratedMoietyFragment) -> dict[str, Path]:
    sidecars = {}
    if fragment.pdb_path is not None:
        sidecars["pdb"] = fragment.pdb_path
    if fragment.sdf_path is not None:
        sidecars["sdf"] = fragment.sdf_path
    return sidecars


def _attachment_id(attachment_config: Any, attachment_index: int) -> str:
    name = str(getattr(attachment_config, "name", "") or "").strip()
    return name or f"attachment_{attachment_index:02d}"
