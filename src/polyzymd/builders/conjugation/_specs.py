"""Reaction-product enrichment from prepared modifier sources."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from polyzymd.builders.conjugation._linkage import ReactionProduct
from polyzymd.builders.conjugation.polymer.fragment import (
    GeneratedPolymerFragment,
    PreparedFragment,
)
from polyzymd.builders.conjugation.polymer.moiety import GeneratedMoietyFragment


def reaction_product_from_moiety(
    fragment: GeneratedMoietyFragment,
    plan: ReactionProduct,
    *,
    attachment_config: Any,
    attachment_index: int,
    reaction_name: str,
) -> ReactionProduct:
    """Adapt a one-residue moiety and resolved plan into a build spec."""
    generated_fragment = _generated_fragment_from_moiety_plan(fragment, plan)
    sidecars = _moiety_sidecars(fragment)
    generic_fragment = PreparedFragment.from_generated_fragment(
        generated_fragment,
        source_identity=str(fragment.sdf_path or fragment.pdb_path or fragment.name),
        source_kind="smiles",
        sidecars=sidecars,
        provenance={"residue_name": fragment.residue_name},
        diagnostics=("Prepared SMILES moiety fragment",),
    )
    return plan.model_copy(
        update={
            "attachment_id": _attachment_id(attachment_config, attachment_index),
            "attachment_index": attachment_index,
            "attachment_config": attachment_config,
            "reaction_name": reaction_name,
            "fragment": generic_fragment,
            "source_sidecars": sidecars,
            "attachment_force_field_domain": _attachment_force_field_domain(attachment_config),
            "diagnostics": (*plan.diagnostics, "Resolved moiety reaction product"),
        }
    )


def reaction_product_from_generated_fragment(
    fragment: GeneratedPolymerFragment,
    sdf_path: Path | str | None,
    plan: ReactionProduct,
    *,
    attachment_config: Any,
    attachment_index: int,
    reaction_name: str,
    charged_sdf_path: Path | str | None = None,
    source_kind: Literal["polymer", "smiles", "pdb_fragment"] = "polymer",
    sidecars: dict[str, Path] | None = None,
) -> ReactionProduct:
    """Adapt a generated polymer fragment and resolved plan into a build spec."""
    resolved_sidecars = dict(sidecars or {})
    if sdf_path is not None:
        resolved_sidecars["sdf"] = Path(sdf_path)
        resolved_sidecars["bond_sdf"] = Path(sdf_path)
    if charged_sdf_path is not None:
        resolved_sidecars["charged_sdf"] = Path(charged_sdf_path)
    generic_fragment = (
        fragment
        if isinstance(fragment, PreparedFragment)
        else PreparedFragment.from_generated_fragment(
            fragment,
            source_identity=str(
                resolved_sidecars.get("sdf") or resolved_sidecars.get("pdb") or fragment.name
            ),
            source_kind=source_kind,
            sidecars=resolved_sidecars,
            diagnostics=(f"Prepared {source_kind} fragment",),
        )
    )
    return plan.model_copy(
        update={
            "attachment_id": _attachment_id(attachment_config, attachment_index),
            "attachment_index": attachment_index,
            "attachment_config": attachment_config,
            "reaction_name": reaction_name,
            "fragment": generic_fragment,
            "source_sidecars": resolved_sidecars,
            "attachment_force_field_domain": _attachment_force_field_domain(attachment_config),
            "diagnostics": (*plan.diagnostics, f"Resolved {source_kind} reaction product"),
        }
    )


def _generated_fragment_from_moiety_plan(
    fragment: GeneratedMoietyFragment,
    plan: ReactionProduct,
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


def _attachment_force_field_domain(attachment_config: Any) -> str:
    """Return the configured modifier force-field domain for this product."""
    moiety = getattr(attachment_config, "moiety", None)
    return str(getattr(moiety, "force_field", "") or "").strip()
