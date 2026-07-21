"""Reaction-product enrichment from prepared modifier sources."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from polyzymd.builders.conjugation.polymer.fragment import (
    GeneratedPolymerFragment,
    PreparedFragment,
)
from polyzymd.builders.conjugation.polymer.moiety import GeneratedMoietyFragment


def prepared_fragment_from_moiety(
    fragment: GeneratedMoietyFragment,
) -> PreparedFragment:
    """Prepare a one-residue moiety without reaction-derived selectors."""
    sidecars = _moiety_sidecars(fragment)
    generic_fragment = PreparedFragment(
        atoms=fragment.atoms,
        bonds=fragment.bonds,
        bond_orders=fragment.bond_orders,
        residues=fragment.residues,
        sequence=None,
        name=fragment.name,
        source_identity=str(fragment.sdf_path or fragment.pdb_path or fragment.name),
        source_kind="smiles",
        sidecars=sidecars,
        provenance={"residue_name": fragment.residue_name},
        diagnostics=("Prepared SMILES moiety fragment",),
    )
    return generic_fragment


def prepare_generated_fragment(
    fragment: GeneratedPolymerFragment,
    sdf_path: Path | str | None,
    *,
    charged_sdf_path: Path | str | None = None,
    source_kind: Literal["polymer", "smiles", "pdb_fragment"] = "polymer",
    sidecars: dict[str, Path] | None = None,
    provenance: dict[str, object] | None = None,
    diagnostics: tuple[str, ...] = (),
) -> PreparedFragment:
    """Prepare a generated fragment without reaction-derived selectors."""
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
            provenance=provenance,
            diagnostics=diagnostics or (f"Prepared {source_kind} fragment",),
        )
    )
    return generic_fragment


def prepare_reaction_fragment(fragment: Any) -> PreparedFragment:
    """Normalize a public reaction input into the construction fragment contract."""
    if isinstance(fragment, PreparedFragment):
        return fragment
    if isinstance(fragment, GeneratedMoietyFragment):
        return prepared_fragment_from_moiety(fragment)
    if isinstance(fragment, GeneratedPolymerFragment):
        return prepare_generated_fragment(fragment, None)
    raise TypeError(
        "Reaction fragments must be GeneratedMoietyFragment, "
        "GeneratedPolymerFragment, or PreparedFragment instances"
    )


def _moiety_sidecars(fragment: GeneratedMoietyFragment) -> dict[str, Path]:
    sidecars = {}
    if fragment.pdb_path is not None:
        sidecars["pdb"] = fragment.pdb_path
    if fragment.sdf_path is not None:
        sidecars["sdf"] = fragment.sdf_path
    return sidecars
