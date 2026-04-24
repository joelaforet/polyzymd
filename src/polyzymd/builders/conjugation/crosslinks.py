"""Validation helpers for explicit Pablo crosslink configuration."""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation.contracts import PabloCrosslinkRequirement
from polyzymd.builders.conjugation.linkers import ModifierLinkageSpec, ModifierLinker
from polyzymd.builders.conjugation.polymer_fragment import GeneratedPolymerFragment


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
