"""Pablo residue-definition reaction exploration helpers."""

from __future__ import annotations

import importlib
from collections.abc import Sequence
from typing import Any

from pydantic import BaseModel, Field, field_validator


class PabloProductDefinitionDiagnostic(BaseModel):
    """Compact summary of one Pablo product residue definition."""

    residue_name: str | None = None
    atom_names: tuple[str, ...] = Field(default_factory=tuple)
    leaving_atom_names: tuple[str, ...] = Field(default_factory=tuple)
    linking_bond: str | None = None
    crosslink: str | None = None


class PabloReactionDiagnostic(BaseModel):
    """Diagnostic result for a Pablo ``ResidueDefinition.react`` exploration."""

    success: bool
    pablo_available: bool
    reactant_smarts: tuple[str, ...]
    product_smarts: tuple[str, ...]
    product_residue_names: tuple[str | None, ...]
    product_set_count: int = 0
    products: tuple[tuple[PabloProductDefinitionDiagnostic, ...], ...] = Field(
        default_factory=tuple
    )
    warnings: tuple[str, ...] = Field(default_factory=tuple)
    error_type: str | None = None
    error_message: str | None = None


class PabloReactionRequest(BaseModel):
    """Serializable inputs for Pablo residue-definition reaction exploration."""

    reactant_smarts: tuple[str, ...] = Field(..., min_length=1)
    product_smarts: tuple[str, ...] = Field(..., min_length=1)
    product_residue_names: tuple[str | None, ...] = Field(..., min_length=1)
    product_descriptions: tuple[str | None, ...] | None = None

    @field_validator("reactant_smarts", "product_smarts")
    @classmethod
    def validate_smarts(cls, values: tuple[str, ...]) -> tuple[str, ...]:
        """Reject empty SMARTS entries."""
        normalized = tuple(value.strip() for value in values)
        if any(not value for value in normalized):
            raise ValueError("SMARTS entries must be non-empty")
        return normalized

    @field_validator("product_residue_names")
    @classmethod
    def validate_product_names(cls, values: tuple[str | None, ...]) -> tuple[str | None, ...]:
        """Normalize product residue names without hardcoding chemistry."""
        normalized: list[str | None] = []
        for value in values:
            if value is None:
                normalized.append(None)
                continue
            residue_name = value.strip().upper()
            if not residue_name:
                raise ValueError("Product residue names must be non-empty when provided")
            normalized.append(residue_name)
        return tuple(normalized)


def explore_pablo_residue_reaction(
    reactants: Sequence[Any],
    request: PabloReactionRequest,
    *,
    pablo_module: Any | None = None,
    raise_on_error: bool = False,
) -> PabloReactionDiagnostic:
    """Run Pablo ``ResidueDefinition.react`` and summarize returned products.

    The helper is intentionally diagnostic: it does not claim the returned
    definitions are production-ready force-field templates. It exposes what
    Pablo can produce so PolyzyMD can map those definitions back to exact PDB
    atom identities in later product-state workflows.
    """
    try:
        pablo = pablo_module or importlib.import_module("openff.pablo")
    except Exception as exc:
        if raise_on_error:
            raise
        return PabloReactionDiagnostic(
            success=False,
            pablo_available=False,
            reactant_smarts=request.reactant_smarts,
            product_smarts=request.product_smarts,
            product_residue_names=request.product_residue_names,
            warnings=("OpenFF Pablo is not importable; residue reaction exploration skipped",),
            error_type=type(exc).__name__,
            error_message=str(exc),
        )

    try:
        product_sets = pablo.ResidueDefinition.react(
            reactants=reactants,
            reactant_smarts=request.reactant_smarts,
            product_smarts=request.product_smarts,
            product_residue_names=request.product_residue_names,
            product_descriptions=request.product_descriptions,
        )
    except Exception as exc:
        if raise_on_error:
            raise
        return PabloReactionDiagnostic(
            success=False,
            pablo_available=True,
            reactant_smarts=request.reactant_smarts,
            product_smarts=request.product_smarts,
            product_residue_names=request.product_residue_names,
            warnings=(
                "Pablo imported, but ResidueDefinition.react did not return product definitions",
            ),
            error_type=type(exc).__name__,
            error_message=str(exc),
        )

    products = tuple(
        tuple(_summarize_definition(definition) for definition in row) for row in product_sets
    )
    warnings = []
    if not products:
        warnings.append("Pablo returned no product definition sets")
    if len(request.product_residue_names) != len(request.product_smarts):
        warnings.append("Product residue name count differs from product SMARTS count")
    return PabloReactionDiagnostic(
        success=True,
        pablo_available=True,
        reactant_smarts=request.reactant_smarts,
        product_smarts=request.product_smarts,
        product_residue_names=request.product_residue_names,
        product_set_count=len(products),
        products=products,
        warnings=tuple(warnings),
    )


def _summarize_definition(definition: Any) -> PabloProductDefinitionDiagnostic:
    """Extract stable diagnostics from a Pablo residue definition-like object."""
    atoms = tuple(getattr(definition, "atoms", ()) or ())
    atom_names = tuple(
        name
        for atom in atoms
        if (name := getattr(atom, "name", getattr(atom, "atom_name", None))) is not None
    )
    leaving_atom_names = tuple(
        name
        for atom in atoms
        if getattr(atom, "leaving", False)
        if (name := getattr(atom, "name", getattr(atom, "atom_name", None))) is not None
    )
    return PabloProductDefinitionDiagnostic(
        residue_name=getattr(definition, "residue_name", getattr(definition, "name", None)),
        atom_names=atom_names,
        leaving_atom_names=leaving_atom_names,
        linking_bond=_metadata_repr(getattr(definition, "linking_bond", None)),
        crosslink=_metadata_repr(getattr(definition, "crosslink", None)),
    )


def _metadata_repr(value: Any) -> str | None:
    """Return a compact representation for optional Pablo bond metadata."""
    if value is None:
        return None
    return repr(value)


__all__ = [
    "PabloProductDefinitionDiagnostic",
    "PabloReactionDiagnostic",
    "PabloReactionRequest",
    "explore_pablo_residue_reaction",
]
