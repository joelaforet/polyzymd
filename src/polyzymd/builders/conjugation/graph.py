"""Typed graph edit result models for covalent modification primitives."""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field


class AddedBond(BaseModel):
    """Atom-index record for a bond created by graph surgery."""

    begin_atom_index: int = Field(..., ge=0, description="Begin atom index in the product molecule")
    end_atom_index: int = Field(..., ge=0, description="End atom index in the product molecule")
    order: int = Field(1, ge=1, description="Integer bond order used by the primitive")


class RdkitGraphEditResult(BaseModel):
    """Structured result from an RDKit graph edit primitive."""

    model_config = {"arbitrary_types_allowed": True}

    product_mol: Any = Field(..., exclude=True, description="Product RDKit molecule")
    removed_protein_atom_indices: tuple[int, ...] = Field(
        default_factory=tuple,
        description="Original protein atom indices removed during graph surgery",
    )
    removed_moiety_atom_indices: tuple[int, ...] = Field(
        default_factory=tuple,
        description="Original moiety atom indices removed during graph surgery",
    )
    added_bond: AddedBond = Field(..., description="Product atom indices for the created bond")
    protein_atom_index_map: dict[int, int] = Field(
        default_factory=dict,
        description="Original-to-product atom index map for retained protein atoms",
    )
    moiety_atom_index_map: dict[int, int] = Field(
        default_factory=dict,
        description="Original-to-product atom index map for retained moiety atoms",
    )
    warnings: tuple[str, ...] = Field(
        default_factory=tuple,
        description="Non-fatal warnings emitted by the conservative executor",
    )
