"""Execution context models for explicit conjugation graph edits."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator

from polyzymd.builders.conjugation.graph import AddedBond, RdkitGraphEditResult
from polyzymd.builders.conjugation.nhs_lys import NhsLysGraphEditPlan, NhsReactiveGroup

EXPLICIT_RDKIT_CONTEXT_KEYS = (
    "conjugation_rdkit_execution",
    "conjugation_graph_edit_execution",
    "rdkit_execution",
)


class ExplicitNhsReactiveGroup(BaseModel):
    """Explicit NHS ester atom-index declaration for an RDKit moiety graph."""

    reactive_carbon_index: int = Field(..., ge=0, description="NHS ester acyl carbon index")
    bridging_oxygen_index: int = Field(..., ge=0, description="NHS ester bridging oxygen index")
    leaving_group_atom_indices: tuple[int, ...] = Field(
        ..., min_length=1, description="Moiety atom indices removed as the NHS leaving group"
    )
    retained_atom_indices: tuple[int, ...] | None = Field(
        None, description="Optional moiety atom indices retained after graph surgery"
    )
    candidate_atom_indices: tuple[int, ...] = Field(
        default_factory=tuple, description="Optional NHS search indices used to produce this group"
    )
    evidence: dict[str, Any] = Field(
        default_factory=dict, description="Serializable evidence for the explicit assignment"
    )
    diagnostics: tuple[str, ...] = Field(
        default_factory=tuple, description="Human-readable assignment diagnostics"
    )

    def to_reactive_group(self, moiety_mol: Any) -> NhsReactiveGroup:
        """Convert this explicit declaration into an NHS reactive group model.

        Parameters
        ----------
        moiety_mol : Any
            RDKit molecule-like object used only to infer retained atoms when
            they are not explicitly supplied.

        Returns
        -------
        NhsReactiveGroup
            Reactive group model accepted by the NHS-Lys planning primitive.
        """
        retained = self.retained_atom_indices
        if retained is None:
            atom_count = moiety_mol.GetNumAtoms()
            leaving = set(self.leaving_group_atom_indices)
            retained = tuple(index for index in range(atom_count) if index not in leaving)

        return NhsReactiveGroup(
            reactive_carbon_index=self.reactive_carbon_index,
            bridging_oxygen_index=self.bridging_oxygen_index,
            leaving_group_atom_indices=self.leaving_group_atom_indices,
            retained_atom_indices=retained,
            candidate_atom_indices=self.candidate_atom_indices,
            evidence={"source": "explicit_execution_context", **self.evidence},
            diagnostics=self.diagnostics,
        )


class RdkitGraphEditExecutionRequest(BaseModel):
    """Explicit in-memory request for one RDKit conjugation graph edit.

    RDKit objects are intentionally carried only in memory and are excluded from
    any Pydantic serialization performed on this request.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    protein_mol: Any = Field(..., exclude=True, description="Protein RDKit molecule")
    moiety_mol: Any = Field(..., exclude=True, description="Moiety RDKit molecule")
    protein_topology_atoms: tuple[Any, ...] | None = Field(
        None, exclude=True, description="Optional protein topology atom metadata"
    )
    protein_topology_bonds: tuple[Any, ...] | None = Field(
        None, exclude=True, description="Optional protein topology bond metadata"
    )
    protein_topology_positions: Mapping[int, Sequence[float]] | Sequence[Sequence[float]] | None = (
        Field(None, exclude=True, description="Optional protein atom positions")
    )
    explicit_site_atom_index: int | None = Field(
        None, ge=0, description="Optional explicit protein site atom index"
    )
    explicit_site_hydrogen_indices: tuple[int, ...] | None = Field(
        None, description="Optional explicit site hydrogen atom indices to remove"
    )
    explicit_nhs_group: ExplicitNhsReactiveGroup | NhsReactiveGroup | None = Field(
        None, description="Optional explicit NHS reactive group indices"
    )
    nhs_candidate_atom_indices: tuple[int, ...] | None = Field(
        None, description="Optional atom indices searched during NHS autodetection"
    )
    sanitize: bool = Field(True, description="Attempt RDKit sanitization after graph surgery")

    @field_validator("explicit_nhs_group", mode="before")
    @classmethod
    def normalize_explicit_nhs_group(cls, value: Any) -> Any:
        """Normalize mapping declarations into explicit NHS group models."""
        if value is None or isinstance(value, NhsReactiveGroup | ExplicitNhsReactiveGroup):
            return value
        if isinstance(value, Mapping):
            return ExplicitNhsReactiveGroup.model_validate(value)
        return value


class RdkitGraphEditExecutionSummary(BaseModel):
    """JSON-safe summary for an explicit RDKit conjugation graph edit."""

    attachment: str
    mechanism: str
    site: dict[str, Any]
    moiety: dict[str, Any]
    added_bond: AddedBond
    removed_protein_atom_indices: tuple[int, ...] = Field(default_factory=tuple)
    removed_moiety_atom_indices: tuple[int, ...] = Field(default_factory=tuple)
    removed_atoms_count: int = Field(..., ge=0)
    product_atom_count: int = Field(..., ge=0)
    warnings: tuple[str, ...] = Field(default_factory=tuple)
    topology_unchanged: bool = True


class RdkitGraphEditExecutionResult(BaseModel):
    """In-memory execution result plus a JSON-safe graph-edit summary."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    plan: NhsLysGraphEditPlan
    graph_edit_result: RdkitGraphEditResult = Field(..., exclude=True)
    summary: RdkitGraphEditExecutionSummary


def extract_explicit_rdkit_execution_request(
    context: Mapping[str, Any] | None,
) -> RdkitGraphEditExecutionRequest | None:
    """Extract an explicit RDKit execution request from builder context.

    Parameters
    ----------
    context : Mapping[str, Any] or None
        Optional build context supplied to the covalent modification builder.

    Returns
    -------
    RdkitGraphEditExecutionRequest or None
        Parsed request when explicit RDKit execution keys are present, otherwise
        ``None``.

    Raises
    ------
    ValueError
        If an explicit execution payload is present but cannot be validated.
    """
    if context is None:
        return None

    payload: Any | None = None
    found = False
    for key in EXPLICIT_RDKIT_CONTEXT_KEYS:
        if key in context:
            payload = context[key]
            found = True
            break

    if not found and {"protein_mol", "moiety_mol"} <= set(context):
        payload = context
        found = True

    if not found or payload is None:
        return None
    if isinstance(payload, RdkitGraphEditExecutionRequest):
        return payload

    try:
        return RdkitGraphEditExecutionRequest.model_validate(payload)
    except ValidationError as exc:
        raise ValueError(f"Invalid explicit RDKit execution context: {exc}") from exc
