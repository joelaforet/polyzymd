"""Internal models for covalent modification builders."""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation.diagnostics import ConjugationDiagnosticsReport
from polyzymd.builders.conjugation.execution import (
    RdkitGraphEditExecutionResult,
    RdkitGraphEditExecutionSummary,
)
from polyzymd.builders.conjugation.metadata import ConjugationMetadata


class ConjugationBuildResult(BaseModel):
    """Result returned by the covalent modification builder."""

    model_config = {"arbitrary_types_allowed": True}

    topology: Any
    metadata: ConjugationMetadata
    diagnostics: ConjugationDiagnosticsReport
    graph_edit_results: list[RdkitGraphEditExecutionResult] = Field(
        default_factory=list,
        exclude=True,
        description="In-memory RDKit graph edit results excluded from serialization",
    )
    graph_edit_summaries: list[RdkitGraphEditExecutionSummary] = Field(
        default_factory=list,
        description="JSON-safe summaries for explicit graph edit execution",
    )
