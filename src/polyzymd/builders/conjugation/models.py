"""Models for covalent modification and conjugation builders."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation.diagnostics import ConjugationDiagnosticsReport
from polyzymd.builders.conjugation.execution import (
    RdkitGraphEditExecutionResult,
    RdkitGraphEditExecutionSummary,
)
from polyzymd.builders.conjugation.metadata import ConjugationMetadata


class ConjugateBuildRequest(BaseModel):
    """Lightweight public request shell for conjugate construction."""

    model_config = {"arbitrary_types_allowed": True}

    config: Any | None = Field(
        default=None,
        description="In-memory SimulationConfig or compatible config object.",
    )
    config_path: Path | None = Field(
        default=None,
        description="Path to a YAML config consumed by the legacy workflow.",
    )
    output_dir: Path | None = Field(
        default=None,
        description="Output directory for generated conjugation artifacts.",
    )
    free_polymer_seed: int | None = Field(
        default=None,
        description="Optional seed forwarded to free-polymer placement.",
    )


class ConjugateBuildResult(BaseModel):
    """Lightweight public result shell for future engine outputs."""

    model_config = {"arbitrary_types_allowed": True}

    output_dir: Path | None = None
    legacy_result: Any | None = Field(
        default=None,
        exclude=True,
        description="Delegated legacy result retained in memory during Phase 1.",
    )


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
