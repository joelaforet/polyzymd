"""Internal models for covalent modification builders."""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel

from polyzymd.builders.conjugation.diagnostics import ConjugationDiagnosticsReport
from polyzymd.builders.conjugation.metadata import ConjugationMetadata


class ConjugationBuildResult(BaseModel):
    """Result returned by the covalent modification builder."""

    model_config = {"arbitrary_types_allowed": True}

    topology: Any
    metadata: ConjugationMetadata
    diagnostics: ConjugationDiagnosticsReport
