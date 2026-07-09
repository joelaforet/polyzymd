"""Shared specs-first construction settings and result models."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation._linkage import CrosslinkValidationResult, ResolvedAttachmentPlan
from polyzymd.builders.conjugation.pablo.ingestion import PabloIngestionResult
from polyzymd.builders.conjugation.pablo.parameterization import (
    InterchangeParameterizationResult,
    InterchangeParameterizationSettings,
)
from polyzymd.builders.conjugation.placement import (
    PackmolModifierPlacementResult,
    PackmolModifierPlacementSettings,
)
from polyzymd.builders.conjugation.relaxation import (
    ConjugateRelaxationResult,
    ConjugateRelaxationSettings,
)
from polyzymd.builders.conjugation.structure.pdb import CrosslinkedPdbAssemblyResult


class ModifierConstructionSettings(BaseModel):
    """Settings for specs-first modifier-linking construction."""

    crosslinked_pdb_name: str = "assembled_crosslinked.pdb"
    placement: PackmolModifierPlacementSettings = Field(
        default_factory=PackmolModifierPlacementSettings
    )
    parameterization: InterchangeParameterizationSettings = Field(
        default_factory=InterchangeParameterizationSettings
    )
    relaxation: ConjugateRelaxationSettings = Field(default_factory=ConjugateRelaxationSettings)
    run_relaxation: bool = True


class ModifierConstructionResult(BaseModel):
    """Specs-first construction result and artifact summary."""

    model_config = {"arbitrary_types_allowed": True}

    output_dir: Path
    resolved_plan: ResolvedAttachmentPlan
    crosslink_validation: CrosslinkValidationResult
    placement: PackmolModifierPlacementResult
    assembly: CrosslinkedPdbAssemblyResult
    pablo: PabloIngestionResult
    parameterization: InterchangeParameterizationResult
    relaxation: ConjugateRelaxationResult | None = None
    product_state_pablo_library: Any | None = Field(default=None, exclude=True)
    crosslinked_pdb_path: Path
    validation_report_path: Path | None = None
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)
