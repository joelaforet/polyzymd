"""Public models for conjugate construction."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field


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


class ConjugationResult(BaseModel):
    """Public conjugation orchestration result.

    The current engine delegates to the config-driven system workflow and keeps
    the legacy result in memory while exposing a stable, lightweight summary for
    facade callers.
    """

    model_config = {"arbitrary_types_allowed": True}

    status: str = "completed"
    output_dir: Path | None = None
    config_path: Path | None = None
    crosslinked_conjugate_pdb_path: Path | None = None
    minimized_conjugate_pdb_path: Path | None = None
    equilibrated_conjugate_pdb_path: Path | None = None
    relaxed_conjugate_pdb_path: Path | None = None
    solvated_pdb_path: Path | None = None
    workflow_json_path: Path | None = None
    final_interchange_created: bool | None = None
    artifact_paths: dict[str, Path] = Field(default_factory=dict)
    legacy_result: Any | None = Field(
        default=None,
        exclude=True,
        description="Delegated legacy result retained in memory during the refactor.",
    )

    @classmethod
    def from_legacy_result(
        cls,
        legacy_result: Any,
        *,
        config_path: Path | str | None = None,
    ) -> "ConjugationResult":
        """Create a public result summary from the config workflow result."""
        construction = getattr(legacy_result, "construction", None)
        smoke = getattr(construction, "smoke", None)
        local_minimization = getattr(construction, "local_minimization", None)
        crosslinked_path = _optional_path(getattr(construction, "crosslinked_pdb_path", None))
        minimized_path = _optional_path(
            getattr(smoke, "minimized_pdb_path", None)
            or getattr(local_minimization, "relaxed_pdb_path", None)
        )
        equilibrated_path = _optional_path(getattr(smoke, "equilibrated_pdb_path", None))
        relaxed_path = _optional_path(getattr(legacy_result, "relaxed_conjugate_pdb_path", None))
        solvated_path = _optional_path(getattr(legacy_result, "solvated_pdb_path", None))
        workflow_path = _optional_path(getattr(legacy_result, "workflow_json_path", None))
        artifact_paths = {
            name: path
            for name, path in {
                "crosslinked_conjugate_pdb": crosslinked_path,
                "minimized_conjugate_pdb": minimized_path,
                "equilibrated_conjugate_pdb": equilibrated_path,
                "relaxed_conjugate_pdb": relaxed_path,
                "solvated_pdb": solvated_path,
                "workflow_json": workflow_path,
            }.items()
            if path is not None
        }
        return cls(
            status="completed",
            output_dir=_optional_path(getattr(legacy_result, "output_dir", None)),
            config_path=_optional_path(config_path),
            crosslinked_conjugate_pdb_path=crosslinked_path,
            minimized_conjugate_pdb_path=minimized_path,
            equilibrated_conjugate_pdb_path=equilibrated_path,
            relaxed_conjugate_pdb_path=relaxed_path,
            solvated_pdb_path=solvated_path,
            workflow_json_path=workflow_path,
            final_interchange_created=getattr(legacy_result, "final_interchange_created", None),
            artifact_paths=artifact_paths,
            legacy_result=legacy_result,
        )


ConjugateBuildResult = ConjugationResult


def _optional_path(value: Any) -> Path | None:
    if value is None:
        return None
    return Path(value)
