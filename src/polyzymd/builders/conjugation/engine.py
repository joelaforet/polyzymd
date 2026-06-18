"""Conjugation orchestration boundary."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from polyzymd.builders.conjugation.models import ConjugateBuildRequest
from polyzymd.builders.conjugation.system_workflow import (
    ConjugatedPolymerSystemResult,
    ConjugatedPolymerSystemSettings,
    build_conjugated_polymer_system_from_config,
    build_conjugated_polymer_system_from_config_path,
)


class ConjugationEngine:
    """Phase-1 boundary for conjugate construction orchestration.

    The engine currently delegates to the existing config-driven conjugated
    polymer workflow. Future phases should move orchestration behind this class
    without duplicating behavior in facade functions.
    """

    def __init__(self, *, settings: ConjugatedPolymerSystemSettings | None = None) -> None:
        self.settings = settings

    def build_from_config(
        self,
        config: Any | Path | str,
        *,
        output_dir: Path | str | None = None,
        settings: ConjugatedPolymerSystemSettings | None = None,
        free_polymer_seed: int | None = None,
    ) -> ConjugatedPolymerSystemResult:
        """Delegate config-based construction to the existing workflow."""
        workflow_settings = self.settings if settings is None else settings
        if isinstance(config, (str, Path)):
            return build_conjugated_polymer_system_from_config_path(
                config,
                output_dir=output_dir,
                settings=workflow_settings,
                free_polymer_seed=free_polymer_seed,
            )

        if output_dir is None:
            raise ValueError(
                "output_dir is required when build_from_config() receives an "
                "in-memory SimulationConfig."
            )
        return build_conjugated_polymer_system_from_config(
            config,
            output_dir=output_dir,
            settings=workflow_settings,
            free_polymer_seed=free_polymer_seed,
        )

    def build_from_request(
        self,
        request: ConjugateBuildRequest,
        *,
        settings: ConjugatedPolymerSystemSettings | None = None,
    ) -> ConjugatedPolymerSystemResult:
        """Build from the lightweight public request shell."""
        source = request.config if request.config is not None else request.config_path
        if source is None:
            raise ValueError("ConjugateBuildRequest requires either config or config_path")
        return self.build_from_config(
            source,
            output_dir=request.output_dir,
            settings=settings,
            free_polymer_seed=request.free_polymer_seed,
        )

    def build(self, *_args: Any, **_kwargs: Any) -> ConjugatedPolymerSystemResult:
        """Build from direct engine inputs once behavior is migrated."""
        raise NotImplementedError(
            "Direct ConjugationEngine.build(...) inputs are not implemented in Phase 1. "
            "Use build_from_config(...) to delegate to the existing workflow."
        )
