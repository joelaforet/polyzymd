"""Public facade functions for conjugate construction."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

from polyzymd.builders.conjugation.engine import ConjugationEngine
from polyzymd.builders.conjugation.models import ConjugateBuildRequest, ConjugationResult

if TYPE_CHECKING:
    from polyzymd.builders.conjugation.system_workflow import ConjugatedPolymerSystemSettings


def build_conjugate_from_config(
    config: Any | Path | str,
    *,
    output_dir: Path | str | None = None,
    settings: ConjugatedPolymerSystemSettings | None = None,
    free_polymer_seed: int | None = None,
) -> ConjugationResult:
    """Build a conjugate using the existing config-driven workflow.

    ``config`` may be an in-memory ``SimulationConfig`` or a path to a YAML
    config. This facade delegates through :class:`ConjugationEngine` so public
    orchestration stays centralized while preserving the working workflow.
    """
    engine = ConjugationEngine(settings=settings)
    return engine.build_from_config(
        config,
        output_dir=output_dir,
        free_polymer_seed=free_polymer_seed,
    )


def build_conjugate(
    request: ConjugateBuildRequest | Any | Path | str | None = None,
    **kwargs: Any,
) -> ConjugationResult:
    """Build a conjugate through the public engine wrapper.

    Supported inputs are a ``ConjugateBuildRequest``, a config object, a config
    path, or direct request keywords with ``protein_pdb_path`` plus
    ``attachments``. Raw molecule/topology keyword construction fails explicitly
    in the engine instead of silently falling back to a direct OpenFF path.
    """
    settings = kwargs.pop("settings", None)
    engine = ConjugationEngine(settings=settings)
    return engine.build(request, **kwargs)
