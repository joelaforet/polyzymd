"""Public facade functions for conjugate construction."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from polyzymd.builders.conjugation.engine import ConjugationEngine
from polyzymd.builders.conjugation.system_workflow import (
    ConjugatedPolymerSystemResult,
    ConjugatedPolymerSystemSettings,
)


def build_conjugate_from_config(
    config: Any | Path | str,
    *,
    output_dir: Path | str | None = None,
    settings: ConjugatedPolymerSystemSettings | None = None,
    free_polymer_seed: int | None = None,
) -> ConjugatedPolymerSystemResult:
    """Build a conjugate using the existing config-driven workflow.

    ``config`` may be an in-memory ``SimulationConfig`` or a path to a YAML
    config. Phase 1 intentionally delegates to the legacy workflow entry
    points while establishing a stable public facade.
    """
    engine = ConjugationEngine(settings=settings)
    return engine.build_from_config(
        config,
        output_dir=output_dir,
        free_polymer_seed=free_polymer_seed,
    )


def build_conjugate(*_args: Any, **_kwargs: Any) -> ConjugatedPolymerSystemResult:
    """Build a conjugate from direct construction inputs.

    Direct, non-config conjugate construction is the future engine API. The
    underlying behavior has not been migrated yet, so callers should use
    :func:`build_conjugate_from_config` for Phase 1.
    """
    raise NotImplementedError(
        "Direct build_conjugate(...) inputs are not implemented in the Phase 1 "
        "conjugation facade. Use build_conjugate_from_config(...) to delegate "
        "to the existing config-driven workflow."
    )
