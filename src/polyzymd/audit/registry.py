"""The catalog of physics checks: one entry per check, the single source of truth.

Every check a build runs is registered here with a stable ID, the severity it has
unless a config's ``modeling_exceptions`` accepts it, and a sentence saying which
modelling mistake it catches. The docs page listing every check and the comments
in the config template are generated from this catalog, so they cannot drift from
the code.
"""

from __future__ import annotations

import difflib
from dataclasses import dataclass
from typing import Callable, Literal, Optional

Severity = Literal["FAIL", "WARN", "INFO"]
Scope = Literal["component", "system"]


@dataclass(frozen=True)
class CheckSpec:
    """One check of a built system.

    Attributes:
        id: Stable identifier, used by ``modeling_exceptions`` and in the provenance
            file (``parameters.bonded``, ``geometry.buried_water``, ...).
        scope: ``component`` if the check runs once per component (and an exception
            may name a component), ``system`` if it runs once per build.
        severity: ``FAIL`` stops the build unless the config accepts it; ``WARN``
            and ``INFO`` are recorded only.
        catches: The modelling mistake the check catches, in one sentence.
        component_types: Component types the check applies to; None means all.
    """

    id: str
    scope: Scope
    severity: Severity
    catches: str
    component_types: Optional[tuple[str, ...]] = None

    def applies_to(self, component_type: str) -> bool:
        return self.component_types is None or component_type in self.component_types


_REGISTRY: dict[str, CheckSpec] = {}
_IMPLEMENTATIONS: dict[str, Callable] = {}


def register(spec: CheckSpec) -> Callable[[Callable], Callable]:
    """Register ``spec`` and the function decorated with it as its implementation.

    Raises:
        ValueError: If another check already has the same ID.
    """
    if spec.id in _REGISTRY:
        raise ValueError(f"check {spec.id!r} is registered twice")

    def decorator(function: Callable) -> Callable:
        _REGISTRY[spec.id] = spec
        _IMPLEMENTATIONS[spec.id] = function
        return function

    return decorator


def get(check_id: str) -> CheckSpec:
    """Return the check with ``check_id``.

    Raises:
        KeyError: If there is none, naming the closest registered IDs.
    """
    try:
        return _REGISTRY[check_id]
    except KeyError:
        close = difflib.get_close_matches(check_id, list(_REGISTRY), n=3)
        hint = f"; did you mean {', '.join(close)}?" if close else ""
        raise KeyError(f"no check named {check_id!r}{hint}") from None


def implementation(check_id: str) -> Callable:
    """Return the function that runs the check with ``check_id``."""
    get(check_id)
    return _IMPLEMENTATIONS[check_id]


def catalog() -> list[CheckSpec]:
    """Return every registered check, ordered by ID."""
    return [_REGISTRY[key] for key in sorted(_REGISTRY)]
