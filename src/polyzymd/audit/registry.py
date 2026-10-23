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
        origin: The real failure that motivated the check, if any.
        component_types: Component types the check applies to; None means all.
    """

    id: str
    scope: Scope
    severity: Severity
    catches: str
    origin: str = ""
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


def catalog_markdown() -> str:
    """The docs page listing every check, generated from the registry.

    ``docs/source/reference/physics_checks.md`` is this text; a test keeps it
    current, so the page cannot drift from the code.
    """
    import polyzymd.audit.checks  # noqa: F401  # registers the checks

    lines = [
        "# Physics checks",
        "",
        "Every PolyzyMD build runs each of these checks on the system the requested engine",
        "will run (OpenMM `system.xml`, or GROMACS `.top`/`.gro` with the production `.mdp`),",
        "and records the result in `build_provenance.json` and `build_provenance.md`.",
        "A check whose default is FAIL stops the build; WARN and INFO are recorded only.",
        "",
        "To keep a deliberate departure, accept the failing (or warning) check in the config with a",
        "reason. The build record then lists the accepted finding, its reason and the",
        "atoms it covered:",
        "",
        "```yaml",
        "modeling_exceptions:",
        "  - check: geometry.buried_water",
        "    component: protein        # optional; every component the check runs on if left out",
        '    reason: "Crystal waters W301 and W305 kept deliberately (PDB 4K6G)"',
        "```",
        "",
        "An exception that matches no failing check is reported, so stale entries do not",
        "linger in a published config.",
        "",
        "This page is generated from PolyzyMD's check registry",
        "(`polyzymd.audit.registry.catalog_markdown`).",
        "",
        "| Check | Runs on | Default | Catches | Origin |",
        "|---|---|---|---|---|",
    ]
    for spec in catalog():
        scope = "each component" if spec.scope == "component" else "the system"
        if spec.component_types:
            scope += f" ({', '.join(spec.component_types)})"
        catches = spec.catches.replace("|", "\\|")
        origin = spec.origin.replace("|", "\\|")
        lines.append(f"| `{spec.id}` | {scope} | {spec.severity} | {catches} | {origin} |")
    return "\n".join(lines) + "\n"
