"""Running the registered checks on a realized system, and what they find."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal, Optional, Sequence

from polyzymd.audit import registry
from polyzymd.audit.realized import RealizedSystem

Status = Literal["PASS", "FAIL", "WARN", "INFO"]


@dataclass
class Finding:
    """What one check found for one component (or for the whole system).

    Attributes:
        check: The check's registered ID.
        component: The component's name, or ``system``.
        status: ``PASS``; or the check's severity when it found a problem
            (``FAIL`` or ``WARN``); ``INFO`` for a check that only reports.
        detail: One line a person or an agent can act on.
        examples: A few of the offending terms or atoms.
    """

    check: str
    component: str
    status: Status
    detail: str
    examples: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "check": self.check,
            "component": self.component,
            "status": self.status,
            "detail": self.detail,
            "examples": list(self.examples),
        }


@dataclass
class AuditContext:
    """What every check reads: the realized system, its components, and shared work."""

    realized: RealizedSystem
    components: Sequence[Any]
    #: Work several checks share (such as the reference comparisons), by key.
    cache: dict[str, Any] = field(default_factory=dict)


def run_checks(
    realized: RealizedSystem, components: Sequence[Any], ids: Optional[Sequence[str]] = None
) -> list[Finding]:
    """Run every registered check (or those in ``ids``) and return their findings."""
    import polyzymd.audit.checks  # noqa: F401  # registers the checks

    context = AuditContext(realized, components)
    findings: list[Finding] = []
    for spec in registry.catalog():
        if ids is not None and spec.id not in ids:
            continue
        findings += registry.implementation(spec.id)(context)
    return findings
