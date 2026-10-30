"""Running the registered checks on a realized system, and what they find."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal, Optional, Sequence

from polyzymd.audit import registry
from polyzymd.audit.realized import RealizedSystem

#: ACCEPTED is a FAIL or WARN the config accepts in ``modeling_exceptions``.
Status = Literal["PASS", "FAIL", "WARN", "INFO", "ACCEPTED"]


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
    #: The status a modeling exception accepted (FAIL or WARN), when status is ACCEPTED.
    accepted_from: Optional[str] = None
    #: The reason the config gives for accepting it.
    reason: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        data = {
            "check": self.check,
            "component": self.component,
            "status": self.status,
            "detail": self.detail,
            "examples": list(self.examples),
        }
        if self.accepted_from is not None:
            data["accepted_from"] = self.accepted_from
            data["reason"] = self.reason
        return data


@dataclass
class AuditContext:
    """What every check reads: the realized system, its components, and shared work."""

    realized: RealizedSystem
    components: Sequence[Any]
    #: The SimulationConfig the system was built from; system checks compare with it.
    config: Any = None
    #: Work several checks share (such as the reference comparisons), by key.
    cache: dict[str, Any] = field(default_factory=dict)
    #: What the caller declares for each modified protein, by component name:
    #: ``{"spec": blocks, "parameter_sets": {...}, "base_dir": path}``, read from
    #: the config, independent of the build's own records.
    declared: Optional[dict[str, Any]] = None


def run_checks(
    realized: RealizedSystem,
    components: Sequence[Any],
    ids: Optional[Sequence[str]] = None,
    config: Any = None,
    declared: Optional[dict[str, Any]] = None,
) -> list[Finding]:
    """Run every registered check (or those in ``ids``) and return their findings.

    ``config`` is the SimulationConfig of the build; checks that compare the
    system with the config report INFO "not checked" without it. ``declared``
    holds the modification blocks and parameter sets the config declares for
    each modified protein (see :attr:`AuditContext.declared`).
    """
    import polyzymd.audit.checks  # noqa: F401  # registers the checks

    context = AuditContext(realized, components, config, declared=declared)
    findings: list[Finding] = []
    for spec in registry.catalog():
        if ids is not None and spec.id not in ids:
            continue
        findings += registry.implementation(spec.id)(context)
    return findings
