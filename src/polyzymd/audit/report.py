"""The audit of one build: every check's finding, with the config's accepted exceptions applied.

A failing check is accepted only by a ``modeling_exceptions`` entry that names
it (and, if the entry names one, its component). An accepted finding keeps its
detail and gains the reason and the atoms it covered, so the record shows both
what PolyzyMD would have refused and why the user chose it. An exception that
matches no failing finding is reported, so stale exceptions do not linger in a
published config.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional, Sequence

from polyzymd.audit.findings import Finding, run_checks
from polyzymd.audit.realized import RealizedSystem

#: Statuses that stop a build unless accepted.
BLOCKING = ("FAIL",)


@dataclass
class AcceptedException:
    """A failing or warning finding the config accepts, and what it covered."""

    check: str
    component: str
    status: str
    reason: str
    detail: str
    atoms: list[list[int]] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "check": self.check,
            "component": self.component,
            "status": self.status,
            "reason": self.reason,
            "detail": self.detail,
            "atoms": self.atoms,
        }


@dataclass
class AuditReport:
    """Every finding of one build, with accepted exceptions applied."""

    findings: list[Finding]
    accepted: list[AcceptedException] = field(default_factory=list)

    @property
    def failures(self) -> list[Finding]:
        return [f for f in self.findings if f.status in BLOCKING]

    @property
    def warnings(self) -> list[Finding]:
        return [f for f in self.findings if f.status == "WARN"]

    @property
    def passed(self) -> bool:
        return not self.failures

    def counts(self) -> dict[str, int]:
        statuses = ("PASS", "FAIL", "WARN", "INFO", "ACCEPTED")
        return {s: sum(f.status == s for f in self.findings) for s in statuses}

    def to_dict(self) -> dict[str, Any]:
        return {
            "passed": self.passed,
            "counts": self.counts(),
            "findings": [f.to_dict() for f in self.findings],
            "accepted_exceptions": [a.to_dict() for a in self.accepted],
        }

    def to_markdown(self) -> str:
        verdict = "PASSED" if self.passed else "FAILED"
        counts = self.counts()
        lines = [
            f"## Physics audit: {verdict}",
            "",
            ", ".join(f"{n} {s.lower()}" for s, n in counts.items() if n),
            "",
            "| Status | Check | Component | Detail |",
            "|---|---|---|---|",
        ]
        order = {"FAIL": 0, "ACCEPTED": 1, "WARN": 2, "PASS": 3, "INFO": 4}
        for f in sorted(self.findings, key=lambda f: (order[f.status], f.check, f.component)):
            detail = f.detail.replace("|", "\\|")
            lines.append(f"| {f.status} | `{f.check}` | {f.component} | {detail} |")
        if self.accepted:
            lines += ["", "### Accepted exceptions", ""]
            for a in self.accepted:
                lines.append(f"- `{a.check}` on {a.component} ({a.status} accepted): {a.reason}")
        return "\n".join(lines) + "\n"


def audit(
    realized: RealizedSystem,
    components: Sequence[Any],
    config: Any = None,
    ids: Optional[Sequence[str]] = None,
) -> AuditReport:
    """Run the checks on ``realized`` and apply the config's ``modeling_exceptions``."""
    findings = run_checks(realized, components, ids, config)
    exceptions = list(getattr(config, "modeling_exceptions", []) or [])
    atoms_of = {c.name: _ranges(c) for c in components}
    used = [False] * len(exceptions)
    accepted: list[AcceptedException] = []
    for finding in findings:
        if finding.status not in ("FAIL", "WARN"):
            continue
        for k, exception in enumerate(exceptions):
            if exception.check != finding.check:
                continue
            if exception.component not in (None, finding.component):
                continue
            used[k] = True
            accepted.append(
                AcceptedException(
                    finding.check,
                    finding.component,
                    finding.status,
                    exception.reason,
                    finding.detail,
                    atoms_of.get(finding.component, [[0, realized.n_atoms]]),
                )
            )
            finding.status = "ACCEPTED"
            finding.detail = f"{finding.detail} [accepted: {exception.reason}]"
            break
    for k, exception in enumerate(exceptions):
        if not used[k]:
            where = f" on {exception.component}" if exception.component else ""
            findings.append(
                Finding(
                    exception.check,
                    exception.component or "system",
                    "WARN",
                    f"modeling_exceptions accepts {exception.check}{where}, but that check found "
                    "nothing to accept; remove the entry if it no longer applies",
                )
            )
    return AuditReport(findings, accepted)


def _ranges(component: Any) -> list[list[int]]:
    from polyzymd.core.components import merge_ranges

    return [[r.start, r.stop] for r in merge_ranges(component.molecules)]
