"""Per-component parameter checks: every term comes from the component's declared sources.

The five checks share one comparison of each component with its reference
(:mod:`polyzymd.audit.references`), made once per audit.
"""

from __future__ import annotations

from typing import Any

from polyzymd.audit.findings import AuditContext, Finding
from polyzymd.audit.references import DomainComparison, compare_component
from polyzymd.audit.registry import CheckSpec, register

#: The kinds of term each check reads.
_KINDS = {
    "parameters.bonded": ("bonds", "angles", "torsions"),
    "parameters.nonbonded": ("lj", "exceptions"),
    "parameters.constraints": ("constraints",),
    "charges.source": ("charge",),
    "masses": ("mass",),
}


def _comparisons(context: AuditContext) -> dict[str, list[DomainComparison]]:
    if "references" not in context.cache:
        context.cache["references"] = {
            component.name: compare_component(context.realized, component)
            for component in context.components
        }
    return context.cache["references"]


def _findings(check: str, severity: str, context: AuditContext) -> list[Finding]:
    findings = []
    for component in context.components:
        problems, examples, totals = [], [], [0, 0]
        for domain in _comparisons(context)[component.name]:
            if domain.errors:
                problems += [f"{domain.domain}: {e}" for e in domain.errors]
                continue
            if check == "charges.source" and not domain.charges_recomputed:
                problems.append(
                    f"{domain.domain}: {domain.reference} is not recomputed by the audit, so "
                    "these charges are not checked"
                )
                continue
            for kind in _KINDS[check]:
                result = domain.kinds[kind]
                totals[0] += result.matched
                totals[1] += result.total
                if not result.ok:
                    problems.append(
                        f"{domain.domain}: {result.total - result.matched} of {result.total} "
                        f"{kind} differ from {domain.reference}"
                    )
                    examples += result.examples
        if not problems:
            findings.append(
                Finding(
                    check,
                    component.name,
                    "PASS",
                    f"{totals[0]} {'/'.join(_KINDS[check])} terms match the declared sources",
                )
            )
            continue
        uncheckable = all("not recomputed" in p for p in problems)
        findings.append(
            Finding(
                check,
                component.name,
                "WARN" if uncheckable else severity,  # type: ignore[arg-type]
                "; ".join(problems),
                examples[:3],
            )
        )
    return findings


def _register(spec: CheckSpec) -> None:
    @register(spec)
    def run(context: AuditContext, _spec: CheckSpec = spec) -> list[Finding]:
        return _findings(_spec.id, _spec.severity, context)


_ORIGIN = "1.3 force-field order bug: Sage terms on ff14SB proteins (Paper 1)"

for _spec in (
    CheckSpec(
        "parameters.bonded",
        "component",
        "FAIL",
        "bonds, angles, torsions or impropers not from the component's declared force field",
        _ORIGIN,
    ),
    CheckSpec(
        "parameters.nonbonded",
        "component",
        "FAIL",
        "Lennard-Jones or 1-4 terms not from the component's declared force field",
        _ORIGIN,
    ),
    CheckSpec(
        "parameters.constraints",
        "component",
        "FAIL",
        "constraints missing, extra, or at lengths other than the declared force field's",
        _ORIGIN + " (Sage bond-to-hydrogen lengths on the protein)",
    ),
    CheckSpec(
        "charges.source",
        "component",
        "FAIL",
        "partial charges that differ from the declared source (library, NAGL model, preset file)",
        "LipA production runs that recomputed charges",
    ),
    CheckSpec(
        "masses",
        "component",
        "FAIL",
        "atom masses that differ from the force field's (unintended hydrogen mass repartitioning)",
    ),
):
    _register(_spec)


def comparisons(context: AuditContext) -> Any:
    """The domain comparisons the parameter checks read (for reports and tests)."""
    return _comparisons(context)
