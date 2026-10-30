"""Per-component parameter checks: every term comes from the component's declared sources.

The five checks share one comparison of each component with its reference
(:mod:`polyzymd.audit.references`), made once per audit.
"""

from __future__ import annotations

from typing import Any

from polyzymd.audit.findings import AuditContext, Finding
from polyzymd.audit.modified import build_reference, compare_domains, is_modified
from polyzymd.audit.references import CONSTRAINTS_FROM, DomainComparison, compare_component
from polyzymd.audit.registry import CheckSpec, register

#: The kinds of term each check reads.
_KINDS = {
    "parameters.bonded": ("bonds", "angles", "torsions"),
    "parameters.nonbonded": ("lj",),
    "parameters.one_four": ("exceptions",),
    "parameters.constraints": ("constraints",),
    "charges.source": ("charge",),
    "masses": ("mass",),
}


def _constraints_from(components: Any) -> str:
    """The small-molecule force field whose Constraints section the build lends to the protein."""
    for component in components:
        for domain in component.domains:
            if domain.charges.method == "water_model" and domain.force_field:
                return domain.force_field[0]
    return CONSTRAINTS_FROM


def _comparisons(context: AuditContext) -> dict[str, list[DomainComparison]]:
    if "references" not in context.cache:
        source = _constraints_from(context.components)
        context.cache["references"] = {
            component.name: (
                _modified_comparisons(context, component)
                if is_modified(component)
                else compare_component(context.realized, component, source)
            )
            for component in context.components
        }
    return context.cache["references"]


def modified_reference(context: AuditContext, component: Any) -> Any:
    """The rebuilt reference of a modified protein (built once per audit), or the error."""
    cache = context.cache.setdefault("modified", {})
    if component.name not in cache:
        try:
            cache[component.name] = build_reference(context.realized, component)
        except Exception as error:  # rebuilding runs mBuild, OpenFF and OpenMM
            cache[component.name] = f"could not rebuild it: {type(error).__name__}: {error}"
    return cache[component.name]


def _modified_comparisons(context: AuditContext, component: Any) -> list[DomainComparison]:
    """Compare each domain of a modified protein with the terms its owner predicts."""
    reference = modified_reference(context, component)
    results = []
    for domain in component.domains:
        label = f"{' + '.join(domain.force_field)}, charges {domain.charges.method}"
        if domain.charges.method == "nagl":
            label += f" {domain.charges.model}"
        result = DomainComparison(component.name, domain.name, label, n_molecules=1)
        if isinstance(reference, str):
            result.errors.append(reference)
        results.append(result)
    if isinstance(reference, str):
        return results
    comparisons = compare_domains(context.realized, component, reference)
    for result in results:
        for kind, comparison in comparisons[result.domain].items():
            result.kinds[kind].add(comparison.matched, comparison.total, comparison.mismatches)
    return results


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


for _spec in (
    CheckSpec(
        "parameters.bonded",
        "component",
        "FAIL",
        "bonds, angles, torsions or impropers not from the component's declared force field",
    ),
    CheckSpec(
        "parameters.nonbonded",
        "component",
        "FAIL",
        "Lennard-Jones parameters not from the component's declared force field",
    ),
    CheckSpec(
        "parameters.one_four",
        "component",
        "FAIL",
        "1-2 and 1-3 exclusions or 1-4 pairs not as the owning force field scales them "
        "(GLYCAM06j leaves some glycan 1-4 pairs unscaled; one global scale gets them wrong)",
    ),
    CheckSpec(
        "parameters.constraints",
        "component",
        "FAIL",
        "constraints missing, extra, or at lengths other than the declared force field's",
    ),
    CheckSpec(
        "charges.source",
        "component",
        "FAIL",
        "partial charges that differ from the declared source (library, NAGL model, preset file)",
    ),
    CheckSpec(
        "masses",
        "component",
        "FAIL",
        "atom masses that differ from the force field's (unintended hydrogen mass repartitioning)",
    ),
):
    _register(_spec)


@register(
    CheckSpec(
        "parameters.coverage",
        "system",
        "FAIL",
        "atoms no component owns, molecules other than the components declare, or terms "
        "(bonds, constraints, exceptions) joining two molecules",
    )
)
def coverage(context: AuditContext) -> list[Finding]:
    """Every atom is in one component, every declared molecule is one molecule of the system,
    and no term joins two molecules (a term that does has no declared source)."""
    from collections import Counter

    index = context.realized.index
    n_atoms = context.realized.n_atoms
    problems: list[str] = []
    owned = Counter(i for c in context.components for i in c.atom_indices())
    unowned = [i for i in range(n_atoms) if owned[i] == 0]
    shared = sorted(i for i, n in owned.items() if n > 1)
    if unowned:
        problems.append(f"{len(unowned)} atoms belong to no component (first: atom {unowned[0]})")
    if shared:
        problems.append(f"{len(shared)} atoms belong to two components (first: atom {shared[0]})")
    size = Counter(index.molecule_of)
    for component in context.components:
        wrong = [
            r
            for r in component.molecules
            if len({index.molecule_of[i] for i in r}) != 1
            or size[index.molecule_of[r.start]] != len(r)
        ]
        if wrong:
            problems.append(
                f"{component.name}: {len(wrong)} of its {len(component.molecules)} molecules are "
                f"not single molecules of the system (first at atom {wrong[0].start}): a bond or "
                "constraint joins or splits them"
            )
    spanning = [
        (kind, atoms)
        for terms in index.terms.values()
        for kind, entries in terms.items()
        for atoms, _ in entries
        if len({index.molecule_of[a] for a in atoms}) > 1
    ]
    if spanning:
        problems.append(
            f"{len(spanning)} terms join two molecules (first: {spanning[0][0]} {spanning[0][1]})"
        )
    for component in context.components:
        if not is_modified(component):
            continue
        reference = modified_reference(context, component)
        if isinstance(reference, str):
            continue
        conflicts = reference.conflicts(index)
        if conflicts:
            problems.append(
                f"{component.name}: {len(conflicts)} terms touch two modification sites, so "
                f"no force field owns them (first: {conflicts[0][0]} {conflicts[0][1]})"
            )
    if problems:
        return [Finding("parameters.coverage", "system", "FAIL", "; ".join(problems))]
    return [
        Finding(
            "parameters.coverage",
            "system",
            "PASS",
            f"{n_atoms} atoms in {len(context.components)} components; every term is inside one "
            "declared molecule",
        )
    ]


def comparisons(context: AuditContext) -> Any:
    """The domain comparisons the parameter checks read (for reports and tests)."""
    return _comparisons(context)


@register(
    CheckSpec(
        "parameters.rosemary_alpha",
        "component",
        "WARN",
        "a protein charged whole with NAGL by its force field (OpenFF Rosemary): an alpha "
        "release is not for production, and OpenFF pairs Rosemary with OPC3 water; reported "
        "only for such proteins",
        ("protein",),
    )
)
def rosemary_alpha(context: AuditContext) -> list[Finding]:
    """Warn about a protein that one force field charges whole with NAGL (Rosemary).

    Proteins on library charges (ff14SB) get no finding, so their reports are
    unchanged.
    """
    from polyzymd.config.schema import is_rosemary

    water = None
    if context.config is not None:
        solvent = getattr(context.config, "solvent", None)
        model = getattr(getattr(solvent, "water", None), "model", None)
        water = getattr(model, "value", model)
    findings = []
    for component in context.components:
        if component.type != "protein":
            continue
        protein = next((d for d in component.domains if d.name == "protein"), None)
        if protein is None or protein.charges.method != "nagl":
            continue
        force_fields = list(protein.force_field)
        notes = []
        if is_rosemary(force_fields):
            notes.append(
                f"{' + '.join(force_fields)} is an OpenFF Rosemary alpha release, not for "
                "production"
            )
        if water is not None and str(water).lower() != "opc3":
            notes.append(f"OpenFF pairs Rosemary with OPC3 water; this build uses {water}")
        if notes:
            detail = (
                "; ".join(notes) + ". It parameterizes the whole protein, modifications included."
            )
            findings.append(Finding("parameters.rosemary_alpha", component.name, "WARN", detail))
        else:
            findings.append(
                Finding(
                    "parameters.rosemary_alpha",
                    component.name,
                    "PASS",
                    f"{' + '.join(force_fields)} charges the whole protein with NAGL",
                )
            )
    return findings
