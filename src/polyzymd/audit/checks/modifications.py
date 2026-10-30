"""Checks of covalently modified proteins: site charges, the backbone, what was applied.

They read the reference the parameter checks rebuild once per audit
(:func:`polyzymd.audit.checks.parameters.modified_reference`) and run on every
protein component; a protein without modifications passes them.
"""

from __future__ import annotations

from typing import Any

from polyzymd.audit.checks.parameters import modified_reference
from polyzymd.audit.findings import AuditContext, Finding
from polyzymd.audit.modified import BACKBONE, is_modified
from polyzymd.audit.registry import CheckSpec, register
from polyzymd.audit.terms import compare

#: Charges agree to this, in e (the parameter checks' tolerance).
_CHARGE_TOLERANCE = 1e-5

_PROTEIN = ("protein",)


def _each_modified(context: AuditContext, check: str, run) -> list[Finding]:
    findings = []
    for component in context.components:
        if component.type not in _PROTEIN:
            continue
        if not is_modified(component):
            findings.append(Finding(check, component.name, "PASS", "no modifications"))
            continue
        reference = modified_reference(context, component)
        if isinstance(reference, str):
            findings.append(Finding(check, component.name, "FAIL", reference))
            continue
        findings.append(run(component, reference))
    return findings


@register(
    CheckSpec(
        "charges.site_integer",
        "component",
        "FAIL",
        "a modification site whose charges, with the rest of its residues, do not sum to their "
        "formal charge",
        _PROTEIN,
    )
)
def site_integer(context: AuditContext) -> list[Finding]:
    def run(component: Any, reference: Any) -> Finding:
        bad = [
            f"{site.name} sums to {site.actual_sum:+.6f} e, not {site.target:+.6f}"
            for site in reference.sites
            if abs(site.actual_sum - site.target) > _CHARGE_TOLERANCE
        ]
        unread = []
        for name, formal in reference.site_formal.items():
            if formal is None:
                unread.append(name)
            elif abs(formal[0] - formal[1]) > _CHARGE_TOLERANCE:
                bad.append(
                    f"{name}: its residues' formal charge is {formal[0]:+.0f} in the built "
                    f"chemistry but {formal[1]:+.0f} from the input residues and the recorded "
                    "changes"
                )
        if unread and not bad:
            return Finding(
                "charges.site_integer",
                component.name,
                "WARN",
                "the input PDB could not be read independently, so the formal charges of "
                f"{', '.join(unread)} are not cross-checked",
            )
        if bad:
            return Finding("charges.site_integer", component.name, "FAIL", "; ".join(bad))
        return Finding(
            "charges.site_integer",
            component.name,
            "PASS",
            f"{len(reference.sites)} sites sum to their targets (formal charge minus the "
            "ff14SB charges of the rest of their residues)",
        )

    return _each_modified(context, "charges.site_integer", run)


@register(
    CheckSpec(
        "charges.site_shift",
        "component",
        "FAIL",
        "a site whose NAGL charges need a uniform shift above the bound to reach its integer "
        "charge (the capped model and the site disagree)",
        _PROTEIN,
    )
)
def site_shift(context: AuditContext) -> list[Finding]:
    def run(component: Any, reference: Any) -> Finding:
        shifted = [site for site in reference.sites if site.shift_per_atom is not None]
        bad = [
            f"{site.name} needs {site.shift_per_atom:+.4f} e per atom (bound {site.bound})"
            for site in shifted
            if abs(site.shift_per_atom) > site.bound
        ]
        if bad:
            return Finding("charges.site_shift", component.name, "FAIL", "; ".join(bad))
        largest = max((abs(site.shift_per_atom) for site in shifted), default=0.0)
        return Finding(
            "charges.site_shift",
            component.name,
            "PASS",
            f"{len(shifted)} NAGL sites; largest shift {largest:.4f} e per atom",
        )

    return _each_modified(context, "charges.site_shift", run)


@register(
    CheckSpec(
        "modifications.applied",
        "component",
        "FAIL",
        "a planned modification missing from the built system: a bond it forms, an atom it "
        "adds or removes, or an atom renamed",
        _PROTEIN,
    )
)
def applied(context: AuditContext) -> list[Finding]:
    def run(component: Any, reference: Any) -> Finding:
        index = context.realized.index
        atoms = set(reference.domain_of)
        built = {
            (min(a, b), max(a, b)) for a in atoms for b in index.bonded_to.get(a, ()) if b in atoms
        }
        missing = sorted(reference.bonds - built)
        extra = sorted(built - reference.bonds)
        names = _built_names(context.realized, sorted(atoms))
        # Attached residues are skipped: readers rename atoms of residues they
        # know (OpenMM's PDB reader gives ACE and NME their standard names).
        renamed = [
            (i, expected)
            for i, expected in zip(sorted(atoms), reference.names)
            if names is not None
            and i not in reference.attached_atoms
            and names[i][0] != _as_read(names[i][1], expected.rsplit(":", 1)[1])
        ]
        problems = []
        if missing:
            problems.append(
                f"{len(missing)} planned bonds are missing (first: "
                f"{_pair(reference, missing[0])})"
            )
        if extra:
            problems.append(
                f"{len(extra)} bonds are not planned (first: {_pair(reference, extra[0])})"
            )
        if renamed:
            problems.append(
                f"{len(renamed)} atoms are not the planned atoms (first: atom {renamed[0][0]} "
                f"should be {renamed[0][1]})"
            )
        if problems:
            return Finding("modifications.applied", component.name, "FAIL", "; ".join(problems))
        return Finding(
            "modifications.applied",
            component.name,
            "PASS",
            f"{len(reference.bonds)} bonds and {len(atoms)} atoms are those of the input "
            "structure with the modifications applied",
        )

    return _each_modified(context, "modifications.applied", run)


@register(
    CheckSpec(
        "parameters.modified_backbone",
        "component",
        "FAIL",
        "the backbone (N, H, CA, HA, C, O) of a modified residue not on ff14SB: its charges "
        "and terms must match its unmodified neighbours'",
        _PROTEIN,
    )
)
def modified_backbone(context: AuditContext) -> list[Finding]:
    def run(component: Any, reference: Any) -> Finding:
        index = context.realized.index
        backbone = set(reference.backbone)
        problems, examples = [], []
        wrong = [
            i
            for i, charge in reference.backbone.items()
            if abs(index.particle[i][0] - charge) > _CHARGE_TOLERANCE
        ]
        if wrong:
            problems.append(f"{len(wrong)} backbone charges are not ff14SB's")
            examples += [
                f"{reference.names[i - min(reference.domain_of)]}: "
                f"{index.particle[i][0]:+.4f} vs {reference.backbone[i]:+.4f}"
                for i in wrong[:3]
            ]
        atoms = list(reference.domain_of)
        actual = index.extract(
            atoms,
            owned=lambda term: bool(backbone & set(term))
            and all(reference.domain_of.get(i) == "protein" for i in term),
            terms_of=atoms,
        )
        expected = reference.expected["protein"]
        for kind in ("bonds", "angles", "torsions"):
            touching = {k: v for k, v in expected[kind].items() if backbone & set(k)}
            comparison = compare(kind, actual[kind], touching)
            if not comparison.ok:
                problems.append(
                    f"{comparison.total - comparison.matched} of {comparison.total} {kind} "
                    "of modified-residue backbones are not ff14SB's"
                )
                examples += [repr(m) for m in comparison.mismatches]
        offset = min(reference.domain_of)
        joined = [
            i
            for i, domain in reference.domain_of.items()
            if domain != "protein"
            and reference.names[i - offset].rsplit(":", 1)[1] in BACKBONE
            and i not in reference.designed_backbone
            and i not in reference.attached_atoms
        ]
        if joined:
            problems.append(
                f"{len(joined)} backbone atoms are in a site although no change touched them "
                f"(first: {reference.names[joined[0] - offset]})"
            )
        if problems:
            return Finding(
                "parameters.modified_backbone",
                component.name,
                "FAIL",
                "; ".join(problems),
                examples[:3],
            )
        return Finding(
            "parameters.modified_backbone",
            component.name,
            "PASS",
            f"{len(backbone)} backbone atoms of modified residues carry ff14SB charges and terms",
        )

    return _each_modified(context, "parameters.modified_backbone", run)


@register(
    CheckSpec(
        "parameters.glycan_force_field",
        "component",
        "WARN",
        "a glycan not parameterized with GLYCAM06j, the recommended glycan force field "
        "(general force fields give stiffer glycosidic angles and other torsion minima)",
        _PROTEIN,
    )
)
def glycan_force_field(context: AuditContext) -> list[Finding]:
    def run(component: Any, reference: Any) -> Finding:
        others = [
            f"{name} on {' + '.join(reference.site_force_fields[name])}"
            for name, site_type in reference.site_types.items()
            if site_type == "glycan"
            and not any("glycam" in f.lower() for f in reference.site_force_fields[name])
        ]
        if others:
            return Finding(
                "parameters.glycan_force_field",
                component.name,
                "WARN",
                "glycans not on GLYCAM06j: " + ", ".join(others) + ". GLYCAM06j is "
                "recommended for glycans; see the glycans how-to.",
            )
        n = sum(t == "glycan" for t in reference.site_types.values())
        return Finding(
            "parameters.glycan_force_field",
            component.name,
            "PASS",
            f"{n} glycans on GLYCAM06j" if n else "no glycans",
        )

    return _each_modified(context, "parameters.glycan_force_field", run)


@register(
    CheckSpec(
        "modifications.protonation_state",
        "component",
        "WARN",
        "a mutation that leaves a side chain in a state other than the usual one at pH 7 "
        "(a CCD component's neutral Asp or Glu, or doubly protonated His)",
        _PROTEIN,
    )
)
def protonation_state(context: AuditContext) -> list[Finding]:
    def run(component: Any, reference: Any) -> Finding:
        if reference.protonation_warnings:
            return Finding(
                "modifications.protonation_state",
                component.name,
                "WARN",
                "; ".join(reference.protonation_warnings),
            )
        return Finding(
            "modifications.protonation_state",
            component.name,
            "PASS",
            "every mutated ionizable side chain is in its usual state at pH 7",
        )

    return _each_modified(context, "modifications.protonation_state", run)


def _built_names(realized: Any, atoms: list[int]) -> dict[int, tuple[str, str]] | None:
    """The (atom name, residue name) of each of ``atoms`` in the realized topology."""
    if realized.topology is None:
        return None
    names = {atom.index: (atom.name, atom.residue.name) for atom in realized.topology.atoms()}
    return {i: names[i] for i in atoms}


def _as_read(residue_name: str, atom_name: str) -> str:
    """``atom_name`` as OpenMM's PDB reader names it in a ``residue_name`` residue.

    The reader renames atoms of residues it knows (an NME cap's CCD ``HN2``
    becomes ``H2``), so a planned name is compared after the same renaming.
    """
    from openmm.app import PDBFile

    PDBFile._loadNameReplacementTables()
    return PDBFile._atomNameReplacements.get(residue_name, {}).get(atom_name, atom_name)


def _pair(reference: Any, pair: tuple[int, int]) -> str:
    offset = min(reference.domain_of)
    return f"{reference.names[pair[0] - offset]}-{reference.names[pair[1] - offset]}"


@register(
    CheckSpec(
        "modifications.sites",
        "component",
        "FAIL",
        "a modification site whose atoms are not those its changes define (side chain from CB "
        "of the anchor, touched backbone atoms, added atoms, attached residues), or that "
        "differ from the component record's",
        _PROTEIN,
    )
)
def sites(context: AuditContext) -> list[Finding]:
    def run(component: Any, reference: Any) -> Finding:
        offset = min(reference.domain_of)
        in_sites = {i for i, d in reference.domain_of.items() if d != "protein"}
        problems = []
        outside = sorted(in_sites - reference.rule_atoms)
        missing = sorted(reference.rule_atoms - in_sites)
        if outside:
            problems.append(
                f"{len(outside)} site atoms are not atoms any change defines (first: "
                f"{reference.names[outside[0] - offset]})"
            )
        if missing:
            problems.append(
                f"{len(missing)} atoms the changes define are in no site (first: "
                f"{reference.names[missing[0] - offset]})"
            )
        if reference.neutral_termini:
            problems.append(
                "a neutral terminus is not supported in PolyzyMD 1.4 (termini are simulated "
                "charged at physiological pH; cap a truncated chain with ACE or NME instead): "
                + ", ".join(reference.neutral_termini)
            )
        for domain in component.domains:
            recorded = {i for r in domain.atoms for i in r}
            derived = {i for i, d in reference.domain_of.items() if d == domain.name}
            if recorded != derived:
                problems.append(
                    f"the component record's {domain.name} domain holds {len(recorded)} atoms, "
                    f"the rebuilt protein's {len(derived)} ({len(recorded ^ derived)} differ)"
                )
        if problems:
            return Finding("modifications.sites", component.name, "FAIL", "; ".join(problems))
        return Finding(
            "modifications.sites",
            component.name,
            "PASS",
            f"{len(in_sites)} site atoms in {len(component.domains) - 1} sites, as the changes "
            "define them",
        )

    return _each_modified(context, "modifications.sites", run)


@register(
    CheckSpec(
        "modifications.declared",
        "component",
        "FAIL",
        "a built protein whose modifications or parameter sets differ from the config's "
        "mutations:, modifications: and parameter sets",
        _PROTEIN,
    )
)
def declared(context: AuditContext) -> list[Finding]:
    from polyzymd.audit.modified import declared_differences

    def run(component: Any, reference: Any) -> Finding:
        given = (context.declared or {}).get(component.name)
        if given is None:
            return Finding(
                "modifications.declared",
                component.name,
                "INFO",
                "no declared modification blocks were given, so the build is not compared "
                "with them",
            )
        problems = declared_differences(reference, component, given)
        if problems:
            return Finding(
                "modifications.declared", component.name, "FAIL", "; ".join(problems[:5])
            )
        return Finding(
            "modifications.declared",
            component.name,
            "PASS",
            f"{len(reference.records)} changes and their parameter sets are the declared ones",
        )

    return _each_modified(context, "modifications.declared", run)


@register(
    CheckSpec(
        "modifications.chemistry",
        "component",
        "FAIL",
        "a modified protein whose chemistry (elements, bond orders, formal charges) differs "
        "from an independent reading of its structure with openff-pablo",
        _PROTEIN,
    )
)
def chemistry(context: AuditContext) -> list[Finding]:
    def run(component: Any, reference: Any) -> Finding:
        if reference.chemistry:
            return Finding(
                "modifications.chemistry", component.name, "FAIL", "; ".join(reference.chemistry)
            )
        return Finding(
            "modifications.chemistry",
            component.name,
            "PASS",
            f"{len(reference.domain_of)} atoms, their bonds and formal charges match an "
            "independent reading of the structure",
        )

    return _each_modified(context, "modifications.chemistry", run)
