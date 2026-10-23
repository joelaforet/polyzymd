"""Restraint checks: the system restrains the atoms, form and strength the config names.

Atom identity is checked with a different library from the one that chose the
atoms: the builder resolves selections with MDTraj on its topology; the audit
selects again with MDAnalysis on the built PDB, inside the component's segment
ID, and requires the same single atom.
"""

from __future__ import annotations

import re
from typing import Any, Optional

from polyzymd.audit.findings import AuditContext, Finding
from polyzymd.audit.registry import CheckSpec, register

#: Name prefix PolyzyMD gives its restraint forces (polyzymd.core.restraints).
RESTRAINT_FORCE_PREFIX = "polyzymd_restraint:"

#: The energy expression of each restraint type, as PolyzyMD writes it.
EXPRESSIONS = {
    "flat_bottom": "step(r - r0) * 0.5 * k * (r - r0)^2",
    "upper_wall": "step(r - r0) * 0.5 * k * (r - r0)^2",
    "harmonic": "0.5 * k * (r - r0)^2",
    "lower_wall": "step(r0 - r) * 0.5 * k * (r0 - r)^2",
}


def _restraint_forces(context: AuditContext) -> dict[str, tuple[Any, list]]:
    """{restraint name: (force, [(atom1, atom2, {param: value})])} from the realized system."""
    system = context.realized.system
    found = {}
    for number in context.realized.index.restraint_forces:
        force = system.getForce(number)
        name = force.getName()
        if not name.startswith(RESTRAINT_FORCE_PREFIX):
            continue
        params = [force.getPerBondParameterName(k) for k in range(force.getNumPerBondParameters())]
        bonds = []
        for b in range(force.getNumBonds()):
            a1, a2, values = force.getBondParameters(b)
            bonds.append((a1, a2, dict(zip(params, values))))
        found[name[len(RESTRAINT_FORCE_PREFIX) :]] = (force, bonds)
    return found


def _enabled(context: AuditContext) -> list[Any]:
    return [r for r in context.config.restraints if r.enabled]


@register(
    CheckSpec(
        "restraint.parameters",
        "system",
        "FAIL",
        "a restraint whose form, distance or force constant differs from the config, or a "
        "restraint missing or extra",
    )
)
def parameters(context: AuditContext) -> list[Finding]:
    if context.config is None:
        return [Finding("restraint.parameters", "system", "INFO", "not checked: no config")]
    wanted = {r.name: r for r in _enabled(context)}
    found = _restraint_forces(context)
    unnamed = [
        n
        for n in context.realized.index.restraint_forces
        if not context.realized.system.getForce(n).getName().startswith(RESTRAINT_FORCE_PREFIX)
    ]
    if not wanted and not found:
        return []
    problems = []
    for name in sorted(set(wanted) - set(found)):
        problems.append(f"{name}: in the config but not in the system")
    for name in sorted(set(found) - set(wanted)):
        problems.append(f"{name}: in the system but not an enabled restraint of the config")
    for name in sorted(set(wanted) & set(found)):
        restraint, (force, bonds) = wanted[name], found[name]
        expression = EXPRESSIONS[restraint.type.value]
        if force.getEnergyFunction() != expression:
            problems.append(
                f"{name}: energy {force.getEnergyFunction()!r}, expected {expression!r}"
            )
        if len(bonds) != 1:
            problems.append(f"{name}: {len(bonds)} atom pairs, expected 1")
            continue
        values = bonds[0][2]
        r0, k = values.get("r0"), values.get("k")
        if r0 is None or abs(r0 - restraint.distance / 10.0) > 1e-9:
            problems.append(f"{name}: r0 {r0} nm, config {restraint.distance} Å")
        if k is None or abs(k - restraint.force_constant) > 1e-6 * max(
            1.0, restraint.force_constant
        ):
            problems.append(f"{name}: k {k}, config {restraint.force_constant} kJ/mol/nm²")
    if unnamed:
        problems.append(
            f"{len(unnamed)} restraint-type forces carry no PolyzyMD restraint name, so nothing "
            "says what they restrain"
        )
    if problems:
        return [Finding("restraint.parameters", "system", "FAIL", "; ".join(problems))]
    return [
        Finding(
            "restraint.parameters",
            "system",
            "PASS",
            f"{len(found)} restraints with the config's form, distance and force constant",
        )
    ]


def _mdanalysis_selection(selection: str) -> str:
    """Translate PolyzyMD's selection words to MDAnalysis's (pdbindex is the PDB serial)."""
    text = re.sub(r"\bpdbindex\s+(\d+)", r"id \1", selection)
    text = re.sub(r"\bchain(?:id)?\s+(\S+)", r"chainID \1", text)
    return text


def _pdb(context: AuditContext) -> Optional[str]:
    return next((path for path in context.realized.files if path.endswith(".pdb")), None)


@register(
    CheckSpec(
        "restraint.atom_identity",
        "system",
        "FAIL",
        "a restrained atom that is not the atom its selection names in its component",
        "pdbindex off-by-one in Paper 1 restraint selections",
    )
)
def atom_identity(context: AuditContext) -> list[Finding]:
    if context.config is None or not _enabled(context):
        return []
    pdb = _pdb(context)
    if pdb is None:
        return [
            Finding("restraint.atom_identity", "system", "WARN", "no built PDB to select atoms in")
        ]
    import MDAnalysis as mda

    universe = mda.Universe(pdb)
    segid = {c.name: c.segid for c in context.components}
    found = _restraint_forces(context)
    problems, verified = [], []
    for restraint in _enabled(context):
        if restraint.name not in found or len(found[restraint.name][1]) != 1:
            continue  # restraint.parameters reports it
        atoms = found[restraint.name][1][0][:2]
        for spec, atom in zip((restraint.atom1, restraint.atom2), atoms):
            query = f"segid {segid.get(spec.component, '?')} and ({_mdanalysis_selection(spec.selection)})"
            try:
                chosen = universe.select_atoms(query)
            except Exception as error:  # MDAnalysis raises several types for a bad selection
                problems.append(f"{restraint.name}: MDAnalysis cannot read {query!r}: {error}")
                continue
            if len(chosen) != 1 or int(chosen.indices[0]) != atom:
                names = ", ".join(
                    f"{a.resname}{a.resid}:{a.name} (atom {a.index})" for a in chosen[:3]
                )
                problems.append(
                    f"{restraint.name}: the system restrains atom {atom}, but {query!r} selects "
                    f"{len(chosen)} atom(s){': ' + names if names else ''}"
                )
            else:
                a = chosen[0]
                verified.append(f"{spec.component} {a.resname}{a.resid}:{a.name}")
    if problems:
        return [Finding("restraint.atom_identity", "system", "FAIL", "; ".join(problems))]
    return [
        Finding(
            "restraint.atom_identity",
            "system",
            "PASS",
            "restrained atoms re-selected with MDAnalysis: " + ", ".join(verified),
        )
    ]
