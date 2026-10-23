"""System checks: the build holds what the config asked for, as a whole.

These compare the realized system with the config itself (not with the
builder's records of it): the components and their declared sources, counts,
net charge, salt or ionic strength, the constraint policy and the nonbonded
settings.
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Optional

from polyzymd.audit.findings import AuditContext, Finding
from polyzymd.audit.registry import CheckSpec, register

#: Molar concentration of pure water, as the builder sizes salt and ionic strength.
WATER_MOLARITY = 55.5
#: Mass of one TIP3P water (Da), the unit of the builder's water-equivalents.
WATER_MASS = 15.999 + 2 * 1.008


def _no_config(check: str) -> list[Finding]:
    return [Finding(check, "system", "INFO", "not checked: the audit was given no config")]


def _of_type(context: AuditContext, type_: str) -> list[Any]:
    return [c for c in context.components if c.type == type_]


def _domain(component: Any, name: str) -> Optional[Any]:
    return next((d for d in component.domains if d.name == name), None)


def _charges(spec: Any) -> tuple[str, Optional[str]]:
    return spec.charges.method, getattr(spec.charges, "model", None)


@register(
    CheckSpec(
        "components.declared",
        "system",
        "FAIL",
        "a component, its segment ID, force field or charge source that differs from the config",
        "parameter checks compare with the builder's records; this ties the records to the config",
    )
)
def declared(context: AuditContext) -> list[Finding]:
    """Every solute and the solvent appear as components with the config's sources."""
    config = context.config
    if config is None:
        return _no_config("components.declared")
    problems: list[str] = []
    records = {c.name: c for c in context.components}
    for solute in config.solutes:
        record = records.get(solute.name)
        if record is None:
            problems.append(f"solute {solute.name!r} has no component in the build")
            continue
        if (record.type, record.segid) != (solute.type, solute.segid):
            problems.append(
                f"{solute.name}: built as {record.type}/{record.segid}, config says "
                f"{solute.type}/{solute.segid}"
            )
        for domain in record.domains:
            if list(domain.force_field) != list(solute.force_field):
                problems.append(
                    f"{solute.name}/{domain.name}: force field {list(domain.force_field)}, "
                    f"config says {solute.force_field}"
                )
            if (domain.charges.method, domain.charges.model) != _charges(solute):
                problems.append(
                    f"{solute.name}/{domain.name}: charges {domain.charges.method} "
                    f"{domain.charges.model or ''}, config says {_charges(solute)}"
                )
    names = {s.name for s in config.solutes} | {"solvent", "ions"}
    for extra in sorted(set(records) - names):
        problems.append(f"component {extra!r} is in the build but not in the config")
    solvent = records.get("solvent")
    if solvent is None:
        problems.append("the build has no solvent component")
    else:
        for component in [solvent] + _of_type(context, "ions"):
            for domain in component.domains:
                if list(domain.force_field) != list(config.solvent.force_field):
                    problems.append(
                        f"{component.name}/{domain.name}: force field {list(domain.force_field)}, "
                        f"config solvent says {config.solvent.force_field}"
                    )
        water = _domain(solvent, "water")
        model = config.solvent.water.model.value
        if water is None or (water.charges.method, water.charges.model) != ("water_model", model):
            problems.append(f"water is not built as {model}")
        for cosolvent in config.solvent.co_solvents:
            domain = _domain(solvent, f"co_solvent:{cosolvent.name}")
            if domain is None:
                problems.append(f"co-solvent {cosolvent.name!r} is not in the build")
            elif (domain.charges.method, domain.charges.model) != _charges(cosolvent):
                problems.append(
                    f"co-solvent {cosolvent.name}: charges {domain.charges.method} "
                    f"{domain.charges.model or ''}, config says {_charges(cosolvent)}"
                )
    if problems:
        return [Finding("components.declared", "system", "FAIL", "; ".join(problems))]
    return [
        Finding(
            "components.declared",
            "system",
            "PASS",
            f"{len(config.solutes)} solutes and the solvent are built with the config's force "
            "fields, charge sources and segment IDs",
        )
    ]


@register(
    CheckSpec(
        "composition.counts",
        "system",
        "FAIL",
        "a molecule or chain count that differs from the config",
    )
)
def counts(context: AuditContext) -> list[Finding]:
    config = context.config
    if config is None:
        return _no_config("composition.counts")
    problems, notes = [], []
    records = {c.type: c for c in context.components}
    if config.polymers is not None:
        polymer = records.get("polymer")
        n = polymer.n_molecules if polymer else 0
        notes.append(f"{n} polymer chains")
        if n != config.polymers.count:
            problems.append(f"{n} polymer chains, config count {config.polymers.count}")
    if config.substrate is not None:
        ligand = records.get("small_molecule")
        n = ligand.n_molecules if ligand else 0
        notes.append(f"{n} {config.substrate.name}")
        if n != 1:
            problems.append(f"{n} copies of {config.substrate.name}, config places one")
    solvent = records.get("solvent")
    for cosolvent in config.solvent.co_solvents:
        domain = _domain(solvent, f"co_solvent:{cosolvent.name}") if solvent else None
        n = len(domain.atoms) if domain else 0
        notes.append(f"{n} {cosolvent.name}")
        if cosolvent.count is not None and n != cosolvent.count:
            problems.append(f"{n} {cosolvent.name}, config count {cosolvent.count}")
    if solvent is not None and _domain(solvent, "water") is not None:
        notes.append(f"{len(_domain(solvent, 'water').atoms)} waters")
    if problems:
        return [Finding("composition.counts", "system", "FAIL", "; ".join(problems))]
    return [
        Finding("composition.counts", "system", "PASS", ", ".join(notes) or "no counts to check")
    ]


@register(
    CheckSpec(
        "composition.polymer_sequences",
        "system",
        "INFO",
        "the polymer sequences used, and the monomer fractions they reach",
        component_types=("polymer",),
    )
)
def polymer_sequences(context: AuditContext) -> list[Finding]:
    polymers = _of_type(context, "polymer")
    if not polymers:
        return []
    chains = Counter()
    for domain in polymers[0].domains:
        chains[domain.name.split(":", 1)[-1]] += len(domain.atoms)
    monomers = Counter()
    for sequence, n in chains.items():
        for label in sequence:
            monomers[label] += n
    total = sum(monomers.values()) or 1
    reached = ", ".join(f"{label} {monomers[label] / total:.3f}" for label in sorted(monomers))
    configured = ""
    if context.config is not None and context.config.polymers is not None:
        configured = "; configured " + ", ".join(
            f"{m.label} {m.probability:g}" for m in context.config.polymers.monomers
        )
    detail = f"{len(chains)} sequences over {sum(chains.values())} chains; monomer fractions {reached}{configured}"
    return [Finding("composition.polymer_sequences", polymers[0].name, "INFO", detail)]


def _mass(context: AuditContext, atoms: Any) -> float:
    system = context.realized.system
    from openmm import unit

    return sum(system.getParticleMass(i).value_in_unit(unit.dalton) for i in atoms)


@register(
    CheckSpec(
        "composition.polymer_mass_fraction",
        "system",
        "INFO",
        "the polymer weight fraction of polymer plus solvent the build reached",
        "20 wt% decision for the Paper 1 rerun",
        component_types=("polymer",),
    )
)
def polymer_mass_fraction(context: AuditContext) -> list[Finding]:
    polymers = _of_type(context, "polymer")
    if not polymers:
        return []
    polymer = _mass(context, polymers[0].atom_indices())
    solvent = sum(
        _mass(context, c.atom_indices())
        for c in context.components
        if c.type in ("solvent", "ions")
    )
    fraction = polymer / (polymer + solvent) if polymer + solvent else 0.0
    return [
        Finding(
            "composition.polymer_mass_fraction",
            polymers[0].name,
            "INFO",
            f"{100 * fraction:.2f} wt% polymer of polymer + solvent ({polymer / 1000:.1f} kDa)",
        )
    ]


def _charge_of(context: AuditContext, atoms: Any) -> float:
    particle = context.realized.index.particle
    return sum(particle[i][0] for i in atoms)


@register(
    CheckSpec(
        "ions.neutral",
        "system",
        "FAIL",
        "a net charge left on a system the config asked to neutralize",
    )
)
def neutral(context: AuditContext) -> list[Finding]:
    total = _charge_of(context, range(context.realized.n_atoms))
    neutralize = True if context.config is None else context.config.solvent.ions.neutralize
    if not neutralize:
        return [
            Finding(
                "ions.neutral",
                "system",
                "INFO",
                f"net charge {total:+.4f} e (not neutralized, as configured)",
            )
        ]
    if abs(total) > 1e-3:
        return [
            Finding(
                "ions.neutral", "system", "FAIL", f"net charge {total:+.4f} e with neutralize on"
            )
        ]
    return [Finding("ions.neutral", "system", "PASS", f"net charge {total:+.6f} e")]


def _ion_census(context: AuditContext) -> tuple[int, int, float, float]:
    """(Na+, Cl-, water-equivalents, sum of n z^2 over charged co-solvents)."""
    ions = _of_type(context, "ions")
    na = sum(len(d.atoms) for c in ions for d in c.domains if d.name == "Na+")
    cl = sum(len(d.atoms) for c in ions for d in c.domains if d.name == "Cl-")
    solvent_atoms = [
        i for c in context.components if c.type in ("solvent", "ions") for i in c.atom_indices()
    ]
    water_equivalents = _mass(context, solvent_atoms) / WATER_MASS
    weight = 0.0
    for component in _of_type(context, "solvent"):
        for domain in component.domains:
            if domain.name.startswith("co_solvent:") and domain.atoms:
                z = round(_charge_of(context, domain.atoms[0]))
                weight += z * z * len(domain.atoms)
    return na, cl, water_equivalents, weight


@register(
    CheckSpec(
        "ions.ionic_strength",
        "system",
        "FAIL",
        "an ionic strength more than one ion away from the configured target",
        "nacl_concentration meaning change; buffer matched by ionic strength (Paper 1 rerun)",
    )
)
def ionic_strength(context: AuditContext) -> list[Finding]:
    config = context.config
    if config is None or config.solvent.ions.ionic_strength is None:
        return []
    target = config.solvent.ions.ionic_strength
    na, cl, water_eq, weight = _ion_census(context)
    litres = water_eq / WATER_MOLARITY
    reached = 0.5 * (na + cl + weight) / litres
    one_ion = 0.5 / litres
    solute_charge = _charge_of(
        context,
        [
            i
            for c in context.components
            if c.type not in ("solvent", "ions")
            for i in c.atom_indices()
        ],
    )
    floor = 0.5 * (abs(round(solute_charge)) + weight) / litres
    detail = (
        f"I = {reached:.4f} M from {na} Na+, {cl} Cl- in {water_eq:.0f} water-equivalents "
        f"(target {target} M, one ion = {one_ion:.4f} M)"
    )
    if abs(reached - target) <= one_ion + 1e-12:
        return [Finding("ions.ionic_strength", "system", "PASS", detail)]
    if floor > target:
        return [
            Finding(
                "ions.ionic_strength",
                "system",
                "WARN",
                detail + f"; neutralizing alone needs {floor:.4f} M, above the target",
            )
        ]
    return [Finding("ions.ionic_strength", "system", "FAIL", detail)]


@register(
    CheckSpec(
        "ions.salt",
        "system",
        "FAIL",
        "NaCl pairs more than one pair away from the configured concentration",
        "nacl_concentration meaning change",
    )
)
def salt(context: AuditContext) -> list[Finding]:
    config = context.config
    if config is None or config.solvent.ions.ionic_strength is not None:
        return []
    target = config.solvent.ions.nacl_concentration
    na, cl, water_eq, _ = _ion_census(context)
    expected = round(water_eq * target / WATER_MOLARITY)
    pairs = min(na, cl)
    detail = (
        f"{pairs} NaCl pairs ({na} Na+, {cl} Cl-) for {target} M in {water_eq:.0f} "
        f"water-equivalents; expected {expected}"
    )
    status = "PASS" if abs(pairs - expected) <= 1 else "FAIL"
    return [Finding("ions.salt", "system", status, detail)]


@register(
    CheckSpec(
        "constraints.hbonds",
        "system",
        "FAIL",
        "a bond to hydrogen left flexible at a time step that needs it constrained",
    )
)
def hbonds(context: AuditContext) -> list[Finding]:
    realized = context.realized
    elements = realized.elements
    index = realized.index
    flexible = [
        atoms
        for terms in index.terms.values()
        for atoms, _ in terms.get("bonds", [])
        if "H" in (elements[atoms[0]], elements[atoms[1]])
    ]
    constrained = Counter()
    for terms in index.terms.values():
        for (a, b), _ in terms.get("constraints", []):
            for h, other in ((a, b), (b, a)):
                if elements[h] == "H" and elements[other] != "H":
                    constrained[h] += 1
    free = [i for i, e in enumerate(elements) if e == "H" and constrained[i] == 0]
    step = realized.settings.get("timestep_fs")
    if step is None and context.config is not None:
        step = context.config.simulation_phases.production.time_step
    detail = (
        f"{len(flexible)} flexible bonds to hydrogen, {len(free)} hydrogens without a "
        f"constraint; time step {step} fs"
    )
    if (flexible or free) and (step is None or step > 1.0):
        return [
            Finding(
                "constraints.hbonds",
                "system",
                "FAIL",
                detail,
                [repr(x) for x in (flexible or free)[:3]],
            )
        ]
    return [Finding("constraints.hbonds", "system", "PASS", detail)]


@register(
    CheckSpec(
        "nonbonded.settings",
        "system",
        "FAIL",
        "a cutoff, switching distance, PME or dispersion correction other than the force field's",
    )
)
def nonbonded(context: AuditContext) -> list[Finding]:
    """The realized nonbonded settings are those the force field files define."""
    config = context.config
    if config is None:
        return _no_config("nonbonded.settings")
    from openff.toolkit import ForceField
    from openff.units import unit

    force_field = ForceField(*config.solvent.force_field)
    vdw = force_field["vdW"]
    cutoff = vdw.cutoff.m_as(unit.nanometer)
    switch = cutoff - vdw.switch_width.m_as(unit.nanometer) if vdw.switch_width.m else None
    dispersion = getattr(vdw, "long_range_dispersion", "isotropic") == "isotropic"
    settings = context.realized.settings
    expected = {
        "nonbonded_method": "pme",
        "cutoff_nm": cutoff,
        "switch_nm": switch,
        "dispersion_correction": dispersion,
    }
    wrong = []
    for key, value in expected.items():
        have = settings.get(key)
        same = (
            abs(have - value) < 1e-6
            if isinstance(value, float) and isinstance(have, float)
            else have == value
        )
        if not same:
            wrong.append(f"{key} {have} (force field: {value})")
    detail = ", ".join(f"{k} {settings.get(k)}" for k in expected)
    if wrong:
        return [Finding("nonbonded.settings", "system", "FAIL", "; ".join(wrong))]
    return [Finding("nonbonded.settings", "system", "PASS", detail)]
