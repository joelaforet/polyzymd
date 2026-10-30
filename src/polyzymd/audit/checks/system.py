"""System checks: the build holds what the config asked for, as a whole.

These compare the realized system with the config itself (not with the
builder's records of it): the components and their declared sources, counts,
net charge, salt or ionic strength, the constraint policy and the nonbonded
settings.
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Optional

import numpy as np

from polyzymd.audit.findings import AuditContext, Finding
from polyzymd.audit.registry import CheckSpec, register
from polyzymd.core.components import is_attachment_segid

#: Molar concentration of pure water, as the builder sizes salt and ionic strength.
WATER_MOLARITY = 55.5
#: Mass of one water molecule (Da), the unit of the builder's water-equivalents
#: (the same for every water model it builds: TIP3P and OPC3 have the same atoms).
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
            expected = _declared_source(solute, domain)
            if expected is None:
                problems.append(
                    f"{solute.name}/{domain.name}: a {_site_type(domain)} site, but the config "
                    "has no modification_parameters for it"
                )
                continue
            force_field, charges = expected
            if list(domain.force_field) != list(force_field):
                problems.append(
                    f"{solute.name}/{domain.name}: force field {list(domain.force_field)}, "
                    f"config says {list(force_field)}"
                )
            if (domain.charges.method, domain.charges.model) != charges:
                problems.append(
                    f"{solute.name}/{domain.name}: charges {domain.charges.method} "
                    f"{domain.charges.model or ''}, config says {charges}"
                )
    names = {s.name for s in config.solutes} | {"solvent", "ions"}
    for extra in sorted(set(records) - names):
        problems.append(f"component {extra!r} is in the build but not in the config")
    solvent = records.get("solvent")
    if solvent is None:
        problems.append("the build has no solvent component")
    else:
        from polyzymd.config.schema import water_force_fields

        model = config.solvent.water.model.value
        # Water and ions take the water model's own file (opc3.offxml) after
        # the solvent force field; co-solvents take the solvent force field.
        water_files = water_force_fields(config.solvent.force_field, model)
        for component in [solvent] + _of_type(context, "ions"):
            for domain in component.domains:
                expected = (
                    list(config.solvent.force_field)
                    if domain.name.startswith("co_solvent:")
                    else water_files
                )
                if list(domain.force_field) != expected:
                    problems.append(
                        f"{component.name}/{domain.name}: force field {list(domain.force_field)}, "
                        f"config solvent and water model say {expected}"
                    )
        water = _domain(solvent, "water")
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
            elif domain.charges.method == "preset":
                from pathlib import Path

                shipped = f"{cosolvent.name.lower().strip()}.sdf"
                if Path(domain.charges.source_file or "").name != shipped:
                    problems.append(
                        f"co-solvent {cosolvent.name}: preset charges from "
                        f"{domain.charges.source_file}, not the shipped {shipped}"
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


def _grafts(context: AuditContext) -> list[tuple[str, list[int]]]:
    """Return each polymer chain grafted onto a component: its segment and atoms."""
    grafts = []
    for component in context.components:
        for segid, ranges in sorted(getattr(component, "segments", {}).items()):
            if is_attachment_segid(segid, kinds="P"):
                grafts.append((segid, [i for r in ranges for i in r]))
    return grafts


def _graft_sequences(context: AuditContext) -> list[tuple[str, str]]:
    """Return (segment, drawn sequence) of each grafted chain the build recorded."""
    return [
        (chain["segment"], chain["sequence"])
        for component in context.components
        for chain in component.provenance.get("grafted_chains", [])
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
    grafted = ""
    if _grafts(context):
        grafted = "; grafted " + (
            ", ".join(f"{name} {sequence}" for name, sequence in _graft_sequences(context))
            or ", ".join(segid for segid, _ in _grafts(context))
        )
    if not polymers:
        if grafted:
            return [Finding("composition.polymer_sequences", "system", "INFO", grafted[2:])]
        return []
    chains = Counter()
    for domain in polymers[0].domains:
        chains[domain.name.split(":", 1)[-1]] += len(domain.atoms)
    config_polymers = context.config.polymers if context.config is not None else None
    # A dynamic chain's sequence starts with its head cap and ends with its tail cap.
    capped = (
        config_polymers is not None and getattr(config_polymers, "end_groups", None) is not None
    )
    monomers = Counter()
    ends = {"head": Counter(), "tail": Counter()}
    for sequence, n in chains.items():
        if capped and len(sequence) >= 3:
            ends["head"][sequence[0]] += n
            ends["tail"][sequence[-1]] += n
            sequence = sequence[1:-1]
        for label in sequence:
            monomers[label] += n
    total = sum(monomers.values()) or 1
    reached = ", ".join(f"{label} {monomers[label] / total:.3f}" for label in sorted(monomers))
    configured = ""
    if config_polymers is not None:
        configured = "; configured " + ", ".join(
            f"{m.label} {m.probability:g}" for m in config_polymers.monomers
        )
    caps = ""
    if capped:
        # Caps that are a monomer's end form carry that monomer's label; with
        # them counted, the fractions are those of every monomer of the chains.
        labels = {m.label for m in config_polymers.monomers}
        units = Counter(monomers)
        for counter in ends.values():
            for label, count in counter.items():
                if label in labels:
                    units[label] += count
        if units != monomers:
            whole = sum(units.values()) or 1
            caps += "; with end units " + ", ".join(
                f"{label} {units[label] / whole:.3f}" for label in sorted(units)
            )
        caps += "; caps " + "; ".join(
            f"{end} " + ", ".join(f"{label} {count}" for label, count in sorted(counter.items()))
            for end, counter in ends.items()
        )
    detail = (
        f"{len(chains)} sequences over {sum(chains.values())} {'soluble ' if grafted else ''}"
        f"chains; monomer fractions "
        f"{reached}{configured}{caps}{grafted}"
    )
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
        component_types=("polymer",),
    )
)
def polymer_mass_fraction(context: AuditContext) -> list[Finding]:
    polymers = _of_type(context, "polymer")
    grafts = _grafts(context)
    if not polymers and not grafts:
        return []
    soluble = _mass(context, polymers[0].atom_indices()) if polymers else 0.0
    grafted = sum(_mass(context, atoms) for _, atoms in grafts)
    polymer = soluble + grafted
    solvent = sum(
        _mass(context, c.atom_indices())
        for c in context.components
        if c.type in ("solvent", "ions")
    )
    fraction = polymer / (polymer + solvent) if polymer + solvent else 0.0
    return [
        Finding(
            "composition.polymer_mass_fraction",
            polymers[0].name if polymers else "system",
            "INFO",
            f"{100 * fraction:.2f} wt% polymer of polymer + solvent ({polymer / 1000:.1f} kDa"
            + (
                f"; {soluble / 1000:.1f} kDa soluble, {grafted / 1000:.1f} kDa grafted in "
                f"{len(grafts)} chain(s))"
                if grafts
                else ")"
            ),
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


#: Daltons in one nm^3 of a 1 g/mL fluid.
DALTON_PER_NM3_AT_1G_PER_ML = 602.214076


def _water_equivalents(context: AuditContext) -> float:
    """The solvent the salt and ionic strength refer to, as the builder sizes it.

    The builder fills the box to its target density: (box volume x density -
    solute mass) / water mass. Co-solvents given as a count or concentration
    are added on top, so they are not part of it.
    """
    box = np.asarray(context.realized.box_nm, dtype=float)
    volume = abs(float(np.linalg.det(box)))
    density = 1.0 if context.config is None else context.config.solvent.box.target_density
    solute = _mass(
        context,
        [
            i
            for c in context.components
            if c.type not in ("solvent", "ions")
            for i in c.atom_indices()
        ],
    )
    return (volume * density * DALTON_PER_NM3_AT_1G_PER_ML - solute) / WATER_MASS


def _co_solvent_charge(context: AuditContext) -> int:
    """Net charge of the charged co-solvent molecules, such as SDS anions."""
    return sum(
        round(_charge_of(context, domain.atoms[0])) * len(domain.atoms)
        for component in _of_type(context, "solvent")
        for domain in component.domains
        if domain.name.startswith("co_solvent:") and domain.atoms
    )


def _ion_census(context: AuditContext) -> tuple[int, int, float, float, float]:
    """(Na+, Cl-, water-equivalents, sum n z^2 of charged co-solvents, charge of all but ions)."""
    ions = _of_type(context, "ions")
    na = sum(len(d.atoms) for c in ions for d in c.domains if d.name == "Na+")
    cl = sum(len(d.atoms) for c in ions for d in c.domains if d.name == "Cl-")
    weight = 0.0
    for component in _of_type(context, "solvent"):
        for domain in component.domains:
            if domain.name.startswith("co_solvent:") and domain.atoms:
                z = round(_charge_of(context, domain.atoms[0]))
                weight += z * z * len(domain.atoms)
    others = _charge_of(
        context, [i for c in context.components if c.type != "ions" for i in c.atom_indices()]
    )
    return na, cl, _water_equivalents(context), weight, others


@register(
    CheckSpec(
        "ions.ionic_strength",
        "system",
        "FAIL",
        "an ionic strength more than one ion away from the configured target",
    )
)
def ionic_strength(context: AuditContext) -> list[Finding]:
    config = context.config
    if config is None:
        return _no_config("ions.ionic_strength")
    if config.solvent.ions.ionic_strength is None:
        return []
    target = config.solvent.ions.ionic_strength
    na, cl, water_eq, weight, others = _ion_census(context)
    litres = water_eq / WATER_MOLARITY
    reached = 0.5 * (na + cl + weight) / litres
    one_ion = 0.5 / litres
    # The least the build can reach: the ions that neutralize everything else.
    floor = 0.5 * (abs(round(others)) + weight) / litres
    detail = (
        f"I = {reached:.4f} M from {na} Na+, {cl} Cl- in {water_eq:.0f} water-equivalents "
        f"(target {target} M, one ion = {one_ion:.4f} M)"
    )
    if abs(reached - target) <= one_ion + 1e-12:
        return [Finding("ions.ionic_strength", "system", "PASS", detail)]
    if floor > target and reached <= floor + one_ion + 1e-12:
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
    )
)
def salt(context: AuditContext) -> list[Finding]:
    config = context.config
    if config is None:
        return _no_config("ions.salt")
    if config.solvent.ions.ionic_strength is not None:
        return []
    target = config.solvent.ions.nacl_concentration
    na, cl, water_eq, _, others = _ion_census(context)
    expected = round(water_eq * target / WATER_MOLARITY)
    # A charged co-solvent (SDS) comes with its own counter-ions; they are not salt.
    co_solvent = _co_solvent_charge(context)
    counter = abs(co_solvent)
    # Salt pairs: the ions beyond those counter-ions and those that neutralize the rest.
    neutralizing = abs(round(others - co_solvent)) if config.solvent.ions.neutralize else 0
    pairs = (na + cl - counter - neutralizing) // 2
    of_co_solvents = f"{counter} counter-ions of charged co-solvents, " if counter else ""
    detail = (
        f"{pairs} NaCl pairs ({na} Na+, {cl} Cl-: {of_co_solvents}{neutralizing} neutralizing) "
        f"for {target} M in {water_eq:.0f} water-equivalents; expected {expected}"
    )
    status = "PASS" if abs(pairs - expected) <= 1 else "FAIL"
    return [Finding("ions.salt", "system", status, detail)]


@register(
    CheckSpec(
        "constraints.hbonds",
        "system",
        "FAIL",
        "a bond to hydrogen left flexible at a time step that needs it constrained, or a step "
        "above 2 fs without hydrogen mass repartitioning",
    )
)
def hbonds(context: AuditContext) -> list[Finding]:
    from openmm import unit

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
    steps = [realized.settings.get("timestep_fs")]
    if context.config is not None:
        phases = context.config.simulation_phases
        steps.append(phases.production.time_step)
        steps += [stage.time_step or 2.0 for stage in phases.equilibration_stages or []]
    step = max((s for s in steps if s is not None), default=None)
    masses = [
        realized.system.getParticleMass(i).value_in_unit(unit.dalton)
        for i, e in enumerate(elements)
        if e == "H"
    ]
    lightest = min(masses, default=None)
    detail = (
        f"{len(flexible)} flexible bonds to hydrogen, {len(free)} hydrogens without a "
        f"constraint; longest time step {step} fs; lightest hydrogen {lightest} Da"
    )
    if (flexible or free) and (step is None or step > 1.0):
        examples = [repr(x) for x in (flexible or free)[:3]]
        return [Finding("constraints.hbonds", "system", "FAIL", detail, examples)]
    if step is not None and step > 2.0 + 1e-9 and lightest is not None and lightest < 2.0:
        return [
            Finding(
                "constraints.hbonds",
                "system",
                "FAIL",
                detail + "; a step above 2 fs needs hydrogen mass repartitioning",
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

    from polyzymd.config.schema import water_force_fields

    force_field = ForceField(
        *water_force_fields(config.solvent.force_field, config.solvent.water.model)
    )
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


@register(
    CheckSpec(
        "composition.density",
        "system",
        "INFO",
        "the density the box reached, against the config's target",
    )
)
def density(context: AuditContext) -> list[Finding]:
    box = np.asarray(context.realized.box_nm, dtype=float)
    volume = abs(float(np.linalg.det(box)))
    total = _mass(context, range(context.realized.n_atoms))
    reached = total / (volume * DALTON_PER_NM3_AT_1G_PER_ML)
    detail = f"{reached:.4f} g/mL in {volume:.1f} nm^3"
    if context.config is not None:
        detail += (
            f" (target {context.config.solvent.box.target_density} g/mL; co-solvents given "
            "as a count or concentration are added on top of it)"
        )
    return [Finding("composition.density", "system", "INFO", detail)]


def _site_type(domain: Any) -> str | None:
    """The site type of a modification-site domain, or None for any other domain."""
    reference = domain.reference or {}
    return reference.get("site_type") if reference.get("site") else None


def _declared_source(solute: Any, domain: Any) -> tuple[list[str], tuple] | None:
    """The force field and charges the config names for ``domain`` of ``solute``.

    A modification site takes its type's set from ``modification_parameters``;
    every other domain takes the solute's own force field and charges.
    """
    site_type = _site_type(domain)
    if site_type is None:
        return list(solute.force_field), _charges(solute)
    sets = getattr(solute, "modification_parameters", None)
    chosen = getattr(sets, site_type, None) if sets is not None else None
    if chosen is None:
        return None
    return list(chosen.force_field), (chosen.charges.method, chosen.charges.model)

