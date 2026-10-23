"""Engine checks: the run parameters the engine reads are the config's.

GROMACS takes the time step, length, temperature, pressure coupling, bond
constraints and position-restraint switches from its ``.mdp`` files, not from
the topology, so they are read back and compared with the config. An OpenMM
run takes them from the config at run time; there the run-time system check
compares each production segment's System with the audited build.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from polyzymd.audit.findings import AuditContext, Finding
from polyzymd.audit.realized import mdp_defines, parse_mdp
from polyzymd.audit.registry import CheckSpec, register

#: The ``define`` that switches on each position-restraint group in a GROMACS build.
#: Written here, not imported from the exporter, so the check does not share its code.
POSRES_DEFINES = {
    "protein_heavy": "POSRES_PROTEIN",
    "protein_backbone": "POSRES_PROTEIN",
    "protein_calpha": "POSRES_PROTEIN",
    "ligand_heavy": "POSRES_LIGAND",
    "polymer_heavy": "POSRES_POLYMER",
}
#: 1 atm in bar, the unit GROMACS takes pressure in.
BAR_PER_ATM = 1.01325


def _close(a: Any, b: float, tolerance: float = 1e-6) -> bool:
    try:
        return abs(float(a) - b) <= tolerance * max(1.0, abs(b))
    except (TypeError, ValueError):
        return False


def _phase_problems(
    what: str,
    options: dict[str, str],
    *,
    step_fs: float,
    ns: float | None,
    kelvin: float | None,
    npt: bool,
    atm: float,
    restrained: set[str],
) -> list[str]:
    problems = []
    if not _close(options.get("dt"), step_fs / 1000.0):
        problems.append(f"{what}: dt {options.get('dt')} ps, config {step_fs} fs")
    if ns is not None:
        expected = int(ns * 1e6 / step_fs)
        steps = int(options.get("nsteps", -1))
        if abs(steps - expected) > 1:
            problems.append(f"{what}: nsteps {steps}, config {ns} ns at {step_fs} fs = {expected}")
    if kelvin is not None and "annealing" not in options:
        temperatures = options.get("ref-t", "").split()
        if not temperatures or any(not _close(t, kelvin) for t in temperatures):
            problems.append(f"{what}: ref_t {options.get('ref-t')}, config {kelvin} K")
    coupled = options.get("pcoupl", "no").lower() != "no"
    if npt and not coupled:
        problems.append(f"{what}: the config says NPT but pcoupl is no")
    if not npt and coupled:
        problems.append(f"{what}: the config says NVT but pcoupl is {options.get('pcoupl')}")
    if npt and coupled and not _close(options.get("ref-p"), atm * BAR_PER_ATM, 1e-4):
        problems.append(f"{what}: ref_p {options.get('ref-p')} bar, config {atm} atm")
    if options.get("constraints", "none").lower() != "h-bonds":
        problems.append(f"{what}: constraints {options.get('constraints')}, PolyzyMD runs h-bonds")
    switched = set(mdp_defines(options.get("define", "")))
    wanted = {POSRES_DEFINES[g] for g in restrained if g in POSRES_DEFINES}
    if wanted - switched:
        problems.append(
            f"{what}: restrains {sorted(restrained)} but define lacks {sorted(wanted - switched)}"
        )
    stray = {d for d in switched if d.startswith("POSRES") and d != "POSRES_EM"} - wanted
    if stray:
        problems.append(
            f"{what}: define switches on {sorted(stray)}, which the config does not restrain"
        )
    return problems


@register(
    CheckSpec(
        "engine.realized",
        "system",
        "FAIL",
        "GROMACS run parameters (.mdp) other than the config's: time step, length, temperature, "
        "pressure coupling, bond constraints or position-restraint switches",
        "testers on GROMACS",
    )
)
def realized(context: AuditContext) -> list[Finding]:
    if context.realized.engine != "gromacs":
        return [
            Finding(
                "engine.realized",
                "system",
                "INFO",
                "OpenMM takes its run parameters from the config at run time; each production "
                "segment's System is compared with this build when it starts",
            )
        ]
    config = context.config
    if config is None:
        return [Finding("engine.realized", "system", "INFO", "not checked: no config")]
    prod_mdp = next((Path(p) for p in context.realized.files if p.endswith("prod.mdp")), None)
    if prod_mdp is None:
        return [Finding("engine.realized", "system", "FAIL", "the GROMACS build has no prod.mdp")]
    phases = config.simulation_phases
    production = phases.production
    problems = _phase_problems(
        "prod.mdp",
        parse_mdp(prod_mdp),
        step_fs=production.time_step,
        ns=production.duration,
        kelvin=config.thermodynamics.temperature,
        npt=production.ensemble.value == "NPT",
        atm=config.thermodynamics.pressure,
        restrained=set(),
    )
    stages = phases.equilibration_stages or []
    files = sorted(prod_mdp.parent.glob("eq_*.mdp"))
    if len(files) != len(stages):
        problems.append(
            f"{len(files)} equilibration .mdp files for {len(stages)} configured stages"
        )
    for path, stage in zip(files, stages):
        ramp = stage.temperature is None
        problems += _phase_problems(
            path.name,
            parse_mdp(path),
            step_fs=stage.time_step or 2.0,
            ns=None if ramp else stage.resolved_duration,
            kelvin=None if ramp else stage.temperature,
            npt=stage.ensemble.value == "NPT",
            atm=config.thermodynamics.pressure,
            restrained={p.group for p in stage.position_restraints},
        )
    if problems:
        return [Finding("engine.realized", "system", "FAIL", "; ".join(problems))]
    return [
        Finding(
            "engine.realized",
            "system",
            "PASS",
            f"prod.mdp and {len(files)} equilibration .mdp files carry the config's time steps, "
            "lengths, temperatures, pressure coupling, h-bond constraints and restraint switches",
        )
    ]
