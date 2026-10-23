"""Engine checks: the run parameters the engine reads are the config's.

GROMACS takes the dynamics (integrator, thermostat and its time constant,
temperature or annealing schedule, pressure coupling), the time step and
length, bond constraints, velocity generation and the position-restraint
switches from its ``.mdp`` files, not from the topology, so they are read back
and compared with the config. An OpenMM run takes them from the config at run
time; there the run-time system check compares each production segment's
System with the audited build.
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
#: The minimization restraints on protein and ligand heavy atoms (``freeze_solute``).
MINIMIZATION_DEFINE = "POSRES_EM"
#: The GROMACS integrator and tcoupl that run each PolyzyMD thermostat.
THERMOSTATS = {
    "LangevinMiddle": ("sd", "no"),
    "Langevin": ("sd", "no"),
    "NoseHoover": ("md", "nose-hoover"),
    "Andersen": ("md-vv", "andersen"),
}
#: 1 atm in bar, the unit GROMACS takes pressure in.
BAR_PER_ATM = 1.01325


def _close(a: Any, b: float, tolerance: float = 1e-6) -> bool:
    try:
        return abs(float(a) - b) <= tolerance * max(1.0, abs(b))
    except (TypeError, ValueError):
        return False


def _yes(options: dict[str, str], key: str) -> bool:
    return options.get(key, "no").lower() == "yes"


def _restraint_problems(what: str, options: dict[str, str], wanted: set[str]) -> list[str]:
    switched = {d for d in mdp_defines(options.get("define", "")) if d.startswith("POSRES")}
    problems = []
    if wanted - switched:
        problems.append(f"{what}: define lacks {sorted(wanted - switched)}")
    if switched - wanted:
        problems.append(
            f"{what}: define switches on {sorted(switched - wanted)}, "
            "which the config does not restrain here"
        )
    return problems


def _phase_problems(
    what: str,
    options: dict[str, str],
    phase: Any,
    *,
    kelvin: float,
    atm: float,
    first: bool,
    restrained: set[str],
) -> list[str]:
    """Compare one dynamics .mdp with its phase (an equilibration stage or production)."""
    step_fs = phase.time_step or 2.0
    problems = []
    if not _close(options.get("dt"), step_fs / 1000.0):
        problems.append(f"{what}: dt {options.get('dt')} ps, config {step_fs} fs")
    ns = getattr(phase, "resolved_duration", None) or phase.duration
    expected = round(ns * 1e6 / step_fs)
    try:
        steps = int(options.get("nsteps", ""))
    except ValueError:
        steps = None
    if steps is None or abs(steps - expected) > 1:
        problems.append(f"{what}: nsteps {steps}, config {ns} ns at {step_fs} fs = {expected}")

    thermostat = phase.thermostat.value if phase.thermostat else "LangevinMiddle"
    integrator, tcoupl = THERMOSTATS.get(thermostat, THERMOSTATS["LangevinMiddle"])
    if options.get("integrator", "md").lower() != integrator:
        problems.append(
            f"{what}: integrator {options.get('integrator', 'md')}, "
            f"thermostat {thermostat} runs as {integrator}"
        )
    if options.get("tcoupl", "no").lower() != tcoupl:
        problems.append(
            f"{what}: tcoupl {options.get('tcoupl', 'no')}, thermostat {thermostat} runs as {tcoupl}"
        )
    tau = phase.thermostat_timescale or 1.0
    taus = options.get("tau-t", "").split()
    if not taus or any(not _close(t, tau) for t in taus):
        problems.append(f"{what}: tau_t {options.get('tau-t')} ps, config {tau} ps")

    ramp = getattr(phase, "is_temperature_ramping", False)
    final = phase.get_final_temperature() if ramp else kelvin
    temperatures = options.get("ref-t", "").split()
    if not temperatures or any(not _close(t, final) for t in temperatures):
        problems.append(f"{what}: ref_t {options.get('ref-t')}, config {final} K")
    annealing = options.get("annealing", "no").split()
    if not ramp and any(a.lower() != "no" for a in annealing):
        problems.append(f"{what}: annealing {options.get('annealing')} on a constant temperature")
    if ramp:
        problems += _annealing_problems(what, options, phase, final, steps, step_fs)

    coupled = options.get("pcoupl", "no").lower() != "no"
    npt = phase.ensemble.value == "NPT"
    if npt and not coupled:
        problems.append(f"{what}: the config says NPT but pcoupl is no")
    if not npt and coupled:
        problems.append(f"{what}: the config says NVT but pcoupl is {options.get('pcoupl')}")
    if npt and coupled and not _close(options.get("ref-p"), atm * BAR_PER_ATM, 1e-4):
        problems.append(f"{what}: ref_p {options.get('ref-p')} bar, config {atm} atm")
    if options.get("constraints", "none").lower() != "h-bonds":
        problems.append(f"{what}: constraints {options.get('constraints')}, PolyzyMD runs h-bonds")

    # Velocities are drawn once, at the start of the first dynamics; every later
    # phase continues from the one before.
    if _yes(options, "gen-vel") != first:
        problems.append(
            f"{what}: gen_vel {options.get('gen-vel', 'no')}; velocities are drawn only "
            "at the start of the first equilibration stage"
        )
    if _yes(options, "continuation") == first:
        problems.append(
            f"{what}: continuation {options.get('continuation', 'no')}; every phase "
            "but the first continues from the one before"
        )
    wanted = {POSRES_DEFINES[g] for g in restrained if g in POSRES_DEFINES}
    return problems + _restraint_problems(what, options, wanted)


def _annealing_problems(
    what: str, options: dict[str, str], stage: Any, final: float, steps: int | None, step_fs: float
) -> list[str]:
    """A ramp runs as one annealing schedule from the start to the final temperature."""
    start = stage.get_start_temperature()
    try:
        times = [float(t) for t in options.get("annealing-time", "").split()]
        temps = [float(t) for t in options.get("annealing-temp", "").split()]
    except ValueError:
        times, temps = [], []
    if options.get("annealing", "no").lower() != "single" or not temps or len(times) != len(temps):
        return [f"{what}: the config ramps {start} -> {final} K but the .mdp has no annealing"]
    problems = []
    if not (_close(temps[0], start) and _close(temps[-1], final)):
        problems.append(
            f"{what}: annealing {temps[0]} -> {temps[-1]} K, config {start} -> {final} K"
        )
    if any(t < start - 1e-6 or t > final + 1e-6 for t in temps) or temps != sorted(temps):
        problems.append(f"{what}: annealing temperatures leave the ramp {start} -> {final} K")
    if steps is not None and times[-1] > steps * step_fs / 1000.0 + 1e-6:
        problems.append(f"{what}: annealing runs past the end of the stage")
    return problems


def _minimization_problems(path: Path, config: Any, components) -> list[str]:
    options = parse_mdp(path)
    problems = []
    if options.get("integrator", "").lower() != "steep":
        problems.append(
            f"em.mdp: integrator {options.get('integrator')}, PolyzyMD minimizes with steep"
        )
    solute = any(getattr(c, "type", None) in ("protein", "small_molecule") for c in components)
    frozen = config.simulation_phases.minimization.freeze_solute and solute
    return problems + _restraint_problems(
        "em.mdp", options, {MINIMIZATION_DEFINE} if frozen else set()
    )


@register(
    CheckSpec(
        "engine.realized",
        "system",
        "FAIL",
        "GROMACS run parameters (.mdp) other than the config's: integrator, thermostat, time "
        "step, length, temperature or annealing, pressure coupling, bond constraints, velocity "
        "generation or position-restraint switches",
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
    folder = prod_mdp.parent
    phases = config.simulation_phases
    kelvin = config.thermodynamics.temperature
    atm = config.thermodynamics.pressure
    stages = phases.equilibration_stages or []

    problems = []
    em_mdp = folder / "em.mdp"
    if em_mdp.is_file():
        problems += _minimization_problems(em_mdp, config, context.components)
    else:
        problems.append("the GROMACS build has no em.mdp")
    # Stage files are named eq_NN_<stage name>.mdp in config order.
    expected = [f"eq_{k + 1:02d}_{stage.name}.mdp" for k, stage in enumerate(stages)]
    found = sorted(p.name for p in folder.glob("eq_*.mdp"))
    if found != sorted(expected):
        problems.append(f"equilibration .mdp files {found}, the config's stages give {expected}")
    for k, (name, stage) in enumerate(zip(expected, stages)):
        if (folder / name).is_file():
            problems += _phase_problems(
                name,
                parse_mdp(folder / name),
                stage,
                kelvin=stage.temperature or kelvin,
                atm=atm,
                first=k == 0,
                restrained={p.group for p in stage.position_restraints},
            )
    # Production continues from equilibration, or from minimization when there
    # is none, so it never draws velocities.
    problems += _phase_problems(
        "prod.mdp",
        parse_mdp(prod_mdp),
        phases.production,
        kelvin=kelvin,
        atm=atm,
        first=False,
        restrained=set(),
    )
    if problems:
        return [Finding("engine.realized", "system", "FAIL", "; ".join(problems))]
    return [
        Finding(
            "engine.realized",
            "system",
            "PASS",
            f"em.mdp, prod.mdp and {len(stages)} equilibration .mdp files carry the config's "
            "dynamics, time steps, lengths, temperatures, pressure coupling, h-bond constraints, "
            "velocity generation and restraint switches",
        )
    ]
