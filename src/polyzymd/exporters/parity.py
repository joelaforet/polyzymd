"""Manual OpenMM/GROMACS parity evidence helpers."""

from __future__ import annotations

import copy
import json
import subprocess
import tempfile
from pathlib import Path
from typing import Any

DEFAULT_GROMACS_ENERGY_TERMS = (
    "Bond",
    "Angle",
    "Proper-Dih.",
    "LJ-14",
    "Coulomb-14",
    "LJ-(SR)",
    "Coulomb-(SR)",
    "Coul.-recip.",
    "Potential",
)


def write_openmm_gromacs_parity_report(
    *,
    openmm_system: Any,
    openmm_positions: Any,
    gromacs_edr_path: Path | str,
    output_path: Path | str,
    openmm_periodic_box_vectors: Any | None = None,
    gmx_command: str = "gmx",
    energy_terms: tuple[str, ...] = DEFAULT_GROMACS_ENERGY_TERMS,
    platform_name: str = "Reference",
) -> Path:
    """Write an identical-position OpenMM/GROMACS energy parity report.

    This helper is intended for manual production evidence. It reads energies
    from a GROMACS rerun ``.edr`` file and evaluates the OpenMM System at the
    same positions used to create that rerun. Metrics are normalized by particle
    count so large-system reports remain easy to compare across reruns.

    Parameters
    ----------
    openmm_system : Any
        Authoritative OpenMM System.
    openmm_positions : Any
        Positions matching the GROMACS rerun coordinate file.
    gromacs_edr_path : Path or str
        GROMACS energy file from ``mdrun -rerun``.
    output_path : Path or str
        Destination JSON report.
    openmm_periodic_box_vectors : Any or None, optional
        Periodic box vectors matching the rerun coordinates, by default ``None``.
    gmx_command : str, optional
        GROMACS executable, by default ``"gmx"``.
    energy_terms : tuple of str, optional
        GROMACS energy terms to extract.
    platform_name : str, optional
        OpenMM platform for deterministic single-point evaluation, by default ``"Reference"``.

    Returns
    -------
    Path
        Written JSON report path.
    """

    gromacs_terms = extract_gromacs_energy_terms(
        gromacs_edr_path,
        energy_terms=energy_terms,
        gmx_command=gmx_command,
    )
    openmm_terms = evaluate_openmm_force_groups(
        openmm_system,
        openmm_positions,
        periodic_box_vectors=openmm_periodic_box_vectors,
        platform_name=platform_name,
    )
    particle_count = int(openmm_system.getNumParticles())
    total_delta = openmm_terms["total_potential_kj_mol"] - gromacs_terms.get("Potential", 0.0)
    bonded_14_gromacs = sum(
        gromacs_terms.get(term, 0.0)
        for term in ("Bond", "Angle", "Proper-Dih.", "Improper-Dih.", "LJ-14", "Coulomb-14")
    )
    report = {
        "schema_version": 1,
        "gromacs_edr_path": str(gromacs_edr_path),
        "particle_count": particle_count,
        "openmm": openmm_terms,
        "gromacs": gromacs_terms,
        "normalized_metrics": {
            "total_potential_delta_kj_mol": total_delta,
            "total_potential_delta_kj_mol_per_particle": total_delta / particle_count,
            "gromacs_explicit_bonded_14_kj_mol": bonded_14_gromacs,
            "gromacs_explicit_bonded_14_kj_mol_per_particle": bonded_14_gromacs / particle_count,
            "gromacs_pme_reciprocal_kj_mol": gromacs_terms.get("Coul.-recip."),
            "gromacs_pme_reciprocal_kj_mol_per_particle": (
                gromacs_terms.get("Coul.-recip.", 0.0) / particle_count
            ),
        },
    }
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return path


def evaluate_openmm_force_groups(
    openmm_system: Any,
    openmm_positions: Any,
    *,
    periodic_box_vectors: Any | None = None,
    platform_name: str = "Reference",
) -> dict[str, Any]:
    """Evaluate OpenMM total potential and per-force-class energies."""

    from openmm import LangevinIntegrator, Platform, unit
    from openmm.app import Simulation, Topology

    system = copy.deepcopy(openmm_system)
    force_groups: dict[int, str] = {}
    for group, force in enumerate(system.getForces()):
        force.setForceGroup(group)
        force_groups[group] = force.__class__.__name__
    integrator = LangevinIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.001 * unit.picoseconds
    )
    platform = Platform.getPlatformByName(platform_name)
    topology = Topology()
    chain = topology.addChain("A")
    residue = topology.addResidue("DUM", chain)
    for index in range(system.getNumParticles()):
        topology.addAtom(f"X{index}", None, residue)
    simulation = Simulation(topology, system, integrator, platform)
    if periodic_box_vectors is not None:
        simulation.context.setPeriodicBoxVectors(*periodic_box_vectors)
    simulation.context.setPositions(openmm_positions)
    total = _state_energy_kj_mol(simulation.context.getState(getEnergy=True))
    by_force: dict[str, float] = {}
    for group, name in force_groups.items():
        energy = _state_energy_kj_mol(simulation.context.getState(getEnergy=True, groups={group}))
        by_force[name] = by_force.get(name, 0.0) + energy
    return {"total_potential_kj_mol": total, "by_force_class_kj_mol": by_force}


def extract_gromacs_energy_terms(
    edr_path: Path | str,
    *,
    energy_terms: tuple[str, ...] = DEFAULT_GROMACS_ENERGY_TERMS,
    gmx_command: str = "gmx",
) -> dict[str, float]:
    """Extract final-frame GROMACS energy terms from an ``.edr`` file."""

    with tempfile.TemporaryDirectory() as tmpdir:
        xvg_path = Path(tmpdir) / "energies.xvg"
        selection = "\n".join(energy_terms) + "\n0\n"
        subprocess.run(
            [gmx_command, "energy", "-f", str(edr_path), "-o", str(xvg_path)],
            input=selection,
            text=True,
            capture_output=True,
            check=True,
        )
        return parse_gromacs_energy_xvg(xvg_path, energy_terms)


def parse_gromacs_energy_xvg(path: Path | str, energy_terms: tuple[str, ...]) -> dict[str, float]:
    """Parse the final data row from a GROMACS energy XVG file."""

    data_rows = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "@")):
            continue
        data_rows.append([float(value) for value in stripped.split()])
    if not data_rows:
        raise ValueError(f"GROMACS energy XVG contains no data rows: {path}")
    final = data_rows[-1]
    values = final[1:]
    if len(values) < len(energy_terms):
        raise ValueError(
            f"GROMACS energy XVG contains {len(values)} values for {len(energy_terms)} terms"
        )
    return {term: values[index] for index, term in enumerate(energy_terms)}


def _state_energy_kj_mol(state: Any) -> float:
    """Return an OpenMM State potential energy in kJ/mol."""

    from openmm import unit

    return float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
