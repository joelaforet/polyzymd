"""The system an engine will actually run, in one representation for every check.

A build is checked against the files of the engine its config requests:

- **OpenMM**: the serialized ``system.xml`` and the built PDB;
- **GROMACS**: the ``.top`` (with its ``.itp`` includes) and ``.gro``, read by
  OpenMM's ``GromacsTopFile`` with the settings of the production ``.mdp``,
  because GROMACS applies bond constraints, cutoffs, PME, switching and the
  dispersion correction through run parameters rather than the topology.

Either way the result is an OpenMM ``System`` plus coordinates and run settings,
so every check is written once and reads both engines the same way.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import Any, Literal, Optional

import numpy as np

from polyzymd.audit.terms import ASK_FOR_SUPPORT, SystemIndex

Engine = Literal["openmm", "gromacs"]


@dataclass
class RealizedSystem:
    """An engine's system, its coordinates and the run settings it will use.

    Attributes:
        engine: ``openmm`` or ``gromacs``.
        system: The OpenMM System holding every force-field term.
        topology: The OpenMM Topology (atoms, residues, bonds).
        positions_nm: Atom positions, shape (n_atoms, 3), in nm.
        box_nm: Periodic box vectors as rows, shape (3, 3), in nm.
        settings: Run settings that are not in the System's forces, such as
            ``constraints`` (``none``, ``h-bonds``, ``all-bonds``, ``h-angles``)
            and ``timestep_fs``; empty entries mean the files do not say.
        files: SHA-256 of every file the system was read from, by path.
    """

    engine: Engine
    system: Any
    topology: Any
    positions_nm: np.ndarray
    box_nm: np.ndarray
    settings: dict[str, Any] = field(default_factory=dict)
    files: dict[str, str] = field(default_factory=dict)

    @cached_property
    def index(self) -> SystemIndex:
        """Every term of the system, indexed once."""
        return SystemIndex.from_system(self.system)

    @property
    def n_atoms(self) -> int:
        return self.system.getNumParticles()

    @property
    def elements(self) -> list[str]:
        return [a.element.symbol if a.element is not None else "" for a in self.topology.atoms()]

    @classmethod
    def from_openmm(
        cls, system_xml: Path, structure_pdb: Path, timestep_fs: Optional[float] = None
    ) -> "RealizedSystem":
        """Read an OpenMM build: ``system.xml`` and its PDB.

        Bond constraints are part of an OpenMM System, so ``settings['constraints']``
        is derived from it: ``h-bonds`` when the constraints are exactly the bonds
        to hydrogen (plus rigid water), ``none`` when there are none, ``custom``
        otherwise.
        """
        import openmm
        from openmm import app, unit

        system_xml, structure_pdb = Path(system_xml), Path(structure_pdb)
        system = openmm.XmlSerializer.deserialize(system_xml.read_text())
        pdb = app.PDBFile(str(structure_pdb))
        _require_same_atoms(system, pdb.topology.getNumAtoms(), structure_pdb)
        box = pdb.topology.getPeriodicBoxVectors()
        if box is None:
            box = system.getDefaultPeriodicBoxVectors()
        realized = cls(
            engine="openmm",
            system=system,
            topology=pdb.topology,
            positions_nm=np.array(pdb.positions.value_in_unit(unit.nanometer)),
            box_nm=np.array([[v.value_in_unit(unit.nanometer) for v in row] for row in box]),
            files={str(p): _sha256(p) for p in (system_xml, structure_pdb)},
        )
        realized.settings = _nonbonded_settings(system)
        realized.settings["constraints"] = _constraint_policy(realized)
        realized.settings["timestep_fs"] = timestep_fs
        return realized

    @classmethod
    def from_gromacs(cls, top: Path, gro: Path, mdp: Path) -> "RealizedSystem":
        """Read a GROMACS build: ``.top`` and ``.gro`` under the run settings of ``mdp``.

        Position restraints (``#ifdef POSRES`` blocks) are not defined, as in
        production.

        Only settings an OpenMM System reproduces exactly are accepted: PME
        electrostatics, a plain or potential-switched van der Waals cutoff equal
        to the Coulomb cutoff, and no or ``EnerPres`` dispersion correction.
        Unset keys take GROMACS's defaults.

        Raises:
            ValueError: If the ``.mdp`` asks for anything else (for example
                reaction-field electrostatics, the potential-shift modifier, or
                different Coulomb and van der Waals cutoffs), or an unknown
                constraint setting.
        """
        from openmm import app, unit

        top, gro, mdp = Path(top), Path(gro), Path(mdp)
        options = {key: value.lower() for key, value in parse_mdp(mdp).items()}
        constraints_key = options.get("constraints", "none")
        constraint_map = {
            "none": None,
            "h-bonds": app.HBonds,
            "all-bonds": app.AllBonds,
            "h-angles": app.HAngles,
        }
        if constraints_key not in constraint_map:
            raise ValueError(
                f"{mdp}: constraints = {constraints_key} is not understood. {ASK_FOR_SUPPORT}"
            )
        rvdw = float(options.get("rvdw", 1.0))
        rcoulomb = float(options.get("rcoulomb", 1.0))
        if abs(rvdw - rcoulomb) > 1e-9:
            raise ValueError(
                f"{mdp}: rvdw ({rvdw}) and rcoulomb ({rcoulomb}) differ. {ASK_FOR_SUPPORT}"
            )
        # GROMACS's default coulombtype is a plain cutoff; OpenMM has no exact match.
        coulombtype = options.get("coulombtype", "cut-off")
        if coulombtype != "pme":
            raise ValueError(
                f"{mdp}: coulombtype = {coulombtype} has no exact OpenMM equivalent; "
                f"the audit reads PME only. {ASK_FOR_SUPPORT}"
            )
        vdwtype = options.get("vdwtype", "cut-off").replace("cutoff", "cut-off")
        modifier = options.get("vdw-modifier", "potential-shift")
        if vdwtype != "cut-off" or modifier not in ("potential-switch", "none"):
            raise ValueError(
                f"{mdp}: vdwtype = {vdwtype} with vdw-modifier = {modifier} has no exact "
                f"OpenMM equivalent; use vdw-modifier = Potential-switch or None. {ASK_FOR_SUPPORT}"
            )
        dispcorr = options.get("dispcorr", "no")
        if dispcorr not in ("no", "enerpres"):
            raise ValueError(
                f"{mdp}: DispCorr = {dispcorr} has no exact OpenMM equivalent "
                "(OpenMM corrects energy and pressure together); use EnerPres or no. "
                f"{ASK_FOR_SUPPORT}"
            )

        # The production run's preprocessor defines (define = -DPOSRES -DFLEXIBLE ...)
        # choose blocks of the topology, so they are read as grompp would.
        defines = mdp_defines(parse_mdp(mdp).get("define", ""))
        gro_file = app.GromacsGroFile(str(gro))
        topology_file = app.GromacsTopFile(
            str(top), periodicBoxVectors=gro_file.getPeriodicBoxVectors(), defines=defines
        )
        kwargs: dict[str, Any] = {
            "nonbondedMethod": app.PME,
            "nonbondedCutoff": rcoulomb * unit.nanometer,
            "constraints": constraint_map[constraints_key],
            # GROMACS keeps water rigid (SETTLE) unless the run defines FLEXIBLE.
            "rigidWater": "FLEXIBLE" not in defines,
            "ewaldErrorTolerance": float(options.get("ewald-rtol", 1e-5)),
        }
        if modifier == "potential-switch":
            kwargs["switchDistance"] = float(options.get("rvdw-switch", 0.0)) * unit.nanometer
        system = topology_file.createSystem(**kwargs)
        for force in system.getForces():
            if force.__class__.__name__ == "NonbondedForce":
                force.setUseDispersionCorrection(dispcorr == "enerpres")
        files = {str(p): _sha256(p) for p in (top, gro, mdp)}
        for include in _includes(top):
            files[str(include)] = _sha256(include)
        _require_same_atoms(system, len(gro_file.positions), gro)
        realized = cls(
            engine="gromacs",
            system=system,
            topology=topology_file.topology,
            positions_nm=np.array(gro_file.positions.value_in_unit(unit.nanometer)),
            box_nm=np.array(
                [
                    [v.value_in_unit(unit.nanometer) for v in row]
                    for row in gro_file.getPeriodicBoxVectors()
                ]
            ),
            files=files,
        )
        realized.settings = _nonbonded_settings(system)
        realized.settings["constraints"] = constraints_key
        dt_ps = options.get("dt")
        realized.settings["timestep_fs"] = float(dt_ps) * 1000 if dt_ps else None
        return realized


def mdp_defines(text: str) -> dict[str, str]:
    """Return the preprocessor defines of an ``.mdp`` ``define`` line (``-DPOSRES -DX=1``)."""
    defines = {}
    for token in text.split():
        if token.startswith("-D") and len(token) > 2:
            name, _, value = token[2:].partition("=")
            defines[name] = value or "1"
    return defines


def parse_mdp(path: Path) -> dict[str, str]:
    """Return the settings of a GROMACS ``.mdp`` file, keys lower-cased with ``_`` as ``-``."""
    options: dict[str, str] = {}
    for line in Path(path).read_text().splitlines():
        line = line.split(";", 1)[0].strip()
        if "=" not in line:
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        options[key.lower().replace("_", "-")] = value
    return options


def _includes(top: Path) -> list[Path]:
    """Return the files a ``.top`` includes with paths relative to it, recursively."""
    found: list[Path] = []
    pending = [top]
    while pending:
        current = pending.pop()
        for match in re.finditer(r'^\s*#include\s+"([^"]+)"', current.read_text(), re.MULTILINE):
            path = (current.parent / match.group(1)).resolve()
            if path.is_file() and path not in found:
                found.append(path)
                pending.append(path)
    return found


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _nonbonded_settings(system: Any) -> dict[str, Any]:
    """Return the nonbonded method, cutoff, switching and dispersion settings of ``system``."""
    import openmm
    from openmm import unit

    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            names = {
                openmm.NonbondedForce.NoCutoff: "no-cutoff",
                openmm.NonbondedForce.CutoffNonPeriodic: "cutoff-nonperiodic",
                openmm.NonbondedForce.CutoffPeriodic: "cutoff-periodic",
                openmm.NonbondedForce.Ewald: "ewald",
                openmm.NonbondedForce.PME: "pme",
                openmm.NonbondedForce.LJPME: "ljpme",
            }
            return {
                "nonbonded_method": names.get(force.getNonbondedMethod(), "unknown"),
                "cutoff_nm": force.getCutoffDistance().value_in_unit(unit.nanometer),
                "switch_nm": (
                    force.getSwitchingDistance().value_in_unit(unit.nanometer)
                    if force.getUseSwitchingFunction()
                    else None
                ),
                "dispersion_correction": force.getUseDispersionCorrection(),
                "ewald_tolerance": force.getEwaldErrorTolerance(),
            }
    return {}


def _require_same_atoms(system: Any, n_structure: int, path: Path) -> None:
    if system.getNumParticles() != n_structure:
        raise ValueError(
            f"{path} has {n_structure} atoms but the system has {system.getNumParticles()}"
        )


#: Longest constrained bond to hydrogen (S-H is about 0.134 nm). A constraint to
#: hydrogen longer than this spans an angle, as HAngles adds.
_MAX_H_BOND_NM = 0.145


def _constraint_policy(realized: RealizedSystem) -> str:
    """Name the constraint policy an OpenMM System applies.

    ``h-bonds`` when every constraint is a bond to hydrogen or holds a water
    rigid, ``none`` when there are none, ``custom`` otherwise (heavy-atom
    constraints, or angle constraints such as HAngles adds). Whether every bond
    to hydrogen is constrained is the constraints check's job, which reads the
    bonds from the force field rather than the structure.
    """
    from openmm import unit

    system = realized.system
    if system.getNumConstraints() == 0:
        return "none"
    elements = realized.elements
    atoms = list(realized.topology.atoms())
    for k in range(system.getNumConstraints()):
        a, b, distance = system.getConstraintParameters(k)
        pair = (elements[a], elements[b])
        if "H" not in pair:
            return "custom"
        if pair == ("H", "H"):
            residue = atoms[a].residue
            in_water = residue is atoms[b].residue and sorted(
                x.element.symbol for x in residue.atoms() if x.element is not None
            ) == ["H", "H", "O"]
            if not in_water:
                return "custom"
        elif distance.value_in_unit(unit.nanometer) > _MAX_H_BOND_NM:
            return "custom"
    return "h-bonds"
