"""Write the bond-complete topology that analyses load.

PolyzyMD writes two files when it builds a system. ``solvated_system.pdb``
is for viewers: it carries atom names and coordinates, but its fixed columns
stop at serial 99,999, above which OpenMM writes serials in hex and
MDAnalysis refuses the CONECT records. ``system.xml`` is the OpenMM System:
every particle, bond and constraint, but no names, residues or elements.
Neither alone is a faithful topology for analysis.

``system.prmtop`` is built from both. It holds every atom, residue, element,
mass, charge and bond with no column widths and no atom limit, and MDAnalysis
reads it natively, taking bonds from the file rather than guessing them. The
force-field parameters it also carries are never read by an analysis.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

LOGGER = logging.getLogger(__name__)

ANALYSIS_TOPOLOGY_NAME = "system.prmtop"
VIEWER_TOPOLOGY_NAME = "solvated_system.pdb"
SYSTEM_XML_NAME = "system.xml"


def write_analysis_topology(topology: Any, system: Any, positions: Any, path: Path) -> Path | None:
    """Write ``system.prmtop`` from an OpenMM topology, system and positions.

    Returns the written path, or ``None`` after logging a warning when ParmEd
    cannot convert the system. A build never fails for want of an analysis
    topology; the analysis loader then falls back to the PDB.
    """
    try:
        import parmed
        from openmm import unit

        structure = parmed.openmm.load_topology(topology, system, xyz=positions)
        _type_constrained_bonds(structure, system, unit, parmed)
        path.parent.mkdir(parents=True, exist_ok=True)
        structure.save(str(path), format="amber", overwrite=True)
    except Exception as exc:  # ParmEd raises a mix of its own and builtin errors
        LOGGER.warning(
            "Could not write %s (%s: %s). Analyses will read the PDB topology instead.",
            path,
            type(exc).__name__,
            exc,
        )
        return None
    LOGGER.info("Wrote analysis topology %s", path)
    return path


def _type_constrained_bonds(structure: Any, system: Any, unit: Any, parmed: Any) -> int:
    """Give bonds that exist only as constraints a zero-stiffness type.

    A bond constrained by ``HBonds`` or rigid water has no harmonic term, so
    ParmEd leaves its type empty and the Amber writer refuses the structure.
    The type written here has zero force constant and the constraint length as
    its equilibrium distance. Analyses read the bond, never the parameters.
    """
    lengths: dict[frozenset[int], float] = {}
    for index in range(system.getNumConstraints()):
        first, second, distance = system.getConstraintParameters(index)
        lengths[frozenset((first, second))] = distance.value_in_unit(unit.angstrom)
    typed = 0
    for bond in structure.bonds:
        if bond.type is not None:
            continue
        length = lengths.get(frozenset((bond.atom1.idx, bond.atom2.idx)))
        if length is None:
            length = bond.measure() or 1.0
        bond.type = parmed.BondType(0.0, length, list=structure.bond_types)
        structure.bond_types.append(bond.type)
        typed += 1
    return typed


def rebuild_analysis_topology(working_dir: Path, *, overwrite: bool = False) -> Path | None:
    """Write ``system.prmtop`` for a run built before the file existed.

    Reads ``solvated_system.pdb`` with OpenMM's own PDB reader, which accepts
    the hex serials it writes, and ``system.xml`` for the bonds. Returns the
    existing file untouched unless ``overwrite`` is set.
    """
    output = working_dir / ANALYSIS_TOPOLOGY_NAME
    if output.exists() and not overwrite:
        return output
    pdb_path = working_dir / VIEWER_TOPOLOGY_NAME
    system_path = working_dir / SYSTEM_XML_NAME
    missing = [str(path) for path in (pdb_path, system_path) if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Cannot rebuild {output}; missing {', '.join(missing)}")

    from openmm import XmlSerializer
    from openmm.app import PDBFile

    pdb = PDBFile(str(pdb_path))
    system = XmlSerializer.deserialize(system_path.read_text())
    if pdb.topology.getNumAtoms() != system.getNumParticles():
        raise ValueError(
            f"{pdb_path} has {pdb.topology.getNumAtoms()} atoms but {system_path} has "
            f"{system.getNumParticles()} particles; they are not from the same build"
        )
    return write_analysis_topology(pdb.topology, system, pdb.positions, output)
