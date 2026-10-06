"""Read a GROMACS run whose ``prod.tpr`` MDAnalysis cannot parse.

MDAnalysis reads a run input (``.tpr``) only up to the file version it
knows: MDAnalysis 2.10 reads GROMACS 2025 files (tpx 137) but not GROMACS
2026 files (tpx 138). For such a run, :func:`universe_from_gromacs_top`
builds the universe from the run's GROMACS topology (``.top``) instead,
parsed by MDAnalysis's ``ITPParser``, and lays it out as MDAnalysis's
``TPRParser`` lays out a TPR: the same atoms, residue numbers, segments,
chain IDs, charges, masses, elements, bonds, angles and dihedrals. An
analysis therefore measures the same thing whichever GROMACS version wrote
the run.

:func:`apply_build_chain_ids` then gives either kind of GROMACS universe
the chain IDs PolyzyMD assigned when it built the system (A protein, B
substrate, C polymer), read from the build's ``solvated_system.pdb``.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe

#: Particle types GROMACS allows in the ``ptype`` column of ``[ atomtypes ]``.
_PTYPES = frozenset({"A", "S", "V", "D", "B"})
_INCLUDE = re.compile(r'^\s*#include\s+"([^"]+)"')
_SECTION = re.compile(r"^\s*\[\s*([A-Za-z_]+)\s*\]")


def tpr_unsupported(error: BaseException) -> bool:
    """Return whether ``error`` is MDAnalysis refusing a TPR file version it cannot read.

    ``MDAnalysis.Universe`` reports it as a ``ValueError`` raised while
    handling the ``NotImplementedError`` of the TPR parser.
    """
    return isinstance(error, ValueError) and isinstance(error.__context__, NotImplementedError)


def system_prefix(config: Any) -> str:
    """Return the prefix of the files PolyzyMD writes for ``config``'s GROMACS run.

    The enzyme name, then the polymer type prefix when polymers are enabled,
    joined with ``_``; ``"system"`` when there is neither.
    """
    parts: list[str] = []
    enzyme_name = getattr(getattr(config, "enzyme", None), "name", None)
    if isinstance(enzyme_name, str) and enzyme_name:
        parts.append(enzyme_name)
    polymers = getattr(config, "polymers", None)
    polymer_prefix = getattr(polymers, "type_prefix", None)
    if getattr(polymers, "enabled", False) is True and isinstance(polymer_prefix, str):
        parts.append(polymer_prefix)
    return "_".join(parts) if parts else "system"


def topology_name(config: Any) -> str:
    """Return the file name of the run's ``.top``: ``gromacs.analysis_topology``, or ``<prefix>.top``."""
    chosen = getattr(getattr(config, "gromacs", None), "analysis_topology", None)
    return chosen if isinstance(chosen, str) and chosen else f"{system_prefix(config)}.top"


def run_input_files(working_dir: str | Path, config: Any) -> list[Path]:
    """Return the GROMACS input files of a run that exist, by the names PolyzyMD writes.

    ``prod.tpr``, ``em.mdp``, each ``eq_NN_<stage>.mdp``, ``prod.mdp``, the
    topology of :func:`topology_name` and every file it includes with
    ``#include "..."`` that resolves in ``working_dir``. Other files in the
    folder, such as a backup ``.top``, are never returned.
    """
    folder = Path(working_dir)
    stages = getattr(getattr(config, "simulation_phases", None), "equilibration_stages", None)
    names = ["prod.tpr", "em.mdp"]
    names += [f"eq_{i:02d}_{stage.name}.mdp" for i, stage in enumerate(stages or [], start=1)]
    names.append("prod.mdp")
    files = [folder / name for name in names if (folder / name).is_file()]
    pending = [folder / topology_name(config)]
    while pending:
        path = pending.pop(0)
        if path in files or not path.is_file():
            continue
        files.append(path)
        for line in path.read_text(errors="replace").splitlines():
            include = _INCLUDE.match(line)
            if include:
                pending.append(path.parent / include.group(1))
    return files


def gromacs_topology_file(directory: str | Path, name: str | None = None) -> Path:
    """Return the GROMACS topology (``.top``) in ``directory``.

    PolyzyMD writes one, ``<system prefix>.top``, beside ``prod.tpr``. With
    ``name`` (from :func:`topology_name`), that file is used when it exists,
    so another ``.top`` in the folder does not matter. Without it, or when it
    is missing, the folder must hold exactly one ``.top``.

    Raises
    ------
    ProtocolError
        When ``directory`` has no ``.top`` file or more than one.
    """
    from polyzymd.analyses.exceptions import ProtocolError

    if name and (Path(directory) / name).is_file():
        return Path(directory) / name
    candidates = sorted(Path(directory).glob("*.top"))
    if len(candidates) == 1:
        return candidates[0]
    found = ", ".join(path.name for path in candidates) or "none"
    raise ProtocolError(
        f"This MDAnalysis cannot read {Path(directory) / 'prod.tpr'}, and PolyzyMD reads "
        f"the run's GROMACS topology instead, but {Path(directory)} needs exactly one "
        f".top file (found: {found}).",
        hint="Keep the <prefix>.top and its .itp files that PolyzyMD wrote beside prod.tpr, "
        "set gromacs.analysis_topology in the config to the .top to use, "
        "or compile prod.tpr with a GROMACS version this MDAnalysis reads.",
    )


def atomic_numbers(top_file: str | Path) -> dict[str, int]:
    """Return the atomic number of each atom type in the ``[ atomtypes ]`` of ``top_file``.

    The ``[ atomtypes ]`` sections are read from ``top_file`` and from every
    file it includes with ``#include "..."``, resolved against the including
    file's folder; an include that does not resolve there, such as a
    force-field folder in the GROMACS installation, is skipped. A line reads
    ``name [bonded type] [atomic number] mass charge ptype V W``; a type whose
    line has no atomic number is left out.
    """
    numbers: dict[str, int] = {}
    seen: set[Path] = set()

    def read(path: Path) -> None:
        path = path.resolve()
        if path in seen or not path.is_file():
            return
        seen.add(path)
        section = None
        for raw in path.read_text().splitlines():
            line = raw.split(";", 1)[0].strip()
            if not line:
                continue
            include = _INCLUDE.match(line)
            if include:
                read(path.parent / include.group(1))
                continue
            header = _SECTION.match(line)
            if header:
                section = header.group(1).lower()
                continue
            if section != "atomtypes" or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 6 or fields[-3] not in _PTYPES:
                continue
            before = fields[:-3]
            if len(before) == 5:
                number = before[2]
            elif len(before) == 4 and re.fullmatch(r"\d+", before[1]):
                number = before[1]
            else:
                continue
            numbers[before[0]] = int(number)

    read(Path(top_file))
    return numbers


def universe_from_gromacs_top(
    top_file: str | Path, trajectories: str | Path | Sequence[str | Path]
) -> Universe:
    """Return the universe of ``trajectories`` with the topology of ``top_file``, laid out as a TPR.

    MDAnalysis's ``ITPParser`` reads the atoms, charges, masses, bonds
    (with constraints, and the two O-H bonds of each SETTLE water), angles
    and dihedrals of ``top_file``. They are then arranged as ``TPRParser``
    arranges a TPR: residues numbered 1, 2, ... through the whole system;
    one segment ``seg_<i>_<moltype>`` per block of consecutive molecules of
    one type, as grompp makes from the ``[ molecules ]`` lines; chain IDs
    and ``moltypes`` named after the molecule type; atom ``ids`` and
    ``molnums`` counted from 0; charges and masses in single precision; and
    elements from the atomic numbers in ``[ atomtypes ]``
    (:func:`atomic_numbers`).
    """
    import MDAnalysis as mda
    from MDAnalysis.core import topologyattrs as attrs
    from MDAnalysis.core.topology import Topology
    from MDAnalysis.guesser.tables import Z2SYMB

    with warnings.catch_warnings():
        # ITPParser warns that it guesses no elements; they come from the atomic numbers.
        warnings.simplefilter("ignore", UserWarning)
        warnings.simplefilter("ignore", DeprecationWarning)
        itp = mda.Universe(str(top_file), topology_format="ITP", infer_system=True)
    atoms = itp.atoms
    n_atoms = len(atoms)

    residue_of_atom = atoms.resindices
    residues = itp.residues
    molnums = np.asarray(residues.molnums) - int(np.min(residues.molnums))
    moltypes = np.asarray(residues.moltypes, dtype=object)
    # A block starts at every residue whose molecule type differs from that of the
    # molecule before it.
    starts_molecule = np.r_[True, molnums[1:] != molnums[:-1]]
    new_block = starts_molecule & np.r_[True, moltypes[1:] != moltypes[:-1]]
    residue_block = np.cumsum(new_block) - 1
    block_moltypes = moltypes[new_block]
    segids = np.asarray(
        [f"seg_{i}_{moltype}" for i, moltype in enumerate(block_moltypes)], dtype=object
    )

    numbers = atomic_numbers(top_file)
    elements = np.asarray([Z2SYMB.get(numbers.get(t, 0), "") for t in atoms.types], dtype=object)
    atom_moltypes = moltypes[residue_of_atom]
    chain_ids = np.asarray(
        [m[14:] if m.startswith("Protein_chain_") else m for m in atom_moltypes], dtype=object
    )

    topology = Topology(
        n_atoms,
        len(residues),
        len(segids),
        attrs=[
            attrs.Atomids(np.arange(n_atoms, dtype=np.int32)),
            attrs.Atomnames(np.asarray(atoms.names, dtype=object)),
            attrs.Atomtypes(np.asarray(atoms.types, dtype=object)),
            attrs.Charges(np.asarray(atoms.charges, dtype=np.float32)),
            attrs.Masses(np.asarray(atoms.masses, dtype=np.float32)),
            attrs.Resids(np.arange(1, len(residues) + 1, dtype=np.int32)),
            attrs.Resnums(np.arange(1, len(residues) + 1, dtype=np.int32)),
            attrs.Resnames(np.asarray(residues.resnames, dtype=object)),
            attrs.Moltypes(moltypes),
            attrs.Molnums(molnums.astype(np.int32)),
            attrs.Segids(segids),
            attrs.ChainIDs(chain_ids),
            attrs.Bonds([tuple(b) for b in itp.bonds.indices]),
            attrs.Angles([tuple(a) for a in itp.angles.indices]),
            attrs.Dihedrals([tuple(d) for d in itp.dihedrals.indices]),
            attrs.Impropers([tuple(i) for i in itp.impropers.indices]),
        ],
        atom_resindex=np.asarray(residue_of_atom),
        residue_segindex=residue_block,
    )
    if any(elements):
        topology.add_TopologyAttr(attrs.Elements(elements))
    if isinstance(trajectories, (str, Path)):
        coordinates: Any = str(trajectories)
    else:
        coordinates = [str(path) for path in trajectories]
        if len(coordinates) == 1:
            coordinates = coordinates[0]
    return mda.Universe(topology, coordinates)


#: The build's system PDB, which carries PolyzyMD's chain IDs.
BUILD_PDB = "solvated_system.pdb"


def build_pdb_file(topology_file: str | Path) -> Path | None:
    """Return the ``solvated_system.pdb`` of the run whose topology is ``topology_file``.

    PolyzyMD writes it in the replicate directory, the parent of the
    ``gromacs`` folder; the ``gromacs`` folder itself is searched first.
    """
    folder = Path(topology_file).parent
    for candidate in (folder / BUILD_PDB, folder.parent / BUILD_PDB):
        if candidate.is_file():
            return candidate
    return None


def apply_build_chain_ids(universe: Any, pdb_file: str | Path | None) -> dict[str, Any]:
    """Give a universe the chain IDs PolyzyMD assigned when it built the system.

    A TPR or ``.top`` names chains after molecule types (``MOL0``, ...), and
    an OpenMM ``system.prmtop`` has no chain IDs, while PolyzyMD's build puts
    the protein on chain A, the substrate on B and the polymers on C. The chain
    IDs of ``pdb_file`` replace those of ``universe`` when it has the same
    number of atoms with the same residue names in order (water and ion atom
    names differ between the two files). Otherwise, or when ``pdb_file``
    cannot be read, the universe is left unchanged.

    Only the residue name and chain ID columns of the ``ATOM`` and ``HETATM``
    records of the first model are read. Serial numbers and ``CONECT``
    records are not, so a PDB from OpenMM above 99,999 atoms, with hex
    serials, is read like any other.

    Returns
    -------
    dict[str, Any]
        ``applied``, the ``source`` path and, when nothing was applied, the
        ``reason``; also stored as ``universe._polyzymd_chain_ids``.
    """
    metadata: dict[str, Any]
    if pdb_file is None:
        metadata = {"applied": False, "source": None, "reason": f"no {BUILD_PDB}"}
    else:
        resnames: list[str] = []
        chains: list[str] = []
        try:
            with open(pdb_file) as handle:
                for line in handle:
                    if line.startswith(("ATOM  ", "HETATM")):
                        resnames.append(line[17:21].strip())
                        chains.append(line[21:22].strip())
                    elif line.startswith("ENDMDL") or line[:6].strip() == "END":
                        break
        except (OSError, UnicodeDecodeError) as error:
            reason = f"cannot read {BUILD_PDB}: {error}"
        else:
            if len(chains) != len(universe.atoms):
                reason = f"{len(chains)} atoms in {BUILD_PDB} for {len(universe.atoms)}"
            elif not np.array_equal(np.asarray(resnames), universe.atoms.resnames.astype(str)):
                reason = f"residue names of {BUILD_PDB} differ from the topology's"
            elif not any(chains):
                reason = f"{BUILD_PDB} has no chain IDs"
            else:
                reason = None
        if reason is None:
            chain_ids = np.asarray(chains, dtype=object)
            if hasattr(universe.atoms, "chainIDs"):
                universe.atoms.chainIDs = chain_ids
            else:
                universe.add_TopologyAttr("chainIDs", chain_ids)
            metadata = {"applied": True, "source": str(pdb_file)}
        else:
            metadata = {"applied": False, "source": str(pdb_file), "reason": reason}
    universe._polyzymd_chain_ids = metadata
    return metadata
