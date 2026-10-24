"""Build dynamic polymer chains with mBuild from star-sited SMILES.

A monomer is given as its repeat unit, with ``*`` on the two atoms that bond to
the neighbouring repeat units (``[*]CC([*])(C)C(=O)OC`` for methyl
methacrylate). The first ``*`` points to the head of the chain and the second
to its tail. A cap (end group) has one ``*`` (``[*][H]``, ``[*]C``, or a whole
monomer that kept its double bond, ``[*]C=C(C)C(=O)OC``). Each end draws its
cap from its own options, so a chain's sequence is its head cap's label, the
labels of its repeat units, and its tail cap's label (``HAABH``).

For each sequence the generator:

1. assembles the chain's chemistry (atoms, bond orders, formal charges and one
   residue per repeat unit) from those SMILES with RDKit;
2. builds the chain's coordinates with mBuild's ``Polymer`` (the
   ``biopolymers-v0.1.0`` tag of github.com/joelaforet/mbuild) from the same
   repeat units and caps, and maps them onto the chemistry by a graph match;
3. assigns partial charges (NAGL, AM1-BCC or Espaloma), relaxes the chain in
   vacuum with its force field, and checks it for rings pierced by a bond and
   for atoms closer than 1.0 Angstrom, building again with new monomer
   conformers when a check fails;
4. writes the charged chain as the per-sequence SDF with a metadata sidecar,
   the same file cached mode reads.
"""

from __future__ import annotations

import json
import logging
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

if TYPE_CHECKING:
    from openff.toolkit import Molecule as OFFMolecule

logger = logging.getLogger(__name__)

DYNAMIC_CACHE_SCHEMA_VERSION = 5

#: The mBuild release dynamic polymers are built with (pixi.toml pins it).
MBUILD_FORK_URL = "https://github.com/joelaforet/mbuild"
MBUILD_FORK_TAG = "biopolymers-v0.1.0"

#: The notebook that helps a user write star-sited SMILES.
STAR_SITE_NOTEBOOK = "https://polyzymd.readthedocs.io/en/latest/tutorials/star_sited_monomers.html"

#: Length (nm) of the bond mBuild places between two repeat units, or a repeat
#: unit and a cap, before the chain is relaxed.
_LINK_SEPARATION_NM = 0.15

#: Two atoms that are neither bonded nor share a bonded neighbour must be
#: farther apart than this after relaxation (half the Packmol tolerance).
MIN_CONTACT_ANGSTROM = 1.0


class PolymerGenerationError(Exception):
    """Raised when a polymer chain cannot be built."""


# ---------------------------------------------------------------------------
# Star-sited SMILES
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StarSitedUnit:
    """A repeat unit or cap read from star-sited SMILES.

    Attributes
    ----------
    smiles : str
        The SMILES as given.
    mol : rdkit.Chem.Mol
        The unit with explicit hydrogens, each ``*`` replaced by a hydrogen
        (the site hydrogen) that the bond to the next unit replaces.
    anchors : tuple[int, ...]
        Atom index of the atom that carried each ``*``, in SMILES order.
    site_hydrogens : tuple[int, ...]
        Atom index of the site hydrogen of each ``*``, in the same order.
    """

    smiles: str
    mol: Any
    anchors: tuple[int, ...]
    site_hydrogens: tuple[int, ...]


def parse_star_smiles(smiles: str, n_sites: int, what: str) -> StarSitedUnit:
    """Read star-sited SMILES with exactly ``n_sites`` ``*`` atoms.

    Parameters
    ----------
    smiles : str
        SMILES with ``*`` (or ``[*]``) on each linking atom.
    n_sites : int
        2 for a repeat unit, 1 for a cap.
    what : str
        How error messages name the unit (``monomer 'A'``).

    Raises
    ------
    ValueError
        If RDKit cannot read the SMILES, the number of ``*`` is wrong, a ``*``
        is not singly bonded to exactly one real atom, or both ``*`` of a
        repeat unit sit on one atom.
    """
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    try:
        mol = Chem.MolFromSmiles(smiles)
    finally:
        RDLogger.EnableLog("rdApp.*")
    if mol is None:
        raise ValueError(
            f"{what}: RDKit cannot read the SMILES {smiles!r}.\nfix: correct the SMILES."
        )
    stars = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
    kind = "repeat unit" if n_sites == 2 else "cap"
    if len(stars) != n_sites:
        example = "[*]CC([*])(C)C(=O)OC" if n_sites == 2 else "[*][H] or [*]C"
        raise ValueError(
            f"{what}: the SMILES {smiles!r} has {len(stars)} '*' sites; a {kind} needs "
            f"exactly {n_sites}.\n"
            f"fix: put '*' on the atom(s) that bond to the neighbouring unit, for example "
            f"{example}. The notebook at {STAR_SITE_NOTEBOOK} helps you find them."
        )
    anchors = []
    for star in stars:
        if star.GetFormalCharge() or star.GetNumExplicitHs():
            raise ValueError(
                f"{what}: a '*' in {smiles!r} carries a charge or hydrogens; the '*' is only "
                "the bond to the next unit.\nfix: write each site as [*] or *, and put any "
                "charge on the atom that has it."
            )
        bonds = star.GetBonds()
        if len(bonds) != 1 or bonds[0].GetBondType() != Chem.BondType.SINGLE:
            raise ValueError(
                f"{what}: each '*' in {smiles!r} must be bonded to one atom by a single bond.\n"
                "fix: write the repeat unit as it is in the chain (the C=C of a vinyl monomer "
                "becomes C-C), with '*' attached to each linking atom."
            )
        neighbour = star.GetNeighbors()[0]
        if neighbour.GetAtomicNum() == 0:
            raise ValueError(
                f"{what}: two '*' in {smiles!r} are bonded to each other.\n"
                "fix: correct the SMILES."
            )
        anchors.append(neighbour.GetIdx())
    if n_sites == 2 and anchors[0] == anchors[1]:
        raise ValueError(
            f"{what}: both '*' of {smiles!r} are on one atom.\n"
            "fix: give a repeat unit of at least two backbone atoms, for example [*]CC[*] for "
            "polyethylene."
        )

    editable = Chem.RWMol(mol)
    for star in stars:
        hydrogen = editable.GetAtomWithIdx(star.GetIdx())
        hydrogen.SetAtomicNum(1)
        hydrogen.SetIsotope(0)
        hydrogen.SetAtomMapNum(0)
        hydrogen.SetNoImplicit(True)
    site_hydrogens = tuple(star.GetIdx() for star in stars)
    try:
        Chem.SanitizeMol(editable)
    except Exception as exc:  # RDKit raises several sanitization exception types
        raise ValueError(
            f"{what}: the SMILES {smiles!r} is not a valid molecule once each '*' is a "
            f"hydrogen ({exc}).\nfix: check valences and charges."
        ) from exc
    with_h = Chem.AddHs(editable.GetMol())
    return StarSitedUnit(smiles, with_h, tuple(anchors), site_hydrogens)


# ---------------------------------------------------------------------------
# Chemistry of one chain
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MonomerUnit:
    """A configured monomer: its label, name, residue name and repeat unit."""

    label: str
    name: str
    residue_name: str
    unit: StarSitedUnit


@dataclass(frozen=True)
class CapUnit:
    """A configured cap: its label, optional name and residue name, and its unit.

    A cap with a ``residue_name`` is a residue of its own; without one, its
    atoms join the residue of the repeat unit it caps.
    """

    label: str
    name: str | None
    residue_name: str | None
    unit: StarSitedUnit


def assemble_chain(
    sequence: str,
    monomers: dict[str, MonomerUnit],
    head: StarSitedUnit,
    tail: StarSitedUnit,
    head_residue: bool = False,
    tail_residue: bool = False,
) -> tuple[Any, list[int]]:
    """Return the RDKit molecule of a capped chain and each atom's residue number.

    ``sequence`` holds the labels of the repeat units only. Atoms are ordered
    residue by residue. Each repeat unit is one residue. The head cap is
    residue 1 when ``head_residue`` is true, otherwise its atoms join the first
    repeat unit's residue; likewise the tail cap and the last repeat unit.
    Bond orders and formal charges are those of the SMILES; each bond between
    units replaces the two site hydrogens.
    """
    from rdkit import Chem

    chain = Chem.RWMol()
    residues: list[int] = []

    def add_unit(unit: StarSitedUnit, residue: int) -> dict[int, int]:
        drop = set(unit.site_hydrogens)
        index = {}
        for atom in unit.mol.GetAtoms():
            if atom.GetIdx() in drop:
                continue
            new = Chem.Atom(atom.GetAtomicNum())
            new.SetFormalCharge(atom.GetFormalCharge())
            new.SetNoImplicit(True)
            index[atom.GetIdx()] = chain.AddAtom(new)
            residues.append(residue)
        for bond in unit.mol.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a in index and b in index:
                chain.AddBond(index[a], index[b], bond.GetBondType())
        return index

    offset = 1 if head_residue else 0
    head_index = add_unit(head, 1)
    previous_tail = head_index[head.anchors[0]]
    for position, label in enumerate(sequence, start=1):
        unit = monomers[label].unit
        index = add_unit(unit, position + offset)
        chain.AddBond(previous_tail, index[unit.anchors[0]], Chem.BondType.SINGLE)
        previous_tail = index[unit.anchors[1]]
    tail_index = add_unit(tail, len(sequence) + offset + (1 if tail_residue else 0))
    chain.AddBond(previous_tail, tail_index[tail.anchors[0]], Chem.BondType.SINGLE)

    mol = chain.GetMol()
    Chem.SanitizeMol(mol)
    return mol, residues


# ---------------------------------------------------------------------------
# Coordinates from mBuild
# ---------------------------------------------------------------------------


def mbuild_provenance() -> dict[str, Any]:
    """Return the version and source of the installed mBuild.

    The commit comes from pip's ``direct_url.json`` when mBuild was installed
    from git (as pixi.toml does); it is None for other installs.
    """
    import importlib.metadata as metadata

    import mbuild

    record: dict[str, Any] = {
        "package": "mbuild",
        "version": mbuild.__version__,
        "url": None,
        "commit": None,
        "requested_revision": None,
    }
    try:
        direct = metadata.distribution("mbuild").read_text("direct_url.json")
    except metadata.PackageNotFoundError:
        direct = None
    if direct:
        payload = json.loads(direct)
        vcs = payload.get("vcs_info", {})
        record.update(
            url=payload.get("url"),
            commit=vcs.get("commit_id"),
            requested_revision=vcs.get("requested_revision"),
        )
    return record


def mbuild_summary() -> str:
    """Return one line naming the installed mBuild, its source and commit."""
    return builder_summary(mbuild_provenance())


def builder_summary(record: dict[str, Any]) -> str:
    """Return one line naming the mBuild of a builder record, its source and commit."""
    source = ", ".join(
        str(record[key]) for key in ("url", "requested_revision", "commit") if record.get(key)
    )
    return f"mbuild {record.get('version')}" + (f" ({source})" if source else "")


#: The reference of mBuild, cited when a build generated polymer chains with it.
MBUILD = (
    "mBuild: Klein, C., Sallai, J., Jones, T. J., Iacovella, C. R., McCabe, C. and Cummings, P. T. "
    "(2016). A hierarchical, component based approach to screening properties of soft "
    "matter. In Foundations of Molecular Modeling and Simulation (Molecular Modeling and "
    "Simulation), 79-92. Springer Singapore. doi:10.1007/978-981-10-1128-3_5"
)


def require_mbuild_polymer() -> None:
    """Raise ImportError when the installed mBuild cannot build these polymers."""
    try:
        from mbuild.polymer import Polymer
    except ImportError as exc:
        raise ImportError(
            "dynamic polymers need mBuild from "
            f"{MBUILD_FORK_URL} (tag {MBUILD_FORK_TAG}), which the PolyzyMD pixi "
            "environments install.\nfix: run polyzymd in the pixi build environment, or "
            f"pip install 'mbuild @ git+{MBUILD_FORK_URL}@{MBUILD_FORK_TAG}'."
        ) from exc
    import inspect

    if "bond_tag" not in inspect.signature(Polymer.add_end_groups).parameters:
        raise ImportError(
            "the installed mBuild has no star-site (tagged) end groups; dynamic polymers need "
            f"mBuild from {MBUILD_FORK_URL} (tag {MBUILD_FORK_TAG}).\nfix: pip install "
            f"'mbuild @ git+{MBUILD_FORK_URL}@{MBUILD_FORK_TAG}', or use the pixi environment."
        )


def _tagged_compound(unit: StarSitedUnit, tags: Sequence[str], seed: int) -> Any:
    """Return an mBuild compound of ``unit`` with its anchors tagged ``tags``.

    The site hydrogens stay in the compound: ``Polymer.add_monomer`` and
    ``add_end_groups`` remove one hydrogen of each tagged atom and put a port
    in its place.
    """
    from mbuild.conversion import from_rdkit

    compound = from_rdkit(unit.mol, smiles_seed=seed)
    particles = list(compound.particles())
    for anchor, tag in zip(unit.anchors, tags):
        particles[anchor].particle_tag = tag
    return compound


def _mbuild_chain(
    sequence: str,
    monomers: dict[str, MonomerUnit],
    head: StarSitedUnit,
    tail: StarSitedUnit,
    seed: int,
) -> tuple[list[int], list[tuple[int, int]], Any]:
    """Build the chain with mBuild's ``Polymer``.

    Returns the atomic number of each particle, its bonds as index pairs, and
    its coordinates in Angstrom. ``Polymer.build`` maps the sorted unique
    characters of the sequence onto the monomers in the order they were
    added, so only the labels the sequence uses are added, sorted.
    """
    import numpy as np
    from mbuild.polymer import Polymer

    polymer = Polymer()
    for label in sorted(set(sequence)):
        compound = _tagged_compound(monomers[label].unit, ("head", "tail"), seed)
        polymer.add_monomer(
            compound, head_tag="head", tail_tag="tail", separation=_LINK_SEPARATION_NM
        )
    for end, cap in (("head", head), ("tail", tail)):
        polymer.add_end_groups(
            _tagged_compound(cap, ("cap",), seed),
            bond_tag="cap",
            separation=_LINK_SEPARATION_NM,
            label=end,
            duplicate=False,
        )
    polymer.build(n=1, sequence=sequence, add_hydrogens=False)

    particles = list(polymer.particles())
    position = {id(p): i for i, p in enumerate(particles)}
    elements = [p.element.atomic_number for p in particles]
    bonds = [(position[id(a)], position[id(b)]) for a, b in polymer.bonds()]
    xyz = np.array([p.pos for p in particles]) * 10.0
    return elements, bonds, xyz


def _map_coordinates(
    chain: Any, elements: list[int], bonds: list[tuple[int, int]], xyz: Any
) -> Any:
    """Return the coordinates of ``chain``'s atoms from mBuild's particles.

    The two graphs are matched on elements and connectivity only, because
    mBuild's ``Polymer.build`` leaves the bonds it makes between units without
    a bond order and its RDKit conversion drops formal charges; the chemistry
    is ``chain``'s. Symmetric atoms (the hydrogens of a methyl) may swap, which
    changes nothing.

    Raises
    ------
    PolymerGenerationError
        If the mBuild chain is not the assembled chain.
    """
    import numpy as np
    from rdkit import Chem

    def skeleton(atomic_numbers, pairs):
        editable = Chem.RWMol()
        for number in atomic_numbers:
            atom = Chem.Atom(int(number))
            atom.SetNoImplicit(True)
            editable.AddAtom(atom)
        for a, b in pairs:
            editable.AddBond(int(a), int(b), Chem.BondType.SINGLE)
        mol = editable.GetMol()
        mol.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(mol)
        return mol

    query = skeleton(
        [atom.GetAtomicNum() for atom in chain.GetAtoms()],
        [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in chain.GetBonds()],
    )
    target = skeleton(elements, bonds)
    if query.GetNumAtoms() != target.GetNumAtoms() or query.GetNumBonds() != target.GetNumBonds():
        raise PolymerGenerationError(
            f"mBuild built {target.GetNumAtoms()} atoms and {target.GetNumBonds()} bonds; "
            f"the chain has {query.GetNumAtoms()} atoms and {query.GetNumBonds()} bonds."
        )
    match = target.GetSubstructMatch(query)
    if len(match) != query.GetNumAtoms():
        raise PolymerGenerationError("the chain mBuild built does not match its SMILES.")
    return np.asarray(xyz)[list(match)]


# ---------------------------------------------------------------------------
# Geometry checks
# ---------------------------------------------------------------------------


def find_pierced_rings(mol: Any, xyz: Any) -> list[tuple[tuple[int, ...], tuple[int, int]]]:
    """Return each (ring, bond) where the bond passes through the ring.

    A bond pierces a ring when its two atoms lie on opposite sides of the
    ring's best-fit plane and the point where it crosses the plane lies inside
    the ring polygon. Bonds that share an atom with the ring are skipped.
    """
    import numpy as np

    xyz = np.asarray(xyz, dtype=float)
    found = []
    for ring in mol.GetRingInfo().AtomRings():
        points = xyz[list(ring)]
        centre = points.mean(axis=0)
        _, _, vt = np.linalg.svd(points - centre)
        u, v, normal = vt[0], vt[1], vt[2]
        polygon = np.column_stack(((points - centre) @ u, (points - centre) @ v))
        ring_atoms = set(ring)
        for bond in mol.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a in ring_atoms or b in ring_atoms:
                continue
            da = (xyz[a] - centre) @ normal
            db = (xyz[b] - centre) @ normal
            if da * db >= 0:
                continue
            cross = xyz[a] + da / (da - db) * (xyz[b] - xyz[a]) - centre
            if _inside(polygon, (cross @ u, cross @ v)):
                found.append((tuple(ring), (a, b)))
    return found


def _inside(polygon: Any, point: tuple[float, float]) -> bool:
    """Return whether a 2D point lies inside a polygon (ray casting)."""
    x, y = point
    inside = False
    n = len(polygon)
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
            inside = not inside
    return inside


def closest_contact(mol: Any, xyz: Any) -> tuple[float, tuple[int, int]] | None:
    """Return the shortest distance (Angstrom) between atoms more than two bonds apart."""
    import numpy as np
    from rdkit import Chem

    xyz = np.asarray(xyz, dtype=float)
    path = Chem.GetDistanceMatrix(mol)
    diff = xyz[:, None, :] - xyz[None, :, :]
    distance = np.sqrt((diff**2).sum(-1))
    distance[path <= 2] = np.inf
    i, j = np.unravel_index(np.argmin(distance), distance.shape)
    if not np.isfinite(distance[i, j]):
        return None
    return float(distance[i, j]), (int(i), int(j))


# ---------------------------------------------------------------------------
# Relaxation
# ---------------------------------------------------------------------------


def relax_in_vacuum(molecule: "OFFMolecule", force_field: Sequence[str]) -> None:
    """Minimize ``molecule``'s conformer in vacuum with ``force_field`` and its charges.

    The chain mBuild builds has bonds of the port length and side chains placed
    without regard to each other; a short minimization with the force field
    the build uses gives the packing and the checks a physical structure.
    Replaces conformer 0.
    """
    import numpy as np
    import openmm
    from openff.interchange import Interchange
    from openff.toolkit import ForceField
    from openff.units import unit
    from openmm import unit as omm_unit

    interchange = Interchange.from_smirnoff(
        ForceField(*force_field),
        molecule.to_topology(),
        charge_from_molecules=[molecule],
    )
    system = interchange.to_openmm_system()
    integrator = openmm.VerletIntegrator(0.001)
    platform = openmm.Platform.getPlatformByName("CPU")
    context = openmm.Context(system, integrator, platform)
    context.setPositions(molecule.conformers[0].m_as(unit.nanometer))
    openmm.LocalEnergyMinimizer.minimize(context, 10.0, 2000)
    positions = context.getState(getPositions=True).getPositions(asNumpy=True)
    xyz = np.asarray(positions.value_in_unit(omm_unit.angstrom))
    molecule._conformers = [unit.Quantity(xyz, unit.angstrom)]


# ---------------------------------------------------------------------------
# SDF writing
# ---------------------------------------------------------------------------


def molecule_to_sdf(path: Path, molecule: "OFFMolecule") -> None:
    """Write ``molecule`` with its charges and per-atom metadata to an SDF.

    The metadata goes into the ``metadata`` property, which
    :func:`polyzymd.utils.topology.topology_from_sdf` unpacks onto the atoms.
    """
    from copy import deepcopy

    packaged = deepcopy(molecule)
    packaged.properties["metadata"] = {
        atom.molecule_atom_index: dict(atom.metadata) for atom in packaged.atoms
    }
    packaged.to_file(str(path), file_format="SDF")


# ---------------------------------------------------------------------------
# The generator
# ---------------------------------------------------------------------------


class PolymerGenerator:
    """Generate charged polymer chains with mBuild and cache them as SDF files.

    Parameters
    ----------
    monomers : dict[str, MonomerUnit]
        The configured monomers by sequence label.
    end_groups : dict[str, dict[str, CapUnit]]
        The cap options of the ``head`` and of the ``tail``, by label.
    cache_directory : Path
        Directory of the charged SDFs and their metadata sidecars.
    force_field : sequence of str
        Force field files of the polymer; the chain is relaxed with them.
    max_retries : int
        Builds tried, each with new monomer conformers, before a chain whose
        rings are pierced or whose atoms overlap is refused.
    charger_type : str
        nagl, am1bcc or espaloma.
    type_prefix : str
        Filename prefix of every chain (``{type_prefix}_seq=..._charged.sdf``).
    nagl_model : str | None
        NAGL model for ``charger_type`` nagl.
    """

    def __init__(
        self,
        monomers: dict[str, MonomerUnit],
        end_groups: dict[str, dict[str, CapUnit]],
        cache_directory: Path,
        force_field: Sequence[str],
        max_retries: int = 10,
        charger_type: str = "nagl",
        type_prefix: str = "polymer",
        nagl_model: str | None = None,
    ):
        self.monomers = monomers
        self.end_groups = end_groups
        self.cache_directory = Path(cache_directory)
        self.force_field = [str(name) for name in force_field]
        self.max_retries = max_retries
        self.charger_type = charger_type.lower()
        self.type_prefix = type_prefix
        if nagl_model is None and self.charger_type == "nagl":
            from polyzymd.utils.charging import DEFAULT_NAGL_MODEL

            nagl_model = DEFAULT_NAGL_MODEL
        self.nagl_model = nagl_model if self.charger_type == "nagl" else None

    @classmethod
    def from_smiles(
        cls,
        monomers: Sequence[tuple[str, str, str, str]],
        head: Sequence[dict[str, Any]],
        tail: Sequence[dict[str, Any]],
        **kwargs: Any,
    ) -> "PolymerGenerator":
        """Create a generator from monomers and the cap options of each end.

        ``monomers`` holds (label, name, residue_name, smiles) tuples; ``head``
        and ``tail`` hold cap options as mappings with ``label``, ``smiles``
        and optionally ``name`` and ``residue_name`` (as in the config).
        """
        units = {
            label: MonomerUnit(
                label, name, residue_name, parse_star_smiles(smiles, 2, f"monomer {label!r}")
            )
            for label, name, residue_name, smiles in monomers
        }
        caps = {
            end: {
                option["label"]: CapUnit(
                    option["label"],
                    option.get("name"),
                    option.get("residue_name"),
                    parse_star_smiles(option["smiles"], 1, f"end_groups.{end} {option['label']!r}"),
                )
                for option in options
            }
            for end, options in (("head", head), ("tail", tail))
        }
        return cls(units, caps, **kwargs)

    def split(self, sequence: str) -> tuple[CapUnit, str, CapUnit]:
        """Return the head cap, the repeat-unit labels and the tail cap of ``sequence``.

        Raises
        ------
        PolymerGenerationError
            If the sequence is too short or names a cap or monomer label that
            the config does not have.
        """
        if len(sequence) < 3:
            raise PolymerGenerationError(
                f"sequence {sequence!r} needs a head cap, at least one repeat unit and a tail cap."
            )
        head, middle, tail = sequence[0], sequence[1:-1], sequence[-1]
        problems = []
        if head not in self.end_groups["head"]:
            problems.append(f"head cap {head!r}")
        if tail not in self.end_groups["tail"]:
            problems.append(f"tail cap {tail!r}")
        unknown = sorted(set(middle) - set(self.monomers))
        if unknown:
            problems.append(f"monomers {unknown}")
        if problems:
            raise PolymerGenerationError(
                f"sequence {sequence!r} uses labels the config does not have: "
                + ", ".join(problems)
                + "."
            )
        return self.end_groups["head"][head], middle, self.end_groups["tail"][tail]

    # -- files ---------------------------------------------------------------

    def sdf_path(self, sequence: str) -> Path:
        """Return the charged SDF path of ``sequence`` in the cache directory."""
        name = f"{self.type_prefix}_seq={sequence}_{len(sequence)}-mer_charged.sdf"
        return self.cache_directory / name

    @staticmethod
    def metadata_path(sdf_path: Path) -> Path:
        """Return the metadata sidecar of a charged SDF."""
        return sdf_path.with_name(f"{sdf_path.name}.metadata.json")

    def cache_metadata(self, sequence: str) -> dict[str, Any]:
        """Return what a cached chain of ``sequence`` must have been built from."""
        builder = mbuild_provenance()
        head, middle, tail = self.split(sequence)
        return {
            "schema_version": DYNAMIC_CACHE_SCHEMA_VERSION,
            "sequence": sequence,
            "monomers": {
                label: {
                    "name": monomer.name,
                    "residue_name": monomer.residue_name,
                    "smiles": monomer.unit.smiles,
                }
                for label, monomer in sorted(self.monomers.items())
                if label in middle
            },
            "end_groups": {
                end: {
                    "label": cap.label,
                    "name": cap.name,
                    "residue_name": cap.residue_name,
                    "smiles": cap.unit.smiles,
                }
                for end, cap in (("head", head), ("tail", tail))
            },
            "charger_type": self.charger_type,
            "nagl_model": self.nagl_model,
            "relaxation_force_field": self.force_field,
            "builder": {
                "package": builder["package"],
                "version": builder["version"],
                "commit": builder["commit"],
            },
        }

    def cache_metadata_matches(self, sdf_path: Path, sequence: str) -> bool:
        """Return whether ``sdf_path`` has a sidecar matching this generator's inputs."""
        path = self.metadata_path(sdf_path)
        if not path.exists():
            logger.warning("Ignoring dynamic polymer cache without metadata sidecar: %s", sdf_path)
            return False
        try:
            cached = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            logger.warning("Ignoring dynamic polymer cache with invalid metadata: %s", path)
            return False
        if not isinstance(cached, dict):
            logger.warning("Ignoring dynamic polymer cache with non-object metadata: %s", path)
            return False
        expected = self.cache_metadata(sequence)
        if cached != expected:
            fields = sorted(
                key for key in set(cached) | set(expected) if cached.get(key) != expected.get(key)
            )
            logger.warning(
                "Ignoring dynamic polymer cache with mismatched metadata fields %s: %s",
                fields,
                path,
            )
            return False
        return True

    # -- building ------------------------------------------------------------

    def _charger(self):
        from polyzymd.utils.charging import get_charger

        if self.charger_type == "nagl":
            return get_charger("nagl", model_name=self.nagl_model)
        return get_charger(self.charger_type)

    def build_molecule(self, sequence: str) -> "OFFMolecule":
        """Build, charge, relax and check one chain; return the OpenFF molecule."""
        import numpy as np
        from openff.toolkit import Molecule
        from openff.units import unit

        head_cap, middle, tail_cap = self.split(sequence)
        require_mbuild_polymer()
        head, tail = head_cap.unit, tail_cap.unit
        chain, residues = assemble_chain(
            middle,
            self.monomers,
            head,
            tail,
            head_residue=bool(head_cap.residue_name),
            tail_residue=bool(tail_cap.residue_name),
        )
        names = self._residue_labels(head_cap, middle, tail_cap)

        problems = []
        base_seed = zlib.crc32(sequence.encode("utf-8")) % (2**31 - 1)
        for attempt in range(self.max_retries):
            seed = base_seed + attempt
            state = np.random.get_state()
            np.random.seed(seed)
            try:
                elements, bonds, xyz = _mbuild_chain(middle, self.monomers, head, tail, seed)
            finally:
                np.random.set_state(state)
            coordinates = _map_coordinates(chain, elements, bonds, xyz)

            molecule = Molecule.from_rdkit(
                chain, allow_undefined_stereo=True, hydrogens_are_explicit=True
            )
            molecule._conformers = [unit.Quantity(coordinates, unit.angstrom)]
            self._label_residues(molecule, names, residues)
            molecule = self._charger().charge_molecule(molecule)
            relax_in_vacuum(molecule, self.force_field)

            relaxed = molecule.conformers[0].m_as(unit.angstrom)
            pierced = find_pierced_rings(chain, relaxed)
            contact = closest_contact(chain, relaxed)
            if not pierced and (contact is None or contact[0] >= MIN_CONTACT_ANGSTROM):
                return molecule
            if pierced:
                ring, bond = pierced[0]
                problems.append(f"bond {bond} pierces ring {ring}")
            else:
                problems.append(f"atoms {contact[1]} are {contact[0]:.2f} A apart")
            logger.info("Polymer %s attempt %d refused: %s", sequence, attempt + 1, problems[-1])
        raise PolymerGenerationError(
            f"polymer {sequence!r}: no build out of {self.max_retries} passed the geometry "
            f"checks (last: {problems[-1]}).\nfix: raise max_retries, or check the monomer "
            "SMILES in the notebook."
        )

    def _residue_labels(self, head: CapUnit, middle: str, tail: CapUnit) -> list[tuple[str, str]]:
        """Return the (residue_name, name) of each residue of a chain, in order."""
        names = []
        if head.residue_name:
            names.append((head.residue_name, head.name or head.label))
        names += [
            (self.monomers[label].residue_name, self.monomers[label].name) for label in middle
        ]
        if tail.residue_name:
            names.append((tail.residue_name, tail.name or tail.label))
        return names

    def _label_residues(
        self, molecule: "OFFMolecule", names: list[tuple[str, str]], residues: list[int]
    ) -> None:
        """Give each atom the residue metadata the rest of PolyzyMD reads from SDFs."""
        for atom, residue in zip(molecule.atoms, residues):
            residue_name, name = names[residue - 1]
            atom.metadata.update(
                {
                    "residue_name": residue_name,
                    "residue_number": residue,
                    "insertion_code": " ",
                    "chain_id": "1",
                    "extended_name": name,
                }
            )

    def generate_polymer(self, sequence: str, force_regenerate: bool = False) -> "OFFMolecule":
        """Return the charged chain of ``sequence``, from the cache or built now.

        The returned molecule is read back from the SDF just written, so a
        chain built now and the same chain read from the cache are identical.
        """
        from polyzymd.utils import get_largest_offmol, topology_from_sdf

        sdf_path = self.sdf_path(sequence)
        if (
            sdf_path.exists()
            and not force_regenerate
            and self.cache_metadata_matches(sdf_path, sequence)
        ):
            logger.info("Loading cached polymer from %s", sdf_path)
            return get_largest_offmol(topology_from_sdf(sdf_path))

        logger.info("Building polymer %s with mBuild", sequence)
        try:
            molecule = self.build_molecule(sequence)
        except PolymerGenerationError:
            raise
        except Exception as exc:
            raise PolymerGenerationError(
                f"polymer {sequence!r} could not be built: {type(exc).__name__}: {exc}"
            ) from exc
        self.cache_directory.mkdir(parents=True, exist_ok=True)
        molecule_to_sdf(sdf_path, molecule)
        self.metadata_path(sdf_path).write_text(
            json.dumps(self.cache_metadata(sequence), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        logger.info("Saved charged SDF: %s", sdf_path)
        return get_largest_offmol(topology_from_sdf(sdf_path))
