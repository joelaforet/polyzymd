"""GLYCAM06j parameters for glycan sites.

GLYCAM06j ships as an OpenMM force field (``amber14/GLYCAM_06j-1.xml``), not
as a SMIRNOFF file, and it leaves selected 1-4 pairs unscaled (a ``<Script>``
in the XML). OpenFF Interchange holds one 1-4 scale per collection and drops
NonbondedForce exceptions on import, so GLYCAM terms never pass through a
force-field file in Interchange. They are built natively and copied in:

1. For each glycan site, a small OpenMM topology is cut out of the protein:
   the linking residue renamed to its GLYCAM linkage template (ASN -> NLN,
   SER -> OLS, THR -> OLT), the glycan residues, and the neighbouring peptide
   groups as ACE and NME caps. ``ForceField(protein.ff14SB.xml,
   GLYCAM_06j-1.xml)`` parameterizes it with an explicit template for every
   residue (mannose and glucose epimers share a bonding graph).
2. Every bonded term touching a site atom, and the charge and Lennard-Jones
   parameters of every site atom, are copied into the protein's Interchange
   in place of the placeholder terms the site was typed with. A term touching
   a site atom reaches at most the previous residue's C and the next
   residue's N, which the caps reproduce. The linkage templates' backbone
   charges equal ff14SB's, which is checked.
3. Two things an Interchange cannot hold are stored on the site and written
   at export: the 1-4 scale of every pair touching a site atom, read from
   the native system (Coulomb and Lennard-Jones separately) and kept as a
   per-pair override with the final charges and Lennard-Jones parameters;
   and GLYCAM's atom-type masses of the site atoms (an Interchange takes
   masses from the elements). :func:`apply_to_openmm` writes both into an
   OpenMM System; :func:`apply_to_gromacs` writes the pairs as ``[ pairs ]``
   function type 2 lines (``i j 2 fudgeQQ qi qj sigma epsilon``), because
   function type 1 cannot scale the Coulomb part per pair, and the masses
   into ``[ atoms ]``.
"""

from __future__ import annotations

import functools
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

#: The parameter-set name of GLYCAM06j.
GLYCAM_FORCE_FIELD = "GLYCAM_06j-1.xml"

#: The OpenMM force-field files that hold GLYCAM06j and its protein side.
OPENMM_FILES = ("amber14/protein.ff14SB.xml", "amber14/GLYCAM_06j-1.xml")

#: GLYCAM's template for a residue that carries a glycan.
LINKAGE_TEMPLATES = {"ASN": "NLN", "SER": "OLS", "THR": "OLT"}

#: The force field that types glycan atoms before their GLYCAM terms replace them.
PLACEHOLDER_FORCE_FIELD = ("openff-2.3.0.offxml",)

#: ff14SB's 1-4 scales, used where a pair's scale cannot be read (a zero charge or epsilon).
_DEFAULT_SCALES = (1 / 1.2, 0.5)

_BACKBONE = ("N", "H", "CA", "HA", "C", "O")


@dataclass
class PairOverride:
    """The 1-4 interaction of one pair of atoms, in protein atom indices."""

    i: int
    j: int
    charge_scale: float
    charge_i: float  # e
    charge_j: float  # e
    sigma: float  # nm
    epsilon: float  # kJ/mol, scale included

    @property
    def charge_product(self) -> float:
        return self.charge_scale * self.charge_i * self.charge_j


@dataclass
class GlycamSite:
    """What one glycan site received from GLYCAM06j."""

    name: str
    atom_indices: list[int]
    linkage: str  # e.g. "A:ASN60 -> NLN"
    residue_templates: list[str]
    n_pairs: int = 0
    n_unscaled_pairs: int = 0
    pairs: list[PairOverride] = field(default_factory=list, repr=False)
    masses: dict[int, float] = field(default_factory=dict, repr=False)  # Da

    def to_dict(self) -> dict:
        data = asdict(self)
        data.pop("pairs")
        data.pop("masses")
        return data


@dataclass
class _NativeTerms:
    charges: dict[int, float]
    lennard_jones: dict[int, tuple[float, float]]  # sigma nm, epsilon kJ/mol
    backbone_charges: dict[int, float]
    bonds: list[tuple[int, int, float, float]]
    angles: list[tuple[int, int, int, float, float]]
    torsions: list[tuple[int, int, int, int, int, float, float]]
    scales: dict[tuple[int, int], tuple[float, float]]
    masses: dict[int, float]


#: The protein atom a glycan bonds to, by linkage residue (N- or O-glycosylation).
GLYCAN_LINK_ATOMS = {"ASN": "ND2", "SER": "OG", "THR": "OG1"}


def glycan_attachment_defaults(
    resname: str,
    atom_name: str | None,
    fragment_atom_name: str | None,
    fragment_resnum: int | None,
    pdb_path: str | Path | None,
) -> tuple[str | None, str, int | None]:
    """Fill the atoms a glycan attachment bonds, where not given.

    - ``atom_name``: the linkage atom of the residue, Asn ND2, Ser OG or Thr
      OG1 (None on any other residue).
    - ``fragment_atom_name``: C1, the anomeric carbon of the reducing-end sugar.
    - ``fragment_resnum`` (when the fragment atom is C1): the first residue of
      the glycan file with a C1. GLYCAM-Web and GlyGen files start with the
      ROH aglycone, which has none, so this is the reducing-end sugar.
    """
    if atom_name is None:
        atom_name = GLYCAN_LINK_ATOMS.get(resname.upper())
    if fragment_atom_name is None:
        fragment_atom_name = "C1"
    if fragment_resnum is None and fragment_atom_name == "C1" and pdb_path is not None:
        path = Path(pdb_path)
        if path.is_file():
            for line in path.read_text().splitlines():
                if line.startswith(("ATOM  ", "HETATM")) and line[12:16].strip() == "C1":
                    fragment_resnum = int(line[22:26])
                    break
    return atom_name, fragment_atom_name, fragment_resnum


def is_glycam(force_field) -> bool:
    """Whether a parameter set's force-field list names GLYCAM."""
    return any("glycam" in str(name).lower() for name in force_field)


def glycam_residue_names() -> set[str]:
    """The residue templates GLYCAM06j defines."""
    return set(_glycam_template_bonds())


@functools.cache
def _glycam_template_bonds() -> dict[str, frozenset[frozenset[str]]]:
    """Each GLYCAM06j residue template's bonds, as pairs of atom names."""
    import xml.etree.ElementTree as ET

    from openmm.app import forcefield

    path = Path(forcefield.__file__).parent / "data" / OPENMM_FILES[1]
    return {
        residue.get("name"): frozenset(
            frozenset((bond.get("atomName1"), bond.get("atomName2")))
            for bond in residue.findall("Bond")
        )
        for residue in ET.parse(path).getroot().iter("Residue")
    }


def check_glycan_pdb(path: str | Path) -> None:
    """Refuse a glycan PDB file that cannot bond to a protein as GLYCAM expects."""
    text = Path(path).read_text().splitlines()
    residues = {line[17:20].strip() for line in text if line.startswith(("ATOM", "HETATM"))}
    if not any(line.startswith("CONECT") for line in text):
        raise ValueError(
            f"{path} has no CONECT records. A glycan's bonds are read from CONECT records; "
            "GLYCAM PDB files from GlyGen or GLYCAM-Web often lack them. Add them (for "
            "example by loading the file in PyMOL and saving it with CONECT records) and "
            "check every bond."
        )
    if "ROH" not in residues:
        raise ValueError(
            f"{path} has no ROH residue at the reducing end. The bond to the protein replaces "
            "the anomeric hydroxyl (ROH), so without it the anomeric carbon has nothing to "
            "leave. Use the GLYCAM file of the glycan with its reducing-end ROH."
        )
    missing = _missing_template_bonds(text)
    if missing:
        raise ValueError(
            f"{path} lacks CONECT records for bonds its GLYCAM06j residues have: "
            + ", ".join(missing[:10])
            + (f" (and {len(missing) - 10} more)" if len(missing) > 10 else "")
            + ". Add every bond as a CONECT record."
        )


def _missing_template_bonds(lines: list[str]) -> list[str]:
    """Bonds of each GLYCAM06j residue in the file that no CONECT record gives."""
    atoms = {}
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            residue = (line[17:20].strip(), line[21], line[22:27].strip())
            atoms[int(line[6:11])] = (residue, line[12:16].strip())
    bonded: set[tuple] = set()
    for line in lines:
        if line.startswith("CONECT"):
            serials = [int(line[i : i + 5]) for i in range(6, len(line.rstrip()), 5)]
            first = atoms.get(serials[0])
            for other in serials[1:]:
                second = atoms.get(other)
                if first and second and first[0] == second[0]:
                    bonded.add((first[0], frozenset((first[1], second[1]))))
    templates = _glycam_template_bonds()
    missing = []
    for residue in dict.fromkeys(residue for residue, _ in atoms.values()):
        names = {name for key, name in atoms.values() if key == residue}
        for bond in sorted(templates.get(residue[0], ()), key=sorted):
            if bond <= names and (residue, bond) not in bonded:
                missing.append(f"{residue[0]}{residue[2]} {'-'.join(sorted(bond))}")
    return missing


def glycam_site_terms(
    molecule, site, library: dict, keys, *, check_backbone: bool = True
) -> tuple[_NativeTerms, GlycamSite]:
    """Parameterize one glycan site natively with GLYCAM06j.

    ``molecule`` is the OpenFF Molecule of the whole protein, ``site`` its
    :class:`~polyzymd.modifying_proteins.topology.ModificationSite`,
    ``library`` the ff14SB library charge of every protein atom by key.
    ``check_backbone`` checks that GLYCAM's linkage-residue backbone charges
    are ff14SB's; a Rosemary protein (whose backbone carries NAGL charges)
    turns it off.
    """
    import openmm
    from openmm import app

    site_atoms = set(site.atom_indices)
    linking = site.residues[0]
    by_residue: dict[tuple, list[int]] = {}
    for index, key in enumerate(keys):
        by_residue.setdefault((key.chain_id, key.resnum, key.icode), []).append(index)
    neighbours = [set() for _ in range(molecule.n_atoms)]
    for bond in molecule.bonds:
        neighbours[bond.atom1_index].add(bond.atom2_index)
        neighbours[bond.atom2_index].add(bond.atom1_index)

    linking_atoms = by_residue[linking]
    name_of = {index: keys[index].atom_name for index in range(molecule.n_atoms)}
    parent = molecule.atom(linking_atoms[0]).metadata["residue_name"]
    label = f"{linking[0]}:{parent}{linking[1]}{linking[2]}"
    if parent not in LINKAGE_TEMPLATES:
        raise ValueError(
            f"GLYCAM06j links glycans to ASN, SER or THR; site {site.name} is on {label}."
        )
    atom_by_name = {name_of[index]: index for index in linking_atoms}
    previous_c = [j for j in neighbours[atom_by_name["N"]] if name_of[j] == "C"]
    next_n = [j for j in neighbours[atom_by_name["C"]] if name_of[j] == "N"]
    if not previous_c or not next_n:
        raise ValueError(
            f"Refused: glycan {site.name} on {label}, the first or last residue of its "
            f"chain. GLYCAM06j's linkage residue ({LINKAGE_TEMPLATES[parent]}) is "
            "parameterized with a neighbouring residue on each side, and a chain-terminal "
            "residue lacks one. File an issue at "
            "https://github.com/joelaforet/polyzymd/issues if you meant to do this and "
            "want this capability."
        )
    glycan_residues = sorted(
        {(keys[i].chain_id, keys[i].resnum, keys[i].icode) for i in site_atoms} - {linking},
        key=lambda residue: min(by_residue[residue]),
    )
    known = glycam_residue_names()
    unknown = sorted(
        {
            molecule.atom(by_residue[residue][0]).metadata["residue_name"]
            for residue in glycan_residues
        }
        - known
    )
    if unknown:
        raise ValueError(
            f"Site {site.name} has residues GLYCAM06j does not define: {', '.join(unknown)}. "
            "GLYCAM06j needs glycans named as GLYCAM residues (GLYCAM-Web or GlyGen files)."
        )

    # The small topology: ACE caps the previous residue, NME the next.
    top = app.Topology()
    chain = top.addChain("A")
    to_local: dict[int, object] = {}
    templates = {}
    previous_residue = (
        keys[previous_c[0]].chain_id,
        keys[previous_c[0]].resnum,
        keys[previous_c[0]].icode,
    )
    next_residue = (keys[next_n[0]].chain_id, keys[next_n[0]].resnum, keys[next_n[0]].icode)

    cap_methyls: set[int] = set()

    def add_cap(residue_key, kept_names, cap_name):
        residue = top.addResidue(cap_name, chain)
        templates[residue] = cap_name
        kept = {name_of[i]: i for i in by_residue[residue_key] if name_of[i] in kept_names}
        for name, index in kept.items():
            element = app.Element.getByAtomicNumber(molecule.atom(index).atomic_number)
            to_local[index] = top.addAtom("CH3" if name == "CA" else name, element, residue)
        methyl = to_local[kept["CA"]]
        cap_methyls.add(kept["CA"])
        for n in range(3):
            hydrogen = top.addAtom(f"HH3{n + 1}", app.element.hydrogen, residue)
            top.addBond(methyl, hydrogen)
        if cap_name == "NME" and "H" not in kept:
            hydrogen = top.addAtom("H", app.element.hydrogen, residue)
            top.addBond(to_local[kept["N"]], hydrogen)

    add_cap(previous_residue, {"CA", "C", "O"}, "ACE")
    linkage_template = LINKAGE_TEMPLATES[parent]
    for residue_key, template in [(linking, linkage_template)] + [
        (residue, molecule.atom(by_residue[residue][0]).metadata["residue_name"])
        for residue in glycan_residues
    ]:
        residue = top.addResidue(template, chain)
        templates[residue] = template
        for index in by_residue[residue_key]:
            element = app.Element.getByAtomicNumber(molecule.atom(index).atomic_number)
            to_local[index] = top.addAtom(name_of[index], element, residue)
    add_cap(next_residue, {"N", "H", "CA"}, "NME")
    for bond in molecule.bonds:
        if bond.atom1_index in to_local and bond.atom2_index in to_local:
            top.addBond(to_local[bond.atom1_index], to_local[bond.atom2_index])

    system = app.ForceField(*OPENMM_FILES).createSystem(
        top,
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        rigidWater=False,
        removeCMMotion=False,
        residueTemplates=templates,
    )
    # Cap methyls stand in for the neighbours' CA; no term touching the site
    # may reach them, so they map to nothing.
    global_of = {
        local.index: index for index, local in to_local.items() if index not in cap_methyls
    }
    site_local = {to_local[index].index for index in site_atoms}
    local_neighbours: dict[int, set[int]] = {}
    for a, b in top.bonds():
        local_neighbours.setdefault(a.index, set()).add(b.index)
        local_neighbours.setdefault(b.index, set()).add(a.index)

    def touching(atoms) -> bool:
        return any(atom in site_local for atom in atoms)

    def mapped(atoms) -> tuple[int, ...]:
        try:
            return tuple(global_of[atom] for atom in atoms)
        except KeyError:
            raise ValueError(
                f"A GLYCAM06j term of site {site.name} reaches a cap atom; the site is too "
                "close to the end of its chain."
            ) from None

    terms = _NativeTerms({}, {}, {}, [], [], [], {}, {})
    for local, index in global_of.items():
        if index in site_atoms:
            terms.masses[index] = system.getParticleMass(local)._value
    for force in system.getForces():
        if isinstance(force, openmm.HarmonicBondForce):
            for n in range(force.getNumBonds()):
                i, j, length, k = force.getBondParameters(n)
                if touching((i, j)):
                    terms.bonds.append((*mapped((i, j)), length._value, k._value))
        elif isinstance(force, openmm.HarmonicAngleForce):
            for n in range(force.getNumAngles()):
                i, j, k_, angle, k = force.getAngleParameters(n)
                if touching((i, j, k_)):
                    terms.angles.append((*mapped((i, j, k_)), angle._value, k._value))
        elif isinstance(force, openmm.PeriodicTorsionForce):
            for n in range(force.getNumTorsions()):
                i, j, k_, l_, periodicity, phase, k = force.getTorsionParameters(n)
                if touching((i, j, k_, l_)):
                    terms.torsions.append(
                        (*mapped((i, j, k_, l_)), periodicity, phase._value, k._value)
                    )
        elif isinstance(force, openmm.NonbondedForce):
            for local, index in global_of.items():
                charge, sigma, epsilon = force.getParticleParameters(local)
                if index in site_atoms:
                    terms.charges[index] = charge._value
                    terms.lennard_jones[index] = (sigma._value, epsilon._value)
                elif index in linking_atoms and name_of[index] in _BACKBONE:
                    terms.backbone_charges[index] = charge._value
            for n in range(force.getNumExceptions()):
                i, j, product, sigma, epsilon = force.getExceptionParameters(n)
                if not touching((i, j)) or not _is_one_four(i, j, local_neighbours):
                    continue
                qi, _, ei = force.getParticleParameters(i)
                qj, _, ej = force.getParticleParameters(j)
                charge_scale = (
                    product._value / (qi._value * qj._value)
                    if abs(qi._value * qj._value) > 1e-12
                    else _DEFAULT_SCALES[0]
                )
                lj_scale = (
                    epsilon._value / np.sqrt(ei._value * ej._value)
                    if ei._value * ej._value > 1e-12
                    else _DEFAULT_SCALES[1]
                )
                pair = tuple(sorted(mapped((i, j))))
                terms.scales[pair] = (charge_scale, lj_scale)

    for index, charge in terms.backbone_charges.items():
        if check_backbone and abs(charge - library[keys[index]]) > 1e-4:
            raise ValueError(
                f"GLYCAM06j's {linkage_template} backbone charge of {name_of[index]} "
                f"({charge:+.4f}) differs from ff14SB's ({library[keys[index]]:+.4f})."
            )
    record = GlycamSite(
        name=site.name,
        atom_indices=sorted(site_atoms),
        linkage=f"{label} -> {linkage_template}",
        residue_templates=[linkage_template]
        + [molecule.atom(by_residue[r][0]).metadata["residue_name"] for r in glycan_residues],
    )
    return terms, record


def inject(interchange, terms: _NativeTerms, record: GlycamSite) -> None:
    """Replace the placeholder terms of a site's atoms with GLYCAM's, in place.

    Fills ``record.pairs`` with the 1-4 override of every pair touching the
    site, from the final charges and Lennard-Jones parameters.
    """
    from openff.interchange.components.potentials import Potential
    from openff.interchange.models import (
        AngleKey,
        BondKey,
        PotentialKey,
        ProperTorsionKey,
        TopologyKey,
    )
    from openff.units import unit

    site = set(record.atom_indices)
    tag = f"glycam::{record.name}::"

    def replace(collection_name, new_terms):
        collection = interchange[collection_name]
        removed = [key for key in collection.key_map if set(key.atom_indices) & site]
        before = {tuple(sorted(key.atom_indices)) for key in removed}
        for key in removed:
            del collection.key_map[key]
        for topology_key, parameters in new_terms:
            potential_key = PotentialKey(
                id=tag + "-".join(map(str, topology_key.atom_indices)),
                mult=getattr(topology_key, "mult", None),
                associated_handler=collection_name,
            )
            collection.key_map[topology_key] = potential_key
            collection.potentials[potential_key] = Potential(parameters=parameters)
        used = set(collection.key_map.values())
        collection.potentials = {k: v for k, v in collection.potentials.items() if k in used}
        return before

    nm, kj = unit.nanometer, unit.kilojoule_per_mole
    before = replace(
        "Bonds",
        [
            (BondKey(atom_indices=(i, j)), {"length": r0 * nm, "k": k * kj / nm**2})
            for i, j, r0, k in terms.bonds
        ],
    )
    _same_terms("bonds", before, terms.bonds, 2, record.name)
    before = replace(
        "Angles",
        [
            (
                AngleKey(atom_indices=(i, j, k_)),
                {"angle": a * unit.radian, "k": k * kj / unit.radian**2},
            )
            for i, j, k_, a, k in terms.angles
        ],
    )
    _same_terms("angles", before, terms.angles, 3, record.name)
    # GLYCAM's impropers go in with the propers: both are periodic torsions
    # over the atoms in the order given, and Interchange's GROMACS export
    # drops an improper whose atoms are not in its own central-atom-first
    # order. The placeholder's impropers on the site are removed.
    multiplicity: dict[tuple, int] = {}
    torsions = []
    for i, j, k_, l_, periodicity, phase, k in terms.torsions:
        atoms = (i, j, k_, l_)
        mult = multiplicity.get(atoms, 0)
        multiplicity[atoms] = mult + 1
        torsions.append(
            (
                ProperTorsionKey(atom_indices=atoms, mult=mult),
                {
                    "k": k * kj,
                    "periodicity": periodicity * unit.dimensionless,
                    "phase": phase * unit.radian,
                    "idivf": 1.0 * unit.dimensionless,
                },
            )
        )
    replace("ProperTorsions", torsions)
    if "ImproperTorsions" in interchange.collections:
        replace("ImproperTorsions", [])

    vdw = interchange["vdW"]
    for index in site:
        sigma, epsilon = terms.lennard_jones[index]
        potential_key = PotentialKey(id=f"{tag}vdW-{index}", associated_handler="vdW")
        vdw.key_map[TopologyKey(atom_indices=(index,))] = potential_key
        vdw.potentials[potential_key] = Potential(
            parameters={"sigma": sigma * nm, "epsilon": epsilon * kj}
        )
    used = set(vdw.key_map.values())
    vdw.potentials = {k: v for k, v in vdw.potentials.items() if k in used}

    if "Constraints" in interchange.collections:
        lengths = {frozenset((i, j)): r0 for i, j, r0, _ in terms.bonds}
        constraints = interchange["Constraints"]
        for key in [k for k in constraints.key_map if set(k.atom_indices) & site]:
            potential_key = PotentialKey(
                id=tag + "constraint-" + "-".join(map(str, key.atom_indices)),
                associated_handler="Bonds",
            )
            constraints.key_map[key] = potential_key
            constraints.potentials[potential_key] = Potential(
                parameters={"distance": lengths[frozenset(key.atom_indices)] * nm}
            )
        used = set(constraints.key_map.values())
        constraints.potentials = {k: v for k, v in constraints.potentials.items() if k in used}

    charges = {
        key.atom_indices[0]: charge.m_as(unit.elementary_charge)
        for key, charge in interchange["Electrostatics"].charges.items()
    }
    lennard_jones = {
        key.atom_indices[0]: (
            vdw.potentials[potential_key].parameters["sigma"].m_as(nm),
            vdw.potentials[potential_key].parameters["epsilon"].m_as(kj),
        )
        for key, potential_key in vdw.key_map.items()
    }
    record.masses = dict(terms.masses)
    record.pairs = []
    for (i, j), (charge_scale, lj_scale) in sorted(terms.scales.items()):
        (si, ei), (sj, ej) = lennard_jones[i], lennard_jones[j]
        record.pairs.append(
            PairOverride(
                i=i,
                j=j,
                charge_scale=float(charge_scale),
                charge_i=float(charges[i]),
                charge_j=float(charges[j]),
                sigma=float(0.5 * (si + sj)),
                epsilon=float(lj_scale * np.sqrt(ei * ej)),
            )
        )
    record.n_pairs = len(record.pairs)
    record.n_unscaled_pairs = int(
        sum(
            abs(scale[0] - 1.0) < 1e-6 and abs(scale[1] - 1.0) < 1e-6
            for scale in terms.scales.values()
        )
    )


def apply_to_openmm(system, sites, offset: int = 0) -> int:
    """Write the GLYCAM06j sites' 1-4 pairs and masses into ``system``.

    ``sites`` are :class:`GlycamSite` records; ``offset`` is the index of the
    protein's first atom in ``system``. The pairs become NonbondedForce
    exceptions; every pair must find its exception. Returns the number of
    exceptions changed.
    """
    import openmm

    for site in sites:
        for index, mass in site.masses.items():
            system.setParticleMass(index + offset, mass)
    pairs = [pair for site in sites for pair in site.pairs]
    forces = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)]
    if len(forces) != 1:
        raise ValueError("Expected one NonbondedForce (combine_nonbonded_forces=True).")
    force = forces[0]
    wanted = {(p.i + offset, p.j + offset): p for p in pairs}
    changed = 0
    for n in range(force.getNumExceptions()):
        i, j, *_ = force.getExceptionParameters(n)
        pair = wanted.get((min(i, j), max(i, j)))
        if pair is not None:
            force.setExceptionParameters(n, i, j, pair.charge_product, pair.sigma, pair.epsilon)
            changed += 1
    if changed != len(wanted):
        raise ValueError(
            f"Found {changed} of {len(wanted)} GLYCAM06j 1-4 exceptions in the system."
        )
    return changed


def apply_to_gromacs(path: str | Path, moleculetype: str, sites) -> int:
    """Write the GLYCAM06j sites' 1-4 pairs and masses into a GROMACS file.

    ``path`` is the ``.top`` or ``.itp`` holding ``[ moleculetype ]``
    ``moleculetype``; atom indices are 1-based within that molecule. Each
    site pair's ``[ pairs ]`` line becomes ``i j 2 fudgeQQ qi qj sigma
    epsilon``, and each site atom's ``[ atoms ]`` mass becomes GLYCAM's.
    Preprocessor lines are left alone. Every pair and atom must be found.
    Returns the number of pair lines rewritten.
    """
    path = Path(path)
    wanted = {(p.i + 1, p.j + 1): p for site in sites for p in site.pairs}
    masses = {index + 1: mass for site in sites for index, mass in site.masses.items()}
    out, section, current = [], None, None
    changed = weighed = 0
    for line in path.read_text().splitlines(keepends=True):
        stripped = line.split(";")[0].strip()
        header = re.match(r"\[\s*(\w+)\s*\]", stripped)
        if header:
            section = header.group(1).lower()
        elif not stripped or stripped.startswith("#"):
            pass
        elif section == "moleculetype":
            current = stripped.split()[0]
        elif current == moleculetype and section == "pairs":
            fields = stripped.split()
            a, b = int(fields[0]), int(fields[1])
            pair = wanted.get((min(a, b), max(a, b)))
            if pair is not None:
                line = (
                    f"{a:7d} {b:7d}  2  {pair.charge_scale:.10f} {pair.charge_i:.10f} "
                    f"{pair.charge_j:.10f} {pair.sigma:.10e} {pair.epsilon:.10e}\n"
                )
                changed += 1
        elif current == moleculetype and section == "atoms":
            fields = stripped.split()
            mass = masses.get(int(fields[0]))
            if mass is not None:
                fields[7] = f"{mass:.6f}"
                line = " ".join(fields) + "\n"
                weighed += 1
        out.append(line)
    if changed != len(wanted) or weighed != len(masses):
        raise ValueError(
            f"Found {changed} of {len(wanted)} GLYCAM06j 1-4 pairs and {weighed} of "
            f"{len(masses)} atoms in {path} ({moleculetype})."
        )
    path.write_text("".join(out))
    return changed


def apply_to_gromacs_topology(
    top_path: str | Path, sites, molecule_ranges, first_molecule: int = 0
) -> int:
    """Write the GLYCAM06j sites into the molecule types of a GROMACS topology.

    ``molecule_ranges`` are the ``[start, stop)`` protein atom ranges of the
    protein's molecules, which are entries ``first_molecule``,
    ``first_molecule + 1``, ... of ``[ molecules ]``. Each site's pairs and
    masses go to the molecule type of the molecule holding it, in indices
    local to that molecule. Molecules that share a molecule type must carry
    the same overrides, written once. Returns the number of pair lines written.
    """
    from dataclasses import replace

    top_path = Path(top_path)
    files = [top_path] + [
        top_path.parent / name
        for name in re.findall(r'^\s*#include\s+"([^"]+)"', top_path.read_text(), re.MULTILINE)
        if (top_path.parent / name).is_file()
    ]
    instances: list[str] = []
    section = None
    for line in top_path.read_text().splitlines():
        stripped = line.split(";")[0].strip()
        header = re.match(r"\[\s*(\w+)\s*\]", stripped)
        if header:
            section = header.group(1).lower()
        elif section == "molecules" and stripped and not stripped.startswith("#"):
            name, count = stripped.split()[:2]
            instances += [name] * int(count)
    per_type: dict[str, list] = {}
    for number, (start, stop) in enumerate(molecule_ranges):
        local = []
        for site in sites:
            inside = [i for i in site.atom_indices if start <= i < stop]
            if not inside:
                continue
            if len(inside) != len(site.atom_indices):
                raise ValueError(f"GLYCAM06j site {site.name} spans two molecules.")
            local.append(
                replace(
                    site,
                    atom_indices=[i - start for i in site.atom_indices],
                    pairs=[replace(q, i=q.i - start, j=q.j - start) for q in site.pairs],
                    masses={i - start: m for i, m in site.masses.items()},
                )
            )
        if not local:
            continue
        name = instances[first_molecule + number]
        key = sorted((tuple(s.atom_indices), len(s.pairs)) for s in local)
        if name in per_type:
            if per_type[name][0] != key:
                raise ValueError(
                    f"Molecules sharing GROMACS molecule type {name} carry different GLYCAM06j "
                    "sites; give each its own molecule type."
                )
            continue
        per_type[name] = [key, local]
    written = 0
    for name, (_, local) in per_type.items():
        pattern = re.compile(
            r"^\s*\[\s*moleculetype\s*\]\s*\n(?:\s*;.*\n)*\s*" + re.escape(name) + r"\s",
            re.MULTILINE,
        )
        holder = next((f for f in files if pattern.search(f.read_text())), None)
        if holder is None:
            raise ValueError(f"No [ moleculetype ] {name} in {top_path} or its includes.")
        written += apply_to_gromacs(holder, name, local)
    return written


def _same_terms(kind, before, new_terms, n_atoms, name) -> None:
    after = {tuple(sorted(term[:n_atoms])) for term in new_terms}
    if before != after:
        raise ValueError(
            f"GLYCAM06j and the placeholder force field disagree on which {kind} touch "
            f"site {name}; the site's bonding differs from GLYCAM's templates."
        )


def _is_one_four(i: int, j: int, neighbours: dict[int, set[int]]) -> bool:
    first = neighbours.get(i, set())
    second = {k for n in first for k in neighbours.get(n, set())}
    third = {k for n in second for k in neighbours.get(n, set())}
    return j in third and j not in second and j not in first and j != i
