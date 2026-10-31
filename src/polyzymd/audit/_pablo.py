"""Read a modified protein's chemistry from its PDB file with openff-pablo, for the audit.

The build perceives a modified protein's chemistry (bond orders, formal
charges) through mBuild's RDKit export (``to_openff``). The audit perceives it
again here, from the PDB file and residue definitions, so a perception error
in the build does not pass unseen:

- every residue no modification touched is read with its Chemical Component
  Dictionary definition (``STD_CCD_CACHE``, protonation variants included);
- the residues a modification touched (an anchor, the residues attached to
  it, a non-canonical mutant), with every residue bonded to them other than
  by a peptide bond, are read as one residue: their atoms are written under
  the anchor's chain and number with a four-letter residue name (Pablo reads
  those, and they never clash with the CCD), attached atoms renamed so names
  stay unique, and its definition is built from those mBuild residues' atoms,
  bonds, bond orders and formal charges. (Pablo allows a residue one
  crosslink; a glycan's first sugar bonds both the protein and the next
  sugar, so crosslinks cannot describe it.)

The touched residues' definitions therefore come from mBuild's residue
chemistry, but not through the RDKit export the build uses.

openff-pablo is not released on PyPI or conda-forge (PolyzyMD pins its git
tag), and it is slated to be absorbed into the OpenFF toolkit no earlier than
January 2027. Every Pablo call is in this module, so switching to the
toolkit's reader means changing only this file.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Any

#: Bond orders of mBuild bonds as RDKit bond types (filled on first use).
_BOND_TYPES: dict[float, Any] = {}


def perceive(protein: Any) -> Any:
    """Return the OpenFF Molecule Pablo reads for ``protein``, in ``atom_keys()`` order.

    Raises:
        ValueError: If Pablo does not read one molecule with the protein's atoms,
            or a modification joins two chain residues (not supported yet).
    """
    from openff.pablo import STD_CCD_CACHE, topology_from_pdb

    keys = [tuple(key) for key in protein.atom_keys()]
    residues = _residues(protein)
    residue_of = {
        particle: key for key, residue in residues.items() for particle in residue.particles()
    }
    groups = _groups(protein, residue_of, residues)
    definitions, renamed, written_as = {}, {}, {}
    for k, (anchor, members, links, in_chain) in enumerate(groups):
        name = f"PZ{k:02X}"
        definition, names = _definition(residues, anchor, members, links, name, in_chain)
        definitions[name] = [definition]
        renamed.update(names)
        for residue in members:
            written_as[residue] = (name, anchor)
    library = STD_CCD_CACHE.with_(definitions) if definitions else STD_CCD_CACHE

    with tempfile.TemporaryDirectory(prefix="polyzymd-pablo-") as folder:
        path = Path(folder) / "protein.pdb"
        protein.save_pdb(path)
        _rewrite(path, written_as, renamed)
        topology = topology_from_pdb(str(path), residue_library=library)
    return _one_molecule(topology, keys, {new: old for old, new in renamed.items()})


def _one_molecule(topology: Any, keys: list[tuple], original: dict) -> Any:
    """Every atom and bond Pablo read, as one Molecule in ``keys`` order.

    A protein may be several molecules (chains no bond joins, such as the two
    chains of a protomer); each atom is placed by its chain, residue number,
    insertion code and name.
    """
    from openff.toolkit import Molecule

    position = {key: k for k, key in enumerate(keys)}
    placed: dict[int, tuple[Any, Any]] = {}
    where: dict[tuple[int, int], int] = {}
    for m, molecule in enumerate(topology.molecules):
        for index, atom in enumerate(molecule.atoms):
            meta = atom.metadata
            read = (
                meta["chain_id"],
                int(meta["residue_number"]),
                str(meta.get("insertion_code", "")).strip(),
                atom.name,
            )
            key = original.get(read, read)
            if key not in position or position[key] in placed:
                raise ValueError(f"Pablo reads an atom the protein does not have: {read}")
            placed[position[key]] = (atom, key)
            where[(m, index)] = position[key]
    if len(placed) != len(keys):
        raise ValueError(f"Pablo reads {len(placed)} atoms, the protein has {len(keys)}")
    combined = Molecule()
    for k in range(len(keys)):
        atom, key = placed[k]
        new = combined.add_atom(
            atom.atomic_number, atom.formal_charge, atom.is_aromatic, name=key[3]
        )
        combined.atom(new).metadata.update(
            {"chain_id": key[0], "residue_number": str(key[1]), "insertion_code": key[2]}
        )
    for m, molecule in enumerate(topology.molecules):
        for bond in molecule.bonds:
            combined.add_bond(
                where[(m, bond.atom1_index)],
                where[(m, bond.atom2_index)],
                bond.bond_order,
                bond.is_aromatic,
            )
    return combined


def _groups(protein: Any, residue_of: dict, residues: dict) -> list[tuple]:
    """The residues of each modification, read as one residue.

    Returns ``(anchor, members, links, in_chain)`` per modification: the
    chain residue it is anchored on (else its first residue), its residues,
    the bonds between them as ``(particle, particle, bond order)``, and
    whether the anchor is a chain residue.
    """
    touched, attached = set(), set()
    for record in protein.records:
        touched.add((record.site.chain_id, record.site.resnum, record.site.icode))
        for chain_id, resnum, icode, _ in record.new_residues:
            touched.add((chain_id, resnum, icode))
            attached.add((chain_id, resnum, icode))

    def chain_residue(key: tuple) -> bool:
        # An attached amino acid (a CCD component) has N, CA and C too; mBuild
        # numbers it after the chain, so the records, not the atoms, say what it is.
        return key not in attached and _in_chain(residues[key])

    parent = {residue: residue for residue in residues}

    def find(residue):
        while parent[residue] != residue:
            parent[residue] = parent[parent[residue]]
            residue = parent[residue]
        return residue

    links = []
    for p1, p2, data in protein.mbuild.bonds(return_bond_order=True):
        r1, r2 = residue_of[p1], residue_of[p2]
        if r1 == r2 or not ({r1, r2} & touched):
            continue
        if {p1.name, p2.name} == {"C", "N"} and chain_residue(r1) and chain_residue(r2):
            continue  # a peptide bond
        links.append((p1, p2, float(data["bond_order"])))
        parent[find(r1)] = find(r2)
    members: dict[tuple, list[tuple]] = {}
    for residue in touched | {residue_of[p] for p1, p2, _ in links for p in (p1, p2)}:
        members.setdefault(find(residue), []).append(residue)
    groups = []
    for root, group in sorted(members.items()):
        chain_residues = [r for r in group if chain_residue(r)]
        if len(chain_residues) > 1:
            raise ValueError(
                "a modification bonds two chain residues ("
                + ", ".join(f"{c}:{n}{i}" for c, n, i in sorted(chain_residues))
                + "); reading it independently is not supported yet"
            )
        anchor = chain_residues[0] if chain_residues else sorted(group)[0]
        in_group = set(group)
        group_links = [link for link in links if residue_of[link[0]] in in_group]
        groups.append((anchor, sorted(group), group_links, bool(chain_residues)))
    return groups


def _in_chain(residue: Any) -> bool:
    return {"N", "CA", "C"} <= {particle.name for particle in residue.particles()}


def _residues(protein: Any) -> dict[tuple, Any]:
    """Every mBuild residue of ``protein``, by (chain, number, insertion code)."""
    out = {}
    for chain in protein.mbuild.chains:
        for residue in protein.mbuild.residues(chain.chain_id):
            out[(chain.chain_id, residue.resnum, residue.icode)] = residue
    return out


def _touched(protein: Any) -> set[tuple]:
    """Residues a modification record changed or brought."""
    touched, attached = set(), set()
    for record in protein.records:
        touched.add((record.site.chain_id, record.site.resnum, record.site.icode))
        for chain_id, resnum, icode, _ in record.new_residues:
            touched.add((chain_id, resnum, icode))
            attached.add((chain_id, resnum, icode))

    def chain_residue(key: tuple) -> bool:
        # An attached amino acid (a CCD component) has N, CA and C too; mBuild
        # numbers it after the chain, so the records, not the atoms, say what it is.
        return key not in attached and _in_chain(residues[key])

    return touched


def _definition(
    residues: dict, anchor: tuple, members: list, links: list, name: str, in_chain: bool
):
    """A Pablo definition of one modification's residues read as one residue.

    Returns the definition and the new key of every atom written under
    another residue or name, by original key.
    """
    from openff.pablo import ResidueDefinition
    from openff.pablo.chem import PEPTIDE_BOND
    from openff.toolkit import Molecule
    from rdkit import Chem

    if not _BOND_TYPES:
        _BOND_TYPES.update(
            {
                1.0: Chem.BondType.SINGLE,
                2.0: Chem.BondType.DOUBLE,
                3.0: Chem.BondType.TRIPLE,
                1.5: Chem.BondType.AROMATIC,
            }
        )
    editable, index, names, renamed = Chem.RWMol(), {}, [], {}
    taken = {particle.name for particle in residues[anchor].particles()}
    counter = 0
    for residue_key in [anchor] + [r for r in members if r != anchor]:
        residue = residues[residue_key]
        charges = residue.atom_formal_charges
        for particle in residue.particles():
            atom = Chem.Atom(particle.element.symbol)
            atom.SetFormalCharge(charges.get(particle.name, 0))
            atom.SetNoImplicit(True)
            index[particle] = editable.AddAtom(atom)
            atom_name = particle.name
            if residue_key != anchor:
                counter += 1
                atom_name = f"Z{counter:03d}"
                while atom_name in taken:
                    counter += 1
                    atom_name = f"Z{counter:03d}"
                taken.add(atom_name)
                renamed[(*residue_key, particle.name)] = (*anchor, atom_name)
            names.append((atom_name, False))
        for p1, p2, data in residue.bonds(return_bond_order=True):
            editable.AddBond(index[p1], index[p2], _BOND_TYPES[float(data["bond_order"])])
    for p1, p2, order in links:
        editable.AddBond(index[p1], index[p2], _BOND_TYPES[order])
    by_name = {atom_name: k for k, (atom_name, _) in enumerate(names)}

    def add(symbol: str, to: str, atom_name: str) -> str:
        while atom_name in taken:
            atom_name += "X"
        taken.add(atom_name)
        atom = Chem.Atom(symbol)
        atom.SetNoImplicit(True)
        new = editable.AddAtom(atom)
        editable.AddBond(by_name[to], new, Chem.BondType.SINGLE)
        names.append((atom_name, True))
        by_name[atom_name] = new
        return atom_name

    if in_chain:
        # The peptide bonds' leaving atoms, where the atom has room for one: a
        # labelled N-terminal N or an amidated C-terminal C already has three bonds.
        if editable.GetAtomWithIdx(by_name["N"]).GetDegree() < 3:
            add("H", "N", "H2")
        if editable.GetAtomWithIdx(by_name["C"]).GetDegree() < 3:
            add("H", add("O", "C", "OXT"), "HXT")
    rdmol = editable.GetMol()
    Chem.SanitizeMol(rdmol)
    molecule = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True, hydrogens_are_explicit=True)
    for atom, (atom_name, leaving) in zip(molecule.atoms, names):
        atom.name = atom_name
        atom.metadata["leaving_atom"] = leaving
    definition = ResidueDefinition.from_molecule(
        molecule, residue_name=name, linking_bond=PEPTIDE_BOND if in_chain else None
    )
    return definition, renamed


def _rewrite(path: Path, written_as: dict, renamed: dict) -> None:
    """Write each modification's atoms as one residue: its name, number and atom names.

    Pablo reads a residue from consecutive lines, so the attached residues'
    lines move to just after their anchor's. Serial numbers (and with them
    the CONECT records) are kept.
    """
    head, atoms, tail = [], [], []
    for line in path.read_text().splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM", "TER")):
            atoms.append(line)
        elif atoms:
            tail.append(line)
        else:
            head.append(line)
    moved: dict[tuple, list[str]] = {}
    kept = []
    for line in atoms:
        if line.startswith("TER"):
            kept.append(line)
            continue
        residue = (line[21].strip(), int(line[22:26]), line[26].strip())
        if residue in written_as:
            name, anchor = written_as[residue]
            chain_id, resnum, icode = anchor
            atom_name = line[12:16].strip()
            new_name = renamed.get((*residue, atom_name), (None,) * 4)[3] or atom_name
            atom_field = f" {new_name:<3}" if len(new_name) < 4 else new_name
            line = (
                "ATOM  "
                + line[6:12]
                + atom_field
                + line[16]
                + name
                + f"{chain_id or ' '}{resnum:>4}{icode or ' '}"
                + line[27:]
            )
            if residue != anchor:
                moved.setdefault(anchor, []).append(line)
                continue
        kept.append(line)
    out = []
    for k, line in enumerate(kept):
        out.append(line)
        if line.startswith("TER"):
            continue
        residue = (line[21].strip(), int(line[22:26]), line[26].strip())
        following = kept[k + 1] if k + 1 < len(kept) else None
        same = (
            following is not None
            and not following.startswith("TER")
            and (following[21].strip(), int(following[22:26]), following[26].strip()) == residue
        )
        if residue in moved and not same:
            out.extend(moved.pop(residue))
    for lines in moved.values():  # an anchor that is itself an attached residue
        out.extend(lines)
    path.write_text("".join(head + out + tail))
