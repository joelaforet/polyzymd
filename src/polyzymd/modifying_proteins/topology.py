"""Which atoms each modification owns, and which segment each atom is in.

A *site* is the set of atoms that take the modification's force field and
charges instead of the protein's. For an attachment it is the anchor
residue's side chain from CB outward plus every attached atom. When the
anchor atom is a backbone atom (an N-terminal label), that atom and its
hydrogens join instead, and the side chain stays with the protein. For a
non-canonical residue (a mutation to a CCD component other than the 20
canonical ones, or to a SMILES side chain) it is the side chain. A later
modification that bonds to a site's atoms grows that site. The backbone
(N, H, CA, HA, C, O) of a modified residue stays with the protein.

A *segment* groups atoms for selections and PSF segment IDs: protein
residues, including modified and mutated ones, are ``PROT``; the residues
an attachment adds are a segment named after the attachment.
"""

from __future__ import annotations

from dataclasses import dataclass

from polyzymd.modifying_proteins.api import AtomKey, ModifiedProtein
from polyzymd.modifying_proteins.sites import CANONICAL

#: Backbone atoms that keep the protein's parameters, termini included.
BACKBONE = frozenset({"N", "H", "CA", "HA", "HA2", "HA3", "C", "O", "H1", "H2", "H3", "OXT", "HXT"})

#: Hydrogens that move with each backbone heavy atom when it is the anchor.
_BACKBONE_HYDROGENS = {
    "N": ("H", "H1", "H2", "H3"),
    "CA": ("HA", "HA2", "HA3"),
    "OXT": ("HXT",),
}

#: Protonation changes that give an ff14SB residue variant (LYN, ASH, GLH,
#: HID/HIE/HIP, CYM). mBuild moves the Asp and Glu proton at OD2 and OE2 only.
_TEMPLATE_PROTONATION = {
    ("LYS", "NZ"),
    ("ASP", "OD2"),
    ("GLU", "OE2"),
    ("HIS", "ND1"),
    ("HIS", "NE2"),
    ("CYS", "SG"),
}

#: CCD codes of the D forms of the canonical residues. Like ``stereo="D"``,
#: a mutation to one of them keeps ff14SB, whose parameters do not depend on
#: handedness.
D_CANONICAL = frozenset(
    "DAL DAR DSG DAS DCY DGN DGL DHI DIL DLE DLY MED DPN DPR DSN DTH DTR DTY DVA".split()
)

#: Name prefix of a site on a non-canonical residue.
RESIDUE_SITE_PREFIX = "R"

PROTEIN_SEGMENT = "PROT"


@dataclass(frozen=True)
class ModificationSite:
    """The atoms one modification owns in the final structure."""

    name: str
    type: str  # glycan, polymer, generic or noncanonical: the first modification's
    types: frozenset[str]  # every modification's type, when several share the site
    residues: tuple[tuple[str, int, str], ...]  # anchor residues (chain, number, icode)
    input_residues: tuple[str, ...]  # their names in the input PDB, for backbone charges
    atoms: tuple[AtomKey, ...]
    atom_indices: tuple[int, ...]  # in the order save_pdb writes the atoms


def modification_sites(protein: ModifiedProtein) -> list[ModificationSite]:
    """Return every site of ``protein`` against its current atoms."""
    keys = protein.atom_keys()
    present = set(keys)
    residue_atoms: dict[tuple[str, int, str], list[AtomKey]] = {}
    for key in keys:
        residue_atoms.setdefault((key.chain_id, key.resnum, key.icode), []).append(key)

    names: list[str] = []
    types: dict[str, str] = {}
    anchors: dict[str, list[tuple[str, int, str]]] = {}
    # Residues whose side chain belongs to the site. A residue bonded only
    # through a backbone atom keeps its side chain on the protein.
    side_chains: dict[str, set[tuple[str, int, str]]] = {}
    all_types: dict[str, set[str]] = {}
    explicit: dict[str, set[AtomKey]] = {}
    fragment_residues: dict[tuple[str, int, str], str] = {}
    # Residues protonated outside the ff14SB variants; they become sites of
    # their own unless an attachment claims them.
    pending: list[tuple[str, int, str]] = []
    n_residue_sites = 0

    def owner(residue: tuple[str, int, str]) -> str | None:
        if residue in fragment_residues:
            return fragment_residues[residue]
        return next((name for name in names if residue in anchors[name]), None)

    def rename(old: AtomKey, new: AtomKey) -> None:
        for atoms in explicit.values():
            if old in atoms:
                atoms.discard(old)
                atoms.add(new)

    for record in protein.records:
        residue = (record.site.chain_id, record.site.resnum, record.site.icode)
        for old, new in record.renamed:
            rename(old, new)
        site = owner(residue)
        if record.kind == "mutate":
            if site is not None:
                _drop(site, names, types, anchors, explicit)
                side_chains.pop(site), all_types.pop(site)
            canonical = record.residue_after in CANONICAL | D_CANONICAL
            if not canonical or _is_side_chain_fragment(record):
                n_residue_sites += 1
                name = f"{RESIDUE_SITE_PREFIX}{n_residue_sites:03d}"
                _new(name, "noncanonical", residue, names, types, anchors, explicit)
                side_chains[name], all_types[name] = {residue}, {"noncanonical"}
            continue
        if record.kind in ("protonate", "deprotonate"):
            atom = record.arguments["atom_name"]
            template = (record.residue_before, atom) in _TEMPLATE_PROTONATION
            if site is not None:
                explicit[site].update(record.added)
            elif not template and residue not in pending:
                pending.append(residue)
            continue
        if residue in pending:
            pending.remove(residue)
        # An attachment joins the site its anchor already belongs to.
        if site is None:
            site = record.name
            _new(site, record.type, residue, names, types, anchors, explicit)
            side_chains[site], all_types[site] = set(), set()
        all_types[site].add(record.type)
        explicit[site].update(record.added)
        anchor_atom = record.arguments["atom_name"]
        if residue not in fragment_residues:
            if residue not in anchors[site]:
                anchors[site].append(residue)
            if anchor_atom in BACKBONE:
                touched = (anchor_atom, *_BACKBONE_HYDROGENS.get(anchor_atom, ()))
                explicit[site].update(AtomKey(*residue, name) for name in touched)
            else:
                side_chains[site].add(residue)
        for chain_id, resnum, icode, _ in record.new_residues:
            fragment_residues[(chain_id, resnum, icode)] = site

    for residue in pending:
        n_residue_sites += 1
        name = f"{RESIDUE_SITE_PREFIX}{n_residue_sites:03d}"
        _new(name, "noncanonical", residue, names, types, anchors, explicit)
        side_chains[name], all_types[name] = {residue}, {"noncanonical"}

    sites = []
    order = {key: index for index, key in enumerate(keys)}
    for name in names:
        atoms = {atom for atom in explicit[name] if atom in present}
        for residue in side_chains[name]:
            atoms.update(a for a in residue_atoms.get(residue, ()) if a.atom_name not in BACKBONE)
        # Every atom of a residue the attachment created, including atoms a
        # reaction moved there from the protein (the thiol H of a maleimide).
        for residue, owner_name in fragment_residues.items():
            if owner_name == name:
                atoms.update(residue_atoms.get(residue, ()))
        ordered = sorted(atoms, key=order.__getitem__)
        sites.append(
            ModificationSite(
                name=name,
                type=types[name],
                types=frozenset(all_types[name]),
                residues=tuple(anchors[name]),
                input_residues=tuple(protein.input_residues[r] for r in anchors[name]),
                atoms=tuple(ordered),
                atom_indices=tuple(order[atom] for atom in ordered),
            )
        )
    return sites


def segments(protein: ModifiedProtein) -> list[str]:
    """Return each atom's segment, in the order ``save_pdb`` writes the atoms."""
    segment_of = {}
    for record in protein.records:
        for chain_id, resnum, icode, _ in record.new_residues:
            segment_of[(chain_id, resnum, icode)] = record.name
    return [
        segment_of.get((key.chain_id, key.resnum, key.icode), PROTEIN_SEGMENT)
        for key in protein.atom_keys()
    ]


def _is_side_chain_fragment(record) -> bool:
    """True for a mutation to a SMILES side chain, which is never an ff14SB residue."""
    target = record.arguments["to"]
    return isinstance(target, dict) and "smiles" in target


def _new(name, site_type, residue, names, types, anchors, explicit) -> None:
    names.append(name)
    types[name] = site_type
    anchors[name] = [residue]
    explicit[name] = set()


def _drop(name, names, types, anchors, explicit) -> None:
    names.remove(name)
    del types[name], anchors[name], explicit[name]
