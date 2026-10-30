"""Modify a protein through stable PolyzyMD names, recording every change.

The functions here call the mBuild biopolymers fork. Their names and
arguments belong to PolyzyMD, so code written against them keeps working
when mBuild's own API changes. Every call records what it did to the
atoms (added, removed and renamed atoms, the change in formal charge),
which is what parameterization and the build audit read later.

Mutations come first: ``mutate`` refuses a protein that already has a
protonation change or an attachment.
"""

from __future__ import annotations

import hashlib
import weakref
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, NamedTuple

from polyzymd.modifying_proteins.sites import Site, resolve, resolve_text

#: Attachment types; each has its own segment prefix and checks.
ATTACHMENT_TYPES = ("glycan", "polymer", "generic")

#: Name prefix of each attachment type. PSF segment IDs hold at most four
#: characters, so names are a letter and three digits (``G001``).
_NAME_PREFIX = {"glycan": "G", "polymer": "P", "generic": "M"}

# Source of each fragment built here, for records and YAML export.
_FRAGMENT_SOURCES: weakref.WeakKeyDictionary = weakref.WeakKeyDictionary()


class AtomKey(NamedTuple):
    """An atom named by chain, residue number, insertion code and atom name."""

    chain_id: str
    resnum: int
    icode: str
    atom_name: str

    def __str__(self) -> str:
        return f"{self.chain_id}:{self.resnum}{self.icode}:{self.atom_name}"


@dataclass(frozen=True)
class ModificationRecord:
    """What one call did: its arguments and its effect on the atoms."""

    kind: str
    site: Site
    arguments: dict
    residue_before: str
    residue_after: str
    added: tuple[AtomKey, ...]
    removed: tuple[AtomKey, ...]
    renamed: tuple[tuple[AtomKey, AtomKey], ...]
    new_residues: tuple[tuple[str, int, str, str], ...]
    formal_charge_change: int
    plan_text: str
    type: str | None = None
    name: str | None = None


@dataclass
class ModificationPlan:
    """What a modification will do, before it is applied."""

    kind: str
    site: Site
    arguments: dict
    text: str
    reaction: str | None = None
    formed: list[tuple[str, str, int]] = field(default_factory=list)
    leaving: list[str] = field(default_factory=list)
    charges: dict[str, int] = field(default_factory=dict)
    type: str | None = None
    name: str | None = None
    _mbuild_plan: Any = field(default=None, repr=False)

    def __str__(self) -> str:
        return self.text

    def draw(self, **kwargs):
        """Draw the site, the fragment and the product (mBuild's drawing)."""
        return self._mbuild_plan.draw(**kwargs)


class ModifiedProtein:
    """A protein and the ordered record of the modifications applied to it."""

    def __init__(self, protein, pdb_path: Path, ccd_paths: list[Path], download: bool):
        self.mbuild = protein
        self.records: list[ModificationRecord] = []
        self.source = {"pdb_path": str(pdb_path), "sha256": _sha256(pdb_path)}
        self._ccd_paths = list(ccd_paths)
        self._download = download
        self.input_residues = {
            (chain_id, residue.resnum, residue.icode): residue.name
            for chain_id, residue in _residues(protein)
        }

    def atom_keys(self) -> list[AtomKey]:
        """Return every atom's key, in the order ``save_pdb`` writes them."""
        keys = _snapshot(self.mbuild)
        return [keys[atom][0] for atom in self.mbuild.particles()]

    def save_pdb(self, path: str | Path) -> None:
        """Write the modified protein as a PDB file with CONECT records."""
        self.mbuild.save_pdb(str(path), overwrite=True)

    def _use_ccd_file(self, path: Path) -> None:
        from mbuild.biopolymers import CCDLibrary

        if path.parent not in self._ccd_paths:
            self._ccd_paths.insert(0, path.parent)
            self.mbuild.library = CCDLibrary(paths=self._ccd_paths, download=self._download)


def load_protein(
    pdb_path: str | Path, *, ccd_paths: list[str | Path] = (), download: bool = False
) -> ModifiedProtein:
    """Load a fully protonated protein PDB file.

    ``ccd_paths`` are folders of ``{CODE}.cif`` component files used for
    non-canonical residues; ``download`` lets mBuild fetch unknown codes
    from the RCSB.
    """
    from mbuild.biopolymers import CCDLibrary, Protein

    pdb_path = Path(pdb_path)
    paths = [Path(p) for p in ccd_paths]
    library = CCDLibrary(paths=paths, download=download)
    return ModifiedProtein(Protein(str(pdb_path), library=library), pdb_path, paths, download)


def fragment_from_smiles(smiles: str, resname: str):
    """Build a fragment from SMILES; ``*`` marks the atom that bonds."""
    from mbuild.biopolymers import prepare_fragment

    fragment = prepare_fragment(smiles, resname)
    _FRAGMENT_SOURCES[fragment] = {"smiles": smiles, "resname": resname}
    return fragment


def fragment_from_pdb(pdb_path: str | Path, charge: int = 0):
    """Read a fragment from a PDB file with CONECT records (e.g. a GLYCAM glycan)."""
    from mbuild.biopolymers import fragment_from_pdb as mbuild_fragment_from_pdb

    fragment = mbuild_fragment_from_pdb(str(pdb_path), charge=charge)
    source = {"pdb_path": str(pdb_path), "sha256": _sha256(Path(pdb_path))}
    if charge:
        source["charge"] = charge
    _FRAGMENT_SOURCES[fragment] = source
    return fragment


def fragment_from_ccd(code: str, link_atom: str, protein: ModifiedProtein | None = None):
    """Build a fragment from a CCD component, bonding through ``link_atom``."""
    from mbuild.biopolymers import fragment_from_ccd as mbuild_fragment_from_ccd

    library = protein.mbuild.library if protein is not None else None
    fragment = mbuild_fragment_from_ccd(code, link_atom, library=library)
    _FRAGMENT_SOURCES[fragment] = {"ccd": code, "link_atom": link_atom}
    return fragment


def fragment_source(fragment) -> dict | None:
    """Return how a fragment built here was made, or None."""
    return _FRAGMENT_SOURCES.get(fragment)


# ---------------------------------------------------------------------------
# Plans
# ---------------------------------------------------------------------------


def plan_mutate(
    protein: ModifiedProtein,
    site: str | Site,
    to,
    *,
    stereo: str | None = None,
    relax: bool = True,
    minimize: bool = True,
    platform: str | None = None,
) -> ModificationPlan:
    """Plan replacing the side chain at ``site``.

    ``to`` is a residue code (canonical, ``"PHE"``, or any CCD component,
    ``"4II"``), ``{"smiles": "*Cc1ccccc1", "resname": "PHX"}`` for a side
    chain whose ``*`` bonds to CA, or ``{"ccd_file": "path/XYZ.cif"}``.
    """
    _check_mutations_first(protein)
    site, residue = _site(protein, site)
    target, arguments = _mutation_target(protein, to)
    if stereo is not None:
        arguments["stereo"] = stereo
    plan = protein.mbuild.plan_mutate(
        residue.resnum,
        target,
        chain_id=site.chain_id,
        icode=site.icode,
        resname=to.get("resname") if isinstance(to, dict) else None,
        stereo=stereo,
        relax=relax,
        minimize=minimize,
        platform=platform,
    )
    return _view("mutate", site, arguments, plan)


def plan_protonate(protein: ModifiedProtein, site: str | Site, atom_name: str):
    """Plan adding a proton at ``atom_name`` of the residue at ``site``."""
    return _plan_protonation(protein, site, atom_name, "protonate")


def plan_deprotonate(protein: ModifiedProtein, site: str | Site, atom_name: str):
    """Plan removing a proton from ``atom_name`` of the residue at ``site``."""
    return _plan_protonation(protein, site, atom_name, "deprotonate")


def _plan_protonation(protein, site, atom_name, kind):
    site, residue = _site(protein, site)
    method = getattr(protein.mbuild, f"plan_{kind}")
    plan = method(residue.resnum, atom_name, chain_id=site.chain_id, icode=site.icode)
    if not (plan.leaving if kind == "deprotonate" else plan.new_hydrogens):
        reason = getattr(plan, "_state", {}).get("noop", "")
        raise ValueError(f"{kind} {site} {atom_name} would change nothing. {reason}".strip())
    return _view(kind, site, {"atom_name": atom_name}, plan)


def plan_attach(
    protein: ModifiedProtein,
    fragment,
    *,
    type: str,
    site: str | Site,
    atom_name: str,
    reaction: str | None = None,
    fragment_atom_name: str | None = None,
    fragment_resnum: int | None = None,
    fragment_resname: str | None = None,
    leaving_atom_names: str | list[str] | None = None,
    fragment_leaving_atom_names: str | list[str] | None = None,
    merge: bool = False,
    bond_order: int = 1,
    name: str | None = None,
    parameters: str | None = None,
    relax: bool = True,
    minimize: bool = True,
    platform: str | None = None,
) -> ModificationPlan:
    """Plan bonding ``fragment`` to ``atom_name`` of the residue at ``site``.

    ``type`` is ``"glycan"``, ``"polymer"`` or ``"generic"``. ``fragment``
    is a fragment built here, or its source written as data
    (``{"smiles": ..., "resname": ...}``, ``{"pdb_path": ...}`` or
    ``{"ccd": ..., "link_atom": ...}``). ``parameters`` names the parameter
    set a config declares for this attachment; it is recorded, not used here.
    The other arguments are mBuild's.
    """
    if type not in ATTACHMENT_TYPES:
        raise ValueError(f"Attachment type {type!r}; use one of {', '.join(ATTACHMENT_TYPES)}.")
    site, residue = _site(protein, site)
    if type == "glycan":
        source = fragment if isinstance(fragment, dict) else fragment_source(fragment)
        if source and "pdb_path" in source:
            from polyzymd.modifying_proteins.glycam import check_glycan_pdb

            check_glycan_pdb(source["pdb_path"])
    if isinstance(fragment, dict):
        fragment = _fragment_from_source(protein, fragment)
    source = fragment_source(fragment)
    arguments = {"atom_name": atom_name, "fragment": source}
    optional = {
        "reaction": reaction,
        "fragment_atom_name": fragment_atom_name,
        "fragment_resnum": fragment_resnum,
        "fragment_resname": fragment_resname,
        "leaving_atom_names": leaving_atom_names,
        "fragment_leaving_atom_names": fragment_leaving_atom_names,
    }
    arguments.update({key: value for key, value in optional.items() if value is not None})
    if merge:
        arguments["merge"] = True
    if bond_order != 1:
        arguments["bond_order"] = bond_order
    if parameters is not None:
        arguments["parameters"] = parameters
    plan = protein.mbuild.plan_attach(
        fragment,
        fragment_atom_name,
        resnum=residue.resnum,
        atom_name=atom_name,
        chain_id=site.chain_id,
        icode=site.icode,
        fragment_resnum=fragment_resnum,
        fragment_resname=fragment_resname,
        bond_order=bond_order,
        relax=relax,
        leaving_atom_names=leaving_atom_names,
        fragment_leaving_atom_names=fragment_leaving_atom_names,
        reaction=reaction,
        merge=merge,
        platform=platform,
        minimize=minimize,
    )
    view = _view("attach", site, arguments, plan)
    view.type = type
    view.name = name
    return view


# ---------------------------------------------------------------------------
# Applying
# ---------------------------------------------------------------------------


def apply(protein: ModifiedProtein, plan: ModificationPlan) -> ModificationRecord:
    """Apply a plan and append the record of what it did."""
    if plan._mbuild_plan.applied:
        raise ValueError("This plan was already applied; plan the modification again.")
    if plan.kind == "mutate":
        _check_mutations_first(protein)
    name = plan.name
    if plan.kind == "attach":
        taken = {record.name for record in protein.records}
        if name is None:
            name = _next_name(_NAME_PREFIX[plan.type], taken)
        elif name in taken:
            raise ValueError(f"A modification named {name!r} already exists.")
        elif not (name.isalnum() and len(name) <= 4):
            raise ValueError(
                f"Attachment name {name!r} must be one to four letters or digits: it "
                "becomes the attachment's PSF segment ID."
            )
    residue_before = resolve(protein.mbuild, plan.site).name
    before = _snapshot(protein.mbuild)
    charge_before = _chemical_charge(protein.mbuild)
    protein.mbuild.apply(plan._mbuild_plan)
    after = _snapshot(protein.mbuild)
    residue_after = protein.mbuild.get_residue(
        plan.site.resnum, chain_id=plan.site.chain_id, icode=plan.site.icode
    ).name
    added = [atom for atom in after if atom not in before]
    old_residues = {(key.chain_id, key.resnum, key.icode) for key, _ in before.values()}
    new_residues = sorted(
        {
            (key.chain_id, key.resnum, key.icode, resname)
            for key, resname in (after[atom] for atom in added)
            if (key.chain_id, key.resnum, key.icode) not in old_residues
        }
    )
    record = ModificationRecord(
        kind=plan.kind,
        site=plan.site,
        arguments=plan.arguments,
        residue_before=residue_before,
        residue_after=residue_after,
        added=tuple(after[atom][0] for atom in added),
        removed=tuple(before[atom][0] for atom in before if atom not in after),
        renamed=tuple(
            (before[atom][0], after[atom][0])
            for atom in before
            if atom in after and before[atom][0] != after[atom][0]
        ),
        new_residues=tuple(new_residues),
        formal_charge_change=_chemical_charge(protein.mbuild) - charge_before,
        plan_text=plan.text,
        type=plan.type,
        name=name,
    )
    protein.records.append(record)
    return record


def mutate(protein: ModifiedProtein, site: str | Site, to, **options) -> ModificationRecord:
    """Replace the side chain at ``site`` (see ``plan_mutate``)."""
    return apply(protein, plan_mutate(protein, site, to, **options))


def protonate(protein: ModifiedProtein, site: str | Site, atom_name: str) -> ModificationRecord:
    """Add a proton at ``atom_name`` of the residue at ``site``."""
    return apply(protein, plan_protonate(protein, site, atom_name))


def deprotonate(protein: ModifiedProtein, site: str | Site, atom_name: str) -> ModificationRecord:
    """Remove a proton from ``atom_name`` of the residue at ``site``."""
    return apply(protein, plan_deprotonate(protein, site, atom_name))


def attach(protein: ModifiedProtein, fragment, **options) -> ModificationRecord:
    """Bond a fragment to the protein (see ``plan_attach``)."""
    return apply(protein, plan_attach(protein, fragment, **options))


def relax_all(
    protein: ModifiedProtein,
    *,
    side_chains: bool = True,
    tolerance: float = 10.0,
    platform: str | None = None,
) -> None:
    """Relax every modified residue and attached fragment together.

    Place each modification with ``relax=False``, then call this once: one
    minimization with the protein backbone fixed (mBuild's
    ``relax_fragments``), as for a glycoprotein with many glycans.
    """
    touched = set()
    for record in protein.records:
        touched.add((record.site.chain_id, record.site.resnum, record.site.icode))
        touched.update((key.chain_id, key.resnum, key.icode) for key in record.added)
    residues = [
        residue
        for chain_id, residue in _residues(protein.mbuild)
        if (chain_id, residue.resnum, residue.icode) in touched
    ]
    if residues:
        protein.mbuild.relax_fragments(
            residues=residues, side_chains=side_chains, tolerance=tolerance, platform=platform
        )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _site(protein: ModifiedProtein, site: str | Site):
    """Return the Site and the mBuild residue it names, checked against the structure."""
    if isinstance(site, Site):
        return site, resolve(protein.mbuild, site)
    return resolve_text(protein.mbuild, site)


def _check_mutations_first(protein: ModifiedProtein) -> None:
    if any(record.kind != "mutate" for record in protein.records):
        raise ValueError(
            "Mutations come before protonation changes and attachments. Apply every "
            "mutation first."
        )


def _next_name(prefix: str, taken: set) -> str:
    number = 1
    while f"{prefix}{number:03d}" in taken:
        number += 1
    if number > 999:
        raise ValueError(f"More than 999 attachments of one type ({prefix}); name them.")
    return f"{prefix}{number:03d}"


def _mutation_target(protein: ModifiedProtein, to) -> tuple[Any, dict]:
    if isinstance(to, str):
        return to.upper(), {"to": to.upper()}
    if isinstance(to, dict) and "smiles" in to:
        from mbuild.biopolymers import prepare_fragment

        if "resname" not in to:
            raise ValueError("A SMILES side chain needs a `resname` (three characters or fewer).")
        return prepare_fragment(to["smiles"], to["resname"]), {
            "to": {"smiles": to["smiles"], "resname": to["resname"]}
        }
    if isinstance(to, dict) and "ccd_file" in to:
        path = Path(to["ccd_file"])
        protein._use_ccd_file(path)
        return path.stem.upper(), {"to": {"ccd_file": str(path), "sha256": _sha256(path)}}
    raise ValueError(
        f"Cannot mutate to {to!r}; give a residue code ('PHE', '4II'), "
        "{'smiles': ..., 'resname': ...} or {'ccd_file': ...}."
    )


def _fragment_from_source(protein: ModifiedProtein, source: dict):
    if "smiles" in source:
        return fragment_from_smiles(source["smiles"], source["resname"])
    if "pdb_path" in source:
        return fragment_from_pdb(source["pdb_path"], charge=source.get("charge", 0))
    if "ccd" in source:
        return fragment_from_ccd(source["ccd"], source["link_atom"], protein)
    if "ccd_file" in source:
        # A local CCD definition (a component the bundled library lacks), named
        # by its file; the protein's library reads the file's folder.
        path = Path(source["ccd_file"])
        protein._use_ccd_file(path)
        fragment = fragment_from_ccd(path.stem.upper(), source["link_atom"], protein)
        _FRAGMENT_SOURCES[fragment] = {
            "ccd_file": str(path),
            "sha256": _sha256(path),
            "link_atom": source["link_atom"],
        }
        return fragment
    raise ValueError(
        f"Cannot build a fragment from {source!r}; give smiles + resname, pdb_path, "
        "ccd + link_atom, or ccd_file + link_atom."
    )


def _view(kind: str, site: Site, arguments: dict, plan) -> ModificationPlan:
    return ModificationPlan(
        kind=kind,
        site=site,
        arguments=arguments,
        text=str(plan),
        reaction=plan.reaction,
        formed=[(_label(a), _label(b), order) for a, b, order in plan.formed],
        leaving=[_label(leaving) for _, leaving in plan.leaving],
        charges={_label(atom): charge for atom, charge in plan.charges.items()},
        _mbuild_plan=plan,
    )


def _label(atom) -> str:
    from mbuild.biopolymers import Residue

    parent = atom.parent
    if isinstance(parent, Residue):
        return f"{parent.name}{parent.resnum}{parent.icode}:{atom.name}"
    return atom.name


def _residues(protein):
    for chain in protein.chains:
        for residue in protein.residues(chain.chain_id):
            yield chain.chain_id, residue


def _snapshot(protein) -> dict:
    """Map each atom particle to its key and residue name."""
    keys = {}
    for chain_id, residue in _residues(protein):
        for atom in residue.particles():
            keys[atom] = (AtomKey(chain_id, residue.resnum, residue.icode, atom.name), residue.name)
    return keys


def _chemical_charge(protein) -> int:
    """Return the protein's net formal charge, summed over its chemistry.

    mBuild's ``net_formal_charge`` sums each residue's template charges, and
    a residue keeps them when an attachment removes a charged atom: a
    C-terminal amide formed by removing OXT still reports OXT's -1 there.
    The chemistry ``to_rdkit`` builds is correct, so it is the source.
    """
    return sum(atom.GetFormalCharge() for atom in protein.to_rdkit().GetAtoms())


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
