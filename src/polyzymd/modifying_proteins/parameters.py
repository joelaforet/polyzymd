"""Force-field parameters and charges for a modified protein.

The protein keeps its declared force field, one with residue library charges
(ff14SB's OpenFF port by default). Each modification site takes the parameter
set of its type (glycan, polymer, generic, noncanonical):

- Atoms outside every site keep the protein's parameters and library
  charges.
- The backbone of a modified residue (N, H, CA, HA, C, O) keeps the protein's
  parameters and the library charges of the residue as it was before any
  site-creating change (the input residue, or the canonical residue a
  mutation made), so it matches its neighbours. These charges come from the
  declared protein force field run on that residue, never from another
  copy of the force field (Joe, 2026-10-30). A force field that charges
  proteins with NAGL (OpenFF Rosemary) has no library charges and
  parameterizes the whole protein instead (:func:`parameterize_whole`).
- Site atoms (the side chain from CB outward plus everything attached) take
  the parameter set's force field. Their charges come from NAGL on a capped
  model of the site, shifted uniformly over the site atoms so that the site
  and the backbone of its residues sum to their formal charge.
- A change that touches a backbone atom (an attachment to the N-terminal
  amine, a neutral N-terminus, a protonated C-terminus) puts that atom and
  its hydrogens into the site; at the C-terminus the whole carboxyl and CA.
- A terminus changed by protonation alone (a neutral N-terminus, a protonated
  C-terminus) is refused in 1.4 (Joe, 2026-10-28); cap a truncated chain with
  ACE or NME instead.
- When a site touches a terminal residue's terminus (an N-terminal label, a
  C-terminal amide), that
  residue's untouched atoms, and the backbone charges of its site target,
  take the protein force field's main-chain template for the residue (ASN,
  not NASN or CASN): the reference caps that terminus with ACE or NME, as
  Amber does (Joe, 2026-10-28). An untouched terminus keeps its terminal
  template.
- A bonded term belongs to the protein only when all its atoms are outside
  every site; a term touching a site atom takes that site's force field
  (CA-CB, the angles at CA involving CB and chi1 are therefore the site's).
  A term touching sites of two force fields takes the first force field in
  the order of the parameter-set names.
- Bonds to hydrogen are constrained everywhere, with or without sites, by one
  force field's Constraints handler (``constraints_from``), each at the bond
  length of the force field that owns the bond.
- A glycan site whose parameter set is GLYCAM06j takes GLYCAM's own charges,
  terms and per-pair 1-4 scales instead (see :mod:`.glycam`).

The protein's terms come from a reference built by reloading the input PDB
and replaying only the changes that keep every residue an ff14SB residue:
canonical mutations and protonation changes to an ff14SB variant (LYN, ASH,
...) on residues outside every site. Residues outside the sites are then
identical in the reference and the modified protein, and each site residue
keeps the residue it had before the changes that made it a site.

Status: this ff14SB + Sage/NAGL site method is interim. It stays until a
benchmark against OpenFF Rosemary (issue #307) decides whether it is still
needed. The capped-model machinery belongs to this method.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from polyzymd.modifying_proteins import glycam
from polyzymd.modifying_proteins.api import ModifiedProtein, attach, load_protein, mutate
from polyzymd.modifying_proteins.api import deprotonate as _deprotonate
from polyzymd.modifying_proteins.api import protonate as _protonate
from polyzymd.modifying_proteins.sites import CANONICAL
from polyzymd.modifying_proteins.topology import (
    _TEMPLATE_PROTONATION,
    ModificationSite,
    modification_sites,
    segments,
)

#: The protein force field the protein atoms keep.
PROTEIN_FORCE_FIELD = ("ff14sb_off_impropers_0.0.4.offxml",)

#: Collections whose terms are split between the protein and the sites.
_SPLIT_COLLECTIONS = ("Bonds", "Angles", "ProperTorsions", "ImproperTorsions", "vdW", "Constraints")

#: Collections whose atom coverage is the same in every force field. Impropers
#: are defined differently by ff14SB and OpenFF, so only their ownership is
#: checked.
_COVERED_COLLECTIONS = ("Bonds", "Angles", "ProperTorsions", "vdW", "Constraints")

#: Guard-rail bound on the uniform charge shift of a site, in e per atom.
DEFAULT_MAX_CHARGE_SHIFT = 0.005

#: Force field whose Constraints handler (bonds to hydrogen) applies everywhere.
CONSTRAINTS_FROM = "openff-2.3.0.offxml"

#: Handlers that are not needed for library charges alone.
_VALENCE_HANDLERS = ("Bonds", "Angles", "ProperTorsions", "ImproperTorsions", "Constraints", "vdW")


@dataclass(frozen=True)
class ParameterSet:
    """The force field and charge model of one kind of site."""

    force_field: tuple[str, ...]
    charges: Mapping[str, str]

    def __post_init__(self) -> None:
        object.__setattr__(self, "force_field", tuple(self.force_field))
        if not self.force_field:
            raise ValueError("A parameter set needs at least one force field file.")
        if self.is_glycam:
            if (
                self.force_field != (glycam.GLYCAM_FORCE_FIELD,)
                or self.charges.get("method") != "library"
            ):
                raise ValueError(
                    f"A GLYCAM06j parameter set is force_field=['{glycam.GLYCAM_FORCE_FIELD}'] "
                    "with charges={'method': 'library'}: GLYCAM's residue templates carry "
                    "both."
                )
            return
        if self.charges.get("method") != "nagl" or not self.charges.get("model"):
            raise ValueError(
                "Site charges come from NAGL on a capped model: give "
                "charges={'method': 'nagl', 'model': '<model .pt file>'}."
            )

    @property
    def is_glycam(self) -> bool:
        return glycam.is_glycam(self.force_field)


@dataclass
class SiteParameters:
    """What one site received: its atoms, force field and charges."""

    name: str
    type: str
    atom_indices: list[int]
    backbone_residues: list[str]
    force_field: list[str]
    charges: dict
    capped_model_atoms: int | None = None  # None for GLYCAM06j sites
    nagl_sum: float | None = None
    target: float | None = None
    shift_total: float | None = None
    shift_per_atom: float | None = None


@dataclass
class ParameterizationRecord:
    """Plain data describing how a modified protein was parameterized."""

    protein_force_field: list[str]
    max_charge_shift: float
    total_charge: float
    segments: list[str]
    sites: list[SiteParameters] = field(default_factory=list)
    glycam: list[glycam.GlycamSite] = field(default_factory=list)
    #: [start, stop) of each molecule (connected component) of the protein, in
    #: atom order: chains joined by a disulfide or a modification are one.
    molecules: list[list[int]] = field(default_factory=list)

    @property
    def glycam_pairs(self) -> list[glycam.PairOverride]:
        """The 1-4 overrides of every GLYCAM06j site, in protein atom indices."""
        return [pair for site in self.glycam for pair in site.pairs]

    def to_dict(self) -> dict:
        data = asdict(self)
        data["glycam"] = [site.to_dict() for site in self.glycam]
        return data


@dataclass
class ParameterizedProtein:
    """The parameterized protein: its Interchange, record and export methods.

    Export through :meth:`to_openmm` and :meth:`to_gromacs`: they write each
    GLYCAM06j glycan's own 1-4 scales, which an Interchange cannot hold.
    :attr:`interchange_without_glycam_pairs` is the bare Interchange, whose
    own exports scale every 1-4 pair uniformly; :attr:`interchange` returns
    it only when there is no GLYCAM06j site. For tuple unpacking,
    ``interchange, record = parameterize(...)`` works the same way.
    """

    interchange_without_glycam_pairs: object
    record: ParameterizationRecord

    @property
    def glycam_pairs(self) -> list[glycam.PairOverride]:
        return self.record.glycam_pairs

    @property
    def interchange(self):
        if self.record.glycam:
            raise ValueError(
                "This protein has GLYCAM06j glycans, whose 1-4 pairs an Interchange "
                "cannot hold: export with to_openmm() / to_gromacs(), or take "
                "interchange_without_glycam_pairs and apply glycam_pairs yourself "
                "(glycam.apply_to_openmm / apply_to_gromacs)."
            )
        return self.interchange_without_glycam_pairs

    def __iter__(self):
        return iter((self.interchange, self.record))

    def to_openmm(self, **kwargs):
        """Return the OpenMM System, GLYCAM06j 1-4 pairs and masses included.

        ``kwargs`` go to ``Interchange.to_openmm``; nonbonded forces are
        always combined, because the pairs are NonbondedForce exceptions.
        """
        kwargs["combine_nonbonded_forces"] = True
        if self.record.glycam and kwargs.get("hydrogen_mass") not in (None, 1.007947):
            raise ValueError(
                "Hydrogen mass repartitioning with GLYCAM06j glycans is not "
                "supported yet. If you need it, open an issue at "
                "https://github.com/joelaforet/polyzymd/issues"
            )
        system = self.interchange_without_glycam_pairs.to_openmm(**kwargs)
        if self.record.glycam:
            glycam.apply_to_openmm(system, self.record.glycam)
        return system

    def to_gromacs(self, prefix: str, **kwargs) -> None:
        """Write ``{prefix}.top`` and ``{prefix}.gro``, GLYCAM06j pairs and masses included.

        ``kwargs`` go to ``Interchange.to_gromacs`` (``monolithic=False``
        writes the protein's molecule type to its own ``.itp``).
        """
        self.interchange_without_glycam_pairs.to_gromacs(prefix=prefix, **kwargs)
        if not self.record.glycam:
            return
        glycam.apply_to_gromacs_topology(
            Path(f"{prefix}.top"), self.record.glycam, self.molecule_ranges
        )

    @property
    def molecule_ranges(self) -> list[list[int]]:
        """``[start, stop)`` of each molecule of the protein, in atom order."""
        return self.record.molecules or [[0, len(self.record.segments)]]


def parameterize(
    protein: ModifiedProtein,
    parameter_sets: Mapping[str, ParameterSet],
    *,
    max_charge_shift: float = DEFAULT_MAX_CHARGE_SHIFT,
    protein_force_field: Sequence[str] = PROTEIN_FORCE_FIELD,
    constraints_from: str = CONSTRAINTS_FROM,
):
    """Return the parameterized protein (:class:`ParameterizedProtein`).

    ``parameter_sets`` maps each site type (``glycan``, ``polymer``,
    ``generic``, ``noncanonical``) to its :class:`ParameterSet`; a type with
    a site and no set raises. ``max_charge_shift`` bounds the uniform shift
    (e per atom) that brings each site to its integer target; a larger shift
    raises, because it means the capped model and the site disagree (the
    guard-rail default is 0.005 e per atom; a config can set another).
    ``constraints_from`` names the force field whose bonds-to-hydrogen
    Constraints handler is used everywhere (the builder passes its
    small-molecule force field).
    """
    from openff.toolkit import ForceField, Topology
    from openff.units import unit

    sites = modification_sites(protein)
    site_sets = resolve_parameter_sets(sites, parameter_sets)
    misused = sorted(
        site.name for site in sites if site_sets[site.name].is_glycam and site.type != "glycan"
    )
    if misused:
        raise ValueError(
            f"GLYCAM06j parameterizes glycan sites only; {', '.join(misused)} are not glycans."
        )

    molecule = to_openff(protein)
    reference = to_openff(reference_protein(protein, sites))
    library = library_charges(reference, protein_force_field)
    keys = _keys(molecule)
    reference_keys = _keys(reference)

    in_site = {index for site in sites for index in site.atom_indices}
    outside = [index for index in range(molecule.n_atoms) if index not in in_site]
    absent = [str(keys[index]) for index in outside if keys[index] not in library]
    if absent:
        raise ValueError(
            "These protein atoms have no counterpart in the protein force field's reference: "
            + ", ".join(absent[:10])
        )

    reference_residue = {
        (key.chain_id, key.resnum, key.icode): atom.metadata["residue_name"]
        for key, atom in zip(reference_keys, reference.atoms)
    }
    charges = np.array([library.get(key, 0.0) for key in keys])
    record = ParameterizationRecord(
        protein_force_field=list(protein_force_field),
        max_charge_shift=max_charge_shift,
        total_charge=0.0,
        segments=segments(protein),
    )
    native = []
    for site in sites:
        parameter_set = site_sets[site.name]
        if parameter_set.is_glycam:
            terms, glycam_site = glycam.glycam_site_terms(molecule, site, library, keys)
            native.append((terms, glycam_site))
            site_charges, info = [terms.charges[i] for i in site.atom_indices], {}
        else:
            site_charges, info = site_charges_from_capped_model(
                molecule, site, library, parameter_set.charges["model"], max_charge_shift
            )
        charges[list(site.atom_indices)] = site_charges
        record.sites.append(
            SiteParameters(
                name=site.name,
                type=site.type,
                atom_indices=list(site.atom_indices),
                backbone_residues=[reference_residue[residue] for residue in site.residues],
                force_field=list(parameter_set.force_field),
                charges=dict(parameter_set.charges),
                **info,
            )
        )
    formal = sum(atom.formal_charge.m_as(unit.elementary_charge) for atom in molecule.atoms)
    if not np.isclose(charges.sum(), formal, atol=1e-6):
        raise ValueError(f"Charges sum to {charges.sum():.6f}, not the formal charge {formal}.")
    record.total_charge = float(charges.sum())
    molecule.partial_charges = charges * unit.elementary_charge

    # One OpenFF Molecule per connected component (chain, or chains joined by
    # a bond); each is a contiguous range, so topology indices are the protein's.
    molecules = split_molecules(molecule)
    record.molecules = [[m.start, m.stop] for m in molecule_ranges(molecule)]
    topology = Topology.from_molecules(molecules)
    positions = molecule.conformers[0] if molecule.conformers else None
    constraints = ForceField(constraints_from)["Constraints"]
    protein_ff = _with_constraints(ForceField(*protein_force_field), constraints)
    # One group per site force field, in the order of the parameter-set names.
    groups: dict[tuple[str, ...], set[int]] = {}
    for site_type in sorted(parameter_sets):
        for site in sites:
            if site.type == site_type:
                parameter_set = site_sets[site.name]
                # GLYCAM sites are typed by a placeholder whose terms are then replaced.
                force_field = (
                    glycam.PLACEHOLDER_FORCE_FIELD
                    if parameter_set.is_glycam
                    else parameter_set.force_field
                )
                groups.setdefault(force_field, set()).update(site.atom_indices)
    if not groups:
        interchange = protein_ff.create_interchange(
            topology, charge_from_molecules=charge_templates(molecules)
        )
        interchange.positions = positions
        return ParameterizedProtein(interchange, record)

    site_interchanges = [
        (
            force_field,
            _with_constraints(ForceField(*force_field), constraints).create_interchange(
                topology, charge_from_molecules=charge_templates(molecules)
            ),
            atoms,
        )
        for force_field, atoms in groups.items()
    ]
    reference_interchange = protein_ff.create_interchange(
        Topology.from_molecules(split_molecules(reference))
    )
    position = {keys[index]: index for index in outside}
    mapping = {index: position[key] for index, key in enumerate(reference_keys) if key in position}
    interchange = _splice(site_interchanges, reference_interchange, mapping, set(outside))
    for terms, glycam_site in native:
        glycam.inject(interchange, terms, glycam_site)
        record.glycam.append(glycam_site)
    _check_charges(interchange, charges)
    interchange.positions = positions
    return ParameterizedProtein(interchange, record)


def resolve_parameter_sets(
    sites: Sequence[ModificationSite],
    parameter_sets: Mapping[str, ParameterSet],
) -> dict[str, ParameterSet]:
    """Return the parameter set of each site, by site name.

    A site can merge modifications of several types (a click onto a
    non-canonical residue is ``noncanonical`` and ``generic``). Their
    parameter sets must be identical; different ones raise, naming both,
    rather than one being picked.
    """
    resolved = {}
    for site in sites:
        names = sorted(site.types)
        missing = [name for name in names if name not in parameter_sets]
        if missing:
            raise ValueError(
                f"No parameter set for the site type(s) {', '.join(missing)} (site {site.name})."
            )
        chosen = parameter_sets[names[0]]
        for name in names[1:]:
            if parameter_sets[name] != chosen:
                raise ValueError(
                    f"Site {site.name} merges a {names[0]} and a {name} modification with "
                    f"different parameter sets ({_describe(chosen)} vs "
                    f"{_describe(parameter_sets[name])}). Give both types the same set."
                )
        resolved[site.name] = chosen
    return resolved


def _describe(parameter_set: ParameterSet) -> str:
    charges = ", ".join(f"{k}={v}" for k, v in sorted(parameter_set.charges.items()))
    return f"{' + '.join(parameter_set.force_field)} with {charges}"


def to_openff(protein: ModifiedProtein):
    """Return the OpenFF Molecule of ``protein``, atoms in ``save_pdb`` order.

    The chemistry (bond orders, formal charges, hydrogens) is mBuild's; each
    atom carries its chain, residue number, insertion code and name.
    """
    from openff.toolkit import Molecule

    rdmol = protein.mbuild.to_rdkit()
    molecule = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True, hydrogens_are_explicit=True)
    for atom, rdatom in zip(molecule.atoms, rdmol.GetAtoms()):
        info = rdatom.GetPDBResidueInfo()
        atom.name = info.GetName().strip()
        atom.metadata["chain_id"] = info.GetChainId().strip()
        atom.metadata["residue_number"] = str(info.GetResidueNumber())
        atom.metadata["insertion_code"] = info.GetInsertionCode().strip()
        atom.metadata["residue_name"] = info.GetResidueName().strip()
    expected = [tuple(key) for key in protein.atom_keys()]
    if [tuple(key) for key in _keys(molecule)] != expected:
        raise ValueError("mBuild's RDKit export does not follow the protein's atom order.")
    return molecule


def reference_protein(
    protein: ModifiedProtein, sites: Sequence[ModificationSite] | None = None
) -> ModifiedProtein:
    """Rebuild ``protein`` with only the changes that keep ff14SB residues.

    Canonical mutations are replayed; protonation changes to an ff14SB
    variant are replayed on residues outside every site. Everything else
    (attachments, non-canonical mutations, other protonation changes) is
    left out, so each site residue keeps the residue it had before.
    """
    if sites is None:
        sites = modification_sites(protein)
    anchors = {residue for site in sites for residue in site.residues}
    reference = load_protein(
        protein.source["pdb_path"], ccd_paths=protein._ccd_paths, download=protein._download
    )
    for record in protein.records:
        residue = (record.site.chain_id, record.site.resnum, record.site.icode)
        if record.kind == "mutate":
            # A canonical mutation is replayed even on a site residue: a
            # residue mutated and then modified keeps the mutant's backbone
            # charges (CYS for S13C followed by a maleimide; Joe, 2026-10-28).
            to = record.arguments["to"]
            if isinstance(to, str) and record.residue_after in CANONICAL:
                options = {"relax": False, "minimize": False}
                if "stereo" in record.arguments:
                    options["stereo"] = record.arguments["stereo"]
                mutate(reference, record.site, to, **options)
        # A site residue keeps the ff14SB charges of the residue as it was in
        # the input PDB (LYS for a deprotonated, then acylated lysine).
        elif record.kind in ("protonate", "deprotonate") and residue not in anchors:
            atom = record.arguments["atom_name"]
            if (record.residue_before, atom) in _TEMPLATE_PROTONATION:
                verb = _protonate if record.kind == "protonate" else _deprotonate
                verb(reference, record.site, atom)
    termini = modified_termini(reference, sites)
    _refuse_neutral_termini(protein, termini)
    for residue, terminus in termini.items():
        _cap(reference, residue, terminus)
    return reference


def molecule_ranges(molecule) -> list[range]:
    """The atom range of each connected component of ``molecule``, in atom order.

    PolyzyMD keeps a protein's chains (with what is bonded to them) in atom
    order, so every component is a contiguous range; a bond between chains
    that are not next to each other would break that, and raises.
    """
    import networkx as nx

    ranges = []
    for component in nx.connected_components(molecule.to_networkx()):
        start, stop = min(component), max(component) + 1
        if len(component) != stop - start:
            raise ValueError(
                f"Atoms {start}-{stop - 1} hold more than one molecule: a bond joins chains "
                "that are not next to each other in the PDB. Order the chains so that bonded "
                "chains are adjacent."
            )
        ranges.append(range(start, stop))
    return sorted(ranges, key=lambda r: r.start)


def split_molecules(molecule) -> list:
    """Return one OpenFF Molecule per connected component of ``molecule``, in atom order.

    Names, metadata, charges and the first conformer carry over.
    """
    from openff.toolkit import Molecule

    ranges = molecule_ranges(molecule)
    if len(ranges) == 1:
        return [molecule]
    pieces = []
    for atoms in ranges:
        piece = Molecule()
        for index in atoms:
            atom = molecule.atom(index)
            piece.add_atom(
                atom.atomic_number,
                atom.formal_charge,
                atom.is_aromatic,
                stereochemistry=atom.stereochemistry,
                name=atom.name,
                metadata=dict(atom.metadata),
            )
        for bond in molecule.bonds:
            if bond.atom1_index in atoms:
                piece.add_bond(
                    bond.atom1_index - atoms.start,
                    bond.atom2_index - atoms.start,
                    bond.bond_order,
                    bond.is_aromatic,
                    stereochemistry=bond.stereochemistry,
                )
        if molecule.partial_charges is not None:
            piece.partial_charges = molecule.partial_charges[atoms.start : atoms.stop]
        if molecule.conformers:
            piece.add_conformer(molecule.conformers[0][atoms.start : atoms.stop])
        pieces.append(piece)
    return pieces


#: Anchor atoms of an attachment that labels each terminus.
_TERMINAL_ANCHORS = {"N": frozenset({"N"}), "C": frozenset({"C", "O", "OXT"})}

NEUTRAL_TERMINUS_MESSAGE = (
    "PolyzyMD 1.4 does not model a neutral N-terminus (NH2) or C-terminus (COOH): "
    "ff14SB has no template for them, and protein termini are normally simulated charged "
    "at physiological pH. If your chain is truncated, cap it instead: attach ACE at the "
    "N-terminus or NME at the C-terminus. If you need a neutral terminus (e.g. a low- or "
    "high-pH study), please open an issue: https://github.com/joelaforet/polyzymd/issues"
)


def _refuse_neutral_termini(protein: ModifiedProtein, termini: Mapping) -> None:
    """Raise when a terminus is changed by protonation alone (Joe, 2026-10-28)."""
    for residue, terminus in termini.items():
        labelled = any(
            record.kind == "attach"
            and (record.site.chain_id, record.site.resnum, record.site.icode) == residue
            and record.arguments["atom_name"] in _TERMINAL_ANCHORS[terminus]
            for record in protein.records
        )
        if not labelled:
            chain_id, resnum, icode = residue
            raise ValueError(
                f"{NEUTRAL_TERMINUS_MESSAGE} (chain {chain_id} residue {resnum}{icode}, "
                f"{terminus}-terminus)"
            )


#: Backbone atoms whose presence in a site means the C-terminus is modified.
_C_TERMINAL_GROUP = frozenset({"C", "O", "OXT", "HXT"})


def modified_termini(reference: ModifiedProtein, sites: Sequence[ModificationSite]) -> dict:
    """Return the terminal residues whose terminus a site touches.

    A residue's N-terminus is modified when it is an N-terminal residue in
    ``reference`` (it has H1, H2 or H3) and a site holds its N; its
    C-terminus when it is a C-terminal residue in ``reference`` (it has OXT)
    and a site holds its C, O, OXT or HXT. The value is "N" or "C".
    """
    names: dict[tuple[str, int, str], set[str]] = {}
    for key in reference.atom_keys():
        names.setdefault((key.chain_id, key.resnum, key.icode), set()).add(key.atom_name)
    termini = {}
    for site in sites:
        for key in site.atoms:
            residue = (key.chain_id, key.resnum, key.icode)
            present = names.get(residue, set())
            if key.atom_name == "N" and present & {"H1", "H2", "H3"}:
                termini[residue] = "N"
            elif key.atom_name in _C_TERMINAL_GROUP and "OXT" in present:
                termini[residue] = "C"
    return termini


def _cap(reference: ModifiedProtein, residue: tuple[str, int, str], terminus: str) -> None:
    """Cap a terminus of ``reference`` with ACE or NME (Amber's convention).

    The capped residue is then a main-chain residue to the protein force
    field, so its library charges are the main-chain template's.
    """
    from polyzymd.modifying_proteins.sites import Site

    chain_id, resnum, icode = residue
    name = reference.mbuild.get_residue(resnum, chain_id=chain_id, icode=icode).name
    site = Site(chain_id=chain_id, code=name, resnum=resnum, icode=icode)
    if terminus == "N":
        _deprotonate(reference, site, "N")
        fragment = {"smiles": "*C(C)=O", "resname": "ACE"}
        attach(reference, fragment, type="generic", site=site, atom_name="N", relax=False)
    else:
        fragment = {"smiles": "*NC", "resname": "NME"}
        attach(
            reference,
            fragment,
            type="generic",
            site=site,
            atom_name="C",
            leaving_atom_names=["OXT"],
            relax=False,
        )


def charge_templates(molecules: Sequence) -> list:
    """One molecule per distinct chemistry, for ``charge_from_molecules``.

    Interchange refuses two isomorphic templates; copies of one chain (the
    protomers of a homotrimer) are identical in atom order, so the first copy
    stands for all. Copies whose charges differ raise.
    """
    import numpy as np

    templates: dict[tuple, Any] = {}
    for molecule in molecules:
        key = (
            tuple(atom.atomic_number for atom in molecule.atoms),
            tuple(atom.formal_charge.m for atom in molecule.atoms),
            tuple(sorted((b.atom1_index, b.atom2_index, b.bond_order) for b in molecule.bonds)),
        )
        first = templates.setdefault(key, molecule)
        if first is not molecule and molecule.partial_charges is not None:
            if not np.allclose(
                first.partial_charges.m, molecule.partial_charges.m, atol=1e-6
            ):
                raise ValueError(
                    "Two copies of one chain carry different charges; give them the same "
                    "modifications or make them differ in chemistry."
                )
    return list(templates.values())


def library_charges(molecule, force_field: Sequence[str] = PROTEIN_FORCE_FIELD) -> dict:
    """Return the protein force field's library charge of every atom, by key."""
    from openff.toolkit import ForceField, Topology
    from openff.units import unit

    charges_only = ForceField(*force_field)
    for handler in _VALENCE_HANDLERS:
        if handler in charges_only.registered_parameter_handlers:
            charges_only.deregister_parameter_handler(handler)
    interchange = charges_only.create_interchange(
        Topology.from_molecules(split_molecules(molecule))
    )
    values = np.zeros(molecule.n_atoms)
    for key, charge in interchange["Electrostatics"].charges.items():
        values[key.atom_indices[0]] = charge.m_as(unit.elementary_charge)
    return dict(zip(_keys(molecule), values))


def capped_model(molecule, core: set[int]):
    """Cut ``core`` and two bonds of neighbours out of ``molecule``, capped with H.

    Returns the capped Molecule and, in its atom order, the index in
    ``molecule`` of every atom cut out (the caps follow). The core may span
    any number of residues. A cut bond of order above one raises.
    """
    from openff.toolkit import Molecule
    from openff.units import unit

    neighbours = [set() for _ in range(molecule.n_atoms)]
    for bond in molecule.bonds:
        neighbours[bond.atom1_index].add(bond.atom2_index)
        neighbours[bond.atom2_index].add(bond.atom1_index)
    kept = set(core)
    for _ in range(2):
        kept |= {other for index in kept for other in neighbours[index]}
    kept_indices = sorted(kept)
    local_index = {index: order for order, index in enumerate(kept_indices)}

    local = Molecule()
    for index in kept_indices:
        atom = molecule.atom(index)
        local.add_atom(
            atom.atomic_number,
            atom.formal_charge.m_as(unit.elementary_charge),
            atom.is_aromatic,
            name=atom.name,
        )
    for bond in molecule.bonds:
        ends = [index for index in (bond.atom1_index, bond.atom2_index) if index in kept]
        if len(ends) == 2:
            local.add_bond(
                local_index[bond.atom1_index],
                local_index[bond.atom2_index],
                bond.bond_order,
                bond.is_aromatic,
            )
        elif len(ends) == 1:
            if bond.bond_order != 1:
                raise ValueError(
                    f"The capped model would cut a bond of order {bond.bond_order} at "
                    f"{_label(molecule.atom(ends[0]))}; widen the site."
                )
            cap = local.add_atom(1, 0, False, name="HC")
            local.add_bond(local_index[ends[0]], cap, 1, False)
    return local, kept_indices


def site_charges_from_capped_model(
    molecule, site: ModificationSite, library: Mapping, nagl_model: str, max_charge_shift: float
):
    """Return the site's charges and a summary of how they were made.

    NAGL charges the capped model of the site's residues and attached atoms.
    The site atoms keep their NAGL charges plus one uniform shift, so the site
    and the backbone of its residues sum to their formal charge with the
    backbone on the protein's library charges.
    """
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
    from openff.units import unit

    site_atoms = list(site.atom_indices)
    site_set = set(site_atoms)
    residues = set(site.residues)
    keys = _keys(molecule)
    backbone = [
        index
        for index, key in enumerate(keys)
        if (key.chain_id, key.resnum, key.icode) in residues and index not in site_set
    ]
    local, kept = capped_model(molecule, site_set | set(backbone))
    NAGLToolkitWrapper().assign_partial_charges(local, partial_charge_method=nagl_model)
    local_charges = local.partial_charges.m_as(unit.elementary_charge)
    order = {index: position for position, index in enumerate(kept)}
    nagl = np.array([local_charges[order[index]] for index in site_atoms])

    formal = sum(
        molecule.atom(index).formal_charge.m_as(unit.elementary_charge)
        for index in site_atoms + backbone
    )
    target = formal - sum(library[keys[index]] for index in backbone)
    shift_total = target - nagl.sum()
    shift = shift_total / len(site_atoms)
    if abs(shift) > max_charge_shift:
        raise ValueError(
            f"Site {site.name} needs a charge shift of {shift:+.4f} e per atom to reach "
            f"{target:+.4f} e (NAGL gives {nagl.sum():+.4f} e over {len(site_atoms)} atoms), "
            f"more than the bound of {max_charge_shift} e. Check the site's formal charges "
            "and protonation states."
        )
    info = {
        "capped_model_atoms": local.n_atoms,
        "nagl_sum": float(nagl.sum()),
        "target": float(target),
        "shift_total": float(shift_total),
        "shift_per_atom": float(shift),
    }
    return nagl + shift, info


def _with_constraints(force_field, constraints):
    """Return ``force_field`` with ``constraints`` registered if it has none."""
    if "Constraints" not in force_field.registered_parameter_handlers:
        force_field.register_parameter_handler(constraints)
    return force_field


def _splice(site_interchanges, reference, mapping: dict[int, int], outside: set[int]):
    """Combine the protein's and each site force field's terms in one Interchange.

    ``site_interchanges`` holds (force field, Interchange, site atoms) per
    site force field, in priority order, each an Interchange of the whole
    modified protein. A term wholly outside every site comes from
    ``reference`` (the protein force field), mapped through ``mapping``; any
    other term comes from the first group whose site atoms it touches.
    Potential ids are prefixed with their source, so two force fields that
    share a parameter id keep their own values. The first Interchange is
    returned, holding the result.
    """
    target = site_interchanges[0][1]
    for name in ("vdW", "Electrostatics"):
        for attr in ("scale_12", "scale_13", "scale_14", "scale_15"):
            for _, interchange, _ in site_interchanges:
                if not np.isclose(
                    getattr(interchange[name], attr), getattr(reference[name], attr), atol=1e-6
                ):
                    raise ValueError(
                        f"The protein and site force fields differ in {name} {attr}; "
                        "they cannot share one system."
                    )
    for _, interchange, _ in site_interchanges:
        if interchange["vdW"].mixing_rule != reference["vdW"].mixing_rule:
            raise ValueError("The protein and site force fields use different mixing rules.")

    for name in _SPLIT_COLLECTIONS:
        sources = [ic for _, ic, _ in site_interchanges] + [reference]
        present = [name in ic.collections for ic in sources]
        if not any(present):
            continue
        if not all(present):
            lacking = [
                "+".join(ff) for (ff, ic, _) in site_interchanges if name not in ic.collections
            ]
            raise ValueError(
                f"Force field(s) {', '.join(lacking) or 'of the protein'} have no {name} "
                "handler; every force field in one system must define it."
            )
        key_map, potentials = {}, {}
        for force_field, interchange, atoms in site_interchanges:
            label = "+".join(force_field)
            earlier = set()
            for ff, _, a in site_interchanges:
                if ff == force_field:
                    break
                earlier |= a
            collection = interchange[name]
            for key, potential_key in collection.key_map.items():
                touched = set(key.atom_indices)
                if not touched & atoms or touched & earlier:
                    continue
                new_potential = potential_key.model_copy(
                    update={"id": f"site:{label}::{potential_key.id}"}
                )
                key_map[key] = new_potential
                potentials[new_potential] = collection.potentials[potential_key]
        source = reference[name]
        for key, potential_key in source.key_map.items():
            mapped = tuple(mapping.get(i) for i in key.atom_indices)
            if None in mapped or not set(mapped) <= outside:
                continue
            new_potential = potential_key.model_copy(update={"id": "protein::" + potential_key.id})
            key_map[key.model_copy(update={"atom_indices": mapped})] = new_potential
            potentials[new_potential] = source.potentials[potential_key]

        collection = target[name]
        original = {tuple(sorted(key.atom_indices)) for key in collection.key_map}
        collection.key_map.clear()
        collection.key_map.update(key_map)
        collection.potentials.clear()
        collection.potentials.update(potentials)
        for key, potential_key in key_map.items():
            from_protein = potential_key.id.startswith("protein::")
            if from_protein != (set(key.atom_indices) <= outside):
                raise ValueError(f"A {name} term on atoms {key.atom_indices} has the wrong owner.")
        if name in _COVERED_COLLECTIONS:
            final = {tuple(sorted(key.atom_indices)) for key in key_map}
            if final != original:
                missing, extra = len(original - final), len(final - original)
                raise ValueError(
                    f"The combined {name} terms do not cover the protein "
                    f"({missing} missing, {extra} extra)."
                )
    return target


def _check_charges(interchange, expected: np.ndarray) -> None:
    from openff.units import unit

    written = np.zeros(len(expected))
    for key, charge in interchange["Electrostatics"].charges.items():
        written[key.atom_indices[0]] = charge.m_as(unit.elementary_charge)
    if not np.allclose(written, expected, atol=1e-10):
        raise ValueError("A force field overwrote the assigned charges.")


def _keys(molecule):
    from polyzymd.modifying_proteins.api import AtomKey

    return [
        AtomKey(
            atom.metadata["chain_id"],
            int(atom.metadata["residue_number"]),
            atom.metadata["insertion_code"],
            atom.name,
        )
        for atom in molecule.atoms
    ]


def _label(atom) -> str:
    data = atom.metadata
    return f"{data['chain_id']}:{data['residue_name']}{data['residue_number']} {atom.name}"
