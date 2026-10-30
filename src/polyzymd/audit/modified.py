"""The expected terms of a covalently modified protein, rebuilt from its inputs.

A modified protein is one component with several domains: the protein
(ff14SB, library charges) and one domain per modification site (a parameter
set: an OpenFF force field with NAGL charges, or GLYCAM06j). Every term of the
realized system has one owning domain by the boundary rule of
:mod:`polyzymd.core.components`: the domain all its atoms belong to, else the
site it touches; a term touching two sites has no owner.

The expected value of each term comes from its owner, computed here from the
component record's reference alone:

- the protein is rebuilt from the input PDB file by replaying the modification
  blocks (with ``polyzymd.modifying_proteins``, which reads the chemistry from
  mBuild, and without relaxation; parameters do not depend on coordinates);
- **protein** terms and library charges: an OpenFF Interchange of the
  protein with only ff14SB residues (canonical mutations and ff14SB-variant
  protonations replayed), built here with the declared protein force field.
  A modified terminus is capped with ACE or NME in that reference, so it
  takes the force field's main-chain charges; no other source of protein
  charges is read;
- **OpenFF site** terms: an Interchange of the whole modified protein with the
  site's force field; site charges: NAGL with the declared model on a capped
  model built here (the site, two bonds of neighbours, every cut bond capped
  with H), shifted uniformly so that the site and the rest of its residues
  sum to their formal charge, the rest of its residues taking ff14SB library
  charges;
- **GLYCAM06j site** terms, charges and 1-4 pairs: OpenMM's own
  ``ForceField(protein.ff14SB.xml, GLYCAM_06j-1.xml)`` on the protein with its
  glycans only (other modifications left out), linkage residues renamed
  NLN/OLS/OLT;
- **1-4 pairs** of the protein and of OpenFF sites: the owner's 1-4 scales
  applied to the expected charges and Lennard-Jones parameters (1-2 and 1-3
  pairs are exclusions).

What is shared with the build: reading the modified protein's chemistry
(``to_openff``), deriving its sites (``modification_sites``) and replaying
the ff14SB-only reference (``reference_protein``). The values (force-field
assignment, charges, the shift and every 1-4 pair) are computed here.

Status: this module checks the interim ff14SB + Sage/NAGL method alongside
the whole-protein path. It may change after the OpenFF Rosemary benchmark
(issue #307).
"""

from __future__ import annotations

import collections
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np

from polyzymd.audit.molecules import reference_system
from polyzymd.audit.terms import TERM_KINDS, SystemIndex, compare
from polyzymd.utils.mbuild_logs import quiet_mbuild

LOGGER = logging.getLogger(__name__)

#: Atoms that are not a residue's side chain (the side chain runs from CB outward).
BACKBONE = frozenset({"N", "H", "CA", "HA", "HA2", "HA3", "C", "O", "H1", "H2", "H3", "OXT", "HXT"})

#: The atoms a change to a backbone atom brings into a site:
#: the touched atom and its hydrogens; a change to the C-terminal carboxylate
#: takes the whole carboxyl and CA, which share ff14SB's terminal charge.
TOUCHED_BACKBONE = {
    "N": ("N", "H", "H1", "H2", "H3"),
    "CA": ("CA", "HA", "HA2", "HA3"),
    "OXT": ("OXT", "HXT", "O", "C", "CA", "HA", "HA2", "HA3"),
}

#: Protonation changes that give an ff14SB residue variant (LYN, ASH, ...) and so no site.
FF14SB_VARIANTS = frozenset(
    {
        ("LYS", "NZ"),
        ("ASP", "OD1"),
        ("ASP", "OD2"),
        ("GLU", "OE1"),
        ("GLU", "OE2"),
        ("HIS", "ND1"),
        ("HIS", "NE2"),
        ("CYS", "SG"),
    }
)

#: Term kinds compared per owned term, beside the per-atom charge, LJ and mass.
_VALENCE = ("bonds", "angles", "torsions", "constraints")

#: Side chains a CCD component brings in a state other than the usual one at pH 7:
#: residue, the atom that marks it, the deprotonation that gives the usual state.
_UNUSUAL = {
    "ASP": ("HD2", "deprotonate", "OD2"),
    "GLU": ("HE2", "deprotonate", "OE2"),
}


@dataclass
class SiteCharges:
    """How one site's charges were expected: raw sum, target, shift."""

    name: str
    type: str
    atoms: list[int]  # realized indices
    glycam: bool
    target: float
    actual_sum: float
    nagl_sum: Optional[float] = None
    shift_per_atom: Optional[float] = None
    bound: Optional[float] = None


@dataclass
class ModifiedReference:
    """Everything the checks of a modified protein compare with."""

    component: str
    domain_of: dict[int, str]
    expected: dict[str, dict[str, dict]]
    sites: list[SiteCharges]
    #: Backbone atoms of modified residues, realized index -> expected ff14SB charge.
    backbone: dict[int, float]
    #: Expected bonds of the component, as sorted realized index pairs.
    bonds: set[tuple[int, int]]
    #: Expected atom names in realized order (chain:resnum:name), for the applied check.
    names: list[str]
    site_types: dict[str, str]
    site_force_fields: dict[str, tuple[str, ...]]
    #: Atoms the site rule puts in a site, from the modification records alone.
    rule_atoms: set[int] = field(default_factory=set)
    #: Backbone atoms the rule puts in a site (a change touched them).
    designed_backbone: set[int] = field(default_factory=set)
    #: Atoms a change added or an attachment brought (not protein residue atoms).
    attached_atoms: set[int] = field(default_factory=set)
    #: Per site: (formal charge of its residues read from the chemistry, the
    #: formal charge of the input residues plus the recorded changes), or None.
    site_formal: dict[str, Optional[tuple[float, float]]] = field(default_factory=dict)
    #: The rebuilt protein's modification records and the component's recipe.
    records: list = field(default_factory=list)
    recipe: dict = field(default_factory=dict)
    #: How the chemistry the build used differs from the one read independently.
    chemistry: list[str] = field(default_factory=list)
    #: Termini changed by protonation alone, which the ff14SB path refuses.
    neutral_termini: list[str] = field(default_factory=list)
    protonation_warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)

    def owner(self, atoms: tuple[int, ...]) -> Optional[str]:
        """The domain a term over ``atoms`` (realized indices) belongs to, or None."""
        domains = [self.domain_of.get(i) for i in atoms]
        if any(d is None for d in domains):
            return None
        sites = {d for d in domains if d != "protein"}
        if len(sites) > 1:
            return None
        return sites.pop() if sites else "protein"

    def conflicts(self, index: SystemIndex) -> list[tuple[str, tuple[int, ...]]]:
        """Terms of the realized system that touch two sites of this component."""
        found = []
        molecules = {index.molecule_of[i] for i in self.domain_of}
        for molecule in molecules:
            for kind, entries in index.terms.get(molecule, {}).items():
                for atoms, _ in entries:
                    sites = {self.domain_of.get(i) for i in atoms} - {"protein", None}
                    if len(sites) > 1:
                        found.append((kind, atoms))
        return found


#: Capping residues of ff14SB that a GLYCAM06j template also matches.
_PROTEIN_CAPS = frozenset({"ACE", "NME", "NHE"})


def is_modified(component: Any) -> bool:
    """Whether ``component`` is a modified protein (its domains carry a rebuild recipe)."""
    return any("modified_protein" in (domain.reference or {}) for domain in component.domains)


def build_reference(realized: Any, component: Any) -> ModifiedReference:
    """Rebuild ``component``'s modified protein and compute every expected term.

    The rebuild places the modifications without relaxation, by design (only
    the chemistry and the parameters are compared); mBuild's warnings about the
    close contacts that leaves describe the rebuild, not the built structure,
    so they are muted.
    """
    LOGGER.info(
        "The audit rebuilds the modified protein without relaxation (by design); mBuild's "
        "placement warnings from that rebuild are muted."
    )
    with quiet_mbuild():
        return _build_reference(realized, component)


def _build_reference(realized: Any, component: Any) -> ModifiedReference:
    from polyzymd.modifying_proteins.parameters import reference_protein, to_openff
    from polyzymd.modifying_proteins.topology import modification_sites

    recipe = next(
        domain.reference["modified_protein"]
        for domain in component.domains
        if "modified_protein" in (domain.reference or {})
    )
    offset = min(component.atom_indices())
    protein = rebuild(recipe)
    # The chemistry is read again, independently, from the PDB file; the
    # expected values below use that reading.
    built_molecule = to_openff(protein)
    try:
        from polyzymd.audit import _pablo

        molecule = _pablo.perceive(protein)
        chemistry = chemistry_differences(molecule, built_molecule, protein.atom_keys())
    except Exception as error:  # Pablo raises many types; any means no independent reading
        molecule = built_molecule
        chemistry = [f"could not read the chemistry independently: {type(error).__name__}: {error}"]
    keys = [tuple(key) for key in protein.atom_keys()]
    n_atoms = len(keys)
    if n_atoms != component.n_atoms:
        raise ValueError(
            f"rebuilding from {recipe['pdb']} and the modification blocks gives {n_atoms} "
            f"atoms, not the component's {component.n_atoms}"
        )
    elements = realized.elements[offset : offset + n_atoms]
    expected_elements = [atom.symbol for atom in molecule.atoms]
    if elements != expected_elements:
        raise ValueError(
            "the component's atoms are not the rebuilt protein's atoms, in order "
            f"(first difference at atom {_first_difference(elements, expected_elements) + offset})"
        )
    real = lambda local: local + offset  # noqa: E731
    if recipe.get("protein_charges"):
        ref = _whole_protein_reference(recipe, component, protein, molecule, keys, real)
        ref.chemistry = chemistry
        return ref
    from polyzymd.config.schema import protein_charge_handlers

    if not protein_charge_handlers(tuple(recipe["protein_force_field"]))[1]:
        raise ValueError(
            f"the protein declares library charges, but {', '.join(recipe['protein_force_field'])} "
            "has no protein library charges"
        )

    sites = modification_sites(protein)
    site_by_name = {site.name: site for site in sites}
    domain_names = [domain.name for domain in component.domains]
    domain_of = {real(i): "protein" for i in range(n_atoms)}
    for site in sites:
        if site.name not in domain_names:
            raise ValueError(
                f"the rebuilt protein has a site {site.name} the component does not list"
            )
        for local in site.atom_indices:
            domain_of[real(local)] = site.name

    sets = recipe["parameter_sets"]
    site_types = {site.name: site.type for site in sites}
    site_sets = {site.name: _site_set(protein, site, sets) for site in sites}
    site_force_fields = {name: tuple(s["force_field"]) for name, s in site_sets.items()}

    # The protein: the declared protein force field on the reference, the
    # input protein with only the residues that force field has templates for.
    reference = to_openff(reference_protein(protein, sites))
    reference_keys = [_key(atom) for atom in reference.atoms]
    position = {key: i for i, key in enumerate(keys)}
    input_residues = set(protein.input_residues)
    # Only atoms of the input's residues: a reference may also hold terminal caps.
    to_real_from_reference = {
        k: real(position[key])
        for k, key in enumerate(reference_keys)
        if key in position and key[:3] in input_residues
    }
    protein_system = reference_system(
        recipe["protein_force_field"],
        reference.to_topology(),
        constraints_from=recipe["constraints_from"],
    )
    protein_index = SystemIndex.from_system(protein_system)
    library = {
        to_real_from_reference[k]: protein_index.particle[k][0] for k in to_real_from_reference
    }
    # A terminus a site touches takes the main-chain template (ASN, not NASN
    # or CASN): the reference caps it with ACE or NME, so the declared protein
    # force field charges it as in the middle of a chain. The termini are
    # found here from the input PDB, and each must be capped in the reference.
    termini = modified_termini(recipe["pdb"], keys, domain_of, offset)
    uncapped = uncapped_termini(reference, termini)
    if uncapped:
        raise ValueError(
            "the reference does not cap these modified termini with ACE or NME: "
            + ", ".join(uncapped)
        )

    expected: dict[str, dict[str, dict]] = {
        name: {kind: {} for kind in TERM_KINDS} for name in domain_names
    }
    per_atom: dict[int, tuple[float, tuple[float, float], float]] = {}
    for k, i in to_real_from_reference.items():
        if domain_of[i] == "protein":
            q, lj, mass = protein_index.particle[k]
            per_atom[i] = (library[i], lj, mass)
    _owned_terms(expected, protein_index, to_real_from_reference, domain_of, "protein")

    # Charges of every site.
    site_charges: list[SiteCharges] = []
    charges = {i: per_atom[i][0] for i in per_atom}
    residue_of = {real(i): keys[i][:3] for i in range(n_atoms)}
    atoms_of_residue: dict[tuple, list[int]] = collections.defaultdict(list)
    for i, residue in residue_of.items():
        atoms_of_residue[residue].append(i)
    formal = {
        real(i): atom.formal_charge.m_as("elementary_charge")
        for i, atom in enumerate(molecule.atoms)
    }

    glycam_sites = [s for s in sites if _is_glycam(site_sets[s.name])]
    native = None
    if glycam_sites:
        native = _native_glycam(recipe, protein, glycam_sites, keys, real)
    for site in sites:
        atoms = [real(i) for i in site.atom_indices]
        residues = {residue_of[i] for i in atoms}
        rest = [i for r in residues for i in atoms_of_residue[r] if domain_of[i] != site.name]
        missing = [i for i in rest if i not in library]
        if missing:
            raise ValueError(
                f"site {site.name}: atoms of its residues outside the site have no ff14SB "
                f"library charge ({', '.join(_label(keys[i - offset]) for i in missing[:5])})"
            )
        target = sum(formal[i] for r in residues for i in atoms_of_residue[r]) - sum(
            library[i] for i in rest
        )
        actual_sum = float(sum(realized.index.particle[i][0] for i in atoms))
        if _is_glycam(site_sets[site.name]):
            for i in atoms:
                charges[i] = native["particle"][i][0]
            site_charges.append(
                SiteCharges(site.name, site.type, atoms, True, float(target), actual_sum)
            )
            continue
        nagl = capped_model_charges(
            molecule, list(site.atom_indices), [i - offset for i in rest], site_sets[site.name]
        )
        raw = float(sum(nagl.values()))
        shift = (target - raw) / len(atoms)
        for local, value in nagl.items():
            charges[real(local)] = value + shift
        site_charges.append(
            SiteCharges(
                site.name,
                site.type,
                atoms,
                False,
                float(target),
                actual_sum,
                nagl_sum=raw,
                shift_per_atom=shift,
                bound=recipe["max_charge_shift"],
            )
        )

    # OpenFF site terms: the site's force field on the whole modified protein.
    from openff.toolkit import Molecule
    from openff.units import unit

    charged = Molecule(molecule)
    charged.partial_charges = (
        np.array([charges.get(real(i), 0.0) for i in range(n_atoms)]) * unit.elementary_charge
    )
    to_real_local = {i: real(i) for i in range(n_atoms)}
    scales: dict[str, tuple[float, float]] = {"protein": _scales(recipe["protein_force_field"])}
    by_force_field: dict[tuple[str, ...], list[str]] = collections.defaultdict(list)
    for site in sites:
        if not _is_glycam(site_sets[site.name]):
            by_force_field[site_force_fields[site.name]].append(site.name)
    for force_field, names in by_force_field.items():
        system = reference_system(
            force_field,
            charged.to_topology(),
            [charged],
            constraints_from=recipe["constraints_from"],
        )
        index = SystemIndex.from_system(system)
        for name in names:
            scales[name] = _scales(force_field)
            for i in site_by_name[name].atom_indices:
                q, lj, mass = index.particle[i]
                per_atom[real(i)] = (charges[real(i)], lj, mass)
            _owned_terms(expected, index, to_real_local, domain_of, name)

    # GLYCAM06j site terms.
    if native is not None:
        for site in glycam_sites:
            for local in site.atom_indices:
                per_atom[real(local)] = native["particle"][real(local)]
            _owned_terms(expected, native["index"], native["to_real"], domain_of, site.name)

    for i, (q, lj, mass) in per_atom.items():
        domain = domain_of[i]
        expected[domain]["charge"][i] = charges.get(i, q)
        expected[domain]["lj"][i] = lj
        expected[domain]["mass"][i] = mass

    bonds = {tuple(sorted((real(b.atom1_index), real(b.atom2_index)))) for b in molecule.bonds}
    ref = ModifiedReference(
        component=component.name,
        domain_of=domain_of,
        expected=expected,
        sites=site_charges,
        backbone={},
        bonds=bonds,
        names=[_label(key) for key in keys],
        site_types=site_types,
        site_force_fields=site_force_fields,
    )
    _expected_pairs(ref, per_atom, charges, scales, native)
    for site in sites:
        for residue in site.residues:
            for i in atoms_of_residue[residue]:
                if keys[i - offset][3] in BACKBONE and domain_of[i] == "protein":
                    ref.backbone[i] = library[i]
    ref.rule_atoms, ref.designed_backbone, ref.attached_atoms = site_rule(protein, keys, offset)
    ref.records = list(protein.records)
    ref.neutral_termini = neutral_termini(termini, ref.records)
    ref.chemistry = chemistry
    ref.recipe = recipe
    input_formal = _input_formal_charges(recipe["pdb"])
    for site in sites:
        residues = {residue_of[i] for i in (real(k) for k in site.atom_indices)}
        chemistry = float(sum(formal[i] for r in residues for i in atoms_of_residue[r]))
        if input_formal is None:
            ref.site_formal[site.name] = None
            continue
        recorded = sum(input_formal.get(r, 0.0) for r in residues) + sum(
            record.formal_charge_change
            for record in protein.records
            if (record.site.chain_id, record.site.resnum, record.site.icode) in residues
            or any(tuple(n[:3]) in residues for n in record.new_residues)
        )
        ref.site_formal[site.name] = (chemistry, float(recorded))
    ref.protonation_warnings = protonation_warnings(protein)
    return ref


def _whole_protein_reference(
    recipe: dict, component: Any, protein: Any, molecule: Any, keys: list, real
) -> ModifiedReference:
    """The reference of a protein one force field parameterizes whole (Rosemary alpha).

    Every atom, modifications included, belongs to the protein domain: its
    terms come from the declared force field on ``molecule`` (the chemistry
    read independently with openff-pablo), and its charges from NAGL with the
    declared model on each of its molecules, both computed here.
    """
    from openff.toolkit import Topology
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper

    from polyzymd.modifying_proteins.parameters import charge_templates, split_molecules

    n_atoms = len(keys)
    # NAGL brings each molecule to its own formal charge: one call per
    # connected component (chain, or chains joined by a bond), in atom order.
    charged = []
    for piece in split_molecules(molecule):
        NAGLToolkitWrapper().assign_partial_charges(
            piece, partial_charge_method=recipe["protein_charges"]["model"]
        )
        charged.append(piece)
    system = reference_system(
        recipe["protein_force_field"],
        Topology.from_molecules(charged),
        charge_templates(charged),
        constraints_from=recipe["constraints_from"],
    )
    index = SystemIndex.from_system(system)
    to_real = {i: real(i) for i in range(n_atoms)}
    domain_of = {real(i): "protein" for i in range(n_atoms)}
    expected: dict[str, dict[str, dict]] = {"protein": {kind: {} for kind in TERM_KINDS}}
    _owned_terms(expected, index, to_real, domain_of, "protein")
    per_atom = {real(i): index.particle[i] for i in range(n_atoms)}
    charges = {i: per_atom[i][0] for i in per_atom}
    for i, (q, lj, mass) in per_atom.items():
        expected["protein"]["charge"][i] = q
        expected["protein"]["lj"][i] = lj
        expected["protein"]["mass"][i] = mass
    ref = ModifiedReference(
        component=component.name,
        domain_of=domain_of,
        expected=expected,
        sites=[],
        backbone={},
        bonds={tuple(sorted((real(b.atom1_index), real(b.atom2_index)))) for b in molecule.bonds},
        names=[_label(key) for key in keys],
        site_types={},
        site_force_fields={},
    )
    _expected_pairs(
        ref, per_atom, charges, {"protein": _scales(recipe["protein_force_field"])}, None
    )
    # One domain by design: no site atoms (rule_atoms stays empty). The records
    # are compared with the config's blocks by modifications.declared.
    ref.records = list(protein.records)
    ref.recipe = recipe
    # A neutral terminus (NH2, COOH) is allowed on Rosemary, which needs no
    # residue templates, so ref.neutral_termini stays empty.
    ref.protonation_warnings = protonation_warnings(protein)
    return ref


def chemistry_differences(perceived: Any, built: Any, keys: list) -> list[str]:
    """Where the build's chemistry differs from the independent reading, atom by atom.

    Compares elements, formal charges, which atoms are bonded, which bonds are
    aromatic (RDKit's perception, the same on both sides) and the orders of
    the other bonds, so two Kekule forms of one aromatic ring agree. Returns
    at most a few differences of each kind, with counts.
    """
    label = [str(key) for key in keys]
    problems = []
    elements = [
        f"{label[i]} is {b.symbol}, read as {a.symbol}"
        for i, (a, b) in enumerate(zip(perceived.atoms, built.atoms))
        if a.atomic_number != b.atomic_number
    ]
    charges = [
        f"{label[i]} has formal charge {b.formal_charge.m:+g}, read as {a.formal_charge.m:+g}"
        for i, (a, b) in enumerate(zip(perceived.atoms, built.atoms))
        if a.formal_charge.m != b.formal_charge.m
    ]

    def bonds(molecule: Any) -> dict:
        """Each bond: "aromatic", or its order."""
        out = {}
        for bond in molecule.to_rdkit().GetBonds():
            i, j = sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
            out[(i, j)] = "aromatic" if bond.GetIsAromatic() else bond.GetBondTypeAsDouble()
        return out

    read, used = bonds(perceived), bonds(built)
    bond_problems = []
    for pair in sorted(set(read) | set(used)):
        a, b = read.get(pair), used.get(pair)
        if a != b:
            what = "missing" if b is None else "extra" if a is None else f"{b}, read as {a}"
            bond_problems.append(f"bond {label[pair[0]]}-{label[pair[1]]}: {what}")
    if perceived.n_atoms != built.n_atoms:
        problems.append(f"{built.n_atoms} atoms built, {perceived.n_atoms} read")
    for kind, found in (
        ("elements", elements),
        ("formal charges", charges),
        ("bonds", bond_problems),
    ):
        if found:
            problems.append(f"{len(found)} {kind} differ (first: {'; '.join(found[:2])})")
    return problems


def site_rule(
    protein: Any, keys: list[tuple], offset: int = 0
) -> tuple[set[int], set[int], set[int]]:
    """The atoms the site rule puts in a site, read from the modification records alone.

    A site atom is a side-chain atom (CB outward) of a residue an attachment,
    a non-canonical mutation or a non-variant protonation change anchors; a
    backbone atom a change touched, with the atoms :data:`TOUCHED_BACKBONE`
    names (the residue's side chain then stays out); an atom a change added;
    or an atom of a residue an attachment brought. Returns those atoms, the
    backbone atoms among them, and the added or attached atoms, in realized
    indices.
    """
    from polyzymd.modifying_proteins.sites import CANONICAL, canonical_name

    residue_atoms: dict[tuple, list[int]] = collections.defaultdict(list)
    index_of = {}
    for i, key in enumerate(keys):
        residue_atoms[tuple(key[:3])].append(i)
        index_of[tuple(key)] = i

    def side_chain(residue: tuple) -> set[int]:
        return {i + offset for i in residue_atoms.get(residue, ()) if keys[i][3] not in BACKBONE}

    def touched(residue: tuple, atom: str) -> set[int]:
        names = TOUCHED_BACKBONE.get(atom, (atom,))
        return {i + offset for i in residue_atoms.get(residue, ()) if keys[i][3] in names}

    rule: set[int] = set()
    designed: set[int] = set()
    brought: set[int] = set()
    attached: set[tuple] = set()
    for record in protein.records:
        residue = (record.site.chain_id, record.site.resnum, record.site.icode)
        added = {index_of[tuple(a)] + offset for a in record.added if tuple(a) in index_of}
        atom = record.arguments.get("atom_name")
        brought |= added
        if record.kind == "mutate":
            to = record.arguments.get("to")
            if not (isinstance(to, str) and canonical_name(record.residue_after) in CANONICAL):
                rule |= side_chain(residue)
            continue
        if residue in attached:
            rule |= added
        elif atom in BACKBONE:
            backbone = touched(residue, atom)
            rule |= backbone | added
            designed |= backbone
        elif record.kind == "attach" or (record.residue_before, atom) not in FF14SB_VARIANTS:
            rule |= side_chain(residue) | added
        for chain_id, resnum, icode, _ in record.new_residues:
            attached.add((chain_id, resnum, icode))
            atoms = {i + offset for i in residue_atoms.get((chain_id, resnum, icode), ())}
            rule |= atoms
            brought |= atoms
    return rule, designed, brought


def declared_differences(reference: ModifiedReference, component: Any, declared: dict) -> list[str]:
    """How the declared modification blocks and parameter sets differ from what was built.

    ``declared`` comes from the caller (the config), independent of the build's
    own records: ``{"entries": [...], "parameter_sets": {type: {...}}, "base_dir":
    path}``, where ``entries`` is the config's own expansion, one entry per site
    (``polyzymd.config.modifications.declared_modifications``). A ``"spec"`` of
    blocks is still read, expanded here.
    """
    base = Path(declared.get("base_dir", "."))
    if "entries" in declared:
        wanted = [_declared_entry(item, base) for item in declared["entries"]]
    else:
        from polyzymd.modifying_proteins.spec import expand

        mutations, modifications = expand(declared["spec"], base)
        wanted = [_declared_mutation(entry, base) for entry in mutations]
        for entry in modifications:
            ((verb, arguments),) = entry.items()
            wanted.append(_declared_change(verb, arguments, base))
    built = [_built_change(record) for record in reference.records]
    problems = []
    for k in range(max(len(wanted), len(built))):
        w = wanted[k] if k < len(wanted) else None
        b = built[k] if k < len(built) else None
        if w != b:
            problems.append(f"change {k + 1}: declared {_show(w)}, built {_show(b)}")
    sets = declared.get("parameter_sets", {})
    for domain in component.domains[1:]:
        site_type = (domain.reference or {}).get("site_type")
        declared_set = sets.get(site_type)
        built_set = {
            "force_field": list(domain.force_field),
            "charges": {"method": domain.charges.method}
            | ({"model": domain.charges.model} if domain.charges.model else {}),
        }
        if declared_set is None or (
            list(declared_set.get("force_field", [])) != built_set["force_field"]
            or dict(declared_set.get("charges", {})) != built_set["charges"]
        ):
            problems.append(
                f"site {domain.name} ({site_type}): declared parameters {declared_set}, "
                f"built with {built_set}"
            )
    return problems


def _declared_entry(item: dict, base: Path) -> tuple:
    """One expanded declared modification, as :func:`_built_change` describes a record."""
    if item["kind"] == "mutate":
        return _declared_mutation(item, base)
    return _declared_change(item["kind"], item, base)


def _declared_mutation(entry: dict, base: Path) -> tuple:
    return ("mutate", _parse(entry["site"]), _target(entry["to"], base), entry.get("stereo"))


def _declared_change(verb: str, arguments: dict, base: Path) -> tuple:
    if verb != "attach":
        return (verb, _parse(arguments["site"]), arguments["atom_name"])
    return (
        "attach",
        _parse(arguments["site"]),
        arguments.get("type"),
        arguments.get("atom_name"),
        arguments.get("reaction"),
        _fragment(arguments.get("fragment") or {}, base),
        *_attachment_details(arguments),
    )


def _attachment_details(arguments: dict) -> tuple:
    """The fields of an attachment that decide which atoms bond and leave."""

    def names(value: Any) -> Any:
        return tuple([value] if isinstance(value, str) else value) if value else None

    return (
        arguments.get("fragment_atom_name"),
        arguments.get("fragment_resnum"),
        arguments.get("fragment_resname"),
        names(arguments.get("leaving_atom_names")),
        names(arguments.get("fragment_leaving_atom_names")),
        arguments.get("attach_at"),
        bool(arguments.get("merge", False)),
        int(arguments.get("bond_order", 1)),
    )


def _built_change(record: Any) -> tuple:
    site = (record.site.chain_id, record.site.resnum, record.site.icode)
    arguments = record.arguments
    if record.kind == "mutate":
        return ("mutate", site, _target(arguments["to"], Path(".")), arguments.get("stereo"))
    if record.kind != "attach":
        return (record.kind, site, arguments["atom_name"])
    return (
        "attach",
        site,
        record.type,
        arguments.get("atom_name"),
        arguments.get("reaction"),
        _fragment(arguments.get("fragment") or {}, Path(".")),
        *_attachment_details(arguments),
    )


def _target(to: Any, base: Path) -> Any:
    if isinstance(to, str):
        return to.upper()
    if "ccd_file" in to:
        return ("ccd_file", to.get("sha256") or _sha256(base / to["ccd_file"]))
    return ("smiles", to.get("smiles"), to.get("resname"))


#: What identifies a grafted chain: its monomers, caps, length, sequence block and seed.
#: Its drawn sequence and replicate are left out, since each replicate draws anew.
_CHAIN_KEYS = ("monomers", "end_groups", "length", "sequence", "random_seed")


def _fragment(fragment: dict, base: Path) -> tuple:
    if "polymer" in fragment:
        import json

        chain = {key: fragment["polymer"].get(key) for key in _CHAIN_KEYS}
        return ("polymer", json.dumps(chain, sort_keys=True, default=str))
    if "pdb_path" in fragment:
        return ("pdb", fragment.get("sha256") or _sha256(base / fragment["pdb_path"]))
    if "ccd" in fragment:
        return ("ccd", fragment["ccd"], fragment.get("link_atom"))
    if "ccd_file" in fragment:
        digest = fragment.get("sha256") or _sha256(base / fragment["ccd_file"])
        return ("ccd_file", digest, fragment.get("link_atom"))
    return ("smiles", fragment.get("smiles"), fragment.get("resname"))


def _sha256(path: Path) -> str:
    import hashlib

    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _show(change: Optional[tuple]) -> str:
    return "nothing" if change is None else repr(change)


def modified_termini(pdb: str, keys: list[tuple], domain_of: dict, offset: int = 0) -> dict:
    """Terminal residues of the input whose terminus a site holds.

    A residue's N-terminus is modified when the input residue is N-terminal
    (it has H1, H2 or H3) and a site holds its N; its C-terminus when the
    input residue is C-terminal (it has OXT) and a site holds its C, O, OXT
    or HXT. The input PDB is read here, line by line.
    Returns ``{residue: "N" or "C"}``.
    """
    names: dict[tuple, set[str]] = collections.defaultdict(set)
    for line in Path(pdb).read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            residue = (line[21].strip(), int(line[22:26]), line[26].strip())
            names[residue].add(line[12:16].strip())
    termini: dict[tuple, str] = {}
    for i, key in enumerate(keys):
        if domain_of.get(i + offset, "protein") == "protein":
            continue
        residue, name = tuple(key[:3]), key[3]
        present = names.get(residue, set())
        if name == "N" and present & {"H1", "H2", "H3"}:
            termini[residue] = "N"
        elif name in ("C", "O", "OXT", "HXT") and "OXT" in present:
            termini[residue] = "C"
    return termini


#: Anchor atoms of an attachment that labels each terminus (a cap, a label).
TERMINAL_ANCHORS = {"N": frozenset({"N"}), "C": frozenset({"C", "O", "OXT"})}


def neutral_termini(termini: dict, records: list) -> list[str]:
    """Modified termini that no attachment labels: a terminus changed by protonation alone.

    These are refused because ff14SB has no templates for a neutral
    N-terminus (NH2) or C-terminus (COOH). If you need one, open an issue at
    https://github.com/joelaforet/polyzymd/issues
    """
    found = []
    for residue, terminus in termini.items():
        labelled = any(
            record.kind == "attach"
            and (record.site.chain_id, record.site.resnum, record.site.icode) == residue
            and record.arguments.get("atom_name") in TERMINAL_ANCHORS[terminus]
            for record in records
        )
        if not labelled:
            found.append(
                f"chain {residue[0]} residue {residue[1]}{residue[2]} ({terminus}-terminus)"
            )
    return found


#: The cap residue and the atom it bonds to, for each terminus.
_CAPS = {"N": ("ACE", "N"), "C": ("NME", "C")}


def uncapped_termini(reference: Any, termini: dict) -> list[str]:
    """Modified termini (``{residue: "N" or "C"}``) that ``reference`` does not cap.

    A terminus is capped when an atom of an ACE (N-terminus) or NME
    (C-terminus) residue is bonded to the residue's N or C.
    """
    found = []
    for residue, terminus in termini.items():
        cap, anchor = _CAPS[terminus]
        capped = any(
            _key(atom) == (*residue, anchor) and other.metadata["residue_name"] == cap
            for bond in reference.bonds
            for atom, other in ((bond.atom1, bond.atom2), (bond.atom2, bond.atom1))
        )
        if not capped:
            found.append(
                f"chain {residue[0]} residue {residue[1]}{residue[2]} ({terminus}-terminus)"
            )
    return found


def _input_formal_charges(pdb: str) -> Optional[dict[tuple, float]]:
    """Formal charge of each residue of the input PDB, read by OpenFF's own PDB reader."""
    try:
        from openff.toolkit import Topology

        topology = Topology.from_pdb(str(pdb))
    except Exception:  # a structure OpenFF cannot read; the check then says so
        return None
    out: dict[tuple, float] = collections.defaultdict(float)
    for atom in topology.atoms:
        meta = atom.metadata
        residue = (
            meta["chain_id"],
            int(meta["residue_number"]),
            str(meta.get("insertion_code", "")).strip(),
        )
        out[residue] += atom.formal_charge.m_as("elementary_charge")
    return out


def rebuild(recipe: dict) -> Any:
    """Rebuild the modified protein from the input PDB and the modification blocks."""
    import hashlib

    from polyzymd.modifying_proteins.api import load_protein
    from polyzymd.modifying_proteins.spec import apply_spec

    path = Path(recipe["pdb"])
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if recipe.get("sha256") and digest != recipe["sha256"]:
        raise ValueError(f"{path} changed since the build (its SHA-256 differs)")
    protein = load_protein(path, ccd_paths=recipe.get("ccd_paths", []), download=recipe["download"])
    apply_spec(protein, recipe["spec"], base_dir=recipe["base_dir"], relax="none")
    return protein


def capped_model_charges(
    molecule: Any, site: list[int], residues: list[int], parameter_set: dict
) -> dict[int, float]:
    """NAGL charges of the ``site`` atoms of ``molecule``, computed on a capped model.

    The model holds the site, the rest of its residues (``residues``) and every
    atom within two bonds of them; each bond from the model to an atom outside
    it is replaced by a bond to a hydrogen.
    """
    core = sorted(set(site) | set(residues))
    from openff.toolkit import Molecule
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper

    neighbours: dict[int, set[int]] = collections.defaultdict(set)
    orders: dict[tuple[int, int], int] = {}
    for bond in molecule.bonds:
        a, b = bond.atom1_index, bond.atom2_index
        neighbours[a].add(b)
        neighbours[b].add(a)
        orders[(min(a, b), max(a, b))] = bond.bond_order
    kept = set(core)
    frontier = set(core)
    for _ in range(2):
        frontier = {j for i in frontier for j in neighbours[i]} - kept
        kept |= frontier
    model = Molecule()
    local: dict[int, int] = {}
    for i in sorted(kept):
        atom = molecule.atom(i)
        local[i] = model.add_atom(
            atom.atomic_number, atom.formal_charge, atom.is_aromatic, stereochemistry=None
        )
    for (a, b), order in orders.items():
        if a in kept and b in kept:
            aromatic = molecule.get_bond_between(a, b).is_aromatic
            model.add_bond(local[a], local[b], order, aromatic)
    for i in sorted(kept):
        for j in neighbours[i] - kept:
            if orders[(min(i, j), max(i, j))] != 1:
                raise ValueError(
                    f"the capped model of a site would cut a bond of order "
                    f"{orders[(min(i, j), max(i, j))]} between atoms {i} and {j}"
                )
            hydrogen = model.add_atom(1, 0, False)
            model.add_bond(local[i], hydrogen, 1, False)
    NAGLToolkitWrapper().assign_partial_charges(
        model, partial_charge_method=parameter_set["charges"]["model"]
    )
    values = model.partial_charges.m_as("elementary_charge")
    return {i: float(values[local[i]]) for i in site}


def protonation_warnings(protein: Any) -> list[str]:
    """Side chains a mutation left in a state other than the usual one at pH 7."""
    keys = protein.atom_keys()
    atoms_of: dict[tuple, set[str]] = collections.defaultdict(set)
    for key in keys:
        atoms_of[(key.chain_id, key.resnum, key.icode)].add(key.atom_name)
    warnings = []
    for record in protein.records:
        if record.kind != "mutate":
            continue
        residue = (record.site.chain_id, record.site.resnum, record.site.icode)
        names = atoms_of.get(residue, set())
        site = f"{residue[0]}:{record.residue_after}{residue[1]}{residue[2]}"
        if record.residue_after in _UNUSUAL:
            marker, verb, atom = _UNUSUAL[record.residue_after]
            if marker in names:
                warnings.append(
                    f"{site} is the neutral acid (as its CCD component defines it); for the "
                    f"charged form add `- {verb}: {{site: {site}, atom_name: {atom}}}`"
                )
        elif record.residue_after == "HIS" and {"HD1", "HE2"} <= names:
            warnings.append(
                f"{site} is doubly protonated (HIP, charged, as its CCD component defines it); "
                f"for a neutral tautomer add `- deprotonate: {{site: {site}, atom_name: NE2}}` "
                "(HID) or ND1 (HIE)"
            )
    return warnings


def compare_domains(realized: Any, component: Any, ref: ModifiedReference) -> dict[str, dict]:
    """Compare each domain's owned terms of ``realized`` with the expected ones.

    Returns, per domain, a :class:`~polyzymd.audit.terms.Comparison` per kind.
    """
    index = realized.index
    component_atoms = list(ref.domain_of)
    results = {}
    for domain in [d.name for d in component.domains]:
        atoms = [i for i, d in ref.domain_of.items() if d == domain]
        actual = index.extract(
            atoms,
            owned=lambda term, domain=domain: ref.owner(term) == domain,
            terms_of=component_atoms,
        )
        results[domain] = {
            kind: compare(kind, actual[kind], ref.expected[domain][kind]) for kind in TERM_KINDS
        }
    return results


def _owned_terms(
    expected: dict,
    index: SystemIndex,
    to_real: dict[int, int],
    domain_of: dict[int, str],
    domain: str,
) -> None:
    """Add to ``expected[domain]`` the valence terms of ``index`` that ``domain`` owns."""

    def owned(term: tuple[int, ...]) -> bool:
        if not all(i in to_real for i in term):
            return False
        domains = {domain_of.get(to_real[i]) for i in term}
        if None in domains:
            return False
        sites = domains - {"protein"}
        if len(sites) > 1:
            return False
        return (sites.pop() if sites else "protein") == domain

    atoms = list(to_real)
    terms = index.extract(atoms, to_real, owned=owned)
    for kind in _VALENCE:
        expected[domain][kind].update(terms[kind])


def _expected_pairs(ref: ModifiedReference, per_atom, charges, scales, native) -> None:
    """Add every expected 1-2, 1-3 (excluded) and 1-4 pair to its owner's terms."""
    neighbours: dict[int, set[int]] = collections.defaultdict(set)
    for a, b in ref.bonds:
        neighbours[a].add(b)
        neighbours[b].add(a)
    for start in ref.domain_of:
        distance = {start: 0}
        frontier = [start]
        for depth in (1, 2, 3):
            frontier = [j for i in frontier for j in neighbours[i] if j not in distance]
            for j in frontier:
                distance.setdefault(j, depth)
        for j, d in distance.items():
            if j <= start or d == 0:
                continue
            pair = (start, j)
            owner = ref.owner(pair)
            if owner is None:
                continue
            if d < 3:
                value = (0.0, 0.0, 0.0)
            elif native is not None and owner in native["sites"]:
                value = native["exceptions"].get(pair)
                if value is None:
                    continue
            else:
                coulomb, lj = scales[owner]
                qi, (si, ei), _ = per_atom[start]
                qj, (sj, ej), _ = per_atom[j]
                epsilon = lj * float(np.sqrt(ei * ej))
                sigma = (si + sj) / 2 if epsilon != 0 else 0.0
                value = (coulomb * charges[start] * charges[j], sigma, epsilon)
            ref.expected[owner]["exceptions"][pair] = value


def _native_glycam(recipe: dict, protein: Any, glycam_sites, keys, real) -> dict:
    """OpenMM's ff14SB + GLYCAM06j on the protein with only its glycans, mapped to the build."""
    import openmm
    from openmm import app

    from polyzymd.modifying_proteins import glycam
    from polyzymd.modifying_proteins.api import load_protein
    from polyzymd.modifying_proteins.parameters import to_openff
    from polyzymd.modifying_proteins.spec import apply_spec
    from polyzymd.modifying_proteins.topology import modification_sites

    glycans_only = _glycans_only(recipe["spec"], protein, glycam_sites)
    replay = load_protein(
        Path(recipe["pdb"]), ccd_paths=recipe.get("ccd_paths", []), download=recipe["download"]
    )
    apply_spec(replay, glycans_only, base_dir=recipe["base_dir"], relax="none")
    replay_keys = [tuple(key) for key in replay.atom_keys()]
    replay_sites = {site.residues[0]: site for site in modification_sites(replay)}

    # Glycan residues are numbered after the chain, so they map through their site.
    residue_map: dict[tuple, tuple] = {}
    for site in glycam_sites:
        anchor = site.residues[0]
        built = sorted({keys[i][:3] for i in site.atom_indices} - {anchor}, key=lambda r: r[1])
        replayed_site = replay_sites.get(anchor)
        if replayed_site is None:
            raise ValueError(f"replaying the glycans alone lost site {site.name}")
        replayed = sorted(
            {replay_keys[i][:3] for i in replayed_site.atom_indices} - {anchor},
            key=lambda r: r[1],
        )
        if len(built) != len(replayed):
            raise ValueError(f"replaying the glycans alone changed site {site.name}")
        residue_map.update(zip(replayed, built))
    position = {key: i for i, key in enumerate(keys)}
    to_real = {}
    for k, key in enumerate(replay_keys):
        mapped = (*residue_map.get(key[:3], key[:3]), key[3])
        if mapped in position:
            to_real[k] = real(position[mapped])

    molecule = to_openff(replay)
    topology = molecule.to_topology().to_openmm()
    names = glycam.glycam_residue_names()
    residue_name = {
        key[:3]: molecule.atom(k).metadata["residue_name"] for k, key in enumerate(replay_keys)
    }
    linkages = {
        site.residues[0]: glycam.LINKAGE_TEMPLATES.get(residue_name[site.residues[0]])
        for site in modification_sites(replay)
        if site.type == "glycan"
    }
    templates = {}
    for residue, first in zip(topology.residues(), _first_atoms(topology)):
        meta = molecule.atom(first).metadata
        key = (meta["chain_id"], int(meta["residue_number"]), meta.get("insertion_code", ""))
        if key in linkages and linkages[key]:
            residue.name = linkages[key]
            templates[residue] = residue.name
        elif residue.name in names:
            templates[residue] = residue.name
        elif residue.name in _PROTEIN_CAPS:
            # ff14SB's ACE/NME also match GLYCAM's ACX/NMX; the protein's caps are ff14SB's.
            templates[residue] = residue.name
    system = app.ForceField(*glycam.OPENMM_FILES).createSystem(
        topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
        removeCMMotion=False,
        residueTemplates=templates,
    )
    index = SystemIndex.from_system(system)
    nonbonded = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
    exceptions = {}
    for n in range(nonbonded.getNumExceptions()):
        i, j, product, sigma, epsilon = nonbonded.getExceptionParameters(n)
        if i in to_real and j in to_real:
            pair = tuple(sorted((to_real[i], to_real[j])))
            eps = epsilon._value
            exceptions[pair] = (product._value, sigma._value if eps != 0 else 0.0, eps)
    return {
        "index": index,
        "to_real": to_real,
        "particle": {to_real[k]: index.particle[k] for k in to_real},
        "exceptions": exceptions,
        "sites": {site.name for site in glycam_sites},
    }


def _glycans_only(spec: dict, protein: Any, glycam_sites) -> dict:
    """The modification blocks with only the glycans and what keeps residues ff14SB's."""
    from polyzymd.modifying_proteins.sites import CANONICAL

    glycan_anchors = {site.residues[0] for site in glycam_sites}
    other_anchors = set()
    for record in protein.records:
        if record.kind == "attach" and record.type != "glycan":
            other_anchors.add((record.site.chain_id, record.site.resnum, record.site.icode))
    mutations = [
        entry
        for entry in spec.get("mutations", [])
        if isinstance(entry.get("to"), str) and entry["to"].upper() in CANONICAL
    ]
    modifications = []
    for entry in spec.get("modifications", []):
        ((verb, arguments),) = entry.items()
        site = _parse(arguments["site"])
        if verb == "attach":
            if arguments.get("type") == "glycan" and site in glycan_anchors:
                modifications.append(entry)
        elif site not in other_anchors:
            modifications.append(entry)
    out = {}
    if mutations:
        out["mutations"] = mutations
    out["modifications"] = modifications
    return out


def _parse(text: str) -> tuple:
    from polyzymd.modifying_proteins.sites import parse_site

    site = parse_site(text)
    return (site.chain_id, site.resnum, site.icode)


def _first_atoms(topology: Any) -> list[int]:
    return [next(iter(residue.atoms())).index for residue in topology.residues()]


def _site_set(protein: Any, site: Any, sets: dict) -> dict:
    """The parameter set a site's type names (all its merged types must agree)."""
    types = {site.type}
    for record in protein.records:
        if record.kind == "attach" and set(record.added) & set(site.atoms):
            types.add(record.type)
    chosen = [sets.get(t) for t in sorted(types)]
    if any(s is None for s in chosen):
        raise ValueError(f"site {site.name}: no parameter set for {', '.join(sorted(types))}")
    if any(s != chosen[0] for s in chosen):
        raise ValueError(f"site {site.name}: its types {sorted(types)} have different sets")
    return chosen[0]


def _is_glycam(parameter_set: dict) -> bool:
    return any("glycam" in str(name).lower() for name in parameter_set["force_field"])


def _scales(force_field) -> tuple[float, float]:
    """The (Coulomb, Lennard-Jones) 1-4 scales of an OpenFF force field."""
    from openff.toolkit import ForceField

    loaded = ForceField(*force_field)
    return float(loaded["Electrostatics"].scale14), float(loaded["vdW"].scale14)


def _key(atom: Any) -> tuple:
    meta = atom.metadata
    return (
        meta["chain_id"],
        int(meta["residue_number"]),
        meta.get("insertion_code", ""),
        atom.name,
    )


def _label(key: tuple) -> str:
    return f"{key[0]}:{key[1]}{key[2]}:{key[3]}"


def _first_difference(a: list, b: list) -> int:
    return next((k for k, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))
