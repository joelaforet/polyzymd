"""Rebuild each domain of a component from its inputs alone, and compare it with the build.

A domain's reference is made with OpenFF directly from the force field files and
charge source the component record names, and from the inputs the record points
to (the protein PDB, an SDF, a SMILES). Every term of the realized system inside
the domain is then compared with the reference, term by term.

Copies of one molecule (thousands of waters, dozens of chains of one sequence)
are compared once: molecules whose terms are identical, in their own atom order,
form one group, and one member of each group is compared with the reference.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Any, Optional

from polyzymd.audit.molecules import (
    graph_hash,
    map_molecule,
    molecule_graph,
    reference_system,
    system_graph,
)
from polyzymd.audit.realized import RealizedSystem
from polyzymd.audit.terms import TERM_KINDS, compare

#: The Constraints section lent to force fields without one (the Amber ports):
#: bonds to hydrogen at the force field's own bond lengths, and rigid water.
CONSTRAINTS_FROM = "openff-2.3.0.offxml"


@dataclass
class KindResult:
    """Agreement of one kind of term over a domain."""

    matched: int = 0
    total: int = 0
    examples: list[str] = field(default_factory=list)

    def add(self, matched: int, total: int, examples: list, copies: int = 1) -> None:
        self.matched += matched * copies
        self.total += total * copies
        for example in examples:
            if len(self.examples) < 3:
                self.examples.append(repr(example))

    @property
    def ok(self) -> bool:
        return self.matched == self.total


@dataclass
class DomainComparison:
    """How the realized terms of one domain compare with its reference."""

    component: str
    domain: str
    reference: str
    n_molecules: int = 0
    kinds: dict[str, KindResult] = field(
        default_factory=lambda: {k: KindResult() for k in TERM_KINDS}
    )
    #: Problems that stopped the comparison (no reference, molecules not the reference).
    errors: list[str] = field(default_factory=list)
    #: False when the charge method cannot be recomputed here (am1bcc, espaloma):
    #: the reference then carries the build's charges, and only the other terms are checked.
    charges_recomputed: bool = True


def compare_component(realized: RealizedSystem, component: Any) -> list[DomainComparison]:
    """Compare every domain of ``component`` with its reference."""
    return [_compare_domain(realized, component, domain) for domain in component.domains]


def _compare_domain(realized: RealizedSystem, component: Any, domain: Any) -> DomainComparison:
    ref = domain.reference or {}
    label = f"{' + '.join(domain.force_field)}, charges {_charge_label(domain.charges)}"
    result = DomainComparison(component.name, domain.name, label)
    try:
        if "pdb" in ref:
            _compare_whole(realized, domain, ref, result)
        elif "sdf" in ref or "smiles" in ref:
            _compare_copies(realized, domain, ref, result)
        else:
            result.errors.append("the component record says nothing to rebuild this domain from")
    except Exception as error:  # OpenFF raises many types; any of them means no reference
        result.errors.append(f"could not build the reference: {type(error).__name__}: {error}")
    return result


def _charge_label(charges: Any) -> str:
    if charges.method == "nagl":
        return f"nagl {charges.model}"
    if charges.method == "preset" and charges.source_file:
        return f"preset from {charges.source_file}"
    return charges.method


def _compare_whole(
    realized: RealizedSystem, domain: Any, ref: dict, result: DomainComparison
) -> None:
    """A domain read from a PDB (the protein): its atoms in order are the reference's atoms."""
    from openff.toolkit import Topology

    templates = None
    if ref.get("custom_substructures"):
        import json
        from pathlib import Path

        templates = json.loads(Path(ref["custom_substructures"]).read_text())
    if templates:
        topology = Topology.from_pdb(ref["pdb"], _custom_substructures=templates)
    else:
        topology = Topology.from_pdb(ref["pdb"])
    atoms = [i for r in domain.atoms for i in r]
    result.n_molecules = len(domain.atoms)
    elements = realized.elements
    expected = [atom.symbol for atom in topology.atoms]
    if len(expected) != len(atoms) or [elements[i] for i in atoms] != expected:
        result.errors.append(
            f"its {len(atoms)} atoms are not the {len(expected)} atoms of {ref['pdb']}, in order"
        )
        return
    system = reference_system(domain.force_field, topology, constraints_from=CONSTRAINTS_FROM)
    expected_terms = _all_terms(system)
    actual = realized.index.extract(atoms, {s: k for k, s in enumerate(atoms)})
    for kind in TERM_KINDS:
        comparison = compare(kind, actual[kind], expected_terms[kind])
        result.kinds[kind].add(comparison.matched, comparison.total, comparison.mismatches)


def _compare_copies(
    realized: RealizedSystem, domain: Any, ref: dict, result: DomainComparison
) -> None:
    """A domain of copies of one molecule: compare one copy of each distinct set of terms."""
    molecule = _reference_molecule(ref)
    graph = molecule_graph(molecule)
    expected_hash = graph_hash(graph)
    result.n_molecules = len(domain.atoms)
    groups: dict[str, list[Any]] = {}
    index = realized.index
    for atom_range in domain.atoms:
        local = index.extract(atom_range, {i: i - atom_range.start for i in atom_range})
        key = hashlib.sha256(
            repr(sorted((k, sorted(v.items())) for k, v in local.items())).encode()
        )
        groups.setdefault(key.hexdigest(), []).append(atom_range)
    charged = _charged(molecule, domain.charges, result)
    expected_terms = None
    if charged is not None:
        expected_terms = _all_terms(_molecule_system(domain, charged))
    elements = realized.elements
    for members in groups.values():
        atoms = list(members[0])
        sub_graph = system_graph(index, atoms, elements)
        mapping = map_molecule(molecule, sub_graph, atoms, identity_only=True)
        if mapping is None and graph_hash(sub_graph) == expected_hash:
            mapping = map_molecule(molecule, sub_graph, atoms)
        if mapping is None:
            result.errors.append(
                f"{len(members)} molecule(s) starting at atom {atoms[0]} are not the reference "
                f"molecule ({molecule.to_smiles(explicit_hydrogens=False)})"
            )
            continue
        actual = index.extract(atoms, mapping)
        terms = expected_terms
        if terms is None:
            # The charges cannot be recomputed here: give the reference the build's charges.
            from openff.toolkit import Molecule
            from openff.units import unit

            copy = Molecule(molecule)
            copy.partial_charges = [actual["charge"][k] for k in range(copy.n_atoms)] * unit.e
            terms = _all_terms(_molecule_system(domain, copy))
        for kind in TERM_KINDS:
            comparison = compare(kind, actual[kind], terms[kind])
            result.kinds[kind].add(
                comparison.matched, comparison.total, comparison.mismatches, copies=len(members)
            )


def _reference_molecule(ref: dict) -> Any:
    from openff.toolkit import Molecule

    if "smiles" in ref:
        return Molecule.from_smiles(ref["smiles"])
    loaded = Molecule.from_file(ref["sdf"], allow_undefined_stereo=True)
    molecules = loaded if isinstance(loaded, list) else [loaded]
    return molecules[ref.get("conformer_index", 0)]


def _charged(molecule: Any, charges: Any, result: DomainComparison) -> Optional[Any]:
    """Return the reference molecule charged by the declared source, or None for library
    charges (the force field assigns them), or when the method cannot be recomputed."""
    from openff.toolkit import Molecule

    method = charges.method
    if method in ("library", "water_model"):
        return _uncharged(molecule)
    if method == "preset":
        if molecule.partial_charges is None:
            raise ValueError("charges preset, but the reference file holds no partial charges")
        return Molecule(molecule)
    if method == "nagl":
        from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper

        copy = Molecule(molecule)
        NAGLToolkitWrapper().assign_partial_charges(copy, partial_charge_method=charges.model)
        return copy
    result.charges_recomputed = False
    return None


def _uncharged(molecule: Any) -> Any:
    from openff.toolkit import Molecule

    copy = Molecule(molecule)
    copy.partial_charges = None
    return copy


def _molecule_system(domain: Any, molecule: Any) -> Any:
    charge_from = [molecule] if molecule.partial_charges is not None else []
    return reference_system(
        domain.force_field, molecule.to_topology(), charge_from, constraints_from=CONSTRAINTS_FROM
    )


def _all_terms(system: Any) -> dict[str, dict]:
    from polyzymd.audit.terms import SystemIndex

    return SystemIndex.from_system(system).extract(range(system.getNumParticles()))
