"""Match atoms of a built system to reference molecules, and parameterize references.

References are built with OpenFF directly, from one force field at a time, so
the audit never shares a code path with the builder it checks.
"""

from __future__ import annotations

from typing import Any, Iterable, Optional, Sequence

import networkx as nx

from polyzymd.audit.terms import SystemIndex


def system_graph(index: SystemIndex, atoms: Iterable[int], elements: Sequence[str]) -> nx.Graph:
    """Return the bond graph of ``atoms``, from bonds and constraints, with elements on nodes.

    A constraint between two hydrogens (the H-H distance of rigid water) is not a bond.
    """
    atoms = set(atoms)
    graph = nx.Graph()
    for i in atoms:
        graph.add_node(i, element=elements[i])
    for i in atoms:
        for j in index.bonded_to.get(i, ()):
            if j in atoms and not (elements[i] == "H" and elements[j] == "H"):
                graph.add_edge(i, j)
    return graph


def molecule_graph(molecule: Any) -> nx.Graph:
    """Return the bond graph of an OpenFF Molecule, with elements on nodes."""
    graph = nx.Graph()
    for i, atom in enumerate(molecule.atoms):
        graph.add_node(i, element=atom.symbol)
    for bond in molecule.bonds:
        graph.add_edge(bond.atom1_index, bond.atom2_index)
    return graph


def graph_hash(graph: nx.Graph) -> str:
    """Return a hash two graphs share when they could be isomorphic."""
    return nx.weisfeiler_lehman_graph_hash(graph, node_attr="element")


def map_molecule(
    molecule: Any, sub_graph: nx.Graph, atoms: Sequence[int], identity_only: bool = False
) -> Optional[dict[int, int]]:
    """Return {system index: molecule index}, or None if ``atoms`` are not ``molecule``.

    PolyzyMD keeps the atom order of its inputs, so the order-preserving map is
    tried first; otherwise an element-preserving graph isomorphism is searched,
    unless ``identity_only``.
    """
    reference = molecule_graph(molecule)
    if len(atoms) != reference.number_of_nodes():
        return None
    ordered = sorted(atoms)
    identity = {s: r for r, s in enumerate(ordered)}
    if all(
        sub_graph.nodes[s]["element"] == reference.nodes[r]["element"] for s, r in identity.items()
    ):
        if {tuple(sorted((identity[a], identity[b]))) for a, b in sub_graph.edges} == {
            tuple(sorted(edge)) for edge in reference.edges
        }:
            return identity
    if identity_only:
        return None
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        sub_graph, reference, node_match=lambda a, b: a["element"] == b["element"]
    )
    return next(matcher.isomorphisms_iter(), None)


def molecule_maps(molecule: Any, sub_graph: nx.Graph, limit: int = 64) -> list[dict[int, int]]:
    """Return up to ``limit`` element-preserving maps of ``sub_graph`` onto ``molecule``.

    A graph without bond orders has symmetries the chemistry does not (the two
    oxygens of a carboxylate, say), so a caller picks among these the map whose
    terms agree best, instead of trusting the first.
    """
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        sub_graph, molecule_graph(molecule), node_match=lambda a, b: a["element"] == b["element"]
    )
    maps = []
    for mapping in matcher.isomorphisms_iter():
        maps.append(mapping)
        if len(maps) >= limit:
            break
    return maps


def reference_system(
    force_field_files: Sequence[str],
    topology: Any,
    charge_from: Sequence[Any] = (),
    constraints_from: Optional[str] = None,
) -> Any:
    """Return the OpenMM System OpenFF builds for ``topology`` from the named force fields only.

    ``constraints_from`` lends its Constraints section to a force field without one
    (the Amber ports), so bonds to hydrogen are constrained as in a PolyzyMD build.
    """
    from openff.toolkit import ForceField

    force_field = ForceField(*force_field_files)
    if constraints_from and "Constraints" not in force_field.registered_parameter_handlers:
        force_field.register_parameter_handler(ForceField(constraints_from)["Constraints"])
    interchange = force_field.create_interchange(topology, charge_from_molecules=list(charge_from))
    return interchange.to_openmm(combine_nonbonded_forces=True)
