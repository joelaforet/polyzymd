"""Tests for offline CGSmiles and mBuild polymer authoring."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import nbformat
import networkx as nx
import pytest
import yaml

from polyzymd.builders.polymer import PolymerBuilder
from polyzymd.builders.polymer_authoring import (
    build_and_export_polymer,
    convert_charge_write_sdf,
    emit_provided_molecules_yaml_snippet,
    place_templates_on_cg_beads,
    resolve_cgsmiles_graph,
    validate_cgsmiles_graph,
    validate_graph_equivalence,
    validate_resolution,
)


def _methyl_template(name: str, port_vectors: list[tuple[float, float, float]]):
    """Create a tiny mBuild methane-like template with labelled ports."""
    import mbuild as mb

    template = mb.Compound(name=name)
    carbon = mb.Compound(name="C", element="C", pos=[0.0, 0.0, 0.0])
    template.add(carbon, label="C")
    for index, position in enumerate(
        ([0.11, 0.0, 0.0], [-0.04, 0.10, 0.0], [-0.04, -0.05, 0.09], [-0.04, -0.05, -0.09]),
        start=1,
    ):
        hydrogen = mb.Compound(name=f"H{index}", element="H", pos=position)
        template.add(hydrogen, label=f"H{index}")
        template.add_bond((carbon, hydrogen), bond_order=1)
    for index, vector in enumerate(port_vectors, start=1):
        template.add(
            mb.Port(anchor=carbon, orientation=vector, separation=0.08), label=f"port_{index}"
        )
    return template


def _star_graph() -> nx.Graph:
    """Return a three-arm coarse star graph with explicit cap nodes."""
    graph = nx.Graph()
    graph.add_node(0, fragname="CORE", position=(0.0, 0.0, 0.0))
    graph.add_node(1, fragname="ARM", position=(0.55, 0.0, 0.0))
    graph.add_node(2, fragname="ARM", position=(-0.28, 0.48, 0.0))
    graph.add_node(3, fragname="ARM", position=(-0.28, -0.48, 0.0))
    graph.add_edges_from([(0, 1, {"order": 1}), (0, 2, {"order": 1}), (0, 3, {"order": 1})])
    return graph


def _authoring_inputs() -> dict[str, object]:
    """Return compact authoring inputs for a tiny carbon star."""
    return {
        "fragments": "{#CORE=[$br]C([$br])[$br],#ARM=[$br]C}",
        "aa_templates": {
            "CORE": _methyl_template("CORE", [(1, 0, 0), (-0.5, 0.86, 0), (-0.5, -0.86, 0)]),
            "ARM": _methyl_template("ARM", [(-1, 0, 0)]),
        },
        "atom_maps": {"CORE": {0: "C"}, "ARM": {0: "C"}},
        "port_maps": {"CORE": {"$br1": ["port_1", "port_2", "port_3"]}, "ARM": {"$br1": "port_1"}},
        "leaving_maps": {"CORE": {"$br1": ["H1", "H2", "H3"]}, "ARM": {"$br1": "H1"}},
    }


def _single_edge_graph(order: float = 1.0, **edge_data: object) -> nx.Graph:
    """Return a two-node coarse graph with one edge."""
    graph = nx.Graph()
    graph.add_node(0, fragname="A", position=(0.0, 0.0, 0.0))
    graph.add_node(1, fragname="A", position=(0.3, 0.0, 0.0))
    graph.add_edge(0, 1, order=order, **edge_data)
    return graph


def _single_interfragment_atom_graph(order: float = 1.0, **edge_data: object) -> nx.Graph:
    """Return a resolved atom graph with one interfragment bond."""
    atom_graph = nx.Graph()
    atom_graph.add_node(0, fragid=[0], mapping=[("A", 0)], element="C")
    atom_graph.add_node(1, fragid=[1], mapping=[("A", 0)], element="C")
    atom_graph.add_edge(0, 1, bonding=("$a", "$a"), order=order, **edge_data)
    return atom_graph


def test_validate_rejects_cycles() -> None:
    """Cycle graphs are outside the Phase15 tree-only contract."""
    graph = nx.Graph()
    for node in range(3):
        graph.add_node(node, fragname="A", position=(float(node), 0.0, 0.0))
    graph.add_edges_from([(0, 1, {"order": 1}), (1, 2, {"order": 1}), (2, 0, {"order": 1})])
    with pytest.raises(ValueError, match="tree graphs only"):
        validate_cgsmiles_graph(graph)


@pytest.mark.parametrize("order", [1.0, 1.5, 2.0, 3.0])
def test_supported_coarse_bond_orders_are_allowed(order: float) -> None:
    """Coarse edge order should mean bond order, not multiplicity."""
    graph = _single_edge_graph(order, aromatic=True) if order == 1.5 else _single_edge_graph(order)
    atom_graph = (
        _single_interfragment_atom_graph(order, aromatic=True)
        if order == 1.5
        else _single_interfragment_atom_graph(order)
    )
    validate_cgsmiles_graph(graph)
    validate_resolution(graph, atom_graph)


@pytest.mark.parametrize("order", [0.0, -1.0, 4.0, 1.25])
def test_unsupported_coarse_bond_orders_rejected(order: float) -> None:
    """Unsupported coarse edge bond orders should fail before authoring."""
    graph = _single_edge_graph(order)
    with pytest.raises(ValueError, match="Unsupported"):
        validate_cgsmiles_graph(graph)


def test_conflicting_aromatic_metadata_rejected() -> None:
    """Aromatic metadata must agree with the normalized 1.5 bond order."""
    with pytest.raises(ValueError, match="Conflicting aromatic metadata"):
        validate_cgsmiles_graph(_single_edge_graph(2.0, aromatic=True))
    with pytest.raises(ValueError, match="Conflicting aromatic metadata"):
        validate_resolution(
            _single_edge_graph(1.5, aromatic=True),
            _single_interfragment_atom_graph(1.5, aromatic=False),
        )


def test_duplicate_interfragment_records_rejected_with_pim1_guidance() -> None:
    """One coarse edge may resolve to exactly one atomistic interfragment bond."""
    graph = nx.Graph()
    graph.add_node(0, fragname="A", position=(0.0, 0.0, 0.0))
    graph.add_node(1, fragname="A", position=(0.3, 0.0, 0.0))
    graph.add_edge(0, 1, order=1)
    atom_graph = nx.Graph()
    for node, fragid in enumerate((0, 0, 1, 1)):
        atom_graph.add_node(node, fragid=[fragid], mapping=[("A", node % 2)], element="C")
    atom_graph.add_edge(0, 2, bonding=("$a", "$a"), order=1)
    atom_graph.add_edge(1, 3, bonding=("$b", "$b"), order=1)

    with pytest.raises(ValueError, match="PIM-1-like doubly linked or multi-point"):
        validate_resolution(graph, atom_graph)


def test_interfragment_order_mismatch_rejected() -> None:
    """Resolved atomistic bond order must match the corresponding CG edge order."""
    with pytest.raises(ValueError, match="bond order does not match CG edge order"):
        validate_resolution(_single_edge_graph(2.0), _single_interfragment_atom_graph(1.0))


def test_missing_maps_rejected_before_arbitrary_leaving_removal() -> None:
    """Users must supply explicit atom, port, and leaving atom maps."""
    graph = _star_graph()
    inputs = _authoring_inputs()
    _, atom_graph = resolve_cgsmiles_graph(inputs["fragments"], graph)
    with pytest.raises(ValueError, match="Missing explicit leaving atom map"):
        build_and_export_polymer(
            cg_graph=graph,
            fragments=inputs["fragments"],
            aa_templates=inputs["aa_templates"],
            atom_maps=inputs["atom_maps"],
            port_maps=inputs["port_maps"],
            leaving_maps={"CORE": {"$br": "H1"}},
            output_sdf=Path("unused.sdf"),
        )
    assert atom_graph.number_of_nodes() > 0


def test_missing_atom_map_rejected() -> None:
    """Explicit local atom maps are required for every bonded CGSmiles atom."""
    graph = _star_graph()
    inputs = _authoring_inputs()
    with pytest.raises(ValueError, match="Missing atom map"):
        build_and_export_polymer(
            cg_graph=graph,
            fragments=inputs["fragments"],
            aa_templates=inputs["aa_templates"],
            atom_maps={"CORE": {}, "ARM": {0: "C"}},
            port_maps=inputs["port_maps"],
            leaving_maps=inputs["leaving_maps"],
            output_sdf=Path("unused.sdf"),
        )


def test_duplicate_or_missing_ports_rejected() -> None:
    """Every mBuild Port must be available once for each consumed edge."""
    graph = _star_graph()
    inputs = _authoring_inputs()
    with pytest.raises(ValueError, match="No unused mBuild Port"):
        build_and_export_polymer(
            cg_graph=graph,
            fragments=inputs["fragments"],
            aa_templates=inputs["aa_templates"],
            atom_maps=inputs["atom_maps"],
            port_maps={"CORE": {"$br1": ["port_1"]}, "ARM": {"$br1": "port_1"}},
            leaving_maps=inputs["leaving_maps"],
            output_sdf=Path("unused.sdf"),
        )


def test_exact_label_required() -> None:
    """Port and leaving maps should not use prefix fallback labels."""
    graph = _star_graph()
    inputs = _authoring_inputs()
    with pytest.raises(ValueError, match="Missing port map"):
        build_and_export_polymer(
            cg_graph=graph,
            fragments=inputs["fragments"],
            aa_templates=inputs["aa_templates"],
            atom_maps=inputs["atom_maps"],
            port_maps={"CORE": {"$br": ["port_1", "port_2", "port_3"]}, "ARM": {"$br": "port_1"}},
            leaving_maps=inputs["leaving_maps"],
            output_sdf=Path("unused.sdf"),
        )


def test_leaving_alternatives_are_consumed_and_exhaustion_rejected() -> None:
    """Leaving sequence entries should be consumed once per occurrence."""
    graph = _star_graph()
    inputs = _authoring_inputs()
    with pytest.raises(ValueError, match="No explicit leaving atom remains"):
        build_and_export_polymer(
            cg_graph=graph,
            fragments=inputs["fragments"],
            aa_templates=inputs["aa_templates"],
            atom_maps=inputs["atom_maps"],
            port_maps=inputs["port_maps"],
            leaving_maps={"CORE": {"$br1": ["H1"]}, "ARM": {"$br1": "H1"}},
            output_sdf=Path("unused.sdf"),
        )


@pytest.mark.conjugation_stack
def test_element_mismatch_rejected(tmp_path: Path) -> None:
    """Graph equivalence should catch template element mismatches."""
    graph = _star_graph()
    inputs = _authoring_inputs()
    bad_templates = dict(inputs["aa_templates"])
    bad_templates["ARM"].labels["C"].element = "O"
    with pytest.raises(ValueError, match="not equivalent"):
        build_and_export_polymer(
            cg_graph=graph,
            fragments=inputs["fragments"],
            aa_templates=bad_templates,
            atom_maps=inputs["atom_maps"],
            port_maps=inputs["port_maps"],
            leaving_maps=inputs["leaving_maps"],
            output_sdf=tmp_path / "bad.sdf",
        )


def test_charge_aromatic_and_stereo_mismatch_rejected() -> None:
    """Graph equivalence should compare supplied charge and stereo metadata."""
    graph = _star_graph()
    inputs = _authoring_inputs()
    _, atom_graph = resolve_cgsmiles_graph(inputs["fragments"], graph)
    compound = place_templates_on_cg_beads(
        graph,
        inputs["aa_templates"],
        atom_graph,
        atom_maps=inputs["atom_maps"],
        port_maps=inputs["port_maps"],
        leaving_maps=inputs["leaving_maps"],
    )
    first_node = next(iter(atom_graph.nodes))
    atom_graph.nodes[first_node]["formal_charge"] = 1
    atom_graph.nodes[first_node]["aromatic"] = True
    atom_graph.nodes[first_node]["stereo"] = "R"

    with pytest.raises(ValueError, match="not equivalent"):
        validate_graph_equivalence(atom_graph, compound)


def test_force_overlap_bonds_preserve_resolved_order() -> None:
    """Created mBuild interfragment bonds should carry resolved bond orders."""
    graph = _star_graph()
    inputs = _authoring_inputs()
    _, atom_graph = resolve_cgsmiles_graph(inputs["fragments"], graph)
    for _, _, data in graph.edges(data=True):
        data["order"] = 2.0
    for atom_u, atom_v, data in atom_graph.edges(data=True):
        fragids_u = tuple(atom_graph.nodes[atom_u].get("fragid", ()))
        fragids_v = tuple(atom_graph.nodes[atom_v].get("fragid", ()))
        if fragids_u != fragids_v:
            data["order"] = 2.0
    compound = place_templates_on_cg_beads(
        graph,
        inputs["aa_templates"],
        atom_graph,
        atom_maps=inputs["atom_maps"],
        port_maps=inputs["port_maps"],
        leaving_maps=inputs["leaving_maps"],
    )

    created_orders = [
        data.get("order")
        for _, _, data in compound.bond_graph.edges(data=True)
        if data.get("order") == 2.0 and data.get("bond_order") == 2.0
    ]
    assert len(created_orders) == 3


def test_stereo_aliases_are_normalized_and_conflicts_rejected() -> None:
    """Node and edge stereo aliases should normalize to canonical stereo."""
    import mbuild as mb

    atom_graph = _single_interfragment_atom_graph(1.0, stereo="E", bond_stereo="E")
    for node in atom_graph.nodes:
        atom_graph.nodes[node]["stereo"] = "R"
        atom_graph.nodes[node]["chirality"] = "R"

    compound = mb.Compound(name="stereo_alias")
    atom_a = mb.Compound(name="C", element="C", pos=[0.0, 0.0, 0.0])
    atom_b = mb.Compound(name="C", element="C", pos=[0.2, 0.0, 0.0])
    atom_a.stereo = "R"
    atom_b.stereo = "R"
    compound.add(atom_a)
    compound.add(atom_b)
    compound.add_bond((atom_a, atom_b), bond_order=1.0)
    compound.bond_graph[atom_a][atom_b]["stereo"] = "E"
    validate_graph_equivalence(atom_graph, compound)

    atom_graph.nodes[0]["chirality"] = "S"
    with pytest.raises(ValueError, match="Conflicting stereo/chirality"):
        validate_graph_equivalence(atom_graph, compound)

    atom_graph.nodes[0]["chirality"] = "R"
    atom_graph.edges[0, 1]["bond_stereo"] = "Z"
    with pytest.raises(ValueError, match="Conflicting stereo/bond_stereo"):
        validate_graph_equivalence(atom_graph, compound)


def test_clash_threshold_is_configurable() -> None:
    """Authoring should reject nonbonded clashes at the requested threshold."""
    graph = _star_graph()
    inputs = _authoring_inputs()
    _, atom_graph = resolve_cgsmiles_graph(inputs["fragments"], graph)
    with pytest.raises(ValueError, match="severe nonbonded clash"):
        place_templates_on_cg_beads(
            graph,
            inputs["aa_templates"],
            atom_graph,
            atom_maps=inputs["atom_maps"],
            port_maps=inputs["port_maps"],
            leaving_maps=inputs["leaving_maps"],
            clash_threshold_nm=2.0,
        )


@pytest.mark.conjugation_stack
def test_tiny_three_arm_star_exports_valid_sdf_and_yaml(tmp_path: Path) -> None:
    """A degree-three star should resolve, stitch ports, charge, and validate."""
    graph = _star_graph()
    inputs = _authoring_inputs()
    result = build_and_export_polymer(
        cg_graph=graph,
        fragments=inputs["fragments"],
        aa_templates=inputs["aa_templates"],
        atom_maps=inputs["atom_maps"],
        port_maps=inputs["port_maps"],
        leaving_maps=inputs["leaving_maps"],
        output_sdf=tmp_path / "star.sdf",
        charge_method="gasteiger",
        name="star",
        yaml_base_path=tmp_path,
    )

    assert result["sdf_path"].exists()
    assert result["molecule"].n_atoms == 14
    snippet = yaml.safe_load(result["yaml_snippet"])
    assert snippet["polymers"]["generation_mode"] == "provided"
    assert snippet["polymers"]["provided_molecules"][0]["entries"][0]["sdf_path"] == "star.sdf"


def test_atomic_sdf_export_cleans_up_failed_temp(monkeypatch, tmp_path: Path) -> None:
    """Failed SDF export should remove same-directory temporary files."""
    graph = _star_graph()
    inputs = _authoring_inputs()
    _, atom_graph = resolve_cgsmiles_graph(inputs["fragments"], graph)
    compound = place_templates_on_cg_beads(
        graph,
        inputs["aa_templates"],
        atom_graph,
        atom_maps=inputs["atom_maps"],
        port_maps=inputs["port_maps"],
        leaving_maps=inputs["leaving_maps"],
    )

    def fail_to_file(self, path, file_format):
        Path(path).write_text("partial", encoding="utf-8")
        raise OSError("mock writer failure")

    from openff.toolkit import Molecule

    monkeypatch.setattr(Molecule, "to_file", fail_to_file)
    with pytest.raises(OSError, match="mock writer failure"):
        convert_charge_write_sdf(compound, tmp_path / "failed.sdf", atom_graph=atom_graph)
    assert not (tmp_path / "failed.sdf").exists()
    assert not list(tmp_path.glob(".*.tmp.sdf"))


def test_provided_molecules_yaml_schema_shape() -> None:
    """The emitted YAML should match the provided_molecules shape."""
    snippet = yaml.safe_load(emit_provided_molecules_yaml_snippet("generated/foo.sdf", name="foo"))
    assert snippet == {
        "polymers": {
            "enabled": True,
            "generation_mode": "provided",
            "provided_molecules": [
                {"name": "foo", "entries": [{"sdf_path": "generated/foo.sdf", "count": 1}]}
            ],
        }
    }


def test_provided_molecules_yaml_unwrapped_option() -> None:
    """The explicit unwrapped option should support merging into generated modes."""
    snippet = yaml.safe_load(
        emit_provided_molecules_yaml_snippet("generated/foo.sdf", name="foo", wrap_polymers=False)
    )
    assert snippet == {
        "provided_molecules": [
            {"name": "foo", "entries": [{"sdf_path": "generated/foo.sdf", "count": 1}]}
        ]
    }


def test_polymer_builder_provided_mode_skips_sequence_generation(monkeypatch) -> None:
    """Provided mode should return only provided molecule pools without dummy fields."""

    class FakeMolecule:
        n_atoms = 3

    sentinel = FakeMolecule()
    builder = PolymerBuilder(generation_mode="provided", provided_molecules=[object()])

    monkeypatch.setattr(builder, "_build_provided_molecules", lambda seed: ([sentinel], [2]))
    molecules, counts = builder.build(count=0, seed=10)

    assert molecules == [sentinel]
    assert counts == [2]
    assert builder.sequence_counts == Counter()
    assert builder.get_packing_info() == ([sentinel], [2])


def test_polymer_builder_provided_mode_validate_after_build(monkeypatch) -> None:
    """Provided mode validation should use packing molecules after a successful build."""

    class FakeMolecule:
        n_atoms = 3

    builder = PolymerBuilder(generation_mode="provided", provided_molecules=[object()])
    monkeypatch.setattr(builder, "_build_provided_molecules", lambda seed: ([FakeMolecule()], [1]))

    builder.build(count=0, seed=10)

    assert builder.validate() is True


def test_polymer_builder_provided_mode_validate_rejects_empty_state() -> None:
    """Provided mode validation should reject a builder that has not built pools."""
    builder = PolymerBuilder(generation_mode="provided", provided_molecules=[object()])

    with pytest.raises(RuntimeError, match="No provided molecules loaded"):
        builder.validate()


def test_polymer_builder_provided_mode_validate_rejects_empty_molecule(monkeypatch) -> None:
    """Provided mode validation should reject malformed provided molecule entries."""

    class EmptyMolecule:
        n_atoms = 0

    builder = PolymerBuilder(generation_mode="provided", provided_molecules=[object()])
    monkeypatch.setattr(builder, "_build_provided_molecules", lambda seed: ([EmptyMolecule()], [1]))
    builder.build(count=0, seed=10)

    with pytest.raises(ValueError, match="has no atoms"):
        builder.validate()


def test_notebook_resource_is_clean_and_compact() -> None:
    """The canonical notebook should be valid, clean, and user-edit minimal."""
    notebook_path = (
        Path(__file__).parents[2]
        / "src"
        / "polyzymd"
        / "templates"
        / "notebooks"
        / "cgsmiles_polymer_scaffold.ipynb"
    )
    notebook = nbformat.read(notebook_path, as_version=4)
    user_edit_cells = [
        cell
        for cell in notebook.cells
        if "polyzymd-user-edit" in cell.get("metadata", {}).get("tags", [])
    ]
    assert len(notebook.cells) == 5
    assert len(user_edit_cells) == 3
    assert all(
        cell.get("execution_count") is None for cell in notebook.cells if cell.cell_type == "code"
    )
    assert all(not cell.get("outputs") for cell in notebook.cells if cell.cell_type == "code")
