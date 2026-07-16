"""Offline mBuild and CGSmiles polymer authoring helpers."""

from __future__ import annotations

import os
import tempfile
from collections import Counter, deque
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import yaml

SUPPORTED_BOND_ORDERS = (1.0, 1.5, 2.0, 3.0)
BOND_ORDER_TOLERANCE = 1.0e-3


def path_to_cgsmiles_graph(coarse: Any) -> nx.Graph:
    """Convert an mBuild Path-like coarse object to a CGSmiles graph.

    Parameters
    ----------
    coarse : object
        mBuild coarse object exposing ``bond_graph`` plus either ``coordinates``
        or per-node positions. Nodes must resolve to a stable ``fragname``.

    Returns
    -------
    networkx.Graph
        Tree graph whose nodes contain ``fragname``, ``cg_name``, and
        ``position`` fields and whose edges contain integer ``order`` fields.
    """
    bond_graph = getattr(coarse, "bond_graph", None)
    if bond_graph is None:
        raise TypeError("Coarse object must expose a bond_graph")

    graph = nx.Graph()
    coordinates = getattr(coarse, "coordinates", None)
    beads = getattr(coarse, "beads", None)
    for node, data in bond_graph.nodes(data=True):
        node_id = int(node)
        cg_name = _coarse_node_name(node_id, data, beads)
        fragname = str(
            data.get("fragname") or (cg_name[1:] if cg_name.startswith("_") else cg_name)
        )
        position = _coarse_node_position(node_id, data, coordinates)
        graph.add_node(node_id, fragname=fragname, cg_name=cg_name, position=position)

    for u_node, v_node, data in bond_graph.edges(data=True):
        order = _normalize_bond_order(data, source=f"CG edge {(int(u_node), int(v_node))}")
        graph.add_edge(int(u_node), int(v_node), **{**dict(data), "order": order})

    validate_cgsmiles_graph(graph)
    return graph


def validate_cgsmiles_graph(cg_graph: nx.Graph) -> None:
    """Validate the coarse graph contract used by Phase15 authoring.

    Parameters
    ----------
    cg_graph : networkx.Graph
        CGSmiles coarse graph with integer node IDs and fragment names.
    """
    if cg_graph.number_of_nodes() == 0:
        raise ValueError("CGSmiles coarse graph must contain at least one node")
    duplicates = [item for item, count in Counter(list(cg_graph.nodes)).items() if count > 1]
    if duplicates:
        raise ValueError(f"Duplicate CG node identifiers are not allowed: {duplicates}")
    if not nx.is_tree(cg_graph):
        raise ValueError(
            "Phase15 CGSmiles authoring supports tree graphs only; cycles are rejected"
        )
    for node, data in cg_graph.nodes(data=True):
        if "fragname" not in data or not str(data["fragname"]).strip():
            raise ValueError(f"CG node {node!r} is missing a fragname")
        position = np.asarray(data.get("position", np.zeros(3)), dtype=float)
        if position.shape != (3,) or not np.isfinite(position).all():
            raise ValueError(f"CG node {node!r} has invalid coordinates")
    for u_node, v_node, data in cg_graph.edges(data=True):
        if u_node not in cg_graph or v_node not in cg_graph:
            raise ValueError(f"CG edge references a missing node: {(u_node, v_node)}")
        data["order"] = _normalize_bond_order(data, source=f"CG edge {(u_node, v_node)}")


def resolve_cgsmiles_graph(fragments: str, cg_graph: nx.Graph) -> tuple[nx.Graph, nx.Graph]:
    """Resolve a CGSmiles fragment string on a coarse graph.

    Parameters
    ----------
    fragments : str
        CGSmiles fragment definition string, for example ``"{#A=...}"``.
    cg_graph : networkx.Graph
        Validated coarse graph.

    Returns
    -------
    tuple[networkx.Graph, networkx.Graph]
        Resolved meta graph and atom graph from CGSmiles.
    """
    from cgsmiles import MoleculeResolver

    validate_cgsmiles_graph(cg_graph)
    cg_edge_orders = {
        tuple(sorted((int(u_node), int(v_node)))): _normalize_bond_order(
            data, source=f"CG edge {(u_node, v_node)}"
        )
        for u_node, v_node, data in cg_graph.edges(data=True)
    }
    resolver_graph = cg_graph.copy()
    for _, _, data in resolver_graph.edges(data=True):
        data["order"] = 1
    resolver = MoleculeResolver.from_graph(
        fragments, resolver_graph, last_all_atom=True, legacy=True
    )
    meta_graph, atom_graph = resolver.resolve_all()
    _apply_coarse_bond_orders_to_atom_graph(atom_graph, cg_edge_orders)
    validate_resolution(cg_graph, atom_graph)
    return meta_graph, atom_graph


def validate_resolution(cg_graph: nx.Graph, atom_graph: nx.Graph) -> set[tuple[int, int]]:
    """Validate that CGSmiles resolved every coarse edge exactly once.

    Parameters
    ----------
    cg_graph : networkx.Graph
        Coarse graph supplied to CGSmiles.
    atom_graph : networkx.Graph
        Resolved all-atom CGSmiles graph.

    Returns
    -------
    set[tuple[int, int]]
        Sorted coarse edge pairs represented by interfragment atomistic bonds.
    """
    cg_edge_orders = {
        tuple(sorted((int(u_node), int(v_node)))): _normalize_bond_order(
            data, source=f"CG edge {(u_node, v_node)}"
        )
        for u_node, v_node, data in cg_graph.edges(data=True)
    }
    records = interfragment_bond_records(atom_graph)
    _validate_one_interfragment_record_per_edge(cg_edge_orders, records)
    resolved_edges = {tuple(sorted((record["cg_u"], record["cg_v"]))) for record in records}
    missing = set(cg_edge_orders) - resolved_edges
    extra = resolved_edges - set(cg_edge_orders)
    if missing or extra:
        raise ValueError(f"Resolution mismatch. Missing={sorted(missing)}, extra={sorted(extra)}")
    return resolved_edges


def interfragment_bond_records(atom_graph: nx.Graph) -> list[dict[str, Any]]:
    """Return atomistic interfragment bond records from a CGSmiles atom graph.

    Parameters
    ----------
    atom_graph : networkx.Graph
        Resolved CGSmiles atom graph.

    Returns
    -------
    list[dict[str, object]]
        Records with coarse node IDs, fragment names, local atom indices,
        CGSmiles bonding labels, and bond orders.
    """
    records: list[dict[str, Any]] = []
    for atom_u, atom_v, data in atom_graph.edges(data=True):
        fragids_u = tuple(atom_graph.nodes[atom_u].get("fragid", ()))
        fragids_v = tuple(atom_graph.nodes[atom_v].get("fragid", ()))
        if fragids_u == fragids_v:
            continue
        if len(fragids_u) != 1 or len(fragids_v) != 1:
            raise ValueError("Phase15 requires one unsquashed CGSmiles fragment per CG bead")
        fragname_u, local_u = _fragment_local_atom(atom_graph, atom_u)
        fragname_v, local_v = _fragment_local_atom(atom_graph, atom_v)
        records.append(
            {
                "cg_u": int(fragids_u[0]),
                "cg_v": int(fragids_v[0]),
                "fragname_u": fragname_u,
                "fragname_v": fragname_v,
                "local_u": local_u,
                "local_v": local_v,
                "bonding": tuple(data.get("bonding", ())),
                "order": _normalize_bond_order(
                    data, source=f"resolved interfragment bond {(atom_u, atom_v)}"
                ),
            }
        )
    return records


def place_templates_on_cg_beads(
    cg_graph: nx.Graph,
    aa_templates: Mapping[str, Any],
    atom_graph: nx.Graph,
    *,
    atom_maps: Mapping[str, Mapping[int | str, str]],
    port_maps: Mapping[str, Mapping[str, str | Sequence[str]]],
    leaving_maps: Mapping[str, Mapping[str, str | Sequence[str]]],
    name: str = "backmapped_polymer",
    clash_threshold_nm: float = 0.10,
) -> Any:
    """Backmap and stitch atomistic mBuild templates on coarse CG beads.

    Parameters
    ----------
    cg_graph : networkx.Graph
        Tree coarse graph with node positions.
    aa_templates : mapping[str, object]
        mBuild atomistic templates keyed by CGSmiles fragment name.
    atom_graph : networkx.Graph
        Resolved CGSmiles atom graph.
    atom_maps : mapping[str, mapping[int | str, str]]
        Explicit maps from CGSmiles local atom indices to mBuild particle labels
        or names for each fragment.
    port_maps : mapping[str, mapping[str, str or sequence[str]]]
        Explicit maps from CGSmiles bonding labels to mBuild Port labels. Use
        a list of labels when one CGSmiles descriptor is intentionally reused.
    leaving_maps : mapping[str, mapping[str, str | sequence[str]]]
        Explicit named atoms to remove for each consumed bonding label.
    name : str, optional
        Output compound name, by default ``"backmapped_polymer"``.
    clash_threshold_nm : float, optional
        Minimum allowed nonbonded distance in nm, by default 0.10.

    Returns
    -------
    mbuild.Compound
        Atomistic mBuild compound assembled with mBuild ports and
        ``force_overlap(..., add_bond=True)``.
    """
    import mbuild as mb

    validate_cgsmiles_graph(cg_graph)
    _validate_authoring_maps(cg_graph, aa_templates, atom_graph, atom_maps, port_maps, leaving_maps)

    system = mb.Compound(name=name)
    instances = _clone_template_instances(cg_graph, aa_templates, atom_maps)
    for cg_node in sorted(instances):
        data = instances[cg_node]
        system.add(data["compound"], label=f"{data['fragname']}[{cg_node}]")

    used_ports: set[tuple[int, str]] = set()
    created_bonds: list[tuple[Any, Any, float]] = []
    for record in _ordered_tree_bond_records(cg_graph, interfragment_bond_records(atom_graph)):
        _consume_bond_record_with_ports(record, instances, port_maps, leaving_maps, used_ports)
        particle_u = instances[record["cg_u"]]["local_to_particle"][record["local_u"]]
        particle_v = instances[record["cg_v"]]["local_to_particle"][record["local_v"]]
        created_bonds.append((particle_u, particle_v, float(record["order"])))

    _normalize_created_bond_orders(system, created_bonds)
    _validate_mbuild_geometry(system, clash_threshold_nm=clash_threshold_nm)
    return system


def validate_graph_equivalence(
    atom_graph: nx.Graph, compound: Any, molecule: Any | None = None
) -> None:
    """Validate CGSmiles, mBuild, and optionally OpenFF graph equivalence.

    Parameters
    ----------
    atom_graph : networkx.Graph
        Resolved CGSmiles atom graph.
    compound : object
        Atomistic mBuild compound produced by this module.
    molecule : object or None, optional
        Optional OpenFF molecule converted from the mBuild compound.
    """
    mbuild_graph = _molecular_graph_from_mbuild(compound)
    cg_graph = _molecular_graph_from_cgsmiles(atom_graph)
    if not _graphs_match(cg_graph, mbuild_graph):
        raise ValueError("Resolved CGSmiles and mBuild atom graphs are not equivalent")
    if molecule is not None:
        openff_graph = _molecular_graph_from_openff(molecule)
        if not _graphs_match(cg_graph, openff_graph):
            raise ValueError("Resolved CGSmiles and OpenFF atom graphs are not equivalent")


def convert_charge_write_sdf(
    compound: Any,
    output_path: Path | str,
    *,
    charge_method: str = "gasteiger",
    atom_graph: nx.Graph | None = None,
) -> Any:
    """Convert an mBuild compound to OpenFF, charge, write SDF, and validate it.

    Parameters
    ----------
    compound : object
        Atomistic mBuild compound.
    output_path : pathlib.Path or str
        Destination SDF path.
    charge_method : str, optional
        OpenFF Toolkit partial charge method or ``"existing"``, by default
        ``"gasteiger"``.
    atom_graph : networkx.Graph or None, optional
        Optional CGSmiles atom graph for equivalence validation.

    Returns
    -------
    openff.toolkit.Molecule
        Charged OpenFF molecule loaded back through the provided molecule
        validator.
    """
    from polyzymd.builders.conjugation.polymer.mbuild import from_mbuild
    from polyzymd.builders.conjugation.polymer.provided_molecules import (
        load_validated_provided_molecule_sdf,
    )

    path = Path(output_path).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    if atom_graph is not None:
        validate_graph_equivalence(atom_graph, compound)
    molecule = from_mbuild(compound)
    if atom_graph is not None:
        validate_graph_equivalence(atom_graph, compound, molecule)
    _assign_partial_charges(molecule, charge_method)
    temp_path = _same_dir_temp_path(path, suffix=".sdf")
    try:
        molecule.to_file(temp_path, file_format="SDF")
        validated = load_validated_provided_molecule_sdf(temp_path)
        os.replace(temp_path, path)
    except Exception:
        temp_path.unlink(missing_ok=True)
        raise
    return validated


def build_and_export_polymer(
    *,
    cg_graph: nx.Graph,
    fragments: str,
    aa_templates: Mapping[str, Any],
    atom_maps: Mapping[str, Mapping[int | str, str]],
    port_maps: Mapping[str, Mapping[str, str | Sequence[str]]],
    leaving_maps: Mapping[str, Mapping[str, str | Sequence[str]]],
    output_sdf: Path | str,
    charge_method: str = "gasteiger",
    name: str = "authored_polymer",
    clash_threshold_nm: float = 0.10,
    yaml_base_path: Path | str | None = None,
) -> dict[str, Any]:
    """Run the compact Phase15 authoring workflow in one call.

    Parameters
    ----------
    cg_graph : networkx.Graph
        Coarse tree graph.
    fragments : str
        CGSmiles fragment definition string.
    aa_templates : mapping[str, object]
        Atomistic mBuild templates keyed by fragment name.
    atom_maps : mapping[str, mapping[int | str, str]]
        Explicit local atom maps.
    port_maps : mapping[str, mapping[str, str or sequence[str]]]
        Explicit bonding-label to mBuild Port maps.
    leaving_maps : mapping[str, mapping[str, str | sequence[str]]]
        Explicit leaving atom maps.
    output_sdf : pathlib.Path or str
        Destination charged SDF.
    charge_method : str, optional
        OpenFF partial charge method, by default ``"gasteiger"``.
    name : str, optional
        Molecule and compound name, by default ``"authored_polymer"``.
    clash_threshold_nm : float, optional
        Minimum allowed nonbonded distance in nm, by default 0.10.
    yaml_base_path : pathlib.Path or str or None, optional
        Base directory used to make the emitted SDF path relative in YAML.

    Returns
    -------
    dict[str, object]
        Result dictionary containing the coarse graph, resolved graphs,
        mBuild compound, OpenFF molecule, SDF path, and YAML snippet.
    """
    _, atom_graph = resolve_cgsmiles_graph(fragments, cg_graph)
    compound = place_templates_on_cg_beads(
        cg_graph,
        aa_templates,
        atom_graph,
        atom_maps=atom_maps,
        port_maps=port_maps,
        leaving_maps=leaving_maps,
        name=name,
        clash_threshold_nm=clash_threshold_nm,
    )
    molecule = convert_charge_write_sdf(
        compound,
        output_sdf,
        charge_method=charge_method,
        atom_graph=atom_graph,
    )
    snippet = emit_provided_molecules_yaml_snippet(output_sdf, name=name, base_path=yaml_base_path)
    snippet_path = Path(output_sdf).with_suffix(".provided_molecules.yaml")
    _atomic_write_text(snippet_path, snippet)
    return {
        "cg_graph": cg_graph,
        "atom_graph": atom_graph,
        "compound": compound,
        "molecule": molecule,
        "sdf_path": Path(output_sdf).resolve(),
        "yaml_snippet": snippet,
        "yaml_path": snippet_path.resolve(),
    }


def emit_provided_molecules_yaml_snippet(
    sdf_path: Path | str,
    *,
    name: str,
    count: int = 1,
    base_path: Path | str | None = None,
    wrap_polymers: bool = True,
) -> str:
    """Emit a configuration snippet for a generated provided molecule SDF.

    Parameters
    ----------
    sdf_path : pathlib.Path or str
        Path to the charged SDF file.
    name : str
        Provided molecule pool name.
    count : int, optional
        Number of copies to pack, by default 1.
    base_path : pathlib.Path or str or None, optional
        Base directory used to emit a relative SDF path, by default ``None``.
    wrap_polymers : bool, optional
        If ``True``, emit a full pasteable ``polymers`` block for provided-only
        mode. If ``False``, emit only the lower-level ``provided_molecules``
        list for merging into an existing dynamic or fragments block.

    Returns
    -------
    str
        YAML snippet suitable under ``polymers`` in ``config.yaml``.
    """
    if count < 1:
        raise ValueError("Provided molecule count must be positive")
    yaml_sdf_path = _display_sdf_path(Path(sdf_path), base_path)
    provided = [
        {
            "name": name,
            "entries": [{"sdf_path": yaml_sdf_path, "count": int(count)}],
        }
    ]
    payload = (
        {
            "polymers": {
                "enabled": True,
                "generation_mode": "provided",
                "provided_molecules": provided,
            }
        }
        if wrap_polymers
        else {"provided_molecules": provided}
    )
    return yaml.safe_dump(payload, sort_keys=False)


def _same_dir_temp_path(path: Path, *, suffix: str) -> Path:
    """Create a same-directory temporary path for atomic replacement."""
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temp_name = tempfile.mkstemp(
        prefix=f".{path.stem}.",
        suffix=f".tmp{suffix}",
        dir=path.parent,
    )
    os.close(handle)
    return Path(temp_name)


def _atomic_write_text(path: Path, text: str) -> None:
    """Write text atomically using a same-directory temporary file."""
    temp_path = _same_dir_temp_path(path, suffix=path.suffix)
    try:
        temp_path.write_text(text, encoding="utf-8")
        os.replace(temp_path, path)
    except Exception:
        temp_path.unlink(missing_ok=True)
        raise


def _display_sdf_path(sdf_path: Path, base_path: Path | str | None) -> str:
    """Return an SDF path for YAML, relative to a config base when requested."""
    if base_path is None:
        return str(sdf_path)
    try:
        return str(sdf_path.resolve().relative_to(Path(base_path).expanduser().resolve()))
    except ValueError:
        return str(sdf_path)


def _coarse_node_name(node_id: int, data: Mapping[str, Any], beads: Any) -> str:
    """Return a stable coarse bead name."""
    if data.get("name"):
        return str(data["name"])
    if beads is not None:
        return str(beads[node_id])
    raise ValueError(f"CG node {node_id!r} is missing a name or fragname")


def _phase15_single_bond_message(edge: tuple[Any, Any]) -> str:
    """Return the Phase15 unsupported-topology message."""
    return (
        "Phase15 supports exactly one atomistic interfragment bond per coarse edge. "
        f"CG edge {edge} is unsupported; PIM-1-like doubly linked or multi-point "
        "residues and cycles should be supplied as charged SDFs through "
        "provided_molecules."
    )


def _validate_one_interfragment_record_per_edge(
    cg_edge_orders: Mapping[tuple[int, int], float], records: Sequence[Mapping[str, Any]]
) -> None:
    """Reject duplicate, missing, extra, or order-mismatched interfragment records."""
    counts = Counter(
        tuple(sorted((int(record["cg_u"]), int(record["cg_v"])))) for record in records
    )
    missing = sorted(edge for edge in cg_edge_orders if counts[edge] == 0)
    duplicate = sorted(edge for edge, count in counts.items() if count > 1)
    extra = sorted(edge for edge in counts if edge not in cg_edge_orders)
    if missing or duplicate or extra:
        raise ValueError(
            _phase15_single_bond_message(("duplicate/missing", "interfragment"))
            + f" Missing={missing}, duplicate={duplicate}, extra={extra}."
        )
    mismatches = []
    for record in records:
        edge = tuple(sorted((int(record["cg_u"]), int(record["cg_v"]))))
        expected = cg_edge_orders[edge]
        actual = _normalize_bond_order(record, source=f"resolved interfragment record {edge}")
        if not _orders_match(expected, actual):
            mismatches.append((edge, expected, actual))
    if mismatches:
        raise ValueError(
            "Resolved interfragment bond order does not match CG edge order: " f"{mismatches}"
        )


def _apply_coarse_bond_orders_to_atom_graph(
    atom_graph: nx.Graph, cg_edge_orders: Mapping[tuple[int, int], float]
) -> None:
    """Copy normalized coarse bond orders onto resolved interfragment bonds."""
    for atom_u, atom_v, data in atom_graph.edges(data=True):
        fragids_u = tuple(atom_graph.nodes[atom_u].get("fragid", ()))
        fragids_v = tuple(atom_graph.nodes[atom_v].get("fragid", ()))
        if fragids_u == fragids_v:
            continue
        if len(fragids_u) != 1 or len(fragids_v) != 1:
            continue
        edge = tuple(sorted((int(fragids_u[0]), int(fragids_v[0]))))
        if edge in cg_edge_orders:
            data["order"] = cg_edge_orders[edge]


def _coarse_node_position(node_id: int, data: Mapping[str, Any], coordinates: Any) -> np.ndarray:
    """Return a finite coarse bead position."""
    if "position" in data:
        position = np.asarray(data["position"], dtype=float)
    elif "pos" in data:
        position = np.asarray(data["pos"], dtype=float)
    elif coordinates is not None:
        position = np.asarray(coordinates[node_id], dtype=float)
    else:
        position = np.zeros(3, dtype=float)
    if position.shape != (3,) or not np.isfinite(position).all():
        raise ValueError(f"CG node {node_id!r} has invalid coordinates")
    return position


def _normalize_bond_order(data: Mapping[str, Any], *, source: str) -> float:
    """Return a supported normalized bond order from edge metadata."""
    aromatic = _normalized_bool_alias(data, "aromatic", "is_aromatic", source=source)
    raw_order = _first_present(data, "order", "bond_order")
    try:
        order = float(1.0 if raw_order is None else raw_order)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Unsupported bond order {raw_order!r} from {source}") from exc
    if order <= 0.0:
        raise ValueError(f"Unsupported zero or negative bond order {order!r} from {source}")
    normalized = _nearest_supported_order(order, source=source)
    if aromatic is not None and bool(aromatic) != _orders_match(normalized, 1.5):
        raise ValueError(
            f"Conflicting aromatic metadata for {source}: order={normalized}, "
            f"aromatic={aromatic}"
        )
    return normalized


def _nearest_supported_order(order: float, *, source: str) -> float:
    """Return the supported bond order matching a numeric order."""
    for supported in SUPPORTED_BOND_ORDERS:
        if _orders_match(order, supported):
            return supported
    raise ValueError(
        f"Unsupported bond order {order!r} from {source}; supported orders are "
        "1.0, 1.5 aromatic, 2.0, and 3.0"
    )


def _orders_match(left: float, right: float) -> bool:
    """Return whether two numeric bond orders match within tolerance."""
    return abs(float(left) - float(right)) < BOND_ORDER_TOLERANCE


def _normalized_bool_alias(
    data: Mapping[str, Any], primary: str, alias: str, *, source: str
) -> bool | None:
    """Return a normalized bool alias pair or raise on conflicts."""
    primary_present = primary in data and data[primary] is not None
    alias_present = alias in data and data[alias] is not None
    if not primary_present and not alias_present:
        return None
    primary_value = bool(data[primary]) if primary_present else None
    alias_value = bool(data[alias]) if alias_present else None
    if primary_present and alias_present and primary_value != alias_value:
        raise ValueError(
            f"Conflicting {primary}/{alias} metadata for {source}: "
            f"{primary}={data[primary]!r}, {alias}={data[alias]!r}"
        )
    return primary_value if primary_present else alias_value


def _normalized_text_alias(
    data: Mapping[str, Any], primary: str, alias: str, *, source: str
) -> str | None:
    """Return a normalized text alias pair or raise on conflicts."""
    primary_value = data.get(primary)
    alias_value = data.get(alias)
    primary_text = str(primary_value) if primary_value is not None else None
    alias_text = str(alias_value) if alias_value is not None else None
    if primary_text is not None and alias_text is not None and primary_text != alias_text:
        raise ValueError(
            f"Conflicting {primary}/{alias} metadata for {source}: "
            f"{primary}={primary_value!r}, {alias}={alias_value!r}"
        )
    return primary_text if primary_text is not None else alias_text


def _resolved_cg_edge_pairs(atom_graph: nx.Graph) -> set[tuple[int, int]]:
    """Return coarse pairs represented by interfragment atomistic bonds."""
    pairs = set()
    for atom_u, atom_v in atom_graph.edges:
        fragids_u = tuple(atom_graph.nodes[atom_u].get("fragid", ()))
        fragids_v = tuple(atom_graph.nodes[atom_v].get("fragid", ()))
        if fragids_u == fragids_v:
            continue
        for fragid_u in fragids_u:
            for fragid_v in fragids_v:
                pairs.add(tuple(sorted((int(fragid_u), int(fragid_v)))))
    return pairs


def _fragment_local_atom(atom_graph: nx.Graph, atom_node: Any) -> tuple[str, int]:
    """Return fragment name and local atom index for a CGSmiles atom node."""
    mapping = atom_graph.nodes[atom_node].get("mapping")
    if not mapping:
        raise ValueError(f"Atom node {atom_node!r} lacks a CGSmiles mapping")
    fragname, local_index = mapping[0]
    return str(fragname), int(local_index)


def _validate_authoring_maps(
    cg_graph: nx.Graph,
    aa_templates: Mapping[str, Any],
    atom_graph: nx.Graph,
    atom_maps: Mapping[str, Mapping[int | str, str]],
    port_maps: Mapping[str, Mapping[str, str | Sequence[str]]],
    leaving_maps: Mapping[str, Mapping[str, str | Sequence[str]]],
) -> None:
    """Validate explicit atom, port, and leaving maps."""
    fragment_names = {str(data["fragname"]) for _, data in cg_graph.nodes(data=True)}
    missing_templates = fragment_names - set(aa_templates)
    if missing_templates:
        raise KeyError(f"Missing mBuild templates for fragments: {sorted(missing_templates)}")
    for fragname in sorted(fragment_names):
        if fragname not in atom_maps:
            raise ValueError(f"Missing explicit atom map for fragment {fragname!r}")
        if fragname not in port_maps:
            raise ValueError(f"Missing explicit port map for fragment {fragname!r}")
        if fragname not in leaving_maps:
            raise ValueError(f"Missing explicit leaving atom map for fragment {fragname!r}")

    for record in interfragment_bond_records(atom_graph):
        labels = tuple(record["bonding"])
        if len(labels) != 2:
            raise ValueError("Every interfragment bond must carry two CGSmiles bonding labels")
        for side, label in (("u", labels[0]), ("v", labels[1])):
            fragname = record[f"fragname_{side}"]
            local_index = record[f"local_{side}"]
            if (
                local_index not in atom_maps[fragname]
                and str(local_index) not in atom_maps[fragname]
            ):
                raise ValueError(
                    f"Missing atom map for fragment {fragname!r} local atom {local_index}"
                )
            _lookup_label_map(port_maps[fragname], label, "port", fragname)
            _lookup_label_map(leaving_maps[fragname], label, "leaving atom", fragname)


def _clone_template_instances(
    cg_graph: nx.Graph,
    aa_templates: Mapping[str, Any],
    atom_maps: Mapping[str, Mapping[int | str, str]],
) -> dict[int, dict[str, Any]]:
    """Clone and place one template instance on each coarse bead."""
    import mbuild as mb

    instances: dict[int, dict[str, Any]] = {}
    for cg_node, data in cg_graph.nodes(data=True):
        fragname = str(data["fragname"])
        instance = mb.clone(aa_templates[fragname])
        instance.name = f"{fragname}_{cg_node}"
        instance.translate_to(np.asarray(data["position"], dtype=float))
        local_to_particle = {
            int(local_index): _find_particle(instance, atom_name)
            for local_index, atom_name in _normalized_atom_map(atom_maps[fragname]).items()
        }
        for local_index, particle in local_to_particle.items():
            particle.cg_node = int(cg_node)
            particle.local_atom_index = int(local_index)
            particle.atom_identity = f"{cg_node}:{local_index}"
        instances[int(cg_node)] = {
            "compound": instance,
            "fragname": fragname,
            "local_to_particle": local_to_particle,
        }
    return instances


def _normalized_atom_map(atom_map: Mapping[int | str, str]) -> dict[int, str]:
    """Return an integer-keyed atom map and reject duplicates."""
    normalized = {int(key): value for key, value in atom_map.items()}
    if len(normalized) != len(atom_map):
        raise ValueError("Duplicate local atom map indices are not allowed")
    if len(set(normalized.values())) != len(normalized):
        raise ValueError("Duplicate mBuild atom names in an atom map are not allowed")
    return normalized


def _ordered_tree_bond_records(
    cg_graph: nx.Graph, records: Sequence[Mapping[str, Any]]
) -> list[Mapping[str, Any]]:
    """Return interfragment records ordered from a tree root outward."""
    cg_edge_orders = {
        tuple(sorted((int(u_node), int(v_node)))): _normalize_bond_order(
            data, source=f"CG edge {(u_node, v_node)}"
        )
        for u_node, v_node, data in cg_graph.edges(data=True)
    }
    _validate_one_interfragment_record_per_edge(cg_edge_orders, records)
    record_by_edge = {
        tuple(sorted((int(record["cg_u"]), int(record["cg_v"])))): record for record in records
    }
    root = min(int(node) for node in cg_graph.nodes)
    ordered = []
    seen = {root}
    queue: deque[int] = deque([root])
    while queue:
        parent = queue.popleft()
        for child in sorted(
            int(node) for node in cg_graph.neighbors(parent) if int(node) not in seen
        ):
            seen.add(child)
            queue.append(child)
            record = dict(record_by_edge[tuple(sorted((parent, child)))])
            if record["cg_u"] != parent:
                record = _swap_record_sides(record)
            ordered.append(record)
    return ordered


def _swap_record_sides(record: dict[str, Any]) -> dict[str, Any]:
    """Swap u/v fields so u is the already-placed parent."""
    swapped = dict(record)
    for suffix in ("cg", "fragname", "local"):
        swapped[f"{suffix}_u"], swapped[f"{suffix}_v"] = (
            record[f"{suffix}_v"],
            record[f"{suffix}_u"],
        )
    bonding = tuple(record["bonding"])
    swapped["bonding"] = (bonding[1], bonding[0])
    return swapped


def _consume_bond_record_with_ports(
    record: Mapping[str, Any],
    instances: Mapping[int, Mapping[str, Any]],
    port_maps: Mapping[str, Mapping[str, str | Sequence[str]]],
    leaving_maps: Mapping[str, Mapping[str, str | Sequence[str]]],
    used_ports: set[tuple[int, str]],
) -> None:
    """Remove explicit leaving atoms and stitch one tree edge with mBuild ports."""
    import mbuild as mb

    labels = tuple(record["bonding"])
    left = instances[int(record["cg_u"])]
    right = instances[int(record["cg_v"])]
    left_port_label = _select_port_label(
        port_maps[str(left["fragname"])],
        labels[0],
        int(record["cg_u"]),
        used_ports,
        left["fragname"],
    )
    right_port_label = _select_port_label(
        port_maps[str(right["fragname"])],
        labels[1],
        int(record["cg_v"]),
        used_ports,
        right["fragname"],
    )
    for cg_node, port_label in (
        (record["cg_u"], left_port_label),
        (record["cg_v"], right_port_label),
    ):
        key = (int(cg_node), str(port_label))
        if key in used_ports:
            raise ValueError(
                f"mBuild port {port_label!r} on CG node {cg_node} was used more than once"
            )
        used_ports.add(key)

    _remove_leaving_atoms(
        left["compound"], leaving_maps[str(left["fragname"])], labels[0], used_ports
    )
    _remove_leaving_atoms(
        right["compound"], leaving_maps[str(right["fragname"])], labels[1], used_ports
    )
    mb.force_overlap(
        move_this=right["compound"],
        from_positions=_find_port(right["compound"], right_port_label),
        to_positions=_find_port(left["compound"], left_port_label),
        add_bond=True,
    )


def _lookup_label_map(label_map: Mapping[str, Any], label: str, kind: str, fragname: str) -> Any:
    """Return a map value for an exact CGSmiles label."""
    label = str(label)
    if label in label_map:
        return label_map[label]
    raise ValueError(f"Missing {kind} map for CGSmiles label {label!r} on fragment {fragname!r}")


def _select_port_label(
    port_map: Mapping[str, str | Sequence[str]],
    label: str,
    cg_node: int,
    used_ports: set[tuple[int, str]],
    fragname: str,
) -> str:
    """Select one explicit unused mBuild Port label for a CGSmiles label."""
    value = _lookup_label_map(port_map, label, "port", str(fragname))
    candidates = (value,) if isinstance(value, str) else tuple(value)
    for candidate in candidates:
        port_label = str(candidate)
        if (cg_node, port_label) not in used_ports:
            return port_label
    raise ValueError(f"No unused mBuild Port remains for CGSmiles label {label!r}")


def _remove_leaving_atoms(
    compound: Any,
    leaving_map: Mapping[str, str | Sequence[str]],
    label: str,
    used_ports: set[tuple[int, str]],
) -> None:
    """Remove explicitly named leaving atoms for one consumed port."""
    value = _lookup_label_map(leaving_map, label, "leaving atom", getattr(compound, "name", ""))
    atom_names = (value,) if isinstance(value, str) else tuple(value)
    if not atom_names:
        raise ValueError("Leaving atom map entries must name at least one atom")
    for atom_name in atom_names:
        key = (id(compound), str(label), str(atom_name))
        if key in used_ports:
            continue
        try:
            particle = _find_particle(compound, atom_name)
        except ValueError:
            continue
        compound.remove(particle)
        used_ports.add(key)
        return
    raise ValueError(f"No explicit leaving atom remains for CGSmiles label {label!r}")


def _find_particle(compound: Any, atom_name: str) -> Any:
    """Find one mBuild particle by label or name."""
    particles = tuple(compound.particles())
    if atom_name in getattr(compound, "labels", {}):
        labelled = compound.labels[atom_name]
        if labelled in particles:
            return labelled
    matches = [particle for particle in particles if getattr(particle, "name", None) == atom_name]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one mBuild atom named or labelled {atom_name!r}")
    return matches[0]


def _find_port(compound: Any, port_label: str) -> Any:
    """Find one mBuild Port by label."""
    port = getattr(compound, "labels", {}).get(port_label)
    if port is None:
        raise ValueError(f"Missing mBuild Port labelled {port_label!r}")
    if not hasattr(port, "anchor"):
        raise TypeError(f"Label {port_label!r} does not reference an mBuild Port")
    return port


def _normalize_created_bond_orders(
    compound: Any, created_bonds: Sequence[tuple[Any, Any, float]]
) -> None:
    """Set bond orders only for bonds created by force_overlap."""
    bond_graph = getattr(compound, "bond_graph", None)
    if bond_graph is None:
        return
    for atom_u, atom_v, order in created_bonds:
        try:
            bond_graph[atom_u][atom_v]["bond_order"] = float(order)
            bond_graph[atom_u][atom_v]["order"] = float(order)
        except KeyError:
            compound.add_bond((atom_u, atom_v), bond_order=float(order))


def _validate_mbuild_geometry(compound: Any, *, clash_threshold_nm: float) -> None:
    """Validate finite coordinates and reject severe nonbonded clashes."""
    if clash_threshold_nm < 0.0:
        raise ValueError("Clash threshold must be non-negative")
    all_particles = [
        particle
        for particle in compound.particles()
        if not bool(getattr(particle, "port_particle", False))
    ]
    particles = [
        particle for particle in all_particles if _particle_symbol(particle).upper() != "H"
    ]
    coordinates = np.asarray([particle.pos for particle in all_particles], dtype=float)
    if coordinates.shape != (len(all_particles), 3) or not np.isfinite(coordinates).all():
        raise ValueError("Backmapped mBuild compound contains invalid coordinates")
    bonded = {frozenset(pair) for pair in compound.bonds()}
    for index, particle_i in enumerate(particles):
        for particle_j in particles[index + 1 :]:
            if frozenset((particle_i, particle_j)) in bonded:
                continue
            if _same_immediate_parent(particle_i, particle_j):
                continue
            if np.linalg.norm(particle_i.pos - particle_j.pos) < clash_threshold_nm:
                raise ValueError("Backmapped mBuild compound contains a severe nonbonded clash")


def _same_immediate_parent(particle_i: Any, particle_j: Any) -> bool:
    """Return whether two particles belong to the same template instance."""
    parent_i = getattr(particle_i, "parent", None)
    parent_j = getattr(particle_j, "parent", None)
    return parent_i is not None and parent_i is parent_j


def _assign_partial_charges(molecule: Any, charge_method: str) -> None:
    """Assign or validate partial charges on an OpenFF molecule."""
    import numpy as local_np

    if charge_method == "existing":
        if molecule.partial_charges is None:
            raise ValueError("charge_method='existing' requires complete input partial charges")
    else:
        molecule.assign_partial_charges(charge_method)
    charges = molecule.partial_charges
    if charges is None or len(charges) != molecule.n_atoms:
        raise ValueError("OpenFF molecule has incomplete partial charges after charge assignment")
    if not local_np.isfinite(charges.m_as("elementary_charge")).all():
        raise ValueError("OpenFF molecule has non-finite partial charges after charge assignment")


def _molecular_graph_from_cgsmiles(atom_graph: nx.Graph) -> nx.Graph:
    """Return a labelled molecular graph from resolved CGSmiles."""
    graph = nx.Graph()
    for node, data in atom_graph.nodes(data=True):
        graph.add_node(
            node,
            element=str(data.get("element", "")).capitalize(),
            **_node_chemistry(data, source=f"CGSmiles atom {node}"),
        )
    for atom_u, atom_v, data in atom_graph.edges(data=True):
        graph.add_edge(
            atom_u,
            atom_v,
            order=_normalize_bond_order(data, source=f"CGSmiles bond {(atom_u, atom_v)}"),
            **_edge_chemistry(data, source=f"CGSmiles bond {(atom_u, atom_v)}"),
        )
    return graph


def _molecular_graph_from_mbuild(compound: Any) -> nx.Graph:
    """Return a labelled molecular graph from an mBuild compound."""
    particles = [
        particle
        for particle in compound.particles()
        if not bool(getattr(particle, "port_particle", False))
    ]
    index_by_particle = {particle: index for index, particle in enumerate(particles)}
    graph = nx.Graph()
    for index, particle in enumerate(particles):
        graph.add_node(
            index,
            element=_particle_symbol(particle).capitalize(),
            **_particle_chemistry(particle),
        )
    for atom_u, atom_v in compound.bonds():
        if atom_u in index_by_particle and atom_v in index_by_particle:
            graph.add_edge(
                index_by_particle[atom_u],
                index_by_particle[atom_v],
                order=_mbuild_graph_bond_order(compound, atom_u, atom_v),
                **_edge_chemistry(
                    _mbuild_graph_bond_data(compound, atom_u, atom_v),
                    source=f"mBuild bond {(index_by_particle[atom_u], index_by_particle[atom_v])}",
                ),
            )
    return graph


def _molecular_graph_from_openff(molecule: Any) -> nx.Graph:
    """Return a labelled molecular graph from an OpenFF molecule."""
    graph = nx.Graph()
    for atom in molecule.atoms:
        graph.add_node(
            atom.molecule_atom_index,
            element=atom.symbol.capitalize(),
            formal_charge=int(atom.formal_charge.m_as("elementary_charge")),
            **_openff_atom_chemistry(atom),
        )
    for bond in molecule.bonds:
        graph.add_edge(
            bond.atom1_index,
            bond.atom2_index,
            order=float(bond.bond_order),
            **_openff_bond_chemistry(bond),
        )
    return graph


def _graphs_match(left: nx.Graph, right: nx.Graph) -> bool:
    """Return whether two labelled molecular graphs are isomorphic."""
    return nx.is_isomorphic(
        left,
        right,
        node_match=_node_attrs_match,
        edge_match=_edge_attrs_match,
    )


def _node_attrs_match(left: Mapping[str, Any], right: Mapping[str, Any]) -> bool:
    """Return whether node chemistry attributes match exactly where provided."""
    keys = set(left) | set(right)
    return all(left.get(key) == right.get(key) for key in keys)


def _edge_attrs_match(left: Mapping[str, Any], right: Mapping[str, Any]) -> bool:
    """Return whether edge chemistry attributes match exactly where provided."""
    if not _orders_match(float(left.get("order", 1)), float(right.get("order", 1))):
        return False
    keys = (set(left) | set(right)) - {"order"}
    return all(left.get(key) == right.get(key) for key in keys)


def _node_chemistry(data: Mapping[str, Any], *, source: str) -> dict[str, Any]:
    """Return node chemistry metadata that source data explicitly provides."""
    attrs: dict[str, Any] = {}
    charge = _first_present(data, "formal_charge", "formalcharge")
    attrs["formal_charge"] = int(charge or 0)
    aromatic = _normalized_bool_alias(data, "aromatic", "is_aromatic", source=source)
    if aromatic is not None and bool(aromatic):
        attrs["aromatic"] = bool(aromatic)
    stereo = _normalized_text_alias(data, "stereo", "chirality", source=source)
    if stereo is not None:
        attrs["stereo"] = stereo
    return attrs


def _particle_chemistry(particle: Any) -> dict[str, Any]:
    """Return mBuild particle chemistry metadata explicitly present on a particle."""
    attrs: dict[str, Any] = {}
    attrs["formal_charge"] = int(getattr(particle, "formal_charge", 0) or 0)
    if hasattr(particle, "aromatic"):
        attrs["aromatic"] = bool(particle.aromatic)
    stereo = _normalized_particle_text_alias(particle, "stereo", "chirality")
    if stereo is not None:
        attrs["stereo"] = stereo
    return attrs


def _openff_atom_chemistry(atom: Any) -> dict[str, Any]:
    """Return optional OpenFF atom chemistry metadata."""
    attrs: dict[str, Any] = {}
    if getattr(atom, "is_aromatic", None):
        attrs["aromatic"] = bool(atom.is_aromatic)
    if getattr(atom, "stereochemistry", None) is not None:
        attrs["stereo"] = str(atom.stereochemistry)
    return attrs


def _edge_chemistry(data: Mapping[str, Any], *, source: str) -> dict[str, Any]:
    """Return edge chemistry metadata that source data explicitly provides."""
    attrs: dict[str, Any] = {}
    stereo = _normalized_text_alias(data, "stereo", "bond_stereo", source=source)
    if stereo is not None:
        attrs["stereo"] = stereo
    return attrs


def _openff_bond_chemistry(bond: Any) -> dict[str, Any]:
    """Return optional OpenFF bond chemistry metadata."""
    attrs: dict[str, Any] = {}
    if getattr(bond, "stereochemistry", None) is not None:
        attrs["stereo"] = str(bond.stereochemistry)
    return attrs


def _first_present(data: Mapping[str, Any], *keys: str) -> Any | None:
    """Return the first present non-None mapping value."""
    for key in keys:
        if key in data and data[key] is not None:
            return data[key]
    return None


def _normalized_particle_text_alias(particle: Any, primary: str, alias: str) -> str | None:
    """Return a normalized particle text alias pair or raise on conflicts."""
    primary_present = hasattr(particle, primary)
    alias_present = hasattr(particle, alias)
    primary_value = getattr(particle, primary) if primary_present else None
    alias_value = getattr(particle, alias) if alias_present else None
    primary_text = str(primary_value) if primary_value is not None else None
    alias_text = str(alias_value) if alias_value is not None else None
    if primary_text is not None and alias_text is not None and primary_text != alias_text:
        raise ValueError(
            f"Conflicting {primary}/{alias} metadata for particle "
            f"{getattr(particle, 'name', '<unnamed>')}: {primary}={primary_value!r}, "
            f"{alias}={alias_value!r}"
        )
    return primary_text if primary_text is not None else alias_text


def _particle_symbol(particle: Any) -> str:
    """Return an element symbol from an mBuild particle."""
    element = getattr(particle, "element", None)
    if isinstance(element, str) and element.strip():
        return element.strip()
    symbol = getattr(element, "symbol", None)
    if symbol:
        return str(symbol)
    return str(getattr(particle, "name", ""))[:1]


def _mbuild_graph_bond_order(compound: Any, atom_u: Any, atom_v: Any) -> float:
    """Return stored mBuild graph bond order."""
    data = _mbuild_graph_bond_data(compound, atom_u, atom_v)
    return _normalize_bond_order(
        data, source=f"mBuild bond {(_particle_symbol(atom_u), _particle_symbol(atom_v))}"
    )


def _mbuild_graph_bond_data(compound: Any, atom_u: Any, atom_v: Any) -> Mapping[str, Any]:
    """Return stored mBuild graph bond metadata."""
    try:
        return compound.bond_graph[atom_u][atom_v]
    except (AttributeError, KeyError, TypeError):
        return {}
