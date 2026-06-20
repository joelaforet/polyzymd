"""Generic local product-state NAGL patch charges for conjugation junctions."""

from __future__ import annotations

import math
import os
from collections import deque
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation.pablo.charge_records import AtomPartialChargeRecord

DEFAULT_PATCH_NAGL_MODEL = "openff-gnn-am1bcc-0.1.0-rc.3.pt"
PATCH_NAGL_MODEL_ENV = "POLYZYMD_CONJUGATE_PATCH_NAGL_MODEL"
PATCH_SOURCE = "production:product-state-local-nagl-charge-bridge"


class LocalChargePatchError(ValueError):
    """Raised when local product-state patch metadata is insufficient."""


class ConjugationChargeConfig(BaseModel):
    """Configuration for local product-state patch charging."""

    patch_radius_bonds: int = Field(default=1, ge=1)
    nagl_model: str | None = None


@dataclass(frozen=True)
class _PatchAtom:
    """Atom metadata used to build and map a local charge patch."""

    key: tuple[str, int | None, str]
    element: str
    atom_name: str
    residue_name: str
    residue_number: int | None
    chain_id: str
    insertion_code: str
    source_atom: Any
    is_real: bool = True


def build_local_product_charge_patch_records(
    spec: Any,
    *,
    product_atoms: Sequence[Any] = (),
    config: ConjugationChargeConfig | None = None,
) -> tuple[tuple[AtomPartialChargeRecord, ...], str]:
    """Build atom charge records from a generic product-state local patch.

    Parameters
    ----------
    spec : Any
        Attachment build spec containing a resolved plan and generated fragment.
    product_atoms : sequence of Any, optional
        Product PDB atoms used to recover source-residue atoms when available,
        by default ``()``.
    config : ConjugationChargeConfig or None, optional
        Local patch charge settings, by default ``None``.

    Returns
    -------
    tuple of tuple and str
        Atom-level charge records for real product atoms and the NAGL model name.
    """
    charge_config = _resolve_charge_config(spec, config)
    graph = _build_product_graph(spec, product_atoms=tuple(product_atoms))
    selected = _select_patch_keys(graph, radius=charge_config.patch_radius_bonds)
    molecule, map_numbers = _build_off_patch_molecule(graph, selected)
    model_name = charge_config.nagl_model or os.environ.get(
        PATCH_NAGL_MODEL_ENV, DEFAULT_PATCH_NAGL_MODEL
    )
    charged = _charge_with_nagl(molecule, model_name=model_name)
    charges = _partial_charges(charged)
    records = []
    for key, map_number in map_numbers.items():
        atom = graph.atoms[key]
        if not atom.is_real:
            continue
        records.append(
            AtomPartialChargeRecord(
                chain_id=atom.chain_id,
                residue_name=atom.residue_name,
                residue_number=atom.residue_number,
                insertion_code=atom.insertion_code,
                atom_name=atom.atom_name,
                charge_e=charges[map_number],
                source=f"{PATCH_SOURCE}:{model_name}:radius-{charge_config.patch_radius_bonds}",
                source_role="local_nagl_patch",
            )
        )
    if not records:
        raise LocalChargePatchError("Local product-state patch did not map any real atoms")
    return tuple(records), model_name


@dataclass
class _PatchGraph:
    """Product-state graph around one covalent conjugation junction."""

    atoms: dict[tuple[str, int | None, str], _PatchAtom]
    bonds: dict[tuple[str, int | None, str], dict[tuple[str, int | None, str], int]]
    roots: tuple[tuple[str, int | None, str], tuple[str, int | None, str]]


def _resolve_charge_config(
    spec: Any, config: ConjugationChargeConfig | None
) -> ConjugationChargeConfig:
    """Resolve local patch configuration from explicit config or spec metadata."""
    if config is not None:
        return config
    raw = getattr(spec, "charge_config", None) or getattr(
        getattr(spec, "attachment_config", None), "charge", None
    )
    if raw is None:
        return ConjugationChargeConfig()
    if isinstance(raw, ConjugationChargeConfig):
        return raw
    patch_radius = getattr(raw, "patch_radius_bonds", None)
    nagl_model = getattr(raw, "nagl_model", None) or getattr(raw, "model_name", None)
    updates = {}
    if patch_radius is not None:
        updates["patch_radius_bonds"] = int(patch_radius)
    if nagl_model is not None:
        updates["nagl_model"] = str(nagl_model)
    return ConjugationChargeConfig(**updates)


def _build_product_graph(spec: Any, *, product_atoms: tuple[Any, ...]) -> _PatchGraph:
    """Build a product-state graph from resolved plan and fragment metadata."""
    plan = getattr(spec, "resolved_plan", None)
    fragment = getattr(spec, "generated_fragment", None)
    if plan is None:
        raise LocalChargePatchError("Local charge patch requires a resolved attachment plan")
    if fragment is None:
        raise LocalChargePatchError("Local charge patch requires generated fragment metadata")
    protein_link = getattr(plan, "protein_link_atom", None)
    modifier_link = getattr(plan, "modifier_link_atom", None)
    if protein_link is None or modifier_link is None:
        raise LocalChargePatchError("Local charge patch requires both resolved link atoms")

    atoms: dict[tuple[str, int | None, str], _PatchAtom] = {}
    bonds: dict[tuple[str, int | None, str], dict[tuple[str, int | None, str], int]] = {}
    roots = (_atom_key(protein_link, role="protein"), _atom_key(modifier_link, role="modifier"))
    leaving = _leaving_keys(plan, fragment)

    for atom in product_atoms:
        if _same_product_residue(
            atom, protein_link, getattr(plan, "protein_product_residue_name", None)
        ):
            key = _atom_key(atom, role="protein")
            if key not in leaving:
                atoms[key] = _patch_atom(atom, key=key, role="protein")
    if roots[0] not in atoms:
        atoms[roots[0]] = _patch_atom(
            protein_link,
            key=roots[0],
            role="protein",
            residue_name=getattr(plan, "protein_product_residue_name", None),
        )

    for atom in tuple(getattr(fragment, "atoms", ()) or ()):
        key = _atom_key(atom, role="modifier")
        if key not in leaving:
            atoms[key] = _patch_atom(
                atom,
                key=key,
                role="modifier",
                residue_name=_mapped_residue_name(spec, atom),
                residue_number=_mapped_residue_number(spec, atom),
                chain_id=_mapped_chain_id(spec, atom),
            )
    if roots[1] not in atoms:
        atoms[roots[1]] = _patch_atom(modifier_link, key=roots[1], role="modifier")

    for left, right, order in _fragment_bonds(fragment):
        left_key = _atom_key(left, role="modifier")
        right_key = _atom_key(right, role="modifier")
        if left_key in atoms and right_key in atoms:
            _add_bond(bonds, left_key, right_key, int(round(float(order))))
    _add_bond(
        bonds,
        roots[0],
        roots[1],
        int(
            round(
                float(getattr(getattr(plan, "pablo_crosslink_requirement", None), "bond_order", 1))
            )
        ),
    )
    _validate_roots_connected(atoms, bonds, roots)
    return _PatchGraph(atoms=atoms, bonds=bonds, roots=roots)


def _select_patch_keys(graph: _PatchGraph, *, radius: int) -> set[tuple[str, int | None, str]]:
    """Select a graph neighborhood around both link atoms."""
    selected: set[tuple[str, int | None, str]] = set()
    queue = deque((root, 0) for root in graph.roots)
    while queue:
        key, distance = queue.popleft()
        if key in selected:
            continue
        selected.add(key)
        if distance >= radius:
            continue
        for neighbor in graph.bonds.get(key, {}):
            queue.append((neighbor, distance + 1))
    return selected


def _build_off_patch_molecule(
    graph: _PatchGraph, selected: set[tuple[str, int | None, str]]
) -> tuple[Any, dict[tuple[str, int | None, str], int]]:
    """Build an OpenFF patch molecule with explicit generic caps."""
    from openff.toolkit import Molecule
    from rdkit import Chem

    rwmol = Chem.RWMol()
    rd_indices: dict[tuple[str, int | None, str], int] = {}
    map_numbers: dict[tuple[str, int | None, str], int] = {}
    for map_number, key in enumerate(sorted(selected), start=1):
        atom = graph.atoms[key]
        rd_atom = Chem.Atom(atom.element or "C")
        rd_atom.SetAtomMapNum(map_number)
        rd_indices[key] = rwmol.AddAtom(rd_atom)
        map_numbers[key] = map_number - 1
    for left in sorted(selected):
        for right, order in graph.bonds.get(left, {}).items():
            if right not in selected or rd_indices[left] >= rd_indices[right]:
                continue
            rwmol.AddBond(rd_indices[left], rd_indices[right], _rdkit_bond_type(order))
    _add_boundary_hydrogen_caps(rwmol, graph, selected, rd_indices)
    mol = rwmol.GetMol()
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol, catchErrors=True)
    mol = Chem.AddHs(mol, addCoords=True)
    offmol = Molecule.from_rdkit(mol, allow_undefined_stereo=True, hydrogens_are_explicit=True)
    for rd_atom, off_atom in zip(mol.GetAtoms(), offmol.atoms, strict=True):
        map_number = rd_atom.GetAtomMapNum()
        if map_number:
            off_atom.metadata["atom_map"] = map_number
    return offmol, map_numbers


def _add_boundary_hydrogen_caps(
    rwmol: Any, graph: _PatchGraph, selected: set[Any], rd_indices: dict[Any, int]
) -> None:
    """Cap selected atoms that have omitted neighbors with hydrogens."""
    from rdkit import Chem

    for key in sorted(selected):
        omitted = [neighbor for neighbor in graph.bonds.get(key, {}) if neighbor not in selected]
        for _neighbor in omitted:
            cap_index = rwmol.AddAtom(Chem.Atom("H"))
            rwmol.AddBond(rd_indices[key], cap_index, Chem.BondType.SINGLE)


def _charge_with_nagl(molecule: Any, *, model_name: str) -> Any:
    """Assign NAGL charges to a local patch molecule."""
    from polyzymd.utils.charging import NAGLCharger

    charger = NAGLCharger(model_name=model_name)
    charged = charger.charge_molecule(molecule)
    if getattr(charged, "partial_charges", None) is None:
        raise LocalChargePatchError(f"NAGL model {model_name} did not return patch charges")
    return charged


def _partial_charges(molecule: Any) -> dict[int, float]:
    """Return charged patch values by zero-based real-atom map index."""
    quantity = getattr(molecule, "partial_charges", None)
    if quantity is None:
        raise LocalChargePatchError("Charged patch molecule is missing partial charges")
    if hasattr(quantity, "m_as"):
        values = tuple(float(value) for value in quantity.m_as("elementary_charge"))
    else:
        values = tuple(float(value) for value in quantity)
    charges: dict[int, float] = {}
    for atom, charge in zip(tuple(molecule.atoms), values, strict=True):
        map_number = (getattr(atom, "metadata", None) or {}).get("atom_map")
        if not map_number:
            continue
        if not math.isfinite(charge):
            raise LocalChargePatchError("NAGL returned a non-finite local patch charge")
        charges[int(map_number) - 1] = charge
    return charges


def _fragment_bonds(fragment: Any) -> tuple[tuple[Any, Any, float], ...]:
    """Return generated-fragment bonds with bond-order metadata."""
    atoms = tuple(getattr(fragment, "atoms", ()) or ())
    lookup = _fragment_atom_lookup(atoms)
    bonds: list[tuple[Any, Any, float]] = []
    seen: set[frozenset[Any]] = set()
    for left, right, order in tuple(getattr(fragment, "bond_orders", ()) or ()):
        left_atom = _resolve_fragment_atom(left, lookup)
        right_atom = _resolve_fragment_atom(right, lookup)
        if left_atom is None or right_atom is None:
            continue
        seen.add(frozenset((id(left_atom), id(right_atom))))
        bonds.append((left_atom, right_atom, float(order)))
    for left, right in tuple(getattr(fragment, "bonds", ()) or ()):
        left_atom = _resolve_fragment_atom(left, lookup)
        right_atom = _resolve_fragment_atom(right, lookup)
        if left_atom is None or right_atom is None:
            continue
        key = frozenset((id(left_atom), id(right_atom)))
        if key not in seen:
            bonds.append((left_atom, right_atom, 1.0))
    return tuple(bonds)


def _fragment_atom_lookup(atoms: tuple[Any, ...]) -> dict[Any, Any]:
    """Build lookup keys for generated-fragment atoms."""
    lookup: dict[Any, Any] = {}
    for atom in atoms:
        for value in (
            getattr(atom, "atom_index", None),
            getattr(atom, "serial", None),
            getattr(atom, "atom_name", None),
            getattr(atom, "name", None),
        ):
            if value not in (None, ""):
                lookup[value] = atom
                lookup[str(value)] = atom
    return lookup


def _resolve_fragment_atom(value: Any, lookup: Mapping[Any, Any]) -> Any | None:
    """Resolve a bond endpoint to a generated-fragment atom."""
    return lookup.get(value) or lookup.get(str(value))


def _leaving_keys(plan: Any, fragment: Any) -> set[tuple[str, int | None, str]]:
    """Return leaving atom keys from resolved plan and fragment metadata."""
    keys = {_atom_key(atom, role="protein") for atom in getattr(plan, "protein_leaving_atoms", ())}
    keys.update(
        _atom_key(atom, role="modifier") for atom in getattr(plan, "modifier_leaving_atoms", ())
    )
    leaving_names = {
        str(name).strip().upper() for name in getattr(fragment, "leaving_atom_names", ())
    }
    leaving_indices = {int(index) for index in getattr(fragment, "leaving_atom_indices", ())}
    leaving_serials = {int(serial) for serial in getattr(fragment, "leaving_atom_serials", ())}
    for atom in tuple(getattr(fragment, "atoms", ()) or ()):
        if (
            str(getattr(atom, "atom_name", "")).strip().upper() in leaving_names
            or getattr(atom, "atom_index", None) in leaving_indices
            or getattr(atom, "serial", None) in leaving_serials
        ):
            keys.add(_atom_key(atom, role="modifier"))
    return keys


def _same_product_residue(atom: Any, link_atom: Any, product_residue_name: str | None) -> bool:
    """Return whether an atom belongs to the product-side protein residue."""
    return (
        str(getattr(atom, "chain_id", "") or "").strip()
        == str(getattr(link_atom, "chain_id", "") or "").strip()
        and getattr(atom, "residue_number", None) == getattr(link_atom, "residue_number", None)
        and str(getattr(atom, "insertion_code", "") or "").strip()
        == str(getattr(link_atom, "insertion_code", "") or "").strip()
        and (
            product_residue_name is None
            or str(getattr(atom, "residue_name", "") or "").strip().upper()
            == str(product_residue_name).strip().upper()
        )
    )


def _patch_atom(
    atom: Any,
    *,
    key: tuple[str, int | None, str],
    role: str,
    residue_name: str | None = None,
    residue_number: int | None = None,
    chain_id: str | None = None,
) -> _PatchAtom:
    """Convert PDB-like metadata to a patch atom."""
    return _PatchAtom(
        key=key,
        element=_atom_element(atom),
        atom_name=str(getattr(atom, "atom_name", "") or getattr(atom, "name", "")).strip(),
        residue_name=(residue_name or str(getattr(atom, "residue_name", "") or "")).strip().upper(),
        residue_number=(
            residue_number if residue_number is not None else getattr(atom, "residue_number", None)
        ),
        chain_id=(chain_id if chain_id is not None else getattr(atom, "chain_id", "") or ""),
        insertion_code=str(getattr(atom, "insertion_code", "") or "").strip(),
        source_atom=atom,
        is_real=role in {"protein", "modifier"},
    )


def _atom_key(atom: Any, *, role: str) -> tuple[str, int | None, str]:
    """Return a stable atom key within one patch graph."""
    serial = getattr(atom, "serial", None)
    index = getattr(atom, "atom_index", None)
    name = str(getattr(atom, "atom_name", "") or getattr(atom, "name", "")).strip().upper()
    residue = getattr(atom, "residue_number", None)
    return (
        role,
        (
            int(index if index is not None else serial)
            if index is not None or serial is not None
            else residue
        ),
        name,
    )


def _atom_element(atom: Any) -> str:
    """Return an element symbol from PDB-like atom metadata."""
    element = str(getattr(atom, "element", "") or "").strip()
    if element:
        return element.capitalize()
    name = str(getattr(atom, "atom_name", "") or getattr(atom, "name", "")).strip()
    letters = "".join(char for char in name if char.isalpha())
    if not letters:
        raise LocalChargePatchError(f"Could not infer element for patch atom {name!r}")
    return letters[0].upper()


def _mapped_residue_name(spec: Any, atom: Any) -> str:
    """Return product residue name for a generated-fragment atom."""
    mapping = _mapping_for_atom(spec, atom)
    return str(mapping.get("target_residue_name", getattr(atom, "residue_name", "")))


def _mapped_residue_number(spec: Any, atom: Any) -> int | None:
    """Return product residue number for a generated-fragment atom."""
    mapping = _mapping_for_atom(spec, atom)
    value = mapping.get("target_residue_number", getattr(atom, "residue_number", None))
    return int(value) if value not in (None, "") else None


def _mapped_chain_id(spec: Any, atom: Any) -> str:
    """Return product chain ID for a generated-fragment atom."""
    mapping = _mapping_for_atom(spec, atom)
    return str(mapping.get("target_chain", getattr(atom, "chain_id", "C") or "C"))


def _mapping_for_atom(spec: Any, atom: Any) -> Mapping[str, Any]:
    """Return product residue mapping for a generated-fragment atom."""
    mappings = getattr(spec, "product_residue_mappings", {}) or {}
    source_number = str(getattr(atom, "residue_number", "") or "")
    return (
        mappings.get(source_number)
        or mappings.get(f"{source_number}{getattr(atom, 'insertion_code', '')}")
        or {}
    )


def _add_bond(bonds: dict[Any, dict[Any, int]], left: Any, right: Any, order: int) -> None:
    """Add an undirected graph bond."""
    if order < 1 or order > 3:
        raise LocalChargePatchError(f"Unsupported local patch bond order {order!r}")
    bonds.setdefault(left, {})[right] = order
    bonds.setdefault(right, {})[left] = order


def _validate_roots_connected(
    atoms: Mapping[Any, Any], bonds: Mapping[Any, Any], roots: tuple[Any, Any]
) -> None:
    """Validate that link atoms are present and connected by the crosslink bond."""
    missing = [root for root in roots if root not in atoms]
    if missing:
        raise LocalChargePatchError(f"Local charge patch roots are missing from graph: {missing}")
    if roots[1] not in bonds.get(roots[0], {}):
        raise LocalChargePatchError(
            "Local charge patch graph is missing the product crosslink bond"
        )


def _rdkit_bond_type(order: int) -> Any:
    """Return an RDKit bond type from integer order metadata."""
    from rdkit import Chem

    if order == 1:
        return Chem.BondType.SINGLE
    if order == 2:
        return Chem.BondType.DOUBLE
    if order == 3:
        return Chem.BondType.TRIPLE
    raise LocalChargePatchError(f"Unsupported local patch bond order {order!r}")
