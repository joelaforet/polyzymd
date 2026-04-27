"""Direct OpenFF linkage bridge for the v1 conjugation deliverable."""

from __future__ import annotations

import json
import math
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation.contracts import (
    ResolvedAttachmentPlan,
    parse_pdb_atom_records,
    placed_fragment_from_resolved_plan,
)
from polyzymd.builders.conjugation.parameterization import build_formal_charge_smoke_template
from polyzymd.builders.conjugation.pdb_assembly import (
    CrosslinkedPdbAssemblyOptions,
    CrosslinkedPdbAssemblyResult,
    PdbAtomRecord,
    PlacedPolymerFragment,
    write_crosslinked_pdb,
)
from polyzymd.builders.conjugation.polymer_fragment import GeneratedPolymerFragment

_ATOM_RECORD_PREFIXES = ("ATOM", "HETATM")


class DirectLinkBondProof(BaseModel):
    """Serializable proof that the explicit protein-polymer link bond exists."""

    protein_atom_serial: int
    modifier_atom_serial: int
    protein_atom_name: str
    modifier_atom_name: str
    protein_residue_name: str
    modifier_residue_name: str
    protein_chain_id: str
    modifier_chain_id: str
    distance_angstrom: float
    conect_present: bool


class DirectOpenFFLinkageResult(BaseModel):
    """Result from direct explicit protein-polymer linkage construction.

    The OpenFF and RDKit objects are intentionally excluded from JSON output so
    that the result can always be serialized as a lightweight build summary.
    """

    model_config = {"arbitrary_types_allowed": True}

    output_dir: Path
    linked_pdb_path: Path
    summary_json_path: Path
    assembly: CrosslinkedPdbAssemblyResult
    linked_bond: DirectLinkBondProof
    finite_coordinates: bool
    atom_count: int
    removed_atom_count: int
    removed_atom_names: tuple[str, ...]
    topology: Any | None = Field(default=None, exclude=True)
    molecule: Any | None = Field(default=None, exclude=True)
    charge_template: Any | None = Field(default=None, exclude=True)
    warnings: tuple[str, ...] = Field(default_factory=tuple)
    limitations: tuple[str, ...] = Field(default_factory=tuple)

    def save_summary(self) -> Path:
        """Write the JSON-safe result summary to disk.

        Returns
        -------
        pathlib.Path
            Path to the written summary JSON file.
        """
        self.summary_json_path.write_text(
            json.dumps(self.model_dump(mode="json"), indent=2) + "\n",
            encoding="utf-8",
        )
        return self.summary_json_path


def build_direct_openff_linkage(
    *,
    protein_pdb_path: Path | str,
    modifier: GeneratedPolymerFragment | PlacedPolymerFragment,
    resolved_plan: ResolvedAttachmentPlan,
    output_dir: Path | str,
    linked_pdb_name: str = "direct_linked_conjugate.pdb",
    summary_json_name: str = "direct_openff_linkage.json",
    build_openff_topology: bool = True,
) -> DirectOpenFFLinkageResult:
    """Build a direct explicit linkage artifact and optional OpenFF topology.

    Parameters
    ----------
    protein_pdb_path : pathlib.Path or str
        Protein PDB containing the resolved protein link atom.
    modifier : GeneratedPolymerFragment or PlacedPolymerFragment
        Placed or generated polymer fragment whose coordinates are already in
        the protein reference frame.
    resolved_plan : ResolvedAttachmentPlan
        Explicit atom-level protein-polymer linkage plan.
    output_dir : pathlib.Path or str
        Directory for linked PDB and summary artifacts.
    linked_pdb_name : str, optional
        Linked PDB artifact filename, by default ``"direct_linked_conjugate.pdb"``.
    summary_json_name : str, optional
        JSON summary filename, by default ``"direct_openff_linkage.json"``.
    build_openff_topology : bool, optional
        Attempt a lazy RDKit/OpenFF conversion of the linked component graph, by
        default ``True``.

    Returns
    -------
    DirectOpenFFLinkageResult
        Linkage proof, artifacts, and optional OpenFF objects.
    """
    artifact_dir = Path(output_dir)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    linked_pdb_path = artifact_dir / linked_pdb_name
    summary_json_path = artifact_dir / summary_json_name

    placed_modifier = placed_fragment_from_resolved_plan(
        _as_placed_fragment(modifier),
        resolved_plan,
    )
    assembly = write_crosslinked_pdb(
        protein_pdb_path,
        placed_modifier,
        resolved_plan.to_pdb_linkage_attachment(),
        linked_pdb_path,
        CrosslinkedPdbAssemblyOptions(),
    )
    atoms = parse_pdb_atom_records(linked_pdb_path)
    finite_coordinates = _coordinates_are_finite(atoms)
    if not finite_coordinates:
        raise ValueError("Direct linkage bridge produced non-finite coordinates")

    conect_bonds = _conect_bonds_from_pdb(linked_pdb_path)
    linked_bond = _build_link_bond_proof(atoms, assembly.added_conect_pair, conect_bonds)
    topology = None
    molecule = None
    charge_template = None
    warnings: list[str] = []
    limitations: list[str] = []

    if build_openff_topology:
        topology, molecule, openff_warnings, openff_limitations = _try_build_openff_topology(
            protein_pdb_path,
            placed_modifier,
            resolved_plan,
            atoms,
        )
        warnings.extend(openff_warnings)
        limitations.extend(openff_limitations)
        if molecule is not None:
            charge_template = build_formal_charge_smoke_template(molecule)
            limitations.append(
                "Direct smoke charge template uses formal charges only and is not a "
                "production-quality charge model"
            )
    else:
        limitations.append("OpenFF topology construction was disabled by the caller")

    if topology is None:
        limitations.append(
            "No production-ready OpenFF topology was emitted; the linked PDB and bond proof are "
            "the authoritative v1 direct-linkage artifacts"
        )

    result = DirectOpenFFLinkageResult(
        output_dir=artifact_dir,
        linked_pdb_path=linked_pdb_path,
        summary_json_path=summary_json_path,
        assembly=assembly,
        linked_bond=linked_bond,
        finite_coordinates=finite_coordinates,
        atom_count=len(atoms),
        removed_atom_count=assembly.removed_protein_atom_count
        + assembly.removed_polymer_atom_count,
        removed_atom_names=assembly.removed_atom_names,
        topology=topology,
        molecule=molecule,
        charge_template=charge_template,
        warnings=(*assembly.warnings, *warnings),
        limitations=tuple(dict.fromkeys(limitations)),
    )
    result.save_summary()
    return result


def _as_placed_fragment(
    modifier: GeneratedPolymerFragment | PlacedPolymerFragment,
) -> PlacedPolymerFragment:
    """Return a placed fragment from supported modifier inputs."""
    if isinstance(modifier, GeneratedPolymerFragment):
        return modifier.to_placed_fragment()
    return modifier


def _coordinates_are_finite(atoms: Sequence[PdbAtomRecord]) -> bool:
    """Return whether all parsed atom coordinates are finite."""
    return all(
        math.isfinite(atom.x) and math.isfinite(atom.y) and math.isfinite(atom.z) for atom in atoms
    )


def _build_link_bond_proof(
    atoms: Sequence[PdbAtomRecord],
    added_conect_pair: tuple[int, int],
    conect_bonds: set[tuple[int, int]],
) -> DirectLinkBondProof:
    """Build a serializable proof object for the linked CONECT pair."""
    serial_to_atom = {atom.serial: atom for atom in atoms if atom.serial is not None}
    protein_serial, modifier_serial = added_conect_pair
    protein_atom = serial_to_atom[protein_serial]
    modifier_atom = serial_to_atom[modifier_serial]
    distance = math.dist(
        (protein_atom.x, protein_atom.y, protein_atom.z),
        (modifier_atom.x, modifier_atom.y, modifier_atom.z),
    )
    return DirectLinkBondProof(
        protein_atom_serial=protein_serial,
        modifier_atom_serial=modifier_serial,
        protein_atom_name=protein_atom.atom_name,
        modifier_atom_name=modifier_atom.atom_name,
        protein_residue_name=protein_atom.residue_name,
        modifier_residue_name=modifier_atom.residue_name,
        protein_chain_id=protein_atom.chain_id,
        modifier_chain_id=modifier_atom.chain_id,
        distance_angstrom=distance,
        conect_present=tuple(sorted(added_conect_pair)) in conect_bonds,
    )


def _try_build_openff_topology(
    protein_pdb_path: Path | str,
    placed_modifier: PlacedPolymerFragment,
    resolved_plan: ResolvedAttachmentPlan,
    atoms: Sequence[PdbAtomRecord],
) -> tuple[Any | None, Any | None, tuple[str, ...], tuple[str, ...]]:
    """Attempt lazy component-graph conversion of the explicit linkage."""
    try:
        from openff.toolkit import Molecule, Topology
        from openff.toolkit.utils.exceptions import OpenFFToolkitException
        from rdkit import Chem
    except ImportError as exc:
        return (
            None,
            None,
            (f"OpenFF/RDKit import unavailable: {exc}",),
            ("Install the conjugation pixi environment to materialize OpenFF objects",),
        )

    known_conversion_errors = (Chem.rdchem.MolSanitizeException, OpenFFToolkitException)
    try:
        molecule = _build_component_graph_openff_molecule(
            protein_pdb_path,
            placed_modifier,
            resolved_plan,
            atoms,
            chem_module=Chem,
            molecule_cls=Molecule,
        )
        topology = Topology.from_molecules([molecule])
    except known_conversion_errors as exc:
        return _openff_conversion_failure_result(exc)
    except (KeyError, IndexError, TypeError, ValueError) as exc:
        raise RuntimeError(
            "Unexpected direct OpenFF topology construction failure while assembling the "
            "linked component graph"
        ) from exc

    return (
        topology,
        molecule,
        (),
        (
            "Direct v1 OpenFF topology uses RDKit component graph surgery with formal-charge "
            "smoke charges; production workflows should replace this with curated residue "
            "templates and charges",
        ),
    )


def _openff_conversion_failure_result(
    exc: Exception,
) -> tuple[Any | None, Any | None, tuple[str, ...], tuple[str, ...]]:
    """Return the best-effort result for documented toolkit conversion failures.

    Parameters
    ----------
    exc : Exception
        Known RDKit or OpenFF toolkit conversion exception.

    Returns
    -------
    tuple[Any | None, Any | None, tuple[str, ...], tuple[str, ...]]
        Empty topology objects plus warning and limitation messages.
    """
    return (
        None,
        None,
        (f"OpenFF component-graph conversion failed: {exc}",),
        (
            "The direct bridge could not materialize a chemically connected OpenFF topology. "
            "Check product-state leaving atoms, polymer connectivity, and custom residue "
            "coverage before solvation.",
        ),
    )


def _build_component_graph_openff_molecule(
    protein_pdb_path: Path | str,
    placed_modifier: PlacedPolymerFragment,
    resolved_plan: ResolvedAttachmentPlan,
    linked_atoms: Sequence[PdbAtomRecord],
    *,
    chem_module: Any,
    molecule_cls: Any,
) -> Any:
    """Build a connected OpenFF molecule from valid component graphs.

    The linked PDB remains the coordinate and identity artifact. The chemistry
    graph is built from the original protein PDB and explicit polymer fragment
    connectivity so RDKit does not infer spurious proximity bonds after Packmol
    placement.
    """
    protein_atoms = parse_pdb_atom_records(protein_pdb_path)
    protein_mol = chem_module.MolFromPDBFile(
        str(protein_pdb_path),
        sanitize=True,
        removeHs=False,
        proximityBonding=True,
    )
    if protein_mol is None:
        raise ValueError(f"RDKit could not load protein PDB {protein_pdb_path}")
    if protein_mol.GetNumAtoms() != len(protein_atoms):
        raise ValueError(
            "RDKit/PDB atom-count mismatch while building direct topology: "
            f"RDKit has {protein_mol.GetNumAtoms()} atoms but PDB has {len(protein_atoms)} atoms"
        )

    protein_product, protein_old_to_new = _protein_product_rdkit_mol(
        protein_mol,
        resolved_plan,
        chem_module=chem_module,
    )
    modifier_product, modifier_old_to_new = _modifier_product_rdkit_mol(
        placed_modifier,
        chem_module=chem_module,
    )
    protein_link_index = _required_mapped_index(
        resolved_plan.protein_link_atom.atom_index,
        protein_old_to_new,
        "protein link atom",
    )
    modifier_link_index = _required_mapped_index(
        resolved_plan.modifier_link_atom.atom_index,
        modifier_old_to_new,
        "modifier link atom",
    )

    offset = protein_product.GetNumAtoms()
    combined = chem_module.CombineMols(protein_product, modifier_product)
    combined_rw = chem_module.RWMol(combined)
    combined_rw.AddBond(
        protein_link_index,
        offset + modifier_link_index,
        chem_module.BondType.SINGLE,
    )
    combined_rw.GetAtomWithIdx(protein_link_index).SetNumRadicalElectrons(0)
    combined_rw.GetAtomWithIdx(offset + modifier_link_index).SetNumRadicalElectrons(0)
    _normalize_product_formal_charges(combined_rw, protein_link_index=protein_link_index)
    combined_mol = combined_rw.GetMol()
    chem_module.SanitizeMol(combined_mol)

    if combined_mol.GetNumAtoms() != len(linked_atoms):
        raise ValueError(
            "Component graph atom count does not match linked PDB order: "
            f"graph has {combined_mol.GetNumAtoms()} atoms but linked PDB has {len(linked_atoms)}"
        )

    combined_mol = _with_linked_pdb_conformer(combined_mol, linked_atoms, chem_module=chem_module)
    molecule = molecule_cls.from_rdkit(
        combined_mol,
        allow_undefined_stereo=True,
        hydrogens_are_explicit=True,
    )
    _apply_linked_pdb_metadata(molecule, linked_atoms)
    return molecule


def _protein_product_rdkit_mol(
    protein_mol: Any,
    resolved_plan: ResolvedAttachmentPlan,
    *,
    chem_module: Any,
) -> tuple[Any, dict[int, int]]:
    """Return protein RDKit product graph and original-to-product atom map."""
    removed_indices = {
        atom.atom_index
        for atom in resolved_plan.protein_leaving_atoms
        if atom.atom_index is not None
    }
    protein_rw = chem_module.RWMol(protein_mol)
    for atom_index in sorted(removed_indices, reverse=True):
        protein_rw.RemoveAtom(atom_index)

    old_to_new = _old_to_new_index_map(protein_mol.GetNumAtoms(), removed_indices)
    link_index = _required_mapped_index(
        resolved_plan.protein_link_atom.atom_index,
        old_to_new,
        "protein link atom",
    )
    _normalize_product_formal_charges(protein_rw, protein_link_index=link_index)
    product = protein_rw.GetMol()
    chem_module.SanitizeMol(product)
    return product, old_to_new


def _modifier_product_rdkit_mol(
    placed_modifier: PlacedPolymerFragment,
    *,
    chem_module: Any,
) -> tuple[Any, dict[int, int]]:
    """Return modifier RDKit product graph and original-to-product atom map."""
    removed_keys = {
        _atom_identity(atom)
        for atom in placed_modifier.atoms
        if _is_removed_modifier_atom(atom, placed_modifier)
    }
    kept_atoms = [
        atom for atom in placed_modifier.atoms if _atom_identity(atom) not in removed_keys
    ]
    old_to_new = {
        atom.atom_index: index
        for index, atom in enumerate(kept_atoms)
        if atom.atom_index is not None
    }

    rw_mol = chem_module.RWMol()
    for atom in kept_atoms:
        element = _element_for_rdkit(atom)
        rd_atom = chem_module.Atom(element)
        rd_atom.SetNoImplicit(True)
        formal_charge = _pdb_formal_charge(atom.charge)
        if formal_charge is not None:
            rd_atom.SetFormalCharge(formal_charge)
        _set_rdkit_pdb_info(rd_atom, atom, chem_module=chem_module)
        rw_mol.AddAtom(rd_atom)

    kept_keys = {_atom_identity(atom) for atom in kept_atoms}
    for atom_1_ref, atom_2_ref, bond_order in _modifier_bond_order_records(placed_modifier):
        atom_1 = _resolve_modifier_ref(placed_modifier, atom_1_ref)
        atom_2 = _resolve_modifier_ref(placed_modifier, atom_2_ref)
        if atom_1 is None or atom_2 is None:
            continue
        if _atom_identity(atom_1) not in kept_keys or _atom_identity(atom_2) not in kept_keys:
            continue
        idx_1 = _required_mapped_index(atom_1.atom_index, old_to_new, "modifier bond atom")
        idx_2 = _required_mapped_index(atom_2.atom_index, old_to_new, "modifier bond atom")
        if rw_mol.GetBondBetweenAtoms(idx_1, idx_2) is None:
            bond_type = _rdkit_bond_type(bond_order, chem_module=chem_module)
            if bond_type == chem_module.BondType.AROMATIC:
                rw_mol.GetAtomWithIdx(idx_1).SetIsAromatic(True)
                rw_mol.GetAtomWithIdx(idx_2).SetIsAromatic(True)
            rw_mol.AddBond(idx_1, idx_2, bond_type)

    _normalize_product_formal_charges(rw_mol)
    product = rw_mol.GetMol()
    chem_module.SanitizeMol(product)
    return product, old_to_new


def _with_linked_pdb_conformer(
    rdkit_mol: Any,
    linked_atoms: Sequence[PdbAtomRecord],
    *,
    chem_module: Any,
) -> Any:
    """Return a copy of an RDKit molecule with linked-PDB coordinates."""
    rw_mol = chem_module.RWMol(rdkit_mol)
    rw_mol.RemoveAllConformers()
    mol = rw_mol.GetMol()
    conformer = chem_module.Conformer(len(linked_atoms))
    for index, atom in enumerate(linked_atoms):
        conformer.SetAtomPosition(index, (atom.x, atom.y, atom.z))
    mol.AddConformer(conformer)
    return mol


def _apply_linked_pdb_metadata(molecule: Any, linked_atoms: Sequence[PdbAtomRecord]) -> None:
    """Patch OpenFF atom names and PDB metadata from the linked artifact."""
    if int(molecule.n_atoms) != len(linked_atoms):
        raise ValueError(
            "OpenFF molecule atom count does not match linked PDB metadata: "
            f"molecule has {molecule.n_atoms} atoms but PDB has {len(linked_atoms)} atoms"
        )
    for off_atom, pdb_atom in zip(molecule.atoms, linked_atoms, strict=True):
        off_atom.name = pdb_atom.atom_name
        off_atom.metadata["chain_id"] = pdb_atom.chain_id
        off_atom.metadata["residue_name"] = pdb_atom.residue_name
        off_atom.metadata["residue_number"] = str(pdb_atom.residue_number)
        off_atom.metadata["insertion_code"] = pdb_atom.insertion_code
        if pdb_atom.serial is not None:
            off_atom.metadata["serial"] = str(pdb_atom.serial)


def _normalize_product_formal_charges(
    rw_mol: Any,
    *,
    protein_link_index: int | None = None,
) -> None:
    """Apply conservative product-state formal charge fixes before sanitization."""
    for atom in rw_mol.GetAtoms():
        if atom.GetSymbol() == "N" and atom.GetDegree() == 4 and atom.GetFormalCharge() == 0:
            atom.SetFormalCharge(1)
        if _is_terminal_oxoanion_radical(atom):
            atom.SetFormalCharge(-1)
            atom.SetNumRadicalElectrons(0)
    if protein_link_index is None:
        return
    link_atom = rw_mol.GetAtomWithIdx(protein_link_index)
    if link_atom.GetSymbol() == "N":
        link_atom.SetFormalCharge(0)


def _is_terminal_oxoanion_radical(atom: Any) -> bool:
    """Return whether an atom is a terminal sulfonate or phosphate oxoanion radical."""
    if (
        atom.GetSymbol() != "O"
        or atom.GetFormalCharge() != 0
        or atom.GetNumRadicalElectrons() == 0
        or atom.GetDegree() != 1
    ):
        return False
    neighbor = atom.GetNeighbors()[0]
    if neighbor.GetSymbol() not in {"S", "P"}:
        return False
    bond = atom.GetOwningMol().GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
    return bond is not None and abs(bond.GetBondTypeAsDouble() - 1.0) < 0.01


def _old_to_new_index_map(atom_count: int, removed_indices: set[int]) -> dict[int, int]:
    """Return an original atom index to retained atom index map."""
    old_to_new: dict[int, int] = {}
    cursor = 0
    for old_index in range(atom_count):
        if old_index in removed_indices:
            continue
        old_to_new[old_index] = cursor
        cursor += 1
    return old_to_new


def _required_mapped_index(
    atom_index: int | None,
    index_map: dict[int, int],
    label: str,
) -> int:
    """Return a mapped atom index or raise a contextual error."""
    if atom_index is None or atom_index not in index_map:
        raise ValueError(f"Cannot map {label} into the product graph")
    return index_map[atom_index]


def _atom_identity(atom: PdbAtomRecord) -> tuple[int | None, int | None, str, int, str]:
    """Return a stable source identity for a PDB atom record."""
    return (
        atom.atom_index,
        atom.serial,
        atom.atom_name,
        atom.residue_number,
        atom.insertion_code,
    )


def _is_removed_modifier_atom(atom: PdbAtomRecord, fragment: PlacedPolymerFragment) -> bool:
    """Return whether a modifier atom is selected as a product leaving atom."""
    return (
        (atom.serial is not None and atom.serial in fragment.leaving_atom_serials)
        or (atom.atom_index is not None and atom.atom_index in fragment.leaving_atom_indices)
        or atom.atom_name.upper() in {name.upper() for name in fragment.leaving_atom_names}
    )


def _resolve_modifier_ref(
    fragment: PlacedPolymerFragment,
    atom_ref: int | str,
) -> PdbAtomRecord | None:
    """Resolve a modifier bond endpoint by serial, index, or atom name."""
    if isinstance(atom_ref, str):
        matches = [atom for atom in fragment.atoms if atom.atom_name.upper() == atom_ref.upper()]
    else:
        serial_matches = [atom for atom in fragment.atoms if atom.serial == atom_ref]
        index_matches = [atom for atom in fragment.atoms if atom.atom_index == atom_ref]
        matches = serial_matches or index_matches
    unique: list[PdbAtomRecord] = []
    seen: set[tuple[int | None, int | None, str, int, str]] = set()
    for atom in matches:
        identity = _atom_identity(atom)
        if identity in seen:
            continue
        seen.add(identity)
        unique.append(atom)
    if len(unique) != 1:
        return None
    return unique[0]


def _modifier_bond_order_records(
    fragment: PlacedPolymerFragment,
) -> tuple[tuple[int | str, int | str, float], ...]:
    """Return modifier bond records with explicit or default bond orders."""
    if fragment.bond_orders:
        return fragment.bond_orders
    return tuple((atom_1, atom_2, 1.0) for atom_1, atom_2 in fragment.bonds)


def _rdkit_bond_type(bond_order: float, *, chem_module: Any) -> Any:
    """Map a numeric bond order to an RDKit bond type."""
    if abs(bond_order - 1.5) < 0.01:
        return chem_module.BondType.AROMATIC
    if abs(bond_order - 3.0) < 0.01:
        return chem_module.BondType.TRIPLE
    if abs(bond_order - 2.0) < 0.01:
        return chem_module.BondType.DOUBLE
    return chem_module.BondType.SINGLE


def _set_rdkit_pdb_info(rd_atom: Any, atom: PdbAtomRecord, *, chem_module: Any) -> None:
    """Attach PDB residue metadata to an RDKit atom."""
    info = chem_module.AtomPDBResidueInfo()
    info.SetName(_format_rdkit_atom_name(atom.atom_name))
    info.SetResidueName(atom.residue_name[:3])
    info.SetResidueNumber(atom.residue_number)
    info.SetChainId(atom.chain_id[:1])
    info.SetInsertionCode(atom.insertion_code[:1])
    info.SetIsHeteroAtom(atom.record_name == "HETATM")
    if atom.serial is not None:
        info.SetSerialNumber(atom.serial)
    rd_atom.SetMonomerInfo(info)


def _format_rdkit_atom_name(atom_name: str) -> str:
    """Format an atom name for RDKit PDB residue metadata."""
    stripped = atom_name.strip()
    if len(stripped) < 4:
        return f" {stripped:<3}"
    return f"{stripped:<4}"[:4]


def _pdb_formal_charge(charge: str) -> int | None:
    """Parse a PDB charge field into an integer formal charge."""
    normalized = charge.strip()
    if not normalized:
        return None
    if len(normalized) == 2 and normalized[0].isdigit() and normalized[1] in {"+", "-"}:
        sign = 1 if normalized[1] == "+" else -1
        return sign * int(normalized[0])
    if len(normalized) == 2 and normalized[0] in {"+", "-"} and normalized[1].isdigit():
        sign = 1 if normalized[0] == "+" else -1
        return sign * int(normalized[1])
    if normalized in {"+", "-"}:
        return 1 if normalized == "+" else -1
    return None


def _element_for_rdkit(atom: PdbAtomRecord) -> str:
    """Return a conservative element symbol for RDKit atom construction."""
    if atom.element:
        return atom.element.strip().capitalize()
    stripped = "".join(char for char in atom.atom_name if char.isalpha())
    if not stripped:
        return "C"
    candidate = stripped[0].upper()
    return "H" if candidate == "H" else candidate.capitalize()


def _conect_bonds_from_pdb(path: Path) -> set[tuple[int, int]]:
    """Parse unique serial-based CONECT bonds from a PDB artifact."""
    bonds: set[tuple[int, int]] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("CONECT"):
                continue
            source = _parse_conect_serial(line[6:11])
            if source is None:
                continue
            for start in range(11, len(line), 5):
                target = _parse_conect_serial(line[start : start + 5])
                if target is None or target == source:
                    continue
                bonds.add(tuple(sorted((source, target))))
    return bonds


def _parse_conect_serial(value: str) -> int | None:
    """Parse a CONECT serial field."""
    stripped = value.strip()
    if not stripped:
        return None
    try:
        return int(stripped)
    except ValueError:
        return None
