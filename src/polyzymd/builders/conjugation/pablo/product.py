"""Product-state Pablo residue definitions for conjugation POC products."""

from __future__ import annotations

import importlib
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any, NamedTuple

from pydantic import BaseModel, ConfigDict, Field

from polyzymd.builders.conjugation._linkage import (
    PabloCrosslinkRequirement,
    parse_pdb_atom_records,
)
from polyzymd.builders.conjugation.pablo.sdf import (
    charged_sdf_formal_charges,
    explicit_bond_order,
    read_sdf_molecules,
    sdf_bond_orders_by_atom_index,
    select_sdf_molecule,
    validate_fragment_atom_indices,
)
from polyzymd.builders.conjugation.structure.parsing import parse_pdb_conect_pairs
from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord


class ProductStatePabloDefinitionSummary(BaseModel):
    """JSON-safe summary of a generated Pablo residue definition."""

    residue_name: str
    chain_id: str | None = None
    residue_number: int | None = None
    atom_names: tuple[str, ...] = Field(default_factory=tuple)
    leaving_atom_names: tuple[str, ...] = Field(default_factory=tuple)
    bond_count: int = 0
    linking_bond: tuple[str, str] | None = None
    crosslink: tuple[str, str] | None = None


class ProductStatePabloLibrary(BaseModel):
    """Pablo cache augmented with product-state residue definitions.

    The ``charge_templates`` field is an optional cache of OpenFF molecule-level
    templates when such templates are available. It is not equivalent to Pablo
    product-state residue provenance; final conjugate parameterization relies on
    ``residue_library`` and ``definitions`` for residue-template coverage.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    residue_library: Any = Field(exclude=True)
    definitions: tuple[Any, ...] = Field(default_factory=tuple, exclude=True)
    charge_templates: tuple[Any, ...] = Field(default_factory=tuple, exclude=True)
    residue_partial_charges: tuple[Any, ...] = Field(default_factory=tuple, exclude=True)
    charge_bridge_report: Any | None = Field(default=None, exclude=True)
    charge_bridge_report_path: Path | None = None
    summaries: tuple[ProductStatePabloDefinitionSummary, ...] = Field(default_factory=tuple)
    crosslink_requirement: PabloCrosslinkRequirement
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)

    @property
    def residue_names(self) -> tuple[str, ...]:
        """Return generated residue names without duplicates."""
        seen: set[str] = set()
        names: list[str] = []
        for summary in self.summaries:
            if summary.residue_name not in seen:
                seen.add(summary.residue_name)
                names.append(summary.residue_name)
        return tuple(names)


class _PolymerExternalBond(NamedTuple):
    """One product-state bond connecting two chain-C residues."""

    left_key: tuple[str, str, int, str]
    right_key: tuple[str, str, int, str]
    left_atom: str
    right_atom: str
    order: int


class _PolymerLinkPlan(NamedTuple):
    """Pablo metadata assigned to one chain-C residue definition."""

    linking_bond: tuple[str, str, int] | None
    crosslink: tuple[str, str] | None
    extra_leaving_bonds: tuple[tuple[str, str], ...]
    atom_name_aliases: tuple[tuple[str, str], ...]


def build_product_state_pablo_library(
    product_pdb: Path | str,
    source_protein_pdb: Path | str | None = None,
    polymer_sdf: Path | str | None = None,
    generated_fragment: Any | None = None,
    resolved_plan: Any | None = None,
    product_residue_mappings: Mapping[str, Mapping[str, Any]] | None = None,
    *,
    charged_polymer_sdf: Path | str | None = None,
    pablo_module: Any | None = None,
) -> ProductStatePabloLibrary:
    """Build POC product-state Pablo residue definitions for an emitted PDB.

    The helper describes the already-modified product state directly. It does
    not call Pablo ``with_crosslink()``, because that API rejects empty leaving
    groups for product PDBs. Instead, it adds reciprocal ``crosslink`` metadata
    to the modified residue definitions and marks the reactant-state leaving
    atoms as absent leaving atoms on those definitions.
    """
    pablo = pablo_module or importlib.import_module("openff.pablo")
    residue_mod = _pablo_residue_module(pablo)
    atom_cls = getattr(pablo, "AtomDefinition", None) or residue_mod.AtomDefinition
    bond_cls = getattr(pablo, "BondDefinition", None) or residue_mod.BondDefinition
    residue_cls = getattr(pablo, "ResidueDefinition", None) or residue_mod.ResidueDefinition

    product_path = Path(product_pdb)
    source_path = Path(source_protein_pdb) if source_protein_pdb is not None else None
    product_atoms = parse_pdb_atom_records(product_path)
    source_atoms = parse_pdb_atom_records(source_path) if source_path is not None else []
    conect_pairs = parse_pdb_conect_pairs(product_path)
    requirement = _coerce_requirement(resolved_plan)
    diagnostics: list[str] = []

    protein_key = _locate_residue_key(
        product_atoms,
        residue_name=requirement.residues[0],
        atom_name=requirement.linking_atoms[0],
        resolved_atom=getattr(resolved_plan, "protein_link_atom", None),
    )
    modifier_key = _locate_residue_key(
        product_atoms,
        residue_name=requirement.residues[1],
        atom_name=requirement.linking_atoms[1],
        resolved_atom=getattr(resolved_plan, "modifier_link_atom", None),
    )

    definitions: list[Any] = []
    summaries: list[ProductStatePabloDefinitionSummary] = []

    protein_definition = _build_protein_product_definition(
        pablo,
        atom_cls,
        bond_cls,
        residue_cls,
        product_atoms=_atoms_for_key(product_atoms, protein_key),
        source_atoms=source_atoms,
        requirement=requirement,
        source_residue_name=_source_protein_residue_name(resolved_plan, source_atoms),
    )
    definitions.append(protein_definition)
    summaries.append(_summarize_definition(protein_definition, protein_key))

    polymer_residue_keys = _product_residue_keys_from_mappings(
        product_atoms,
        product_residue_mappings,
    )
    polymer_atoms = _scoped_polymer_product_atoms(product_atoms, polymer_residue_keys)
    source_residue_aliases = _source_residue_aliases(
        product_atoms,
        product_residue_mappings,
    )
    if polymer_residue_keys is not None:
        diagnostics.append(
            "Product-state definitions were scoped to one attachment's emitted chain-C residues"
        )

    fragment_bonds = _fragment_bonds(generated_fragment)
    sdf_bonds = _polymer_sdf_bonds(polymer_sdf, generated_fragment)
    if sdf_bonds:
        fragment_bonds = (*fragment_bonds, *sdf_bonds)
    product_bonds, external_polymer_bonds = _product_bonds_and_links_by_residue(
        polymer_atoms,
        conect_pairs,
        fragment_bonds=fragment_bonds,
        source_residue_aliases=source_residue_aliases,
    )
    product_formal_charges = _product_formal_charges(
        polymer_atoms,
        generated_fragment,
        charged_polymer_sdf=charged_polymer_sdf,
        source_residue_aliases=source_residue_aliases,
    )
    polymer_link_plans, polymer_link_diagnostics = _plan_polymer_external_links(
        external_polymer_bonds,
        reserved_crosslink_keys={modifier_key},
    )
    diagnostics.extend(polymer_link_diagnostics)
    if polymer_sdf is not None:
        diagnostics.append(
            f"Polymer SDF sidecar recorded for product-state library: {Path(polymer_sdf)}"
        )
    if sdf_bonds:
        diagnostics.append("Polymer SDF bond orders were applied to product-state definitions")
    if not fragment_bonds:
        diagnostics.append("No generated-fragment bond graph was available; using PDB CONECT bonds")

    for key, residue_atoms in _product_polymer_residues(polymer_atoms):
        is_modified = key == modifier_key
        link_plan = polymer_link_plans.get(key, _EMPTY_POLYMER_LINK_PLAN)
        definition = _build_pdb_residue_definition(
            atom_cls,
            bond_cls,
            residue_cls,
            residue_atoms=residue_atoms,
            bonds=product_bonds.get(key, ()),
            formal_charges=product_formal_charges.get(key, {}),
            linking_bond=link_plan.linking_bond,
            leaving_atoms=_leaving_atoms_for_residue(
                resolved_plan,
                side="modifier" if is_modified else "none",
            ),
            crosslink=(
                (requirement.linking_atoms[1], requirement.linking_atoms[0])
                if is_modified
                else link_plan.crosslink
            ),
            extra_leaving_bonds=link_plan.extra_leaving_bonds,
            atom_name_aliases=dict(link_plan.atom_name_aliases),
            bond_order=int(requirement.bond_order),
            description="PolyzyMD product-state polymer residue definition",
        )
        definitions.append(definition)
        summaries.append(_summarize_definition(definition, key))

    _validate_no_whole_polymer_collapse(summaries)
    residue_library = pablo.STD_CCD_CACHE.with_(definitions)
    diagnostics.append(
        "Product-state residue definitions were added directly; Pablo with_crosslink() was not used"
    )
    return ProductStatePabloLibrary(
        residue_library=residue_library,
        definitions=tuple(definitions),
        summaries=tuple(summaries),
        crosslink_requirement=requirement,
        diagnostics=tuple(diagnostics),
    )


def build_product_state_pablo_library_for_specs(
    product_pdb: Path | str,
    source_protein_pdb: Path | str | None,
    specs: Iterable[Any],
    *,
    pablo_module: Any | None = None,
) -> Any:
    """Build product-state Pablo residue definitions for resolved attachment specs."""
    attachment_specs = tuple(specs)
    if not attachment_specs:
        raise ValueError("Product-state Pablo library generation requires at least one spec")

    libraries = []
    for spec in attachment_specs:
        generated_fragment = _spec_generated_fragment(spec)
        sdf_path = _spec_sdf_path(spec)
        charged_sdf_path = _spec_charged_sdf_path(spec)
        _validate_spec_sidecars(spec, sdf_path=sdf_path, generated_fragment=generated_fragment)
        libraries.append(
            build_product_state_pablo_library(
                product_pdb=product_pdb,
                source_protein_pdb=source_protein_pdb,
                polymer_sdf=sdf_path,
                charged_polymer_sdf=charged_sdf_path,
                generated_fragment=generated_fragment,
                resolved_plan=getattr(spec, "resolved_plan", None),
                product_residue_mappings=_spec_product_residue_mappings(spec),
                pablo_module=pablo_module,
            )
        )

    if len(libraries) == 1:
        return libraries[0]

    pablo = pablo_module or importlib.import_module("openff.pablo")
    definitions = tuple(
        definition for library in libraries for definition in getattr(library, "definitions", ())
    )
    summaries = tuple(
        summary for library in libraries for summary in getattr(library, "summaries", ())
    )
    diagnostics = tuple(
        diagnostic for library in libraries for diagnostic in getattr(library, "diagnostics", ())
    )
    return ProductStatePabloLibrary(
        residue_library=pablo.STD_CCD_CACHE.with_(definitions),
        definitions=definitions,
        summaries=summaries,
        crosslink_requirement=attachment_specs[0].resolved_plan.pablo_crosslink_requirement,
        diagnostics=diagnostics,
    )


def _spec_generated_fragment(spec: Any) -> Any:
    """Return a product-library-compatible fragment adapter from a resolved spec."""
    generated_fragment = getattr(spec, "generated_fragment", None)
    if generated_fragment is not None:
        return generated_fragment
    fragment = getattr(spec, "fragment", None)
    to_generated = getattr(fragment, "to_generated_polymer_fragment", None)
    if callable(to_generated):
        return to_generated()
    return fragment


def _spec_sdf_path(spec: Any) -> Path | None:
    """Return the bond-order SDF sidecar recorded on a resolved spec, if present."""
    for sidecars in (
        getattr(spec, "source_sidecars", None),
        getattr(getattr(spec, "fragment", None), "sidecars", None),
    ):
        if sidecars and sidecars.get("bond_sdf") is not None:
            return Path(sidecars["bond_sdf"])
        if sidecars and sidecars.get("sdf") is not None:
            return Path(sidecars["sdf"])
    return None


def _spec_charged_sdf_path(spec: Any) -> Path | None:
    """Return the production charged SDF sidecar recorded on a resolved spec."""
    for sidecars in (
        getattr(spec, "source_sidecars", None),
        getattr(getattr(spec, "fragment", None), "sidecars", None),
    ):
        if sidecars and sidecars.get("charged_sdf") is not None:
            return Path(sidecars["charged_sdf"])
    return None


def _spec_product_residue_mappings(spec: Any) -> Mapping[str, Mapping[str, Any]] | None:
    """Return assembly residue mappings recorded for one attachment spec."""
    mappings = getattr(spec, "product_residue_mappings", None)
    return mappings or None


def _validate_spec_sidecars(
    spec: Any,
    *,
    sdf_path: Path | None,
    generated_fragment: Any,
) -> None:
    """Fail before Pablo work when a spec lacks required bond-order sidecars."""
    if sdf_path is not None:
        if not sdf_path.exists():
            raise ValueError(f"Product-state Pablo SDF sidecar does not exist: {sdf_path}")
        return

    fragment = getattr(spec, "fragment", None)
    if getattr(fragment, "source_kind", None) != "polymer":
        return
    if not _fragment_has_complete_bond_orders(generated_fragment):
        attachment_id = getattr(spec, "attachment_id", "<unknown>")
        raise ValueError(
            "Product-state Pablo library generation for polymer attachment "
            f"{attachment_id!r} requires an SDF sidecar or complete generated-fragment "
            "bond orders. Regenerate the polymer PDB/SDF pair and preserve the 'sdf' "
            "source sidecar on the attachment spec."
        )


def _fragment_has_complete_bond_orders(generated_fragment: Any) -> bool:
    """Return whether every generated-fragment bond has an explicit bond order."""
    bonds = {
        frozenset((atom1, atom2)) for atom1, atom2 in getattr(generated_fragment, "bonds", ()) or ()
    }
    if not bonds:
        return True
    bond_orders = {
        frozenset((entry[0], entry[1]))
        for entry in getattr(generated_fragment, "bond_orders", ()) or ()
        if len(entry) >= 3
    }
    return bonds.issubset(bond_orders)


def _pablo_residue_module(pablo_module: Any) -> Any:
    """Return Pablo's residue module with lazy import fallback."""
    residue_mod = getattr(pablo_module, "residue", None)
    if residue_mod is not None:
        return residue_mod
    return importlib.import_module("openff.pablo.residue")


def _coerce_requirement(resolved_plan: Any | None) -> PabloCrosslinkRequirement:
    """Extract the reactant-state Pablo requirement from a resolved plan."""
    requirement = getattr(resolved_plan, "pablo_crosslink_requirement", None)
    if requirement is None:
        raise ValueError("resolved_plan.pablo_crosslink_requirement is required")
    if isinstance(requirement, PabloCrosslinkRequirement):
        return requirement
    return PabloCrosslinkRequirement.model_validate(requirement)


def _residue_key(atom: PdbAtomRecord) -> tuple[str, str, int, str]:
    """Return the exact PDB residue key used by the product definitions."""
    return (atom.chain_id, atom.residue_name, atom.residue_number, atom.insertion_code)


def _locate_residue_key(
    atoms: list[PdbAtomRecord],
    *,
    residue_name: str,
    atom_name: str,
    resolved_atom: Any | None,
) -> tuple[str, str, int, str]:
    """Find the concrete product-PDB residue containing a crosslink endpoint."""
    candidates = [
        atom
        for atom in atoms
        if atom.residue_name.upper() == residue_name.upper()
        and atom.atom_name.upper() == atom_name.upper()
    ]
    if resolved_atom is not None:
        resolved_number = getattr(resolved_atom, "residue_number", None)
        resolved_chain = getattr(resolved_atom, "chain_id", None)
        scoped = [
            atom
            for atom in candidates
            if (resolved_chain in (None, "") or atom.chain_id == resolved_chain)
            and (resolved_number is None or atom.residue_number == resolved_number)
        ]
        if scoped:
            candidates = scoped
    if len(candidates) != 1:
        raise ValueError(
            "Expected exactly one product-state crosslink endpoint "
            f"{residue_name}:{atom_name}, found {len(candidates)}"
        )
    return _residue_key(candidates[0])


def _atoms_for_key(
    atoms: list[PdbAtomRecord],
    key: tuple[str, str, int, str],
) -> tuple[PdbAtomRecord, ...]:
    """Return product atoms for a residue key."""
    return tuple(atom for atom in atoms if _residue_key(atom) == key)


def _source_protein_residue_name(
    resolved_plan: Any | None, source_atoms: list[PdbAtomRecord]
) -> str:
    """Return the source protein residue name used to seed LYX definitions."""
    selector = getattr(getattr(resolved_plan, "contract", None), "protein_endpoint", None)
    selector = getattr(selector, "selector", None)
    if selector is not None:
        return str(getattr(selector, "residue_name", "LYS") or "LYS").upper()
    link_atom = getattr(resolved_plan, "protein_link_atom", None)
    if link_atom is not None:
        for atom in source_atoms:
            if atom.chain_id == getattr(
                link_atom, "chain_id", None
            ) and atom.residue_number == getattr(link_atom, "residue_number", None):
                return atom.residue_name.upper()
    return "LYS"


def _build_protein_product_definition(
    pablo: Any,
    atom_cls: Any,
    bond_cls: Any,
    residue_cls: Any,
    *,
    product_atoms: tuple[PdbAtomRecord, ...],
    source_atoms: list[PdbAtomRecord],
    requirement: PabloCrosslinkRequirement,
    source_residue_name: str,
) -> Any:
    """Build the modified protein residue definition from the source residue template."""
    product_names = {atom.atom_name for atom in product_atoms}
    leaving_names = set(requirement.leaving_atoms[0])
    source_template = _select_source_template(
        pablo,
        source_residue_name,
        required_names=product_names | leaving_names,
    )
    peptide_leaving_names = _template_linking_leaving_names(source_template)
    source_symbols = {atom.atom_name: atom.element for atom in source_atoms}
    atoms = []
    retained_names = product_names | leaving_names | peptide_leaving_names
    for atom in source_template.atoms:
        if atom.name not in retained_names:
            continue
        symbol = atom.symbol or source_symbols.get(atom.name) or _guess_element(atom.name)
        charge = _protein_product_atom_charge(
            atom,
            requirement=requirement,
            product_names=product_names,
        )
        atoms.append(
            atom.replace(
                symbol=symbol,
                leaving=atom.leaving or atom.name in leaving_names,
                charge=charge,
            )
        )
    existing = {atom.name for atom in atoms}
    for missing_name in sorted(retained_names - existing):
        atoms.append(
            atom_cls.with_defaults(
                missing_name,
                source_symbols.get(missing_name) or _guess_element(missing_name),
                leaving=missing_name in leaving_names or missing_name in peptide_leaving_names,
            )
        )
    bonds = tuple(
        bond
        for bond in source_template.bonds
        if bond.atom1 in {atom.name for atom in atoms}
        and bond.atom2 in {atom.name for atom in atoms}
    )
    bonds = (
        *bonds,
        *(
            bond_cls.with_defaults(requirement.linking_atoms[0], leaving_name, order=1)
            for leaving_name in sorted(leaving_names)
            if not any(leaving_name in (bond.atom1, bond.atom2) for bond in bonds)
        ),
    )
    bonds = _coalesce_leaving_bonds(
        bond_cls,
        bonds,
        linking_atom=requirement.linking_atoms[0],
        leaving_names=leaving_names,
    )
    crosslink = bond_cls.with_defaults(
        requirement.linking_atoms[0],
        requirement.linking_atoms[1],
        order=int(requirement.bond_order),
    )
    return residue_cls(
        residue_name=requirement.residues[0],
        description=f"PolyzyMD product-state {requirement.residues[0]} from {source_residue_name}",
        linking_bond=_protein_linking_bond(bond_cls, source_template, product_names),
        crosslink=crosslink,
        atoms=tuple(atoms),
        bonds=bonds,
        virtual_sites=tuple(getattr(source_template, "virtual_sites", ())),
    )


def _protein_product_atom_charge(
    atom: Any,
    *,
    requirement: PabloCrosslinkRequirement,
    product_names: set[str],
) -> int:
    """Return formal charge for a product-state protein atom definition."""
    charge = int(getattr(atom, "charge", 0) or 0)
    if atom.name != requirement.linking_atoms[0]:
        return charge
    if (atom.symbol or _guess_element(atom.name)).upper() != "N":
        return charge
    removed_hydrogens = [
        name for name in requirement.leaving_atoms[0] if _guess_element(name).upper() == "H"
    ]
    if charge > 0 and removed_hydrogens and not set(removed_hydrogens).issubset(product_names):
        return max(0, charge - 1)
    return charge


def _select_source_template(pablo: Any, residue_name: str, *, required_names: set[str]) -> Any:
    """Select the Pablo standard residue definition that best covers atom names."""
    definitions = tuple(pablo.STD_CCD_CACHE[residue_name])
    if not definitions:
        raise ValueError(f"Pablo STD_CCD_CACHE has no definitions for {residue_name}")
    return max(
        definitions,
        key=lambda definition: len(required_names & {atom.name for atom in definition.atoms}),
    )


def _template_linking_leaving_names(source_template: Any) -> set[str]:
    """Return absent peptide-neighbor leaving atoms from a Pablo protein template."""
    if getattr(source_template, "linking_bond", None) is None:
        return set()
    return {
        atom.name
        for atom in getattr(source_template, "atoms", ())
        if getattr(atom, "leaving", False)
    }


def _fragment_bonds(generated_fragment: Any | None) -> tuple[tuple[Any, Any, int], ...]:
    """Return generated-fragment bonds with product-mappable atom descriptors."""
    if generated_fragment is None:
        return ()
    atom_lookup = _fragment_atom_lookup(generated_fragment)
    unique_atom_names = _unique_fragment_atom_names(generated_fragment)
    order_lookup: dict[frozenset[Any], int] = {}
    for entry in getattr(generated_fragment, "bond_orders", ()) or ():
        if len(entry) >= 3:
            order_lookup[frozenset((entry[0], entry[1]))] = explicit_bond_order(
                entry[2], source="generated polymer fragment"
            )
    bonds = []
    for atom1, atom2 in getattr(generated_fragment, "bonds", ()) or ():
        resolved1 = (
            _fragment_atom_descriptor(atom_lookup.get(atom1), unique_atom_names=unique_atom_names)
            or atom1
        )
        resolved2 = (
            _fragment_atom_descriptor(atom_lookup.get(atom2), unique_atom_names=unique_atom_names)
            or atom2
        )
        bonds.append((resolved1, resolved2, order_lookup.get(frozenset((atom1, atom2)), 1)))
    return tuple(bonds)


def _polymer_sdf_bonds(
    polymer_sdf: Path | str | None,
    generated_fragment: Any | None,
) -> tuple[tuple[Any, Any, int], ...]:
    """Return SDF bond orders mapped to product-PDB atom descriptors."""
    if polymer_sdf is None:
        return ()
    if generated_fragment is None:
        raise ValueError(
            "polymer_sdf was provided but generated_fragment is unavailable; product-state "
            "Pablo definitions need the generated fragment atom indices to map SDF bond orders "
            "onto product PDB atom names."
        )
    sdf_path = Path(polymer_sdf)
    if not sdf_path.exists():
        raise ValueError(f"Polymer SDF sidecar does not exist: {sdf_path}")

    molecules = read_sdf_molecules(sdf_path, source_label="Polymer SDF")
    fragment_atoms = tuple(getattr(generated_fragment, "atoms", ()) or ())
    atoms_by_index = validate_fragment_atom_indices(
        fragment_atoms,
        source_label="Generated polymer fragment atoms for SDF bond orders",
    )
    if not atoms_by_index:
        raise ValueError(
            "Generated polymer fragment atoms do not carry atom_index values, so SDF bond "
            "orders cannot be mapped onto product-state Pablo definitions."
        )
    mol = select_sdf_molecule(
        molecules,
        expected_atoms=len(fragment_atoms),
        sdf_path=sdf_path,
        source_label="Polymer SDF sidecar",
        prefer_most_bonds=True,
    )
    sdf_bonds = sdf_bond_orders_by_atom_index(
        mol,
        fragment_atoms=fragment_atoms,
        sdf_path=sdf_path,
        source_label="Polymer SDF sidecar",
    )
    unique_atom_names = _unique_fragment_atom_names(generated_fragment)
    bonds: list[tuple[Any, Any, int]] = []
    for begin_index, end_index, order in sdf_bonds:
        atom1 = _fragment_atom_descriptor(
            atoms_by_index.get(begin_index), unique_atom_names=unique_atom_names
        )
        atom2 = _fragment_atom_descriptor(
            atoms_by_index.get(end_index), unique_atom_names=unique_atom_names
        )
        if atom1 is None or atom2 is None:
            raise ValueError(
                "Polymer SDF bond endpoint could not be mapped to a generated-fragment atom: "
                f"SDF atoms {begin_index + 1}-{end_index + 1} in "
                f"{sdf_path}. Regenerate the polymer PDB/SDF pair from the same source molecule."
            )
        bonds.append((atom1, atom2, order))
    return tuple(bonds)


def _product_residue_keys_from_mappings(
    product_atoms: list[PdbAtomRecord],
    product_residue_mappings: Mapping[str, Mapping[str, Any]] | None,
) -> set[tuple[str, str, int, str]] | None:
    """Return concrete product residue keys emitted for one attachment."""
    if not product_residue_mappings:
        return None
    keys: set[tuple[str, str, int, str]] = set()
    for mapping in product_residue_mappings.values():
        target_chain = str(mapping.get("target_chain", "C"))
        target_number = _parse_int(mapping.get("target_residue_number"))
        if target_number is None:
            continue
        keys.update(
            _residue_key(atom)
            for atom in product_atoms
            if atom.chain_id == target_chain and atom.residue_number == target_number
        )
    return keys


def _scoped_polymer_product_atoms(
    product_atoms: list[PdbAtomRecord],
    polymer_residue_keys: set[tuple[str, str, int, str]] | None,
) -> list[PdbAtomRecord]:
    """Return atoms with chain-C residues limited to one attachment when known."""
    if polymer_residue_keys is None:
        return product_atoms
    return [
        atom
        for atom in product_atoms
        if atom.chain_id != "C" or _residue_key(atom) in polymer_residue_keys
    ]


def _source_residue_aliases(
    product_atoms: list[PdbAtomRecord],
    product_residue_mappings: Mapping[str, Mapping[str, Any]] | None,
) -> dict[tuple[str, int, str], tuple[int, str]]:
    """Map emitted product residues back to source fragment residue numbers."""
    if not product_residue_mappings:
        return {}
    aliases: dict[tuple[str, int, str], tuple[int, str]] = {}
    for source_key, mapping in product_residue_mappings.items():
        source_number = _parse_source_residue_number(source_key, mapping)
        if source_number is None:
            continue
        target_chain = str(mapping.get("target_chain", "C"))
        target_number = _parse_int(mapping.get("target_residue_number"))
        if target_number is None:
            continue
        insertion_code = ""
        for atom in product_atoms:
            if atom.chain_id == target_chain and atom.residue_number == target_number:
                insertion_code = atom.insertion_code
                break
        aliases[(target_chain, target_number, insertion_code)] = (source_number, "")
    return aliases


def _parse_source_residue_number(
    source_key: str,
    mapping: Mapping[str, Any],
) -> int | None:
    """Return the source fragment residue number from an assembly mapping."""
    mapped = _parse_int(mapping.get("source_residue_number"))
    if mapped is not None:
        return mapped
    digits = "".join(char for char in str(source_key) if char.isdigit())
    return _parse_int(digits)


def _product_formal_charges(
    product_atoms: list[PdbAtomRecord],
    generated_fragment: Any | None,
    *,
    charged_polymer_sdf: Path | str | None = None,
    source_residue_aliases: Mapping[tuple[str, int, str], tuple[int, str]] | None = None,
) -> dict[tuple[str, str, int, str], dict[str, int]]:
    """Map non-zero generated-fragment formal charges onto product PDB residues."""
    if generated_fragment is None:
        return {}
    charged_formal_charges = _charged_sdf_formal_charges(
        charged_polymer_sdf,
        generated_fragment=generated_fragment,
    )
    unique_atom_names = _unique_fragment_atom_names(generated_fragment)
    product_lookup = _product_lookup(
        product_atoms,
        source_residue_aliases=source_residue_aliases,
    )
    charges: dict[tuple[str, str, int, str], dict[str, int]] = defaultdict(dict)
    for fallback_index, atom in enumerate(getattr(generated_fragment, "atoms", ()) or ()):
        raw_atom_index = getattr(atom, "atom_index", None)
        atom_index = int(raw_atom_index) if raw_atom_index is not None else fallback_index
        charge = charged_formal_charges.get(atom_index, getattr(atom, "formal_charge", None))
        if charge in (None, 0):
            continue
        descriptor = _fragment_atom_descriptor(atom, unique_atom_names=unique_atom_names)
        if descriptor is None:
            continue
        product_atom = product_lookup.get(descriptor)
        if product_atom is None:
            continue
        charges[_residue_key(product_atom)][product_atom.atom_name] = int(charge)
    return {key: dict(value) for key, value in charges.items()}


def _charged_sdf_formal_charges(
    charged_polymer_sdf: Path | str | None,
    *,
    generated_fragment: Any,
) -> dict[int, int]:
    """Return non-zero formal charges from a validated charged SDF sidecar."""
    if charged_polymer_sdf is None:
        return {}
    fragment_atoms = tuple(getattr(generated_fragment, "atoms", ()) or ())
    return charged_sdf_formal_charges(charged_polymer_sdf, fragment_atoms=fragment_atoms)


def _explicit_bond_order(value: Any, *, source: str) -> int:
    """Return a supported integer bond order or raise an actionable error."""
    return explicit_bond_order(value, source=source)


def _fragment_atom_lookup(generated_fragment: Any) -> dict[Any, Any]:
    """Build lookups for generated-fragment atom references."""
    atoms = tuple(getattr(generated_fragment, "atoms", ()) or ())
    atom_names = [getattr(atom, "atom_name", None) for atom in atoms]
    unique_atom_names = {name for name in atom_names if name and atom_names.count(name) == 1}
    lookup: dict[Any, Any] = {}
    for atom in atoms:
        serial = getattr(atom, "serial", None)
        atom_index = getattr(atom, "atom_index", None)
        atom_name = getattr(atom, "atom_name", None)
        if serial is not None:
            lookup.setdefault(serial, atom)
        if atom_index is not None:
            lookup.setdefault(atom_index, atom)
        if atom_name in unique_atom_names:
            lookup.setdefault(atom_name, atom)
    return lookup


def _unique_fragment_atom_names(generated_fragment: Any | None) -> set[str]:
    """Return atom names that uniquely identify generated-fragment atoms."""
    if generated_fragment is None:
        return set()
    atom_names = [
        getattr(atom, "atom_name", None) for atom in getattr(generated_fragment, "atoms", ()) or ()
    ]
    return {name for name in atom_names if name and atom_names.count(name) == 1}


def _fragment_atom_descriptor(
    atom: Any | None,
    *,
    unique_atom_names: set[str] | None = None,
) -> str | tuple[int | None, str] | None:
    """Return a product-PDB lookup descriptor for a generated-fragment atom."""
    if atom is None:
        return None
    atom_name = getattr(atom, "atom_name", None)
    if not atom_name:
        return None
    if unique_atom_names is not None and atom_name in unique_atom_names:
        return atom_name
    return (getattr(atom, "residue_number", None), atom_name)


def _protein_linking_bond(
    bond_cls: Any, source_template: Any, product_names: set[str]
) -> Any | None:
    """Return the peptide linking bond for a modified protein residue."""
    linking_bond = getattr(source_template, "linking_bond", None)
    if linking_bond is not None:
        return linking_bond
    if {"N", "C"}.issubset(product_names):
        return bond_cls.with_defaults("C", "N", order=1)
    return None


def _product_bonds_and_links_by_residue(
    product_atoms: list[PdbAtomRecord],
    conect_pairs: tuple[tuple[int, int], ...],
    *,
    fragment_bonds: tuple[tuple[Any, Any, int], ...],
    source_residue_aliases: Mapping[tuple[str, int, str], tuple[int, str]] | None = None,
) -> tuple[
    dict[tuple[str, str, int, str], tuple[tuple[str, str, int], ...]],
    tuple[_PolymerExternalBond, ...],
]:
    """Collect intra-residue bonds and chain-C inter-residue product bonds."""
    serial_to_atom = {atom.serial: atom for atom in product_atoms if atom.serial is not None}
    by_residue: dict[tuple[str, str, int, str], dict[tuple[str, str], int]] = defaultdict(dict)
    external_bonds: set[_PolymerExternalBond] = set()
    for serial1, serial2 in conect_pairs:
        atom1 = serial_to_atom.get(serial1)
        atom2 = serial_to_atom.get(serial2)
        if atom1 is None or atom2 is None:
            continue
        if _residue_key(atom1) == _residue_key(atom2):
            key = _residue_key(atom1)
            by_residue[key].setdefault(_ordered_bond_key(atom1.atom_name, atom2.atom_name), 1)
        else:
            _record_external_polymer_bond(external_bonds, atom1, atom2, 1)

    product_lookup = _product_lookup(
        product_atoms,
        source_residue_aliases=source_residue_aliases,
    )
    for raw1, raw2, order in fragment_bonds:
        atom1 = product_lookup.get(raw1)
        atom2 = product_lookup.get(raw2)
        if atom1 is None or atom2 is None:
            continue
        if _residue_key(atom1) == _residue_key(atom2):
            by_residue[_residue_key(atom1)][
                _ordered_bond_key(atom1.atom_name, atom2.atom_name)
            ] = order
        else:
            _record_external_polymer_bond(external_bonds, atom1, atom2, order)

    return {
        key: tuple((atom1, atom2, order) for (atom1, atom2), order in sorted(bonds.items()))
        for key, bonds in by_residue.items()
    }, tuple(
        sorted(
            external_bonds,
            key=lambda bond: (
                bond.left_key[0],
                bond.left_key[2],
                bond.left_key[3],
                bond.right_key[2],
                bond.right_key[3],
                bond.left_atom,
                bond.right_atom,
            ),
        )
    )


def _record_external_polymer_bond(
    external_bonds: set[_PolymerExternalBond],
    atom1: PdbAtomRecord,
    atom2: PdbAtomRecord,
    order: int,
) -> None:
    """Record a chain-C inter-residue bond without collapsing monomers."""
    if atom1.chain_id != "C" or atom2.chain_id != "C" or atom1.chain_id != atom2.chain_id:
        return
    left, right = sorted(
        (atom1, atom2), key=lambda atom: (atom.residue_number, atom.insertion_code)
    )
    if left.residue_number == right.residue_number and left.insertion_code == right.insertion_code:
        return
    external_bonds.add(
        _PolymerExternalBond(
            left_key=_residue_key(left),
            right_key=_residue_key(right),
            left_atom=left.atom_name,
            right_atom=right.atom_name,
            order=int(order),
        )
    )


_POLYMER_LINK_EXIT_ATOM = "POU"
_POLYMER_LINK_ENTRY_ATOM = "PIN"
_EMPTY_POLYMER_LINK_PLAN = _PolymerLinkPlan(None, None, (), ())


def _plan_polymer_external_links(
    external_bonds: tuple[_PolymerExternalBond, ...],
    *,
    reserved_crosslink_keys: set[tuple[str, str, int, str]],
) -> tuple[dict[tuple[str, str, int, str], _PolymerLinkPlan], list[str]]:
    """Assign chain-C external bonds to uniform Pablo polymer link metadata."""
    link_bonds: dict[tuple[str, str, int, str], tuple[str, str, int]] = {}
    crosslinks: dict[tuple[str, str, int, str], tuple[str, str]] = {}
    extra_leaving_bonds: dict[tuple[str, str, int, str], list[tuple[str, str]]] = defaultdict(list)
    aliases: dict[tuple[str, str, int, str], dict[str, str]] = defaultdict(dict)
    diagnostics: list[str] = []
    leaving_index = 1

    polymer_link_edges = _polymer_link_edge_keys(external_bonds, diagnostics)
    linking_edges = []
    for bond in external_bonds:
        if _polymer_external_edge_key(bond) in polymer_link_edges:
            linking_edges.append(bond)
            continue
        if bond.left_key in reserved_crosslink_keys or bond.right_key in reserved_crosslink_keys:
            diagnostics.append(
                "Nonsequential product polymer bond touches a residue whose Pablo crosslink "
                f"slot is reserved: {bond.left_key}:{bond.left_atom}-{bond.right_key}:{bond.right_atom}"
            )
            continue
        if bond.left_key in crosslinks or bond.right_key in crosslinks:
            diagnostics.append(
                "Multiple product polymer crosslinks requested for one residue; "
                f"skipping {bond.left_key}:{bond.left_atom}-{bond.right_key}:{bond.right_atom}"
            )
            continue
        crosslinks[bond.left_key] = (bond.left_atom, bond.right_atom)
        crosslinks[bond.right_key] = (bond.right_atom, bond.left_atom)
        leaving_index = _add_link_leaving_bonds(
            extra_leaving_bonds,
            bond,
            leaving_index,
            left_prefix="X",
            right_prefix="Y",
        )

    for bond in _orient_polymer_external_bonds(tuple(linking_edges), diagnostics):
        link = (_POLYMER_LINK_EXIT_ATOM, _POLYMER_LINK_ENTRY_ATOM, bond.order)
        if not _set_uniform_polymer_link(link_bonds, bond.left_key, link):
            diagnostics.append(
                "Conflicting product polymer link order for residue " f"{bond.left_key}"
            )
            continue
        if not _set_uniform_polymer_link(link_bonds, bond.right_key, link):
            diagnostics.append(
                "Conflicting product polymer link order for residue " f"{bond.right_key}"
            )
            continue
        aliases[bond.left_key][bond.left_atom] = _POLYMER_LINK_EXIT_ATOM
        aliases[bond.right_key][bond.right_atom] = _POLYMER_LINK_ENTRY_ATOM
        leaving_index = _add_link_leaving_bonds(
            extra_leaving_bonds,
            _PolymerExternalBond(
                left_key=bond.left_key,
                right_key=bond.right_key,
                left_atom=_POLYMER_LINK_EXIT_ATOM,
                right_atom=_POLYMER_LINK_ENTRY_ATOM,
                order=bond.order,
            ),
            leaving_index,
            left_prefix="L",
            right_prefix="R",
        )

    keys = set(link_bonds) | set(crosslinks) | set(extra_leaving_bonds) | set(aliases)
    return {
        key: _PolymerLinkPlan(
            linking_bond=link_bonds.get(key),
            crosslink=crosslinks.get(key),
            extra_leaving_bonds=tuple(extra_leaving_bonds.get(key, ())),
            atom_name_aliases=tuple(sorted(aliases.get(key, {}).items())),
        )
        for key in keys
    }, diagnostics


def _polymer_link_edge_keys(
    external_bonds: tuple[_PolymerExternalBond, ...],
    diagnostics: list[str],
) -> set[frozenset[tuple[str, str, int, str]]]:
    """Identify polymer backbone links from chain-C connectivity, not residue numbers."""
    adjacency: dict[
        tuple[str, str, int, str],
        list[tuple[tuple[str, str, int, str], str, str, int]],
    ] = defaultdict(list)
    bonds_by_component_edge: dict[frozenset[tuple[str, str, int, str]], _PolymerExternalBond] = {}
    for bond in external_bonds:
        adjacency[bond.left_key].append(
            (bond.right_key, bond.left_atom, bond.right_atom, bond.order)
        )
        adjacency[bond.right_key].append(
            (bond.left_key, bond.right_atom, bond.left_atom, bond.order)
        )
        bonds_by_component_edge[_polymer_external_edge_key(bond)] = bond

    link_edges: set[frozenset[tuple[str, str, int, str]]] = set()
    unvisited = set(adjacency)
    while unvisited:
        component = _polymer_component(next(iter(unvisited)), adjacency)
        unvisited.difference_update(component)
        component_edges = {
            edge_key for edge_key in bonds_by_component_edge if edge_key.issubset(component)
        }
        if not component_edges:
            continue
        if all(len(adjacency[key]) <= 2 for key in component):
            link_edges.update(component_edges)
            if any(
                abs(bond.left_key[2] - bond.right_key[2]) != 1
                for edge_key, bond in bonds_by_component_edge.items()
                if edge_key in component_edges
            ):
                diagnostics.append(
                    "Product polymer links were assigned from chain-C connectivity rather "
                    "than residue-number adjacency"
                )
            continue

        diagnostics.append(
            "Branched product polymer connectivity detected; residue-number-adjacent "
            "bonds were treated as polymer links"
        )
        link_edges.update(
            edge_key
            for edge_key, bond in bonds_by_component_edge.items()
            if edge_key in component_edges and abs(bond.left_key[2] - bond.right_key[2]) == 1
        )
    return link_edges


def _polymer_external_edge_key(
    bond: _PolymerExternalBond,
) -> frozenset[tuple[str, str, int, str]]:
    """Return an orientation-independent key for an external polymer bond."""
    return frozenset((bond.left_key, bond.right_key))


def _orient_polymer_external_bonds(
    external_bonds: tuple[_PolymerExternalBond, ...],
    diagnostics: list[str],
) -> tuple[_PolymerExternalBond, ...]:
    """Orient polymer residue bonds by graph connectivity, not residue number."""
    adjacency: dict[
        tuple[str, str, int, str],
        list[tuple[tuple[str, str, int, str], str, str, int]],
    ] = defaultdict(list)
    for bond in external_bonds:
        adjacency[bond.left_key].append(
            (bond.right_key, bond.left_atom, bond.right_atom, bond.order)
        )
        adjacency[bond.right_key].append(
            (bond.left_key, bond.right_atom, bond.left_atom, bond.order)
        )

    for key, neighbors in adjacency.items():
        if len(neighbors) > 2:
            diagnostics.append(
                "Product polymer residue has more than two inter-residue bonds; "
                f"Pablo polymer-link metadata may be incomplete for {key}"
            )

    oriented: list[_PolymerExternalBond] = []
    visited_edges: set[frozenset[tuple[str, str, int, str]]] = set()
    unvisited = set(adjacency)
    while unvisited:
        component = _polymer_component(next(iter(unvisited)), adjacency)
        unvisited.difference_update(component)
        endpoints = sorted(
            (key for key in component if len(adjacency[key]) == 1),
            key=_polymer_residue_sort_key,
        )
        start = endpoints[0] if endpoints else sorted(component, key=_polymer_residue_sort_key)[0]
        previous: tuple[str, str, int, str] | None = None
        current = start
        while True:
            next_entries = [
                entry
                for entry in adjacency[current]
                if entry[0] != previous and frozenset((current, entry[0])) not in visited_edges
            ]
            if not next_entries:
                break
            neighbor, current_atom, neighbor_atom, order = sorted(
                next_entries, key=lambda item: _polymer_residue_sort_key(item[0])
            )[0]
            visited_edges.add(frozenset((current, neighbor)))
            oriented.append(
                _PolymerExternalBond(
                    left_key=current,
                    right_key=neighbor,
                    left_atom=current_atom,
                    right_atom=neighbor_atom,
                    order=order,
                )
            )
            previous, current = current, neighbor
            if current == start:
                break
    return tuple(oriented)


def _polymer_residue_sort_key(key: tuple[str, str, int, str]) -> tuple[int, str, str, str]:
    """Sort chain-C residue keys by emitted residue order."""
    return (key[2], key[3], key[0], key[1])


def _polymer_component(
    start: tuple[str, str, int, str],
    adjacency: dict[
        tuple[str, str, int, str],
        list[tuple[tuple[str, str, int, str], str, str, int]],
    ],
) -> set[tuple[str, str, int, str]]:
    """Return the connected component containing ``start``."""
    seen: set[tuple[str, str, int, str]] = set()
    stack = [start]
    while stack:
        key = stack.pop()
        if key in seen:
            continue
        seen.add(key)
        stack.extend(neighbor for neighbor, *_rest in adjacency[key] if neighbor not in seen)
    return seen


def _set_uniform_polymer_link(
    link_bonds: dict[tuple[str, str, int, str], tuple[str, str, int]],
    key: tuple[str, str, int, str],
    link: tuple[str, str, int],
) -> bool:
    """Set or validate the uniform polymer link for one residue."""
    existing = link_bonds.get(key)
    if existing is not None:
        return existing == link
    link_bonds[key] = link
    return True


def _add_link_leaving_bonds(
    extra_leaving_bonds: dict[tuple[str, str, int, str], list[tuple[str, str]]],
    bond: _PolymerExternalBond,
    leaving_index: int,
    *,
    left_prefix: str,
    right_prefix: str,
) -> int:
    """Add absent leaving atoms needed for Pablo product-state link matching."""
    left_leaver = f"{left_prefix}{leaving_index % 1000:03d}"
    right_leaver = f"{right_prefix}{leaving_index % 1000:03d}"
    extra_leaving_bonds[bond.left_key].append((bond.left_atom, left_leaver))
    extra_leaving_bonds[bond.right_key].append((bond.right_atom, right_leaver))
    return leaving_index + 1


def _product_lookup(
    product_atoms: list[PdbAtomRecord],
    *,
    source_residue_aliases: Mapping[tuple[str, int, str], tuple[int, str]] | None = None,
) -> dict[Any, PdbAtomRecord]:
    """Build flexible product atom lookups for generated-fragment bond references."""
    lookup: dict[Any, PdbAtomRecord] = {}
    atom_name_counts = Counter(atom.atom_name for atom in product_atoms)
    for atom in product_atoms:
        if atom.serial is not None:
            lookup.setdefault(atom.serial, atom)
        if atom.atom_index is not None:
            lookup.setdefault(atom.atom_index, atom)
        if atom_name_counts[atom.atom_name] == 1:
            lookup.setdefault(atom.atom_name, atom)
        lookup.setdefault((atom.residue_number, atom.atom_name), atom)
        lookup.setdefault((atom.chain_id, atom.residue_number, atom.atom_name), atom)
        if source_residue_aliases:
            alias = source_residue_aliases.get(
                (atom.chain_id, atom.residue_number, atom.insertion_code)
            )
            if alias is not None:
                source_number, source_insertion = alias
                lookup.setdefault((source_number, atom.atom_name), atom)
                lookup.setdefault((source_number, source_insertion, atom.atom_name), atom)
                lookup.setdefault((atom.chain_id, source_number, atom.atom_name), atom)
    return lookup


def _ordered_bond(atom1: str, atom2: str, order: int) -> tuple[str, str, int]:
    """Return a stable unordered bond tuple with order metadata."""
    left, right = sorted((atom1, atom2))
    return (left, right, order)


def _ordered_bond_key(atom1: str, atom2: str) -> tuple[str, str]:
    """Return a stable unordered bond key without order metadata."""
    left, right = sorted((atom1, atom2))
    return (left, right)


def _product_polymer_residues(
    product_atoms: list[PdbAtomRecord],
) -> tuple[tuple[tuple[str, str, int, str], tuple[PdbAtomRecord, ...]], ...]:
    """Return chain-C product residues without collapsing the polymer."""
    grouped: dict[tuple[str, str, int, str], list[PdbAtomRecord]] = defaultdict(list)
    for atom in product_atoms:
        if atom.chain_id == "C":
            grouped[_residue_key(atom)].append(atom)
    return tuple((key, tuple(atoms)) for key, atoms in grouped.items())


def _leaving_atoms_for_residue(
    resolved_plan: Any | None, *, side: str
) -> tuple[PdbAtomRecord, ...]:
    """Return reactant-state leaving atoms for one product residue side."""
    if side == "modifier":
        return tuple(getattr(resolved_plan, "modifier_leaving_atoms", ()) or ())
    if side == "protein":
        return tuple(getattr(resolved_plan, "protein_leaving_atoms", ()) or ())
    return ()


def _build_pdb_residue_definition(
    atom_cls: Any,
    bond_cls: Any,
    residue_cls: Any,
    *,
    residue_atoms: tuple[PdbAtomRecord, ...],
    bonds: tuple[tuple[str, str, int], ...],
    formal_charges: Mapping[str, int],
    linking_bond: tuple[str, str, int] | None,
    leaving_atoms: tuple[PdbAtomRecord, ...],
    crosslink: tuple[str, str] | None,
    extra_leaving_bonds: tuple[tuple[str, str], ...],
    atom_name_aliases: dict[str, str],
    bond_order: int,
    description: str,
) -> Any:
    """Build one non-protein product residue definition from PDB atoms and bonds."""
    residue_name = residue_atoms[0].residue_name
    atoms = []
    for atom in residue_atoms:
        canonical_name = atom_name_aliases.get(atom.atom_name, atom.atom_name)
        synonyms = (atom.atom_name,) if canonical_name != atom.atom_name else ()
        atoms.append(
            atom_cls.with_defaults(
                canonical_name,
                atom.element or _guess_element(atom.atom_name),
                synonyms=synonyms,
                charge=_formal_charge_for_atom(atom, formal_charges),
            )
        )
    existing = {atom.name for atom in atoms}
    bond_defs = [
        bond_cls.with_defaults(
            atom_name_aliases.get(atom1, atom1),
            atom_name_aliases.get(atom2, atom2),
            order=order,
        )
        for atom1, atom2, order in bonds
    ]
    for leaving_atom in leaving_atoms:
        if leaving_atom.atom_name not in existing:
            atoms.append(
                atom_cls.with_defaults(
                    leaving_atom.atom_name,
                    leaving_atom.element or _guess_element(leaving_atom.atom_name),
                    leaving=True,
                    charge=_formal_charge_for_atom(leaving_atom, formal_charges),
                )
            )
            existing.add(leaving_atom.atom_name)
        else:
            atoms = [
                atom.replace(leaving=True) if atom.name == leaving_atom.atom_name else atom
                for atom in atoms
            ]
        link_atom = crosslink[0] if crosslink is not None else None
        if link_atom and not any(
            leaving_atom.atom_name in (bond.atom1, bond.atom2) for bond in bond_defs
        ):
            bond_defs.append(
                bond_cls.with_defaults(link_atom, leaving_atom.atom_name, order=bond_order)
            )
    for link_atom, leaving_name in extra_leaving_bonds:
        if leaving_name not in existing:
            atoms.append(atom_cls.with_defaults(leaving_name, "H", leaving=True))
            existing.add(leaving_name)
        bond_defs.append(bond_cls.with_defaults(link_atom, leaving_name, order=1))
    crosslink_def = (
        bond_cls.with_defaults(crosslink[0], crosslink[1], order=bond_order)
        if crosslink is not None
        else None
    )
    linking_bond_def = (
        bond_cls.with_defaults(linking_bond[0], linking_bond[1], order=linking_bond[2])
        if linking_bond is not None
        else None
    )
    if crosslink is not None:
        bond_defs = list(
            _coalesce_leaving_bonds(
                bond_cls,
                tuple(bond_defs),
                linking_atom=crosslink[0],
                leaving_names={atom.name for atom in atoms if getattr(atom, "leaving", False)},
            )
        )
    return residue_cls(
        residue_name=residue_name,
        description=description,
        linking_bond=linking_bond_def,
        crosslink=crosslink_def,
        atoms=tuple(atoms),
        bonds=tuple(bond_defs),
        virtual_sites=(),
    )


def _coalesce_leaving_bonds(
    bond_cls: Any,
    bonds: tuple[Any, ...],
    *,
    linking_atom: str,
    leaving_names: set[str],
) -> tuple[Any, ...]:
    """Make multi-atom leaving groups a single fragment from the linking atom.

    Pablo permits at most one leaving atom directly bonded to a crosslinking
    atom. Some product-state mechanisms remove multiple hydrogens from that
    atom; for matching purposes they must be represented as one absent leaving
    fragment rather than multiple direct leaving neighbors.
    """
    direct_leavers = sorted(
        bond.atom2 if bond.atom1 == linking_atom else bond.atom1
        for bond in bonds
        if linking_atom in (bond.atom1, bond.atom2)
        and (bond.atom2 if bond.atom1 == linking_atom else bond.atom1) in leaving_names
    )
    if len(direct_leavers) <= 1:
        return bonds

    retained = direct_leavers[0]
    rewritten = [
        bond
        for bond in bonds
        if not (
            linking_atom in (bond.atom1, bond.atom2)
            and (bond.atom2 if bond.atom1 == linking_atom else bond.atom1) in direct_leavers[1:]
        )
    ]
    previous = retained
    for leaving_name in direct_leavers[1:]:
        rewritten.append(bond_cls.with_defaults(previous, leaving_name, order=1))
        previous = leaving_name
    return tuple(rewritten)


def _summarize_definition(
    definition: Any,
    key: tuple[str, str, int, str],
) -> ProductStatePabloDefinitionSummary:
    """Build a JSON-safe summary for one Pablo definition."""
    crosslink = getattr(definition, "crosslink", None)
    linking_bond = getattr(definition, "linking_bond", None)
    return ProductStatePabloDefinitionSummary(
        residue_name=str(getattr(definition, "residue_name", key[1])),
        chain_id=key[0],
        residue_number=key[2],
        atom_names=tuple(atom.name for atom in getattr(definition, "atoms", ())),
        leaving_atom_names=tuple(
            atom.name
            for atom in getattr(definition, "atoms", ())
            if getattr(atom, "leaving", False)
        ),
        bond_count=len(tuple(getattr(definition, "bonds", ()))),
        linking_bond=(linking_bond.atom1, linking_bond.atom2) if linking_bond is not None else None,
        crosslink=(crosslink.atom1, crosslink.atom2) if crosslink is not None else None,
    )


def _validate_no_whole_polymer_collapse(
    summaries: list[ProductStatePabloDefinitionSummary],
) -> None:
    """Reject the old anti-pattern of describing all chain-C atoms as one POLY residue.

    Direct SMILES-moiety requests can legitimately contain a single chain-C
    residue, such as one GlcNAc/NAG unit for an N-glycosylation validation example.
    """
    if any(summary.residue_name == "POLY" for summary in summaries):
        raise ValueError("Product-state Pablo definitions must not collapse chain C to POLY")


def _parse_int(value: Any) -> int | None:
    """Parse an integer value or return ``None``."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _parse_formal_charge(charge: str | None) -> int:
    """Parse the fixed-width PDB formal-charge field."""
    if not charge:
        return 0
    text = charge.strip()
    if len(text) == 2 and text[0].isdigit() and text[1] in {"+", "-"}:
        magnitude = int(text[0])
        return magnitude if text[1] == "+" else -magnitude
    return 0


def _formal_charge_for_atom(atom: Any, formal_charges: Mapping[str, int]) -> int:
    """Return source-mapped SDF formal charge, falling back to the PDB charge field."""
    atom_name = getattr(atom, "atom_name", "")
    if atom_name in formal_charges:
        return int(formal_charges[atom_name])
    return _parse_formal_charge(getattr(atom, "charge", None))


def _guess_element(atom_name: str) -> str:
    """Guess an element symbol from a PDB atom name."""
    letters = "".join(char for char in atom_name.strip() if char.isalpha())
    if not letters:
        return "C"
    first = letters[0].upper()
    if len(letters) >= 2 and letters[:2].title() in {"Cl", "Br", "Na", "Mg", "Ca", "Zn"}:
        return letters[:2].title()
    return first


__all__ = [
    "ProductStatePabloDefinitionSummary",
    "ProductStatePabloLibrary",
    "build_product_state_pablo_library",
    "build_product_state_pablo_library_for_specs",
]
