"""N-glycosylation reaction-template adapter."""

from __future__ import annotations

from pathlib import Path
from typing import Any, ClassVar

from pydantic import BaseModel, Field, field_validator

from polyzymd.builders.conjugation._linkage import (
    ExplicitLinkageContract,
    LinkageBond,
    PdbAtomSelector,
    ReactiveEndpoint,
    ResolvedAttachmentPlan,
    parse_pdb_atom_records,
    resolve_explicit_linkage_contract,
)
from polyzymd.builders.conjugation.polymer import GeneratedMoietyFragment
from polyzymd.builders.conjugation.reactions.base import (
    ReactionContext,
    ReactionResult,
    ReactionTemplate,
)
from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord


class GlycanAnomericGroup(BaseModel):
    """Resolved reducing-end anomeric carbon and leaving hydroxyl atoms."""

    reactive_carbon_index: int = Field(..., ge=0, description="Anomeric carbon atom index")
    hydroxyl_oxygen_index: int = Field(..., ge=0, description="Anomeric hydroxyl oxygen atom index")
    ring_oxygen_index: int = Field(..., ge=0, description="Ring oxygen bonded to C1")
    leaving_atom_indices: tuple[int, ...] = Field(
        ..., min_length=1, description="Anomeric hydroxyl atoms removed during linkage"
    )
    candidate_atom_indices: tuple[int, ...] = Field(
        default_factory=tuple, description="SMARTS/connectivity candidates considered"
    )
    evidence: dict[str, Any] = Field(default_factory=dict)


class NGlycosylationReactionSettings(BaseModel):
    """Settings owned by the built-in N-glycosylation reaction template."""

    source_site_residue_name: str = Field(
        "ASN",
        min_length=1,
        max_length=4,
        description="Input protein residue expected at the attachment site",
    )
    product_site_residue_name: str = Field(
        "ASX",
        min_length=1,
        max_length=4,
        description="Product residue name for the linked asparagine site",
    )
    product_moiety_residue_name: str | None = Field(
        None,
        min_length=1,
        max_length=4,
        description="Optional product residue name for the glycan moiety; defaults to fragment",
    )
    target_atom_name: str = Field(
        "ND2",
        min_length=1,
        max_length=4,
        description="Asparagine atom that forms the N-glycosidic bond",
    )
    bond_order: int | float = Field(1, description="Order of the product C-N linkage bond")
    target_bond_length_angstrom: float = Field(
        1.45,
        gt=0,
        description="Default target length for the Asn ND2 to anomeric-carbon bond",
    )

    @field_validator(
        "source_site_residue_name",
        "product_site_residue_name",
        "product_moiety_residue_name",
        "target_atom_name",
    )
    @classmethod
    def normalize_pdb_label(cls, value: str | None) -> str | None:
        """Normalize PDB residue and atom labels used by the adapter."""
        return value.strip().upper() if value is not None else None

    def product_residue_names(self, moiety_residue_name: str) -> tuple[str, str]:
        """Return protein-site and glycan product residue names in Pablo order."""
        return (
            self.product_site_residue_name,
            self.product_moiety_residue_name or moiety_residue_name.strip().upper(),
        )


class NGlycosylationReaction(ReactionTemplate):
    """Reaction template for N-glycosylation of Asn ND2 by a reducing-end glycan."""

    name: ClassVar[str] = "n_glycosylation"
    aliases: ClassVar[tuple[str, ...]] = ("n_glycan_asn", "asn_n_glycosylation")
    description: ClassVar[str] = "N-glycosylation of asparagine ND2 by a glycan C1."
    Settings: ClassVar[type[NGlycosylationReactionSettings]] = NGlycosylationReactionSettings
    coordinate_backend_mechanism: ClassVar[str] = "n_glycosylation"
    mapped_reaction_smarts: ClassVar[str] = (
        "[N:1]([H:2]).[C:3]([O:4][H:5])([O:6])>>[N:1][C:3]([O:6])"
    )
    atom_role_metadata: ClassVar[tuple[dict[str, Any], ...]] = (
        {"map_number": 1, "participant": "site", "role": "linking", "label": "asn_nd2"},
        {"map_number": 2, "participant": "site", "role": "leaving", "label": "nd2_hydrogen"},
        {
            "map_number": 3,
            "participant": "moiety",
            "role": "linking",
            "label": "anomeric_carbon",
        },
        {
            "map_number": 4,
            "participant": "moiety",
            "role": "leaving",
            "label": "anomeric_hydroxyl_oxygen",
        },
        {
            "map_number": 5,
            "participant": "moiety",
            "role": "leaving",
            "label": "anomeric_hydroxyl_hydrogen",
        },
        {
            "map_number": 6,
            "participant": "moiety",
            "role": "geometry_anchor",
            "label": "ring_oxygen",
        },
    )

    @classmethod
    def default_settings(cls) -> NGlycosylationReactionSettings:
        """Return a fresh settings model with template-owned defaults."""
        return cls.Settings()

    @classmethod
    def settings_from_attachment(cls, attachment: Any) -> NGlycosylationReactionSettings:
        """Resolve N-glycosylation settings from a config attachment plus defaults."""
        defaults = cls.default_settings()
        site = getattr(attachment, "site", None)
        mechanism = getattr(attachment, "mechanism", None)
        product_residues = getattr(mechanism, "product_residues", None)
        bond = getattr(mechanism, "bond", None)
        return cls.Settings(
            source_site_residue_name=_coalesce_text(
                _field(site, "residue_name"), defaults.source_site_residue_name
            ),
            product_site_residue_name=_coalesce_text(
                _field(product_residues, "site"), defaults.product_site_residue_name
            ),
            product_moiety_residue_name=_coalesce_optional_text(
                _field(product_residues, "moiety"), defaults.product_moiety_residue_name
            ),
            target_atom_name=_coalesce_text(
                _field(site, "atom_name"),
                _field(bond, "site_atom"),
                defaults.target_atom_name,
            ),
            bond_order=_field_if_set(bond, "order", defaults.bond_order),
            target_bond_length_angstrom=_field_if_set(
                bond,
                "target_bond_length_angstrom",
                defaults.target_bond_length_angstrom,
            ),
        )

    @classmethod
    def reaction_smarts(cls) -> str:
        """Return the atom-mapped reaction SMARTS used for diagnostics."""
        return cls.mapped_reaction_smarts

    @classmethod
    def role_metadata(cls) -> tuple[dict[str, Any], ...]:
        """Return atom-map role metadata for diagnostics and generic preflights."""
        return tuple(dict(role) for role in cls.atom_role_metadata)

    @classmethod
    def build_role_model(cls) -> Any:
        """Build the generic atom-mapped role model for N-glycosylation diagnostics."""
        from polyzymd.builders.conjugation.reactions._roles import (
            AtomMappedReaction,
            AtomRoleSpec,
            ReactionParticipant,
        )

        return AtomMappedReaction.from_reaction_smarts(
            name=cls.coordinate_backend_mechanism,
            reaction_smarts=cls.reaction_smarts(),
            participants=(
                ReactionParticipant(name="asparagine site", role="site", reactant_index=0),
                ReactionParticipant(name="reducing-end glycan", role="moiety", reactant_index=1),
            ),
            atom_roles=tuple(AtomRoleSpec.model_validate(role) for role in cls.role_metadata()),
            description=cls.description,
        )

    @classmethod
    def detect_anomeric_group(cls, moiety_fragment: GeneratedMoietyFragment) -> GlycanAnomericGroup:
        """Detect the glycan anomeric carbon and hydroxyl leaving atoms."""
        mol = _rdkit_mol_from_moiety_fragment(moiety_fragment)
        return detect_glycan_anomeric_group(mol)

    @classmethod
    def build_contract(
        cls,
        site_config: Any,
        moiety_fragment: GeneratedMoietyFragment,
        *,
        protein_pdb_path: Path | str | None = None,
        settings: NGlycosylationReactionSettings | None = None,
    ) -> ExplicitLinkageContract:
        """Build a generic linkage contract for an Asn-glycan attachment."""
        resolved = settings or cls.default_settings()
        site_selector = _site_selector_from_config(site_config, settings=resolved)
        group = cls.detect_anomeric_group(moiety_fragment)
        moiety_atoms = _moiety_pdb_atoms(moiety_fragment)
        atoms_by_index = {atom.atom_index: atom for atom in moiety_atoms}
        reactive_atom = atoms_by_index.get(group.reactive_carbon_index)
        if reactive_atom is None:
            raise ValueError(
                "Detected glycan anomeric carbon could not be mapped to the moiety fragment"
            )
        leaving_atoms = _atoms_for_indices(moiety_atoms, group.leaving_atom_indices)
        protein_leaving_selectors: tuple[PdbAtomSelector, ...] = ()
        if protein_pdb_path is not None:
            protein_leaving_selectors = (
                _selector_from_pdb_atom(_resolve_asn_nd2_hydrogen(protein_pdb_path, site_selector)),
            )

        product_site, product_moiety = resolved.product_residue_names(moiety_fragment.residue_name)
        modifier_selector = _selector_from_pdb_atom(reactive_atom)
        return ExplicitLinkageContract(
            protein_endpoint=ReactiveEndpoint(
                participant="protein",
                selector=site_selector,
                product_residue_name=product_site,
                leaving_atom_selectors=protein_leaving_selectors,
            ),
            modifier_endpoint=ReactiveEndpoint(
                participant="modifier",
                selector=modifier_selector,
                product_residue_name=product_moiety,
                leaving_atom_selectors=tuple(
                    _selector_from_pdb_atom(atom) for atom in leaving_atoms
                ),
            ),
            bond=LinkageBond(
                protein_atom_name=resolved.target_atom_name,
                modifier_atom_name=reactive_atom.atom_name,
                bond_order=resolved.bond_order,
                target_bond_length_angstrom=resolved.target_bond_length_angstrom,
            ),
            mechanism_name=cls.name,
        )

    @classmethod
    def resolve_plan(
        cls,
        protein_pdb_path: Path | str,
        site_config: Any,
        moiety_fragment: GeneratedMoietyFragment,
        *,
        settings: NGlycosylationReactionSettings | None = None,
    ) -> ResolvedAttachmentPlan:
        """Resolve an Asn-glycan attachment plan against protein and moiety records."""
        contract = cls.build_contract(
            site_config,
            moiety_fragment,
            protein_pdb_path=protein_pdb_path,
            settings=settings,
        )
        return resolve_explicit_linkage_contract(
            protein_pdb_path,
            _moiety_pdb_atoms(moiety_fragment),
            contract,
        )

    def plan(self, context: ReactionContext) -> ReactionResult:
        """Resolve a plan from a reaction context carrying protein, site, and moiety."""
        site_config = context.metadata.get("site_config") or context.metadata.get("site")
        settings = context.metadata.get("settings")
        if context.protein is None or context.moiety is None or site_config is None:
            raise ValueError(
                "N-glycosylation planning requires context.protein, context.moiety, and "
                "metadata['site_config']"
            )
        plan = self.resolve_plan(
            context.protein,
            site_config,
            context.moiety,
            settings=settings,
        )
        return ReactionResult(plan=plan, metadata={"mechanism": self.name})


def detect_glycan_anomeric_group(mol: Any) -> GlycanAnomericGroup:
    """Detect a reducing-end anomeric carbon with a removable hydroxyl group."""
    from rdkit import Chem

    pattern = Chem.MolFromSmarts("[#6;R](-[#8;R])-[#8]")
    if pattern is None:
        raise ValueError("Internal anomeric-carbon SMARTS could not be parsed")

    smarts_candidates = {match[0] for match in mol.GetSubstructMatches(pattern)}
    candidates: list[GlycanAnomericGroup] = []
    for carbon_index in sorted(smarts_candidates):
        atom = mol.GetAtomWithIdx(carbon_index)
        if atom.GetAtomicNum() != 6 or not atom.IsInRing():
            continue
        oxygen_neighbors = [
            neighbor
            for neighbor in atom.GetNeighbors()
            if neighbor.GetAtomicNum() == 8
            and mol.GetBondBetweenAtoms(carbon_index, neighbor.GetIdx()).GetBondTypeAsDouble()
            == 1.0
        ]
        ring_oxygens = [oxygen for oxygen in oxygen_neighbors if oxygen.IsInRing()]
        hydroxyl_oxygens = [
            oxygen for oxygen in oxygen_neighbors if _oxygen_bound_hydrogen_indices(oxygen)
        ]
        for ring_oxygen in ring_oxygens:
            for hydroxyl_oxygen in hydroxyl_oxygens:
                if ring_oxygen.GetIdx() == hydroxyl_oxygen.GetIdx():
                    continue
                hydroxyl_hydrogens = _oxygen_bound_hydrogen_indices(hydroxyl_oxygen)
                leaving_indices = (hydroxyl_oxygen.GetIdx(), *hydroxyl_hydrogens)
                candidates.append(
                    GlycanAnomericGroup(
                        reactive_carbon_index=carbon_index,
                        hydroxyl_oxygen_index=hydroxyl_oxygen.GetIdx(),
                        ring_oxygen_index=ring_oxygen.GetIdx(),
                        leaving_atom_indices=tuple(sorted(leaving_indices)),
                        candidate_atom_indices=tuple(sorted(smarts_candidates)),
                        evidence={
                            "smarts": "[#6;R](-[#8;R])-[#8]",
                            "oxygen_neighbor_indices": tuple(
                                sorted(oxygen.GetIdx() for oxygen in oxygen_neighbors)
                            ),
                            "hydroxyl_hydrogen_indices": hydroxyl_hydrogens,
                        },
                    )
                )

    if not candidates:
        raise ValueError(
            "No glycan anomeric motif was found. Expected one ring carbon bonded to a ring "
            "oxygen and an exocyclic hydroxyl oxygen."
        )
    unique = {
        (
            candidate.reactive_carbon_index,
            candidate.hydroxyl_oxygen_index,
            candidate.ring_oxygen_index,
        ): candidate
        for candidate in candidates
    }
    if len(unique) > 1:
        assignments = ", ".join(
            f"C{carbon}/O{hydroxyl}/O{ring}" for carbon, hydroxyl, ring in sorted(unique)
        )
        raise ValueError(f"Ambiguous glycan anomeric motif assignments: {assignments}")
    return next(iter(unique.values()))


def _rdkit_mol_from_moiety_fragment(fragment: GeneratedMoietyFragment) -> Any:
    """Return an RDKit molecule for a generated moiety fragment."""
    from rdkit import Chem

    existing = getattr(fragment, "rdkit_mol", None) or getattr(fragment, "mol", None)
    if existing is not None:
        return existing
    if fragment.sdf_path is not None and Path(fragment.sdf_path).exists():
        supplier = Chem.SDMolSupplier(str(fragment.sdf_path), removeHs=False, sanitize=True)
        mols = [mol for mol in supplier if mol is not None]
        if len(mols) == 1:
            return mols[0]

    rw_mol = Chem.RWMol()
    atoms_by_index = {atom.atom_index: atom for atom in fragment.atoms}
    for atom_index in sorted(atoms_by_index):
        fragment_atom = atoms_by_index[atom_index]
        rd_atom = Chem.Atom((fragment_atom.element or fragment_atom.atom_name[0]).strip().title())
        rd_atom.SetFormalCharge(fragment_atom.formal_charge or 0)
        rw_mol.AddAtom(rd_atom)

    rdkit_index_by_fragment_index = {
        fragment_index: rdkit_index
        for rdkit_index, fragment_index in enumerate(sorted(atoms_by_index))
    }
    serial_to_fragment_index = {
        atom.serial: atom.atom_index for atom in fragment.atoms if atom.serial is not None
    }
    for left_serial, right_serial, order in fragment.bond_orders:
        left_index = serial_to_fragment_index.get(int(left_serial), int(left_serial) - 1)
        right_index = serial_to_fragment_index.get(int(right_serial), int(right_serial) - 1)
        rw_mol.AddBond(
            rdkit_index_by_fragment_index[left_index],
            rdkit_index_by_fragment_index[right_index],
            _rdkit_bond_type(order),
        )
    mol = rw_mol.GetMol()
    Chem.SanitizeMol(mol)
    return mol


def _rdkit_bond_type(order: float) -> Any:
    """Map a numeric bond order to an RDKit bond type."""
    from rdkit import Chem

    if order == 1.0:
        return Chem.BondType.SINGLE
    if order == 2.0:
        return Chem.BondType.DOUBLE
    if order == 3.0:
        return Chem.BondType.TRIPLE
    if order == 1.5:
        return Chem.BondType.AROMATIC
    raise ValueError(f"Unsupported moiety bond order for N-glycosylation detection: {order}")


def _site_selector_from_config(
    site_config: Any,
    *,
    settings: NGlycosylationReactionSettings,
) -> PdbAtomSelector:
    """Build the protein endpoint selector, defaulting the omitted atom to ND2."""
    residue_name = _coalesce_text(
        _field(site_config, "residue_name"), settings.source_site_residue_name
    )
    atom_name = _coalesce_text(_field(site_config, "atom_name"), settings.target_atom_name)
    if residue_name != settings.source_site_residue_name:
        raise ValueError(
            "N-glycosylation target residue must be "
            f"{settings.source_site_residue_name}, got {residue_name}"
        )
    if atom_name != settings.target_atom_name:
        raise ValueError(
            "N-glycosylation target atom must be " f"{settings.target_atom_name}, got {atom_name}"
        )
    return PdbAtomSelector(
        chain_id=_coalesce_text(_field(site_config, "chain_id"), "A"),
        residue_name=residue_name,
        residue_number=int(_required_field(site_config, "residue_number")),
        atom_name=atom_name,
        insertion_code=str(_field(site_config, "insertion_code", "") or ""),
        atom_serial=_field(site_config, "atom_serial"),
        atom_index=_field(site_config, "atom_index"),
    )


def _resolve_asn_nd2_hydrogen(
    protein_pdb_path: Path | str,
    selector: PdbAtomSelector,
) -> PdbAtomRecord:
    """Resolve one leaving hydrogen bound to the selected Asn ND2 atom."""
    atoms = parse_pdb_atom_records(protein_pdb_path)
    nd2_matches = [atom for atom in atoms if selector.matches(atom)]
    if len(nd2_matches) != 1:
        raise ValueError(
            "Expected exactly one ASN ND2 atom for N-glycosylation site "
            f"{selector.chain_id}:{selector.residue_number}:{selector.atom_name}, "
            f"found {len(nd2_matches)}"
        )
    nd2_atom = nd2_matches[0]
    hydrogens = _geometry_bound_hydrogens(atoms, nd2_atom)
    if not hydrogens:
        raise ValueError(
            "N-glycosylation requires one explicit Asn ND2 hydrogen to remove; "
            f"found none for {selector.chain_id}:ASN{selector.residue_number}:ND2"
        )
    return hydrogens[0]


def _geometry_bound_hydrogens(
    atoms: tuple[PdbAtomRecord, ...], target_atom: PdbAtomRecord, *, cutoff: float = 1.35
) -> tuple[PdbAtomRecord, ...]:
    """Return same-residue hydrogens within a local N-H distance cutoff."""
    from math import dist

    target_xyz = (target_atom.x, target_atom.y, target_atom.z)
    hydrogens = [
        atom
        for atom in atoms
        if atom.chain_id.upper() == target_atom.chain_id.upper()
        and atom.residue_name.upper() == target_atom.residue_name.upper()
        and atom.residue_number == target_atom.residue_number
        and atom.insertion_code.upper() == target_atom.insertion_code.upper()
        and _is_hydrogen_atom(atom)
        and dist(target_xyz, (atom.x, atom.y, atom.z)) <= cutoff
    ]
    return tuple(
        sorted(
            hydrogens,
            key=lambda atom: (
                atom.atom_index if atom.atom_index is not None else 10**9,
                atom.serial if atom.serial is not None else 10**9,
                atom.atom_name.upper(),
            ),
        )
    )


def _moiety_pdb_atoms(fragment: GeneratedMoietyFragment) -> tuple[PdbAtomRecord, ...]:
    """Convert a generated moiety fragment to PDB atom records."""
    atoms: list[PdbAtomRecord] = []
    for atom in fragment.atoms:
        pdb_atom = atom.to_pdb_atom()
        if not pdb_atom.chain_id:
            pdb_atom = pdb_atom.model_copy(update={"chain_id": "C"})
        atoms.append(pdb_atom)
    return tuple(atoms)


def _atoms_for_indices(
    atoms: tuple[PdbAtomRecord, ...], atom_indices: tuple[int, ...]
) -> tuple[PdbAtomRecord, ...]:
    """Return PDB atoms for resolved fragment atom indices."""
    atoms_by_index = {atom.atom_index: atom for atom in atoms}
    missing = [index for index in atom_indices if index not in atoms_by_index]
    if missing:
        raise ValueError(f"Detected glycan leaving atoms could not be mapped: {missing}")
    return tuple(atoms_by_index[index] for index in atom_indices)


def _selector_from_pdb_atom(atom: PdbAtomRecord) -> PdbAtomSelector:
    """Build an exact selector for a resolved PDB atom."""
    return PdbAtomSelector(
        chain_id=atom.chain_id or "C",
        residue_name=atom.residue_name,
        residue_number=atom.residue_number,
        insertion_code=atom.insertion_code,
        atom_name=atom.atom_name,
        atom_serial=atom.serial,
        atom_index=atom.atom_index,
    )


def _oxygen_bound_hydrogen_indices(oxygen: Any) -> tuple[int, ...]:
    """Return explicit hydrogens bonded to an oxygen atom."""
    return tuple(
        sorted(
            neighbor.GetIdx() for neighbor in oxygen.GetNeighbors() if neighbor.GetAtomicNum() == 1
        )
    )


def _is_hydrogen_atom(atom: PdbAtomRecord) -> bool:
    """Return whether an atom record represents hydrogen."""
    return (atom.element or "").upper() == "H" or atom.atom_name.upper().startswith("H")


def _field(obj: Any, name: str, default: Any = None) -> Any:
    """Return a field from a mapping or object."""
    if obj is None:
        return default
    if isinstance(obj, dict):
        return obj.get(name, default)
    return getattr(obj, name, default)


def _field_if_set(obj: Any, name: str, default: Any = None) -> Any:
    """Return a Pydantic field only when explicitly supplied, otherwise default."""
    if obj is None:
        return default
    if isinstance(obj, dict):
        return obj.get(name, default)
    fields_set = getattr(obj, "model_fields_set", None)
    if fields_set is not None and name not in fields_set:
        return default
    return getattr(obj, name, default)


def _required_field(obj: Any, name: str) -> Any:
    """Return a required field from a mapping or object."""
    value = _field(obj, name)
    if value is None:
        raise ValueError(f"site_config.{name} is required for N-glycosylation")
    return value


def _coalesce_text(*values: str | None) -> str:
    """Return the first non-empty normalized text value."""
    for value in values:
        if value is None:
            continue
        normalized = str(value).strip().upper()
        if normalized:
            return normalized
    raise ValueError("Expected at least one non-empty text value")


def _coalesce_optional_text(*values: str | None) -> str | None:
    """Return the first non-empty normalized text value, or None."""
    for value in values:
        if value is None:
            continue
        normalized = str(value).strip().upper()
        if normalized:
            return normalized
    return None


__all__ = [
    "GlycanAnomericGroup",
    "NGlycosylationReaction",
    "NGlycosylationReactionSettings",
    "detect_glycan_anomeric_group",
]
