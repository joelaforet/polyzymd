"""N-glycosylation reaction-template adapter."""

from __future__ import annotations

from pathlib import Path
from typing import Any, ClassVar

from pydantic import BaseModel, Field, field_validator

from polyzymd.builders.conjugation._linkage import (
    ExplicitLinkageContract,
    LinkageBond,
    PdbAtomSelector,
    ReactionProduct,
    ReactiveEndpoint,
    parse_pdb_atom_records,
    resolve_explicit_linkage_contract,
)
from polyzymd.builders.conjugation._pdb_fragment import PdbFragmentLoadResult
from polyzymd.builders.conjugation.polymer import GeneratedMoietyFragment, GeneratedPolymerFragment
from polyzymd.builders.conjugation.reactions.base import (
    PdbFragmentCompatibilityResult,
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


class PdbFragmentLinkageDiagnostic(BaseModel):
    """Residue-level diagnostic for one inter-residue glycan bond."""

    atom_1_serial: int
    atom_1_name: str
    atom_1_residue: str
    atom_2_serial: int
    atom_2_name: str
    atom_2_residue: str
    plausible_glycosidic: bool


class ResidueResolvedGlycanPdbProfileResult(BaseModel):
    """N-glycosylation profile resolved from a residue-resolved glycan PDB fragment."""

    fragment: GeneratedPolymerFragment
    reducing_c1_atom_index: int
    reducing_c1_serial: int
    hydroxyl_oxygen_atom_index: int
    hydroxyl_hydrogen_atom_index: int
    retained_ring_oxygen_atom_index: int
    retained_ring_oxygen_serial: int
    leaving_atom_indices: tuple[int, int]
    leaving_group_representation: str
    linkage_diagnostics: tuple[PdbFragmentLinkageDiagnostic, ...]


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
    supports_coordinate_assembly: ClassVar[bool] = True
    profile_sidecar_key: ClassVar[str] = "n_glycosylation_profile"
    site_participant_name: ClassVar[str] = "asparagine site"
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
                ReactionParticipant(name=cls.site_participant_name, role="site", reactant_index=0),
                ReactionParticipant(name="reducing-end glycan", role="moiety", reactant_index=1),
            ),
            atom_roles=tuple(AtomRoleSpec.model_validate(role) for role in cls.role_metadata()),
            description=cls.description,
        )

    @classmethod
    def resolve_pdb_fragment_source(
        cls,
        pdb_fragment: PdbFragmentLoadResult,
        attachment: Any,
        *,
        settings: NGlycosylationReactionSettings | None = None,
    ) -> PdbFragmentCompatibilityResult:
        """Validate and resolve a residue-resolved glycan PDB fragment."""
        profile = residue_resolved_glycan_pdb_profile_from_fragment(
            pdb_fragment,
            link_site=_field(_field(attachment, "moiety"), "link_site"),
            leaving_atom_names=tuple(
                _field(_field(_field(attachment, "mechanism"), "leaving_atoms"), "moiety", ()) or ()
            ),
        )
        return PdbFragmentCompatibilityResult(
            fragment=profile.fragment,
            reactive_sequence_index=0,
            reactive_selector={
                "chain_id": "C",
                "residue_name": _reactive_residue_name(profile),
                "residue_number": _reactive_residue_number(profile),
                "atom_name": _reactive_atom_name(profile),
                "atom_serial": profile.reducing_c1_serial,
            },
            sidecar_payload={
                cls.profile_sidecar_key: profile.model_dump(mode="json", exclude={"fragment"})
            },
            diagnostics=(f"Resolved residue-resolved glycan PDB fragment for {cls.name}",),
        )

    @classmethod
    def detect_anomeric_group(cls, moiety_fragment: GeneratedMoietyFragment) -> GlycanAnomericGroup:
        """Detect the glycan anomeric carbon and hydroxyl leaving atoms."""
        mol = _rdkit_mol_from_moiety_fragment(moiety_fragment)
        return detect_glycan_anomeric_group(mol)

    @classmethod
    def site_selector_from_config(
        cls,
        site_config: Any,
        *,
        settings: NGlycosylationReactionSettings,
    ) -> PdbAtomSelector:
        """Resolve and validate this mechanism's protein endpoint."""
        return _site_selector_from_config(
            site_config,
            settings=settings,
            mechanism_name=cls.name,
        )

    @classmethod
    def resolve_protein_leaving_atom(
        cls,
        protein_pdb_path: Path | str,
        selector: PdbAtomSelector,
    ) -> PdbAtomRecord:
        """Resolve the deterministic Asn ND2 hydrogen removed by N-glycosylation."""
        return _resolve_bonded_site_hydrogen(
            protein_pdb_path,
            selector,
            mechanism_name=cls.name,
            preferred_name="HD21",
        )

    @classmethod
    def build_contract(
        cls,
        site_config: Any,
        moiety_fragment: GeneratedMoietyFragment | GeneratedPolymerFragment,
        *,
        protein_pdb_path: Path | str | None = None,
        settings: NGlycosylationReactionSettings | None = None,
    ) -> ExplicitLinkageContract:
        """Build a generic linkage contract for an Asn-glycan attachment."""
        resolved = settings or cls.default_settings()
        site_selector = cls.site_selector_from_config(site_config, settings=resolved)
        if isinstance(moiety_fragment, GeneratedMoietyFragment):
            group = cls.detect_anomeric_group(moiety_fragment)
        else:
            group = _anomeric_group_from_pdb_fragment(moiety_fragment)
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
                _selector_from_pdb_atom(
                    cls.resolve_protein_leaving_atom(protein_pdb_path, site_selector)
                ),
            )

        product_site, product_moiety = resolved.product_residue_names(reactive_atom.residue_name)
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
                allow_external_leaving_residue=not isinstance(
                    moiety_fragment, GeneratedMoietyFragment
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
        moiety_fragment: GeneratedMoietyFragment | GeneratedPolymerFragment,
        *,
        prepared_fragment: Any | None = None,
        attachment_id: str = "",
        attachment_index: int = 1,
        attachment_config: Any = None,
        attachment_force_field_domain: str = "",
        diagnostics: tuple[str, ...] = (),
        settings: NGlycosylationReactionSettings | None = None,
    ) -> ReactionProduct:
        """Resolve an Asn-glycan attachment plan against protein and moiety records."""
        from polyzymd.builders.conjugation._specs import resolve_prepared_reaction_fragment

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
            fragment=resolve_prepared_reaction_fragment(moiety_fragment, prepared_fragment),
            attachment_id=attachment_id,
            attachment_index=attachment_index,
            attachment_config=attachment_config,
            reaction_name=cls.name,
            attachment_force_field_domain=attachment_force_field_domain,
            diagnostics=diagnostics,
        )

    def resolve_attachment(
        self,
        protein_pdb_path: Path | str,
        site_config: Any,
        fragment: GeneratedMoietyFragment | GeneratedPolymerFragment,
        *,
        prepared_fragment: Any | None = None,
        settings: NGlycosylationReactionSettings | None = None,
    ) -> ReactionProduct:
        """Resolve an experimental N-glycosylation generic attachment plan."""
        return self.resolve_plan(
            protein_pdb_path,
            site_config,
            fragment,
            prepared_fragment=prepared_fragment,
            settings=settings,
        )


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


def residue_resolved_glycan_pdb_profile_from_fragment(
    load_result: PdbFragmentLoadResult,
    *,
    link_site: Any | None = None,
    leaving_atom_names: tuple[str, ...] = (),
) -> ResidueResolvedGlycanPdbProfileResult:
    """Validate a residue-resolved glycan PDB fragment for N-glycosylation.

    Parameters
    ----------
    load_result : PdbFragmentLoadResult
        Generic PDB-fragment ingestion result.

    Returns
    -------
    ResidueResolvedGlycanPdbProfileResult
        Mechanism-compatible fragment and glycan-specific diagnostics.

    Raises
    ------
    ValueError
        If the fragment lacks a valid reducing-end anomeric C1 hydroxyl group.
    """
    serial_to_atom = {
        atom.serial: atom for atom in load_result.source_atoms if atom.serial is not None
    }
    atoms_by_index = {atom.atom_index: atom for atom in load_result.source_atoms}
    group, c1_atom, hydroxyl_oxygen, hydroxyl_hydrogen, representation = _resolve_pdb_reducing_end(
        load_result.source_atoms,
        load_result.serial_bonds,
        serial_to_atom,
        link_site=link_site,
        leaving_atom_names=leaving_atom_names,
    )
    diagnostics = _linkage_diagnostics(load_result.serial_bonds, serial_to_atom)
    if not any(item.plausible_glycosidic for item in diagnostics):
        raise ValueError(
            "Residue-resolved glycan PDB fragment lacks plausible inter-residue C-O "
            f"glycosidic bonds in {load_result.source_path}"
        )
    fragment = load_result.to_fragment(
        reactive_atom_serial=c1_atom.serial,
        reactive_atom_index=c1_atom.atom_index,
        leaving_atom_serials=(hydroxyl_oxygen.serial, hydroxyl_hydrogen.serial),
        leaving_atom_indices=(hydroxyl_oxygen.atom_index, hydroxyl_hydrogen.atom_index),
    )
    return ResidueResolvedGlycanPdbProfileResult(
        fragment=fragment,
        reducing_c1_atom_index=c1_atom.atom_index,
        reducing_c1_serial=c1_atom.serial,
        hydroxyl_oxygen_atom_index=hydroxyl_oxygen.atom_index,
        hydroxyl_hydrogen_atom_index=hydroxyl_hydrogen.atom_index,
        retained_ring_oxygen_atom_index=group.ring_oxygen_index,
        retained_ring_oxygen_serial=atoms_by_index[group.ring_oxygen_index].serial,
        leaving_atom_indices=(hydroxyl_oxygen.atom_index, hydroxyl_hydrogen.atom_index),
        leaving_group_representation=representation,
        linkage_diagnostics=diagnostics,
    )


def _reactive_residue_name(profile: ResidueResolvedGlycanPdbProfileResult) -> str:
    """Return the residue name for the glycan reactive atom."""
    atom = profile.fragment.atoms[profile.reducing_c1_atom_index]
    return atom.residue_name


def _reactive_residue_number(profile: ResidueResolvedGlycanPdbProfileResult) -> int:
    """Return the residue number for the glycan reactive atom."""
    atom = profile.fragment.atoms[profile.reducing_c1_atom_index]
    return atom.residue_number


def _reactive_atom_name(profile: ResidueResolvedGlycanPdbProfileResult) -> str:
    """Return the atom name for the graph-detected glycan reactive atom."""
    atom = profile.fragment.atoms[profile.reducing_c1_atom_index]
    return atom.atom_name


def _anomeric_group_from_pdb_fragment(fragment: GeneratedPolymerFragment) -> GlycanAnomericGroup:
    """Detect an anomeric group from a generated fragment graph."""
    group = _detect_pdb_anomeric_group(
        _moiety_pdb_atoms(fragment),
        _fragment_serial_bonds(fragment),
        link_site=_fragment_reactive_selector(fragment),
    )
    leaving_indices = set(group.leaving_atom_indices)
    if leaving_indices != set(fragment.leaving_atom_indices):
        raise ValueError(
            "Residue-resolved glycan PDB N-glycosylation leaving selectors do not match "
            "the graph-detected O/H leaving group"
        )
    return group


def _resolve_roh_leaving_atoms(
    atoms: tuple[PdbAtomRecord, ...],
) -> tuple[PdbAtomRecord, PdbAtomRecord]:
    """Resolve the supported separate-residue hydroxyl-cap O1 and HO1 atoms."""
    roh_atoms = [atom for atom in atoms if atom.residue_name.upper() == "ROH"]
    if not roh_atoms:
        raise ValueError(
            "N-glycosylation PDB fragment profile requires one ROH residue with O1 and HO1 "
            "leaving atoms"
        )
    residues = {(atom.chain_id, atom.residue_number, atom.insertion_code) for atom in roh_atoms}
    if len(residues) != 1:
        raise ValueError(
            "Separate-residue hydroxyl-cap representation requires exactly one ROH residue"
        )
    o1 = [atom for atom in roh_atoms if atom.atom_name.upper() == "O1"]
    ho1 = [atom for atom in roh_atoms if atom.atom_name.upper() == "HO1"]
    if len(o1) != 1 or len(ho1) != 1:
        raise ValueError("ROH leaving group must contain unique atoms named O1 and HO1")
    if _normalized_element(o1[0]) != "O" or _normalized_element(ho1[0]) != "H":
        raise ValueError("ROH leaving group requires O1 element O and HO1 element H")
    return o1[0], ho1[0]


def _resolve_pdb_reducing_end(
    atoms: tuple[PdbAtomRecord, ...],
    bonds: tuple[tuple[int, int], ...],
    serial_to_atom: dict[int, PdbAtomRecord],
    *,
    link_site: Any | None = None,
    leaving_atom_names: tuple[str, ...] = (),
) -> tuple[GlycanAnomericGroup, PdbAtomRecord, PdbAtomRecord, PdbAtomRecord, str]:
    """Resolve a structural reducing-end C1 and hydroxyl leaving group."""
    group = _detect_pdb_anomeric_group(
        atoms,
        bonds,
        link_site=link_site,
        leaving_atom_names=leaving_atom_names,
    )
    atoms_by_index = {atom.atom_index: atom for atom in atoms}
    c1_atom = atoms_by_index[group.reactive_carbon_index]
    hydroxyl_oxygen = atoms_by_index[group.hydroxyl_oxygen_index]
    hydroxyl_hydrogens = [
        atoms_by_index[index]
        for index in group.leaving_atom_indices
        if index != group.hydroxyl_oxygen_index
    ]
    if len(hydroxyl_hydrogens) != 1:
        raise ValueError("Graph-detected glycan leaving group must contain exactly one hydrogen")
    representation = (
        "separate_residue_hydroxyl_cap"
        if _residue_key(c1_atom) != _residue_key(hydroxyl_oxygen)
        else "local_oh"
    )
    return group, c1_atom, hydroxyl_oxygen, hydroxyl_hydrogens[0], representation


def _detect_pdb_anomeric_group(
    atoms: tuple[PdbAtomRecord, ...],
    bonds: tuple[tuple[int, int], ...],
    *,
    link_site: Any | None = None,
    leaving_atom_names: tuple[str, ...] = (),
) -> GlycanAnomericGroup:
    """Detect an anomeric C from explicit graph chemistry, not atom names."""
    serial_to_atom = {atom.serial: atom for atom in atoms if atom.serial is not None}
    neighbors = _serial_neighbors(bonds)
    link_matches = _matching_link_site_atoms(atoms, link_site)
    candidate_atoms = link_matches or atoms
    candidates: list[GlycanAnomericGroup] = []

    for carbon in candidate_atoms:
        if carbon.serial is None or carbon.atom_index is None or _normalized_element(carbon) != "C":
            continue
        oxygen_neighbors = [
            serial_to_atom[serial]
            for serial in neighbors.get(carbon.serial, set())
            if serial in serial_to_atom and _normalized_element(serial_to_atom[serial]) == "O"
        ]
        hydroxyl_pairs = [
            (oxygen, hydrogen)
            for oxygen in oxygen_neighbors
            for hydrogen in _bonded_hydrogens(oxygen, neighbors, serial_to_atom)
            if oxygen.serial is not None
            and hydrogen.serial is not None
            and neighbors.get(oxygen.serial, set()) == {carbon.serial, hydrogen.serial}
            if _leaving_names_match((oxygen, hydrogen), leaving_atom_names)
        ]
        retained_oxygens = [
            oxygen
            for oxygen in oxygen_neighbors
            if not _bonded_hydrogens(oxygen, neighbors, serial_to_atom)
            and _oxygen_connects_into_sugar_graph(oxygen, carbon, neighbors, serial_to_atom)
        ]
        for hydroxyl_oxygen, hydroxyl_hydrogen in hydroxyl_pairs:
            for retained_oxygen in retained_oxygens:
                if retained_oxygen.serial == hydroxyl_oxygen.serial:
                    continue
                candidates.append(
                    GlycanAnomericGroup(
                        reactive_carbon_index=carbon.atom_index,
                        hydroxyl_oxygen_index=hydroxyl_oxygen.atom_index,
                        ring_oxygen_index=retained_oxygen.atom_index,
                        leaving_atom_indices=(
                            hydroxyl_oxygen.atom_index,
                            hydroxyl_hydrogen.atom_index,
                        ),
                        candidate_atom_indices=tuple(
                            sorted(
                                atom.atom_index
                                for atom in candidate_atoms
                                if atom.atom_index is not None and _normalized_element(atom) == "C"
                            )
                        ),
                        evidence={
                            "source": "pdb_conect_graph",
                            "reactive_carbon_serial": carbon.serial,
                            "hydroxyl_oxygen_serial": hydroxyl_oxygen.serial,
                            "hydroxyl_hydrogen_serial": hydroxyl_hydrogen.serial,
                            "retained_ring_oxygen_serial": retained_oxygen.serial,
                        },
                    )
                )

    unique = {
        (
            candidate.reactive_carbon_index,
            candidate.hydroxyl_oxygen_index,
            candidate.ring_oxygen_index,
        ): candidate
        for candidate in candidates
    }
    if not unique:
        message = (
            "No glycan anomeric motif was found. Expected one graph-detected carbon bonded "
            "to a hydroxyl O-H leaving group and a retained ring oxygen."
        )
        if leaving_atom_names:
            message += f" Configured leaving atoms: {', '.join(leaving_atom_names)}."
        raise ValueError(message)
    if len(unique) > 1 and not link_matches:
        assignments = ", ".join(
            f"C{carbon}/O{hydroxyl}/O{ring}" for carbon, hydroxyl, ring in sorted(unique)
        )
        raise ValueError(
            "Ambiguous glycan anomeric motif assignments; configure the existing moiety.link_site "
            f"reactive atom selector to choose one: {assignments}"
        )
    if len(unique) > 1:
        assignments = ", ".join(
            f"C{carbon}/O{hydroxyl}/O{ring}" for carbon, hydroxyl, ring in sorted(unique)
        )
        raise ValueError(f"Configured moiety.link_site still leaves ambiguous motif: {assignments}")
    return next(iter(unique.values()))


def _oxygen_connects_into_sugar_graph(
    oxygen: PdbAtomRecord,
    carbon: PdbAtomRecord,
    neighbors: dict[int, set[int]],
    serial_to_atom: dict[int, PdbAtomRecord],
) -> bool:
    """Return whether an oxygen is retained in a nonterminal ring-like graph."""
    if oxygen.serial is None or carbon.serial is None:
        return False
    retained_neighbors = neighbors.get(oxygen.serial, set()) - {carbon.serial}
    heavy_neighbors = {
        serial
        for serial in retained_neighbors
        if (atom := serial_to_atom.get(serial)) is not None and _normalized_element(atom) != "H"
    }
    if not heavy_neighbors:
        return False
    return _has_alternate_heavy_path_to_carbon(
        starts=heavy_neighbors,
        target=carbon.serial,
        excluded_edge=frozenset((oxygen.serial, carbon.serial)),
        neighbors=neighbors,
        serial_to_atom=serial_to_atom,
    )


def _has_alternate_heavy_path_to_carbon(
    *,
    starts: set[int],
    target: int,
    excluded_edge: frozenset[int],
    neighbors: dict[int, set[int]],
    serial_to_atom: dict[int, PdbAtomRecord],
) -> bool:
    """Return whether a heavy-atom path reaches target without one excluded edge."""
    visited: set[int] = set()
    queue = list(starts)
    while queue:
        current = queue.pop(0)
        if current in visited:
            continue
        visited.add(current)
        if current == target:
            return True
        for neighbor in neighbors.get(current, set()):
            if frozenset((current, neighbor)) == excluded_edge:
                continue
            atom = serial_to_atom.get(neighbor)
            if atom is None or _normalized_element(atom) == "H" or neighbor in visited:
                continue
            queue.append(neighbor)
    return False


def _resolve_local_reducing_hydroxyl(
    atoms: tuple[PdbAtomRecord, ...],
    bonds: tuple[tuple[int, int], ...],
    *,
    link_site: Any | None = None,
    leaving_atom_names: tuple[str, ...] = (),
) -> tuple[PdbAtomRecord, PdbAtomRecord, PdbAtomRecord]:
    """Resolve an ordinary residue-local anomeric C1 hydroxyl from graph connectivity."""
    serial_to_atom = {atom.serial: atom for atom in atoms if atom.serial is not None}
    neighbors = _serial_neighbors(bonds)
    link_matches = _matching_link_site_atoms(atoms, link_site)
    c1_candidates = link_matches or [
        atom
        for atom in atoms
        if atom.atom_name.upper() == "C1" and _normalized_element(atom) == "C"
    ]
    candidates: list[tuple[PdbAtomRecord, PdbAtomRecord, PdbAtomRecord]] = []
    for carbon in c1_candidates:
        if _normalized_element(carbon) != "C":
            continue
        oxygen_neighbors = [
            serial_to_atom[serial]
            for serial in neighbors.get(carbon.serial, set())
            if serial in serial_to_atom and _normalized_element(serial_to_atom[serial]) == "O"
        ]
        residue_local_oxygens = [
            atom for atom in oxygen_neighbors if _residue_key(atom) == _residue_key(carbon)
        ]
        hydroxyl_pairs = [
            (oxygen, hydrogen)
            for oxygen in residue_local_oxygens
            for hydrogen in _bonded_hydrogens(oxygen, neighbors, serial_to_atom)
            if _leaving_names_match((oxygen, hydrogen), leaving_atom_names)
        ]
        ring_like_oxygens = [
            atom
            for atom in residue_local_oxygens
            if not _bonded_hydrogens(atom, neighbors, serial_to_atom)
            and _oxygen_connects_into_sugar_graph(atom, carbon, neighbors, serial_to_atom)
        ]
        for oxygen, hydrogen in hydroxyl_pairs:
            if any(ring_oxygen.serial != oxygen.serial for ring_oxygen in ring_like_oxygens):
                candidates.append((carbon, oxygen, hydrogen))

    if len(candidates) != 1:
        raise ValueError(
            "Residue-resolved glycan PDB fragment requires a unique reducing-end C1 bonded "
            f"to an explicit hydroxyl O/H group; found {len(candidates)} candidates"
        )
    return candidates[0]


def _has_explicit_roh_cap(atoms: tuple[PdbAtomRecord, ...]) -> bool:
    """Return whether the supported separate-residue hydroxyl cap is present."""
    return any(atom.residue_name.upper() == "ROH" for atom in atoms)


def _resolve_reducing_c1(
    bonds: tuple[tuple[int, int], ...],
    serial_to_atom: dict[int, PdbAtomRecord],
    roh_o1: PdbAtomRecord,
    *,
    link_site: Any | None = None,
) -> PdbAtomRecord:
    """Resolve the reducing sugar C1 bonded to ROH O1."""
    link_matches = _matching_link_site_atoms(tuple(serial_to_atom.values()), link_site)
    candidates = []
    for left, right in bonds:
        if roh_o1.serial not in {left, right}:
            continue
        other = right if left == roh_o1.serial else left
        atom = serial_to_atom[other]
        if (
            atom.residue_name.upper() != "ROH"
            and atom.atom_name.upper() == "C1"
            and _normalized_element(atom) == "C"
            and (not link_matches or atom.serial == link_matches[0].serial)
        ):
            candidates.append(atom)
    if len(candidates) != 1:
        raise ValueError(
            "Residue-resolved glycan PDB fragment requires a unique reducing sugar C1 "
            f"bonded to hydroxyl-cap O1; found {len(candidates)} candidates"
        )
    return candidates[0]


def _serial_neighbors(bonds: tuple[tuple[int, int], ...]) -> dict[int, set[int]]:
    """Return an undirected adjacency map keyed by PDB atom serial."""
    neighbors: dict[int, set[int]] = {}
    for left, right in bonds:
        neighbors.setdefault(left, set()).add(right)
        neighbors.setdefault(right, set()).add(left)
    return neighbors


def _bonded_hydrogens(
    atom: PdbAtomRecord,
    neighbors: dict[int, set[int]],
    serial_to_atom: dict[int, PdbAtomRecord],
) -> tuple[PdbAtomRecord, ...]:
    """Return explicit hydrogens directly bonded to an atom."""
    hydrogens = [
        serial_to_atom[serial]
        for serial in neighbors.get(atom.serial, set())
        if serial in serial_to_atom and _is_hydrogen_atom(serial_to_atom[serial])
    ]
    return tuple(sorted(hydrogens, key=lambda item: (item.serial, item.atom_name)))


def _leaving_names_match(
    atoms: tuple[PdbAtomRecord, ...], leaving_atom_names: tuple[str, ...]
) -> bool:
    """Return whether a candidate hydroxyl matches configured leaving atom names."""
    if not leaving_atom_names:
        return True
    expected = sorted(name.strip().upper() for name in leaving_atom_names)
    observed = sorted(atom.atom_name.strip().upper() for atom in atoms)
    return observed == expected


def _matching_link_site_atoms(
    atoms: tuple[PdbAtomRecord, ...], link_site: Any | None
) -> tuple[PdbAtomRecord, ...]:
    """Return atoms matching an optional configured moiety link-site selector."""
    if link_site is None:
        return ()
    selector = PdbAtomSelector(
        chain_id=_coalesce_text(_field(link_site, "chain_id"), "C"),
        residue_name=_coalesce_text(_field(link_site, "residue_name")),
        residue_number=int(_required_selector_field(link_site, "residue_number")),
        atom_name=_coalesce_text(_field(link_site, "atom_name")),
        insertion_code=str(_field(link_site, "insertion_code", "") or ""),
        atom_serial=_field(link_site, "atom_serial"),
        atom_index=_field(link_site, "atom_index"),
    )
    matches = tuple(atom for atom in atoms if selector.matches(atom))
    if len(matches) != 1:
        raise ValueError(
            "Configured moiety.link_site must resolve exactly one glycan atom; "
            f"found {len(matches)}"
        )
    return matches


def _validate_roh_reducing_subgraph(
    bonds: tuple[tuple[int, int], ...],
    *,
    roh_o1: PdbAtomRecord,
    roh_ho1: PdbAtomRecord,
    c1_atom: PdbAtomRecord,
) -> None:
    """Validate the exact ROH leaving-group subgraph accepted for cleavage."""
    if c1_atom.residue_name.upper() == "ROH":
        raise ValueError("Reducing sugar C1 must belong to a non-ROH residue")
    if _normalized_element(c1_atom) != "C":
        raise ValueError("Reducing sugar C1 must have element C")
    bond_set = {frozenset(bond) for bond in bonds}
    required = {
        frozenset((roh_ho1.serial, roh_o1.serial)),
        frozenset((roh_o1.serial, c1_atom.serial)),
    }
    if not required.issubset(bond_set):
        raise ValueError(
            "Separate-residue hydroxyl-cap representation requires exact ROH:HO1-ROH:O1-C1 "
            "reducing subgraph"
        )

    neighbors: dict[int, set[int]] = {
        roh_ho1.serial: set(),
        roh_o1.serial: set(),
        c1_atom.serial: set(),
    }
    for left, right in bonds:
        if left in neighbors:
            neighbors[left].add(right)
        if right in neighbors:
            neighbors[right].add(left)
    if neighbors[roh_ho1.serial] != {roh_o1.serial}:
        raise ValueError("Hydroxyl-cap ROH:HO1 must be bonded only to ROH:O1")
    if neighbors[roh_o1.serial] != {roh_ho1.serial, c1_atom.serial}:
        raise ValueError("Hydroxyl-cap ROH:O1 must be bonded exactly to ROH:HO1 and reducing C1")


def _linkage_diagnostics(
    bonds: tuple[tuple[int, int], ...], serial_to_atom: dict[int, PdbAtomRecord]
) -> tuple[PdbFragmentLinkageDiagnostic, ...]:
    """Build inter-residue glycan linkage diagnostics from graph bonds."""
    diagnostics = []
    for left, right in bonds:
        atom_1 = serial_to_atom[left]
        atom_2 = serial_to_atom[right]
        if _residue_key(atom_1) == _residue_key(atom_2):
            continue
        pair = {_normalized_element(atom_1), _normalized_element(atom_2)}
        diagnostics.append(
            PdbFragmentLinkageDiagnostic(
                atom_1_serial=left,
                atom_1_name=atom_1.atom_name,
                atom_1_residue=_residue_label(atom_1),
                atom_2_serial=right,
                atom_2_name=atom_2.atom_name,
                atom_2_residue=_residue_label(atom_2),
                plausible_glycosidic=pair == {"C", "O"},
            )
        )
    return tuple(diagnostics)


def _residue_key(atom: PdbAtomRecord) -> tuple[str, int, str, str]:
    """Return a residue identity key."""
    return (atom.chain_id, atom.residue_number, atom.insertion_code, atom.residue_name)


def _residue_label(atom: PdbAtomRecord) -> str:
    """Return a compact residue label for diagnostics."""
    return f"{atom.chain_id}:{atom.residue_name}{atom.residue_number}{atom.insertion_code}"


def _normalized_element(atom: PdbAtomRecord) -> str:
    """Return an uppercase element symbol for validation."""
    return atom.element.strip().upper()


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
    mechanism_name: str,
) -> PdbAtomSelector:
    """Build and validate a glycosylation protein endpoint selector."""
    residue_name = _coalesce_text(
        _field(site_config, "residue_name"), settings.source_site_residue_name
    )
    atom_name = _coalesce_text(_field(site_config, "atom_name"), settings.target_atom_name)
    if residue_name != settings.source_site_residue_name:
        raise ValueError(
            f"{mechanism_name} target residue must be "
            f"{settings.source_site_residue_name}, got {residue_name}"
        )
    if atom_name != settings.target_atom_name:
        raise ValueError(
            f"{mechanism_name} target atom must be {settings.target_atom_name}, got {atom_name}"
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


def _resolve_bonded_site_hydrogen(
    protein_pdb_path: Path | str,
    selector: PdbAtomSelector,
    *,
    mechanism_name: str,
    preferred_name: str | None = None,
    residue_library: Any | None = None,
) -> PdbAtomRecord:
    """Resolve one leaving hydrogen from Pablo residue-template connectivity."""
    from polyzymd.builders.conjugation.pablo.residue_library import bonded_hydrogen_names

    atoms = parse_pdb_atom_records(protein_pdb_path)
    site_matches = [atom for atom in atoms if selector.matches(atom)]
    if len(site_matches) != 1:
        raise ValueError(
            f"Expected exactly one {selector.residue_name} {selector.atom_name} atom for "
            f"{mechanism_name} site "
            f"{selector.chain_id}:{selector.residue_number}:{selector.atom_name}, "
            f"found {len(site_matches)}"
        )
    site_atom = site_matches[0]
    template_hydrogen_names = bonded_hydrogen_names(
        selector.residue_name,
        selector.atom_name,
        residue_library=residue_library,
    )
    residue_atoms = [atom for atom in atoms if _same_residue(atom, site_atom)]
    observed = [
        atom
        for atom in residue_atoms
        if atom.atom_name.strip().upper() in template_hydrogen_names and _is_hydrogen_atom(atom)
    ]
    if len(observed) == 1:
        return observed[0]

    preferred = (preferred_name or "").strip().upper()
    if preferred and len(observed) > 1:
        preferred_matches = [
            atom for atom in observed if atom.atom_name.strip().upper() == preferred
        ]
        if len(preferred_matches) == 1:
            return preferred_matches[0]

    if len(observed) != 1:
        raise ValueError(
            f"{mechanism_name} requires exactly one explicit hydrogen bonded to "
            f"{selector.residue_name} {selector.atom_name}"
            + (f", or an unambiguous preferred {preferred} hydrogen; " if preferred else "; ")
            + f"found {len(observed)} for {selector.chain_id}:{selector.residue_name}"
            f"{selector.residue_number}:{selector.atom_name}. "
            f"Template bonded H names: {', '.join(template_hydrogen_names)}. No coordinate "
            "geometry fallback is used."
        )


def _same_residue(left: PdbAtomRecord, right: PdbAtomRecord) -> bool:
    """Return whether two atoms share a PDB residue identity."""
    return (
        left.chain_id.upper() == right.chain_id.upper()
        and left.residue_name.upper() == right.residue_name.upper()
        and left.residue_number == right.residue_number
        and left.insertion_code.upper() == right.insertion_code.upper()
    )


def _moiety_pdb_atoms(
    fragment: GeneratedMoietyFragment | GeneratedPolymerFragment,
) -> tuple[PdbAtomRecord, ...]:
    """Convert a generated moiety fragment to PDB atom records."""
    atoms: list[PdbAtomRecord] = []
    for atom in fragment.atoms:
        pdb_atom = atom.to_pdb_atom()
        if not pdb_atom.chain_id:
            pdb_atom = pdb_atom.model_copy(update={"chain_id": "C"})
        atoms.append(pdb_atom)
    return tuple(atoms)


def _fragment_serial_bonds(fragment: GeneratedPolymerFragment) -> tuple[tuple[int, int], ...]:
    """Return generated-fragment bonds normalized to atom serial pairs."""
    atoms_by_index = {atom.atom_index: atom for atom in fragment.atoms}
    atoms_by_serial = {atom.serial: atom for atom in fragment.atoms if atom.serial is not None}
    atoms_by_name = {atom.atom_name.upper(): atom for atom in fragment.atoms}
    serial_bonds = []
    for left, right in fragment.bonds:
        left_atom = _fragment_bond_atom(left, atoms_by_index, atoms_by_serial, atoms_by_name)
        right_atom = _fragment_bond_atom(right, atoms_by_index, atoms_by_serial, atoms_by_name)
        if left_atom.serial is None or right_atom.serial is None:
            raise ValueError("Generated glycan fragment graph requires serial-numbered atoms")
        serial_bonds.append(tuple(sorted((left_atom.serial, right_atom.serial))))
    if not serial_bonds:
        raise ValueError("Generated glycan fragment graph requires explicit bonds")
    return tuple(sorted(set(serial_bonds)))


def _fragment_bond_atom(
    value: int | str,
    atoms_by_index: dict[int, Any],
    atoms_by_serial: dict[int, Any],
    atoms_by_name: dict[str, Any],
) -> Any:
    """Resolve a generated-fragment bond endpoint to an atom record."""
    if isinstance(value, int):
        if value in atoms_by_serial:
            return atoms_by_serial[value]
        if value in atoms_by_index:
            return atoms_by_index[value]
    else:
        atom = atoms_by_name.get(value.strip().upper())
        if atom is not None:
            return atom
    raise ValueError(f"Generated glycan fragment bond references unknown atom {value!r}")


def _fragment_reactive_selector(fragment: GeneratedPolymerFragment) -> PdbAtomSelector | None:
    """Return the existing reactive selector as a graph-disambiguation hint."""
    matches = [
        atom
        for atom in fragment.atoms
        if (
            fragment.reactive_atom_serial is not None
            and atom.serial == fragment.reactive_atom_serial
        )
        or (
            fragment.reactive_atom_index is not None
            and atom.atom_index == fragment.reactive_atom_index
        )
        or (
            fragment.reactive_atom_name is not None
            and atom.atom_name.upper() == fragment.reactive_atom_name.upper()
        )
    ]
    if len(matches) != 1:
        return None
    atom = matches[0].to_pdb_atom()
    return _selector_from_pdb_atom(atom)


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


def _required_selector_field(obj: Any, name: str) -> Any:
    """Return a required field from an atom selector mapping or object."""
    value = _field(obj, name)
    if value is None:
        raise ValueError(f"moiety.link_site.{name} is required for N-glycosylation")
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
    "PdbFragmentLinkageDiagnostic",
    "ResidueResolvedGlycanPdbProfileResult",
    "NGlycosylationReaction",
    "NGlycosylationReactionSettings",
    "detect_glycan_anomeric_group",
    "residue_resolved_glycan_pdb_profile_from_fragment",
]
