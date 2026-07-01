"""NHS-Lys reaction-template adapter."""

from __future__ import annotations

from pathlib import Path
from typing import Any, ClassVar

from pydantic import BaseModel, Field, field_validator

from polyzymd.builders.conjugation.reactions._rdkit_graph import (
    AddedBond,
    LysineReactiveSite,
    NhsLysAttachmentSite,
    NhsLysGraphEditPlan,
    NhsReactiveGroup,
    RdkitGraphEditResult,
    detect_nhs_reactive_group,
    execute_nhs_lys_amide_rdkit_graph_edit,
    extract_lysine_reactive_site,
    plan_nhs_lys_amide,
)
from polyzymd.builders.conjugation.reactions.base import ReactionTemplate


class NhsLysReactionSettings(BaseModel):
    """Settings owned by the built-in NHS-Lys reaction template."""

    source_site_residue_name: str = Field(
        "LYS",
        min_length=1,
        max_length=4,
        description="Input protein residue expected at the attachment site",
    )
    product_site_residue_name: str = Field(
        "LYX",
        min_length=1,
        max_length=4,
        description="Product residue name for the linked lysine site",
    )
    product_moiety_residue_name: str = Field(
        "NHX",
        min_length=1,
        max_length=4,
        description="Product residue name for the linked NHS-derived moiety residue",
    )
    target_atom_name: str = Field(
        "NZ",
        min_length=1,
        max_length=4,
        description="Protein atom that forms the new amide bond",
    )
    bond_order: int | float = Field(1, description="Order of the product amide linkage bond")
    target_bond_length_angstrom: float = Field(
        1.33,
        gt=0,
        description="Default target length for the lysine NZ to acyl-carbon bond",
    )
    max_nz_hydrogens_to_remove: int = Field(
        2,
        ge=1,
        description="Maximum number of lysine NZ hydrogens removed by the product policy",
    )
    canonicalize_site_hydrogens: bool = Field(
        True,
        description="Hint that source lysine hydrogens should be canonicalized before linking",
    )
    minimize_product_junction: bool = Field(
        True,
        description="Hint that the NHS-Lys junction should be locally minimized after assembly",
    )

    @field_validator(
        "source_site_residue_name",
        "product_site_residue_name",
        "product_moiety_residue_name",
        "target_atom_name",
    )
    @classmethod
    def normalize_pdb_label(cls, value: str) -> str:
        """Normalize PDB residue and atom labels used by the adapter."""
        return value.strip().upper()

    @property
    def product_residue_names(self) -> tuple[str, str]:
        """Return protein-site and moiety product residue names in Pablo order."""
        return (self.product_site_residue_name, self.product_moiety_residue_name)


class NhsLysReaction(ReactionTemplate):
    """Reaction template for NHS ester coupling to lysine NZ."""

    name: ClassVar[str] = "nhs_lys"
    aliases: ClassVar[tuple[str, ...]] = ("nhs_lys_amide",)
    description: ClassVar[str] = "NHS ester coupling to lysine NZ to form an amide."
    Settings: ClassVar[type[NhsLysReactionSettings]] = NhsLysReactionSettings
    coordinate_backend_mechanism: ClassVar[str] = "nhs_lys_amide"
    supports_coordinate_assembly: ClassVar[bool] = True
    mapped_reaction_smarts: ClassVar[str] = (
        "[N:1]([H:2]).[C:3](=[O:4])[O:5][N:6]>>[N:1][C:3](=[O:4])"
    )
    atom_role_metadata: ClassVar[tuple[dict[str, Any], ...]] = (
        {"map_number": 1, "participant": "site", "role": "linking", "label": "lysine_nz"},
        {
            "map_number": 2,
            "participant": "site",
            "role": "leaving",
            "label": "site_hydrogen",
        },
        {
            "map_number": 3,
            "participant": "moiety",
            "role": "linking",
            "label": "nhs_acyl_carbon",
        },
        {
            "map_number": 4,
            "participant": "moiety",
            "role": "geometry_anchor",
            "label": "carbonyl_oxygen",
        },
        {
            "map_number": 5,
            "participant": "moiety",
            "role": "leaving",
            "label": "nhs_bridging_oxygen",
        },
    )

    @classmethod
    def default_settings(cls) -> NhsLysReactionSettings:
        """Return a fresh settings model with template-owned defaults."""
        return cls.Settings()

    @classmethod
    def settings_from_attachment(cls, attachment: Any) -> NhsLysReactionSettings:
        """Resolve NHS-Lys settings from a config attachment plus template defaults."""
        defaults = cls.default_settings()
        site = getattr(attachment, "site", None)
        mechanism = getattr(attachment, "mechanism", None)
        product_residues = getattr(mechanism, "product_residues", None)
        bond = getattr(mechanism, "bond", None)
        return cls.Settings(
            source_site_residue_name=_coalesce_text(
                getattr(site, "residue_name", None), defaults.source_site_residue_name
            ),
            product_site_residue_name=_coalesce_text(
                getattr(product_residues, "site", None), defaults.product_site_residue_name
            ),
            product_moiety_residue_name=_coalesce_text(
                getattr(product_residues, "moiety", None), defaults.product_moiety_residue_name
            ),
            target_atom_name=_coalesce_text(
                getattr(site, "atom_name", None),
                getattr(bond, "site_atom", None),
                defaults.target_atom_name,
            ),
            bond_order=getattr(bond, "order", defaults.bond_order),
            target_bond_length_angstrom=getattr(
                bond,
                "target_bond_length_angstrom",
                defaults.target_bond_length_angstrom,
            ),
            max_nz_hydrogens_to_remove=defaults.max_nz_hydrogens_to_remove,
            canonicalize_site_hydrogens=defaults.canonicalize_site_hydrogens,
            minimize_product_junction=defaults.minimize_product_junction,
        )

    @classmethod
    def product_residue_names(
        cls,
        settings: NhsLysReactionSettings | None = None,
    ) -> tuple[str, str]:
        """Return the product residue names for this reaction template."""
        resolved = settings or cls.default_settings()
        return resolved.product_residue_names

    @classmethod
    def build_contract(
        cls,
        site_config: Any,
        fragment: Any,
        *,
        protein_pdb_path: Path | str | None = None,
        settings: NhsLysReactionSettings | None = None,
    ) -> Any:
        """Build the generic NHS-Lys linkage contract from site and fragment data.

        Parameters
        ----------
        site_config : Any
            Attachment-site configuration for the source lysine.
        fragment : Any
            Generated fragment containing an NHS ester reactive group.
        protein_pdb_path : pathlib.Path, str, or None, optional
            Optional protein PDB used to resolve the leaving lysine hydrogen,
            by default ``None``.
        settings : NhsLysReactionSettings or None, optional
            Resolved reaction settings, by default ``None``.

        Returns
        -------
        Any
            Generic explicit linkage contract equivalent to the NHS-Lys reaction.
        """
        resolved = settings or cls.default_settings()
        linker = _build_nhs_lys_linker(
            target_chain=_coalesce_text(getattr(site_config, "chain_id", None), "A"),
            target_residue_name=_coalesce_text(
                getattr(site_config, "residue_name", None), resolved.source_site_residue_name
            ),
            target_residue_number=_required_int(
                getattr(site_config, "residue_number", None), "site.residue_number"
            ),
            target_insertion_code=str(getattr(site_config, "insertion_code", "") or ""),
            settings=resolved,
        )
        contract = linker.generic_contract(fragment, protein_pdb_path=protein_pdb_path)
        return contract.model_copy(
            update={
                "bond": contract.bond.model_copy(
                    update={
                        "bond_order": resolved.bond_order,
                        "target_bond_length_angstrom": resolved.target_bond_length_angstrom,
                    }
                )
            }
        )

    @classmethod
    def resolve_plan(
        cls,
        protein_pdb_path: Path | str,
        site_config: Any,
        fragment: Any,
        *,
        settings: NhsLysReactionSettings | None = None,
    ) -> Any:
        """Resolve an NHS-Lys attachment plan through the generic contract path."""
        from polyzymd.builders.conjugation._linkage import resolve_explicit_linkage_contract

        contract = cls.build_contract(
            site_config,
            fragment,
            protein_pdb_path=protein_pdb_path,
            settings=settings,
        )
        return resolve_explicit_linkage_contract(protein_pdb_path, fragment, contract)

    def resolve_attachment(
        self,
        protein_pdb_path: Path | str,
        site_config: Any,
        fragment: Any,
        *,
        settings: NhsLysReactionSettings | None = None,
    ) -> Any:
        """Resolve an NHS-Lys attachment plan for the generic construction flow."""
        return self.resolve_plan(
            protein_pdb_path,
            site_config,
            fragment,
            settings=settings,
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
        """Build the generic atom-mapped role model for NHS-Lys diagnostics."""
        from polyzymd.builders.conjugation.reactions._roles import (
            AtomMappedReaction,
            AtomRoleSpec,
            ReactionParticipant,
        )

        return AtomMappedReaction.from_reaction_smarts(
            name=cls.coordinate_backend_mechanism,
            reaction_smarts=cls.reaction_smarts(),
            participants=(
                ReactionParticipant(name="lysine site", role="site", reactant_index=0),
                ReactionParticipant(name="NHS ester moiety", role="moiety", reactant_index=1),
            ),
            atom_roles=tuple(AtomRoleSpec.model_validate(role) for role in cls.role_metadata()),
            description=cls.description,
        )

    @classmethod
    def plan_graph_edit(
        cls,
        site: Any,
        reactive_group: Any,
        *,
        site_hydrogen_indices_to_remove: Any | None = None,
    ) -> Any:
        """Plan an explicit NHS-Lys graph edit."""
        return plan_nhs_lys_amide(
            site,
            reactive_group,
            site_hydrogen_indices_to_remove=site_hydrogen_indices_to_remove,
        )


def _coalesce_text(*values: str | None) -> str:
    """Return the first non-empty text value."""
    for value in values:
        if value is None:
            continue
        normalized = str(value).strip()
        if normalized:
            return normalized
    raise ValueError("Expected at least one non-empty text value")


def _build_nhs_lys_linker(
    *,
    target_residue_number: int,
    target_chain: str,
    target_residue_name: str,
    target_insertion_code: str,
    settings: NhsLysReactionSettings,
) -> Any:
    """Create the private NHS-Lys linker implementation from template settings."""
    from polyzymd.builders.conjugation._linkage import NhsLysModifierLinker

    return NhsLysModifierLinker(
        target_chain=target_chain,
        target_residue_number=target_residue_number,
        target_residue_name=target_residue_name,
        target_insertion_code=target_insertion_code,
        target_atom_name=settings.target_atom_name,
        lysine_target_resname=settings.product_site_residue_name,
        modifier_target_resname=settings.product_moiety_residue_name,
        max_nz_hydrogens_to_remove=settings.max_nz_hydrogens_to_remove,
    )


def _required_int(value: Any, label: str) -> int:
    """Return an integer config value or raise a field-specific error."""
    if value is None:
        raise ValueError(f"{label} is required for NHS-Lys linker creation")
    return int(value)


__all__ = [
    "AddedBond",
    "LysineReactiveSite",
    "NhsLysAttachmentSite",
    "NhsLysGraphEditPlan",
    "NhsLysReaction",
    "NhsLysReactionSettings",
    "NhsReactiveGroup",
    "RdkitGraphEditResult",
    "detect_nhs_reactive_group",
    "execute_nhs_lys_amide_rdkit_graph_edit",
    "extract_lysine_reactive_site",
    "plan_nhs_lys_amide",
]
