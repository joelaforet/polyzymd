"""O-glycosylation reaction template for serine and threonine sites."""

from __future__ import annotations

from pathlib import Path
from typing import Any, ClassVar

from polyzymd.builders.conjugation._linkage import PdbAtomSelector
from polyzymd.builders.conjugation.reactions.n_glycosylation import (
    NGlycosylationReaction,
    NGlycosylationReactionSettings,
    _coalesce_optional_text,
    _coalesce_text,
    _field,
    _field_if_set,
    _resolve_bonded_site_hydrogen,
)
from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord


class OGlycosylationReactionSettings(NGlycosylationReactionSettings):
    """Settings for O-glycosylation at one resolved Ser or Thr hydroxyl."""

    source_site_residue_name: str = "SER"
    product_site_residue_name: str = "OLS"
    target_atom_name: str = "OG"
    target_bond_length_angstrom: float = 1.43


class OGlycosylationReaction(NGlycosylationReaction):
    """O-glycosylation of Ser OG or Thr OG1 by a reducing-end glycan."""

    name: ClassVar[str] = "o_glycosylation"
    aliases: ClassVar[tuple[str, ...]] = ("o_glycan", "ser_thr_o_glycosylation")
    description: ClassVar[str] = "O-glycosylation of serine or threonine by a glycan C1."
    Settings: ClassVar[type[OGlycosylationReactionSettings]] = OGlycosylationReactionSettings
    coordinate_backend_mechanism: ClassVar[str] = "o_glycosylation"
    profile_sidecar_key: ClassVar[str] = "o_glycosylation_profile"
    site_participant_name: ClassVar[str] = "serine or threonine site"
    mapped_reaction_smarts: ClassVar[str] = "[O:1][H:2].[C:3]([O:4][H:5])([O:6])>>[O:1][C:3]([O:6])"
    atom_role_metadata: ClassVar[tuple[dict[str, Any], ...]] = (
        {"map_number": 1, "participant": "site", "role": "linking", "label": "site_oxygen"},
        {"map_number": 2, "participant": "site", "role": "leaving", "label": "site_hydrogen"},
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
    def settings_from_attachment(cls, attachment: Any) -> OGlycosylationReactionSettings:
        """Resolve GLYCAM-compatible O-linked site labels from attachment config."""
        site = _field(attachment, "site")
        residue_name = _coalesce_text(_field(site, "residue_name"), "")
        defaults = cls._site_defaults(residue_name)
        mechanism = _field(attachment, "mechanism")
        product_residues = _field(mechanism, "product_residues")
        bond = _field(mechanism, "bond")
        return cls.Settings(
            source_site_residue_name=residue_name,
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
    def build_contract(cls, site_config: Any, moiety_fragment: Any, **kwargs: Any) -> Any:
        """Build an O-glycosylation contract with site-specific defaults."""
        if kwargs.get("settings") is None:
            kwargs["settings"] = cls._site_defaults(_field(site_config, "residue_name"))
        return super().build_contract(site_config, moiety_fragment, **kwargs)

    @classmethod
    def resolve_protein_leaving_atom(
        cls,
        protein_pdb_path: Path | str,
        selector: PdbAtomSelector,
    ) -> PdbAtomRecord:
        """Resolve the hydroxyl hydrogen removed from Ser OG or Thr OG1."""
        return _resolve_bonded_site_hydrogen(
            protein_pdb_path,
            selector,
            mechanism_name=cls.name,
        )

    @classmethod
    def _site_defaults(cls, residue_name: Any) -> OGlycosylationReactionSettings:
        normalized = str(residue_name or "").strip().upper()
        if normalized == "SER":
            return cls.Settings(
                source_site_residue_name="SER",
                product_site_residue_name="OLS",
                target_atom_name="OG",
            )
        if normalized == "THR":
            return cls.Settings(
                source_site_residue_name="THR",
                product_site_residue_name="OLT",
                target_atom_name="OG1",
            )
        raise ValueError(
            "o_glycosylation requires site.residue_name SER or THR; "
            f"received {normalized or '<missing>'}"
        )


__all__ = ["OGlycosylationReaction", "OGlycosylationReactionSettings"]
