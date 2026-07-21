"""Base interfaces for conjugation reaction templates."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, ClassVar

from pydantic import BaseModel, Field


class PdbFragmentCompatibilityResult(BaseModel):
    """Mechanism-owned result for a compatible PDB-fragment moiety source."""

    model_config = {"arbitrary_types_allowed": True}

    fragment: Any = Field(..., exclude=True)
    reactive_sequence_index: int | None = None
    reactive_selector: dict[str, int | str] | None = None
    sidecar_payload: dict[str, Any] = Field(default_factory=dict)
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)


class ReactionTemplate(ABC):
    """Base class for conjugation reaction mechanisms.

    Concrete mechanisms resolve a protein attachment request and generated
    moiety into a generic atom-level attachment plan. Construction code consumes
    that plan without mechanism-specific product-state branches.
    """

    name: ClassVar[str]
    description: ClassVar[str] = ""
    aliases: ClassVar[tuple[str, ...]] = ()
    supports_coordinate_assembly: ClassVar[bool] = False

    @classmethod
    def identifiers(cls) -> tuple[str, ...]:
        """Return registry identifiers accepted for this template."""
        return (cls.name, *cls.aliases)

    @classmethod
    def resolve_pdb_fragment_source(
        cls,
        pdb_fragment: Any,
        attachment: Any,
        *,
        settings: Any | None = None,
    ) -> PdbFragmentCompatibilityResult:
        """Resolve a compatible residue-resolved PDB-fragment moiety.

        Parameters
        ----------
        pdb_fragment : Any
            Generic loaded PDB-fragment result from the provider.
        attachment : Any
            Attachment configuration selecting this reaction template.
        settings : Any or None, optional
            Mechanism-specific settings, by default ``None``.

        Returns
        -------
        PdbFragmentCompatibilityResult
            Template-owned construction fragment and metadata.

        Raises
        ------
        ValueError
            If this reaction template does not opt in to PDB-fragment input.
        """
        mechanism_name = getattr(getattr(attachment, "mechanism", None), "name", cls.name)
        raise ValueError(
            "attachment.moiety.input_path loaded a generic PDB fragment, but mechanism "
            f"{mechanism_name!r} does not support PDB-fragment moiety sources"
        )

    @abstractmethod
    def resolve_attachment(
        self,
        protein_pdb_path: Any,
        site_config: Any,
        fragment: Any,
        *,
        prepared_fragment: Any | None = None,
        settings: Any | None = None,
    ) -> Any:
        """Resolve a generic attachment plan for construction.

        Parameters
        ----------
        protein_pdb_path : Any
            Protein PDB source used to resolve site atoms.
        site_config : Any
            Attachment-site configuration or compatible object.
        fragment : Any
            Generated moiety or polymer fragment from a provider.
        settings : Any or None, optional
            Mechanism-specific settings, by default ``None``.

        Returns
        -------
        Any
            A mechanism-owned generic resolved attachment plan.
        """
