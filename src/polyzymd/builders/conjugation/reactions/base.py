"""Base interfaces for conjugation reaction templates."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, ClassVar


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

    @abstractmethod
    def resolve_attachment(
        self,
        protein_pdb_path: Any,
        site_config: Any,
        fragment: Any,
        *,
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
