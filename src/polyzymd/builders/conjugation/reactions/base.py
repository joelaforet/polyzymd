"""Base interfaces for conjugation reaction templates."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, ClassVar

from pydantic import BaseModel, Field


class ReactionContext(BaseModel):
    """Placeholder context passed to future reaction templates."""

    model_config = {"arbitrary_types_allowed": True}

    protein: Any | None = None
    moiety: Any | None = None
    metadata: dict[str, Any] = Field(default_factory=dict)


class ReactionResult(BaseModel):
    """Placeholder result returned by future reaction templates."""

    model_config = {"arbitrary_types_allowed": True}

    plan: Any | None = None
    metadata: dict[str, Any] = Field(default_factory=dict)


class ReactionTemplate(ABC):
    """Base class for conjugation reaction mechanisms.

    Concrete mechanisms resolve a protein attachment request and generated
    moiety into a generic atom-level attachment plan. Construction code consumes
    that plan without mechanism-specific product-state branches.
    """

    name: ClassVar[str]
    description: ClassVar[str] = ""
    aliases: ClassVar[tuple[str, ...]] = ()

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

    def plan(self, context: ReactionContext) -> ReactionResult:
        """Resolve a generic plan from a reaction context."""
        site_config = context.metadata.get("site_config") or context.metadata.get("site")
        settings = context.metadata.get("settings")
        if context.protein is None or context.moiety is None or site_config is None:
            raise ValueError(
                "Reaction planning requires context.protein, context.moiety, and "
                "metadata['site_config']"
            )
        plan = self.resolve_attachment(
            context.protein,
            site_config,
            context.moiety,
            settings=settings,
        )
        return ReactionResult(plan=plan, metadata={"mechanism": self.name})
