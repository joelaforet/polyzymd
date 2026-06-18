"""Base interfaces for conjugation reaction templates."""

from __future__ import annotations

from abc import ABC
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
    """Base class for conjugation reaction templates.

    Phase 1 keeps reaction templates as lightweight descriptors. Concrete
    chemistry remains in the legacy implementation until migration.
    """

    name: ClassVar[str]
    description: ClassVar[str] = ""
    aliases: ClassVar[tuple[str, ...]] = ()

    @classmethod
    def identifiers(cls) -> tuple[str, ...]:
        """Return registry identifiers accepted for this template."""
        return (cls.name, *cls.aliases)

    def plan(self, _context: ReactionContext) -> ReactionResult:
        """Plan the reaction once concrete templates are migrated."""
        raise NotImplementedError(
            f"Reaction template {self.name!r} has not been migrated to the new "
            "conjugation reaction API yet."
        )
