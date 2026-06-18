"""NHS-Lys reaction-template placeholder."""

from __future__ import annotations

from importlib import import_module
from types import ModuleType
from typing import ClassVar

from polyzymd.builders.conjugation.reactions.base import (
    ReactionContext,
    ReactionResult,
    ReactionTemplate,
)


class NhsLysReaction(ReactionTemplate):
    """Descriptor for the legacy NHS ester to lysine amide implementation.

    The concrete implementation currently lives in
    ``polyzymd.builders.conjugation.nhs_lys`` and related workflow helpers.
    Future phases should migrate that behavior behind this template.
    """

    name: ClassVar[str] = "nhs_lys"
    aliases: ClassVar[tuple[str, ...]] = ("nhs_lys_amide",)
    description: ClassVar[str] = "NHS ester coupling to lysine NZ to form an amide."
    legacy_module_path: ClassVar[str] = "polyzymd.builders.conjugation.nhs_lys"
    legacy_symbols: ClassVar[tuple[str, ...]] = (
        "LysineReactiveSite",
        "NhsReactiveGroup",
        "NhsLysGraphEditPlan",
        "detect_nhs_reactive_group",
        "execute_nhs_lys_amide_rdkit_graph_edit",
        "extract_lysine_reactive_site",
        "plan_nhs_lys_amide",
    )

    @classmethod
    def load_legacy_module(cls) -> ModuleType:
        """Import and return the current NHS-Lys implementation module."""
        return import_module(cls.legacy_module_path)

    def plan(self, _context: ReactionContext) -> ReactionResult:
        """Raise until the legacy NHS-Lys planner is migrated."""
        raise NotImplementedError(
            "NHS-Lys reaction planning has not been migrated to the new reaction "
            "template API yet. Use the legacy functions in "
            "polyzymd.builders.conjugation.nhs_lys for Phase 1."
        )
