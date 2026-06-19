"""Built-in conjugation reaction-template registry."""

from __future__ import annotations

from polyzymd.builders.conjugation.reactions.base import ReactionTemplate
from polyzymd.builders.conjugation.reactions.n_glycosylation import NGlycosylationReaction
from polyzymd.builders.conjugation.reactions.nhs_lys import NhsLysReaction

_REACTIONS: dict[str, type[ReactionTemplate]] = {
    **dict.fromkeys(NhsLysReaction.identifiers(), NhsLysReaction),
    **dict.fromkeys(NGlycosylationReaction.identifiers(), NGlycosylationReaction),
}


def list_reactions() -> dict[str, type[ReactionTemplate]]:
    """Return the available built-in reaction templates by identifier."""
    return dict(_REACTIONS)


def get_reaction(name: str) -> type[ReactionTemplate]:
    """Return a built-in reaction template class by name or alias."""
    try:
        return _REACTIONS[name]
    except KeyError as exc:
        available = ", ".join(sorted(_REACTIONS))
        raise KeyError(f"Unknown conjugation reaction {name!r}. Available: {available}") from exc
