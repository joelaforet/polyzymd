"""Default native dynamic polymer reaction markers."""

DEFAULT_REACTION_MARKER = "default"


def get_atrp_reaction_paths() -> dict[str, str]:
    """Get default native dynamic reaction markers.

    Returns:
        Dictionary with keys ``initiation``, ``polymerization``, and
        ``termination`` mapped to the literal ``"default"`` marker.
    """
    return {
        "initiation": DEFAULT_REACTION_MARKER,
        "polymerization": DEFAULT_REACTION_MARKER,
        "termination": DEFAULT_REACTION_MARKER,
    }


def get_atrp_initiation_path() -> str:
    """Get the default initiation marker."""
    return DEFAULT_REACTION_MARKER


def get_atrp_polymerization_path() -> str:
    """Get the default polymerization marker."""
    return DEFAULT_REACTION_MARKER


def get_atrp_termination_path() -> str:
    """Get the default termination marker."""
    return DEFAULT_REACTION_MARKER


def is_default_atrp_reaction_set(
    initiation: str,
    polymerization: str,
    termination: str,
) -> bool:
    """Return whether selectors are the native default reaction set.

    Parameters
    ----------
    initiation : str
        Candidate initiation selector.
    polymerization : str
        Candidate polymerization selector.
    termination : str
        Candidate termination selector.

    Returns
    -------
    bool
        ``True`` when all three selectors are the literal default marker.
    """
    return all(
        str(candidate).strip().lower() == DEFAULT_REACTION_MARKER
        for candidate in (initiation, polymerization, termination)
    )
