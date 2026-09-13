"""Specification the analysis scaffold renderer works from."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ScaffoldSpec:
    """What ``polyzymd new-analysis`` renders.

    Parameters
    ----------
    name : str
        Plugin name in snake_case.
    class_name : str
        PascalCase class prefix used for the generated settings model and
        plugin class.
    """

    name: str
    class_name: str

    @property
    def title(self) -> str:
        """Human-readable title derived from the plugin name."""
        return self.name.replace("_", " ").title()
