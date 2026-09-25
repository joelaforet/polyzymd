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
    module : str or None, optional
        Module the generated test imports the plugin from. ``None`` (default)
        means the built-in ``polyzymd.analyses.<name>``; a study analysis is
        imported by its bare ``<name>``.
    """

    name: str
    class_name: str
    module: str | None = None

    @property
    def import_path(self) -> str:
        """Module the generated test imports the plugin from."""
        return self.module or f"polyzymd.analyses.{self.name}"

    @property
    def title(self) -> str:
        """Human-readable title derived from the plugin name."""
        return self.name.replace("_", " ").title()
