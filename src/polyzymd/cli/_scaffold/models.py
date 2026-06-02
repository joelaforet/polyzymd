"""Models used by the analysis scaffold renderer."""

from __future__ import annotations

from dataclasses import dataclass

# Default single-file MDAnalysis-native scaffold style
DEFAULT_STYLE = "simple"
VALID_STYLES = (DEFAULT_STYLE, "dict")
ADVANCED_STYLES = ("dict",)


@dataclass(frozen=True)
class ScaffoldSpec:
    """Rendering specification for a generated analysis plugin.

    Parameters
    ----------
    name : str
        Plugin name in snake_case.
    class_name : str
        PascalCase class prefix used for generated classes.
    style : str
        Scaffold style. ``"simple"`` creates the default single-file
        MDAnalysis-native contributor path. ``"dict"`` creates an advanced
        package scaffold using canonical artifact payloads.
    """

    name: str
    class_name: str
    style: str

    @property
    def title(self) -> str:
        """Return the human-readable title for the plugin name.

        Returns
        -------
        str
            Plugin title derived from ``name``.
        """
        return self.name.replace("_", " ").title()

    @property
    def uses_single_file_layout(self) -> bool:
        """Return whether the scaffold should use a single plugin file.

        Returns
        -------
        bool
            True when ``style`` is the default simple scaffold style.
        """
        return self.style == DEFAULT_STYLE

    @property
    def uses_package_layout(self) -> bool:
        """Return whether the scaffold should generate a plugin package.

        Returns
        -------
        bool
            True when ``style`` requests an advanced package scaffold.
        """
        return self.style in ADVANCED_STYLES
