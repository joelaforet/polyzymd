"""Models used by the analysis scaffold renderer."""

from __future__ import annotations

from dataclasses import dataclass

# Default single-file MDAnalysis-native scaffold style
DEFAULT_STYLE = "simple"
# Single-file observable-contract scaffold: a Settings model plus compute()
CONTRACT_STYLE = "contract"
VALID_STYLES = (DEFAULT_STYLE, "dict", CONTRACT_STYLE)
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
        package scaffold using canonical artifact payloads. ``"contract"``
        creates the single-file observable-contract plugin.
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
            True for the simple and contract scaffold styles.
        """
        return self.style in (DEFAULT_STYLE, CONTRACT_STYLE)

    @property
    def uses_package_layout(self) -> bool:
        """Return whether the scaffold should generate a plugin package.

        Returns
        -------
        bool
            True when ``style`` requests an advanced package scaffold.
        """
        return self.style in ADVANCED_STYLES

    @property
    def uses_contract_layout(self) -> bool:
        """Return whether the scaffold should emit an observable-contract plugin.

        Returns
        -------
        bool
            True when ``style`` requests the contract scaffold.
        """
        return self.style == CONTRACT_STYLE
