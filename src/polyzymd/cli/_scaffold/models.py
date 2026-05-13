"""Models used by the analysis scaffold renderer."""

from __future__ import annotations

from dataclasses import dataclass

VALID_STYLES = ("measurement", "dict", "pydantic")
ADVANCED_STYLES = ("dict", "pydantic")


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
        Scaffold style. ``"measurement"`` creates the simple single-file
        contributor path; ``"dict"`` and ``"pydantic"`` create the advanced
        runner-backed package scaffold.
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
    def uses_pydantic_results(self) -> bool:
        """Return whether the scaffold should generate Pydantic results.

        Returns
        -------
        bool
            True when ``style`` is ``"pydantic"``.
        """
        return self.style == "pydantic"

    @property
    def uses_measurement_api(self) -> bool:
        """Return whether the scaffold should use the measurement API.

        Returns
        -------
        bool
            True when ``style`` is ``"measurement"``.
        """
        return self.style == "measurement"
