"""Runtime validation helpers for simulation configuration references."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def collect_reference_warnings(config: Any) -> list[str]:
    """Return warnings for missing files referenced by a simulation config.

    These checks intentionally live outside the Pydantic schema so lightweight
    configuration parsing can succeed before external structure files are staged.

    Parameters
    ----------
    config : Any
        Simulation configuration or config-like object with PolyzyMD schema
        attributes.

    Returns
    -------
    list[str]
        Human-readable warnings for missing referenced files or directories.
    """

    warnings: list[str] = []
    enzyme = getattr(config, "enzyme", None)
    _check_file(warnings, getattr(enzyme, "pdb_path", None), "enzyme PDB")

    substrate = getattr(config, "substrate", None)
    if substrate is not None:
        _check_file(warnings, getattr(substrate, "sdf_path", None), "substrate SDF")

    polymers = getattr(config, "polymers", None)
    if polymers is not None and bool(getattr(polymers, "enabled", False)):
        _check_polymer_references(warnings, polymers)

    return warnings


def _check_file(warnings: list[str], path_value: Any, label: str) -> None:
    """Append a warning when a referenced file path is missing."""

    if path_value is None:
        return
    path = Path(path_value)
    if not path.is_file():
        warnings.append(f"Missing {label}: {path}")


def _check_directory(warnings: list[str], path_value: Any, label: str) -> Path | None:
    """Append a warning when a referenced directory path is missing."""

    if path_value is None:
        return None
    path = Path(path_value)
    if not path.is_dir():
        warnings.append(f"Missing {label}: {path}")
        return None
    return path


def _check_polymer_references(warnings: list[str], polymers: Any) -> None:
    """Check polymer SDF references."""

    generation_mode = str(getattr(polymers, "generation_mode", "")).lower()
    if generation_mode.endswith("cached"):
        sdf_directory = _check_directory(
            warnings,
            getattr(polymers, "sdf_directory", None),
            "polymer SDF directory",
        )
        if sdf_directory is not None:
            _check_cached_polymer_sdfs(warnings, polymers, sdf_directory)


def _check_cached_polymer_sdfs(warnings: list[str], polymers: Any, sdf_directory: Path) -> None:
    """Warn when a cached-polymer directory has no matching SDF files."""

    type_prefix = getattr(polymers, "type_prefix", None)
    length = getattr(polymers, "length", None)
    if not type_prefix or length is None:
        return

    pattern = f"{type_prefix}_seq=*_{length}-mer_charged.sdf"
    if not any(sdf_directory.glob(pattern)):
        warnings.append(
            "Missing cached polymer SDF files: " f"no files matching {sdf_directory / pattern}"
        )
