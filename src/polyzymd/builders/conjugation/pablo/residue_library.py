"""Pablo residue-library helpers for source-backed conjugation workflows."""

from __future__ import annotations

import copy
import importlib
from typing import Any

from pydantic import BaseModel, Field


class PabloResidueLibraryError(RuntimeError):
    """Raised when a requested Pablo residue-library option is unsupported."""


class PabloResidueLibraryDiagnostic(BaseModel):
    """JSON-safe diagnostic emitted while preparing a Pablo residue library."""

    severity: str
    message: str
    details: dict[str, Any] = Field(default_factory=dict)


class PabloResidueLibraryResult(BaseModel):
    """Prepared Pablo residue library and JSON-safe diagnostics."""

    model_config = {"arbitrary_types_allowed": True}

    residue_library: Any = Field(exclude=True)
    diagnostics: list[PabloResidueLibraryDiagnostic] = Field(default_factory=list)


def build_pablo_residue_library(
    policy: Any | None,
    *,
    pablo_module: Any | None = None,
) -> PabloResidueLibraryResult:
    """Build a Pablo residue library from ``STD_CCD_CACHE`` and policy settings.

    Parameters
    ----------
    policy : Any or None
        CCD/Pablo policy configuration. The helper uses duck typing to avoid a
        hard dependency on the config module.
    pablo_module : Any or None, optional
        Imported ``openff.pablo`` module, by default ``None``. When omitted,
        OpenFF Pablo is imported lazily inside this function.

    Returns
    -------
    PabloResidueLibraryResult
        Pablo residue library plus JSON-safe diagnostics.

    Raises
    ------
    PabloResidueLibraryError
        If Pablo does not expose the required cache or crosslink API.
    """
    if pablo_module is None:
        pablo_module = importlib.import_module("openff.pablo")

    try:
        residue_library = pablo_module.STD_CCD_CACHE
    except AttributeError as exc:
        raise PabloResidueLibraryError("OpenFF Pablo does not expose STD_CCD_CACHE") from exc

    diagnostics: list[PabloResidueLibraryDiagnostic] = []
    residue_library = _apply_lookup_policy(residue_library, policy, diagnostics)

    crosslinks = list(_policy_attr(policy, "crosslinks", ()))
    if crosslinks and not hasattr(residue_library, "with_crosslink"):
        raise PabloResidueLibraryError(
            "The installed OpenFF Pablo cache does not support with_crosslink(), "
            "which is required for configured conjugation crosslinks"
        )

    for index, crosslink in enumerate(crosslinks):
        crosslink_kwargs = _crosslink_kwargs(crosslink)
        try:
            residue_library = residue_library.with_crosslink(**crosslink_kwargs)
        except TypeError as exc:
            raise PabloResidueLibraryError(
                "OpenFF Pablo rejected the configured crosslink arguments; check that the "
                "installed Pablo version supports residues, linking_atoms, leaving_atoms, "
                "and bond_order keyword arguments"
            ) from exc
        diagnostics.append(
            PabloResidueLibraryDiagnostic(
                severity="info",
                message="Configured Pablo crosslink added to residue library",
                details={"index": index, **_json_safe_crosslink_details(crosslink_kwargs)},
            )
        )

    if _policy_attr(policy, "ccd_cache_directory", None) is not None:
        diagnostics.append(
            PabloResidueLibraryDiagnostic(
                severity="warning",
                message="Custom CCD cache directories are recorded but not yet applied",
                details={
                    "ccd_cache_directory": str(_policy_attr(policy, "ccd_cache_directory", None))
                },
            )
        )
    residue_definition_files = tuple(_policy_attr(policy, "residue_definition_files", ()))
    if residue_definition_files:
        files = [str(path) for path in residue_definition_files]
        raise PabloResidueLibraryError(
            "Custom residue definition files were requested, but PolyzyMD does not yet "
            "apply them before Pablo ingestion. Product-state workflows that require "
            "custom residue definitions must stop before ingestion instead of silently "
            f"continuing with the standard CCD cache. Requested files: {files}"
        )

    return PabloResidueLibraryResult(
        residue_library=residue_library,
        diagnostics=diagnostics,
    )


def bonded_hydrogen_names(
    residue_name: str,
    atom_name: str,
    *,
    residue_library: Any | None = None,
    pablo_module: Any | None = None,
) -> tuple[str, ...]:
    """Return Pablo-template hydrogen names bonded to a residue atom.

    Parameters
    ----------
    residue_name : str
        CCD residue name to query.
    atom_name : str
        Residue atom whose bonded hydrogens are requested.
    residue_library : Any or None, optional
        Injectable Pablo residue-library object, by default ``None``.
    pablo_module : Any or None, optional
        Injectable ``openff.pablo`` module, by default ``None``.

    Returns
    -------
    tuple of str
        Sorted hydrogen atom names from authoritative Pablo residue connectivity.
    """
    if residue_library is None:
        if pablo_module is None:
            pablo_module = importlib.import_module("openff.pablo")
        residue_library = pablo_module.STD_CCD_CACHE
    definitions = _residue_definitions(residue_library, residue_name)
    if not definitions:
        raise PabloResidueLibraryError(f"Pablo residue library has no {residue_name} definitions")

    observed_sets = []
    for definition in definitions:
        atom_symbols = {
            str(getattr(atom, "name", "")).upper(): str(getattr(atom, "symbol", "")).upper()
            for atom in getattr(definition, "atoms", ())
        }
        bonded = set()
        for bond in getattr(definition, "bonds", ()):  # noqa: B007
            left = str(getattr(bond, "atom1", "")).upper()
            right = str(getattr(bond, "atom2", "")).upper()
            if left == atom_name.upper() and atom_symbols.get(right) == "H":
                bonded.add(right)
            if right == atom_name.upper() and atom_symbols.get(left) == "H":
                bonded.add(left)
        if bonded:
            observed_sets.append(frozenset(bonded))
    unique_sets = set(observed_sets)
    if len(unique_sets) != 1:
        raise PabloResidueLibraryError(
            f"Pablo residue library has inconsistent {residue_name}:{atom_name} bonded H names: "
            f"{sorted(tuple(sorted(item)) for item in unique_sets)}"
        )
    return tuple(sorted(next(iter(unique_sets))))


def _residue_definitions(residue_library: Any, residue_name: str) -> tuple[Any, ...]:
    """Return residue definitions from a Pablo-like residue library."""
    key = residue_name.upper()
    if hasattr(residue_library, "get"):
        definitions = residue_library.get(key)
    else:
        definitions = residue_library[key]
    if definitions is None:
        return ()
    if isinstance(definitions, tuple):
        return definitions
    if isinstance(definitions, list):
        return tuple(definitions)
    return (definitions,)


def _apply_lookup_policy(
    residue_library: Any,
    policy: Any | None,
    diagnostics: list[PabloResidueLibraryDiagnostic],
) -> Any:
    """Apply CCD lookup policy without mutating the global Pablo cache."""
    lookup_policy = _policy_attr(policy, "lookup_policy", None)
    lookup_value = getattr(lookup_policy, "value", lookup_policy)
    if lookup_value != "auto_download":
        diagnostics.append(
            PabloResidueLibraryDiagnostic(
                severity="info",
                message="Pablo residue library uses the configured CCD lookup policy",
                details={"lookup_policy": lookup_value or "pablo_default"},
            )
        )
        return residue_library

    derived_library = _derive_residue_library(residue_library)
    try:
        derived_library.auto_download = True
    except Exception as exc:  # noqa: BLE001 - third-party cache types vary
        raise PabloResidueLibraryError(
            "The installed OpenFF Pablo cache does not expose a writable auto_download policy"
        ) from exc

    diagnostics.append(
        PabloResidueLibraryDiagnostic(
            severity="info",
            message="Pablo residue library CCD auto-download policy enabled",
            details={"lookup_policy": "auto_download"},
        )
    )
    return derived_library


def _derive_residue_library(residue_library: Any) -> Any:
    """Return a derived cache object so global Pablo state is not modified."""
    if hasattr(residue_library, "with_"):
        try:
            return residue_library.with_({})
        except TypeError:
            try:
                return residue_library.with_()
            except TypeError as exc:
                raise PabloResidueLibraryError(
                    "The installed OpenFF Pablo cache exposes with_(), but not a supported "
                    "copy/derive signature"
                ) from exc
    try:
        return copy.copy(residue_library)
    except Exception as exc:  # noqa: BLE001 - third-party cache types vary
        raise PabloResidueLibraryError(
            "The installed OpenFF Pablo cache cannot be copied before applying lookup policy"
        ) from exc


def _crosslink_kwargs(crosslink: Any) -> dict[str, Any]:
    """Extract Pablo ``with_crosslink`` keyword arguments from a config object."""
    residues = tuple(_policy_attr(crosslink, "residues", ()))
    linking_atoms = tuple(_policy_attr(crosslink, "linking_atoms", ()))
    leaving_atoms = tuple(tuple(group) for group in _policy_attr(crosslink, "leaving_atoms", ()))
    return {
        "residues": residues,
        "linking_atoms": linking_atoms,
        "leaving_atoms": leaving_atoms,
        "bond_order": _policy_attr(crosslink, "bond_order", 1),
    }


def _json_safe_crosslink_details(crosslink_kwargs: dict[str, Any]) -> dict[str, Any]:
    """Return JSON-safe crosslink details for diagnostics."""
    return {
        "residues": list(crosslink_kwargs["residues"]),
        "linking_atoms": list(crosslink_kwargs["linking_atoms"]),
        "leaving_atoms": [list(group) for group in crosslink_kwargs["leaving_atoms"]],
        "bond_order": crosslink_kwargs["bond_order"],
    }


def _policy_attr(policy: Any | None, name: str, default: Any) -> Any:
    """Return a duck-typed policy attribute."""
    if policy is None:
        return default
    return getattr(policy, name, default)
