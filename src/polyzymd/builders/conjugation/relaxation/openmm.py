"""Public OpenMM orchestration for conjugate relaxation and validation."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from polyzymd.builders.conjugation.relaxation._openmm_workflows import (
    run_conjugate_relaxation_workflow,
    run_openmm_validation_workflow,
)
from polyzymd.builders.conjugation.relaxation.models import (
    ConjugateRelaxationResult,
    ConjugateRelaxationSettings,
    OpenMMValidationResult,
    OpenMMValidationSettings,
)


def validate_openmm_product(
    interchange: Any,
    output_dir: Path | str,
    *,
    protein_heavy_atom_indices: tuple[int, ...] | None = None,
    settings: OpenMMValidationSettings | None = None,
    crosslinked_pdb_path: Path | str | None = None,
    attachment_specs: tuple[Any, ...] = (),
) -> OpenMMValidationResult:
    """Run restrained minimization and short OpenMM product validation.

    Parameters
    ----------
    interchange : Any
        OpenFF Interchange-like object exposing OpenMM conversion methods.
    output_dir : pathlib.Path or str
        Directory for validation artifacts.
    protein_heavy_atom_indices : tuple of int or None, optional
        Protein heavy-atom indices to restrain when all-heavy restraints are
        disabled. When all-heavy restraints are enabled, the complete
        non-hydrogen atom set is inferred from the OpenMM topology.
    settings : OpenMMValidationSettings or None, optional
        Validation settings, by default ``None``.
    crosslinked_pdb_path : pathlib.Path, str, or None, optional
        Product PDB used to measure resolved crosslink lengths, by default
        ``None``.
    attachment_specs : tuple of Any, optional
        Resolved attachment build specs used for crosslink-specific diagnostics,
        by default ``()``.

    Returns
    -------
    OpenMMValidationResult
        Energies and artifact paths for the OpenMM validation.
    """
    return run_openmm_validation_workflow(
        interchange,
        output_dir,
        protein_heavy_atom_indices=protein_heavy_atom_indices,
        settings=settings,
        crosslinked_pdb_path=crosslinked_pdb_path,
        attachment_specs=attachment_specs,
    )


def relax_conjugate(
    interchange: Any,
    output_dir: Path | str,
    *,
    product_pdb_path: Path | str,
    attachment_specs: tuple[Any, ...],
    assembly: Any | None = None,
    settings: ConjugateRelaxationSettings | None = None,
) -> ConjugateRelaxationResult:
    """Relax a product-state conjugate with staged minimization then fixed-protein MD.

    Parameters
    ----------
    interchange : Any
        OpenFF Interchange-like object exposing OpenMM conversion methods.
    output_dir : pathlib.Path or str
        Directory for relaxation artifacts.
    product_pdb_path : pathlib.Path or str
        Product PDB carrying assembly serial and residue mapping metadata.
    attachment_specs : tuple of Any
        Resolved attachment build specs used to resolve generic product linkage atoms.
    assembly : Any or None, optional
        Assembly result with ``added_conect_pairs`` and residue mappings, by default ``None``.
    settings : ConjugateRelaxationSettings or None, optional
        Relaxation settings, by default ``None``.

    Returns
    -------
    ConjugateRelaxationResult
        Energies, validation metrics, and artifact paths.
    """
    return run_conjugate_relaxation_workflow(
        interchange,
        output_dir,
        product_pdb_path=product_pdb_path,
        attachment_specs=attachment_specs,
        assembly=assembly,
        settings=settings,
    )
