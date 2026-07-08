"""Public OpenMM orchestration for conjugate relaxation."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from polyzymd.builders.conjugation.relaxation._openmm_workflows import (
    run_conjugate_relaxation_workflow,
)
from polyzymd.builders.conjugation.relaxation.models import (
    ConjugateRelaxationResult,
    ConjugateRelaxationSettings,
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
