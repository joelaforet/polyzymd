"""The component record of a modified protein, for the physics audit and provenance.

A covalently modified protein is one molecule, so it is one component. Its
atoms are written under several segment IDs (``PROT`` for the protein, one
per attachment: ``G001``, ``P001``, ``M001``) and split into domains: the
protein (ff14SB, library charges) and one domain per modification site (its
parameter set). A protein on OpenFF Rosemary alpha is one domain (Rosemary,
NAGL charges on the whole molecule), modifications included. Every domain's
reference names what the audit needs to rebuild the protein from its inputs
alone: the input PDB file, the modification blocks and the parameter sets.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

from polyzymd.core.components import (
    DEFAULT_SEGIDS,
    AtomRange,
    ChargeSource,
    ComponentRecord,
    Domain,
)
from polyzymd.modifying_proteins.api import ModifiedProtein
from polyzymd.modifying_proteins.parameters import (
    CONSTRAINTS_FROM,
    PROTEIN_FORCE_FIELD,
    ParameterizationRecord,
    ParameterSet,
)
from polyzymd.modifying_proteins.spec import to_spec
from polyzymd.modifying_proteins.topology import PROTEIN_SEGMENT


def protein_component(
    protein: ModifiedProtein,
    record: ParameterizationRecord,
    parameter_sets: Mapping[str, ParameterSet],
    *,
    name: str = "protein",
    segid: str = DEFAULT_SEGIDS["protein"],
    start: int = 0,
    base_dir: str | Path = ".",
    constraints_from: str = CONSTRAINTS_FROM,
    placement: str = "fixed",
) -> ComponentRecord:
    """Return the component record of ``protein`` parameterized as ``record`` says.

    ``start`` is the index of the protein's first atom in the built system;
    ``base_dir`` is where relative paths in the modification blocks resolve.
    """
    n_atoms = len(record.segments)
    in_site = {index for site in record.sites for index in site.atom_indices}
    reference = {
        "pdb": protein.source["pdb_path"],
        "sha256": protein.source["sha256"],
        "ccd_paths": [str(path) for path in protein._ccd_paths],
        "download": protein._download,
        # Paths relative to the config's folder, with each file's SHA-256, so the
        # protein can be rebuilt (and re-audited) after the folder moves.
        "spec": to_spec(protein, relative_to=base_dir),
        "files": _input_files(protein, base_dir),
        "base_dir": str(Path(base_dir).resolve()),
        "parameter_sets": {
            site_type: {
                "force_field": list(parameter_set.force_field),
                "charges": dict(parameter_set.charges),
            }
            for site_type, parameter_set in sorted(parameter_sets.items())
        },
        "max_charge_shift": record.max_charge_shift,
        "protein_force_field": list(record.protein_force_field),
        "protein_charges": dict(record.protein_charges) if record.protein_charges else None,
        "constraints_from": constraints_from,
    }
    whole = record.protein_charges
    domains = [
        Domain(
            name="protein",
            atoms=_ranges([index + start for index in range(n_atoms) if index not in in_site]),
            force_field=tuple(record.protein_force_field or PROTEIN_FORCE_FIELD),
            # Rosemary alpha charges the whole protein with NAGL; ff14SB has library charges.
            charges=(
                ChargeSource("nagl", model=whole.get("model")) if whole else ChargeSource("library")
            ),
            reference={"modified_protein": reference, "site": None, "site_type": None},
        )
    ]
    for site in record.sites:
        charges = site.charges
        source = (
            ChargeSource("nagl", model=charges.get("model"))
            if charges.get("method") == "nagl"
            else ChargeSource("library")
        )
        domains.append(
            Domain(
                name=site.name,
                atoms=_ranges([index + start for index in site.atom_indices]),
                force_field=tuple(site.force_field),
                charges=source,
                reference={
                    "modified_protein": reference,
                    "site": site.name,
                    "site_type": site.type,
                },
            )
        )
    segments: dict[str, list[int]] = {}
    for index, segment in enumerate(record.segments):
        if segment != PROTEIN_SEGMENT:
            segments.setdefault(segment, []).append(index + start)
    return ComponentRecord(
        name=name,
        type="protein",
        segid=segid,
        placement=placement,  # type: ignore[arg-type]
        molecules=[
            AtomRange(start + first, start + stop)
            for first, stop in (record.molecules or [[0, n_atoms]])
        ],
        domains=domains,
        provenance={
            "pdb": dict(protein.source),
            "modifications": reference["spec"],
            "parameterization": record.to_dict(),
        },
        segments={key: tuple(_ranges(atoms)) for key, atoms in segments.items()},
    )


def _input_files(protein: ModifiedProtein, base_dir: str | Path) -> dict[str, str]:
    """Each input file of the modifications, relative to ``base_dir``, with its SHA-256."""
    import os

    files: dict[str, str] = {}

    def visit(value) -> None:
        if isinstance(value, dict):
            for key in ("pdb_path", "ccd_file"):
                if key in value and value.get("sha256"):
                    path = os.path.relpath(Path(value[key]).resolve(), Path(base_dir).resolve())
                    files[path] = value["sha256"]
            for item in value.values():
                visit(item)

    for record in protein.records:
        visit(record.arguments)
    pdb = protein.source["pdb_path"]
    files[os.path.relpath(Path(pdb).resolve(), Path(base_dir).resolve())] = protein.source["sha256"]
    return files


def _ranges(indices: Sequence[int]) -> tuple[AtomRange, ...]:
    """Return the contiguous atom ranges of sorted ``indices``."""
    ranges: list[AtomRange] = []
    for index in sorted(indices):
        if ranges and ranges[-1].stop == index:
            ranges[-1] = AtomRange(ranges[-1].start, index + 1)
        else:
            ranges.append(AtomRange(index, index + 1))
    return tuple(ranges)
