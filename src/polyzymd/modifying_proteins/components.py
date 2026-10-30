"""The component record of a modified protein, for the physics audit and provenance.

A covalently modified protein is one molecule, so it is one component. Its
atoms are written under several segment IDs (``PROT`` for the protein, one
per attachment: ``G001``, ``P001``, ``M001``) and split into domains: the
protein (ff14SB, library charges) and one domain per modification site (its
parameter set). Every domain's reference names what the audit needs to
rebuild the protein from its inputs alone: the input PDB file, the
modification blocks and the parameter sets.
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
        "spec": to_spec(protein),
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
        "constraints_from": constraints_from,
    }
    domains = [
        Domain(
            name="protein",
            atoms=_ranges([index + start for index in range(n_atoms) if index not in in_site]),
            force_field=tuple(record.protein_force_field or PROTEIN_FORCE_FIELD),
            charges=ChargeSource("library"),
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
    for index, segid in enumerate(record.segments):
        if segid != PROTEIN_SEGMENT:
            segments.setdefault(segid, []).append(index + start)
    return ComponentRecord(
        name=name,
        type="protein",
        segid=DEFAULT_SEGIDS["protein"],
        placement=placement,  # type: ignore[arg-type]
        molecules=[AtomRange(start, start + n_atoms)],
        domains=domains,
        provenance={
            "pdb": dict(protein.source),
            "modifications": reference["spec"],
            "parameterization": record.to_dict(),
        },
        segments={segid: tuple(_ranges(atoms)) for segid, atoms in segments.items()},
    )


def _ranges(indices: Sequence[int]) -> tuple[AtomRange, ...]:
    """Return the contiguous atom ranges of sorted ``indices``."""
    ranges: list[AtomRange] = []
    for index in sorted(indices):
        if ranges and ranges[-1].stop == index:
            ranges[-1] = AtomRange(ranges[-1].start, index + 1)
        else:
            ranges.append(AtomRange(index, index + 1))
    return tuple(ranges)
