"""Components of a built system: which atoms each one owns, what the config promised for them.

A built system is a list of components (the protein, a ligand, the polymers, the
solvent, the ions). Each component records three things:

- **its atoms**: index ranges in the built system, set by the builder that placed
  them, so nothing has to recover a component's atoms from chain letters later;
- **its contract**: per-atom *domains*, each naming where the atoms' valence terms
  (bonds, angles, torsions, impropers, Lennard-Jones) and their partial charges
  come from. A term whose atoms lie in different domains takes the source of the
  domain listed later (see :data:`BOUNDARY_RULE`);
- **its provenance**: the inputs and settings it was made from.

The builders produce these records and the physics audit (:mod:`polyzymd.audit`)
checks the built system against them, so this module imports neither.
"""

from __future__ import annotations

import hashlib
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Iterator, Literal, Optional, Sequence

#: Default PSF segment ID of each built-in component type. A user-defined
#: component sets its own, which must not collide with these.
DEFAULT_SEGIDS: dict[str, str] = {
    "protein": "PROT",
    "small_molecule": "LIG",
    "polymer": "POLY",
    "solvent": "SOL",
    "ions": "ION",
}
#: PSF segment IDs hold at most four characters.
SEGID_MAX_LENGTH = 4

#: How a term whose atoms lie in different domains is assigned: it belongs to a
#: domain only if every one of its atoms does; otherwise it takes the source of
#: the domain listed later in :attr:`ComponentRecord.domains`. A modified protein
#: lists the protein first and the modification site after it, so every term that
#: touches the site comes from the site's source.
BOUNDARY_RULE = "a term belongs to a domain only if all its atoms do; otherwise the later domain"

ChargeMethod = Literal["library", "nagl", "am1bcc", "espaloma", "preset", "water_model"]


@dataclass(frozen=True)
class ChargeSource:
    """Where a domain's partial charges come from.

    ``method`` is one of:

    - ``library``: the library charges of the domain's force field (ff14SB, ion entries);
    - ``water_model``: the literature charges of the water model named in ``model``;
    - ``nagl``: NAGL, with the model file named in ``model``;
    - ``am1bcc`` / ``espaloma``: computed by that method;
    - ``preset``: read from a file (``source_file``), such as a charged polymer SDF.
    """

    method: ChargeMethod
    model: Optional[str] = None
    source_file: Optional[str] = None


@dataclass(frozen=True)
class AtomRange:
    """Atoms ``start`` up to but excluding ``stop`` of the built system."""

    start: int
    stop: int

    def __len__(self) -> int:
        return self.stop - self.start

    def __iter__(self) -> Iterator[int]:
        return iter(range(self.start, self.stop))


@dataclass(frozen=True)
class Domain:
    """Atoms that share one valence source and one charge source.

    ``reference`` says how to make the domain's molecule from its inputs alone,
    so the audit can rebuild it independently of the build: ``{"pdb": path,
    "custom_substructures": path or None}`` (all its molecules at once, in
    order), ``{"sdf": path, "conformer_index": i}`` or ``{"smiles": text}`` (one
    molecule, of which every atom range of the domain is a copy).
    """

    name: str
    atoms: tuple[AtomRange, ...]
    force_field: tuple[str, ...]
    charges: ChargeSource
    reference: Optional[dict[str, Any]] = None

    def atom_indices(self) -> list[int]:
        return [i for r in self.atoms for i in r]


@dataclass
class ComponentRecord:
    """One component of a built system.

    Attributes:
        name: The handle restraints, exceptions and analysis selections use.
        type: The component type (``protein``, ``small_molecule``, ``polymer``,
            ``solvent``, ``ions`` or a user-defined type).
        segid: The PSF segment ID written for every atom of the component.
        placement: ``fixed`` if it arrived with coordinates, ``packed`` if
            PolyzyMD placed it.
        molecules: Atom range of every molecule, in system order.
        domains: Where the atoms' parameters and charges come from.
        provenance: Inputs, their SHA-256, counts and settings.
    """

    name: str
    type: str
    segid: str
    placement: Literal["fixed", "packed"]
    molecules: list[AtomRange]
    domains: list[Domain]
    provenance: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        validate_segid(self.segid)
        owned = set(self.atom_indices())
        covered = [i for d in self.domains for i in d.atom_indices()]
        if len(covered) != len(set(covered)) or set(covered) != owned:
            raise ValueError(
                f"component {self.name!r}: its domains must cover each of its "
                f"{len(owned)} atoms exactly once (they cover {len(set(covered))}, "
                f"{len(covered) - len(set(covered))} more than once)"
            )

    @property
    def n_atoms(self) -> int:
        return sum(len(r) for r in self.molecules)

    @property
    def n_molecules(self) -> int:
        return len(self.molecules)

    def atom_indices(self) -> list[int]:
        return [i for r in self.molecules for i in r]

    def to_dict(self) -> dict[str, Any]:
        """Return the record as plain data, with atom ranges merged where contiguous."""
        data = asdict(self)
        data["n_atoms"] = self.n_atoms
        data["n_molecules"] = self.n_molecules
        data["molecules"] = [[r.start, r.stop] for r in merge_ranges(self.molecules)]
        for domain, raw in zip(self.domains, data["domains"]):
            raw["atoms"] = [[r.start, r.stop] for r in merge_ranges(domain.atoms)]
        data["boundary_rule"] = BOUNDARY_RULE
        return data


def validate_segid(segid: str, reserved: Iterable[str] = ()) -> None:
    """Raise ValueError if ``segid`` is empty, too long, not plain ASCII, or reserved."""
    if not segid or len(segid) > SEGID_MAX_LENGTH or not segid.isascii() or " " in segid:
        raise ValueError(
            f"segment ID {segid!r} must be 1-{SEGID_MAX_LENGTH} ASCII characters without "
            "spaces (the PSF limit)"
        )
    if segid in set(reserved):
        raise ValueError(f"segment ID {segid!r} is reserved for a built-in component type")


def merge_ranges(ranges: Sequence[AtomRange]) -> list[AtomRange]:
    """Merge adjacent atom ranges, so 94 consecutive polymer chains read as one range."""
    merged: list[AtomRange] = []
    for r in sorted(ranges, key=lambda r: r.start):
        if merged and merged[-1].stop == r.start:
            merged[-1] = AtomRange(merged[-1].start, r.stop)
        else:
            merged.append(r)
    return merged


def molecule_ranges(sizes: Sequence[int], start: int) -> list[AtomRange]:
    """Return consecutive atom ranges for molecules of ``sizes`` beginning at ``start``."""
    ranges = []
    for size in sizes:
        ranges.append(AtomRange(start, start + size))
        start += size
    return ranges


def file_record(path: Optional[Path | str]) -> Optional[dict[str, str]]:
    """Return the path and SHA-256 of an input file, or None if there is none."""
    if path is None:
        return None
    path = Path(path)
    if not path.is_file():
        return {"path": str(path), "sha256": None}
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def check_layout(components: Sequence[ComponentRecord], n_atoms: int) -> None:
    """Raise ValueError unless the components cover atoms 0..n_atoms-1 exactly once,
    and their segment IDs and names are unique."""
    seen: list[int] = sorted(i for c in components for i in c.atom_indices())
    if seen != list(range(n_atoms)):
        missing = sorted(set(range(n_atoms)) - set(seen))
        raise ValueError(
            f"the components cover {len(set(seen))} of {n_atoms} atoms"
            + (f" (first missing atom {missing[0]})" if missing else " (some atoms twice)")
        )
    for attribute in ("name", "segid"):
        values = [getattr(c, attribute) for c in components]
        duplicates = sorted({v for v in values if values.count(v) > 1})
        if duplicates:
            raise ValueError(f"two components share the {attribute} {duplicates[0]!r}")


def atom_segids(components: Sequence[ComponentRecord], n_atoms: int) -> list[str]:
    """Return the segment ID of every atom of the system, in atom order."""
    segids = [""] * n_atoms
    for component in components:
        for i in component.atom_indices():
            segids[i] = component.segid
    return segids


def write_pdb_segids(path: Path | str, segids: Sequence[str]) -> None:
    """Write each atom's segment ID into columns 73-76 of a PDB file, in place.

    MDAnalysis reads these columns as ``segid`` and PyMOL as ``segi``. The
    ATOM and HETATM records must be in system atom order, as OpenMM writes them.

    Raises:
        ValueError: If the file holds a different number of atoms.
    """
    path = Path(path)
    lines = path.read_text().splitlines(keepends=True)
    atom = 0
    for k, line in enumerate(lines):
        if line.startswith(("ATOM  ", "HETATM")):
            if atom >= len(segids):
                raise ValueError(f"{path} holds more than the {len(segids)} atoms given")
            body = line.rstrip("\r\n").ljust(80)
            lines[k] = body[:72] + segids[atom].ljust(4) + body[76:].rstrip() + "\n"
            atom += 1
    if atom != len(segids):
        raise ValueError(f"{path} holds {atom} atoms, not {len(segids)}")
    path.write_text("".join(lines))


def write_index(path: Path | str, components: Sequence[ComponentRecord]) -> Path:
    """Write a GROMACS index file with one group per component, named by the component.

    GROMACS atom numbers count from 1. Each group is also written under its
    segment ID when that differs from the name, so either handle selects it.
    """
    path = Path(path)
    blocks = []
    for component in components:
        numbers = [str(i + 1) for i in component.atom_indices()]
        rows = [" ".join(numbers[k : k + 15]) for k in range(0, len(numbers), 15)]
        for group in dict.fromkeys((component.name, component.segid)):
            blocks.append(f"[ {group} ]\n" + "\n".join(rows) + "\n")
    path.write_text("\n".join(blocks))
    return path
