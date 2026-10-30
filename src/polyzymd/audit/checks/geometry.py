"""Geometry checks: the starting coordinates are a system that can be simulated as intended.

Distances are measured under periodic boundary conditions in the realized box,
with MDAnalysis's neighbour search (triclinic boxes included). Clash thresholds
follow the Packmol tolerance the build used: closer than half of it is a clash
(as the builder's own check counts it); between half and the full tolerance is
reported, since Packmol honours its tolerance only approximately and
minimization resolves such contacts.
"""

from __future__ import annotations

from collections import deque
from typing import Any, Iterable, Optional

import numpy as np

from polyzymd.audit.findings import AuditContext, Finding
from polyzymd.audit.registry import CheckSpec, register

#: Protein heavy atoms within BURIAL_RADIUS_NM of a water oxygen at or above which
#: the water is reported as buried. In correctly built Paper 1 systems (CALB,
#: LipA, RML; 45,000-92,000 atoms) the most buried water had 45-50.
BURIED_NEIGHBOURS = 75
BURIAL_RADIUS_NM = 0.8
#: Heavy-atom distance at which a polymer chain is in contact with the protein.
CONTACT_NM = 0.4


def _heavy(context: AuditContext, atoms: Iterable[int]) -> np.ndarray:
    elements = context.realized.elements
    return np.array([i for i in atoms if elements[i] != "H"], dtype=int)


def _atoms(context: AuditContext, *types: str) -> list[int]:
    return [i for c in context.components if c.type in types for i in c.atom_indices()]


def _box(context: AuditContext) -> np.ndarray:
    """The realized box as MDAnalysis dimensions (Å and degrees)."""
    from MDAnalysis.lib.mdamath import triclinic_box

    vectors = np.asarray(context.realized.box_nm) * 10.0
    return triclinic_box(*vectors).astype(np.float32)


def _positions_angstrom(context: AuditContext) -> np.ndarray:
    return (np.asarray(context.realized.positions_nm) * 10.0).astype(np.float32)


def _pairs(context: AuditContext, a: np.ndarray, b: np.ndarray, cutoff_angstrom: float):
    """(pairs, distances) of atoms of ``a`` and ``b`` closer than the cutoff, periodic."""
    from MDAnalysis.lib.distances import capped_distance

    if len(a) == 0 or len(b) == 0:
        return np.empty((0, 2), dtype=int), np.empty(0)
    positions = _positions_angstrom(context)
    pairs, distances = capped_distance(
        positions[a], positions[b], cutoff_angstrom, box=_box(context), return_distances=True
    )
    return pairs, distances


def _tolerances(context: AuditContext) -> tuple[float, float]:
    """Packmol tolerances (Å) of the solvent and of the polymer packing."""
    config = context.config
    if config is None:
        return 2.0, 2.0
    polymer = config.polymers.packing.tolerance if config.polymers is not None else 2.0
    return config.solvent.box.tolerance, polymer


def _overlap(
    check: str,
    severity: str,
    context: AuditContext,
    a: np.ndarray,
    b: np.ndarray,
    tolerance: float,
    what: str,
) -> list[Finding]:
    pairs, distances = _pairs(context, a, b, tolerance)
    clash = distances < 0.5 * tolerance
    close = int(np.count_nonzero(~clash))
    nearest = f"{distances.min():.2f} Å" if len(distances) else "none within tolerance"
    if clash.any():
        worst = np.argsort(distances)[:3]
        examples = [f"atoms {a[pairs[k, 0]]}-{b[pairs[k, 1]]}: {distances[k]:.2f} Å" for k in worst]
        return [
            Finding(
                check,
                "system",
                severity,  # type: ignore[arg-type]
                f"{int(clash.sum())} {what} pairs closer than {0.5 * tolerance:.2f} Å (half the "
                f"Packmol tolerance {tolerance} Å); nearest {nearest}",
                examples,
            )
        ]
    return [
        Finding(
            check,
            "system",
            "PASS",
            f"no {what} pair closer than {0.5 * tolerance:.2f} Å; {close} between that and "
            f"{tolerance} Å (minimization resolves these); nearest {nearest}",
        )
    ]


@register(
    CheckSpec(
        "geometry.water_overlap",
        "system",
        "FAIL",
        "a solvent heavy atom clashing with a solute heavy atom (closer than half the Packmol "
        "tolerance the build used)",
    )
)
def water_overlap(context: AuditContext) -> list[Finding]:
    tolerance, _ = _tolerances(context)
    solvent = _heavy(context, _atoms(context, "solvent", "ions"))
    solute = _heavy(context, _atoms(context, "protein", "small_molecule", "polymer"))
    return _overlap(
        "geometry.water_overlap", "FAIL", context, solvent, solute, tolerance, "solvent-solute"
    )


@register(
    CheckSpec(
        "geometry.polymer_overlap",
        "system",
        "FAIL",
        "a polymer heavy atom clashing with the protein or ligand",
    )
)
def polymer_overlap(context: AuditContext) -> list[Finding]:
    _, tolerance = _tolerances(context)
    polymer = _heavy(context, _atoms(context, "polymer"))
    if len(polymer) == 0:
        return []
    fixed = _heavy(context, _atoms(context, "protein", "small_molecule"))
    return _overlap(
        "geometry.polymer_overlap", "FAIL", context, polymer, fixed, tolerance, "polymer-solute"
    )


@register(
    CheckSpec(
        "geometry.interchain_overlap",
        "system",
        "WARN",
        "polymer chains clashing with each other (usually resolved by minimization; a sign of a "
        "forced Packmol result)",
    )
)
def interchain_overlap(context: AuditContext) -> list[Finding]:
    _, tolerance = _tolerances(context)
    polymers = [c for c in context.components if c.type == "polymer"]
    if not polymers:
        return []
    chain_of = {}
    for k, molecule in enumerate(polymers[0].molecules):
        for i in molecule:
            chain_of[i] = k
    heavy = _heavy(context, chain_of)
    pairs, distances = _pairs(context, heavy, heavy, 0.5 * tolerance)
    between = [
        (heavy[i], heavy[j], d)
        for (i, j), d in zip(pairs, distances)
        if chain_of[heavy[i]] < chain_of[heavy[j]]
    ]
    if between:
        nearest = min(d for _, _, d in between)
        return [
            Finding(
                "geometry.interchain_overlap",
                polymers[0].name,
                "WARN",
                f"{len(between)} atom pairs of different chains closer than {0.5 * tolerance:.2f} Å; "
                f"nearest {nearest:.2f} Å",
                [
                    f"atoms {a}-{b}: {d:.2f} Å"
                    for a, b, d in sorted(between, key=lambda x: x[2])[:3]
                ],
            )
        ]
    return [
        Finding(
            "geometry.interchain_overlap",
            polymers[0].name,
            "PASS",
            f"no two chains closer than {0.5 * tolerance:.2f} Å",
        )
    ]


@register(
    CheckSpec(
        "geometry.buried_water",
        "system",
        "WARN",
        "water buried inside the protein (may be deliberate, as a crystal water)",
        component_types=("protein",),
    )
)
def buried_water(context: AuditContext) -> list[Finding]:
    protein = _heavy(context, _atoms(context, "protein"))
    if len(protein) == 0:
        return []
    solvent = next((c for c in context.components if c.type == "solvent"), None)
    water = next((d for d in solvent.domains if d.name == "water"), None) if solvent else None
    if water is None:
        return []
    elements = context.realized.elements
    oxygens = np.array([i for r in water.atoms for i in r if elements[i] == "O"], dtype=int)
    pairs, _ = _pairs(context, oxygens, protein, BURIAL_RADIUS_NM * 10.0)
    counts = np.bincount(pairs[:, 0], minlength=len(oxygens)) if len(pairs) else np.zeros(1)
    buried = np.flatnonzero(counts >= BURIED_NEIGHBOURS)
    most = int(counts.max()) if len(counts) else 0
    detail = (
        f"most protein heavy atoms within {BURIAL_RADIUS_NM} nm of a water oxygen: {most} "
        f"(buried at {BURIED_NEIGHBOURS} or more)"
    )
    if len(buried):
        return [
            Finding(
                "geometry.buried_water",
                "system",
                "WARN",
                f"{len(buried)} waters buried in the protein; " + detail,
                [f"water oxygen {oxygens[k]}: {counts[k]} neighbours" for k in buried[:3]],
            )
        ]
    return [Finding("geometry.buried_water", "system", "PASS", detail)]


def _whole(context: AuditContext, atoms: list[int]) -> np.ndarray:
    """Positions (nm) of ``atoms``, each molecule made whole and placed in the periodic
    image nearest the first molecule, so a ligand beside the protein is measured there."""
    positions = np.asarray(context.realized.positions_nm, dtype=float)
    box = np.asarray(context.realized.box_nm, dtype=float)
    inverse = np.linalg.inv(box)
    index = context.realized.index
    members = set(atoms)
    out: dict[int, np.ndarray] = {}
    placed: set[int] = set()
    for start in atoms:
        if start in out:
            continue
        out[start] = positions[start]
        queue = deque([start])
        while queue:
            i = queue.popleft()
            for j in index.bonded_to.get(i, ()):
                if j in members and j not in out:
                    d = positions[j] - positions[i]
                    s = d @ inverse
                    out[j] = out[i] + (s - np.round(s)) @ box
                    queue.append(j)
        if start != atoms[0]:
            piece = [i for i in out if i not in placed]
            centre = np.mean([out[i] for i in piece], axis=0)
            anchor = np.mean([out[i] for i in placed], axis=0)
            s = (centre - anchor) @ inverse
            shift = -np.round(s) @ box
            for i in piece:
                out[i] = out[i] + shift
        placed.update(out)
    return np.array([out[i] for i in atoms])


def _lattice(box: np.ndarray) -> list[np.ndarray]:
    """The 26 lattice vectors to the neighbouring periodic images (nm)."""
    import itertools

    return [
        np.array(shift) @ box
        for shift in itertools.product((-1, 0, 1), repeat=3)
        if shift != (0, 0, 0)
    ]


def _image_distance(points: np.ndarray, box: np.ndarray, cutoff: float) -> Optional[float]:
    """Shortest distance (nm) between ``points`` and their periodic images, if below ``cutoff``."""
    from MDAnalysis.lib.distances import capped_distance

    best = None
    angstrom = (points * 10.0).astype(np.float32)
    for vector in _lattice(box):
        shifted = (points + vector) * 10.0
        _, distances = capped_distance(
            angstrom, shifted.astype(np.float32), cutoff * 10.0, return_distances=True
        )
        if len(distances):
            low = float(distances.min()) / 10.0
            best = low if best is None else min(best, low)
    return best


@register(
    CheckSpec(
        "geometry.periodic_image",
        "system",
        "FAIL",
        "a fixed component closer to its own periodic image than the box padding guarantees",
    )
)
def periodic_image(context: AuditContext) -> list[Finding]:
    """The protein and ligand, made whole, stay 2 x padding from their images."""
    fixed = [int(i) for i in _heavy(context, _atoms(context, "protein", "small_molecule"))]
    if not fixed:
        return []
    config = context.config
    padding = 1.2 if config is None else config.solvent.box.padding
    if config is not None and config.polymers is not None:
        # The cell reserves room for the chains as well (polymers.packing.padding).
        padding += config.polymers.packing.padding
    guaranteed = 2.0 * padding
    box = np.asarray(context.realized.box_nm, dtype=float)
    nearest = _image_distance(_whole(context, fixed), box, guaranteed)
    if nearest is not None and nearest < 0.95 * guaranteed:
        return [
            Finding(
                "geometry.periodic_image",
                "system",
                "FAIL",
                f"the solute comes within {nearest:.2f} nm of its own image; the box padding "
                f"({padding} nm) promises {guaranteed:.2f} nm",
            )
        ]
    shown = f"{nearest:.2f} nm" if nearest is not None else f"more than {guaranteed:.2f} nm"
    return [
        Finding(
            "geometry.periodic_image",
            "system",
            "PASS",
            f"the solute stays {shown} from its own image (padding {padding} nm)",
        )
    ]


@register(
    CheckSpec(
        "geometry.image_bridging",
        "system",
        "FAIL",
        "a polymer chain touching both the protein and the protein's periodic image",
        component_types=("polymer",),
    )
)
def image_bridging(context: AuditContext) -> list[Finding]:
    """No chain, made whole, contacts the protein through two different images."""
    polymers = [c for c in context.components if c.type == "polymer"]
    protein_atoms = [int(i) for i in _heavy(context, _atoms(context, "protein"))]
    if not polymers or not protein_atoms:
        return []
    from MDAnalysis.lib.distances import capped_distance

    box = np.asarray(context.realized.box_nm, dtype=float)
    protein = _whole(context, protein_atoms) * 10.0
    vectors = [np.zeros(3)] + _lattice(box)
    bridging = []
    for k, molecule in enumerate(polymers[0].molecules):
        chain_atoms = [int(i) for i in _heavy(context, molecule)]
        chain = (_whole(context, chain_atoms) * 10.0).astype(np.float32)
        touched = set()
        for n, vector in enumerate(vectors):
            image = (protein + vector * 10.0).astype(np.float32)
            pairs = capped_distance(chain, image, CONTACT_NM * 10.0, return_distances=False)
            if len(pairs):
                touched.add(n)
        if len(touched) > 1:
            bridging.append(k)
    if bridging:
        return [
            Finding(
                "geometry.image_bridging",
                polymers[0].name,
                "FAIL",
                f"{len(bridging)} chains touch the protein and one of its periodic images "
                f"(within {CONTACT_NM} nm), so they bridge copies of the protein",
                [
                    f"chain {k} starting at atom {polymers[0].molecules[k].start}"
                    for k in bridging[:3]
                ],
            )
        ]
    return [
        Finding(
            "geometry.image_bridging",
            polymers[0].name,
            "PASS",
            f"no chain of {polymers[0].n_molecules} touches the protein through two images",
        )
    ]


#: Closest two atoms of a modified protein may come without a bond or angle between
#: them, by the number of hydrogens in the pair. Each limit is shorter than a bond
#: between such atoms (heavy-heavy 1.2-1.5 Å, X-H 0.96-1.1 Å, H-H 0.74 Å) and far
#: inside a van der Waals contact (H...H near 2.4 Å, a hydrogen bond's H...O near
#: 1.8 Å), so only an atom placed through another is caught. A hydrogen placed
#: through the protein is no less wrong than a heavy atom: the minimization before
#: a run must tear it out, and that moves the frozen solute.
MODIFIED_PROTEIN_CLASH_ANGSTROM = {
    "heavy-heavy": 1.0,
    "hydrogen-heavy": 1.0,
    "hydrogen-hydrogen": 0.8,
}


@register(
    CheckSpec(
        "geometry.modified_protein_clash",
        "component",
        "FAIL",
        "atoms of a modified protein without a bond or angle between them closer than "
        "1.0 Å (heavy-heavy, hydrogen-heavy) or 0.8 Å (hydrogen-hydrogen): a fragment or "
        "side chain placed through the protein, which minimization can only fix by "
        "stretching bonds or moving the protein",
        ("protein",),
    )
)
def modified_protein_clash(context: AuditContext) -> list[Finding]:
    from polyzymd.audit.modified import is_modified

    findings = []
    index = context.realized.index
    elements = context.realized.elements
    topology_atoms = list(context.realized.topology.atoms())
    limits = MODIFIED_PROTEIN_CLASH_ANGSTROM
    kinds = ("heavy-heavy", "hydrogen-heavy", "hydrogen-hydrogen")
    for component in context.components:
        if component.type != "protein" or not is_modified(component):
            continue
        atoms = np.asarray(component.atom_indices(), dtype=int)
        pairs, distances = _close_pairs(context, atoms, max(limits.values()))
        clashes = []
        for (a, b), distance in zip(pairs, distances):
            i, j = int(atoms[a]), int(atoms[b])
            if i >= j:
                continue
            kind = kinds[(elements[i] == "H") + (elements[j] == "H")]
            if distance >= limits[kind]:
                continue
            near = index.bonded_to.get(i, set())
            if j in near or any(j in index.bonded_to.get(k, set()) for k in near):
                continue  # bonded, or the ends of an angle
            clashes.append((float(distance), i, j, kind))
        if clashes:
            clashes.sort()

            def label(i: int) -> str:
                atom = topology_atoms[i]
                return f"{atom.residue.name}{atom.residue.id}:{atom.name}"

            counts = ", ".join(
                f"{n} {kind} pairs closer than {limits[kind]} Å"
                for kind in kinds
                if (n := sum(c[3] == kind for c in clashes))
            )
            worst = "; ".join(f"{label(i)}-{label(j)} {d:.2f} Å" for d, i, j, _ in clashes[:5])
            findings.append(
                Finding(
                    "geometry.modified_protein_clash",
                    component.name,
                    "FAIL",
                    f"{counts} without a bond or angle between them (closest: {worst})",
                )
            )
        else:
            findings.append(
                Finding(
                    "geometry.modified_protein_clash",
                    component.name,
                    "PASS",
                    "no atom pair without a bond or angle between them closer than "
                    f"{limits['heavy-heavy']} Å ({limits['hydrogen-hydrogen']} Å for two "
                    "hydrogens)",
                )
            )
    return findings


def _close_pairs(context: AuditContext, atoms: np.ndarray, cutoff_angstrom: float):
    """Pairs of ``atoms`` closer than the cutoff, periodic when the structure has a box."""
    from MDAnalysis.lib.distances import capped_distance

    if len(atoms) == 0:
        return np.empty((0, 2), dtype=int), np.empty(0)
    box = _box(context) if context.realized.box_nm is not None else None
    positions = _positions_angstrom(context)[atoms]
    return capped_distance(positions, positions, cutoff_angstrom, box=box, return_distances=True)
