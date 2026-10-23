"""Read every parameter of an OpenMM System once, and compare sets of terms.

A term is keyed by the atoms it acts on, in the indices of whatever it is being
compared with, and carries its parameters in OpenMM's units (nm, kJ/mol, rad,
elementary charge, dalton).
"""

from __future__ import annotations

import collections
import math
from dataclasses import dataclass, field
from typing import Any, Callable, Iterable, Mapping, Optional

import numpy as np

#: Every kind of term the audit compares, in report order.
TERM_KINDS = ("bonds", "angles", "torsions", "lj", "exceptions", "constraints", "mass", "charge")

#: Forces that change how a system is run but hold no force-field parameters.
RUN_CONTROL_FORCES = frozenset(
    {
        "CMMotionRemover",
        "MonteCarloBarostat",
        "MonteCarloAnisotropicBarostat",
        "MonteCarloFlexibleBarostat",
        "MonteCarloMembraneBarostat",
        "AndersenThermostat",
    }
)

#: Appended to each refusal of something the audit cannot read yet, so a user
#: whose system needs it asks for it rather than works around the check.
ASK_FOR_SUPPORT = (
    "If your system needs it, open an issue at "
    "https://github.com/joelaforet/polyzymd/issues so it can be supported."
)

#: Forces PolyzyMD's restraints are made of. They are set aside, not compared as
#: force-field terms; the restraint checks account for each one.
RESTRAINT_FORCES = frozenset(
    {"CustomBondForce", "CustomCentroidBondForce", "CustomExternalForce", "CustomCompoundBondForce"}
)


@dataclass
class SystemIndex:
    """The terms of an OpenMM System, grouped by the molecule of their first atom.

    Built once per System, so that pulling the terms of one molecule out of a
    100,000-atom system costs only the size of that molecule.
    """

    molecule_of: list[int]
    terms: dict[int, dict[str, list[tuple[tuple[int, ...], Any]]]]
    particle: list[tuple[float, tuple[float, float], float]]
    bonded_to: dict[int, set[int]] = field(default_factory=dict)
    restraint_forces: list[int] = field(default_factory=list)

    @classmethod
    def from_system(cls, system: Any) -> "SystemIndex":
        """Index ``system``; molecules are the connected components of bonds and constraints.

        Raises:
            ValueError: If ``system`` has a force this index cannot read (for
                example a CustomNonbondedForce, CMAP or Ryckaert-Bellemans
                torsions), or more than one NonbondedForce. Skipping such a
                force would let its parameters pass every comparison unread.
        """
        import openmm
        from openmm import unit

        n = system.getNumParticles()
        parent = list(range(n))

        def find(i: int) -> int:
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        raw: dict[str, list[tuple[tuple[int, ...], Any]]] = collections.defaultdict(list)
        bonded_to: dict[int, set[int]] = collections.defaultdict(set)
        for k in range(system.getNumConstraints()):
            a, b, d = system.getConstraintParameters(k)
            raw["constraints"].append(((a, b), d.value_in_unit(unit.nanometer)))
        nonbonded = None
        restraint_forces, unreadable = [], []
        for number, force in enumerate(system.getForces()):
            name = force.__class__.__name__
            if name in RUN_CONTROL_FORCES:
                continue
            if name in RESTRAINT_FORCES:
                restraint_forces.append(number)
            elif isinstance(force, openmm.HarmonicBondForce):
                for k in range(force.getNumBonds()):
                    a, b, r, kk = force.getBondParameters(k)
                    raw["bonds"].append(((a, b), (r._value, kk._value)))
                    parent[find(a)] = find(b)
                    bonded_to[a].add(b)
                    bonded_to[b].add(a)
            elif isinstance(force, openmm.HarmonicAngleForce):
                for k in range(force.getNumAngles()):
                    a, b, c, theta, kk = force.getAngleParameters(k)
                    raw["angles"].append(((a, b, c), (theta._value, kk._value)))
            elif isinstance(force, openmm.PeriodicTorsionForce):
                for k in range(force.getNumTorsions()):
                    a, b, c, d, per, phase, kk = force.getTorsionParameters(k)
                    if abs(kk._value) > 1e-10:
                        raw["torsions"].append(
                            (
                                (a, b, c, d),
                                (per, round(phase._value % (2 * math.pi), 5), round(kk._value, 6)),
                            )
                        )
            elif isinstance(force, openmm.NonbondedForce) and nonbonded is None:
                nonbonded = force
                for k in range(force.getNumExceptions()):
                    a, b, qq, sigma, eps = force.getExceptionParameters(k)
                    raw["exceptions"].append(
                        (
                            (a, b),
                            (qq._value, sigma._value if eps._value != 0 else 0.0, eps._value),
                        )
                    )
            else:
                unreadable.append(f"{name} (force {number})")
        if unreadable:
            raise ValueError(
                "the system has forces the audit cannot read: "
                + ", ".join(unreadable)
                + ". PolyzyMD builds with one combined NonbondedForce and harmonic or "
                "periodic valence forces. " + ASK_FOR_SUPPORT
            )
        if nonbonded is None:
            raise ValueError("the system has no NonbondedForce")
        for (a, b), _ in raw["constraints"]:
            parent[find(a)] = find(b)
            bonded_to[a].add(b)
            bonded_to[b].add(a)
        molecule_of = [find(i) for i in range(n)]
        terms: dict[int, dict[str, list]] = collections.defaultdict(
            lambda: collections.defaultdict(list)
        )
        for kind, entries in raw.items():
            for atoms, value in entries:
                terms[molecule_of[atoms[0]]][kind].append((atoms, value))
        particle = []
        for i in range(n):
            q, sigma, eps = nonbonded.getParticleParameters(i)
            lj = (sigma._value, eps._value) if eps._value != 0 else (0.0, 0.0)
            particle.append((q._value, lj, system.getParticleMass(i).value_in_unit(unit.dalton)))
        return cls(molecule_of, dict(terms), particle, dict(bonded_to), restraint_forces)

    def extract(
        self, atoms: Iterable[int], index_map: Optional[Mapping[int, int]] = None
    ) -> dict[str, dict]:
        """Return every term wholly inside ``atoms``, keyed in ``index_map``'s indices."""
        atoms = set(atoms)
        m: Callable[[int], int] = (lambda i: i) if index_map is None else index_map.__getitem__
        out: dict[str, dict] = {kind: {} for kind in TERM_KINDS}
        torsions: dict[tuple, list] = collections.defaultdict(list)
        for i in atoms:
            q, lj, mass = self.particle[i]
            out["charge"][m(i)] = q
            out["lj"][m(i)] = lj
            out["mass"][m(i)] = mass
        for molecule in {self.molecule_of[i] for i in atoms}:
            for kind, entries in self.terms.get(molecule, {}).items():
                for term_atoms, value in entries:
                    if not atoms.issuperset(term_atoms):
                        continue
                    mapped = tuple(m(i) for i in term_atoms)
                    if kind == "angles":
                        out[kind][mapped if mapped[0] < mapped[2] else mapped[::-1]] = value
                    elif kind == "torsions":
                        # A torsion read backwards is the same torsion; any other
                        # order is a different dihedral, so it keeps its own key.
                        torsions[min(mapped, mapped[::-1])].append(value)
                    else:
                        out[kind][tuple(sorted(mapped))] = value
        out["torsions"] = {key: sorted(values) for key, values in torsions.items()}
        return out


@dataclass
class Comparison:
    """How many terms of one kind agree, and the first few that do not."""

    kind: str
    matched: int
    total: int
    mismatches: list[tuple[Any, Any, Any]]

    @property
    def ok(self) -> bool:
        return self.matched == self.total


def compare(kind: str, actual: Mapping, expected: Mapping) -> Comparison:
    """Compare two term dicts over the union of their keys.

    Charges must agree to 1e-5 e and every other parameter to a relative 1e-4
    (absolute 1e-6); a term present on one side only is a mismatch.
    """
    atol = 1e-5 if kind == "charge" else 1e-6
    keys = set(actual) | set(expected)
    bad = []
    for key in keys:
        a, e = actual.get(key), expected.get(key)
        if a is None or e is None:
            bad.append((key, a, e))
            continue
        a_arr = np.atleast_1d(np.array(a, dtype=float))
        e_arr = np.atleast_1d(np.array(e, dtype=float))
        if a_arr.shape != e_arr.shape or not np.allclose(a_arr, e_arr, rtol=1e-4, atol=atol):
            bad.append((key, a, e))
    return Comparison(kind, len(keys) - len(bad), len(keys), sorted(bad, key=repr)[:3])
