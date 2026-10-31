"""Polymer chains grafted onto a protein.

A graft is a chain built by PolyzyMD's star-sited polymer generator, from
the same monomer and cap definitions a ``type: polymer`` solute uses, and
bonded to the protein with ``attach(type="polymer", ...)``. The chain enters
mBuild as a fragment of residues, one per repeat unit (and per cap that has
its own residue name). ``attach_at`` names the position that bonds,
counted over the whole chain from 1 at the head cap, as in a ``sequence``
block; the reaction stated with the attachment decides which of its atoms
bond and which leave.

The chain is described as plain data::

    {
        "monomers": [{"label": "A", "name": "SBMA", "residue_name": "SBM",
                      "smiles": "[*]CC([*])(C)C(=O)O...", "probability": 1.0}, ...],
        "end_groups": {"head": [{"label": "H", "smiles": "[*][H]", "probability": 1.0}],
                       "tail": [...]},
        "length": 9,                      # counts the two caps, as for a solute
        "sequence": {"fixed": {5: "NHS"}},  # optional, as for a solute: the middle
        "random_seed": 1,                 # optional: one chain for every replicate
        "force_field": ["openff-2.3.0.offxml"],  # relaxes the chain in vacuum
    }

A chain without ``random_seed`` is drawn anew for each replicate
(``"replicate": n``, which ``polyzymd build`` adds). The drawn sequence (head
cap, repeat units, tail cap) is recorded with the fragment, so applying the
record again builds the same chain: the generator seeds the chain's
coordinates from its sequence.
"""

from __future__ import annotations

import hashlib
import random
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any

#: Keys a polymer fragment description may hold.
POLYMER_KEYS = frozenset(
    {
        "monomers",
        "end_groups",
        "length",
        "sequence",
        "random_seed",
        "replicate",
        "force_field",
        "drawn",
    }
)


def chain_seed(basis: int | str, site: str) -> int:
    """Return the seed of the chain grafted at ``site``, derived from ``basis``.

    ``basis`` is the entry's ``random_seed`` (one conjugate for every
    replicate) or ``"replicate <n>"`` (a new draw per replicate). One entry
    with several sites then draws an independent chain at each, and the same
    entry draws the same chains every time it is built.
    """
    digest = hashlib.sha256(f"{basis}:{site}".encode("utf-8")).hexdigest()
    return int(digest[:8], 16)


def cap_labels(polymer: dict[str, Any]) -> dict[str, dict[str, str]]:
    """Return, for each end, every cap label and name mapped to its label."""
    return {
        end: {
            **{o["label"]: o["label"] for o in polymer["end_groups"][end]},
            **{o["name"]: o["label"] for o in polymer["end_groups"][end] if o.get("name")},
        }
        for end in ("head", "tail")
    }


def chain_pattern(polymer: dict[str, Any]) -> list[str | None]:
    """Return the label at each position 1..length of the chain, None where drawn."""
    from polyzymd.config.schema import PolymerSequence

    monomers = polymer["monomers"]
    label_of = {m["label"]: m["label"] for m in monomers}
    label_of.update({m["name"]: m["label"] for m in monomers if m.get("name")})
    return PolymerSequence(**(polymer.get("sequence") or {})).pattern(
        int(polymer["length"]), label_of, cap_labels(polymer)
    )


def draw_chain(polymer: dict[str, Any], site: str | None = None) -> str:
    """Return the chain's sequence, one label per position: head cap, repeat units, tail cap.

    A recorded ``drawn`` sequence is returned as it is. Otherwise the head
    cap, the repeat units and the tail cap are drawn in that order, as for the
    chains of a polymer solute, keeping the positions ``sequence`` pins. The
    draw is seeded with ``random_seed`` when the entry gives one (the same
    chain in every replicate), else with ``replicate`` (a new chain per
    replicate); for a ``site`` the seed is :func:`chain_seed` of that.
    """
    from polyzymd.builders.polymer import draw_units

    if polymer.get("drawn"):
        return polymer["drawn"]
    monomers = polymer["monomers"]
    pattern = chain_pattern(polymer)
    basis: int | str | None = polymer.get("random_seed")
    if basis is None and polymer.get("replicate") is not None:
        basis = f"replicate {polymer['replicate']}"
    drawn_positions = None in pattern
    if drawn_positions and basis is None:
        raise ValueError(
            "a grafted chain with drawn positions needs a seed: the replicate number (given "
            "by polyzymd build) or `random_seed` for one chain in every replicate"
        )
    if basis is None:
        seed = 0
    elif site is not None:
        seed = chain_seed(basis, site)
    elif isinstance(basis, int):
        seed = basis
    else:
        seed = chain_seed(basis, "")

    def cap(end: str) -> str:
        options = polymer["end_groups"][end]
        return random.choices(
            [option["label"] for option in options],
            weights=[option["probability"] for option in options],
            k=1,
        )[0]

    state = random.getstate()
    random.seed(seed)
    try:
        head = pattern[0] or cap("head")
        units = draw_units(
            pattern[1:-1],
            [m["label"] for m in monomers],
            [float(m.get("probability", 0.0)) for m in monomers],
        )
        tail = pattern[-1] or cap("tail")
    finally:
        random.setstate(state)
    return head + units + tail


def anchor_residue(polymer: dict[str, Any], attach_at: int) -> int:
    """Return the residue number, in the fragment, of position ``attach_at`` (1..length).

    Positions 1 and ``length`` are the caps: a cap is attached to only when it
    carries a monomer's label (an end unit), since the bond then forms on that
    unit's atoms.
    """
    length = int(polymer["length"])
    if not 1 <= attach_at <= length:
        raise ValueError(
            f"attach_at {attach_at} is outside the chain's positions 1..{length} "
            f"(the caps are 1 and {length})"
        )
    sequence = draw_chain(polymer)
    head = next(o for o in polymer["end_groups"]["head"] if o["label"] == sequence[0])
    tail = next(o for o in polymer["end_groups"]["tail"] if o["label"] == sequence[-1])
    monomer_labels = {m["label"] for m in polymer["monomers"]}
    for position, end, cap in ((1, "head", head), (length, "tail", tail)):
        if attach_at == position and cap["label"] not in monomer_labels:
            raise ValueError(
                f"attach_at {attach_at} is the {end} cap {cap['label']!r}, not a monomer; "
                f"attach at 2..{length - 1}, or give the {end} a cap that is a monomer's end unit"
            )
    head_residue = 1 if head.get("residue_name") else 0
    if attach_at == 1:
        return 1
    if attach_at == length:
        return (length - 2) + head_residue + (1 if tail.get("residue_name") else 0)
    return (attach_at - 1) + head_residue


def build_chain(polymer: dict[str, Any], sequence: str):
    """Build the chain of ``sequence`` with the star-sited generator (an OpenFF Molecule)."""
    from polyzymd.builders.polymer_generator import PolymerGenerator

    unknown = sorted(set(polymer) - POLYMER_KEYS)
    if unknown:
        raise ValueError(f"unknown polymer fragment keys {unknown}; use {sorted(POLYMER_KEYS)}")
    monomers = [
        (m["label"], m["name"], m.get("residue_name") or m["name"][:3].upper(), m["smiles"])
        for m in polymer["monomers"]
    ]
    with tempfile.TemporaryDirectory(prefix="polyzymd-graft-") as cache:
        generator = PolymerGenerator.from_smiles(
            monomers,
            polymer["end_groups"]["head"],
            polymer["end_groups"]["tail"],
            cache_directory=Path(cache),
            force_field=list(polymer["force_field"]),
            charger_type="nagl",
        )
        return generator.build_molecule(sequence)


def fragment_from_molecule(molecule) -> Any:
    """Return an OpenFF molecule with residue metadata as an mBuild fragment of residues.

    Atoms are named element + count within their residue (C1, C2, H1, ...),
    bond orders and formal charges come from the molecule, so nothing is
    perceived again.
    """
    from mbuild import Compound
    from mbuild.biopolymers.residue import Residue
    from openff.units import unit

    xyz = molecule.conformers[0].m_as(unit.nanometer)
    fragment = Compound(name="PLY")
    residues: dict[int, Any] = {}
    counts: dict[int, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    particles = []
    for atom, position in zip(molecule.atoms, xyz):
        number = int(atom.metadata["residue_number"])
        if number not in residues:
            residues[number] = Residue(
                resname=atom.metadata["residue_name"], resnum=number, hetatm=True
            )
            residues[number].atom_formal_charges = {}
        counts[number][atom.symbol] += 1
        name = f"{atom.symbol}{counts[number][atom.symbol]}"
        particle = Compound(name=name, element=atom.symbol, pos=position)
        residues[number].add(particle)
        particles.append(particle)
        charge = int(atom.formal_charge.m_as(unit.elementary_charge))
        if charge:
            residues[number].atom_formal_charges[name] = charge
    for number in sorted(residues):
        residue = residues[number]
        residue.formal_charge = sum(residue.atom_formal_charges.values())
        fragment.add(residue)
    for bond in molecule.bonds:
        fragment.add_bond(
            (particles[bond.atom1_index], particles[bond.atom2_index]),
            bond_order=float(bond.bond_order),
        )
    return fragment


def fragment_from_polymer(polymer: dict[str, Any]) -> tuple[Any, dict[str, Any]]:
    """Build a grafted chain; return its fragment and its source with the drawn sequence."""
    sequence = draw_chain(polymer)
    fragment = fragment_from_molecule(build_chain(polymer, sequence))
    source = {key: value for key, value in polymer.items() if key != "drawn"}
    source["drawn"] = sequence
    return fragment, {"polymer": source}
