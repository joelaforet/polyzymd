"""Shared setup for the physics audit tests.

The audit compares a build with reference parameters that OpenFF assigns to
each molecule from scratch (``polyzymd.audit.references.reference_system``).
A reference depends only on its inputs: the force field, the molecule, the
charges given and where constraints come from. It does not depend on the
build being audited. Most audit tests break one copy of the same small build
and check that the audit catches it, so they ask for the same references
again; building each takes about a second and a test needs about nine.

Here each reference is built once per distinct set of inputs and every call
gets its own copy, so a test that edits what it is given cannot change what
another test sees. The key holds every input in full (force field names and
file contents, the topology, the charge molecules with their charges, the
constraints source), so two different inputs never share a reference. A call
that raises is not cached. When a reference is first built, its copy is
checked term by term against the original.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

_REFERENCES: dict[str, str] = {}


def _file_digest(name: str) -> str | None:
    """SHA-256 of a force field given as a path; None for one OpenFF finds by name."""
    path = Path(name)
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def _key(force_field_files, topology, charge_from, constraints_from) -> str:
    files = list(force_field_files)
    parts = [
        repr(files),
        repr([_file_digest(name) for name in files]),
        topology.to_json(),
        repr([molecule.to_json() for molecule in charge_from]),
        repr(constraints_from),
        repr(_file_digest(constraints_from) if constraints_from else None),
    ]
    return hashlib.sha256("\0".join(parts).encode()).hexdigest()


@pytest.fixture(autouse=True, scope="module")
def _shared_reference_systems():
    """Build each audit reference once for all audit tests; hand out copies."""
    import openmm

    from polyzymd.audit import references
    from polyzymd.audit.terms import SystemIndex

    build = references.reference_system

    def shared(force_field_files, topology, charge_from=(), constraints_from=None):
        key = _key(force_field_files, topology, charge_from, constraints_from)
        if key not in _REFERENCES:
            system = build(
                force_field_files, topology, charge_from, constraints_from=constraints_from
            )
            xml = openmm.XmlSerializer.serialize(system)
            copy = openmm.XmlSerializer.deserialize(xml)
            assert SystemIndex.from_system(copy) == SystemIndex.from_system(
                system
            ), "an audit reference changed on its way through XML"
            _REFERENCES[key] = xml
        return openmm.XmlSerializer.deserialize(_REFERENCES[key])

    patch = pytest.MonkeyPatch()
    patch.setattr(references, "reference_system", shared)
    yield
    patch.undo()
