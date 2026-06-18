"""Tests for conjugation source-protein preparation helpers."""

from __future__ import annotations

import sys
from types import ModuleType, SimpleNamespace

from polyzymd.builders.conjugation.protein_preparation import (
    ProteinCanonicalizationSettings,
    canonicalize_protein_hydrogens,
)


def test_canonicalize_protein_hydrogens_reports_configured_ph(monkeypatch, tmp_path):
    """The helper should surface the pH used for OpenMM hydrogen regeneration."""
    source = tmp_path / "source.pdb"
    output = tmp_path / "source_canonical.pdb"
    source.write_text("ATOM      1  N   LYS A  23       0.000   0.000   0.000  1.00  0.00           N\nEND\n")

    calls = {}

    class FakePDBFixer:
        def __init__(self, filename):
            calls["filename"] = filename
            self.topology = SimpleNamespace(atoms=lambda: ())
            self.positions = []

        def findMissingResidues(self):
            calls["findMissingResidues"] = True

        def findNonstandardResidues(self):
            calls["findNonstandardResidues"] = True

        def replaceNonstandardResidues(self):
            calls["replaceNonstandardResidues"] = True

        def findMissingAtoms(self):
            calls["findMissingAtoms"] = True

        def addMissingAtoms(self):
            calls["addMissingAtoms"] = True

    class FakeModeller:
        def __init__(self, topology, positions):
            self.topology = SimpleNamespace(
                atoms=lambda: (SimpleNamespace(element=SimpleNamespace(symbol="H")),)
            )
            self.positions = positions

        def delete(self, hydrogens):
            calls["deleted_hydrogens"] = len(hydrogens)

        def addHydrogens(self, force_field, pH):
            calls["force_field"] = force_field.name
            calls["ph"] = pH

    class FakeForceField:
        def __init__(self, name):
            self.name = name

    class FakePDBFile:
        @staticmethod
        def writeFile(topology, positions, handle, keepIds=True):
            calls["keepIds"] = keepIds
            handle.write(
                "ATOM      1  N   LYS A  23       0.000   0.000   0.000  1.00  0.00           N\n"
                "ATOM      2  H   LYS A  23       0.000   0.000   0.000  1.00  0.00           H\n"
                "END\n"
            )

    app_module = ModuleType("openmm.app")
    app_module.ForceField = FakeForceField
    app_module.Modeller = FakeModeller
    app_module.PDBFile = FakePDBFile
    pdbfixer_module = ModuleType("pdbfixer")
    pdbfixer_module.PDBFixer = FakePDBFixer
    monkeypatch.setitem(sys.modules, "openmm", ModuleType("openmm"))
    monkeypatch.setitem(sys.modules, "openmm.app", app_module)
    monkeypatch.setitem(sys.modules, "pdbfixer", pdbfixer_module)

    result = canonicalize_protein_hydrogens(
        source,
        output,
        settings=ProteinCanonicalizationSettings(ph=6.5, force_field_name="test.xml"),
    )

    assert result.output_path == output
    assert result.ph == 6.5
    assert result.force_field_name == "test.xml"
    assert result.hydrogen_count == 1
    assert calls["ph"] == 6.5
    assert calls["force_field"] == "test.xml"
    assert calls["deleted_hydrogens"] == 1
    assert calls["keepIds"] is True
    assert any("pH: 6.5" in diagnostic for diagnostic in result.diagnostics)
