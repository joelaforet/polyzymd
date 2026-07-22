"""Tests for OpenMM handoff artifact helpers."""

from __future__ import annotations

import json
import sys
import types
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import pytest

from polyzymd.builders.openmm_artifacts import (
    SOLVATED_CONJUGATE_PDB,
    SOLVATED_SYSTEM_PDB,
    _reject_disconnected_hetatm,
    _write_solute_audit,
    build_openmm_artifacts,
    ensure_conjugation_pdb_alias,
    resolve_skip_build_artifacts,
)
from polyzymd.config.schema import BuildScope


def _config(*, conjugation_enabled: bool) -> SimpleNamespace:
    """Create a minimal config-like object for helper tests.

    Parameters
    ----------
    conjugation_enabled : bool
        Whether the fake config should enable conjugation routing.

    Returns
    -------
    SimpleNamespace
        Minimal config object consumed by helper functions.
    """
    return SimpleNamespace(conjugation=SimpleNamespace(enabled=conjugation_enabled), restraints=[])


def test_ensure_conjugation_pdb_alias_copies_and_refreshes(tmp_path: Path) -> None:
    """The standard PDB alias should be copied and refreshed from conjugation output."""
    source = tmp_path / SOLVATED_CONJUGATE_PDB
    alias = tmp_path / SOLVATED_SYSTEM_PDB
    source.write_text("fresh conjugation pdb\n", encoding="utf-8")
    alias.write_text("stale pdb\n", encoding="utf-8")

    resolved = ensure_conjugation_pdb_alias(tmp_path)

    assert resolved == alias
    assert alias.read_text(encoding="utf-8") == source.read_text(encoding="utf-8")


def test_resolve_skip_build_conjugation_falls_back_to_route_specific_pdb(
    tmp_path: Path,
) -> None:
    """Conjugation skip-build should copy the route-specific PDB when alias is absent."""
    system_path = tmp_path / "system.xml"
    source = tmp_path / SOLVATED_CONJUGATE_PDB
    system_path.write_text("<System />", encoding="utf-8")
    source.write_text("conjugation pdb\n", encoding="utf-8")

    pdb_path, resolved_system_path = resolve_skip_build_artifacts(
        _config(conjugation_enabled=True),
        tmp_path,
    )

    assert pdb_path == tmp_path / SOLVATED_SYSTEM_PDB
    assert resolved_system_path == system_path
    assert pdb_path.read_text(encoding="utf-8") == "conjugation pdb\n"


def test_resolve_skip_build_conjugation_refreshes_stale_existing_alias(
    tmp_path: Path,
) -> None:
    """Conjugation skip-build should refresh stale standard PDB aliases."""
    system_path = tmp_path / "system.xml"
    source = tmp_path / SOLVATED_CONJUGATE_PDB
    alias = tmp_path / SOLVATED_SYSTEM_PDB
    system_path.write_text("<System />", encoding="utf-8")
    source.write_text("fresh conjugation pdb\n", encoding="utf-8")
    alias.write_text("stale non-conjugation pdb\n", encoding="utf-8")

    pdb_path, resolved_system_path = resolve_skip_build_artifacts(
        _config(conjugation_enabled=True),
        tmp_path,
    )

    assert pdb_path == alias
    assert resolved_system_path == system_path
    assert alias.read_text(encoding="utf-8") == "fresh conjugation pdb\n"


def test_resolve_skip_build_always_requires_system_xml(tmp_path: Path) -> None:
    """Skip-build validation should fail actionably when system.xml is missing."""
    (tmp_path / SOLVATED_SYSTEM_PDB).write_text("pdb\n", encoding="utf-8")

    with pytest.raises(FileNotFoundError, match="system.xml"):
        resolve_skip_build_artifacts(_config(conjugation_enabled=True), tmp_path)


def test_build_openmm_artifacts_routes_conjugation_without_heavy_write(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Conjugation export-only builds should refresh the standard PDB alias."""
    captured: dict[str, object] = {}
    source = tmp_path / SOLVATED_CONJUGATE_PDB
    alias = tmp_path / SOLVATED_SYSTEM_PDB
    alias.write_text("stale pdb\n", encoding="utf-8")
    result = SimpleNamespace(
        system_builder=SimpleNamespace(get_component_info=lambda: {"builder": True}),
        require_final_interchange=lambda: "final-interchange",
        get_component_info=lambda: {"conjugate": True},
    )

    def fake_build_conjugate_from_config(*args: object, **kwargs: object) -> object:
        """Record conjugation facade arguments for assertions.

        Parameters
        ----------
        *args : object
            Positional arguments forwarded by the helper.
        **kwargs : object
            Keyword arguments forwarded by the helper.

        Returns
        -------
        object
            Fake conjugation workflow result.
        """
        captured["args"] = args
        captured["kwargs"] = kwargs
        source.write_text("fresh conjugation pdb\n", encoding="utf-8")
        return result

    class FakeSettings:
        """Fake conjugation settings class with the required flag."""

        def __init__(
            self,
            *,
            create_final_interchange: bool,
            pdb_fragment_output_mode: str = "coordinate_only",
        ) -> None:
            self.create_final_interchange = create_final_interchange
            self.pdb_fragment_output_mode = pdb_fragment_output_mode

        def model_copy(self, *, update: dict[str, object]) -> "FakeSettings":
            """Return a fake settings copy with selected updates."""
            values = {
                "create_final_interchange": self.create_final_interchange,
                "pdb_fragment_output_mode": self.pdb_fragment_output_mode,
            }
            values.update(update)
            return FakeSettings(**values)

    conjugation_module = types.ModuleType("polyzymd.builders.conjugation")
    conjugation_module.build_conjugate_from_config = fake_build_conjugate_from_config
    workflow_module = types.ModuleType("polyzymd.builders.conjugation.system_workflow")
    workflow_module.ConjugatedPolymerSystemSettings = FakeSettings
    native_module = types.ModuleType("polyzymd.builders.conjugation.native_openmm_glycam")
    native_module.native_glycam_enabled = lambda _config: False
    monkeypatch.setitem(sys.modules, "polyzymd.builders.conjugation", conjugation_module)
    monkeypatch.setitem(
        sys.modules, "polyzymd.builders.conjugation.system_workflow", workflow_module
    )
    monkeypatch.setitem(
        sys.modules, "polyzymd.builders.conjugation.native_openmm_glycam", native_module
    )

    artifacts = build_openmm_artifacts(
        sim_config=_config(conjugation_enabled=True),
        working_dir=tmp_path,
        polymer_seed=7,
        write_system=False,
    )

    kwargs = cast(dict[str, Any], captured["kwargs"])
    assert kwargs["output_dir"] == tmp_path
    assert kwargs["free_polymer_seed"] == 7
    assert kwargs["settings"].create_final_interchange is True
    assert artifacts.require_final_interchange() == "final-interchange"
    assert artifacts.get_component_info() == {"conjugate": True}
    assert artifacts.pdb_path == alias
    assert alias.read_text(encoding="utf-8") == "fresh conjugation pdb\n"
    assert not artifacts.system_path.exists()


def test_build_openmm_artifacts_uses_native_glycam_product_state_route(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Native GLYCAM CLI builds should not remain on coordinate-only defaults."""
    captured: dict[str, object] = {}
    source = tmp_path / SOLVATED_CONJUGATE_PDB

    result = SimpleNamespace(
        system_builder=None,
        exact_export_bundle="exact-bundle",
        require_final_interchange=lambda: "exact-bundle",
        get_component_info=lambda: {},
    )

    def fake_build_conjugate_from_config(*args: object, **kwargs: object) -> object:
        """Record native workflow settings and create the route-specific PDB."""
        captured["kwargs"] = kwargs
        source.write_text("native conjugation pdb\n", encoding="utf-8")
        return result

    class FakeSettings:
        """Fake settings class with Pydantic-like copy support."""

        def __init__(
            self,
            *,
            create_final_interchange: bool,
            pdb_fragment_output_mode: str = "coordinate_only",
        ) -> None:
            self.create_final_interchange = create_final_interchange
            self.pdb_fragment_output_mode = pdb_fragment_output_mode

        def model_copy(self, *, update: dict[str, object]) -> "FakeSettings":
            """Return a fake settings copy with selected updates."""
            data = {
                "create_final_interchange": self.create_final_interchange,
                "pdb_fragment_output_mode": self.pdb_fragment_output_mode,
            }
            data.update(update)
            return FakeSettings(**data)

    conjugation_module = types.ModuleType("polyzymd.builders.conjugation")
    conjugation_module.build_conjugate_from_config = fake_build_conjugate_from_config
    workflow_module = types.ModuleType("polyzymd.builders.conjugation.system_workflow")
    workflow_module.ConjugatedPolymerSystemSettings = FakeSettings
    native_module = types.ModuleType("polyzymd.builders.conjugation.native_openmm_glycam")
    native_module.native_glycam_enabled = lambda _config: True
    monkeypatch.setitem(sys.modules, "polyzymd.builders.conjugation", conjugation_module)
    monkeypatch.setitem(
        sys.modules, "polyzymd.builders.conjugation.system_workflow", workflow_module
    )
    monkeypatch.setitem(
        sys.modules, "polyzymd.builders.conjugation.native_openmm_glycam", native_module
    )

    artifacts = build_openmm_artifacts(
        sim_config=_config(conjugation_enabled=True),
        working_dir=tmp_path,
        polymer_seed=11,
        write_system=False,
    )

    settings = cast(dict[str, Any], captured["kwargs"])["settings"]
    assert settings.pdb_fragment_output_mode == "experimental_pablo"
    assert artifacts.require_final_interchange() == "exact-bundle"
    assert artifacts.exact_export_bundle == "exact-bundle"


def test_conjugation_workflow_settings_inherit_openmm_platform() -> None:
    """Config-driven conjugation relaxation should follow the configured OpenMM platform."""
    from polyzymd.builders.conjugation.system_workflow import (
        ConjugatedPolymerSystemSettings,
        _settings_with_config_defaults,
    )

    settings = ConjugatedPolymerSystemSettings()
    config = SimpleNamespace(openmm=SimpleNamespace(platform="CPU"))

    resolved = _settings_with_config_defaults(settings, config)

    assert resolved.relaxation.platform_name == "CPU"


def test_build_openmm_artifacts_routes_standard_system_builder(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Non-conjugation builds should preserve the standard SystemBuilder route."""
    config = _config(conjugation_enabled=False)
    builder = SimpleNamespace(
        build_from_config=lambda **kwargs: ("interchange", kwargs),
        get_component_info=lambda: {"standard": True},
    )
    system_builder_module = types.ModuleType("polyzymd.builders.system_builder")
    system_builder_module.SystemBuilder = SimpleNamespace(from_config=lambda arg: builder)
    monkeypatch.setitem(sys.modules, "polyzymd.builders.system_builder", system_builder_module)

    artifacts = build_openmm_artifacts(
        sim_config=config,
        working_dir=tmp_path,
        polymer_seed=3,
        write_system=False,
    )

    interchange, kwargs = artifacts.require_final_interchange()
    assert interchange == "interchange"
    assert kwargs == {"config": config, "working_dir": tmp_path, "polymer_seed": 3}
    assert artifacts.get_component_info() == {"standard": True}


def test_structure_scope_returns_only_assembled_conjugate(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Structure scope should expose the validated PDB without system handoff artifacts."""
    assembled = tmp_path / "conjugate-construction" / "assembled_crosslinked.pdb"
    assembled.parent.mkdir(parents=True)
    assembled.write_text("END\n", encoding="utf-8")
    captured: dict[str, object] = {}
    result = SimpleNamespace(
        system_builder=None,
        crosslinked_conjugate_pdb_path=assembled,
        workflow_json_path=tmp_path / "conjugated_polymer_system_workflow.json",
    )

    def fake_build(*args: object, **kwargs: object) -> object:
        captured["settings"] = kwargs["settings"]
        return result

    monkeypatch.setattr("polyzymd.builders.conjugation.build_conjugate_from_config", fake_build)

    artifacts = build_openmm_artifacts(
        sim_config=_config(conjugation_enabled=True),
        working_dir=tmp_path,
        polymer_seed=1,
        scope=BuildScope.STRUCTURE,
    )

    assert artifacts.pdb_path == assembled
    assert artifacts.system_path is None
    assert artifacts.builder is None
    assert artifacts.scope is BuildScope.STRUCTURE
    settings = cast(Any, captured["settings"])
    assert settings.build_scope is BuildScope.STRUCTURE
    with pytest.raises(RuntimeError, match="do not contain a final Interchange"):
        artifacts.require_final_interchange()
    assert not (tmp_path / SOLVATED_SYSTEM_PDB).exists()
    assert not (tmp_path / "system.xml").exists()


def test_nonconjugate_structure_scope_fails_before_system_builder(tmp_path: Path) -> None:
    """A bare protein cannot request the conjugation-only structure checkpoint."""
    with pytest.raises(ValueError, match="requires conjugation.enabled"):
        build_openmm_artifacts(
            sim_config=_config(conjugation_enabled=False),
            working_dir=tmp_path,
            polymer_seed=1,
            scope=BuildScope.STRUCTURE,
        )


def test_standard_solute_scope_emits_only_isolated_audited_artifacts(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Solute scope must use the primary-only branch and audit its serialized System."""
    from openmm import System

    calls: list[str] = []
    system = System()
    system.addParticle(12.0)
    system.addParticle(12.0)

    class FakeBuilder:
        interchange = "primary-interchange"

        def build_isolated_primary_from_config(self, config: object) -> object:
            calls.append("isolated")
            return self.interchange

        def build_from_config(self, **kwargs: object) -> object:
            raise AssertionError("full-system build entered for solute scope")

        def save_topology(self, path: Path) -> None:
            path.write_text(
                "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
                "HETATM    2  C1  LIG C   1       1.500   0.000   0.000  1.00  0.00           C\n"
                "CONECT    1    2\nEND\n",
                encoding="utf-8",
            )

        def get_openmm_components(self) -> tuple[None, System, None]:
            return None, system, None

    builder = FakeBuilder()
    module = types.ModuleType("polyzymd.builders.system_builder")
    module.SystemBuilder = SimpleNamespace(from_config=lambda config: builder)
    monkeypatch.setitem(sys.modules, "polyzymd.builders.system_builder", module)

    artifacts = build_openmm_artifacts(
        sim_config=_config(conjugation_enabled=False),
        working_dir=tmp_path,
        polymer_seed=1,
        scope=BuildScope.SOLUTE,
    )

    solute_dir = tmp_path / "solute"
    assert calls == ["isolated"]
    assert artifacts.pdb_path == solute_dir / "solute.pdb"
    assert artifacts.system_path == solute_dir / "system.xml"
    assert {path.name for path in solute_dir.iterdir()} == {
        "solute.pdb",
        "system.xml",
        "openmm_build_audit.json",
    }
    assert not (tmp_path / "system.xml").exists()
    audit = json.loads((solute_dir / "openmm_build_audit.json").read_text(encoding="utf-8"))
    assert audit["pdb_atom_count"] == audit["system_particle_count"] == 2
    assert audit["component_counts"] == {
        "substrate": 0,
        "free_polymers": 0,
        "water": 0,
        "ions": 0,
        "liquids": 0,
    }
    assert audit["barostat_count"] == audit["restraint_force_count"] == 0


def test_native_exact_solute_materializes_authoritative_bundle(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Native solute scope must persist the exact System, PDB, sidecar, and audits."""
    import numpy as np
    from openmm import NonbondedForce, System, Vec3, XmlSerializer, unit
    from openmm.app import Element, Topology

    from polyzymd.exporters.exact_openmm import ExactExportBundle, collect_exact_exception_sidecar

    topology = Topology()
    residue = topology.addResidue("ALA", topology.addChain("A"), "1")
    topology.addAtom("C1", Element.getBySymbol("C"), residue)
    topology.addAtom("C2", Element.getBySymbol("C"), residue)
    topology.addBond(*tuple(topology.atoms()))
    vectors = (Vec3(3, 0, 0), Vec3(0, 2.2, 0), Vec3(0, 0, 2.2)) * unit.nanometer
    topology.setPeriodicBoxVectors(vectors)
    positions = np.asarray([[1.1, 1.1, 1.1], [1.8, 1.1, 1.1]]) * unit.nanometer
    system = System()
    nonbonded = NonbondedForce()
    nonbonded.setNonbondedMethod(NonbondedForce.PME)
    nonbonded.setCutoffDistance(1.0 * unit.nanometer)
    for charge in (-0.2, 0.2):
        system.addParticle(12.0)
        nonbonded.addParticle(charge, 0.3, 0.1)
    system.addForce(nonbonded)
    system.setDefaultPeriodicBoxVectors(*vectors)
    audit = {
        "scope": "solute",
        "preparation_only_warning": "charged PME preparation/export-only; not NPT-ready",
    }
    sidecar = collect_exact_exception_sidecar(topology, system, audit=audit)
    bundle = ExactExportBundle(
        topology=topology,
        system=system,
        positions=positions,
        private_baseline_interchange=object(),
        sidecar=sidecar,
        audit=audit,
    )
    builder = SimpleNamespace(interchange=None)
    result = SimpleNamespace(
        system_builder=builder,
        exact_export_bundle=bundle,
        require_final_interchange=lambda: bundle,
        get_component_info=lambda: {},
    )

    class FakeSettings:
        def __init__(self, **values: object) -> None:
            self.__dict__.update(values)

        def model_copy(self, *, update: dict[str, object]) -> "FakeSettings":
            return FakeSettings(**{**self.__dict__, **update})

    conjugation_module = types.ModuleType("polyzymd.builders.conjugation")
    conjugation_module.build_conjugate_from_config = lambda *args, **kwargs: result
    workflow_module = types.ModuleType("polyzymd.builders.conjugation.system_workflow")
    workflow_module.ConjugatedPolymerSystemSettings = FakeSettings
    native_module = types.ModuleType("polyzymd.builders.conjugation.native_openmm_glycam")
    native_module.native_glycam_enabled = lambda _config: True
    monkeypatch.setitem(sys.modules, "polyzymd.builders.conjugation", conjugation_module)
    monkeypatch.setitem(
        sys.modules, "polyzymd.builders.conjugation.system_workflow", workflow_module
    )
    monkeypatch.setitem(
        sys.modules, "polyzymd.builders.conjugation.native_openmm_glycam", native_module
    )

    artifacts = build_openmm_artifacts(
        sim_config=_config(conjugation_enabled=True),
        working_dir=tmp_path,
        polymer_seed=1,
        scope=BuildScope.SOLUTE,
    )

    solute_dir = tmp_path / "solute"
    assert artifacts.exact_export_bundle is bundle
    assert {path.name for path in solute_dir.iterdir()} == {
        "solute.pdb",
        "system.xml",
        "exact_openmm_exceptions.json",
        "native_openmm_glycam_audit.json",
        "openmm_build_audit.json",
    }
    serialized = XmlSerializer.deserialize((solute_dir / "system.xml").read_text())
    assert serialized.getNumParticles() == topology.getNumAtoms() == sidecar.particle_count
    build_audit = json.loads((solute_dir / "openmm_build_audit.json").read_text())
    assert build_audit["route"] == "native_openmm_glycam"
    assert build_audit["exact_hashes"]["atom_order"] == sidecar.atom_order_hash
    assert build_audit["barostat_count"] == build_audit["restraint_force_count"] == 0


def test_assembled_solute_rejects_disconnected_lig_hetatm(tmp_path: Path) -> None:
    """An unowned assembled LIG must fail even when its residue name is not allowlisted."""
    pdb_path = tmp_path / "assembled.pdb"
    pdb_path.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  C1  LIG C   1       4.000   0.000   0.000  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="disconnected/unowned HETATM"):
        _reject_disconnected_hetatm(pdb_path)

    pdb_path.write_text(pdb_path.read_text(encoding="utf-8") + "CONECT    1    2\n")
    _reject_disconnected_hetatm(pdb_path)


def test_solute_audit_rejects_real_distance_restraint(tmp_path: Path) -> None:
    """Audit restraint counts must detect a force created through the real apply path."""
    from openmm import System, XmlSerializer
    from openmm.app import Element, Topology

    from polyzymd.core.restraints import AtomSelection, RestraintDefinition, RestraintType

    pdb_path = tmp_path / "solute.pdb"
    pdb_path.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  CB  ALA A   1       1.500   0.000   0.000  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    system = System()
    system.addParticle(12.0)
    system.addParticle(12.0)
    topology = Topology()
    residue = topology.addResidue("ALA", topology.addChain("A"), "1")
    topology.addAtom("CA", Element.getBySymbol("C"), residue)
    topology.addAtom("CB", Element.getBySymbol("C"), residue)
    restraint = RestraintDefinition(
        restraint_type=RestraintType.HARMONIC,
        name="audit-test",
        atom1=AtomSelection("index 0"),
        atom2=AtomSelection("index 1"),
    )
    restraint.apply(topology, system)
    assert system.getForce(0).getName() == "PolyzyMD restraint: audit-test"
    system_path = tmp_path / "system.xml"
    system_path.write_text(XmlSerializer.serialize(system), encoding="utf-8")

    with pytest.raises(RuntimeError, match="restraint-like force"):
        _write_solute_audit(
            builder=SimpleNamespace(),
            solute_dir=tmp_path,
            pdb_path=pdb_path,
            system_path=system_path,
        )
