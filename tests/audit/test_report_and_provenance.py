"""Accepted exceptions, the build record and the generated check catalog."""

from __future__ import annotations

from pathlib import Path

import pytest

from polyzymd.audit.findings import Finding

DOCS = Path(__file__).resolve().parents[2] / "docs" / "source" / "reference" / "physics_checks.md"


def test_the_check_catalog_page_is_generated_from_the_registry() -> None:
    """Regenerate with: python -c "from polyzymd.audit.registry import catalog_markdown; ..." """
    from polyzymd.audit.registry import catalog, catalog_markdown

    text = catalog_markdown()
    assert DOCS.read_text() == text, "docs/source/reference/physics_checks.md is out of date"
    for spec in catalog():
        assert f"`{spec.id}`" in text


class _Exception:
    def __init__(self, check, component, reason):
        self.check, self.component, self.reason = check, component, reason


class _Config:
    def __init__(self, exceptions):
        self.modeling_exceptions = exceptions


def _audit_with(findings, exceptions, monkeypatch):
    from polyzymd.audit import report

    monkeypatch.setattr(report, "run_checks", lambda *args: [Finding(*f) for f in findings])

    class _Realized:
        n_atoms = 10

    return report.audit(_Realized(), [], _Config(exceptions))


def test_an_accepted_failure_passes_and_is_recorded_with_its_reason(monkeypatch) -> None:
    findings = [
        ("geometry.buried_water", "system", "WARN", "2 waters buried"),
        ("parameters.bonded", "protein", "FAIL", "Sage terms on the protein"),
    ]
    accepted = _audit_with(
        findings,
        [_Exception("parameters.bonded", "protein", "custom protein force field, see methods")],
        monkeypatch,
    )
    assert accepted.passed
    assert [a.check for a in accepted.accepted] == ["parameters.bonded"]
    assert accepted.accepted[0].reason.startswith("custom protein")
    assert {f.status for f in accepted.findings} == {"WARN", "ACCEPTED"}
    bonded = next(f for f in accepted.findings if f.check == "parameters.bonded")
    assert (bonded.accepted_from, bonded.reason) == (
        "FAIL",
        "custom protein force field, see methods",
    )
    assert bonded.to_dict()["accepted_from"] == "FAIL"


def test_an_exception_for_another_component_does_not_accept(monkeypatch) -> None:
    findings = [("parameters.bonded", "protein", "FAIL", "Sage terms on the protein")]
    report = _audit_with(findings, [_Exception("parameters.bonded", "substrate", "x")], monkeypatch)
    assert not report.passed
    assert any(f.check == "modeling_exceptions.unused" for f in report.findings)


def test_a_stale_exception_is_reported(monkeypatch) -> None:
    report = _audit_with(
        [], [_Exception("geometry.buried_water", None, "crystal waters")], monkeypatch
    )
    assert report.passed
    assert [(f.check, f.status) for f in report.findings] == [
        ("modeling_exceptions.unused", "WARN")
    ]


def test_an_unknown_check_in_the_config_names_the_closest() -> None:
    from pydantic import ValidationError

    from polyzymd.config.schema import ModelingException

    with pytest.raises(ValidationError, match="did you mean geometry.buried_water"):
        ModelingException(check="geometry.buried_waters", reason="x")
    with pytest.raises(ValidationError, match="reason"):
        ModelingException(check="geometry.buried_water", reason="  ")


def test_an_audit_that_cannot_run_is_a_failure_not_a_pass() -> None:
    from polyzymd.simulation.artifact_integrity import _audit

    report = _audit(object(), ["a component"], None)
    assert not report.passed
    assert "could not run" in report.findings[0].detail


def test_an_exception_that_could_never_match_is_refused_by_the_config() -> None:
    """A component a check never reports on is a config error, not a silent no-op."""
    import yaml
    from pydantic import ValidationError

    from polyzymd.config.schema import SimulationConfig

    quickstart = Path(__file__).resolve().parents[2] / "examples" / "quickstart" / "config.yaml"
    data = yaml.safe_load(quickstart.read_text())
    data["solutes"][0]["pdb_path"] = str(quickstart.parent / "trpcage.pdb")
    data["modeling_exceptions"] = [
        {"check": "composition.polymer_sequences", "component": "trpcage", "reason": "x"}
    ]
    with pytest.raises(ValidationError, match="cannot report on 'trpcage'"):
        SimulationConfig.model_validate(data)
    data["modeling_exceptions"] = [
        {"check": "parameters.bonded", "component": "trpcage", "reason": "custom FF, see methods"}
    ]
    assert SimulationConfig.model_validate(data).modeling_exceptions[0].component == "trpcage"


def test_an_unreadable_system_is_an_audit_failure() -> None:
    from polyzymd.simulation.artifact_integrity import _audit

    def unreadable():
        raise ValueError("system.xml is truncated")

    report = _audit(unreadable, ["a component"], None)
    assert not report.passed
    assert "truncated" in report.findings[0].detail


def test_the_build_record_lists_the_methods_to_cite() -> None:
    """Packmol when it packed the build; Machado and Pantano when ionic_strength set the ions."""
    from types import SimpleNamespace

    from polyzymd.simulation.artifact_integrity import methods_to_cite, provenance_markdown

    def config(strength):
        return SimpleNamespace(
            solvent=SimpleNamespace(ions=SimpleNamespace(ionic_strength=strength))
        )

    def others(cited):
        return [c for c in cited if c.startswith(("Martínez", "Machado"))]

    assert others(methods_to_cite(config(None))) == []
    packed = {"packmol_seed": 1, "solvent_packmol_seed": 1}
    cited = others(methods_to_cite(config(None), packed))
    assert len(cited) == 1 and "doi:10.1002/jcc.21224" in cited[0]
    # packmol_seed alone is recorded before anything is packed; it is not evidence.
    assert others(methods_to_cite(config(None), {"packmol_seed": 1})) == []
    cited = others(methods_to_cite(config(0.13), {"polymer_packmol_seed": 2}))
    assert [c.split(" ")[0] for c in cited] == ["Martínez,", "Machado,"]
    assert "doi:10.1021/acs.jctc.9b00953" in cited[1]
    text = provenance_markdown({"components": [], "cite": cited})
    assert "## Methods to cite" in text and "PACKMOL" in text and "Machado" in text
    assert "explanation/physics_audit.html" in text


def _build_config(*, ligand_charges=None, polymer=None, force_field="openff-2.3.0.offxml"):
    """A config-shaped namespace with a protein, an optional ligand and polymer, and water."""
    from types import SimpleNamespace

    def solute(fields, method, model=None, **keys):
        return SimpleNamespace(
            force_field=fields, charges=SimpleNamespace(method=method, model=model), **keys
        )

    solutes = [solute(["ff14sb_off_impropers_0.0.4.offxml"], "library")]
    if ligand_charges:
        solutes.append(solute([force_field], *ligand_charges))
    if polymer is not None:
        solutes.append(solute(["openff-1.3.0.offxml"], "preset", enabled=polymer))
    return SimpleNamespace(
        solutes=solutes,
        solvent=SimpleNamespace(
            force_field=[force_field],
            water=SimpleNamespace(model=SimpleNamespace(value="tip3p")),
            ions=SimpleNamespace(ionic_strength=None),
            co_solvents=[],
        ),
    )


def test_methods_to_cite_name_polyzymd_the_force_fields_water_and_charge_models() -> None:
    from polyzymd import __version__
    from polyzymd.simulation.artifact_integrity import methods_to_cite

    cited = methods_to_cite(
        _build_config(ligand_charges=("nagl", "openff-gnn-am1bcc-1.0.0.pt")),
        {"solvent_packmol_seed": 1},
    )
    assert cited[0].startswith("Laforet, Joseph R., Jr. PolyzyMD") and __version__ in cited[0]
    assert cited[1].startswith("OpenFF Toolkit") and OPENFF_TOOLKIT in cited[1]
    assert cited[2].startswith("OpenFF Interchange") and OPENFF_INTERCHANGE in cited[2]
    assert cited[3].startswith("ff14sb_off_impropers_0.0.4.offxml (the OpenFF SMIRNOFF port")
    assert "doi:10.1021/acs.jctc.5b00255" in cited[3]
    assert cited[4].startswith("openff-2.3.0.offxml (a release of OpenFF Sage 2)")
    assert "doi:10.1021/acs.jctc.3c00039" in cited[4]
    assert cited[5].startswith("openff-2.3.0.offxml (OpenFF Sage 2.3.0)")
    assert "doi:10.1021/acs.jctc.6c00169" in cited[5]
    assert cited[6].startswith("TIP3P water") and "doi:10.1063/1.445869" in cited[6]
    assert cited[7].startswith("AshGC (OpenFF NAGL model openff-gnn-am1bcc-1.0.0.pt)")
    assert "doi:10.1021/acs.jctc.6c00169" in cited[7]
    assert "doi:10.1002/(SICI)1096-987X(20000130)21:2<132::AID-JCC5>3.0.CO;2-P" in cited[8]
    assert "doi:10.1002/jcc.10128" in cited[9]
    assert "doi:10.1002/jcc.21224" in cited[10] and len(cited) == 11


OPENFF_TOOLKIT = "https://github.com/openforcefield/openff-toolkit"
OPENFF_INTERCHANGE = "https://github.com/openforcefield/openff-interchange"


def test_methods_to_cite_name_the_engine_and_the_openff_versions() -> None:
    from polyzymd.simulation.artifact_integrity import methods_to_cite

    record = {
        "engine": "openmm",
        "openmm_version": "8.2.0",
        "openff_toolkit_version": "0.16.9",
        "openff_interchange_version": "0.4.2",
    }
    openmm = methods_to_cite(_build_config(), record=record)
    assert openmm[1].startswith("OpenMM 8.2.0: Eastman, P.")
    assert "doi:10.1021/acs.jpcb.3c06662" in openmm[1]
    assert not any("doi:10.1016/j.softx.2015.06.001" in c for c in openmm)
    assert openmm[2] == (
        f"OpenFF Toolkit 0.16.9: no journal article; cite the release ({OPENFF_TOOLKIT})."
    )
    assert openmm[3].startswith("OpenFF Interchange 0.4.2:") and OPENFF_INTERCHANGE in openmm[3]
    gromacs = methods_to_cite(_build_config(), record={**record, "engine": "gromacs"})
    assert gromacs[1].startswith("GROMACS: Abraham, M. J.")
    assert "doi:10.1016/j.softx.2015.06.001" in gromacs[1]
    assert not any("doi:10.1021/acs.jpcb.3c06662" in c for c in gromacs)
    # Without a record (no engine known) neither engine is cited.
    assert not any(c.startswith(("OpenMM", "GROMACS")) for c in methods_to_cite(_build_config()))


def test_sage_2_3_and_ashgc_cite_their_paper_and_older_files_do_not() -> None:
    from polyzymd.simulation.artifact_integrity import methods_to_cite

    def cites_2026(cited):
        return [c for c in cited if "doi:10.1021/acs.jctc.6c00169" in c]

    for name in ("openff-2.3.0.offxml", "openff_unconstrained-2.3.0-rc2.offxml"):
        assert len(cites_2026(methods_to_cite(_build_config(force_field=name)))) == 1
    sage_2_2 = methods_to_cite(_build_config(force_field="openff-2.2.1.offxml"))
    assert cites_2026(sage_2_2) == []
    assert any("doi:10.1021/acs.jctc.3c00039" in c for c in sage_2_2)
    ashgc = cites_2026(
        methods_to_cite(
            _build_config(
                ligand_charges=("nagl", "openff-gnn-am1bcc-1.0.0.pt"),
                force_field="openff-2.2.1.offxml",
            )
        )
    )
    assert len(ashgc) == 1 and ashgc[0].startswith("AshGC (OpenFF NAGL model ")
    older = methods_to_cite(
        _build_config(ligand_charges=("nagl", "openff-gnn-am1bcc-0.1.0-rc.3.pt"))
    )
    assert not any(c.startswith("AshGC") for c in older)
    assert any(
        c.startswith("OpenFF NAGL charge model openff-gnn-am1bcc-0.1.0-rc.3.pt") for c in older
    )


def _ions_record(*force_field):
    return {
        "components": [
            {"type": "ions", "domains": [{"name": "Na+", "force_field": list(force_field)}]}
        ]
    }


def test_ion_parameters_cite_joung_and_cheatham_when_the_build_has_ions() -> None:
    from polyzymd.simulation.artifact_integrity import methods_to_cite

    def ions(cited):
        return [c for c in cited if "doi:10.1021/jp8001614" in c]

    cited = ions(methods_to_cite(_build_config(), record=_ions_record("openff-2.3.0.offxml")))
    assert cited == [
        "Na+ and Cl- parameters of openff-2.3.0.offxml (TIP3P set): Joung, I. S. and "
        "Cheatham, T. E. (2008). Determination of alkali and halide monovalent ion parameters "
        "for use in explicitly solvated biomolecular simulations. J. Phys. Chem. B 112(30), "
        "9020-9041. doi:10.1021/jp8001614"
    ]
    assert ions(methods_to_cite(_build_config())) == []
    assert ions(methods_to_cite(_build_config(), record=_ions_record("/data/my_ff.offxml"))) == []


def test_openff_ion_parameters_are_the_joung_cheatham_tip3p_set() -> None:
    """The values Joung and Cheatham (2008) give in Table 5 for TIP3P."""
    import xml.etree.ElementTree as ElementTree

    openforcefields = pytest.importorskip("openforcefields")
    offxml = Path(openforcefields.__file__).parent / "offxml"
    table_5 = {"[#11+1:1]": ("1.369", "0.0874393"), "[#17X0-1:1]": ("2.513", "0.035591")}
    for name in ("openff-1.3.0.offxml", "openff-2.0.0.offxml", "openff-2.3.0.offxml"):
        atoms = ElementTree.parse(offxml / name).getroot().find("vdW")
        found = {
            atom.get("smirks"): (
                atom.get("rmin_half").split(" ")[0],
                atom.get("epsilon").split(" ")[0],
            )
            for atom in atoms
            if atom.get("smirks") in table_5
        }
        assert found == table_5, name


def test_methods_to_cite_follow_what_the_build_used() -> None:
    from polyzymd.simulation.artifact_integrity import methods_to_cite

    apo = methods_to_cite(_build_config())
    assert not any("AM1-BCC" in c or "NAGL" in c or "Espaloma" in c for c in apo)
    # One file named by several components is cited once.
    assert sum(c.startswith("openff-2.3.0.offxml (a release of OpenFF Sage 2)") for c in apo) == 1
    espaloma = methods_to_cite(_build_config(ligand_charges=("espaloma",)))
    assert any("doi:10.1021/acs.jpca.4c01287" in c for c in espaloma)
    assert not any("AM1-BCC model" in c for c in espaloma)
    parsley = methods_to_cite(_build_config(force_field="openff_unconstrained-1.3.0.offxml"))
    assert any("doi:10.1021/acs.jctc.1c00571" in c for c in parsley)
    custom = methods_to_cite(_build_config(force_field="/data/my_ff.offxml"))
    assert "my_ff.offxml: PolyzyMD has no reference for this file; cite its source." in custom
    # A disabled polymer entry is not built, so its force field is not cited.
    assert methods_to_cite(_build_config(polymer=False)) == apo
    assert any(
        "doi:10.1021/acs.jctc.1c00571" in c for c in methods_to_cite(_build_config(polymer=True))
    )
