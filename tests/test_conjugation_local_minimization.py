"""Tests for post-crosslink local minimization helpers."""

from __future__ import annotations

import numpy as np

from polyzymd.builders.conjugation._linkage import PabloCrosslinkRequirement
from polyzymd.builders.conjugation._relaxation import (
    LocalLinkageAtomSelector,
    LocalLinkageSelectors,
    LocalMinimizationSettings,
    _build_simulation_with_platform_fallback,
    analyze_crosslink_geometry,
    build_product_state_pablo_policy,
    product_state_pablo_crosslink_requirement,
    write_pdb_with_replaced_coordinates,
)

POC_CROSSLINKED_PDB = (
    "src/polyzymd/builders/conjugation/poc/output/walkthrough-generated-product/"
    "seeded_random_10mer_crosslinked.pdb"
)


def test_analyze_crosslink_geometry_reports_known_poc_clash():
    """The seeded product should report generic product-link metrics."""
    settings = _poc_linkage_settings()

    metrics = analyze_crosslink_geometry(POC_CROSSLINKED_PDB, settings=settings)

    assert len(metrics.linkages) == 1
    linkage = metrics.linkages[0]
    assert linkage.reciprocal_product_link_conect is True
    assert 1.25 <= linkage.protein_modifier_distance_angstrom <= 1.65
    assert linkage.passes is True
    assert metrics.passes is True


def test_default_product_state_pablo_policy_does_not_emit_hz2_hz3():
    """Default local minimization Pablo policy must not use reactant leaving atoms."""
    policy = build_product_state_pablo_policy(
        POC_CROSSLINKED_PDB,
        settings=_poc_linkage_settings(),
    )

    assert policy is not None
    dumped = policy.model_dump(mode="json")
    assert "HZ2" not in str(dumped)
    assert "HZ3" not in str(dumped)
    assert policy.crosslinks[0].leaving_atoms == ((), ())


def _poc_linkage_settings() -> LocalMinimizationSettings:
    """Return generic linkage settings for the bundled POC product."""
    return LocalMinimizationSettings(
        linkages=(
            LocalLinkageSelectors(
                protein_link_selector=LocalLinkageAtomSelector(
                    serial=None,
                    chain_id="A",
                    residue_name="LYX",
                    residue_number=23,
                    atom_name="NZ",
                ),
                modifier_link_selector=LocalLinkageAtomSelector(
                    serial=None,
                    chain_id="C",
                    residue_name="NHX",
                    residue_number=5,
                    atom_name="C047",
                ),
                target_bond_length_angstrom=1.33,
                label="poc_linkage",
            ),
        )
    )


def test_resolved_product_state_crosslink_uses_empty_pablo_leaving_atoms():
    """Resolved reactant leaving atoms are cleared for already-modified product PDBs."""
    reactant_state_requirement = PabloCrosslinkRequirement(
        residues=("LYX", "NHX"),
        linking_atoms=("NZ", "C047"),
        leaving_atoms=(("H11", "H13"), ("O021",)),
        bond_order=1,
    )

    requirement = product_state_pablo_crosslink_requirement(
        POC_CROSSLINKED_PDB,
        pablo_crosslink_requirement=reactant_state_requirement,
    )

    assert requirement.residues == ("LYX", "NHX")
    assert requirement.linking_atoms == ("NZ", "C047")
    assert requirement.leaving_atoms == ((), ())


def test_product_state_crosslink_can_derive_from_resolved_plan():
    """Product-state requirements should accept resolved-plan shaped objects."""

    class ResolvedPlanLike:
        pablo_crosslink_requirement = PabloCrosslinkRequirement(
            residues=("LYX", "NHX"),
            linking_atoms=("NZ", "C047"),
            leaving_atoms=(("H11", "H13"), ("O021",)),
            bond_order=1,
        )

    requirement = product_state_pablo_crosslink_requirement(
        POC_CROSSLINKED_PDB,
        resolved_plan=ResolvedPlanLike(),
    )
    policy = build_product_state_pablo_policy(
        POC_CROSSLINKED_PDB,
        resolved_plan=ResolvedPlanLike(),
    )

    assert requirement.residues == ("LYX", "NHX")
    assert requirement.linking_atoms == ("NZ", "C047")
    assert requirement.leaving_atoms == ((), ())
    assert policy is not None
    assert policy.crosslinks[0].leaving_atoms == ((), ())


def test_write_pdb_with_replaced_coordinates_preserves_noncoordinate_text(tmp_path):
    """Coordinate replacement should preserve atom metadata and CONECT lines."""
    template = tmp_path / "template.pdb"
    template.write_text(
        "ATOM      1  NZ  LYX A  23       1.000   2.000   3.000  1.00  0.00           N1+\n"
        "HETATM    2 C047 NHX C   5       4.000   5.000   6.000  1.00  0.00           C  \n"
        "CONECT    1    2\n"
        "END\n",
        encoding="utf-8",
    )
    output = tmp_path / "relaxed.pdb"

    write_pdb_with_replaced_coordinates(
        template,
        np.asarray([[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]]),
        output,
    )

    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[0][:30] == "ATOM      1  NZ  LYX A  23    "
    assert lines[0][30:54] == "  10.000  20.000  30.000"
    assert lines[0][54:] == "  1.00  0.00           N1+"
    assert lines[1][:30] == "HETATM    2 C047 NHX C   5    "
    assert lines[1][30:54] == "  40.000  50.000  60.000"
    assert lines[2] == "CONECT    1    2"
    assert lines[3] == "END"


def test_local_minimization_platform_creation_falls_back_to_cpu():
    """Context creation failures on CUDA should retry CPU when no platform is requested."""

    class FakePlatform:
        def __init__(self, name: str):
            self._name = name

        def getName(self):
            return self._name

    class FakeOpenMM:
        class Platform:
            @staticmethod
            def getPlatformByName(name: str):
                return FakePlatform(name)

        class VerletIntegrator:
            def __init__(self, timestep):
                self.timestep = timestep

    class FakeOpenMMUnit:
        picoseconds = 1.0

    class FakeOpenMMApp:
        class Simulation:
            def __init__(self, topology, system, integrator, platform):
                if platform.getName() in {"CUDA", "OpenCL"}:
                    raise RuntimeError("CUDA_ERROR_UNSUPPORTED_PTX_VERSION")
                self.platform = platform

    simulation, platform_name = _build_simulation_with_platform_fallback(
        FakeOpenMM,
        FakeOpenMMApp,
        object(),
        object(),
        LocalMinimizationSettings(),
        FakeOpenMMUnit,
    )

    assert platform_name == "CPU"
    assert simulation.platform.getName() == "CPU"
