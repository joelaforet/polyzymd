"""Focused tests for exact OpenMM exception export helpers."""

from __future__ import annotations

import json

import pytest

from polyzymd.exporters.exact_openmm import (
    AtomIdentity,
    ExactConstraintRecord,
    ExactExceptionRecord,
    ExactExceptionSidecar,
    ExactExportBundle,
    ExactExportError,
    ExactNonbondedMetadata,
    ExactParticleRecord,
    ExactTopologyBondRecord,
    authoritative_atom_order_hash,
    authoritative_topology_hash,
    exception_hash,
    gromacs_atom_order_hash,
    gromacs_topology_hash,
    particle_hash,
)
from polyzymd.exporters.gromacs import GromacsExporter, patch_gromacs_topology_with_exact_exceptions


def _sidecar() -> ExactExceptionSidecar:
    """Return a minimal exact exception sidecar for topology patch tests."""
    atom1 = AtomIdentity(index=1, name="C1", residue_name="MOL", residue_id="1", chain_id="C")
    atom2 = AtomIdentity(index=2, name="O1", residue_name="MOL", residue_id="1", chain_id="C")
    atom3 = AtomIdentity(index=3, name="H1", residue_name="MOL", residue_id="1", chain_id="C")
    particles = (
        ExactParticleRecord(index=1, charge_e=0.2, sigma_nm=0.3, epsilon_kj_mol=0.4),
        ExactParticleRecord(index=2, charge_e=-0.5, sigma_nm=0.31, epsilon_kj_mol=0.5),
        ExactParticleRecord(index=3, charge_e=0.1, sigma_nm=0.2, epsilon_kj_mol=0.0),
    )
    exceptions = (
        ExactExceptionRecord(
            exception_index=0,
            i=1,
            j=2,
            charge_product_e2=-0.1,
            sigma_nm=0.305,
            epsilon_kj_mol=0.4472135955,
            category="glycam_unscaled_14",
            atom_i=atom1,
            atom_j=atom2,
        ),
        ExactExceptionRecord(
            exception_index=1,
            i=2,
            j=3,
            charge_product_e2=0.0,
            sigma_nm=0.0,
            epsilon_kj_mol=0.0,
            category="zero_exclusion",
            atom_i=atom2,
            atom_j=atom3,
        ),
    )
    atoms = (atom1, atom2, atom3)
    topology_bonds = (
        ExactTopologyBondRecord(i=1, j=2, atom_i=atom1, atom_j=atom2),
        ExactTopologyBondRecord(i=2, j=3, atom_i=atom2, atom_j=atom3),
    )
    return ExactExceptionSidecar(
        particle_count=3,
        exception_count=2,
        nonzero_exception_count=1,
        zero_exception_count=1,
        constraint_count=1,
        nonbonded_metadata=ExactNonbondedMetadata(
            method="PME",
            method_code=4,
            cutoff_nm=1.0,
            use_switching_function=False,
            use_dispersion_correction=True,
            ewald_error_tolerance=5e-4,
        ),
        atoms=atoms,
        topology_bonds=topology_bonds,
        atom_order_hash=authoritative_atom_order_hash(atoms),
        topology_hash=authoritative_topology_hash(atoms, topology_bonds),
        gromacs_atom_order_hash=gromacs_atom_order_hash(atoms),
        gromacs_topology_hash=gromacs_topology_hash(atoms, topology_bonds),
        exception_hash=exception_hash(exceptions),
        particle_hash=particle_hash(particles),
        particles=particles,
        exceptions=exceptions,
        constraints=(ExactConstraintRecord(constraint_index=0, i=1, j=3, length_nm=0.109),),
    )


def _topology_text(*, pair: str | None = None, extra: list[str] | None = None) -> str:
    """Return a minimal split-ready GROMACS topology."""
    lines = [
        "[ defaults ]",
        "; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ",
        "1 2 yes 1.0 1.0",
        "[ atomtypes ]",
        "[ moleculetype ]",
        "MOL 3",
        "[ atoms ]",
        "1 C 1 MOL C1 1 0.2 12.0",
        "2 O 1 MOL O1 1 -0.5 16.0",
        "3 H 1 MOL H1 1 0.1 1.0",
    ]
    if pair is not None:
        lines.extend(["[ pairs ]", pair])
    lines.extend(extra or ["[ bonds ]", "1 2", "2 3"])
    lines.extend(["[ system ]", "x", "[ molecules ]", "MOL 1"])
    return "\n".join(lines) + "\n"


def _raw_pair() -> str:
    """Return a raw local pair row equivalent to the base sidecar."""
    return "1 2 2 1 0.2 -0.5 0.305 0.4472135955"


def _atom(index: int, name: str = "A", residue: str = "MOL") -> AtomIdentity:
    """Return a compact atom identity."""
    return AtomIdentity(index=index, name=name, residue_name=residue, residue_id="1", chain_id="A")


def _sidecar_from_records(
    particles: tuple[ExactParticleRecord, ...],
    exceptions: tuple[ExactExceptionRecord, ...],
    atoms: tuple[AtomIdentity, ...] | None = None,
    topology_bonds: tuple[ExactTopologyBondRecord, ...] = (),
) -> ExactExceptionSidecar:
    """Build a sidecar from explicit particles and exceptions."""
    default_names = ("C1", "O1")
    resolved_atoms = atoms or tuple(
        AtomIdentity(
            index=particle.index,
            name=default_names[(particle.index - 1) % len(default_names)],
            residue_name="MOL",
            residue_id="1",
            chain_id="A",
        )
        for particle in particles
    )
    return ExactExceptionSidecar(
        particle_count=len(particles),
        exception_count=len(exceptions),
        nonzero_exception_count=sum(1 for record in exceptions if not record.is_zero),
        zero_exception_count=sum(1 for record in exceptions if record.is_zero),
        constraint_count=0,
        nonbonded_metadata=ExactNonbondedMetadata(
            method="PME",
            method_code=4,
            cutoff_nm=1.0,
            use_switching_function=False,
            use_dispersion_correction=True,
            ewald_error_tolerance=5e-4,
        ),
        atoms=resolved_atoms,
        topology_bonds=topology_bonds,
        atom_order_hash=authoritative_atom_order_hash(resolved_atoms),
        topology_hash=authoritative_topology_hash(resolved_atoms, topology_bonds),
        gromacs_atom_order_hash=gromacs_atom_order_hash(resolved_atoms),
        gromacs_topology_hash=gromacs_topology_hash(resolved_atoms, topology_bonds),
        exception_hash=exception_hash(exceptions),
        particle_hash=particle_hash(particles),
        particles=particles,
        exceptions=exceptions,
    )


def _glycoprotein_two_water_sidecar() -> ExactExceptionSidecar:
    """Return MOL0 plus two repeated three-atom waters exact exceptions."""
    particles = tuple(
        ExactParticleRecord(index=index, charge_e=charge, sigma_nm=0.3, epsilon_kj_mol=0.4)
        for index, charge in enumerate((0.2, -0.5, -0.8, 0.4, 0.4, -0.8, 0.4, 0.4), start=1)
    )
    exceptions = [
        ExactExceptionRecord(
            exception_index=0,
            i=1,
            j=2,
            charge_product_e2=-0.1,
            sigma_nm=0.305,
            epsilon_kj_mol=0.4472135955,
            category="glycam_unscaled_14",
            atom_i=_atom(1),
            atom_j=_atom(2),
        )
    ]
    for left, right in ((3, 4), (3, 5), (4, 5), (6, 7), (6, 8), (7, 8)):
        exceptions.append(
            ExactExceptionRecord(
                exception_index=len(exceptions),
                i=left,
                j=right,
                charge_product_e2=0.0,
                sigma_nm=0.0,
                epsilon_kj_mol=0.0,
                category="zero_exclusion",
                atom_i=_atom(left, residue="HOH"),
                atom_j=_atom(right, residue="HOH"),
            )
        )
    atoms = (
        AtomIdentity(index=1, name="C1", residue_name="MOL", residue_id="1", chain_id="A"),
        AtomIdentity(index=2, name="O1", residue_name="MOL", residue_id="1", chain_id="A"),
        AtomIdentity(index=3, name="O", residue_name="HOH", residue_id="1", chain_id="A"),
        AtomIdentity(index=4, name="H1", residue_name="HOH", residue_id="1", chain_id="A"),
        AtomIdentity(index=5, name="H2", residue_name="HOH", residue_id="1", chain_id="A"),
        AtomIdentity(index=6, name="O", residue_name="HOH", residue_id="1", chain_id="A"),
        AtomIdentity(index=7, name="H1", residue_name="HOH", residue_id="1", chain_id="A"),
        AtomIdentity(index=8, name="H2", residue_name="HOH", residue_id="1", chain_id="A"),
    )
    topology_bonds = (ExactTopologyBondRecord(i=1, j=2, atom_i=atoms[0], atom_j=atoms[1]),)
    return _sidecar_from_records(particles, tuple(exceptions), atoms, topology_bonds)


def _glycoprotein_two_water_topology() -> str:
    """Return a topology with one MOL0 and two repeated MOL1 waters."""
    return (
        "\n".join(
            [
                "[ defaults ]",
                "1 2 yes 1.0 1.0",
                "[ moleculetype ]",
                "MOL0 3",
                "[ atoms ]",
                "1 C 1 MOL C1 1 0.2 12.0",
                "2 O 1 MOL O1 1 -0.5 16.0",
                "[ pairs ]",
                _raw_pair(),
                "[ bonds ]",
                "1 2",
                "[ moleculetype ]",
                "MOL1 2",
                "[ atoms ]",
                "1 O 1 HOH O 1 -0.8 16.0",
                "2 H 1 HOH H1 1 0.4 1.0",
                "3 H 1 HOH H2 1 0.4 1.0",
                "[ bonds ]",
                "[ system ]",
                "x",
                "[ molecules ]",
                "MOL0 1",
                "MOL1 2",
            ]
        )
        + "\n"
    )


def test_sidecar_validates_hashes_and_counts() -> None:
    """Sidecar schema should validate deterministic counts and hashes."""
    sidecar = _sidecar()

    assert sidecar.nonzero_exception_count == 1
    assert sidecar.zero_exception_count == 1
    assert sidecar.nonzero_exceptions[0].i == 1


def test_sidecar_rejects_hash_mismatch() -> None:
    """Sidecar schema should fail closed on exception hash mismatch."""
    data = _sidecar().model_dump(mode="python")
    data["exception_hash"] = "bad"

    with pytest.raises(ValueError, match="exception_hash"):
        ExactExceptionSidecar.model_validate(data)


def test_sidecar_rejects_tampered_atom_identity_hash() -> None:
    """Sidecar schema should fail closed when atom records are tampered."""
    sidecar = _sidecar()
    data = sidecar.model_dump(mode="python")
    data["atoms"][1]["name"] = "TAMPERED"

    with pytest.raises(ValueError, match="atom_order_hash"):
        ExactExceptionSidecar.model_validate(data)


def test_sidecar_rejects_tampered_topology_bond_hash() -> None:
    """Sidecar schema should fail closed when topology bond records are tampered."""
    sidecar = _sidecar()
    data = sidecar.model_dump(mode="python")
    data["topology_bonds"][1]["i"] = 1

    with pytest.raises(ValueError, match="topology_hash"):
        ExactExceptionSidecar.model_validate(data)


def test_sidecar_load_rejects_v1_with_migration_message(tmp_path) -> None:
    """Schema v1 sidecars should fail closed with regeneration guidance."""

    path = tmp_path / "exact_openmm_exceptions.json"
    data = _sidecar().model_dump(mode="json")
    data["schema_version"] = 1
    path.write_text(json.dumps(data), encoding="utf-8")

    with pytest.raises(ExactExportError, match="schema v1 is unsupported"):
        ExactExceptionSidecar.load(path)


def test_sidecar_load_rejects_future_schema_with_specific_version(tmp_path) -> None:
    """Unsupported schema errors should name the actual unsupported version."""

    path = tmp_path / "exact_openmm_exceptions.json"
    data = _sidecar().model_dump(mode="json")
    data["schema_version"] = 999
    path.write_text(json.dumps(data), encoding="utf-8")

    with pytest.raises(ExactExportError, match="schema v999 is unsupported"):
        ExactExceptionSidecar.load(path)


def test_bundle_blocks_raw_gromacs_export() -> None:
    """Exact bundle should reject accidental raw Interchange-style export."""
    bundle = ExactExportBundle(
        topology=object(),
        system=object(),
        positions=object(),
        private_baseline_interchange=object(),
        sidecar=_sidecar(),
    )

    with pytest.raises(ExactExportError, match="not a vanilla OpenFF Interchange"):
        bundle.to_gromacs(prefix="raw")


def test_patch_topology_writes_function2_pairs_and_exclusions(tmp_path) -> None:
    """Topology patch should write exact local function-2 pairs and exclusions."""
    top_path = tmp_path / "system.top"
    top_path.write_text(_topology_text(pair=_raw_pair()))

    audit = patch_gromacs_topology_with_exact_exceptions(top_path, _sidecar())
    patched = top_path.read_text()

    assert " no " in patched
    assert "     1      2 2" in patched
    assert "     2      3" in patched
    assert audit["patched_pair_count"] == 1
    assert audit["patched_exclusion_count"] == 2
    assert audit["zero_pairs_in_patched_pairs"] == 0


def test_patch_topology_collapses_multiple_raw_pair_sections(tmp_path) -> None:
    """Topology patch should write local pair sections only where needed."""
    top_path = tmp_path / "system.top"
    top_path.write_text(
        "\n".join(
            [
                "[ defaults ]",
                "1 2 yes 1.0 1.0",
                "[ moleculetype ]",
                "MOL 3",
                "[ atoms ]",
                "1 C 1 MOL C1 1 0.2 12.0",
                "2 O 1 MOL O1 1 -0.5 16.0",
                "3 H 1 MOL H1 1 0.1 1.0",
                "[ pairs ]",
                _raw_pair(),
                "[ bonds ]",
                "1 2",
                "2 3",
                "[ moleculetype ]",
                "SOL 2",
                "[ atoms ]",
                "1 H 1 SOL H1 1 0.0 1.0",
                "[ angles ]",
                "[ system ]",
                "x",
                "[ molecules ]",
                "MOL 1",
                "SOL 0",
            ]
        )
        + "\n"
    )

    audit = patch_gromacs_topology_with_exact_exceptions(top_path, _sidecar())
    patched = top_path.read_text()

    assert patched.count("[ pairs ]") == 1
    assert audit["patched_pair_count"] == 1


def test_patch_topology_preserves_split_molecule_constraints(tmp_path) -> None:
    """Topology patch should not inject global constraints into split molecule types."""
    top_path = tmp_path / "system.top"
    top_path.write_text(_topology_text(pair=_raw_pair()))

    patch_gromacs_topology_with_exact_exceptions(top_path, _sidecar())
    patched = top_path.read_text()

    assert "[ constraints ]" not in patched


def test_patch_topology_writes_first_molecule_zero_exclusions(tmp_path) -> None:
    """Topology patch should prevent first-molecule zero-exclusion leakage."""
    top_path = tmp_path / "system.top"
    top_path.write_text(
        _topology_text(pair=_raw_pair(), extra=["[ exclusions ]", "1 3", "[ bonds ]", "1 2", "2 3"])
    )

    patch_gromacs_topology_with_exact_exceptions(top_path, _sidecar())
    patched = top_path.read_text()

    assert patched.count("[ exclusions ]") == 1
    assert "     2      3" in patched


def test_patch_topology_maps_repeated_water_exclusions_locally(tmp_path) -> None:
    """Repeated waters should keep one local exclusion set that expands exactly."""
    top_path = tmp_path / "system.top"
    top_path.write_text(_glycoprotein_two_water_topology())

    audit = patch_gromacs_topology_with_exact_exceptions(
        top_path, _glycoprotein_two_water_sidecar()
    )
    patched = top_path.read_text()

    assert patched.count("[ pairs ]") == 1
    assert patched.count("[ exclusions ]") == 2
    assert audit["patched_pair_count"] == 1
    assert audit["patched_exclusion_count"] == 7
    assert audit["patched_local_exclusion_row_count"] == 4
    assert "     1      2" in patched
    assert "     1      3" in patched
    assert "     2      3" in patched


def test_patch_topology_rejects_repeated_molecule_semantic_mismatch(tmp_path) -> None:
    """Repeated molecule copies must have identical local exception parameters."""
    top_path = tmp_path / "system.top"
    top_path.write_text(
        "\n".join(
            [
                "[ defaults ]",
                "1 2 yes 1.0 1.0",
                "[ moleculetype ]",
                "MOL 3",
                "[ atoms ]",
                "1 C 1 MOL C1 1 0.2 12.0",
                "2 O 1 MOL O1 1 -0.5 16.0",
                "[ pairs ]",
                _raw_pair(),
                "[ system ]",
                "x",
                "[ molecules ]",
                "MOL 2",
            ]
        )
        + "\n"
    )
    particles = tuple(
        ExactParticleRecord(index=index, charge_e=charge, sigma_nm=0.3, epsilon_kj_mol=0.4)
        for index, charge in enumerate((0.2, -0.5, 0.2, -0.5), start=1)
    )
    exceptions = (
        ExactExceptionRecord(
            exception_index=0,
            i=1,
            j=2,
            charge_product_e2=-0.1,
            sigma_nm=0.305,
            epsilon_kj_mol=0.4472135955,
            category="glycam_unscaled_14",
            atom_i=_atom(1),
            atom_j=_atom(2),
        ),
        ExactExceptionRecord(
            exception_index=1,
            i=3,
            j=4,
            charge_product_e2=-0.2,
            sigma_nm=0.305,
            epsilon_kj_mol=0.4472135955,
            category="glycam_unscaled_14",
            atom_i=_atom(3),
            atom_j=_atom(4),
        ),
    )

    with pytest.raises(RuntimeError, match="inconsistent local exception semantics"):
        patch_gromacs_topology_with_exact_exceptions(
            top_path, _sidecar_from_records(particles, exceptions)
        )


def test_patch_topology_rejects_cross_copy_exception(tmp_path) -> None:
    """Exact exceptions spanning molecule copies must fail closed."""
    top_path = tmp_path / "system.top"
    top_path.write_text(
        "\n".join(
            [
                "[ defaults ]",
                "1 2 yes 1.0 1.0",
                "[ moleculetype ]",
                "MOL 3",
                "[ atoms ]",
                "1 C 1 MOL C1 1 0.2 12.0",
                "2 O 1 MOL O1 1 -0.5 16.0",
                "[ pairs ]",
                _raw_pair(),
                "[ system ]",
                "x",
                "[ molecules ]",
                "MOL 2",
            ]
        )
        + "\n"
    )
    particles = tuple(
        ExactParticleRecord(index=index, charge_e=charge, sigma_nm=0.3, epsilon_kj_mol=0.4)
        for index, charge in enumerate((0.2, -0.5, 0.2, -0.5), start=1)
    )
    exceptions = (
        ExactExceptionRecord(
            exception_index=0,
            i=2,
            j=3,
            charge_product_e2=-0.1,
            sigma_nm=0.305,
            epsilon_kj_mol=0.4472135955,
            category="glycam_unscaled_14",
            atom_i=_atom(2),
            atom_j=_atom(3),
        ),
    )

    with pytest.raises(RuntimeError, match="inter-copy or inter-type"):
        patch_gromacs_topology_with_exact_exceptions(
            top_path, _sidecar_from_records(particles, exceptions)
        )


def test_patch_topology_fails_on_pair_mismatch(tmp_path) -> None:
    """Topology patch should fail closed when raw pairs do not match sidecar pairs."""
    top_path = tmp_path / "system.top"
    top_path.write_text(_topology_text(pair="1 3 2 1 0.2 0.1 0.0 0.0"))

    with pytest.raises(RuntimeError, match="raw pair set"):
        patch_gromacs_topology_with_exact_exceptions(top_path, _sidecar())


def test_patch_topology_preflight_failure_leaves_file_unchanged(tmp_path) -> None:
    """Preflight atom identity mismatches should not write partial patches."""

    top_path = tmp_path / "system.top"
    original = _topology_text(pair=_raw_pair())
    top_path.write_text(original)
    bad_atom = AtomIdentity(index=2, name="BAD", residue_name="MOL", residue_id="1", chain_id="")
    stale = _sidecar().model_copy(
        update={"atoms": (_sidecar().atoms[0], bad_atom, _sidecar().atoms[2])}
    )

    with pytest.raises(RuntimeError, match="atom order/identity mismatch"):
        patch_gromacs_topology_with_exact_exceptions(top_path, stale)

    assert top_path.read_text() == original


def test_patch_topology_reordered_raw_baseline_fails_before_write(tmp_path) -> None:
    """A same-count reordered raw baseline should not validate against OpenMM identity."""

    top_path = tmp_path / "system.top"
    original = _topology_text(pair=_raw_pair()).replace(
        "1 C 1 MOL C1 1 0.2 12.0\n2 O 1 MOL O1 1 -0.5 16.0",
        "1 C 1 MOL O1 1 0.2 12.0\n2 O 1 MOL C1 1 -0.5 16.0",
    )
    top_path.write_text(original)

    with pytest.raises(RuntimeError, match="atom order/identity mismatch"):
        patch_gromacs_topology_with_exact_exceptions(top_path, _sidecar())

    assert top_path.read_text() == original


def test_patch_topology_stale_gromacs_identity_fails_before_write(tmp_path) -> None:
    """Saved output-specific identity should fail closed when reused for another topology."""

    top_path = tmp_path / "system.top"
    original = _topology_text(pair=_raw_pair())
    top_path.write_text(original)
    stale_atom = AtomIdentity(index=1, name="OLD", residue_name="MOL", residue_id="1", chain_id="")
    stale_atoms = (stale_atom, *_sidecar().atoms[1:])
    stale = _sidecar().model_copy(
        update={
            "gromacs_atoms": stale_atoms,
            "gromacs_topology_bonds": _sidecar().topology_bonds,
            "gromacs_atom_order_hash": gromacs_atom_order_hash(stale_atoms),
            "gromacs_topology_hash": gromacs_topology_hash(stale_atoms, _sidecar().topology_bonds),
        }
    )

    with pytest.raises(RuntimeError, match="atom order/identity mismatch"):
        patch_gromacs_topology_with_exact_exceptions(top_path, stale)

    assert top_path.read_text() == original


def test_patch_topology_accepts_water_residue_normalization(tmp_path) -> None:
    """Water aliasing and water residue renumbering are the only accepted normalization."""

    top_path = tmp_path / "system.top"
    top_path.write_text(_glycoprotein_two_water_topology().replace("HOH", "SOL"))
    sidecar = _glycoprotein_two_water_sidecar()
    water_alias_atoms = tuple(
        (
            atom.model_copy(update={"residue_name": "WAT", "residue_id": str(atom.index)})
            if atom.index >= 3
            else atom
        )
        for atom in sidecar.atoms
    )
    sidecar = sidecar.model_copy(update={"atoms": water_alias_atoms})

    audit = patch_gromacs_topology_with_exact_exceptions(top_path, sidecar)

    assert audit["patched_exclusion_count"] == 7


def test_patch_topology_handles_noncontiguous_same_type_rows(tmp_path) -> None:
    """Repeated noncontiguous molecule rows should use cumulative copy ordinals."""

    top_path = tmp_path / "system.top"
    top_path.write_text(
        "\n".join(
            [
                "[ defaults ]",
                "1 2 yes 1.0 1.0",
                "[ moleculetype ]",
                "MOL 3",
                "[ atoms ]",
                "1 C 1 MOL C1 1 0.2 12.0",
                "2 O 1 MOL O1 1 -0.5 16.0",
                "[ pairs ]",
                _raw_pair(),
                "[ moleculetype ]",
                "SOL 2",
                "[ atoms ]",
                "1 H 1 SOL H1 1 0.0 1.0",
                "[ system ]",
                "x",
                "[ molecules ]",
                "MOL 1",
                "SOL 1",
                "MOL 1",
            ]
        )
        + "\n"
    )
    atoms = (
        AtomIdentity(index=1, name="C1", residue_name="MOL", residue_id="1", chain_id=""),
        AtomIdentity(index=2, name="O1", residue_name="MOL", residue_id="1", chain_id=""),
        AtomIdentity(index=3, name="H1", residue_name="SOL", residue_id="1", chain_id=""),
        AtomIdentity(index=4, name="C1", residue_name="MOL", residue_id="1", chain_id=""),
        AtomIdentity(index=5, name="O1", residue_name="MOL", residue_id="1", chain_id=""),
    )
    particles = tuple(
        ExactParticleRecord(index=index, charge_e=charge, sigma_nm=0.3, epsilon_kj_mol=0.4)
        for index, charge in enumerate((0.2, -0.5, 0.0, 0.2, -0.5), start=1)
    )
    exceptions = (
        ExactExceptionRecord(
            exception_index=0,
            i=1,
            j=2,
            charge_product_e2=-0.1,
            sigma_nm=0.305,
            epsilon_kj_mol=0.4472135955,
            category="glycam_unscaled_14",
            atom_i=atoms[0],
            atom_j=atoms[1],
        ),
        ExactExceptionRecord(
            exception_index=1,
            i=4,
            j=5,
            charge_product_e2=-0.1,
            sigma_nm=0.305,
            epsilon_kj_mol=0.4472135955,
            category="glycam_unscaled_14",
            atom_i=atoms[3],
            atom_j=atoms[4],
        ),
    )

    audit = patch_gromacs_topology_with_exact_exceptions(
        top_path, _sidecar_from_records(particles, exceptions, atoms)
    )

    assert audit["patched_pair_count"] == 2


def test_exact_bundle_export_rejects_component_info_before_writes(tmp_path) -> None:
    """Exact route should fail early instead of silently returning empty posres."""

    bundle = ExactExportBundle(
        topology=object(),
        system=object(),
        positions=object(),
        private_baseline_interchange=object(),
        sidecar=_sidecar(),
    )
    exporter = GromacsExporter(bundle, object(), component_info=object())

    with pytest.raises(ValueError, match="does not support component_info"):
        exporter.export(tmp_path, prefix="exact")

    assert not any(tmp_path.iterdir())
