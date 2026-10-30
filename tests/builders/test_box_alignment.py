"""How the build turns the fixed solutes before sizing a rectangular box."""

from __future__ import annotations

import copy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from polyzymd.builders.system_builder import SystemBuilder

QUICKSTART = Path(__file__).resolve().parents[2] / "examples" / "quickstart"


class _Topology:
    """The part of an OpenFF Topology that the alignment reads and writes."""

    def __init__(self, positions_nm: np.ndarray, atomic_numbers: list[int]) -> None:
        from openff.units import Quantity

        self.positions = positions_nm
        self.atoms = [
            SimpleNamespace(atomic_number=n, mass=Quantity(1.008 if n == 1 else 12.011, "dalton"))
            for n in atomic_numbers
        ]

    def get_positions(self):
        from openff.units import Quantity

        return Quantity(self.positions, "nanometer")

    def set_positions(self, quantity) -> None:
        self.positions = np.asarray(quantity.m_as("nanometer"))


def _write_pdb(path: Path, positions_nm: np.ndarray) -> None:
    lines = [
        f"ATOM  {i + 1:5d}  CA  ALA A{i + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
        for i, (x, y, z) in enumerate(positions_nm * 10.0)
    ]
    path.write_text("".join(lines) + "END\n")


def _aligned(tmp_path: Path, chain_seed: int) -> tuple[dict, np.ndarray, Path]:
    """Align a rod-shaped protein carrying a grafted chain drawn from ``chain_seed``."""
    rng = np.random.default_rng(1)
    protein = rng.uniform(-0.5, 0.5, size=(300, 3)) * np.array([0.4, 1.0, 6.0])
    protein = protein @ np.linalg.qr(rng.normal(size=(3, 3)))[0].T
    chain = np.random.default_rng(chain_seed).uniform(-2.0, 2.0, size=(40, 3)) + [2.0, 0.0, 0.0]
    ligand = np.array([[0.0, 0.0, 0.0]])
    positions = np.vstack([protein, chain, ligand])
    pdb = tmp_path / f"modified_protein_{chain_seed}.pdb"
    _write_pdb(pdb, positions[:340])
    builder = SimpleNamespace(
        _combined_topology=_Topology(positions, [6] * len(positions)),
        _modified_protein=SimpleNamespace(
            record=SimpleNamespace(segments=["PROT"] * 300 + ["P001"] * 40), pdb_path=pdb
        ),
        _build_provenance={},
    )
    SystemBuilder._align_for_box(builder, SimpleNamespace(alignment="principal_axes"))
    return builder._build_provenance["box_alignment"], builder._combined_topology.positions, pdb


def test_the_rotation_ignores_grafted_chains_and_the_modified_pdb_is_turned(tmp_path) -> None:
    first, turned, pdb = _aligned(tmp_path, chain_seed=2)
    second, _, _ = _aligned(tmp_path, chain_seed=3)
    assert np.allclose(first["rotation"], second["rotation"])
    assert first["excluded_attachment_atoms"] == 40
    written = np.array(
        [
            [float(line[c : c + 8]) for c in (30, 38, 46)]
            for line in pdb.read_text().splitlines()
            if line.startswith("ATOM")
        ]
    )
    assert np.allclose(written / 10.0, turned[:340], atol=1e-3)


def test_the_dry_run_estimate_uses_the_build_alignment(monkeypatch) -> None:
    """Mass-weighted heavy atoms, as the build turns the protein, not every atom equally."""
    import MDAnalysis as mda
    import yaml

    from polyzymd.cli import main
    from polyzymd.config.schema import SimulationConfig
    from polyzymd.utils.boxvectors import (
        describe_box_plan,
        fixed_solute_alignment,
        plan_rectangular_box,
    )

    data = yaml.safe_load((QUICKSTART / "config.yaml").read_text())
    data["solutes"][0]["pdb_path"] = str(QUICKSTART / "trpcage.pdb")
    data["solvent"]["box"] = {"shape": "rectangular", "padding": {"x": 1.0, "y": 1.2, "z": 1.6}}
    config = SimulationConfig.model_validate(copy.deepcopy(data))
    printed: list[str] = []
    monkeypatch.setattr(main, "colored_echo", lambda text="", **_: printed.append(text))
    main._echo_box_plan(config, phase="build")

    atoms = mda.Universe(str(QUICKSTART / "trpcage.pdb")).atoms
    positions = atoms.positions / 10.0
    numbers = np.array([1 if a.element.upper() == "H" else 6 for a in atoms])
    rotation, center = fixed_solute_alignment(positions, atoms.masses, numbers)
    plan = plan_rectangular_box(
        (positions - center) @ rotation.T + center, (1.0, 1.2, 1.6), margin_nm=0.2
    )
    assert describe_box_plan(plan) in printed[-1]


def test_padding_is_logged_as_numbers() -> None:
    from polyzymd.config.schema import BoxConfig

    assert BoxConfig(shape="rectangular", padding={"x": 1, "y": 1.2, "z": 1.6}).padding_label == (
        "1 / 1.2 / 1.6"
    )
    assert BoxConfig(padding=1.2).padding_label == "1.2"


def test_the_cutoff_reads_modification_force_fields_once(monkeypatch) -> None:
    from polyzymd.config import schema

    cutoffs = {"sage.offxml": 0.9, "polymer.offxml": 1.2}
    calls: list[str] = []

    def cutoff(name: str):
        calls.append(name)
        return cutoffs.get(name)

    monkeypatch.setattr(schema, "_vdw_cutoff_nm", cutoff)
    sets = SimpleNamespace(
        as_dict=lambda: {"polymer": SimpleNamespace(force_field=["polymer.offxml"])}
    )
    config = SimpleNamespace(
        solvent=SimpleNamespace(force_field=["sage.offxml"]),
        solutes=[SimpleNamespace(force_field=["sage.offxml"], modification_parameters=sets)],
    )
    assert schema.SimulationConfig._nonbonded_cutoff_nm(config) == pytest.approx(1.2)
    assert "polymer.offxml" in calls
    # The cache lives at module level, so it is shared across configs.
    monkeypatch.undo()
    assert hasattr(schema._vdw_cutoff_nm, "cache_info")
