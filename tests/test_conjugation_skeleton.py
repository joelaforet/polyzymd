"""Tests for Phase 0-1 conjugation skeleton."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest
from pydantic import ValidationError

from polyzymd.builders.conjugation.exceptions import ConjugationNotImplementedError
from polyzymd.config.schema import (
    CcdLookupPolicy,
    ConjugationCcdPabloPolicyConfig,
    ConjugationChainPolicyConfig,
    ConjugationConfig,
    SimulationConfig,
)


def _minimal_simulation_config_data() -> dict:
    """Build a minimal valid simulation config dictionary.

    Returns
    -------
    dict
        Dictionary suitable for :meth:`SimulationConfig.model_validate`.
    """
    return {
        "name": "conjugation-skeleton",
        "engine": "openmm",
        "enzyme": {"name": "enzyme", "pdb_path": "enzyme.pdb"},
        "thermodynamics": {"temperature": 300.0},
        "simulation_phases": {
            "equilibration_stages": [
                {
                    "name": "eq",
                    "duration": 0.1,
                    "samples": 10,
                    "ensemble": "NVT",
                    "temperature": 300.0,
                }
            ],
            "production": {
                "ensemble": "NPT",
                "duration": 1.0,
                "samples": 10,
                "time_step": 2.0,
                "report_interval": 1000,
                "checkpoint_interval": 5000,
            },
        },
    }


class TestConjugationConfigParsing:
    """Configuration parsing tests for the conjugation skeleton."""

    def test_simulation_config_accepts_absent_conjugation(self):
        """Top-level conjugation remains optional for legacy configs."""
        config = SimulationConfig.model_validate(_minimal_simulation_config_data())
        assert config.conjugation is None

    def test_simulation_config_accepts_disabled_conjugation(self):
        """Explicit disabled conjugation parses without changing defaults."""
        data = _minimal_simulation_config_data()
        data["conjugation"] = {"enabled": False}

        config = SimulationConfig.model_validate(data)

        assert config.conjugation is not None
        assert config.conjugation.enabled is False
        assert config.conjugation.chain_policy.protein_chain == "A"
        assert config.conjugation.chain_policy.moiety_chain == "C"

    def test_enabled_ingest_existing_block_parses(self):
        """Enabled ingest_existing config parses with Pablo policy placeholders."""
        data = _minimal_simulation_config_data()
        data["conjugation"] = {
            "enabled": True,
            "mode": "ingest_existing",
            "source_pdb_path": "prebuilt_conjugate.pdb",
            "ccd_pablo": {"lookup_policy": "offline_cached"},
        }

        config = SimulationConfig.model_validate(data)

        assert config.conjugation is not None
        assert config.conjugation.enabled is True
        assert config.conjugation.mode.value == "ingest_existing"
        assert config.conjugation.source_pdb_path.name == "prebuilt_conjugate.pdb"

    def test_enabled_ingest_existing_defaults_to_enzyme_pdb_path(self):
        """Ingest-existing configs should use enzyme.pdb_path when source is omitted."""
        data = _minimal_simulation_config_data()
        data["conjugation"] = {"enabled": True, "mode": "ingest_existing"}

        config = SimulationConfig.model_validate(data)

        assert config.conjugation is not None
        assert config.conjugation.source_pdb_path == config.enzyme.pdb_path

    def test_enabled_construct_block_parses(self):
        """Enabled construct config parses attachment placeholders."""
        data = _minimal_simulation_config_data()
        data["conjugation"] = {
            "enabled": True,
            "mode": "construct",
            "attachments": [
                {
                    "name": "lys23-peg",
                    "site": {"chain_id": "A", "residue_name": "LYS", "residue_number": 23},
                    "moiety": {"name": "PEG", "smiles": "COCCO"},
                    "mechanism": {"name": "amide"},
                }
            ],
        }

        config = SimulationConfig.model_validate(data)

        assert config.conjugation is not None
        assert config.conjugation.mode.value == "construct"
        assert config.conjugation.attachments[0].moiety.name == "PEG"

    def test_invalid_mode_rejected(self):
        """Unknown conjugation modes fail validation."""
        with pytest.raises(ValidationError):
            ConjugationConfig(enabled=True, mode="unsupported")

    def test_invalid_ccd_policy_rejected(self):
        """Unknown CCD/Pablo policies fail validation."""
        with pytest.raises(ValidationError):
            ConjugationConfig(
                enabled=True,
                ccd_pablo={"lookup_policy": "always_guess"},
            )

    def test_ingest_existing_accepts_attachments(self):
        """Ingest-existing mode can accept additional attachment requests."""
        config = ConjugationConfig(
            enabled=True,
            mode="ingest_existing",
            attachments=[
                {
                    "name": "lys23-peg",
                    "moiety": {"name": "PEG", "smiles": "COCCO"},
                }
            ],
        )

        assert config.attachments[0].name == "lys23-peg"

    def test_construct_accepts_source_pdb_path(self):
        """Construct mode can ingest a prepared base PDB before new attachments."""
        config = ConjugationConfig(
            enabled=True,
            mode="construct",
            source_pdb_path="prebuilt_conjugate.pdb",
        )

        assert config.source_pdb_path == Path("prebuilt_conjugate.pdb")

    def test_ingest_existing_rejects_disabled_pablo_policy(self):
        """Ingest-existing mode requires Pablo validation to be enabled."""
        with pytest.raises(ValidationError, match="ccd_pablo.enabled"):
            ConjugationConfig(
                enabled=True,
                mode="ingest_existing",
                ccd_pablo={"enabled": False},
            )

    def test_source_backed_construct_rejects_disabled_pablo_policy(self):
        """Construct workflows with a source PDB require Pablo ingestion."""
        with pytest.raises(ValidationError, match="ccd_pablo.enabled"):
            ConjugationConfig(
                enabled=True,
                mode="construct",
                source_pdb_path="prebuilt_conjugate.pdb",
                ccd_pablo={"enabled": False},
            )

    def test_pablo_policy_defaults_to_auto_download(self):
        """CCD/Pablo lookup policy should default to auto-download."""
        policy = ConjugationCcdPabloPolicyConfig()

        assert policy.lookup_policy == CcdLookupPolicy.AUTO_DOWNLOAD

    def test_pablo_crosslink_config_validates_shape(self):
        """Crosslink config should require paired residues, atoms, and leaving groups."""
        policy = ConjugationCcdPabloPolicyConfig(
            crosslinks=[
                {
                    "residues": ("LYX", "NHX"),
                    "linking_atoms": ("NZ", "C"),
                    "leaving_atoms": (("HZ1",), ("O1", "O2")),
                    "bond_order": 1,
                }
            ]
        )

        assert policy.crosslinks[0].residues == ("LYX", "NHX")
        assert policy.crosslinks[0].leaving_atoms == (("HZ1",), ("O1", "O2"))

    def test_pablo_crosslink_config_rejects_invalid_shape(self):
        """Invalid crosslink tuple lengths should fail validation."""
        with pytest.raises(ValidationError):
            ConjugationCcdPabloPolicyConfig(
                crosslinks=[
                    {
                        "residues": ("LYX", "NHX"),
                        "linking_atoms": ("NZ",),
                        "leaving_atoms": (("HZ1",), ("O1", "O2")),
                    }
                ]
            )

    def test_chain_policy_defaults(self):
        """Chain policy defaults preserve the PolyzyMD convention."""
        policy = ConjugationChainPolicyConfig()
        assert policy.protein_chain == "A"
        assert policy.substrate_chain == "B"
        assert policy.moiety_chain == "C"
        assert policy.solvent_start_chain == "D"

    def test_chain_policy_rejects_overlap(self):
        """Component and solvent chain IDs must not overlap."""
        with pytest.raises(ValidationError, match="must not overlap"):
            ConjugationChainPolicyConfig(solvent_start_chain="C")

    def test_chain_policy_rejects_nonalpha(self):
        """Chain IDs must be alphabetic single characters."""
        with pytest.raises(ValidationError, match="alphabetic"):
            ConjugationChainPolicyConfig(protein_chain="1")


class TestSystemBuilderConjugationHook:
    """SystemBuilder compatibility tests for the pre-solvation hook."""

    def test_absent_conjugation_returns_combined_topology(self):
        """Legacy configs without conjugation remain unchanged."""
        from polyzymd.builders.system_builder import SystemBuilder

        topology = object()
        builder = SystemBuilder.__new__(SystemBuilder)
        builder._combined_topology = topology

        config = MagicMock()
        config.conjugation = None

        assert builder._apply_conjugation(config) is topology

    def test_disabled_conjugation_returns_combined_topology(self):
        """Disabled conjugation remains a no-op through SystemBuilder."""
        from polyzymd.builders.system_builder import SystemBuilder

        topology = object()
        builder = SystemBuilder.__new__(SystemBuilder)
        builder._combined_topology = topology

        config = MagicMock()
        config.conjugation = ConjugationConfig(enabled=False)

        assert builder._apply_conjugation(config) is topology

    def test_enabled_conjugation_fails_before_solvation(self):
        """Enabled legacy paths fail with a public workflow pointer."""
        from polyzymd.builders.system_builder import SystemBuilder

        builder = SystemBuilder.__new__(SystemBuilder)
        builder._combined_topology = object()
        builder._enzyme_topology = object()
        builder._working_dir = None
        builder._n_enzyme_molecules = 1
        builder._n_substrate_molecules = 0
        builder._n_polymer_chains = 0

        config = MagicMock()
        config.conjugation = ConjugationConfig(enabled=True, mode="construct")

        with pytest.raises(ConjugationNotImplementedError, match="build_conjugate_from_config"):
            builder._apply_conjugation(config)

    def test_pack_polymers_preserves_conjugate_chain_ids(self, monkeypatch):
        """PACKMOL should not overwrite chain IDs on pre-linked conjugates."""
        from polyzymd.builders.system_builder import SystemBuilder
        from polyzymd.utils import packmol as packmol_module

        class FakeAtom:
            def __init__(self, chain_id: str) -> None:
                self.metadata = {"chain_id": chain_id}

        class FakeMolecule:
            def __init__(self, atoms: list[FakeAtom]) -> None:
                self.atoms = atoms

        class FakeTopology:
            def __init__(self, molecules: list[FakeMolecule]) -> None:
                self.molecules = molecules
                self.n_atoms = sum(len(mol.atoms) for mol in molecules)

        protein_atom = FakeAtom("A")
        attached_polymer_atom = FakeAtom("C")
        free_polymer_atom = FakeAtom("X")
        packed_topology = FakeTopology(
            [
                FakeMolecule([protein_atom, attached_polymer_atom]),
                FakeMolecule([free_polymer_atom]),
            ]
        )

        def fake_pack_polymers(**kwargs):
            return packed_topology

        monkeypatch.setattr(packmol_module, "pack_polymers", fake_pack_polymers)

        builder = SystemBuilder.__new__(SystemBuilder)
        builder._combined_topology = object()
        builder._polymer_molecules = [object()]
        builder._polymer_counts = [1]
        builder._preserve_enzyme_chain_ids = True

        assert builder.pack_polymers(box_vectors_nm=[1.0, 1.0, 1.0]) is packed_topology
        assert protein_atom.metadata["chain_id"] == "A"
        assert attached_polymer_atom.metadata["chain_id"] == "C"
        assert free_polymer_atom.metadata["chain_id"] == "X"
