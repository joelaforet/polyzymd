"""Tests for Phase 0-1 conjugation skeleton."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

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


def _minimal_attachment_data(enabled: bool = True) -> dict:
    """Build a minimal valid conjugation attachment dictionary.

    Parameters
    ----------
    enabled : bool, optional
        Whether the attachment is enabled, by default True.

    Returns
    -------
    dict
        Dictionary suitable for :class:`ConjugationConfig` validation.
    """
    return {
        "name": "lys23-peg",
        "enabled": enabled,
        "site": {"chain_id": "A", "residue_name": "LYS", "residue_number": 23},
        "moiety": {
            "name": "PEG",
            "force_field": "openff-2.0.0.offxml",
            "smiles": "COCCO",
        },
        "mechanism": {"name": "amide"},
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

    def test_enabled_conjugation_with_attachment_parses_without_mode(self):
        """Enabled conjugation parses when at least one attachment is enabled."""
        data = _minimal_simulation_config_data()
        data["conjugation"] = {
            "enabled": True,
            "attachments": [_minimal_attachment_data()],
            "ccd_pablo": {"lookup_policy": "offline_cached"},
        }

        config = SimulationConfig.model_validate(data)

        assert config.conjugation is not None
        assert config.conjugation.enabled is True
        assert config.conjugation.attachments[0].moiety.name == "PEG"

    def test_enabled_conjugation_with_no_attachments_rejected(self):
        """Enabled conjugation requires at least one enabled attachment."""
        data = _minimal_simulation_config_data()
        data["conjugation"] = {"enabled": True}

        with pytest.raises(ValidationError, match="at least one enabled attachment"):
            SimulationConfig.model_validate(data)

    def test_enabled_conjugation_with_only_disabled_attachments_rejected(self):
        """Enabled conjugation ignores disabled attachments for the opt-in contract."""
        with pytest.raises(ValidationError, match="at least one enabled attachment"):
            ConjugationConfig(enabled=True, attachments=[_minimal_attachment_data(enabled=False)])

    def test_disabled_attachment_does_not_require_force_field_owner(self):
        """Disabled attachments remain inert without force-field ownership."""
        enabled = _minimal_attachment_data()
        disabled = _minimal_attachment_data(enabled=False)
        disabled["name"] = "disabled"
        disabled["moiety"].pop("force_field")

        config = ConjugationConfig(enabled=True, attachments=[enabled, disabled])

        assert config.attachments[1].moiety.force_field is None

    def test_stale_mode_rejected_as_extra_field(self):
        """Old mode fields are rejected by the public conjugation config."""
        data = _minimal_simulation_config_data()
        data["conjugation"] = {
            "enabled": True,
            "mode": "construct",
            "attachments": [_minimal_attachment_data()],
        }

        with pytest.raises(ValidationError, match="Extra inputs are not permitted"):
            SimulationConfig.model_validate(data)

    def test_stale_structure_path_rejected_as_extra_field(self):
        """Old structure path fields are rejected by the public conjugation config."""
        data = _minimal_simulation_config_data()
        data["conjugation"] = {
            "enabled": True,
            "structure_path": "prebuilt_conjugate.pdb",
            "attachments": [_minimal_attachment_data()],
        }

        with pytest.raises(ValidationError, match="Extra inputs are not permitted"):
            SimulationConfig.model_validate(data)

    def test_stale_moiety_role_rejected_as_extra_field(self):
        """Old moiety role labels are rejected by the public conjugation config."""
        data = _minimal_simulation_config_data()
        attachment = _minimal_attachment_data()
        attachment["moiety"]["role"] = "moiety"
        data["conjugation"] = {
            "enabled": True,
            "attachments": [attachment],
        }

        with pytest.raises(ValidationError, match="Extra inputs are not permitted"):
            SimulationConfig.model_validate(data)

    def test_invalid_ccd_policy_rejected(self):
        """Unknown CCD/Pablo policies fail validation."""
        with pytest.raises(ValidationError):
            ConjugationConfig(
                enabled=True,
                attachments=[_minimal_attachment_data()],
                ccd_pablo={"lookup_policy": "always_guess"},
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
    """SystemBuilder compatibility tests for conjugate chain handling."""

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
