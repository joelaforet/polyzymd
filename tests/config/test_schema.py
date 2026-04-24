"""Test that all public modules can be imported."""

from pathlib import Path

import pytest
from pydantic import ValidationError

from polyzymd.config.schema import SimulationConfig


@pytest.fixture
def minimal_config_data():
    """Provide minimal valid SimulationConfig input data."""
    return {
        "name": "test_simulation",
        "enzyme": {"name": "TestEnzyme", "pdb_path": "test.pdb"},
        "thermodynamics": {"temperature": 300.0},
        "simulation_phases": {
            "equilibration_stages": [
                {
                    "name": "eq1",
                    "duration": 0.1,
                    "temperature": 300.0,
                    "ensemble": "NVT",
                }
            ],
            "production": {
                "ensemble": "NPT",
                "duration": 1.0,
                "samples": 10,
            },
        },
    }


class TestImports:
    """Test basic package imports."""

    def test_import_polyzymd(self):
        """Test main package import."""
        import polyzymd

        assert hasattr(polyzymd, "__version__")

    def test_import_config(self):
        """Test config module import."""
        from polyzymd.config import SimulationConfig

        assert SimulationConfig is not None

    def test_import_conjugation_builder_without_rdkit(self, monkeypatch):
        """Importing the conjugation builder should not require RDKit."""
        import builtins
        import importlib

        import polyzymd.builders  # noqa: F401 - isolate this test to the conjugation builder

        real_import = builtins.__import__

        def block_rdkit_import(name, globals=None, locals=None, fromlist=(), level=0):
            """Block RDKit imports while allowing all other imports."""
            if name == "rdkit" or name.startswith("rdkit."):
                raise ImportError("RDKit intentionally blocked for import test")
            return real_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", block_rdkit_import)
        module = importlib.import_module("polyzymd.builders.conjugation.builder")

        assert module.CovalentModificationBuilder is not None

    def test_import_schema(self):
        """Test schema imports."""
        from polyzymd.config.schema import (
            EnzymeConfig,
            OutputConfig,
            PolymerConfig,
            SimulationPhasesConfig,
            SolventConfig,
            SubstrateConfig,
            ThermodynamicsConfig,
        )

        assert EnzymeConfig is not None
        assert SubstrateConfig is not None
        assert PolymerConfig is not None
        assert SolventConfig is not None
        assert ThermodynamicsConfig is not None
        assert SimulationPhasesConfig is not None
        assert OutputConfig is not None

    def test_version_format(self):
        """Test version string format."""
        import polyzymd

        version = polyzymd.__version__
        # Should be semver format: X.Y.Z
        parts = version.split(".")
        assert len(parts) >= 2, f"Version {version} should have at least major.minor"
        assert parts[0].isdigit(), f"Major version should be numeric: {parts[0]}"
        assert parts[1].isdigit(), f"Minor version should be numeric: {parts[1]}"


class TestConfigValidation:
    """Test configuration validation."""

    def test_enzyme_config_required_fields(self):
        """Test EnzymeConfig requires name and pdb_path."""
        from pydantic import ValidationError

        from polyzymd.config.schema import EnzymeConfig

        # Should fail without required fields
        with pytest.raises(ValidationError):
            EnzymeConfig()

        # Should succeed with required fields
        config = EnzymeConfig(name="TestEnzyme", pdb_path="test.pdb")
        assert config.name == "TestEnzyme"
        assert config.pdb_path == Path("test.pdb")

    def test_monomer_probabilities_sum(self):
        """Test that monomer probabilities must sum to 1.0."""
        from pydantic import ValidationError

        from polyzymd.config.schema import MonomerSpec, PolymerConfig

        # Should fail if probabilities don't sum to 1.0
        with pytest.raises(ValidationError):
            PolymerConfig(
                type_prefix="TEST",
                monomers=[
                    MonomerSpec(label="A", probability=0.5, name="MonA"),
                    MonomerSpec(label="B", probability=0.3, name="MonB"),  # Sum = 0.8
                ],
                length=5,
                count=2,
                sdf_directory=Path("/tmp/test"),  # Required for cached mode
            )

        # Should succeed if probabilities sum to 1.0
        config = PolymerConfig(
            type_prefix="TEST",
            monomers=[
                MonomerSpec(label="A", probability=0.7, name="MonA"),
                MonomerSpec(label="B", probability=0.3, name="MonB"),
            ],
            length=5,
            count=2,
            sdf_directory=Path("/tmp/test"),  # Required for cached mode
        )
        assert len(config.monomers) == 2

    def test_thermodynamics_config(self):
        """Test ThermodynamicsConfig validation and defaults."""
        from pydantic import ValidationError

        from polyzymd.config.schema import ThermodynamicsConfig

        # Temperature is required (no default)
        with pytest.raises(ValidationError):
            ThermodynamicsConfig()

        # Should succeed with temperature provided
        config = ThermodynamicsConfig(temperature=300.0)
        assert config.temperature == 300.0
        assert config.pressure == 1.0  # Default pressure


class TestCoSolventVolumeValidation:
    """Test co-solvent volume fraction / concentration validation."""

    def test_concentration_cosolvent_does_not_crash_validator(self):
        """Concentration-based co-solvents must not crash validate_volume_fractions.

        Regression: sum() over volume_fraction raised TypeError when any
        co-solvent used concentration (volume_fraction=None).
        """
        from polyzymd.config.schema import CoSolventSpec, SolventConfig

        config = SolventConfig(co_solvents=[CoSolventSpec(name="urea", concentration=2.0)])
        assert len(config.co_solvents) == 1

    def test_mixed_volume_fraction_and_concentration_cosolvents(self):
        """Mixed volume_fraction + concentration co-solvents should validate."""
        from polyzymd.config.schema import CoSolventSpec, SolventConfig

        config = SolventConfig(
            co_solvents=[
                CoSolventSpec(name="dmso", volume_fraction=0.3),
                CoSolventSpec(name="urea", concentration=2.0),
            ]
        )
        assert len(config.co_solvents) == 2

    def test_volume_fractions_exceeding_one_rejected(self):
        """Total volume fractions >= 1.0 should still be rejected."""
        from pydantic import ValidationError

        from polyzymd.config.schema import CoSolventSpec, SolventConfig

        with pytest.raises(ValidationError, match="volume fraction must be < 1.0"):
            SolventConfig(
                co_solvents=[
                    CoSolventSpec(name="dmso", volume_fraction=0.6),
                    CoSolventSpec(name="ethanol", volume_fraction=0.5),
                ]
            )

    def test_volume_fractions_below_one_accepted(self):
        """Total volume fractions < 1.0 should pass validation."""
        from polyzymd.config.schema import CoSolventSpec, SolventConfig

        config = SolventConfig(
            co_solvents=[
                CoSolventSpec(name="dmso", volume_fraction=0.3),
                CoSolventSpec(name="ethanol", volume_fraction=0.2),
            ]
        )
        assert len(config.co_solvents) == 2


class TestSimulationPhasesConfig:
    """Test staged equilibration requirements."""

    def test_requires_equilibration_stages(self):
        from pydantic import ValidationError

        from polyzymd.config.schema import SimulationPhaseConfig, SimulationPhasesConfig

        production = SimulationPhaseConfig(
            ensemble="NPT",
            duration=1.0,
            samples=10,
            time_step=2.0,
        )

        with pytest.raises(ValidationError, match="requires 'equilibration_stages'"):
            SimulationPhasesConfig(production=production)

    def test_rejects_empty_equilibration_stages(self):
        from pydantic import ValidationError

        from polyzymd.config.schema import SimulationPhaseConfig, SimulationPhasesConfig

        production = SimulationPhaseConfig(
            ensemble="NPT",
            duration=1.0,
            samples=10,
            time_step=2.0,
        )

        with pytest.raises(ValidationError, match="must contain at least one stage"):
            SimulationPhasesConfig(equilibration_stages=[], production=production)


# ---------------------------------------------------------------------------
# B13 – statepoint export handles concentration-based co-solvents
# ---------------------------------------------------------------------------


class TestStatepointCoSolventExport:
    """to_statepoint must not crash when co-solvents use concentration instead
    of volume_fraction."""

    def _make_config(self, *, volume_fraction=None, concentration=None):
        """Build a minimal SimulationConfig with one co-solvent."""
        from unittest.mock import MagicMock

        config = MagicMock()
        config.enzyme.name = "CALB"
        config.thermodynamics.temperature = 310.0
        config.substrate = None
        config.polymers = None

        cs = MagicMock()
        cs.name = "urea"
        cs.volume_fraction = volume_fraction
        cs.concentration = concentration
        config.solvent.co_solvents = [cs]
        return config

    def test_volume_fraction_exported_as_fraction_key(self):
        """volume_fraction co-solvent produces a _fraction key."""
        from polyzymd.config.schema import SimulationConfig

        _to_statepoint = SimulationConfig.__dict__["to_signac_statepoint"]
        cfg = self._make_config(volume_fraction=0.3)
        sp = _to_statepoint(cfg)
        assert sp["cosolvent_urea_fraction"] == 0.3
        assert "cosolvent_urea_molarity" not in sp

    def test_concentration_exported_as_molarity_key(self):
        """concentration co-solvent produces a _molarity key, not _fraction."""
        from polyzymd.config.schema import SimulationConfig

        _to_statepoint = SimulationConfig.__dict__["to_signac_statepoint"]
        cfg = self._make_config(concentration=2.0)
        sp = _to_statepoint(cfg)
        assert sp["cosolvent_urea_molarity"] == 2.0
        assert "cosolvent_urea_fraction" not in sp

    def test_no_crash_with_none_volume_fraction(self):
        """Previously crashed with TypeError when volume_fraction was None."""
        from polyzymd.config.schema import SimulationConfig

        _to_statepoint = SimulationConfig.__dict__["to_signac_statepoint"]
        cfg = self._make_config(concentration=1.5)
        # Should not raise
        sp = _to_statepoint(cfg)
        assert "cosolvent_urea_molarity" in sp


class TestEngineConfig:
    """Tests for engine configuration fields."""

    def test_default_engine_is_openmm(self, minimal_config_data):
        """Default engine should be openmm for backward compatibility."""
        config = SimulationConfig(**minimal_config_data)
        assert config.engine == "openmm"

    def test_gromacs_engine_parses(self, minimal_config_data):
        """Setting engine to gromacs should work."""
        minimal_config_data["engine"] = "gromacs"
        config = SimulationConfig(**minimal_config_data)
        assert config.engine == "gromacs"

    def test_invalid_engine_rejected(self, minimal_config_data):
        """Invalid engine names should be rejected."""
        minimal_config_data["engine"] = "lammps"
        with pytest.raises(ValidationError):
            SimulationConfig(**minimal_config_data)

    def test_openmm_engine_config_defaults(self, minimal_config_data):
        """OpenMM engine config should have sensible defaults."""
        config = SimulationConfig(**minimal_config_data)
        assert config.openmm.platform == "CUDA"
        assert config.openmm.precision == "mixed"

    def test_gromacs_engine_config_defaults(self, minimal_config_data):
        """GROMACS engine config should have sensible defaults."""
        config = SimulationConfig(**minimal_config_data)
        assert config.gromacs.gmx_binary is None
        assert config.gromacs.grompp_flags == "-maxwarn 1"
        assert config.gromacs.mdrun_flags_equilibration is None
        assert config.gromacs.mdrun_flags_production is None
        assert config.gromacs.command_prefix is None
        assert config.gromacs.mpi_launcher_flags == ""
        assert config.gromacs.env_exports == {}
        assert config.gromacs.setup_commands == []
        assert config.gromacs.ntmpi == 1
        assert config.gromacs.slurm_ntasks is None
        assert config.gromacs.ntomp == 8
        assert config.gromacs.gpu is False
        assert config.gromacs.memory == "16G"

    def test_gromacs_slurm_ntasks_default_none(self, minimal_config_data):
        """GROMACS slurm_ntasks should default to None."""
        config = SimulationConfig(**minimal_config_data)
        assert config.gromacs.slurm_ntasks is None

    def test_gromacs_slurm_ntasks_override(self, minimal_config_data):
        """GROMACS slurm_ntasks should accept explicit override values."""
        minimal_config_data["gromacs"] = {"slurm_ntasks": 16}
        config = SimulationConfig(**minimal_config_data)
        assert config.gromacs.slurm_ntasks == 16

    def test_gromacs_slurm_ntasks_zero_rejected(self, minimal_config_data):
        """GROMACS slurm_ntasks must be a positive integer."""
        minimal_config_data["gromacs"] = {"slurm_ntasks": 0}
        with pytest.raises(ValidationError):
            SimulationConfig(**minimal_config_data)

    def test_gromacs_engine_config_gpu_enabled(self, minimal_config_data):
        """GROMACS config should accept explicit GPU settings."""
        minimal_config_data["gromacs"] = {"gpu": True, "ntomp": 4}
        config = SimulationConfig(**minimal_config_data)
        assert config.gromacs.gpu is True
        assert config.gromacs.ntomp == 4

    def test_gromacs_ntomp_validation(self, minimal_config_data):
        """GROMACS ntomp must be positive."""
        minimal_config_data["gromacs"] = {"ntomp": 0}
        with pytest.raises(ValidationError):
            SimulationConfig(**minimal_config_data)

    def test_gromacs_engine_config_custom(self, minimal_config_data):
        """Custom GROMACS settings should be parsed."""
        minimal_config_data["engine"] = "gromacs"
        minimal_config_data["gromacs"] = {
            "gmx_binary": "gmx_mpi",
            "mdrun_flags": "-ntmpi 1 -ntomp 8",
            "mdrun_flags_equilibration": "-ntomp 4",
            "mdrun_flags_production": "-ntomp 8 -plumed plumed_setup.dat",
            "command_prefix": "singularity exec --rocm --bind /scratch /path/to/gromacs.sif",
            "mpi_launcher_flags": "-genv I_MPI_FABRICS shm:tcp",
            "env_exports": {
                "GMX_GPU_DD_COMMS": "true",
                "GMX_FORCE_UPDATE_DEFAULT_GPU": "true",
            },
            "setup_commands": [
                "source /opt/gromacs/bin/GMXRC",
                "export PATH=$PATH:/opt/plumed/bin",
            ],
        }
        config = SimulationConfig(**minimal_config_data)
        assert config.gromacs.gmx_binary == "gmx_mpi"
        assert config.gromacs.mdrun_flags == "-ntmpi 1 -ntomp 8"
        assert config.gromacs.mdrun_flags_equilibration == "-ntomp 4"
        assert config.gromacs.mdrun_flags_production == "-ntomp 8 -plumed plumed_setup.dat"
        assert (
            config.gromacs.command_prefix
            == "singularity exec --rocm --bind /scratch /path/to/gromacs.sif"
        )
        assert config.gromacs.mpi_launcher_flags == "-genv I_MPI_FABRICS shm:tcp"
        assert config.gromacs.env_exports["GMX_GPU_DD_COMMS"] == "true"
        assert config.gromacs.setup_commands[0] == "source /opt/gromacs/bin/GMXRC"

    def test_old_config_without_engine_field(self, minimal_config_data):
        """Old configs without engine field should default to openmm."""
        minimal_config_data.pop("engine", None)
        config = SimulationConfig(**minimal_config_data)
        assert config.engine == "openmm"


class TestGromacsEngineConfigWarnings:
    """Tests for GROMACS config warning validators."""

    def test_gpu_ntmpi_warning(self, minimal_config_data, caplog):
        """gpu=True + ntmpi>1 should log a warning."""
        import logging

        minimal_config_data["gromacs"] = {"gpu": True, "ntmpi": 4, "gpus": 2}
        with caplog.at_level(logging.WARNING):
            config = SimulationConfig(**minimal_config_data)

        assert any("thread-MPI" in r.message for r in caplog.records)
        assert config.gromacs.gpu is True
        assert config.gromacs.ntmpi == 4

    def test_no_warning_gpu_ntmpi_1(self, minimal_config_data, caplog):
        """gpu=True + ntmpi=1 should not warn."""
        import logging

        minimal_config_data["gromacs"] = {"gpu": True, "ntmpi": 1}
        with caplog.at_level(logging.WARNING):
            SimulationConfig(**minimal_config_data)

        assert not any("thread-MPI" in r.message for r in caplog.records)

    def test_no_warning_cpu_ntmpi_high(self, minimal_config_data, caplog):
        """gpu=False + ntmpi>1 should not warn."""
        import logging

        minimal_config_data["gromacs"] = {"gpu": False, "ntmpi": 4}
        with caplog.at_level(logging.WARNING):
            SimulationConfig(**minimal_config_data)

        assert not any("thread-MPI" in r.message for r in caplog.records)

    def test_mdrun_flags_ntmpi_conflict_warns(self, minimal_config_data, caplog):
        """mdrun_flags with -ntmpi != config ntmpi should warn."""
        import logging

        minimal_config_data["gromacs"] = {"ntmpi": 1, "mdrun_flags": "-ntmpi 4"}
        with caplog.at_level(logging.WARNING):
            SimulationConfig(**minimal_config_data)

        assert any("ntmpi" in r.message and "override" in r.message for r in caplog.records)

    def test_mdrun_flags_ntomp_conflict_warns(self, minimal_config_data, caplog):
        """mdrun_flags with -ntomp != config ntomp should warn."""
        import logging

        minimal_config_data["gromacs"] = {"ntomp": 8, "mdrun_flags": "-ntomp 16"}
        with caplog.at_level(logging.WARNING):
            SimulationConfig(**minimal_config_data)

        assert any("ntomp" in r.message and "override" in r.message for r in caplog.records)

    def test_mdrun_flags_matching_no_warning(self, minimal_config_data, caplog):
        """mdrun_flags matching config fields should not warn."""
        import logging

        minimal_config_data["gromacs"] = {
            "ntmpi": 1,
            "ntomp": 12,
            "mdrun_flags": "-ntmpi 1 -ntomp 12 -nb gpu",
        }
        with caplog.at_level(logging.WARNING):
            SimulationConfig(**minimal_config_data)

        assert not any("override" in r.message for r in caplog.records)

    def test_mdrun_flags_malformed_no_crash(self, minimal_config_data):
        """Malformed mdrun_flags should not crash validation."""
        minimal_config_data["gromacs"] = {"mdrun_flags": "-ntmpi 'unclosed"}
        config = SimulationConfig(**minimal_config_data)
        assert config.gromacs.mdrun_flags == "-ntmpi 'unclosed"

    def test_mdrun_flags_empty_no_warning(self, minimal_config_data, caplog):
        """Empty mdrun_flags should not produce warnings."""
        import logging

        minimal_config_data["gromacs"] = {"mdrun_flags": ""}
        with caplog.at_level(logging.WARNING):
            SimulationConfig(**minimal_config_data)

        assert not any("override" in r.message for r in caplog.records)
