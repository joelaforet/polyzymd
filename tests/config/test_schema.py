"""Test that all public modules can be imported."""

from pathlib import Path

import pytest


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
