"""Test that all public modules can be imported."""

import warnings
from pathlib import Path

import pytest
from pydantic import ValidationError

from polyzymd.config.schema import (
    BoxConfig,
    BoxShape,
    BuildScope,
    ConjugationConfig,
    SimulationConfig,
)


@pytest.fixture
def minimal_config_data():
    """Provide minimal valid SimulationConfig input data."""
    return {
        "name": "test_simulation",
        "engine": "openmm",
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
                "checkpoint_interval": 60.0,
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
        from polyzymd.config import BuildConfig, BuildScope, SimulationConfig

        assert SimulationConfig is not None
        assert BuildConfig is not None
        assert BuildScope is not None

    def test_import_conjugation_public_api_without_rdkit(self, monkeypatch):
        """Importing the conjugation public API should not require RDKit."""
        import builtins
        import importlib

        import polyzymd.builders  # noqa: F401 - isolate this test to conjugation imports

        real_import = builtins.__import__

        def block_rdkit_import(name, globals=None, locals=None, fromlist=(), level=0):
            """Block RDKit imports while allowing all other imports."""
            if name == "rdkit" or name.startswith("rdkit."):
                raise ImportError("RDKit intentionally blocked for import test")
            return real_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", block_rdkit_import)
        module = importlib.import_module("polyzymd.builders.conjugation")

        assert module.build_conjugate_from_config is not None

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


class TestConjugationPlacementConfig:
    """Test the public Packmol timeout policy."""

    def test_packmol_timeout_defaults_to_fifteen_minutes(self):
        """Conjugation placement should be bounded by default."""
        config = ConjugationConfig()

        assert config.placement.timeout_seconds == pytest.approx(900.0)

    def test_packmol_timeout_accepts_seconds_or_null(self):
        """Users should be able to shorten or explicitly disable the timeout."""
        bounded = ConjugationConfig(placement={"timeout_seconds": 30.0})
        unbounded = ConjugationConfig(placement={"timeout_seconds": None})

        assert bounded.placement.timeout_seconds == pytest.approx(30.0)
        assert unbounded.placement.timeout_seconds is None

    def test_packmol_timeout_rejects_nonpositive_values(self):
        """A configured timeout must be positive."""
        with pytest.raises(ValidationError, match="timeout_seconds"):
            ConjugationConfig(placement={"timeout_seconds": 0.0})


class TestBuildConfig:
    """Test the public construction endpoint contract."""

    def test_build_scope_defaults_to_system(self, minimal_config_data):
        """Existing configurations should continue to build complete systems."""
        config = SimulationConfig(**minimal_config_data)

        assert config.build.scope is BuildScope.SYSTEM
        assert config.model_dump(mode="json")["build"] == {"scope": "system"}

    @pytest.mark.parametrize(
        ("value", "expected"),
        [
            ("structure", BuildScope.STRUCTURE),
            ("solute", BuildScope.SOLUTE),
            ("system", BuildScope.SYSTEM),
        ],
    )
    def test_build_scope_accepts_each_public_endpoint(self, minimal_config_data, value, expected):
        """Each documented endpoint should validate from configuration data."""
        minimal_config_data["build"] = {"scope": value}

        config = SimulationConfig(**minimal_config_data)

        assert config.build.scope is expected

    def test_build_scope_rejects_unknown_endpoint(self, minimal_config_data):
        """Typos should fail before construction starts."""
        minimal_config_data["build"] = {"scope": "parameterized"}

        with pytest.raises(ValidationError, match="build.scope"):
            SimulationConfig(**minimal_config_data)

    def test_build_config_rejects_unknown_fields(self, minimal_config_data):
        """Build-only options should not be silently ignored."""
        minimal_config_data["build"] = {"scope": "solute", "through": "parameterization"}

        with pytest.raises(ValidationError, match="build.through"):
            SimulationConfig(**minimal_config_data)


class TestBoxConfig:
    """Test canonical rectangular-box configuration."""

    def test_orthorhombic_is_canonical(self):
        """The explicit rectangular name should validate without migration."""
        config = BoxConfig(shape="orthorhombic")

        assert config.shape is BoxShape.ORTHORHOMBIC
        assert config.model_dump(mode="json")["shape"] == "orthorhombic"

    def test_cube_migrates_with_one_deprecation_warning(self):
        """Legacy cube input should retain behavior while naming it accurately."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            config = BoxConfig(shape="cube")

        deprecations = [item for item in caught if item.category is DeprecationWarning]
        assert len(deprecations) == 1
        assert "orthorhombic" in str(deprecations[0].message)
        assert config.shape is BoxShape.ORTHORHOMBIC


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

    @pytest.mark.parametrize("length", [1, 2])
    def test_dynamic_polymer_length_requires_non_empty_middle_sequence(self, length: int):
        """Dynamic polymer configuration should reject lengths below three."""
        from polyzymd.config.schema import MonomerSpec, PolymerConfig, ReactionConfig

        reactions = ReactionConfig(
            initiation=Path("initiation.rxn"),
            polymerization=Path("polymerization.rxn"),
            termination=Path("termination.rxn"),
        )

        with pytest.raises(ValidationError, match="length.*>= 3"):
            PolymerConfig(
                type_prefix="TEST",
                generation_mode="dynamic",
                monomers=[
                    MonomerSpec(label="A", probability=1.0, name="EGMA", smiles="C=C"),
                ],
                length=length,
                count=1,
                reactions=reactions,
            )

    def test_cached_polymer_length_one_is_accepted(self):
        """Cached polymer configuration should allow single-monomer SDF inputs."""
        from polyzymd.config.schema import MonomerSpec, PolymerConfig

        config = PolymerConfig(
            type_prefix="TEST",
            generation_mode="cached",
            monomers=[MonomerSpec(label="A", probability=1.0, name="EGMA")],
            length=1,
            count=1,
            sdf_directory=Path("/tmp/test"),
        )

        assert config.length == 1

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

    def test_ion_config_defaults_to_neutralization_without_salt(self):
        """IonConfig should neutralize without adding bulk salt by default."""
        from polyzymd.config.schema import IonConfig

        config = IonConfig()
        assert config.neutralize is True
        assert config.nacl_concentration == 0.0

    def test_ion_config_nacl_description_documents_target_concentration(self):
        """IonConfig should document NaCl-equivalent target concentration."""
        from polyzymd.config.schema import IonConfig

        description = IonConfig.model_fields["nacl_concentration"].description

        assert description is not None
        assert "NaCl-equivalent" in description
        assert "target concentration" in description
        assert "pool" not in description
        assert "Additional" not in description

    @pytest.mark.parametrize(
        ("field", "value"),
        [
            ("local_junction_patching", False),
            ("patch_radius_bonds", 2),
            ("preserve_total_charge", False),
            ("integer_tolerance", 1.0e-3),
        ],
    )
    def test_conjugation_rejects_removed_product_charge_knobs(self, field, value):
        """Conjugation config should expose no product charge tuning knobs."""
        from polyzymd.config.schema import ConjugationConfig

        with pytest.raises(ValidationError, match="charge"):
            ConjugationConfig(charge={field: value})


class TestConjugationReactionSmartsConfig:
    """Tests for generic atom-mapped reaction metadata in conjugation config."""

    def test_mechanism_accepts_reaction_smarts_with_atom_roles(self):
        """Mechanism config should validate generic atom-map role metadata."""
        from polyzymd.config.schema import ConjugationMechanismConfig

        mechanism = ConjugationMechanismConfig(
            name="generic_amide",
            reaction_smarts="[N:1].[C:2]>>[N:1]-[C:2]",
            atom_roles=[
                {"map_number": 1, "participant": "site", "role": "linking", "label": "donor"},
                {"map_number": 2, "participant": "moiety", "role": "linking"},
            ],
        )

        assert mechanism.reaction_smarts == "[N:1].[C:2]>>[N:1]-[C:2]"
        assert mechanism.atom_roles[0].map_number == 1
        assert mechanism.atom_roles[0].participant == "site"

    def test_reaction_smarts_requires_reaction_separator(self):
        """Invalid reaction SMARTS should fail before workflow execution."""
        from polyzymd.config.schema import ConjugationMechanismConfig

        with pytest.raises(ValidationError, match="reaction_smarts must contain"):
            ConjugationMechanismConfig(reaction_smarts="[N:1].[C:2]")

    def test_explicit_linkage_accepts_reaction_smarts_as_metadata(self):
        """PDB explicit linkage config may carry mapped SMARTS for generic preflights."""
        from polyzymd.config.schema import ConjugationAttachmentConfig

        attachment = ConjugationAttachmentConfig.model_validate(
            {
                "name": "mapped-explicit-linkage",
                "site": {
                    "chain_id": "A",
                    "residue_name": "AAA",
                    "residue_number": 10,
                    "atom_name": "N1",
                },
                "moiety": {
                    "name": "modifier",
                    "input_path": "modifier.pdb",
                    "link_site": {
                        "chain_id": "C",
                        "residue_name": "BBB",
                        "residue_number": 1,
                        "atom_name": "C1",
                    },
                },
                "mechanism": {
                    "name": "explicit_linkage",
                    "reaction_smarts": "[N:1].[C:2]>>[N:1]-[C:2]",
                    "atom_roles": [
                        {"map_number": 1, "participant": "site", "role": "linking"},
                        {"map_number": 2, "participant": "moiety", "role": "linking"},
                    ],
                    "product_residues": {"site": "AAA", "moiety": "BBB"},
                },
            }
        )

        assert attachment.mechanism.reaction_smarts == "[N:1].[C:2]>>[N:1]-[C:2]"
        assert len(attachment.mechanism.atom_roles) == 2


class TestCoSolventCompositionValidation:
    """Test co-solvent mole fraction and concentration validation."""

    def test_concentration_cosolvent_does_not_crash_validator(self):
        """Concentration-based co-solvents should validate without mole fractions."""
        from polyzymd.config.schema import CoSolventSpec, SolventConfig

        config = SolventConfig(co_solvents=[CoSolventSpec(name="urea", concentration=2.0)])
        assert len(config.co_solvents) == 1

    def test_mixed_mole_fraction_and_concentration_cosolvents(self):
        """Mixed mole_fraction and concentration co-solvents should validate."""
        from polyzymd.config.schema import CoSolventSpec, SolventConfig

        config = SolventConfig(
            co_solvents=[
                CoSolventSpec(name="dmso", mole_fraction=0.3),
                CoSolventSpec(name="urea", concentration=2.0),
            ]
        )
        assert len(config.co_solvents) == 2

    def test_mole_fractions_exceeding_one_rejected(self):
        """Total mole fractions >= 1.0 should be rejected."""
        from pydantic import ValidationError

        from polyzymd.config.schema import CoSolventSpec, SolventConfig

        with pytest.raises(ValidationError, match="mole fraction must be < 1.0"):
            SolventConfig(
                co_solvents=[
                    CoSolventSpec(name="dmso", mole_fraction=0.6),
                    CoSolventSpec(name="ethanol", mole_fraction=0.5),
                ]
            )

    def test_mole_fractions_below_one_accepted(self):
        """Total mole fractions < 1.0 should pass validation."""
        from polyzymd.config.schema import CoSolventSpec, SolventConfig

        config = SolventConfig(
            co_solvents=[
                CoSolventSpec(name="dmso", mole_fraction=0.3),
                CoSolventSpec(name="ethanol", mole_fraction=0.2),
            ]
        )
        assert len(config.co_solvents) == 2

    def test_exactly_one_composition_method_required(self):
        """Each co-solvent should have exactly one composition method."""
        from polyzymd.config.schema import CoSolventSpec

        with pytest.raises(ValidationError, match="Must specify either 'mole_fraction'"):
            CoSolventSpec(name="dmso")

        with pytest.raises(ValidationError, match="Cannot specify both 'mole_fraction'"):
            CoSolventSpec(name="dmso", mole_fraction=0.1, concentration=1.0)


class TestRunDirectoryNaming:
    """Test centralized run-directory naming helpers."""

    def test_output_config_format_directory_name_accepts_solvent_tokens(self):
        """Output formatter should support solvent template placeholders directly."""
        from polyzymd.config.schema import OutputConfig

        output = OutputConfig(
            naming_template=(
                "{enzyme}_{primary_solvent}_{cosolvent_composition}_{solvent_composition}_"
                "r{replicate}"
            ),
        )

        assert (
            output.format_directory_name(
                enzyme="LipA",
                substrate="apo",
                polymer_type="none",
                temperature=310.0,
                replicate=2,
                primary_solvent="water_tip3p",
                cosolvent_composition="dmso_30molpct",
                solvent_composition="water_tip3p_dmso_30molpct",
            )
            == "LipA_water_tip3p_dmso_30molpct_water_tip3p_dmso_30molpct_r2"
        )
        assert (
            output.format_directory_name(
                enzyme="LipA",
                substrate="apo",
                polymer_type="none",
                temperature=310.0,
                replicate=2,
                primary_solvent="water_tip3p",
                cosolvent_composition="dmso_30molpct",
            )
            == "LipA_water_tip3p_dmso_30molpct_water_tip3p_dmso_30molpct_r2"
        )

    def test_format_run_directory_name_uses_output_formatter(self, minimal_config_data):
        """Simulation formatter should route centralized values through OutputConfig."""
        from polyzymd.config.schema import SimulationConfig

        config = SimulationConfig(**minimal_config_data)
        calls = []
        original = config.output.format_directory_name

        def spy_format_directory_name(**values):
            calls.append(values)
            return original(**values)

        object.__setattr__(config.output, "format_directory_name", spy_format_directory_name)

        assert config.format_run_directory_name(4) == "TestEnzyme_apo_none_1ns_300K_run4"
        assert calls == [config._run_directory_template_values(4)]

    def test_solvent_token_formatting(self, minimal_config_data):
        """Solvent placeholders should render safe composition tokens."""
        minimal_config_data["solvent"] = {
            "primary": {"type": "water", "model": "tip3p"},
            "co_solvents": [
                {"name": "urea", "concentration": 2.5},
                {
                    "name": "tert butanol",
                    "smiles": "CC(C)(C)O",
                    "density": 0.78,
                    "mole_fraction": 0.125,
                },
                {"name": "dmso", "mole_fraction": 0.30},
            ],
        }
        minimal_config_data["output"] = {
            "naming_template": "{primary_solvent}_{cosolvent_composition}_{solvent_composition}"
        }

        config = SimulationConfig(**minimal_config_data)

        assert config.format_run_directory_name() == (
            "water_tip3p_dmso_30molpct_tert_butanol_12p5molpct_urea_2p5M_"
            "water_tip3p_dmso_30molpct_tert_butanol_12p5molpct_urea_2p5M"
        )

    def test_high_mole_fraction_solvent_token_does_not_crash(self, minimal_config_data):
        """Student DMSO case should format using mole-fraction tokens."""
        minimal_config_data["solvent"] = {
            "primary": {"type": "water", "model": "tip3p"},
            "co_solvents": [{"name": "dmso", "mole_fraction": 0.99}],
        }
        minimal_config_data["output"] = {"naming_template": "{cosolvent_composition}"}

        config = SimulationConfig(**minimal_config_data)

        assert config.format_run_directory_name() == "dmso_99molpct"
        assert config._cosolvent_composition_token() == "dmso_99molpct"

    def test_absent_cosolvent_uses_none_and_primary_only_composition(self, minimal_config_data):
        """No co-solvents should produce none and primary-only solvent composition."""
        minimal_config_data["output"] = {
            "naming_template": "{primary_solvent}_{cosolvent_composition}_{solvent_composition}"
        }

        config = SimulationConfig(**minimal_config_data)

        assert config.format_run_directory_name() == "water_tip3p_none_water_tip3p"

    def test_non_water_primary_solvent_is_sanitized(self, minimal_config_data):
        """Non-water primary solvent names should normalize unsafe punctuation."""
        minimal_config_data["solvent"] = {"primary": {"type": "ACN / MeOH"}}
        minimal_config_data["output"] = {"naming_template": "{primary_solvent}"}

        config = SimulationConfig(**minimal_config_data)

        assert config.format_run_directory_name() == "acn_meoh"

    def test_format_run_directory_name_and_working_directory_match(
        self, minimal_config_data, tmp_path
    ):
        """The public formatter should match get_working_directory().name."""
        minimal_config_data["output"] = {
            "scratch_directory": str(tmp_path),
            "naming_template": "{enzyme}_{temperature}K_run{replicate}_{solvent_composition}",
        }
        config = SimulationConfig(**minimal_config_data)

        run_name = config.format_run_directory_name(3)

        assert run_name == "TestEnzyme_300K_run3_water_tip3p"
        assert config.get_working_directory(3).name == run_name

    def test_discover_replicate_dirs_with_solvent_placeholder_and_middle_replicate(
        self, minimal_config_data, tmp_path
    ):
        """Discovery should work when solvent placeholders are present after replicate."""
        minimal_config_data["output"] = {
            "scratch_directory": str(tmp_path),
            "naming_template": "run{replicate}_{solvent_composition}_{enzyme}",
        }
        config = SimulationConfig(**minimal_config_data)
        first = tmp_path / config.format_run_directory_name(1)
        third = tmp_path / config.format_run_directory_name(3)
        first.mkdir()
        third.mkdir()
        (tmp_path / "runX_water_tip3p_TestEnzyme").mkdir()

        assert config.discover_replicate_dirs() == [(1, first), (3, third)]

    def test_discover_replicate_dirs_requires_replicate_placeholder(
        self, minimal_config_data, tmp_path
    ):
        """Discovery should fail clearly when the template has no replicate token."""
        minimal_config_data["output"] = {
            "scratch_directory": str(tmp_path),
            "naming_template": "{enzyme}_{solvent_composition}",
        }
        config = SimulationConfig(**minimal_config_data)

        with pytest.raises(ValueError, match=r"does not include the \{replicate\} placeholder"):
            config.discover_replicate_dirs()


class TestSimulationPhasesConfig:
    """Test staged equilibration requirements."""

    def test_ignores_deprecated_report_interval(self, caplog):
        """Legacy intervals are ignored so samples controls trajectory cadence."""
        from polyzymd.config.schema import SimulationPhaseConfig

        phase = SimulationPhaseConfig(
            ensemble="NPT",
            duration=1.0,
            samples=10,
            report_interval=50000,
            checkpoint_interval=60.0,
        )

        assert "report_interval" not in phase.model_dump()
        assert "Ignoring deprecated" in caplog.text

    def test_requires_equilibration_stages(self):
        from pydantic import ValidationError

        from polyzymd.config.schema import SimulationPhaseConfig, SimulationPhasesConfig

        production = SimulationPhaseConfig(
            ensemble="NPT",
            duration=1.0,
            samples=10,
            checkpoint_interval=60.0,
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
            checkpoint_interval=60.0,
            time_step=2.0,
        )

        with pytest.raises(ValidationError, match="must contain at least one stage"):
            SimulationPhasesConfig(equilibration_stages=[], production=production)

    def test_rejects_legacy_segments_field(self):
        """Legacy simulation_phases.segments should be rejected."""
        from pydantic import ValidationError

        from polyzymd.config.schema import EquilibrationStageConfig, SimulationPhasesConfig

        with pytest.raises(ValidationError, match="Extra inputs are not permitted"):
            SimulationPhasesConfig(
                equilibration_stages=[
                    EquilibrationStageConfig(
                        name="heating",
                        ensemble="NVT",
                        duration=1.0,
                        samples=10,
                        temperature=300.0,
                    )
                ],
                production={
                    "ensemble": "NPT",
                    "duration": 1.0,
                    "samples": 10,
                    "checkpoint_interval": 60.0,
                },
                segments=[],
            )


class TestEquilibrationTemperatureRamp:
    """Validate increment-based temperature ramps and derived durations."""

    @staticmethod
    def _stage(**overrides):
        from polyzymd.config.schema import EquilibrationStageConfig

        values = {
            "name": "heating",
            "temperature_start": 60.0,
            "temperature_end": 300.0,
            "temperature_increment": 1.0,
            "temperature_interval_steps": 600,
            "time_step": 2.0,
        }
        values.update(overrides)
        return EquilibrationStageConfig(**values)

    def test_derives_duration_without_user_duration(self):
        stage = self._stage()

        assert stage.resolved_duration == pytest.approx(0.288)
        assert stage.temperature_ramp_updates == 240

    def test_shortens_final_increment_to_hit_endpoint(self):
        stage = self._stage(
            temperature_start=100.0,
            temperature_end=112.0,
            temperature_increment=5.0,
            temperature_interval_steps=100,
        )

        assert stage.temperature_ramp_updates == 3
        assert stage.resolved_duration == pytest.approx(0.0006)
        assert stage.temperature_at_step(200, 300) == pytest.approx(110.0)
        assert stage.temperature_at_step(300, 300) == pytest.approx(112.0)

    def test_fractional_increment_does_not_add_spurious_update(self):
        stage = self._stage(
            temperature_start=60.0,
            temperature_end=60.7,
            temperature_increment=0.1,
            temperature_interval_steps=100,
        )

        assert stage.temperature_ramp_updates == 7
        assert stage.resolved_duration == pytest.approx(0.0014)

    def test_serialized_ramp_omits_derived_duration_and_round_trips(self):
        from polyzymd.config.schema import EquilibrationStageConfig

        stage = self._stage()
        serialized = stage.model_dump()

        assert "duration" not in serialized
        restored = EquilibrationStageConfig.model_validate(serialized)
        assert restored.resolved_duration == stage.resolved_duration
        assert restored.temperature_ramp_updates == stage.temperature_ramp_updates

    def test_rejects_duration_with_increment_based_ramp(self):
        with pytest.raises(ValidationError, match="Do not specify 'duration'"):
            self._stage(duration=0.2)

    def test_requires_duration_for_constant_temperature(self):
        from polyzymd.config.schema import EquilibrationStageConfig

        with pytest.raises(ValidationError, match="require 'duration'"):
            EquilibrationStageConfig(name="equilibration", temperature=300.0)

    @pytest.mark.parametrize("end", [60.0, 50.0])
    def test_requires_increasing_temperature(self, end):
        with pytest.raises(ValidationError, match="must be <"):
            self._stage(temperature_end=end)

    def test_migrates_legacy_increment_and_interval(self, caplog):
        from polyzymd.config.schema import EquilibrationStageConfig

        with pytest.warns(DeprecationWarning, match="temperature_interval_steps"):
            stage = EquilibrationStageConfig(
                name="heating",
                duration=0.1,
                temperature_start=60.0,
                temperature_end=300.0,
                temperature_increment=1.0,
                temperature_interval=1200.0,
            )

        assert stage.temperature_increment == pytest.approx(1.0)
        assert stage.temperature_interval_steps == 600
        assert stage.resolved_duration == pytest.approx(0.288)
        assert "converted to temperature_interval_steps=600" in caplog.text

    def test_migrates_duration_only_ramp_with_legacy_defaults(self, caplog):
        from polyzymd.config.schema import EquilibrationStageConfig

        with pytest.warns(DeprecationWarning, match="legacy 1200 fs default"):
            stage = EquilibrationStageConfig(
                name="heating",
                duration=0.3,
                temperature_start=60.0,
                temperature_end=300.0,
            )

        assert stage.temperature_increment == pytest.approx(1.0)
        assert stage.temperature_interval_steps == 600
        assert stage.resolved_duration == pytest.approx(0.288)
        assert "legacy 1200 fs default" in caplog.text

    def test_constant_stage_does_not_emit_ramp_deprecation_warning(self):
        import warnings

        from polyzymd.config.schema import EquilibrationStageConfig

        with warnings.catch_warnings():
            warnings.simplefilter("error", DeprecationWarning)
            stage = EquilibrationStageConfig(
                name="equilibration",
                duration=0.1,
                temperature=300.0,
            )

        assert stage.resolved_duration == pytest.approx(0.1)

    def test_rejects_decimal_rate(self):
        from polyzymd.config.schema import EquilibrationStageConfig

        with pytest.raises(ValidationError, match="not supported"):
            EquilibrationStageConfig(
                name="heating",
                temperature_start=60.0,
                temperature_end=300.0,
                temperature_ramp_rate=0.8333333333333334,
            )

    def test_rejects_old_and_new_intervals_together(self):
        from polyzymd.config.schema import EquilibrationStageConfig

        with pytest.raises(ValidationError, match="Cannot specify both"):
            EquilibrationStageConfig(
                name="heating",
                temperature_start=60.0,
                temperature_end=300.0,
                temperature_increment=1.0,
                temperature_interval=1200.0,
                temperature_interval_steps=600,
            )


# ---------------------------------------------------------------------------
# B13 – statepoint export handles concentration-based co-solvents
# ---------------------------------------------------------------------------


class TestStatepointCoSolventExport:
    """to_statepoint must not crash when co-solvents use concentration instead
    of mole_fraction."""

    def _make_config(self, *, mole_fraction=None, concentration=None):
        """Build a minimal SimulationConfig with one co-solvent."""
        from unittest.mock import MagicMock

        config = MagicMock()
        config.enzyme.name = "CALB"
        config.thermodynamics.temperature = 310.0
        config.substrate = None
        config.polymers = None

        cs = MagicMock()
        cs.name = "urea"
        cs.mole_fraction = mole_fraction
        cs.concentration = concentration
        config.solvent.co_solvents = [cs]
        return config

    def test_mole_fraction_exported_as_mole_fraction_key(self):
        """mole_fraction co-solvent produces a _mole_fraction key."""
        from polyzymd.config.schema import SimulationConfig

        _to_statepoint = SimulationConfig.__dict__["to_signac_statepoint"]
        cfg = self._make_config(mole_fraction=0.3)
        sp = _to_statepoint(cfg)
        assert sp["cosolvent_urea_mole_fraction"] == 0.3
        assert "cosolvent_urea_molarity" not in sp

    def test_concentration_exported_as_molarity_key(self):
        """concentration co-solvent produces a _molarity key, not _fraction."""
        from polyzymd.config.schema import SimulationConfig

        _to_statepoint = SimulationConfig.__dict__["to_signac_statepoint"]
        cfg = self._make_config(concentration=2.0)
        sp = _to_statepoint(cfg)
        assert sp["cosolvent_urea_molarity"] == 2.0
        assert "cosolvent_urea_mole_fraction" not in sp

    def test_no_crash_with_none_mole_fraction(self):
        """Concentration statepoints should allow no mole_fraction."""
        from polyzymd.config.schema import SimulationConfig

        _to_statepoint = SimulationConfig.__dict__["to_signac_statepoint"]
        cfg = self._make_config(concentration=1.5)
        # Should not raise
        sp = _to_statepoint(cfg)
        assert "cosolvent_urea_molarity" in sp


class TestOutputConfig:
    """Tests for output path configuration."""

    def test_defaults_use_projects_as_effective_scratch(self):
        """Default output settings should write scratch data under projects."""
        from polyzymd.config.schema import OutputConfig

        config = OutputConfig()

        assert config.projects_directory == Path(".")
        assert config.scratch_directory is None
        assert config.effective_scratch_directory == Path(".")

    def test_base_directory_field_is_rejected(self):
        """OutputConfig should reject the removed base_directory field."""
        from polyzymd.config.schema import OutputConfig

        with pytest.raises(ValidationError, match="base_directory"):
            OutputConfig(base_directory="/tmp/old-output")

    def test_output_base_directory_yaml_data_is_rejected(self, minimal_config_data):
        """SimulationConfig data should reject stale output.base_directory."""
        minimal_config_data["output"] = {"base_directory": "/tmp/old-output"}

        with pytest.raises(ValidationError, match="base_directory"):
            SimulationConfig(**minimal_config_data)


class TestEngineConfig:
    """Tests for engine configuration fields."""

    def test_explicit_openmm_engine_parses(self, minimal_config_data):
        """Explicit OpenMM engine should be accepted."""
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

    def test_missing_engine_field_is_rejected(self, minimal_config_data):
        """Configs without engine field should fail validation."""
        minimal_config_data.pop("engine", None)
        with pytest.raises(ValidationError):
            SimulationConfig(**minimal_config_data)


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
