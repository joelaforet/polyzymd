"""Tests for conjugation polymer recipes and Polymerist smoke generation."""

from __future__ import annotations

import sys
from pathlib import Path
from types import ModuleType

import pytest
from pydantic import ValidationError

from polyzymd.builders.conjugation.polymer_recipe import (
    PolymerRecipe,
    _polymerist_to_explicit_h_rdkit_mol,
    _write_rdkit_sdf_sidecar,
    generate_polymerist_smoke_polymer,
    sbma_egpma_nhs_recipe,
    sbma_nhs_egpma_acb_recipe,
)
from polyzymd.config.loader import load_config
from polyzymd.config.schema import SimulationConfig


def _minimal_simulation_config_data() -> dict:
    """Build a minimal valid simulation config dictionary."""
    return {
        "name": "conjugation-polymer-recipe",
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


def _recipe_dict() -> dict:
    """Return a realistic SBMA/EGPMA/NHS polymer recipe dictionary."""
    return sbma_egpma_nhs_recipe(length=9, seed=7).model_dump(mode="json")


def test_recipe_validates_probabilities_and_pdb_safe_residue_names():
    """Recipes should reject invalid probabilities and residue names."""
    recipe_data = _recipe_dict()
    recipe_data["monomers"][0]["probability"] = 0.5

    with pytest.raises(ValidationError, match="probabilities must sum"):
        PolymerRecipe.model_validate(recipe_data)

    recipe_data = _recipe_dict()
    recipe_data["monomers"][0]["residue_name"] = "SBMA"

    with pytest.raises(ValidationError, match="at most 3 characters"):
        PolymerRecipe.model_validate(recipe_data)


def test_seeded_sequence_is_deterministic_and_forces_center_nhs():
    """The same seed should generate the same centered NHS reactive residue."""
    recipe = sbma_egpma_nhs_recipe(length=9, seed=123)

    first = recipe.generate_sequence()
    second = recipe.generate_sequence()

    assert first == second
    assert len(first) == 9
    assert first[recipe.effective_reactive_index] == "C"


def test_forced_reactive_index_places_nhs_at_known_residue():
    """A configured reactive index should override the default centered NHS."""
    recipe = sbma_egpma_nhs_recipe(length=8, seed=11, reactive_monomer_index=2)
    sequence = recipe.generate_sequence()

    assert recipe.effective_reactive_index == 2
    assert sequence[2] == "C"


def test_fixed_acb_sequence_maps_to_sbma_nhs_egpma():
    """The v1 deterministic ACB recipe should map to SBMA:NHS:EGPMA."""
    recipe = sbma_nhs_egpma_acb_recipe()

    assert recipe.generate_sequence() == "ACB"
    assert recipe.length == 3
    assert recipe.effective_reactive_index == 1
    assert [recipe.monomer_by_label(label).name for label in recipe.generate_sequence()] == [
        "SBMA",
        "NHS",
        "EGPMA",
    ]


def test_fixed_sequence_validates_length_and_labels():
    """Fixed sequences should reject wrong lengths and undeclared labels."""
    with pytest.raises(ValidationError, match="length must match"):
        sbma_egpma_nhs_recipe(length=3, reactive_monomer_index=1, fixed_sequence="AC")

    with pytest.raises(ValidationError, match="declared monomer labels"):
        sbma_egpma_nhs_recipe(length=3, reactive_monomer_index=1, fixed_sequence="AZB")

    with pytest.raises(ValidationError, match="reactive monomer label"):
        sbma_egpma_nhs_recipe(length=3, reactive_monomer_index=1, fixed_sequence="ABC")


def test_invalid_forced_label_and_index_are_rejected():
    """Forced reactive selectors must refer to declared residues in range."""
    recipe_data = _recipe_dict()
    recipe_data["reactive_monomer_label"] = "Z"

    with pytest.raises(ValidationError, match="Reactive monomer label"):
        PolymerRecipe.model_validate(recipe_data)

    recipe_data = _recipe_dict()
    recipe_data["reactive_monomer_index"] = 99

    with pytest.raises(ValidationError, match="within the polymer length"):
        PolymerRecipe.model_validate(recipe_data)


def test_simulation_config_accepts_yaml_polymer_recipe(tmp_path):
    """Conjugation moieties should parse polymer recipes from YAML."""
    data = _minimal_simulation_config_data()
    data["conjugation"] = {
        "enabled": True,
        "mode": "construct",
        "attachments": [
            {
                "name": "lys23-sbma-egpma-nhs",
                "site": {"chain_id": "A", "residue_name": "LYS", "residue_number": 23},
                "moiety": {"name": "SBMA-EGPMA-NHS", "polymer_recipe": _recipe_dict()},
                "mechanism": {"name": "nhs_lys_amide"},
            }
        ],
    }

    config = SimulationConfig.model_validate(data)

    recipe = config.conjugation.attachments[0].moiety.polymer_recipe
    assert isinstance(recipe, PolymerRecipe)
    assert recipe.monomer_by_label("C").name == "NHS"

    config_path = tmp_path / "config.yaml"
    config_path.write_text("""
name: conjugation-polymer-recipe
engine: openmm
enzyme:
  name: enzyme
  pdb_path: enzyme.pdb
thermodynamics:
  temperature: 300.0
simulation_phases:
  equilibration_stages:
    - name: eq
      duration: 0.1
      samples: 10
      ensemble: NVT
      temperature: 300.0
  production:
    ensemble: NPT
    duration: 1.0
    samples: 10
    time_step: 2.0
    report_interval: 1000
    checkpoint_interval: 5000
conjugation:
  enabled: true
  mode: construct
  attachments:
    - name: lys23-sbma-egpma-nhs
      site:
        chain_id: A
        residue_name: LYS
        residue_number: 23
      moiety:
        name: SBMA-EGPMA-NHS
        recipe:
          name: SBMA-EGPMA-NHS
          length: 9
          seed: 7
          forced_reactive_monomer_label: C
          monomers:
            - label: A
              name: SBMA
              residue_name: SBM
              smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])[N+](C([H])([H])[H])(C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])S(=O)(=O)[O-])C([H])([H])[H]"
              probability: 0.945
            - label: B
              name: EGPMA
              residue_name: EGP
              smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])Oc1c([H])c([H])c([H])c([H])c1[H])C([H])([H])[H]"
              probability: 0.045
            - label: C
              name: NHS
              residue_name: NHS
              smiles: CC(=C)C(=O)ON1C(=O)CCC1=O
              probability: 0.01
      mechanism:
        name: nhs_lys_amide
""".strip())

    loaded = load_config(config_path)
    loaded_recipe = loaded.conjugation.attachments[0].moiety.polymer_recipe
    assert isinstance(loaded_recipe, PolymerRecipe)
    assert loaded_recipe.generate_sequence()[loaded_recipe.effective_reactive_index] == "C"


def test_real_polymerist_generation_smoke(tmp_path):
    """Polymerist should consume the real SBMA/EGPMA/NHS recipe when installed."""
    recipe = sbma_egpma_nhs_recipe(length=3, seed=5, reactive_monomer_index=1)

    try:
        result = generate_polymerist_smoke_polymer(recipe, tmp_path / "polymerist", max_retries=1)
    except Exception as exc:
        pytest.skip(f"Polymerist generation stack unavailable in this environment: {exc}")

    assert result.sequence[1] == "C"
    assert result.monomer_group_path.exists()
    assert result.pdb_path is not None and result.pdb_path.exists()
    assert result.pdb_path.with_suffix(".sdf").exists()
    assert result.object_type
    assert result.atom_count is None or result.atom_count > 0


def test_polymerist_to_rdkit_adds_explicit_hydrogens_and_fixes_nitrogen(
    monkeypatch,
):
    """Polymerist conversion should preserve charge fixes before explicit H addition."""
    calls = {}

    class FakeAtom:
        """Minimal atom exposing the RDKit atom methods used by conversion."""

        def __init__(self):
            self.formal_charge = 0

        def GetSymbol(self):
            """Return the atom element symbol."""
            return "N"

        def GetDegree(self):
            """Return the atom degree."""
            return 4

        def GetFormalCharge(self):
            """Return the current formal charge."""
            return self.formal_charge

        def SetFormalCharge(self, charge):
            """Set the formal charge."""
            self.formal_charge = charge

    class FakeMol:
        """Minimal implicit-hydrogen RDKit-like molecule."""

        def __init__(self, other=None):
            self.atom = FakeAtom()
            self.cache_updated = False
            if other is not None:
                self.atom.formal_charge = other.atom.formal_charge

        def GetAtoms(self):
            """Return atoms for formal-charge adjustment."""
            return (self.atom,)

        def GetNumConformers(self):
            """Return one conformer so AddHs should keep coordinates."""
            return 1

        def UpdatePropertyCache(self, *, strict):
            """Record property-cache update calls before adding hydrogens."""
            self.cache_updated = strict is False

    class FakeExplicitMol:
        """Minimal explicit-hydrogen RDKit-like molecule."""

        def __init__(self):
            self.cache_updated = False

        def UpdatePropertyCache(self, *, strict):
            """Record property-cache update calls."""
            self.cache_updated = strict is False

    class FakePolymer:
        """Minimal Polymerist-like object exposing RDKit conversion."""

        def __init__(self):
            self.mol = FakeMol()

        def to_rdkit(self):
            """Return the fake molecule for conversion."""
            return self.mol

    explicit_mol = FakeExplicitMol()

    def fake_add_hs(mol, *, addCoords):
        """Record AddHs arguments and return an explicit-H molecule."""
        calls["mol"] = mol
        calls["addCoords"] = addCoords
        calls["charge_before_add_hs"] = mol.atom.formal_charge
        calls["cache_updated_before_add_hs"] = mol.cache_updated
        return explicit_mol

    fake_chem = ModuleType("rdkit.Chem")
    fake_chem.AddHs = fake_add_hs
    fake_chem.Mol = FakeMol
    fake_rdkit = ModuleType("rdkit")
    fake_rdkit.Chem = fake_chem
    monkeypatch.setitem(sys.modules, "rdkit", fake_rdkit)
    monkeypatch.setitem(sys.modules, "rdkit.Chem", fake_chem)

    polymer = FakePolymer()
    result = _polymerist_to_explicit_h_rdkit_mol(polymer)

    assert result is explicit_mol
    assert calls["addCoords"] is True
    assert calls["charge_before_add_hs"] == 1
    assert calls["cache_updated_before_add_hs"] is True
    assert explicit_mol.cache_updated is True
    assert calls["mol"] is not polymer.mol
    assert polymer.mol.atom.formal_charge == 0


def test_generate_polymerist_smoke_polymer_returns_explicit_rdkit_sidecar(
    monkeypatch,
    tmp_path,
):
    """Smoke generation should return an explicit-H RDKit mol and SDF sidecar."""
    built_sequences = []

    class FakeFragmentGenerator:
        """Fake fragment generator avoiding Polymerist imports."""

        def __init__(self, *, cache_directory, **_kwargs):
            self.cache_directory = cache_directory

        def load_or_generate(self, **_kwargs):
            """Return a placeholder monomer group."""
            return object()

        def get_cache_path(self, recipe_name):
            """Return the expected monomer-group cache path."""
            return self.cache_directory / f"{recipe_name}.json"

    class FakeAtom:
        """Minimal neutral tetravalent nitrogen atom."""

        def __init__(self):
            self.formal_charge = 0

        def GetSymbol(self):
            """Return the element symbol."""
            return "N"

        def GetDegree(self):
            """Return the atom degree."""
            return 4

        def GetFormalCharge(self):
            """Return the formal charge."""
            return self.formal_charge

        def SetFormalCharge(self, charge):
            """Set the formal charge."""
            self.formal_charge = charge

    class FakeImplicitMol:
        """Minimal RDKit-like molecule returned by Polymerist."""

        def __init__(self, other=None):
            self.atom = FakeAtom()
            self.cache_updated = False
            if other is not None:
                self.atom.formal_charge = other.atom.formal_charge

        def GetAtoms(self):
            """Return atoms for formal-charge adjustment."""
            return (self.atom,)

        def GetNumConformers(self):
            """Return zero conformers."""
            return 0

        def UpdatePropertyCache(self, *, strict):
            """Record property-cache update calls before adding hydrogens."""
            self.cache_updated = strict is False

    class FakeExplicitMol:
        """Minimal explicit-H RDKit-like molecule returned by AddHs."""

        def __init__(self):
            self.cache_updated = False

        def UpdatePropertyCache(self, *, strict):
            """Record cache update behavior."""
            self.cache_updated = strict is False

        def GetNumAtoms(self):
            """Return an explicit-H atom count."""
            return 8

    class FakePolymer:
        """Fake Polymerist polymer with RDKit conversion."""

        n_particles = 2

        def __init__(self):
            self.mol = FakeImplicitMol()

        def to_rdkit(self):
            """Return the implicit-H fake molecule."""
            return self.mol

    class FakePolymerGenerator:
        """Fake polymer generator avoiding Polymerist structure building."""

        def __init__(self, **_kwargs):
            pass

        def _build_polymer_structure(self, *, sequence, **_kwargs):
            """Write a PDB and return a fake polymer object."""
            built_sequences.append(sequence)
            pdb_path = tmp_path / "polymer.pdb"
            pdb_path.write_text("END\n")
            return FakePolymer(), pdb_path

    explicit_mol = FakeExplicitMol()

    def fake_add_hs(mol, *, addCoords):
        """Return the explicit-H molecule."""
        assert addCoords is False
        assert mol.cache_updated is True
        return explicit_mol

    def fake_mol_to_mol_file(mol, path):
        """Write the fake SDF sidecar."""
        assert mol is explicit_mol
        Path(path).write_text("fake sdf\n")

    fake_fragment_module = ModuleType("polyzymd.builders.fragment_generator")
    fake_fragment_module.FragmentGenerator = FakeFragmentGenerator
    fake_polymer_module = ModuleType("polyzymd.builders.polymer_generator")
    fake_polymer_module.PolymerGenerator = FakePolymerGenerator
    fake_reactions_module = ModuleType("polyzymd.data.reactions")
    fake_reactions_module.get_atrp_reaction_paths = lambda: {
        "initiation": tmp_path / "initiation.rxn",
        "polymerization": tmp_path / "polymerization.rxn",
        "termination": tmp_path / "termination.rxn",
    }
    fake_chem = ModuleType("rdkit.Chem")
    fake_chem.AddHs = fake_add_hs
    fake_chem.Mol = FakeImplicitMol
    fake_chem.MolToMolFile = fake_mol_to_mol_file
    fake_rdkit = ModuleType("rdkit")
    fake_rdkit.Chem = fake_chem
    monkeypatch.setitem(sys.modules, "polyzymd.builders.fragment_generator", fake_fragment_module)
    monkeypatch.setitem(sys.modules, "polyzymd.builders.polymer_generator", fake_polymer_module)
    monkeypatch.setitem(sys.modules, "polyzymd.data.reactions", fake_reactions_module)
    monkeypatch.setitem(sys.modules, "rdkit", fake_rdkit)
    monkeypatch.setitem(sys.modules, "rdkit.Chem", fake_chem)

    recipe = sbma_nhs_egpma_acb_recipe()
    result = generate_polymerist_smoke_polymer(recipe, tmp_path / "cache")

    assert built_sequences == ["ACB"]
    assert result.rdkit_mol is explicit_mol
    assert result.sdf_path == tmp_path / "polymer.sdf"
    assert result.sdf_path.exists()
    assert result.atom_count == 8
    assert "rdkit_mol" not in result.model_dump()


def test_polymer_generator_forwards_energy_minimize_flag(monkeypatch, tmp_path):
    """Polymer generation should let callers skip Polymerist minimization."""
    from polyzymd.builders import polymer_generator as polymer_generator_module

    calls = []

    class FakeMonomerGroup:
        def __init__(self, monomers):
            self.monomers = monomers
            self.term_orient = {}

    class FakeSourceGroup:
        monomers = {"A_1-site": object(), "A_2-site": object()}

    def fake_build_linear_polymer(**kwargs):
        calls.append(kwargs)
        return object()

    def fake_write_pdb(path, chain, *, resname_map=None):
        path.write_text("END\n")

    monkeypatch.setattr(polymer_generator_module, "MonomerGroup", FakeMonomerGroup)
    monkeypatch.setattr(
        polymer_generator_module,
        "build_linear_polymer",
        fake_build_linear_polymer,
    )
    monkeypatch.setattr(polymer_generator_module, "mbmol_to_rdmol", lambda chain: object())
    monkeypatch.setattr(polymer_generator_module, "summarize_ring_piercing", lambda mol: {})
    monkeypatch.setattr(polymer_generator_module, "mbmol_to_openmm_pdb", fake_write_pdb)

    generator = polymer_generator_module.PolymerGenerator.__new__(
        polymer_generator_module.PolymerGenerator
    )
    generator.monomer_group = FakeSourceGroup()
    generator.cache_directory = tmp_path
    generator.max_retries = 1

    _, pdb_path = generator._build_polymer_structure(
        sequence="AAA",
        monomer_names={"A": "A"},
        energy_minimize=False,
    )

    assert pdb_path.exists()
    assert calls[0]["energy_minimize"] is False


def test_rdkit_sdf_sidecar_write_failure_surfaces(monkeypatch, tmp_path):
    """SDF sidecar writer failures should fail clearly instead of being swallowed."""

    class FakeMol:
        """Minimal RDKit-like molecule for the writer boundary."""

        def GetAtoms(self):
            """Return fake atoms for formal-charge adjustment."""
            return ()

    def fail_mol_to_mol_file(_mol, _path):
        """Raise a representative writer failure."""
        raise OSError("mock SDF writer failure")

    fake_chem = ModuleType("rdkit.Chem")
    fake_chem.MolToMolFile = fail_mol_to_mol_file
    fake_rdkit = ModuleType("rdkit")
    fake_rdkit.Chem = fake_chem
    monkeypatch.setitem(sys.modules, "rdkit", fake_rdkit)
    monkeypatch.setitem(sys.modules, "rdkit.Chem", fake_chem)

    with pytest.raises(RuntimeError, match="Failed to write required RDKit SDF sidecar"):
        _write_rdkit_sdf_sidecar(FakeMol(), tmp_path / "polymer.sdf")
