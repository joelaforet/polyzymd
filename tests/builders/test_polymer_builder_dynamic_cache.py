"""Tests for PolymerBuilder dynamic cache routing."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from polyzymd.builders.polymer import PolymerBuilder
from polyzymd.builders.polymer_generator import PolymerGenerator

#: Star-sited repeat units of the two test monomers.
SBMA = "[*]CC([*])(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]"
EGPMA = "[*]CC([*])(C)C(=O)OCCOc1ccccc1"

#: One fixed hydrogen cap (label H) on each end.
H_CAP_OPTION = {"label": "H", "smiles": "[*][H]", "probability": 1.0}
H_CAPS = {"head": [H_CAP_OPTION], "tail": [H_CAP_OPTION]}

SAGE = ["openff-2.3.0.offxml"]


def _dynamic_builder(tmp_path: Path, **keys: object) -> PolymerBuilder:
    """A dynamic SBMA/EGPMA builder whose cache is ``tmp_path / 'cache'``.

    Parameters
    ----------
    tmp_path : Path
        Test directory.
    **keys : object
        PolymerBuilder arguments that replace the defaults.

    Returns
    -------
    PolymerBuilder
        Builder in dynamic mode.
    """
    arguments = {
        "characters": ["A", "B"],
        "probabilities": [0.0, 1.0],
        "length": 5,
        "type_prefix": "PEGMA",
        "cache_directory": tmp_path / "cache",
        "generation_mode": "dynamic",
        "monomer_smiles": {"SBMA": SBMA, "EGPMA": EGPMA},
        "monomer_names": {"A": "SBMA", "B": "EGPMA"},
        "residue_names": {"SBMA": "SBM", "EGPMA": "EGP"},
        "end_groups": H_CAPS,
        "force_field": SAGE,
    }
    arguments.update(keys)
    return PolymerBuilder(**arguments)


def _write_metadata(generator: PolymerGenerator, sdf_path: Path, sequence: str) -> None:
    """Write the sidecar ``generator`` expects next to ``sdf_path``.

    Parameters
    ----------
    generator : PolymerGenerator
        Generator whose inputs the sidecar records.
    sdf_path : Path
        The charged SDF.
    sequence : str
        Canonical sequence of the chain.
    """
    generator.metadata_path(sdf_path).write_text(
        json.dumps(generator.cache_metadata(sequence)), encoding="utf-8"
    )


@pytest.mark.parametrize(
    ("missing", "message"),
    [
        ({"end_groups": None}, "requires end_groups"),
        ({"force_field": None}, "requires the polymer force_field"),
        ({"monomer_smiles": None}, "requires monomer_smiles"),
    ],
)
def test_dynamic_builder_requires_its_inputs(tmp_path: Path, missing, message) -> None:
    """Dynamic mode needs the repeat units, the caps and the force field."""
    with pytest.raises(ValueError, match=message):
        _dynamic_builder(tmp_path, **missing)


@pytest.mark.parametrize("length", [1, 2])
def test_dynamic_builder_refuses_chains_without_a_repeat_unit(tmp_path: Path, length: int) -> None:
    """The length counts the two caps, so a dynamic chain is at least 3 long."""
    with pytest.raises(ValueError, match="requires length >= 3"):
        _dynamic_builder(tmp_path, length=length)


def test_dynamic_builder_accepts_one_repeat_unit(tmp_path: Path) -> None:
    assert _dynamic_builder(tmp_path, length=3).length == 3


def test_dynamic_builder_creates_the_generator_from_its_settings(tmp_path: Path) -> None:
    """The generator gets the monomers, residue names, caps, force field and charger."""
    builder = _dynamic_builder(tmp_path, max_retries=3, nagl_model="openff-gnn-am1bcc-1.0.0.pt")
    builder._ensure_generators_initialized()
    generator = builder._polymer_generator

    assert {label: m.residue_name for label, m in generator.monomers.items()} == {
        "A": "SBM",
        "B": "EGP",
    }
    assert generator.monomers["A"].unit.smiles == SBMA
    assert {
        end: [(label, cap.unit.smiles) for label, cap in caps.items()]
        for end, caps in generator.end_groups.items()
    } == {"head": [("H", "[*][H]")], "tail": [("H", "[*][H]")]}
    assert generator.force_field == SAGE
    assert generator.max_retries == 3
    assert generator.nagl_model == "openff-gnn-am1bcc-1.0.0.pt"
    assert generator.cache_directory == tmp_path / "cache"


def test_dynamic_builder_routes_cache_directory_hits_through_generator(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Dynamic mode should not load cache-directory SDFs without metadata validation."""
    cache_path = tmp_path / "PEGMA_seq=HBBBH_5-mer_charged.sdf"
    cache_path.write_text("stale cache", encoding="utf-8")
    sentinel = object()
    builder = _dynamic_builder(tmp_path, cache_directory=tmp_path)

    def fail_load_cache(path: Path) -> object:
        """Fail if PolymerBuilder bypasses dynamic cache validation.

        Parameters
        ----------
        path : Path
            Unexpected cache path.

        Returns
        -------
        object
            Never returned because the test fails first.
        """
        pytest.fail(f"Dynamic builder loaded cache directly: {path}")

    monkeypatch.setattr(builder, "_load_from_sdf", fail_load_cache)
    monkeypatch.setattr(builder, "_generate_polymer", lambda sequence: sentinel)

    assert builder._get_or_create_molecule("HBBBH") is sentinel


def test_dynamic_builder_hands_cache_hits_to_the_generator(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """_generate_polymer asks the generator, which decides whether its cache is valid."""
    builder = _dynamic_builder(tmp_path)
    builder._ensure_generators_initialized()
    sentinel = object()
    calls = []

    def fake_generate(sequence: str) -> object:
        calls.append(sequence)
        return sentinel

    monkeypatch.setattr(builder._polymer_generator, "generate_polymer", fake_generate)

    assert builder._get_or_create_molecule("HBBBH") is sentinel
    assert calls == ["HBBBH"]


def test_dynamic_builder_routes_sdf_directory_generated_artifact_through_generator(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Dynamic mode should not directly load generated-style SDF directory hits."""
    sdf_path = tmp_path / "PEGMA_seq=HBBBH_5-mer_charged.sdf"
    sdf_path.write_text("stale generated artifact", encoding="utf-8")
    sentinel = object()
    builder = _dynamic_builder(tmp_path, sdf_directory=tmp_path)

    def fail_load_cache(path: Path) -> object:
        """Fail if generated-style SDF loading bypasses validation.

        Parameters
        ----------
        path : Path
            Unexpected cache path.

        Returns
        -------
        object
            Never returned because the test fails first.
        """
        pytest.fail(f"Dynamic builder loaded generated-style SDF directly: {path}")

    monkeypatch.setattr(builder, "_load_from_sdf", fail_load_cache)
    monkeypatch.setattr(builder, "_generate_polymer", lambda sequence: sentinel)

    assert builder._get_or_create_molecule("HBBBH") is sentinel


def test_dynamic_builder_loads_generated_sdf_directory_hit_when_metadata_matches(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Dynamic mode should load generated-style SDF hits with matching metadata."""
    sdf_path = tmp_path / "PEGMA_seq=HBBBH_5-mer_charged.sdf"
    sdf_path.write_text("validated generated artifact", encoding="utf-8")
    sentinel = object()
    builder = _dynamic_builder(tmp_path, sdf_directory=tmp_path)
    builder._ensure_generators_initialized()
    _write_metadata(builder._polymer_generator, sdf_path, "HBBBH")

    monkeypatch.setattr(builder, "_load_from_sdf", lambda path: sentinel)
    monkeypatch.setattr(
        builder,
        "_generate_polymer",
        lambda sequence: pytest.fail("Validated generated SDF was not loaded"),
    )

    assert builder._get_or_create_molecule("HBBBH") is sentinel


@pytest.mark.parametrize(
    "stale",
    [
        {"monomer_smiles": {"SBMA": SBMA, "EGPMA": "[*]CC([*])(C)C(=O)OCCO"}},
        {"end_groups": {"head": [{**H_CAP_OPTION, "smiles": "[*]C"}], "tail": [H_CAP_OPTION]}},
        {"residue_names": {"SBMA": "SBM", "EGPMA": "EGM"}},
        {"charger_type": "am1bcc"},
        {"force_field": ["openff-2.2.0.offxml"]},
    ],
    ids=["smiles", "end_groups", "residue_name", "charger", "force_field"],
)
def test_dynamic_builder_regenerates_generated_sdf_directory_hit_with_stale_metadata(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    stale: dict[str, object],
) -> None:
    """Dynamic mode should regenerate generated-style SDF hits with stale metadata."""
    sdf_path = tmp_path / "PEGMA_seq=HBBBH_5-mer_charged.sdf"
    sdf_path.write_text("stale generated artifact", encoding="utf-8")
    old = _dynamic_builder(tmp_path, cache_directory=tmp_path / "stale", **stale)
    old._ensure_generators_initialized()
    _write_metadata(old._polymer_generator, sdf_path, "HBBBH")
    sentinel = object()
    builder = _dynamic_builder(tmp_path, sdf_directory=tmp_path)

    monkeypatch.setattr(
        builder,
        "_load_from_sdf",
        lambda path: pytest.fail("Stale generated SDF was loaded directly"),
    )
    monkeypatch.setattr(builder, "_generate_polymer", lambda sequence: sentinel)

    assert builder._get_or_create_molecule("HBBBH") is sentinel


def test_metadata_of_unused_monomers_does_not_invalidate_the_cache(tmp_path: Path) -> None:
    """A chain's sidecar records only the monomers its sequence uses."""
    sdf_path = tmp_path / "PEGMA_seq=HBBBH_5-mer_charged.sdf"
    sdf_path.write_text("generated artifact", encoding="utf-8")
    old = _dynamic_builder(tmp_path, monomer_smiles={"SBMA": "[*]CC([*])C", "EGPMA": EGPMA})
    old._ensure_generators_initialized()
    _write_metadata(old._polymer_generator, sdf_path, "HBBBH")
    builder = _dynamic_builder(tmp_path, sdf_directory=tmp_path)

    assert builder._dynamic_sdf_metadata_matches(sdf_path, "HBBBH")
    assert set(builder._polymer_generator.cache_metadata("HBBBH")["monomers"]) == {"B"}


@pytest.mark.parametrize("sidecar", [None, "not json", "[1, 2]"])
def test_a_missing_or_unreadable_sidecar_does_not_match(tmp_path: Path, sidecar) -> None:
    """Without a readable sidecar object the generator builds the chain again."""
    builder = _dynamic_builder(tmp_path)
    builder._ensure_generators_initialized()
    generator = builder._polymer_generator
    sdf_path = generator.sdf_path("HBBBH")
    sdf_path.parent.mkdir(parents=True)
    sdf_path.write_text("generated artifact", encoding="utf-8")
    if sidecar is not None:
        generator.metadata_path(sdf_path).write_text(sidecar, encoding="utf-8")

    assert not generator.cache_metadata_matches(sdf_path, "HBBBH")


def test_dynamic_builder_loads_a_user_sdf_as_given(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A user SDF in sdf_directory that is not a cache entry is loaded without a sidecar."""
    sdf_directory = tmp_path / "user"
    sdf_directory.mkdir()
    builder = _dynamic_builder(tmp_path, sdf_directory=sdf_directory)
    user_sdf = builder._get_sdf_path("HBBBH", sdf_directory)
    user_sdf.write_text("user molecule", encoding="utf-8")
    sentinel = object()
    monkeypatch.setattr(builder, "_is_dynamic_cache_like_sdf", lambda path, sequence: False)
    monkeypatch.setattr(builder, "_load_from_sdf", lambda path: sentinel)
    monkeypatch.setattr(
        builder,
        "_generate_polymer",
        lambda sequence: pytest.fail("A user SDF was rebuilt"),
    )

    assert builder._get_or_create_molecule("HBBBH") is sentinel


def test_cached_builder_loads_sdf_directory_hit_without_dynamic_metadata(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Cached mode should continue to load user-provided SDF directory hits."""
    sdf_path = tmp_path / "EGMA_seq=A_1-mer_charged.sdf"
    sdf_path.write_text("curated user structure", encoding="utf-8")
    sentinel = object()
    builder = PolymerBuilder(
        characters=["A"],
        probabilities=[1.0],
        length=1,
        type_prefix="EGMA",
        sdf_directory=tmp_path,
        cache_directory=tmp_path / "cache",
        generation_mode="cached",
    )

    monkeypatch.setattr(builder, "_load_from_sdf", lambda path: sentinel)

    assert builder._get_or_create_molecule("A") is sentinel


def test_cached_mode_finds_the_file_dynamic_mode_writes(tmp_path: Path) -> None:
    """Both modes name a chain's SDF after type_prefix, so cached mode reuses a dynamic cache."""
    dynamic = _dynamic_builder(tmp_path, cache_directory=tmp_path, type_prefix="SE")
    dynamic._ensure_generators_initialized()
    written = dynamic._polymer_generator.sdf_path("HABAH")
    cached = PolymerBuilder(
        characters=["A", "B"],
        probabilities=[0.5, 0.5],
        length=5,
        type_prefix="SE",
        sdf_directory=tmp_path,
        generation_mode="cached",
    )
    assert cached._get_sdf_path("HABAH", tmp_path) == written


def test_missing_cached_sdf_error_names_real_config_keys(tmp_path: Path) -> None:
    """The error for a missing cached chain points at config keys that exist."""
    cached = PolymerBuilder(
        characters=["A"],
        probabilities=[1.0],
        length=5,
        type_prefix="SE",
        sdf_directory=tmp_path,
        generation_mode="cached",
    )
    with pytest.raises(FileNotFoundError) as error:
        cached._get_or_create_molecule("AAAAA")
    assert "allow_generation" not in str(error.value)
    assert "generation_mode: dynamic" in str(error.value)


def test_molecule_source_reports_preset_for_a_user_sdf_in_dynamic_mode(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A user SDF that dynamic mode loads as given keeps its own charges, so they are preset."""
    sdf_directory = tmp_path / "user"
    sdf_directory.mkdir()
    builder = _dynamic_builder(tmp_path, sdf_directory=sdf_directory)
    user_sdf = builder._get_sdf_path("HBBBH", sdf_directory)
    user_sdf.write_text("user molecule", encoding="utf-8")
    monkeypatch.setattr(builder, "_is_dynamic_cache_like_sdf", lambda path, sequence: False)

    assert builder.molecule_source("HBBBH") == {
        "sdf_path": user_sdf,
        "charge_method": "preset",
        "nagl_model": None,
    }

    monkeypatch.setattr(builder, "_is_dynamic_cache_like_sdf", lambda path, sequence: True)
    assert builder.molecule_source("HBBBH")["charge_method"] == builder._charger_type


# ---------------------------------------------------------------------------
# Cached chains that dynamic generation made
# ---------------------------------------------------------------------------

STAR_MMA = "[*]CC([*])(C)C(=O)OC"
STAR_HEMA = "[*]CC([*])(C)C(=O)OCCO"


#: Cap options: a hydrogen or a methyl head, and a hydrogen tail.
CACHED_CAPS = {
    "head": [
        {"label": "H", "smiles": "[*][H]", "probability": 0.5},
        {"label": "M", "smiles": "[*]C", "probability": 0.5, "name": "methyl"},
    ],
    "tail": [{"label": "H", "smiles": "[*][H]", "probability": 1.0}],
}


def _cached_builder(
    tmp_path: Path, mma: str = STAR_MMA, caps: dict | None = None
) -> PolymerBuilder:
    return PolymerBuilder(
        characters=["A", "B"],
        probabilities=[0.5, 0.5],
        length=4,
        type_prefix="MH",
        sdf_directory=tmp_path,
        monomer_smiles={"MMA": mma, "HEMA": STAR_HEMA},
        monomer_names={"A": "MMA", "B": "HEMA"},
        end_groups=caps if caps is not None else CACHED_CAPS,
    )


def _cap_record(label: str, smiles: str, name: str | None = None) -> dict:
    """The per-chain cap record dynamic generation writes in a sidecar."""
    return {"label": label, "name": name, "residue_name": None, "smiles": smiles}


def _write_sidecar(tmp_path: Path, sequence: str = "MABH") -> Path:
    sdf = tmp_path / f"MH_seq={sequence}_{len(sequence)}-mer_charged.sdf"
    sdf.write_text("chain", encoding="utf-8")
    sidecar = sdf.with_name(f"{sdf.name}.metadata.json")
    head = _cap_record("M", "[*]C", "methyl") if sequence[0] == "M" else _cap_record("H", "[*][H]")
    sidecar.write_text(
        json.dumps(
            {
                "schema_version": 5,
                "sequence": sequence,
                "monomers": {
                    "A": {"name": "MMA", "residue_name": "MMA", "smiles": STAR_MMA},
                    "B": {"name": "HEMA", "residue_name": "HEM", "smiles": STAR_HEMA},
                },
                "end_groups": {"head": head, "tail": _cap_record("H", "[*][H]")},
                "builder": {"package": "mbuild", "version": "1.4.0", "commit": "abc"},
            }
        ),
        encoding="utf-8",
    )
    return sdf


def test_cached_mode_refuses_a_chain_built_from_other_smiles(tmp_path: Path) -> None:
    """A config naming chemistry its cached chains were not built from is refused."""
    sdf = _write_sidecar(tmp_path)
    builder = _cached_builder(tmp_path, mma="[*]CC([*])(C)C(=O)OCC")
    with pytest.raises(ValueError, match="built from other chemistry") as error:
        builder._check_cached_chemistry(sdf)
    assert "monomer 'A'" in str(error.value)


def test_cached_mode_refuses_a_chain_with_other_caps(tmp_path: Path) -> None:
    """The chain's own cap (by label) is compared with the config's option of that label."""
    sdf = _write_sidecar(tmp_path, "MABH")
    caps = {
        "head": [
            {"label": "H", "smiles": "[*][H]", "probability": 0.5},
            {"label": "M", "smiles": "[*]CC", "probability": 0.5},
        ],
        "tail": CACHED_CAPS["tail"],
    }
    builder = _cached_builder(tmp_path, caps=caps)
    with pytest.raises(ValueError, match="end_groups.head 'M'") as error:
        builder._check_cached_chemistry(sdf)
    assert "'[*]CC'" in str(error.value) and "'[*]C'" in str(error.value)

    # A chain whose head is the unchanged hydrogen option still matches.
    builder._check_cached_chemistry(_write_sidecar(tmp_path, "HABH"))


def test_cached_mode_accepts_the_chain_caps_of_its_config(tmp_path: Path) -> None:
    """Each chain records the caps it drew; both draws of the config match."""
    builder = _cached_builder(tmp_path, caps={**CACHED_CAPS})
    builder._check_cached_chemistry(_write_sidecar(tmp_path, "MABH"))
    builder._check_cached_chemistry(_write_sidecar(tmp_path, "HBAH"))


def test_cached_mode_accepts_equivalent_smiles_and_records_the_builder(tmp_path: Path) -> None:
    """The same molecule written differently matches; the chain's builder is recorded."""
    sdf = _write_sidecar(tmp_path)
    builder = _cached_builder(tmp_path, mma="COC(=O)C([*])(C)C[*]")
    builder._check_cached_chemistry(sdf)
    assert builder.chain_builders() == [{"package": "mbuild", "version": "1.4.0", "commit": "abc"}]
    assert builder.chain_chemistry()["chains"][sdf.name]["end_groups"] == {
        "head": _cap_record("M", "[*]C", "methyl"),
        "tail": _cap_record("H", "[*][H]"),
    }


def test_cached_mode_warns_when_it_cannot_check_the_chemistry(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    sdf = tmp_path / "MH_seq=MABH_4-mer_charged.sdf"
    sdf.write_text("chain", encoding="utf-8")
    builder = _cached_builder(tmp_path)
    builder._check_cached_chemistry(sdf)
    assert "not checked" in caplog.text
    assert builder.chain_builders() == []


def test_cached_mode_draws_caps_when_the_config_gives_them(tmp_path: Path) -> None:
    """A cached config with end_groups finds chains named with their cap labels."""
    builder = PolymerBuilder(
        characters=["A"],
        probabilities=[1.0],
        length=4,
        type_prefix="MH",
        sdf_directory=tmp_path,
        end_groups={
            "head": [{"label": "H", "smiles": "[*][H]", "probability": 1.0}],
            "tail": [{"label": "M", "smiles": "[*]C", "probability": 1.0}],
        },
    )
    assert builder._draw_sequence() == "HAAM"


def test_cached_mode_refuses_a_chain_whose_cap_the_config_lacks(tmp_path: Path) -> None:
    sdf = tmp_path / "MH_seq=HABH_4-mer_charged.sdf"
    sdf.write_text("chain", encoding="utf-8")
    sdf.with_name(f"{sdf.name}.metadata.json").write_text(
        json.dumps(
            {
                "schema_version": 5,
                "monomers": {},
                "end_groups": {
                    "head": {"label": "H", "smiles": "[*][H]"},
                    "tail": {"label": "H", "smiles": "[*][H]"},
                },
            }
        ),
        encoding="utf-8",
    )
    builder = PolymerBuilder(
        characters=["A", "B"],
        probabilities=[0.5, 0.5],
        length=4,
        type_prefix="MH",
        sdf_directory=tmp_path,
        end_groups={
            "head": [{"label": "M", "smiles": "[*]C", "probability": 1.0}],
            "tail": [{"label": "H", "smiles": "[*][H]", "probability": 1.0}],
        },
    )
    with pytest.raises(ValueError, match="is not one of the config's caps"):
        builder._check_cached_chemistry(sdf)
