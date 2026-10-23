"""A real build: its components, and the system each engine will run read back the same.

Builds the small all-components system once (trp-cage, substrate, two-monomer
polymers, DMSO, water, NaCl), exports it through PolyzyMD's own OpenMM and
GROMACS paths, and reads both back as realized systems.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

pytest.importorskip("openff.interchange")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _support.small_system import build  # noqa: E402

from polyzymd.audit.realized import RealizedSystem, parse_mdp  # noqa: E402
from polyzymd.audit.terms import TERM_KINDS, compare  # noqa: E402
from polyzymd.utils.packmol import packing_order  # noqa: E402

#: Replicate 2 draws BB x3, AA x1, AB x2, which Packmol places as BB, AB, AA:
#: chain order in the system differs from the order the sequences were drawn.
REPLICATE = 2


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    folder = tmp_path_factory.mktemp("small_system")
    builder, config = build(folder, replicate=REPLICATE)
    return builder, config, folder


@pytest.fixture(scope="module")
def engines(built):
    """The small system as OpenMM and GROMACS will run it."""
    import openmm

    from polyzymd.exporters.gromacs import GromacsExporter

    builder, config, folder = built
    interchange = builder.create_interchange()
    xml = folder / "system.xml"
    xml.write_text(
        openmm.XmlSerializer.serialize(interchange.to_openmm(combine_nonbonded_forces=True))
    )
    pdb = folder / "build" / "solvated_system.pdb"
    openmm_side = RealizedSystem.from_openmm(xml, pdb, timestep_fs=2.0)
    if shutil.which("gmx") is None and not Path.home().joinpath(".pixi/bin/gmx").exists():
        pytest.skip("GROMACS is not installed")
    files = GromacsExporter(interchange, config, builder.get_component_info()).export(
        folder / "gromacs"
    )
    gromacs_side = RealizedSystem.from_gromacs(files["top"], files["gro"], files["prod_mdp"])
    return openmm_side, gromacs_side


def test_components_cover_the_system_in_build_order(built) -> None:
    """Every atom belongs to one component, in the order the builder placed them."""
    builder, _, _ = built
    components = builder.components()
    assert [(c.name, c.type, c.segid, c.placement) for c in components] == [
        ("protein", "protein", "PROT", "fixed"),
        ("substrate", "small_molecule", "LIG", "fixed"),
        ("polymer", "polymer", "POLY", "packed"),
        ("solvent", "solvent", "SOL", "packed"),
        ("ions", "ions", "ION", "packed"),
    ]
    assert sum(c.n_atoms for c in components) == builder.solvated_topology.n_atoms
    info = builder.get_component_info()
    sizes = {c.name: c.n_atoms for c in components}
    assert (sizes["protein"], sizes["substrate"], sizes["polymer"]) == (
        info.n_protein_atoms,
        info.n_substrate_atoms,
        info.n_polymer_atoms,
    )


def test_each_component_names_its_parameter_and_charge_sources(built) -> None:
    """The records say what the config asked for, per domain."""
    builder, config, _ = built
    components = {c.name: c for c in builder.components()}
    protein = components["protein"].domains[0]
    assert protein.force_field == (config.force_field.protein,)
    assert protein.charges.method == "library"
    substrate = components["substrate"].domains[0]
    assert substrate.force_field == (config.force_field.small_molecule,)
    assert (substrate.charges.method, substrate.charges.model) == (
        "nagl",
        config.force_field.nagl_model,
    )
    polymer = components["polymer"]
    assert polymer.n_molecules == config.polymers.count
    for domain in polymer.domains:
        sequence = domain.name.split(":")[1]
        assert domain.charges.method == "preset"
        assert domain.charges.source_file.endswith(f"TP_seq={sequence}_2-mer_charged.sdf")
    solvent = {d.name: d for d in components["solvent"].domains}
    assert solvent["water"].charges.method == "water_model"
    assert solvent["co_solvent:dmso"].charges.method == "preset"
    assert solvent["co_solvent:dmso"].charges.source_file.endswith("dmso.sdf")
    assert components["protein"].provenance["pdb"]["sha256"]


def test_each_polymer_domain_holds_its_own_sequence(built) -> None:
    """A domain's chains have the atoms of the SDF it names, so sources are not swapped."""
    from openff.toolkit import Molecule

    builder, _, _ = built
    counts = builder._polymer_counts
    assert packing_order(counts) != list(range(len(counts))), "test needs reordered chains"
    topology = builder.solvated_topology
    elements = [atom.symbol for atom in topology.atoms]
    polymer = next(c for c in builder.components() if c.name == "polymer")
    for domain in polymer.domains:
        reference = Molecule.from_file(domain.charges.source_file)
        expected = sorted(atom.symbol for atom in reference.atoms)
        for chain in domain.atoms:
            assert sorted(elements[i] for i in chain) == expected


def test_both_engines_read_back_the_same_system(engines) -> None:
    """PolyzyMD's OpenMM and GROMACS exports hold the same terms, atom for atom."""
    openmm_side, gromacs_side = engines
    assert openmm_side.n_atoms == gromacs_side.n_atoms
    everyone = range(openmm_side.n_atoms)
    left = openmm_side.index.extract(everyone)
    right = gromacs_side.index.extract(everyone)
    for kind in TERM_KINDS:
        result = compare(kind, left[kind], right[kind])
        assert result.ok, f"{kind}: {result.matched}/{result.total} agree, e.g. {result.mismatches}"


def test_both_engines_run_with_the_same_settings(engines) -> None:
    """Constraints, cutoff, switching and dispersion agree; GROMACS reads them from the .mdp."""
    openmm_side, gromacs_side = engines
    assert openmm_side.settings["constraints"] == gromacs_side.settings["constraints"] == "h-bonds"
    for key in ("nonbonded_method", "dispersion_correction"):
        assert openmm_side.settings[key] == gromacs_side.settings[key], key
    for key in ("cutoff_nm", "switch_nm"):
        assert openmm_side.settings[key] == pytest.approx(gromacs_side.settings[key]), key
    assert gromacs_side.settings["timestep_fs"] == pytest.approx(2.0)
    assert any(path.endswith(".itp") for path in gromacs_side.files)


def test_mdp_settings_are_read_whatever_their_spelling(tmp_path) -> None:
    """Keys are case- and separator-insensitive, and comments are dropped."""
    mdp = tmp_path / "prod.mdp"
    mdp.write_text(
        "; header\nconstraints = h-bonds ; bonds to H\nDispCorr = EnerPres\nrvdw_switch=0.8\n"
    )
    assert parse_mdp(mdp) == {
        "constraints": "h-bonds",
        "dispcorr": "EnerPres",
        "rvdw-switch": "0.8",
    }
