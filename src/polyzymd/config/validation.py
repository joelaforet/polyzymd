"""Runtime validation helpers for simulation configuration references."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def collect_reference_warnings(config: Any) -> list[str]:
    """Return warnings for missing files referenced by a simulation config.

    These checks intentionally live outside the Pydantic schema so lightweight
    configuration parsing can succeed before external structure files are staged.

    Parameters
    ----------
    config : Any
        Simulation configuration or config-like object with PolyzyMD schema
        attributes.

    Returns
    -------
    list[str]
        Human-readable warnings for missing referenced files or directories.
    """

    warnings: list[str] = []
    enzyme = getattr(config, "enzyme", None)
    _check_file(warnings, getattr(enzyme, "pdb_path", None), "enzyme PDB")

    substrate = getattr(config, "substrate", None)
    if substrate is not None:
        _check_file(warnings, getattr(substrate, "sdf_path", None), "substrate SDF")

    polymers = getattr(config, "polymers", None)
    if polymers is not None and bool(getattr(polymers, "enabled", False)):
        _check_polymer_references(warnings, polymers)

    _check_am1bcc(warnings, config)
    return warnings


def _check_am1bcc(warnings: list[str], config: Any) -> None:
    """Warn when the config asks for AM1-BCC charges and no toolkit here computes them."""
    from polyzymd.utils import charging

    substrate = getattr(config, "substrate", None)
    solvent = getattr(config, "solvent", None)
    polymers = getattr(config, "polymers", None)
    methods = [getattr(substrate, "charge_method", None)]
    methods += [cs.charges.method for cs in getattr(solvent, "co_solvents", None) or []]
    if polymers is not None and getattr(polymers, "enabled", False):
        methods.append(getattr(polymers, "charger", None))
    if "am1bcc" in {getattr(m, "value", m) for m in methods} and not charging.am1bcc_available():
        warnings.append(charging.AM1BCC_NEEDS_AMBERTOOLS)


def _check_file(warnings: list[str], path_value: Any, label: str) -> None:
    """Append a warning when a referenced file path is missing."""

    if path_value is None:
        return
    path = Path(path_value)
    if not path.is_file():
        warnings.append(f"Missing {label}: {path}")


def _check_directory(warnings: list[str], path_value: Any, label: str) -> Path | None:
    """Append a warning when a referenced directory path is missing."""

    if path_value is None:
        return None
    path = Path(path_value)
    if not path.is_dir():
        warnings.append(f"Missing {label}: {path}")
        return None
    return path


def _check_polymer_references(warnings: list[str], polymers: Any) -> None:
    """Check the SDF directory of cached polymers."""

    generation_mode = str(getattr(polymers, "generation_mode", "")).lower()
    if generation_mode.endswith("cached"):
        sdf_directory = _check_directory(
            warnings,
            getattr(polymers, "sdf_directory", None),
            "polymer SDF directory",
        )
        if sdf_directory is not None:
            _check_cached_polymer_sdfs(warnings, polymers, sdf_directory)


def _check_cached_polymer_sdfs(warnings: list[str], polymers: Any, sdf_directory: Path) -> None:
    """Warn when a cached-polymer directory has no matching SDF files."""

    type_prefix = getattr(polymers, "type_prefix", None)
    length = getattr(polymers, "length", None)
    if not type_prefix or length is None:
        return

    pattern = f"{type_prefix}_seq=*_{length}-mer_charged.sdf"
    if not any(sdf_directory.glob(pattern)):
        warnings.append(
            "Missing cached polymer SDF files: " f"no files matching {sdf_directory / pattern}"
        )


#: Elements the NAGL charge model handles (the atom features of
#: ``openff-gnn-am1bcc-0.1.0-rc.3.pt``).
NAGL_ELEMENTS = frozenset({"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"})


def _key(config: Any, solute: Any, field: str) -> str:
    """The config key of ``field`` on ``solute``, as the loader names it (``solutes.0.pdb_path``)."""
    index = next(i for i, entry in enumerate(config.solutes) if entry is solute)
    return f"solutes.{index}.{field}"


def require_inputs(config: Any) -> None:
    """Raise ``ValueError`` when an input of a new build is missing or cannot be used.

    Cheap checks of what the build reads: the enzyme PDB has atoms, the
    substrate SDF holds the configured conformer, co-solvent and monomer SMILES
    parse, monomers carry the group the initiation reaction needs, the charge
    method can charge each molecule, and the force fields are installed.
    :meth:`SimulationConfig.require_buildable` calls this.
    """
    pdb = Path(config.enzyme.pdb_path)
    if not pdb.is_file():
        raise ValueError(
            f"protein PDB {pdb} does not exist.\n"
            f"fix: correct {_key(config, config.enzyme, 'pdb_path')}."
        )
    with open(pdb, errors="replace") as stream:
        if not any(line.startswith(("ATOM", "HETATM")) for line in stream):
            raise ValueError(
                f"protein PDB {pdb} has no ATOM or HETATM records.\n"
                "fix: give a PDB file of the prepared protein (polyzymd clean-pdb writes one)."
            )

    from rdkit import Chem

    substrate = config.substrate
    if substrate is not None:
        sdf = Path(substrate.sdf_path)
        if not sdf.is_file():
            raise ValueError(
                f"small-molecule SDF {sdf} does not exist.\n"
                f"fix: correct {_key(config, substrate, 'sdf_path')}."
            )
        molecules = list(Chem.SDMolSupplier(str(sdf), removeHs=False))
        if substrate.conformer_index >= len(molecules):
            raise ValueError(
                f"small-molecule SDF {sdf} holds {len(molecules)} conformer(s), so conformer_index "
                f"{substrate.conformer_index} does not exist.\n"
                f"fix: set conformer_index between 0 and {max(len(molecules) - 1, 0)}."
            )
        molecule = molecules[substrate.conformer_index]
        if molecule is None:
            raise ValueError(
                f"RDKit cannot read conformer {substrate.conformer_index} of {sdf}.\n"
                "fix: write the SDF again from your docking or drawing program."
            )
        _require_chargeable(molecule, substrate.charge_method, f"solute {substrate.name!r}")

    from polyzymd.data.solvent_molecules import is_bundled_solvent, split_counter_ions

    for cosolvent in config.solvent.co_solvents:
        if cosolvent.charges.method == "preset" and is_bundled_solvent(cosolvent.name):
            _require_shipped_molecule(cosolvent)
            continue
        what = f"co-solvent {cosolvent.name!r}"
        smiles, _, _ = split_counter_ions(cosolvent.smiles, cosolvent.name)
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(
                f"{what}: RDKit cannot read the SMILES {cosolvent.smiles!r}.\n"
                "fix: correct the SMILES."
            )

        if molecule.GetNumAtoms() == 1 and molecule.GetAtomWithIdx(0).GetFormalCharge():
            raise ValueError(
                f"{what} is the ion {smiles}; the build adds only Na+ and Cl- ions.\n"
                "fix: remove this co-solvent and use solvent.ions.nacl_concentration."
            )
        _require_chargeable(molecule, cosolvent.charges.method, what)

    polymers = config.polymers
    if polymers is not None and polymers.enabled and polymers.generation_mode.value == "dynamic":
        from polyzymd.builders.polymer_generator import parse_star_smiles, require_mbuild_polymer

        units = [
            (f"monomer {m.label!r}", parse_star_smiles(m.smiles, 2, f"monomer {m.label!r}"))
            for m in polymers.monomers
        ]
        caps = [
            (f"end_groups.{end} {option.label!r}", option)
            for end in ("head", "tail")
            for option in getattr(polymers.end_groups, end)
        ]
        units += [(what, parse_star_smiles(option.smiles, 1, what)) for what, option in caps]
        for what, unit in units:
            _require_chargeable(unit.mol, polymers.charger, what)
        from polyzymd.config.schema import reserved_residue_names

        named = [(f"monomer {m.label!r}", m.residue_name) for m in polymers.monomers]
        named += [(what, option.residue_name) for what, option in caps]
        for what, residue_name in named:
            if residue_name and residue_name.upper() in reserved_residue_names():
                raise ValueError(
                    f"{what} has residue_name '{residue_name}', an amino-acid, nucleic-acid, "
                    "water or ion name, so selections such as 'protein' would include the "
                    "polymer.\n"
                    f"fix: set the residue_name of {what} in {_key(config, polymers, 'monomers')} "
                    "or end_groups to another 3-letter name."
                )
        try:
            require_mbuild_polymer()
        except ImportError as exc:
            raise ValueError(str(exc)) from exc

    from openff.toolkit.typing.engines.smirnoff.forcefield import get_available_force_fields

    available = set(get_available_force_fields())
    named = [(f"solute {s.name!r}", s.force_field) for s in config.solutes]
    named.append(("solvent", config.solvent.force_field))
    for what, files in named:
        for name in files:
            if name not in available and not Path(name).is_file():
                raise ValueError(
                    f"{what}: force_field {name!r} is not an installed force field.\n"
                    "fix: correct the name (for example ff14sb_off_impropers_0.0.4.offxml or "
                    "openff-2.3.0.offxml), or give the path of an .offxml file."
                )

    from polyzymd.utils.charging import installed_nagl_model

    nagl_model = getattr(config, "nagl_model", None)
    if nagl_model is not None and installed_nagl_model(nagl_model) is None:
        from openff.nagl_models import list_available_nagl_models

        installed = ", ".join(sorted(path.name for path in list_available_nagl_models()))
        raise ValueError(
            f"charges model {nagl_model!r} is not an installed NAGL model.\n"
            f"fix: in each charges: {{method: nagl, model: ...}}, name one of {installed}, "
            "or give the path of a model file."
        )


def _require_shipped_molecule(cosolvent: Any) -> None:
    """Raise ``ValueError`` when a preset co-solvent's SMILES is not the molecule PolyzyMD ships."""
    from rdkit import Chem

    from polyzymd.data.solvent_molecules import _SOLVENTS_DIR, split_counter_ions

    library = _SOLVENTS_DIR / f"{cosolvent.name.lower().strip()}.sdf"
    shipped = next(iter(Chem.SDMolSupplier(str(library), removeHs=False)), None)
    smiles, _, _ = split_counter_ions(cosolvent.smiles, cosolvent.name)
    given = Chem.MolFromSmiles(smiles)
    if shipped is None or given is None:
        return
    if Chem.MolToSmiles(Chem.RemoveHs(shipped)) != Chem.MolToSmiles(given):
        raise ValueError(
            f"co-solvent {cosolvent.name!r}: its SMILES {cosolvent.smiles!r} is not the "
            f"{cosolvent.name} PolyzyMD ships charges for, so charges preset cannot apply.\n"
            "fix: remove smiles to use the shipped molecule, or set charges: {method: nagl} "
            "(or am1bcc, espaloma)."
        )


def _require_chargeable(molecule: Any, method: Any, what: str) -> None:
    """Raise ``ValueError`` when ``method`` cannot assign charges to an RDKit molecule."""
    method = getattr(method, "value", method)
    if method == "nagl":
        elements = {atom.GetSymbol() for atom in molecule.GetAtoms()}
        others = sorted(elements - NAGL_ELEMENTS)
        if others:
            raise ValueError(
                f"{what} contains {', '.join(others)}, which the NAGL charge model does not "
                f"cover (it covers {', '.join(sorted(NAGL_ELEMENTS))}).\n"
                "fix: set charges: {method: am1bcc} (needs AmberTools), or leave the molecule "
                "out."
            )
    elif method == "am1bcc":
        from polyzymd.utils import charging

        if not charging.am1bcc_available():
            raise ValueError(
                f"{what}: charges method am1bcc needs AmberTools or OpenEye, and neither is "
                "installed in this environment.\nfix: use charges: {method: nagl}."
            )
    elif method == "espaloma":
        import importlib.util

        if importlib.util.find_spec("espaloma_charge") is None:
            raise ValueError(
                f"{what}: charges method espaloma needs the espaloma-charge package, which is not "
                "installed in this environment.\nfix: use charges: {method: nagl}."
            )
